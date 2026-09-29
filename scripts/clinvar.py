"""Match the genome's SNV calls against ClinVar (GRCh37 SNVs) by chrom + VCF position.

Reads data/genome_23andme.txt and reference/variant_summary.txt.gz.
Writes results/clinvar_hits.tsv and results/clinvar_summary.json.
"""

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
CLINVAR = ROOT / "reference" / "variant_summary.txt.gz"
GENOME = ROOT / "data" / "genome_23andme.txt"
COMP = str.maketrans("ACGT", "TGCA")

COLS = [
    "#AlleleID",
    "Type",
    "Name",
    "GeneSymbol",
    "ClinicalSignificance",
    "RS# (dbSNP)",
    "PhenotypeList",
    "Assembly",
    "Chromosome",
    "ReviewStatus",
    "NumberSubmitters",
    "VariationID",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "OriginSimple",
]

def selected() -> list[str]:
    """Well-known markers reported whatever the genotype (config: selected_markers)."""
    from run_config import get
    return get("selected_markers")

STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
}


def category(sig: str) -> str:
    s = sig.lower()
    if "conflicting" in s:
        return "Conflicting"
    if re.search(r"(^|[/;, ])(likely )?pathogenic", s) and "benign" not in s:
        return "Pathogenic / likely pathogenic"
    if "risk factor" in s or "risk allele" in s:
        return "Risk factor"
    if "drug response" in s:
        return "Drug response"
    if "protective" in s:
        return "Protective"
    if "association" in s or "affects" in s:
        return "Association"
    if "uncertain" in s:
        return "Uncertain significance"
    if "benign" in s:
        return "Benign / likely benign"
    return "Other / not provided"


def load_clinvar() -> pd.DataFrame:
    cv = pd.read_csv(CLINVAR, sep="\t", usecols=COLS, dtype=str, low_memory=False)
    cv = cv[(cv["Assembly"] == "GRCh37") & (cv["Type"] == "single nucleotide variant")]
    cv = cv[cv["PositionVCF"].str.fullmatch(r"\d+", na=False)]
    cv = cv[
        cv["ReferenceAlleleVCF"].isin(list("ACGT"))
        & cv["AlternateAlleleVCF"].isin(list("ACGT"))
    ]
    cv["pos"] = cv["PositionVCF"].astype(int)
    cv = cv.rename(
        columns={
            "Chromosome": "chrom",
            "ReferenceAlleleVCF": "ref",
            "AlternateAlleleVCF": "alt",
        }
    )
    return cv


def main() -> None:
    g = pd.read_csv(
        GENOME,
        sep="\t",
        comment="#",
        names=["rsid", "chrom", "pos", "gt"],
        dtype={"chrom": str},
    )
    g = g[g["gt"] != "--"]
    cv = load_clinvar()
    m = g.merge(cv, on=["chrom", "pos"], how="inner")

    # Build/strand check per site: are the called alleles within ClinVar's ref + alts at that position?
    site_alleles = m.groupby(["chrom", "pos"]).agg(
        gt=("gt", "first"), ref=("ref", "first"), alts=("alt", lambda a: set(a))
    )
    fwd = rc = neither = 0
    bad_sites = set()
    for (chrom, pos), (gt, ref, alts) in zip(site_alleles.index, site_alleles.itertuples(index=False)):
        allowed = alts | {ref}
        if set(gt) <= allowed:
            fwd += 1
        elif set(gt.translate(COMP)) <= allowed:
            rc += 1
            bad_sites.add((chrom, pos))
        else:
            neither += 1
            bad_sites.add((chrom, pos))
    ref_present = int(
        sum(ref in gt for gt, ref, _ in site_alleles.itertuples(index=False))
    )

    rs_match = m[m["RS# (dbSNP)"] != "-1"]
    rs_agree = float((rs_match["rsid"] == "rs" + rs_match["RS# (dbSNP)"]).mean())

    m["alt_copies"] = [gt.count(a) for gt, a in zip(m["gt"], m["alt"])]
    # "c.1601=" style records describe the reference-restoring allele where GRCh37 itself carries the
    # variant (e.g. F5 rs6025); carrying their "alt" means carrying the normal allele, so they are not hits.
    m["ref_restoring"] = m["Name"].str.contains(r"c\.\S*=", na=False) & ~m[
        "Name"
    ].str.contains(">", na=False)
    # Sites whose called alleles don't fit ClinVar's alleles on the forward strand are unreliable matches.
    # "c.=" records are kept out on purpose: their classification describes the transcript's allele, not the
    # GRCh37 base (checked 2026-09-28: 497 of 566 such positions have no record for the GRCh37 allele at all).
    fits = pd.Series([(c, p) not in bad_sites for c, p in zip(m["chrom"], m["pos"])], index=m.index)
    hits = m[(m["alt_copies"] > 0) & ~m["ref_restoring"] & fits].copy()
    hits["category"] = hits["ClinicalSignificance"].map(category)
    hits["stars"] = hits["ReviewStatus"].map(STARS).fillna(0).astype(int)
    fam_sex = (
        (ROOT / "data" / "genome.fam").read_text().split()[4]
    )  # plink-inferred: 1 male, 2 female
    haploid_chroms = ["X", "Y"] if fam_sex == "1" else ["Y"]
    male_haploid = hits["chrom"].isin(haploid_chroms)
    hits["zygosity"] = "heterozygous"
    hits.loc[hits["alt_copies"] == 2, "zygosity"] = "homozygous"
    hits.loc[male_haploid, "zygosity"] = "hemizygous"
    hits.loc[male_haploid & (hits["gt"].str[0] != hits["gt"].str[1]), "zygosity"] = (
        "het call on haploid (likely error)"
    )

    out = hits[
        [
            "category",
            "stars",
            "zygosity",
            "GeneSymbol",
            "chrom",
            "pos",
            "rsid",
            "gt",
            "ref",
            "alt",
            "ClinicalSignificance",
            "ReviewStatus",
            "NumberSubmitters",
            "PhenotypeList",
            "Name",
            "VariationID",
            "#AlleleID",
            "OriginSimple",
        ]
    ]
    out = out.sort_values(
        ["category", "stars", "GeneSymbol"], ascending=[True, False, True]
    )
    out.to_csv(ROOT / "results" / "clinvar_hits.tsv", sep="\t", index=False)

    # Well-known markers, reported whatever the genotype (ref calls included).
    sel = m[m["rsid"].isin(selected())].copy()
    sel["category"] = sel["ClinicalSignificance"].map(category)
    sel["stars"] = sel["ReviewStatus"].map(STARS).fillna(0).astype(int)
    sel[
        [
            "rsid",
            "GeneSymbol",
            "chrom",
            "pos",
            "gt",
            "ref",
            "alt",
            "alt_copies",
            "category",
            "stars",
            "ClinicalSignificance",
            "ReviewStatus",
            "Name",
            "PhenotypeList",
            "VariationID",
        ]
    ].to_csv(ROOT / "results" / "selected_markers.tsv", sep="\t", index=False)

    m["category"] = m["ClinicalSignificance"].map(category)
    m["stars"] = m["ReviewStatus"].map(STARS).fillna(0).astype(int)
    plp = m[m["category"] == "Pathogenic / likely pathogenic"]
    plp_tested = {
        s: int(plp.loc[plp["stars"] >= s, ["chrom", "pos"]].drop_duplicates().shape[0])
        for s in (0, 1, 2, 3)
    }

    summary = {
        "clinvar_file": CLINVAR.name,
        "clinvar_grch37_snv_rows": len(cv),
        "genome_called_snvs": len(g),
        "matched_rows": len(m),
        "matched_sites": len(site_alleles),
        "strand_forward_consistent": fwd,
        "strand_only_if_complemented": rc,
        "strand_inconsistent": neither,
        "sites_carrying_ref_allele": ref_present,
        "rsid_agreement": rs_agree,
        "plp_sites_tested_by_min_stars": plp_tested,
        "hit_rows": len(hits),
        "hit_sites": int(hits[["chrom", "pos"]].drop_duplicates().shape[0]),
        "sites_excluded_strand": len(bad_sites),
        "hits_by_category": hits["category"].value_counts().to_dict(),
        "hits_by_category_stars": {
            f"{c}|{s}": int(n)
            for (c, s), n in hits.groupby(["category", "stars"]).size().items()
        },
        "hits_by_category_zygosity": {
            f"{c}|{z}": int(n)
            for (c, z), n in hits.groupby(["category", "zygosity"]).size().items()
        },
    }
    (ROOT / "results" / "clinvar_summary.json").write_text(
        json.dumps(summary, indent=1)
    )
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
