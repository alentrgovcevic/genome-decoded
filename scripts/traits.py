"""Single-SNP traits and blood group from the raw array, plus two imputed sites (ACTN3, ABO O deletion).

Effect alleles are on the GRCh37 plus strand; sources are in SOURCES (local ClinVar rows plus dbSNP/
SNPedia pages checked on 2026-09-28). Writes results/traits.json.
"""

import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
BCF = str(ROOT / "env" / "bin" / "bcftools")

# rsid, gene, trait, effect allele (plus strand), what the effect allele does, source
TRAITS = [
    (
        "rs4988235",
        "LCT/MCM6",
        "Lactose tolerance",
        "A",
        "lactase persistence (can digest milk as an adult)",
        "ClinVar 'LACTASE PERSISTENCE' (G>A); dbSNP rs4988235",
    ),
    (
        "rs12913832",
        "HERC2",
        "Eye colour",
        "G",
        "blue eyes (SNPedia: GG blue 99% of the time, AA brown 80%)",
        "ClinVar 'BLUE/NONBLUE EYES' (A>G); SNPedia Rs12913832",
    ),
    (
        "rs17822931",
        "ABCC11",
        "Earwax type",
        "T",
        "dry earwax (TT); C is wet",
        "ClinVar c.538G>A (C>T plus strand); SNPedia Rs17822931",
    ),
    (
        "rs762551",
        "CYP1A2",
        "Caffeine metabolism",
        "A",
        "CYP1A2*1F, faster (inducible) caffeine metabolism (SNPedia: AA)",
        "ClinVar c.-9-154C>A; SNPedia Rs762551",
    ),
    (
        "rs671",
        "ALDH2",
        "Alcohol flush",
        "A",
        "ALDH2*2: alcohol flushing / hangover susceptibility",
        "ClinVar c.1510G>A (drug response)",
    ),
    (
        "rs28929474",
        "SERPINA1",
        "Alpha-1 antitrypsin (PiZ)",
        "T",
        "Z allele (ClinVar: pathogenic for alpha-1 antitrypsin deficiency)",
        "ClinVar c.1096G>A p.Glu366Lys, Pathogenic",
    ),
    (
        "rs17580",
        "SERPINA1",
        "Alpha-1 antitrypsin (PiS)",
        "A",
        "S allele (ClinVar: pathogenic, low penetrance, for alpha-1 antitrypsin deficiency)",
        "ClinVar c.863A>T p.Glu288Val, Pathogenic/low penetrance",
    ),
]
TAS2R38 = [
    ("rs713598", "G", "C"),
    ("rs1726866", "G", "A"),
    ("rs10246939", "C", "T"),
]  # (rsid, PAV allele, AVI allele)
ABO_B = ("rs8176746", "T")  # B-tagging p.Leu266Met, plus strand T = B (dbSNP)
ABO_A2 = (
    "rs1053878",
    "A",
)  # A2-tagging p.Pro156Leu, plus strand A = A2 (dbSNP, PMC7906510)


def raw_genotypes() -> dict:
    g = pd.read_csv(
        ROOT / "data" / "genome_23andme.txt",
        sep="\t",
        comment="#",
        names=["rsid", "chrom", "pos", "gt"],
        dtype=str,
    )
    return dict(zip(g["rsid"], g["gt"]))


def imputed(chrom: str, pos: int) -> list[dict]:
    fp = ROOT / "data" / "imputed" / f"chr{chrom}.vcf.gz"
    if not fp.exists():
        return []
    out = subprocess.run(
        [
            BCF,
            "query",
            "-r",
            f"{chrom}:{pos}",
            "-f",
            "%ID\t%POS\t%REF\t%ALT\t%INFO/DR2\t[%GT\t%DS\t%GP]\n",
            str(fp),
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    rows = []
    for line in out.strip().splitlines():
        rid, p, ref, alt, dr2, gt, ds, gp = line.split("\t")
        rows.append(
            {
                "id": rid,
                "pos": int(p),
                "ref": ref,
                "alt": alt,
                "dr2": float(dr2),
                "gt": gt,
                "ds": float(ds),
                "gp": gp,
            }
        )
    return rows


def main() -> None:
    g = raw_genotypes()
    out = {"single": [], "tas2r38": None, "abo": None, "actn3": None}
    for rsid, gene, trait, eff, meaning, src in TRAITS:
        gt = g.get(rsid, "not typed")
        out["single"].append(
            {
                "rsid": rsid,
                "gene": gene,
                "trait": trait,
                "genotype": gt,
                "effect_allele": eff,
                "copies": gt.count(eff) if gt != "not typed" else None,
                "meaning": meaning,
                "source": src,
            }
        )

    gts = {r: g.get(r, "--") for r, _, _ in TAS2R38}
    n_het = sum(v[0] != v[1] for v in gts.values())
    pav = [sum(c == p for c in gts[r]) for r, p, _ in TAS2R38]
    out["tas2r38"] = {
        "genotypes": gts,
        "pav_copies_per_site": pav,
        "n_het": n_het,
        "call": "PAV/AVI (one taster, one non-taster haplotype)"
        if pav == [1, 1, 1]
        else (
            "PAV/PAV"
            if pav == [2, 2, 2]
            else ("AVI/AVI" if pav == [0, 0, 0] else "rare combination")
        ),
    }

    # ACTN3 R577X (rs1815739, chr11:66328095) — not on the chip, imputed. T = X.
    a = [r for r in imputed("11", 66328095)]
    out["actn3"] = a[0] if a else None

    # ABO: B allele from typed rs8176746; O allele (rs8176719 frameshift deletion) from imputation.
    # 1000G codes rs8176719 as REF (no insertion, = O) / ALT (+C, functional A or B), per dbSNP.
    o = [r for r in imputed("9", 136132908) if len(r["alt"]) > len(r["ref"])]
    b_copies = g.get(ABO_B[0], "--").count(ABO_B[1])
    a2 = g.get(ABO_A2[0], "--").count(ABO_A2[1])
    abo = {
        "rs8176746": g.get(ABO_B[0]),
        "b_copies": b_copies,
        "rs1053878": g.get(ABO_A2[0]),
        "a2_copies": a2,
    }
    if o:
        o_copies = round(2 - o[0]["ds"])
        a_copies = 2 - o_copies - b_copies
        alleles = ["A"] * a_copies + ["B"] * b_copies + ["O"] * o_copies
        group = (
            "AB"
            if "A" in alleles and "B" in alleles
            else ("A" if "A" in alleles else ("B" if "B" in alleles else "O"))
        )
        abo.update(
            {
                "rs8176719": o[0],
                "o_copies": o_copies,
                "a_copies": a_copies,
                "alleles": "/".join(alleles),
                "group": group,
            }
        )
    out["abo"] = abo
    (ROOT / "results" / "traits.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
