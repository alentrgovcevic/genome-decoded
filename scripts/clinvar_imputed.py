"""ClinVar records at imputed (not measured) positions where SELF carries the variant allele.

Scans data/imputed/chrN.vcf.gz (Beagle IMP flag set) against ClinVar GRCh37 SNVs whose category is
pathogenic/likely pathogenic, risk factor, drug response or protective. Reports the max genotype
probability, so confident calls can be separated from guesses. Writes results/clinvar_imputed.tsv.
"""

import subprocess
from pathlib import Path

import pandas as pd

from clinvar import STARS, category, load_clinvar

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
BCF = str(ROOT / "env" / "bin" / "bcftools")
KEEP = {"Pathogenic / likely pathogenic", "Risk factor", "Drug response", "Protective"}


def main() -> None:
    cv = load_clinvar()
    cv = cv[
        ~(
            cv["Name"].str.contains(r"c\.\S*=", na=False)
            & ~cv["Name"].str.contains(">", na=False)
        )
    ]
    cv["category"] = cv["ClinicalSignificance"].map(category)
    cv = cv[cv["category"].isin(KEEP)]
    cv["stars"] = cv["ReviewStatus"].map(STARS).fillna(0).astype(int)
    rows = []
    for c in [str(i) for i in range(1, 23)] + ["X"]:
        fp = ROOT / "data" / "imputed" / f"chr{c}.vcf.gz"
        sub = cv[cv["chrom"] == c]
        if not fp.exists() or sub.empty:
            continue
        targets = ROOT / "data" / "imputed" / f"clinvar_targets.chr{c}.tsv"
        sub[["chrom", "pos"]].drop_duplicates().sort_values("pos").to_csv(
            targets, sep="\t", header=False, index=False
        )
        q = subprocess.run(
            [
                BCF,
                "query",
                "-T",
                str(targets),
                "-i",
                "INFO/IMP=1",
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DS\t%GP]\n",
                str(fp),
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        for line in q.strip().splitlines():
            ch, pos, ref, alt, gt, ds, gp = line.split("\t")
            rows.append(
                {
                    "chrom": ch,
                    "pos": int(pos),
                    "ref_imp": ref,
                    "alt_imp": alt,
                    "gt": gt,
                    "ds": float(ds),
                    # haploid X records have no GP; the dosage is then the probability of the ALT allele
                    "max_gp": max(map(float, gp.split(","))) if gp != "." else max(float(ds), 1 - float(ds)),
                }
            )
    imp = pd.DataFrame(rows)
    m = imp.merge(
        cv,
        left_on=["chrom", "pos", "ref_imp", "alt_imp"],
        right_on=["chrom", "pos", "ref", "alt"],
    )
    m["alt_copies"] = (
        m["gt"].str.replace("|", "/").str.split("/").map(lambda a: a.count("1"))
    )
    hits = m[m["alt_copies"] > 0].sort_values(
        ["category", "max_gp"], ascending=[True, False]
    )
    cols = [
        "category",
        "stars",
        "GeneSymbol",
        "chrom",
        "pos",
        "RS# (dbSNP)",
        "gt",
        "alt_copies",
        "max_gp",
        "ds",
        "ClinicalSignificance",
        "PhenotypeList",
        "Name",
        "VariationID",
    ]
    hits[cols].to_csv(ROOT / "results" / "clinvar_imputed.tsv", sep="\t", index=False)
    print(
        f"imputed ClinVar-relevant sites checked: {len(m)}; carrying variant allele: {len(hits)}"
    )
    print(
        hits.groupby(
            ["category", pd.cut(hits["max_gp"], [0, 0.9, 0.99, 1.0001], right=False)],
            observed=True,
        ).size()
    )


if __name__ == "__main__":
    main()
