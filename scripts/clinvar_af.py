"""1000 Genomes allele frequencies for every ClinVar record where SELF carries the variant allele.

Local substitute for gnomAD: querying gnomAD online would send this person's variant positions to an
external service, and the full gnomAD site files are hundreds of GB. 1000G phase 3 (2,504 people) resolves
frequencies down to ~0.02% and is enough to separate rare (< 1%) from common variants.
Reads results/clinvar_hits.tsv (measured) and results/clinvar_imputed.tsv (imputed).
Writes results/clinvar_af.tsv: chrom, pos, ref, alt, af_all, af_eur, source.
"""

import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
P2 = str(ROOT / "env-x86" / "bin" / "plink2")
BCF = str(ROOT / "env" / "bin" / "bcftools")
D = ROOT / "data" / "clinvar_af"


def imputed_alleles(imp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c, sub in imp.groupby("chrom"):
        t = D / f"imp_targets.{c}.tsv"
        sub[["chrom", "pos"]].drop_duplicates().sort_values("pos").to_csv(
            t, sep="\t", header=False, index=False
        )
        q = subprocess.run(
            [
                BCF,
                "query",
                "-T",
                str(t),
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\n",
                str(ROOT / "data" / "imputed" / f"chr{c}.vcf.gz"),
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        rows += [l.split("\t") for l in q.strip().splitlines()]
    return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt"]).astype(
        {"pos": int}
    )


def freqs(sites: pd.DataFrame, keep: Path | None, out: str) -> pd.Series:
    rng = D / "sites.range"
    sites[["chrom", "pos"]].assign(end=sites["pos"]).drop_duplicates().to_csv(
        rng, sep="\t", header=False, index=False
    )
    args = [
        P2,
        "--pfile",
        str(ROOT / "reference" / "panel" / "all_phase3"),
        "vzs",
        "--extract",
        "range",
        str(rng),
        "--set-all-var-ids",
        "@:#:$r:$a",
        "--new-id-max-allele-len",
        "1000",
        "missing",
        "--nonfounders",
        "--freq",
        "--out",
        str(D / out),
        "--silent",
    ]
    if keep:
        args += ["--keep", str(keep)]
    subprocess.run(args, check=True)
    f = pd.read_csv(D / f"{out}.afreq", sep="\t")
    f = f[~f["ALT_FREQS"].astype(str).str.contains(",")]
    return pd.Series(f["ALT_FREQS"].astype(float).to_numpy(), index=f["ID"])


def main() -> None:
    D.mkdir(parents=True, exist_ok=True)
    meas = pd.read_csv(
        ROOT / "results" / "clinvar_hits.tsv", sep="\t", dtype={"chrom": str}
    )
    meas = meas[["chrom", "pos", "ref", "alt"]].assign(source="measured")
    imp = pd.read_csv(
        ROOT / "results" / "clinvar_imputed.tsv", sep="\t", dtype={"chrom": str}
    )
    imp = imputed_alleles(imp).assign(source="imputed")
    sites = pd.concat([meas, imp]).drop_duplicates(["chrom", "pos", "ref", "alt"])
    sites = sites[sites["chrom"].isin([str(i) for i in range(1, 23)] + ["X"])]
    psam = pd.read_csv(ROOT / "reference" / "panel" / "all_phase3.psam", sep="\t")
    eur = D / "eur.keep"
    psam.loc[psam["SuperPop"] == "EUR", ["#IID"]].to_csv(eur, index=False, header=False)
    f_all = freqs(sites, None, "all")
    f_eur = freqs(sites, eur, "eur")
    ids = (
        sites["chrom"]
        + ":"
        + sites["pos"].astype(str)
        + ":"
        + sites["ref"]
        + ":"
        + sites["alt"]
    )
    sites["af_all"] = ids.map(f_all)
    sites["af_eur"] = ids.map(f_eur)
    sites.to_csv(ROOT / "results" / "clinvar_af.tsv", sep="\t", index=False)
    known = sites["af_all"].notna()
    print(
        f"{len(sites)} sites; with 1000G frequency: {known.sum()}; rare (<1% overall): {(sites['af_all'] < 0.01).sum()}; "
        f"absent from 1000G: {(~known).sum()}"
    )


if __name__ == "__main__":
    main()
