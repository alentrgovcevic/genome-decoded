"""Genome-wide chip-artefact check: homozygous calls for rare alleles, observed vs expected.

For every measured autosomal SNV that 1000G also has, take the 1000G European frequency q of the allele
SELF is homozygous for. Under Hardy-Weinberg the expected number of homozygotes for alleles with q < 1% is
sum(q^2) over those sites; the observed count above that is a direct estimate of how many such calls are
probe artefacts (Weedon et al. 2021 BMJ: rare variant calls on consumer arrays are mostly false). Sites are
written to results/qc_rare_hom.tsv so other sections can flag them. Summary: results/qc_rare_hom.json.
"""

import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
P2 = str(ROOT / "env-x86" / "bin" / "plink2")
D = ROOT / "data" / "qc_rare_hom"
RARE = 0.01


def main() -> None:
    D.mkdir(parents=True, exist_ok=True)
    g = pd.read_csv(
        ROOT / "data" / "genome_23andme.txt",
        sep="\t",
        comment="#",
        names=["rsid", "chrom", "pos", "gt"],
        dtype={"chrom": str},
    )
    g = g[
        g["chrom"].isin([str(i) for i in range(1, 23)])
        & g["gt"].str.fullmatch("[ACGT]{2}")
    ]
    rng = D / "sites.range"
    g[["chrom", "pos"]].assign(end=g["pos"]).to_csv(
        rng, sep="\t", header=False, index=False
    )
    subprocess.run(
        [
            P2,
            "--pfile",
            str(ROOT / "reference" / "panel" / "all_phase3"),
            "vzs",
            "--extract",
            "range",
            str(rng),
            "--keep",
            str(ROOT / "reference" / "panel" / "eur_unrelated.keep"),
            "--snps-only",
            "just-acgt",
            "--max-alleles",
            "2",
            "--set-all-var-ids",
            "@:#:$r:$a",
            "--rm-dup",
            "exclude-all",
            "--nonfounders",
            "--freq",
            "--out",
            str(D / "eur"),
            "--silent",
        ],
        check=True,
    )
    f = pd.read_csv(D / "eur.afreq", sep="\t")
    f[["chrom", "pos", "ref", "alt"]] = f["ID"].str.split(":", expand=True)
    f["pos"] = f["pos"].astype(int)
    m = g.merge(f[["chrom", "pos", "ref", "alt", "ALT_FREQS"]], on=["chrom", "pos"])
    hom = m[m["gt"].str[0] == m["gt"].str[1]].copy()
    hom["allele"] = hom["gt"].str[0]
    hom = hom[
        hom["allele"].isin(hom["ref"])
        | (hom["allele"] == hom["alt"])
        | (hom["allele"] == hom["ref"])
    ]
    hom["q"] = [
        a if al == alt else 1 - a
        for al, alt, a in zip(hom["allele"], hom["alt"], hom["ALT_FREQS"])
    ]
    hom = hom[(hom["allele"] == hom["ref"]) | (hom["allele"] == hom["alt"])]
    all_sites = m.copy()
    all_sites["q_ref"] = 1 - all_sites["ALT_FREQS"]
    # Expected homozygotes for a rare allele, over every compared site: each allele with frequency q < 1% contributes q^2.
    exp_alt = float(
        (all_sites.loc[all_sites["ALT_FREQS"] < RARE, "ALT_FREQS"] ** 2).sum()
    )
    exp_ref = float((all_sites.loc[all_sites["q_ref"] < RARE, "q_ref"] ** 2).sum())
    rare = hom[hom["q"] < RARE].copy()
    rare[["rsid", "chrom", "pos", "gt", "q"]].to_csv(
        ROOT / "results" / "qc_rare_hom.tsv", sep="\t", index=False
    )
    res = {
        "sites_compared": int(len(m)),
        "homozygous_calls": int(len(hom)),
        "rare_allele_threshold": RARE,
        "observed_rare_homozygotes": int(len(rare)),
        "expected_rare_homozygotes": round(exp_alt + exp_ref, 2),
        "observed_rare_hom_q0": int((rare["q"] == 0).sum()),
    }
    res["excess_likely_artefacts"] = round(
        res["observed_rare_homozygotes"] - res["expected_rare_homozygotes"], 1
    )
    (ROOT / "results" / "qc_rare_hom.json").write_text(json.dumps(res, indent=1))
    print(res)


if __name__ == "__main__":
    main()
