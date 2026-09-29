"""Imputation accuracy against truth: the 503 1000G Europeans imputed from chip positions only.

scripts/ref_impute.sh cut 503 sequenced Europeans down to SELF's chip sites and imputed them with a panel that
excludes them, exactly as SELF was imputed. At the imputed (non-chip) sites their sequenced genotypes are the truth.
For a random 1-in-STEP sample of imputed biallelic SNVs per chromosome this computes, by minor-allele frequency bin:
  - dosage r^2 (squared correlation of imputed dosage with true allele count, across the 503 people);
  - non-reference concordance (best-guess genotype correct, among genotypes where truth or call is non-reference);
  - accuracy by max genotype probability (calibration of the GP flags used elsewhere in the report).
Replaces the chr22-only masked check, which was dominated by common homozygous sites.
Writes results/impute_benchmark.json.
"""

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
BCF = str(ROOT / "env" / "bin" / "bcftools")
P2 = str(ROOT / "env-x86" / "bin" / "plink2")
D = ROOT / "data" / "refimp"
STEP = 200
MAF_BINS = [0, 0.001, 0.005, 0.01, 0.05, 0.2, 0.5]
GP_BINS = [0, 0.8, 0.9, 0.99, 1.0001]


def chrom_data(c: int, rng: np.random.Generator):
    vcf = D / f"chr{c}.vcf.gz"
    samples = subprocess.run(
        [BCF, "query", "-l", str(vcf)], capture_output=True, text=True, check=True
    ).stdout.split()
    eur = [s for s in samples if s != "SELF"]
    # 1) positions of imputed biallelic SNVs, without genotype data (fast); 2) sample 1 in STEP;
    # 3) one indexed query for dosage and genotype probability at just those sites.
    sites = subprocess.run([BCF, "query", "-i", "INFO/IMP=1 && TYPE=\"snp\" && N_ALT=1", "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
                            "-s", "SELF", str(vcf)], capture_output=True, text=True, check=True).stdout.splitlines()
    pick = [sites[k] for k in np.flatnonzero(rng.random(len(sites)) < 1 / STEP)]
    tgt = D / f"bench_targets.chr{c}.tsv"
    tgt.write_text("".join("\t".join(s.split("\t")[:2]) + "\n" for s in pick))
    q = subprocess.run([BCF, "query", "-T", str(tgt), "-i", "INFO/IMP=1 && TYPE=\"snp\" && N_ALT=1", "-s", ",".join(eur),
                        "-f", "%CHROM:%POS:%REF:%ALT[\t%DS|%GP]\n", str(vcf)], capture_output=True, text=True, check=True).stdout
    rows = [l.split("\t") for l in q.splitlines()]
    ids = [r[0] for r in rows]
    ds = np.array([[float(x.split("|")[0]) for x in r[1:]] for r in rows])
    maxgp = np.array([[max(map(float, x.split("|")[1].split(","))) for x in r[1:]] for r in rows])
    # Truth from the sequenced 1000G genotypes, same people, same variant IDs.
    idf = D / f"bench_ids.chr{c}.txt"
    idf.write_text("\n".join(ids) + "\n")
    out = D / f"bench_truth.chr{c}"
    Path(str(idf).replace(".txt", ".alt")).write_text(
        "\n".join(f"{i}\t{i.split(':')[3]}" for i in ids) + "\n"
    )
    subprocess.run(
        [
            P2,
            "--pfile",
            str(ROOT / "reference" / "panel" / "all_phase3"),
            "vzs",
            "--chr",
            str(c),
            "--keep",
            str(D / "eur503.txt"),
            "--snps-only",
            "just-acgt",
            "--max-alleles",
            "2",
            "--set-all-var-ids",
            "@:#:$r:$a",
            "--rm-dup",
            "exclude-all",
            "--extract",
            str(idf),
            "--export",
            "A",
            "--export-allele",
            str(idf).replace(".txt", ".alt"),
            "--out",
            str(out),
            "--silent",
        ],
        check=True,
    )
    raw = pd.read_csv(f"{out}.raw", sep=r"\s+")
    raw = raw.set_index("IID").loc[eur]
    tcols = {c2.rsplit("_", 1)[0]: c2 for c2 in raw.columns[5:]}
    sel = [k for k, i in enumerate(ids) if i in tcols]
    truth = raw[[tcols[ids[k]] for k in sel]].to_numpy(dtype=float).T
    return ds[sel], maxgp[sel], truth


def main() -> None:
    rng = np.random.default_rng(20260928)
    DS, GP, T = [], [], []
    for c in range(1, 23):
        if not (D / f"chr{c}.vcf.gz").exists():
            continue
        ds, gp, t = chrom_data(c, rng)
        DS.append(ds), GP.append(gp), T.append(t)
        print(f"chr{c}: {len(ds)} sampled imputed sites")
    ds, gp, t = np.vstack(DS), np.vstack(GP), np.vstack(T)
    ok = ~np.isnan(t).any(axis=1)
    ds, gp, t = ds[ok], gp[ok], t[ok]
    af = t.mean(axis=1) / 2
    maf = np.minimum(af, 1 - af)
    best = np.rint(ds)
    res = {
        "n_sites": int(len(ds)),
        "n_people": int(ds.shape[1]),
        "sampling": f"1 in {STEP} imputed SNVs",
        "by_maf": [],
        "by_gp": [],
    }
    for lo, hi in zip(MAF_BINS[:-1], MAF_BINS[1:]):
        m = (maf > lo) & (maf <= hi)
        if m.sum() < 20:
            continue
        r2 = [
            np.corrcoef(a, b)[0, 1] ** 2
            for a, b in zip(ds[m], t[m])
            if a.std() > 0 and b.std() > 0
        ]
        nonref = (t[m] > 0) | (best[m] > 0)
        conc = float((best[m][nonref] == t[m][nonref]).mean()) if nonref.any() else None
        res["by_maf"].append(
            {
                "maf": f"{lo}-{hi}",
                "sites": int(m.sum()),
                "median_dosage_r2": float(np.median(r2)) if r2 else None,
                "nonref_concordance": conc,
            }
        )
    for lo, hi in zip(GP_BINS[:-1], GP_BINS[1:]):
        m = (gp >= lo) & (gp < hi)
        if m.sum():
            res["by_gp"].append(
                {
                    "max_gp": f"{lo}-{min(hi, 1)}",
                    "genotypes": int(m.sum()),
                    "accuracy": float((best[m] == t[m]).mean()),
                }
            )
    (ROOT / "results" / "impute_benchmark.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
