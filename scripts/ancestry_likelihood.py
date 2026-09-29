"""PCA-free check of population similarity: genotype log-likelihood of SELF under each population's
allele frequencies (AADR present-day HO, LD-pruned shared SNPs).

log L(pop) = sum over SNPs of log Binomial(g_self; 2, f_pop), f_pop = (a+0.5)/(n+1). Every class is
subsampled to the same size (SUB_N people, N_REP replicates) because frequency noise from small groups
otherwise dominates the ranking. Reports mean ΔlogL to the best class and how often each class ranks first.
Writes results/membership/aadr_likelihood.json.
"""

import json
import subprocess
from pathlib import Path

import sys

import numpy as np
import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
P1 = str(ROOT / "env" / "bin" / "plink")
D = ROOT / "data" / "aadr"
OUT = ROOT / "results" / "membership"
SUB_N = 10
N_REP = 100


def main() -> None:
    sys.path.insert(0, str(ROOT / "scripts"))
    OUT.mkdir(parents=True, exist_ok=True)
    raw = OUT / "aadr_mp"
    if not Path(f"{raw}.raw").exists():
        subprocess.run(
            [
                P1,
                "--bfile",
                str(D / "mp"),
                "--recode",
                "A",
                "--out",
                str(raw),
                "--silent",
            ],
            check=True,
        )
    df = pd.read_csv(f"{raw}.raw", sep=r"\s+")
    iid = df["IID"].to_numpy()
    G = df.iloc[:, 6:].to_numpy(dtype=float)  # count of the A1 allele, NaN = missing
    pops = pd.read_csv(D / "id_pop.tsv", sep=r"\s+", names=["IID", "pop"]).set_index(
        "IID"
    )["pop"]
    lab = np.array(
        [pops.get(i, "SELF" if i == "SELF" else None) for i in iid], dtype=object
    )
    me = G[iid == "SELF"][0]
    ok = ~np.isnan(me)
    # Pool groups that can't be told apart (see ancestry_membership.py), then compare every class at the same
    # sample size: frequencies from small groups are noisier and cost likelihood, which otherwise makes the
    # ranking follow n (first attempt: Spanish n=173 > Russian n=71 > French n=60 > ...).
    from ancestry_membership import pools

    pooled = {p: g for g, ps in pools().items() for p in ps}
    lab = np.array(
        [pooled.get(x, x) if x is not None else None for x in lab], dtype=object
    )
    counts = pd.Series(lab[(lab != "SELF") & (lab != None)]).value_counts()  # noqa: E711
    classes = sorted(counts[counts >= SUB_N].index)
    x = me[ok]
    het = np.log(np.where(x == 1, 2, 1))
    rng = np.random.default_rng(20260928)
    reps = np.zeros((N_REP, len(classes)))
    for r in range(N_REP):
        for k, c in enumerate(classes):
            idx = rng.choice(np.flatnonzero(lab == c), SUB_N, replace=False)
            g = G[idx][:, ok]
            f = (np.nansum(g, axis=0) + 0.5) / (2 * np.sum(~np.isnan(g), axis=0) + 1)
            reps[r, k] = np.sum(x * np.log(f) + (2 - x) * np.log(1 - f) + het)
    delta = reps - reps.max(axis=1, keepdims=True)
    mean_d = delta.mean(axis=0)
    rank1 = (delta == 0).mean(axis=0)
    order = np.argsort(-mean_d)
    n_snp = int(ok.sum())
    best = classes[order[0]]
    res = {
        "n_snps": n_snp,
        "subsample_n": SUB_N,
        "replicates": N_REP,
        "best": best,
        "ranking": [
            {
                "pop": classes[k],
                "mean_delta_logL": float(mean_d[k]),
                "n_total": int(counts[classes[k]]),
                "rank1_fraction": float(rank1[k]),
            }
            for k in order[:10]
        ],
    }
    (OUT / "aadr_likelihood.json").write_text(json.dumps(res, indent=1))
    print(
        f"{n_snp} SNPs, {SUB_N} per class x {N_REP};",
        "; ".join(
            f"{r['pop']} Δ{r['mean_delta_logL']:.1f} (n={r['n_total']}, rank1 {r['rank1_fraction']:.2f})"
            for r in res["ranking"][:8]
        ),
    )


if __name__ == "__main__":
    main()
