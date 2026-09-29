"""Population membership without projection bias, with a check of how far labels can be trusted.

For each reference set (AADR West Eurasia, the AADR regional zoom subset, 1000G EUR):
  1. PCA with SELF *in* the sample (plink2 --pca), so SELF is not shrunk toward the centre as a projected
     sample would be.
  2. Keep only PCs whose eigenvalue exceeds the Marchenko-Pastur noise edge (1 + sqrt(n/p))^2 × 1.05.
  3. Linear discriminant analysis on those PCs (pooled within-population covariance, equal priors):
     posterior membership probability for SELF per population and per region, plus a chi-square test
     of whether SELF fits its best population at all (Mahalanobis distance, k df).
  4. Leave-one-out: each reference individual is classified using everyone else; the per-population
     accuracy says how distinguishable that label is.
Writes results/membership/<label>.json.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
P2 = str(ROOT / "env-x86" / "bin" / "plink2")
OUT = ROOT / "results" / "membership"
MIN_N = 5
def pools() -> dict[str, list[str]]:
    """Groups the leave-one-out check can't tell apart (config: ancestry.pooled_groups)."""
    from run_config import get
    return get("ancestry.pooled_groups")


def pca(
    bfile: Path, keep: Path, out: Path, n_pcs: int = 10
) -> tuple[pd.DataFrame, np.ndarray, int]:
    subprocess.run(
        [
            P2,
            "--bfile",
            str(bfile),
            "--keep",
            str(keep),
            "--maf",
            "0.01",
            "--pca",
            str(n_pcs),
            "--out",
            str(out),
            "--silent",
        ],
        check=True,
    )
    ev = pd.read_csv(f"{out}.eigenvec", sep="\t").rename(
        columns={"#FID": "FID", "#IID": "IID"}
    )
    val = np.loadtxt(f"{out}.eigenval")
    # SNPs actually used = "N variants remaining after main filters" in the plink2 log (after --maf).
    line = next(l for l in open(f"{out}.log") if "remaining after main filters" in l)
    n_snp = int(line.split()[0])
    return ev, val, n_snp


def lda(X: np.ndarray, y: np.ndarray, classes: list[str]):
    means = np.stack([X[y == c].mean(axis=0) for c in classes])
    resid = np.concatenate([X[y == c] - means[i] for i, c in enumerate(classes)])
    cov = resid.T @ resid / (len(X) - len(classes))
    return means, np.linalg.inv(cov)


def posterior(
    x: np.ndarray, means: np.ndarray, icov: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    d = x - means
    m2 = np.einsum("ij,jk,ik->i", d, icov, d)
    logp = -0.5 * m2
    p = np.exp(logp - logp.max())
    return p / p.sum(), m2


def run(label: str, bfile: Path, keep: Path, pop_of: dict, region_of: dict) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    ev, val, n_snp = pca(bfile, keep, OUT / f"pca_{label}")
    ev["pop"] = ev["IID"].map(pop_of)
    ev.loc[ev["IID"] == "SELF", "pop"] = "SELF"
    ref = ev[ev["IID"] != "SELF"].dropna(subset=["pop"])
    counts = ref["pop"].value_counts()
    classes = sorted(counts[counts >= MIN_N].index)
    ref = ref[ref["pop"].isin(classes)]
    n = len(ref) + 1
    edge = (1 + np.sqrt(n / n_snp)) ** 2
    k = int(max(1, (val > edge * 1.05).sum()))
    pcs = [f"PC{i}" for i in range(1, k + 1)]
    X, y = ref[pcs].to_numpy(), ref["pop"].to_numpy()
    means, icov = lda(X, y, classes)
    me = ev.loc[ev["IID"] == "SELF", pcs].to_numpy()[0]
    post, m2 = posterior(me, means, icov)
    best = int(np.argmax(post))
    fit_p = float(stats.chi2.sf(m2[best], df=k))
    # Leave-one-out accuracy per population.
    correct = {c: 0 for c in classes}
    for i in range(len(X)):
        mask = np.ones(len(X), bool)
        mask[i] = False
        mu, ic = lda(X[mask], y[mask], classes)
        p, _ = posterior(X[i], mu, ic)
        correct[y[i]] += classes[int(np.argmax(p))] == y[i]
    loo = {c: correct[c] / int((y == c).sum()) for c in classes}
    region_post = {}
    for c, p in zip(classes, post):
        r = region_of.get(c, c)
        region_post[r] = region_post.get(r, 0) + float(p)
    order = np.argsort(-post)
    res = {
        "label": label,
        "n_ref": int(len(ref)),
        "n_snps": n_snp,
        "eigenvalues": val.tolist(),
        "mp_edge": float(edge),
        "n_pcs_used": k,
        "self_top": [
            {
                "pop": classes[i],
                "posterior": float(post[i]),
                "n": int(counts[classes[i]]),
                "loo_accuracy": loo[classes[i]],
            }
            for i in order[:8]
        ],
        "self_region_posterior": dict(
            sorted(region_post.items(), key=lambda kv: -kv[1])
        ),
        "fit_p_best": fit_p,
        "loo_accuracy": loo,
        "loo_overall": float(sum(correct.values()) / len(X)),
    }
    (OUT / f"{label}.json").write_text(json.dumps(res, indent=1))
    top = res["self_top"][:4]
    print(
        label,
        f"k={k} PCs (edge {edge:.3f}); LOO overall {res['loo_overall']:.2f}; fit p={fit_p:.3f};",
        ", ".join(
            f"{t['pop']} {t['posterior']:.2f} (LOO {t['loo_accuracy']:.2f}, n={t['n']})"
            for t in top
        ),
    )
    print(
        "   regions:",
        {k2: round(v, 2) for k2, v in res["self_region_posterior"].items()},
    )
    return res


def main() -> None:
    sys.path.insert(0, str(ROOT / "scripts"))
    from run_config import get

    region_of = {p: r for r, ps in get("ancestry.regions").items() for p in ps}
    idp = pd.read_csv(
        ROOT / "data" / "aadr" / "id_pop.tsv", sep=r"\s+", names=["IID", "pop"]
    )
    pop_of = dict(zip(idp["IID"], idp["pop"]))
    d = ROOT / "data" / "aadr"
    # Pooled run: groups the leave-one-out check can't separate are merged into one class.
    pooled = {p: g for g, ps in pools().items() for p in ps}
    pop_pooled = {i: pooled.get(p, p) for i, p in pop_of.items()}
    # A pooled group takes the region of its members (the first member's, if they span regions).
    region_pooled = {**region_of, **{g: region_of.get(ps[0], g) for g, ps in pools().items()}}
    OUT.mkdir(parents=True, exist_ok=True)
    kp = OUT / "aadr_weur_pooled.keep"
    kp.write_text((d / "ref.keep").read_text() + "SELF\tSELF\n")
    run("aadr_weur_pooled", d / "mp", kp, pop_pooled, region_pooled)
    for label, keep in (
        ("aadr_weur", d / "ref.keep"),
        ("aadr_zoom", d / "zoom_clean.keep"),
    ):
        k2 = OUT / f"{label}.keep"
        k2.parent.mkdir(parents=True, exist_ok=True)
        k2.write_text(keep.read_text() + "SELF\tSELF\n")
        run(label, d / "mp", k2, pop_of, region_of)
    psam = pd.read_csv(ROOT / "reference" / "panel" / "all_phase3.psam", sep="\t")
    kg_pop = dict(zip(psam["#IID"], psam["Population"]))
    fam = pd.read_csv(
        ROOT / "data" / "anc" / "mp.fam",
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["FID", "IID"],
    )
    eur = fam[
        fam["IID"].map(psam.set_index("#IID")["SuperPop"]).eq("EUR")
        | fam["IID"].eq("SELF")
    ]
    k3 = OUT / "kg_eur.keep"
    eur.to_csv(k3, sep="\t", header=False, index=False)
    run("kg_eur", ROOT / "data" / "anc" / "mp", k3, kg_pop, {})


if __name__ == "__main__":
    main()
