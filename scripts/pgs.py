"""Polygenic scores (PGS Catalog, GRCh37-harmonised) for SELF vs 1000 Genomes Europeans.

Primary comparison: SELF and the 503 unrelated 1000G Europeans were imputed *together, the same way*
(scripts/ref_impute.sh: everyone cut down to SELF's chip sites, imputed against 1000G minus those 503), so
both sides are imputed dosages and share any shrinkage. Sensitivity: SELF's imputed dosages vs the 503's
sequenced genotypes (the earlier, asymmetric comparison).

Scores use variants present in both datasets, matched as chr:pos:REF:ALT in either orientation.
Strand-ambiguous SNPs are kept: both sides use 1000G's REF/ALT coding, and PGS Catalog harmonised files are
on the GRCh37 forward strand. Reported per score: percentile, bootstrap 95% interval over the reference
(1,000 resamples), z, and the percentile within each European population.
Output: results/pgs/pgs_summary.json. Run from the run directory after scripts/ref_impute.sh.
"""

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
P2 = str(ROOT / "env-x86" / "bin" / "plink2")
D = ROOT / "data" / "pgs"
R = ROOT / "results" / "pgs"
def scores() -> list[str]:
    """config pgs.both plus the list for the inferred sex (pgs.male / pgs.female)."""
    from run_config import get, sex
    p = get("pgs")
    return p.get("both", []) + p.get(sex(), [])
N_BOOT = 1000


def run(*args: str) -> None:
    subprocess.run([P2, *args, "--silent"], check=True)


def build_refimp() -> None:
    """One pgen with SELF + 503 EUR imputed dosages (from data/refimp/)."""
    if (D / "refimp.pgen").exists():
        return
    parts = []
    for c in range(1, 23):
        out = D / f"refimp_chr{c}"
        run(
            "--vcf",
            str(ROOT / "data" / "refimp" / f"chr{c}.vcf.gz"),
            "dosage=DS",
            "--snps-only",
            "just-acgt",
            "--max-alleles",
            "2",
            "--set-all-var-ids",
            "@:#:$r:$a",
            "--rm-dup",
            "exclude-all",
            "--make-pgen",
            "--out",
            str(out),
        )
        parts.append(str(out))
    (D / "refimp_list.txt").write_text("\n".join(parts) + "\n")
    run(
        "--pmerge-list",
        str(D / "refimp_list.txt"),
        "pfile",
        "--make-pgen",
        "--out",
        str(D / "refimp"),
    )


EUR_SEQ = ROOT / "reference" / "panel" / "eur_seq"  # built once by scripts/build_panel.sh


def build_eur_sequenced() -> None:
    if not Path(f"{EUR_SEQ}.pgen").exists():
        raise SystemExit(f"{EUR_SEQ}.pgen missing: run the workflow's panel rule (scripts/build_panel.sh) first")


def pvar_ids(prefix: Path) -> set:
    ids = set()
    with open(f"{prefix}.pvar") as fh:
        for line in fh:
            if not line.startswith("#"):
                ids.add(line.split("\t", 3)[2])
    return ids


def score_file(pgs: str, shared: set) -> dict:
    fp = ROOT / "reference" / "pgs" / f"{pgs}_hmPOS_GRCh37.txt.gz"
    df = pd.read_csv(fp, sep="\t", comment="#", dtype=str, low_memory=False)
    other = (
        df["other_allele"] if "other_allele" in df else df.get("hm_inferOtherAllele")
    )
    if "hm_inferOtherAllele" in df:
        other = other.fillna(df["hm_inferOtherAllele"])
    df = df.assign(
        chr=df["hm_chr"],
        pos=df["hm_pos"],
        ea=df["effect_allele"],
        oa=other,
        w=pd.to_numeric(df["effect_weight"], errors="coerce"),
    ).dropna(subset=["chr", "pos", "ea", "oa", "w"])
    total_w, n_total = df["w"].abs().sum(), len(df)
    df = df[df["ea"].str.fullmatch("[ACGT]") & df["oa"].str.fullmatch("[ACGT]")]
    id1 = df["chr"] + ":" + df["pos"] + ":" + df["oa"] + ":" + df["ea"]
    id2 = df["chr"] + ":" + df["pos"] + ":" + df["ea"] + ":" + df["oa"]
    df = df.assign(
        id=np.where(id1.isin(shared), id1, np.where(id2.isin(shared), id2, ""))
    )
    kept = df[df["id"] != ""].drop_duplicates("id")
    kept[["id", "ea", "w"]].to_csv(
        D / f"{pgs}.score", sep="\t", index=False, header=["ID", "A1", "WEIGHT"]
    )
    return {
        "n_variants_catalog": n_total,
        "n_used": len(kept),
        "weight_fraction_used": float(kept["w"].abs().sum() / total_w)
        if total_w
        else None,
    }


def score(prefix: Path, pgs: str, out: str) -> pd.DataFrame:
    run(
        "--pfile",
        str(prefix),
        "--score",
        str(D / f"{pgs}.score"),
        "1",
        "2",
        "3",
        "header",
        "cols=+scoresums",
        "no-mean-imputation",
        "--out",
        str(R / out),
    )
    s = pd.read_csv(R / f"{out}.sscore", sep="\t").rename(columns={"#IID": "IID"})
    assert s["ALLELE_CT"].nunique() == 1, (
        f"{out}: samples scored on different variant counts"
    )
    return s


def pct(x: float, dist: np.ndarray) -> float:
    return float((dist < x).mean() * 100)


def main() -> None:
    D.mkdir(parents=True, exist_ok=True)
    R.mkdir(parents=True, exist_ok=True)
    build_refimp()
    build_eur_sequenced()
    shared = pvar_ids(D / "refimp") & pvar_ids(EUR_SEQ)
    pops = pd.read_csv(
        ROOT / "reference" / "panel" / "all_phase3.psam", sep="\t"
    ).set_index("#IID")["Population"]
    rng = np.random.default_rng(20260928)
    summary = {
        "shared_variants": len(shared),
        "method": "SELF and reference both imputed identically (ref_impute.sh)",
        "scores": {},
    }
    SCORES = scores()
    for pgs in SCORES:
        meta = json.loads((ROOT / "reference" / "pgs" / f"{pgs}.json").read_text())
        info = score_file(pgs, shared)
        imp = score(D / "refimp", pgs, f"{pgs}_refimp")
        seq = score(EUR_SEQ, pgs, f"{pgs}_eur")
        col = "SCORE1_SUM"
        mine = float(imp.loc[imp["IID"] == "SELF", col].iloc[0])
        ref = imp[imp["IID"] != "SELF"]
        dist = ref[col].to_numpy()
        boots = [pct(mine, rng.choice(dist, len(dist))) for _ in range(N_BOOT)]
        by_pop = {
            p: pct(mine, g[col].to_numpy())
            for p, g in ref.assign(pop=ref["IID"].map(pops)).groupby("pop")
        }
        pub = meta.get("publication") or {}
        summary["scores"][pgs] = {
            **info,
            "trait": meta["trait_reported"],
            "name": meta.get("name"),
            "publication": f"{pub.get('firstauthor')} {pub.get('date_publication', '')[:4]}, {pub.get('journal')}",
            "doi": pub.get("doi"),
            "self_score": mine,
            "ref_mean": float(dist.mean()),
            "ref_sd": float(dist.std(ddof=1)),
            "z": float((mine - dist.mean()) / dist.std(ddof=1)),
            "percentile": pct(mine, dist),
            "percentile_ci95": [
                float(np.percentile(boots, 2.5)),
                float(np.percentile(boots, 97.5)),
            ],
            "percentile_by_population": by_pop,
            "ref_n": int(len(dist)),
            "percentile_vs_sequenced_reference": pct(mine, seq[col].to_numpy()),
            "allele_ct": int(imp["ALLELE_CT"].iloc[0]),
            "ref_hist": np.histogram(dist, bins=30)[0].tolist(),
            "ref_bins": np.histogram(dist, bins=30)[1].tolist(),
        }
        s = summary["scores"][pgs]
        print(
            pgs,
            meta["trait_reported"],
            f"pct {s['percentile']:.0f} [{s['percentile_ci95'][0]:.0f}-{s['percentile_ci95'][1]:.0f}]",
            f"(vs sequenced {s['percentile_vs_sequenced_reference']:.0f}); by pop",
            {k: round(v) for k, v in by_pop.items()},
            f"weight used {info['weight_fraction_used'] * 100:.0f}%",
        )
    k = len(SCORES)
    summary["multiplicity_note"] = {
        "n_scores": k,
        "expected_max_percentile_under_null": 100 * k / (k + 1),
    }
    (R / "pgs_summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
