"""QC of a consumer DNA raw export (see parse_input.py for formats) and conversion to 23andMe format for plink.

Usage: qc.py <raw_file>   (run from the run directory)
Writes data/genome_23andme.txt and results/qc.json. MT rows are counted but not carried forward.
"""

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
RAW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
CHROMS = [str(c) for c in range(1, 23)] + ["X", "Y"]


def load() -> pd.DataFrame:
    from parse_input import load as parse

    if RAW is None:
        raise SystemExit("usage: qc.py <raw_file>")
    df, fmt = parse(RAW)
    load.fmt = fmt
    load.n_mt = int((df["chrom"] == "MT").sum())
    return df[df["chrom"] != "MT"].reset_index(drop=True)


def classify(gt: str) -> str:
    if gt == "--":
        return "no_call"
    if set(gt) & {"D", "I"}:
        return "indel"
    return "hom" if gt[0] == gt[1] else "het"


def infer_sex(x_het: float, y_call: float) -> str:
    """Sex from chrX heterozygosity and chrY call rate; 'unknown' between the thresholds."""
    if x_het < 0.05 and y_call > 0.8:
        return "male"
    if x_het > 0.1 and y_call < 0.5:
        return "female"
    return "unknown"


def main() -> None:
    df = load()
    df["cls"] = df["gt"].map(classify)
    unknown_chroms = sorted(set(df["chrom"]) - set(CHROMS))
    dup_rsid = int(df["rsid"].duplicated().sum())
    dup_pos = int(df.duplicated(["chrom", "pos"]).sum())

    per_chrom = []
    for c in CHROMS:
        sub = df[df["chrom"] == c]
        n = len(sub)
        counts = sub["cls"].value_counts().to_dict()
        called_snv = counts.get("hom", 0) + counts.get("het", 0)
        per_chrom.append(
            {
                "chrom": c,
                "n": n,
                "no_call": counts.get("no_call", 0),
                "indel": counts.get("indel", 0),
                "hom": counts.get("hom", 0),
                "het": counts.get("het", 0),
                "call_rate": 1 - counts.get("no_call", 0) / n if n else None,
                "het_rate": counts.get("het", 0) / called_snv if called_snv else None,
                "min_pos": int(sub["pos"].min()) if n else None,
                "max_pos": int(sub["pos"].max()) if n else None,
            }
        )

    totals = df["cls"].value_counts().to_dict()
    auto = df[df["chrom"].isin(CHROMS[:22])]
    auto_called = auto["cls"].isin(["hom", "het"]).sum()
    gt_counts = df.loc[df["cls"].isin(["hom", "het"]), "gt"].value_counts().to_dict()

    qc = {
        "raw_file": RAW.name,
        "input_format": load.fmt,
        "mt_rows_ignored": load.n_mt,
        "n_variants": len(df),
        "class_totals": {k: int(v) for k, v in totals.items()},
        "call_rate": 1 - totals.get("no_call", 0) / len(df),
        "autosomal_het_rate": float((auto["cls"] == "het").sum() / auto_called),
        "genotype_counts": {k: int(v) for k, v in gt_counts.items()},
        "per_chrom": per_chrom,
        "unknown_chroms": unknown_chroms,
        "duplicate_rsids": dup_rsid,
        "duplicate_positions": dup_pos,
    }
    # Sex from the data, passed to plink explicitly: plink --23file calls anyone with Y calls male, and female
    # arrays carry noise calls on Y probes (one female export had 9.3% of Y sites called).
    x = next(c for c in per_chrom if c["chrom"] == "X")
    y = next(c for c in per_chrom if c["chrom"] == "Y")
    x_het, y_call = (x["het_rate"] or 0.0), (y["call_rate"] or 0.0)
    qc["inferred_sex"] = infer_sex(x_het, y_call)
    qc["sex_evidence"] = {"x_het_rate": x_het, "y_call_rate": y_call}

    out_data = ROOT / "data"
    out_res = ROOT / "results"
    out_data.mkdir(exist_ok=True)
    out_res.mkdir(exist_ok=True)
    (out_data / "sex.txt").write_text(
        {"male": "m", "female": "f"}.get(qc["inferred_sex"], "0") + "\n"
    )

    keep = df[df["cls"] != "indel"].copy()
    if (
        qc["inferred_sex"] == "female"
    ):  # Y calls in a female are probe noise; plink rejects them
        yc = (keep["chrom"] == "Y") & (keep["gt"] != "--")
        qc["female_y_calls_blanked"] = int(yc.sum())
        keep.loc[keep["chrom"] == "Y", "gt"] = "--"
    with open(out_data / "genome_23andme.txt", "w") as fh:
        fh.write("# rsid\tchromosome\tposition\tgenotype\n")
        keep[["rsid", "chrom", "pos", "gt"]].to_csv(
            fh, sep="\t", header=False, index=False
        )
    qc["n_written_23andme"] = len(keep)

    (out_res / "qc.json").write_text(json.dumps(qc, indent=1))
    print(
        json.dumps(
            {k: v for k, v in qc.items() if k not in ("per_chrom", "genotype_counts")},
            indent=1,
        )
    )


if __name__ == "__main__":
    sys.exit(main())
