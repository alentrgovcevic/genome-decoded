"""Read a consumer DNA raw export into one table: rsid, chrom (1-22, X, Y, MT), pos, gt (two letters or "--").

Formats (auto-detected from the first lines):
  myheritage / ftdna  CSV, quoted, header RSID,CHROMOSOME,POSITION,RESULT (FTDNA uses the same layout)
  23andme             tab-separated, '#' comments, columns rsid chromosome position genotype
  ancestrydna         tab-separated, '#' comments, header rsid chromosome position allele1 allele2;
                      chromosomes 23 = X, 24 = Y, 25 = XY (pseudo-autosomal, mapped to X), 26 = MT;
                      no-call allele is '0'
All positions must be GRCh37 (every major consumer export still is); the file header is checked when it says so.
"""

import io
import re
from pathlib import Path

import pandas as pd

ANCESTRY_CHROM = {"23": "X", "24": "Y", "25": "X", "26": "MT"}


def detect(path: Path) -> str:
    head = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for _ in range(40):
            line = fh.readline()
            if not line:
                break
            head.append(line)
    text = "".join(head)
    if "MyHeritage" in text or re.search(
        r'^"?RSID"?,"?CHROMOSOME"?', text, re.M | re.I
    ):
        return "myheritage"
    if (
        re.search(r"^rsid\tchromosome\tposition\tallele1\tallele2", text, re.M | re.I)
        or "AncestryDNA" in text
    ):
        return "ancestrydna"
    if "23andMe" in text or re.search(
        r"^#\s*rsid\tchromosome\tposition\tgenotype", text, re.M
    ):
        return "23andme"
    data = [l for l in head if l.strip() and not l.startswith("#")]
    if data and data[0].count("\t") == 3:
        return "23andme"
    raise ValueError(f"unrecognised raw DNA format: {path.name}")


def build_from_header(path: Path) -> str | None:
    with open(path, encoding="utf-8", errors="replace") as fh:
        for _ in range(40):
            l = fh.readline()
            m = re.search(r"(build ?3[67]|GRCh3[78]|hg1[89])", l, re.I)
            if m:
                return m.group(1)
    return None


def load(path: str | Path) -> tuple[pd.DataFrame, str]:
    path = Path(path)
    fmt = detect(path)
    build = build_from_header(path)
    if build and re.search(r"38|hg38", build, re.I):
        raise ValueError(
            f"{path.name} declares {build}; this pipeline needs GRCh37 coordinates"
        )
    if fmt == "myheritage":
        df = pd.read_csv(path, comment="#", dtype=str)
        df = df.iloc[:, :4]
        df.columns = ["rsid", "chrom", "pos", "gt"]
    elif fmt == "23andme":
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            dtype=str,
            names=["rsid", "chrom", "pos", "gt"],
        )
    else:  # ancestrydna
        lines = [
            l
            for l in open(path, encoding="utf-8", errors="replace")
            if not l.startswith("#")
        ]
        df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
        df.columns = ["rsid", "chrom", "pos", "a1", "a2"]
        df["chrom"] = df["chrom"].map(lambda c: ANCESTRY_CHROM.get(c, c))
        df["gt"] = (df["a1"] + df["a2"]).where(
            (df["a1"] != "0") & (df["a2"] != "0"), "--"
        )
        df = df[["rsid", "chrom", "pos", "gt"]]
    df["chrom"] = (
        df["chrom"]
        .str.upper()
        .str.replace("^CHR", "", regex=True)
        .replace({"23": "X", "24": "Y", "MT": "MT", "M": "MT"})
    )
    df["gt"] = (
        df["gt"]
        .fillna("--")
        .str.upper()
        .str.replace(r"^([ACGTDI])$", r"\1\1", regex=True)
    )  # haploid single letters
    df.loc[~df["gt"].str.fullmatch(r"[ACGTDI]{2}"), "gt"] = "--"
    df = df[pd.to_numeric(df["pos"], errors="coerce").notna()].copy()
    df["pos"] = df["pos"].astype(int)
    return df.reset_index(drop=True), fmt
