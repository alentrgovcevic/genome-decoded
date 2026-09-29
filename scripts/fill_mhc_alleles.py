"""Fill the missing second allele ('0') at homozygous MHC sites of data/genome.bim from 1000G alleles.

In a single-sample plink1 .bim a homozygous site is written '0 C'; HIBAG then cannot match it to a
model SNP's allele pair. Replacing '0' leaves the genotype unchanged (the .bed code still means
homozygous for the listed allele). Writes data/hla/mhc.{bed,bim,fam}.
Needs data/hla/kg_mhc.pvar (plink2 --make-just-pvar of the 1000G MHC region).
"""

import shutil
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
D = ROOT / "data" / "hla"


def main() -> None:
    bim = pd.read_csv(
        ROOT / "data" / "genome.bim",
        sep="\t",
        header=None,
        names=["chr", "id", "cm", "pos", "a1", "a2"],
        dtype={"chr": str, "a1": str, "a2": str},
    )
    pvar = pd.read_csv(
        D / "kg_mhc.pvar",
        sep="\t",
        comment="#",
        header=None,
        usecols=[1, 3, 4],
        names=["pos", "ref", "alt"],
    )
    pvar = pvar[
        pvar["ref"].str.len().eq(1) & pvar["alt"].str.len().eq(1)
    ].drop_duplicates("pos", keep=False)
    alleles = {p: (r, a) for p, r, a in pvar.itertuples(index=False)}
    filled = unresolved = 0
    in_mhc = (bim["chr"] == "6") & bim["pos"].between(25_000_000, 35_000_000)
    for i in bim.index[in_mhc & (bim["a1"] == "0")]:
        pair = alleles.get(bim.at[i, "pos"])
        present = bim.at[i, "a2"]
        if pair and present in pair:
            bim.at[i, "a1"] = pair[1] if present == pair[0] else pair[0]
            filled += 1
        else:
            unresolved += 1
    bim.to_csv(D / "mhc.bim", sep="\t", header=False, index=False)
    for ext in ("bed", "fam"):
        shutil.copy(ROOT / "data" / f"genome.{ext}", D / f"mhc.{ext}")
    print(f"filled {filled}, unresolved {unresolved}")


if __name__ == "__main__":
    main()
