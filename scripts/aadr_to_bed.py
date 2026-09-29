"""Convert selected individuals from an AADR transposed packed genotype file (TGENO) to PLINK 1 .bed/.bim/.fam.

TGENO layout (EIGENSOFT, as found in AADR v66): a 48-byte header "TGENO <nind> <nsnp> <hash> <hash>"
(file size = 48 + nind * rlen, checked), then one record per
individual, each max(48, ceil(nsnp/4)) bytes, four SNPs per byte, most significant bits first.
A 2-bit value is the count of the .snp file's first (reference) allele; 3 = missing.
convertf 8.0.0 from bioconda refuses this layout ("bad ID (SNP): TGENO ...").

Usage: aadr_to_bed.py <prefix> <poplist> <out_prefix>
"""

import math
import os
import sys

import numpy as np
import pandas as pd

# Count of reference allele (0,1,2,3=missing) -> PLINK 2-bit code with A1 = alt, A2 = ref:
# 0 ref copies = hom A1 (0b00), 1 = het (0b10), 2 = hom A2 (0b11), missing = 0b01.
TO_PLINK = np.array([0b00, 0b10, 0b11, 0b01], dtype=np.uint8)


def main(prefix: str, poplist: str, out: str) -> None:
    ind = pd.read_csv(
        f"{prefix}.ind", sep=r"\s+", header=None, names=["id", "sex", "pop"]
    )
    snp = pd.read_csv(
        f"{prefix}.snp",
        sep=r"\s+",
        header=None,
        names=["id", "chr", "cm", "pos", "ref", "alt"],
    )
    pops = set(open(poplist).read().split())
    with open(f"{prefix}.geno", "rb") as fh:
        head = fh.read(48).split(b"\0")[0].split()
        assert head[0] == b"TGENO", head[0]
        nind, nsnp = int(head[1]), int(head[2])
        assert (nind, nsnp) == (len(ind), len(snp)), (nind, nsnp, len(ind), len(snp))
        rlen = max(48, math.ceil(nsnp / 4))
        size = os.path.getsize(f"{prefix}.geno")
        assert size == 48 + nind * rlen, f"unexpected TGENO size {size}"
        # Only Human Origins array genotypes (".HO"): shotgun genomes (".DG", ".SG") form a technology batch
        # that the regional zoom PCA picked out as outliers.
        sel = ind.index[ind["pop"].isin(pops) & ind["id"].str.endswith(".HO")].to_numpy()
        g = np.empty((len(sel), nsnp), dtype=np.uint8)
        for k, i in enumerate(sel):
            fh.seek(48 + rlen * i)
            rec = np.frombuffer(fh.read(rlen), dtype=np.uint8)
            vals = np.stack([(rec >> s) & 3 for s in (6, 4, 2, 0)], axis=1).reshape(-1)[
                :nsnp
            ]
            g[k] = vals
    # SNP-major PLINK .bed: magic, then per SNP ceil(n/4) bytes, sample order low bits first.
    codes = TO_PLINK[g.T]  # nsnp x nsel
    n = codes.shape[1]
    pad = (-n) % 4
    if pad:
        codes = np.concatenate([codes, np.zeros((nsnp, pad), dtype=np.uint8)], axis=1)
    c = codes.reshape(nsnp, -1, 4)
    packed = (
        c[:, :, 0] | (c[:, :, 1] << 2) | (c[:, :, 2] << 4) | (c[:, :, 3] << 6)
    ).astype(np.uint8)
    with open(f"{out}.bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        fh.write(packed.tobytes())
    chrom = snp["chr"].astype(str).replace({"23": "X", "24": "Y", "90": "MT"})
    pd.DataFrame(
        {
            "chr": chrom,
            "id": snp["id"],
            "cm": 0,
            "pos": snp["pos"],
            "a1": snp["alt"],
            "a2": snp["ref"],
        }
    ).to_csv(f"{out}.bim", sep="\t", header=False, index=False)
    sexc = ind.loc[sel, "sex"].map({"M": 1, "F": 2}).fillna(0).astype(int)
    pd.DataFrame(
        {
            "fid": ind.loc[sel, "pop"],
            "iid": ind.loc[sel, "id"],
            "p": 0,
            "m": 0,
            "sex": sexc,
            "ph": -9,
        }
    ).to_csv(f"{out}.fam", sep="\t", header=False, index=False)
    miss = float((g == 3).mean())
    print(
        f"{len(sel)} individuals from {len(pops)} populations, {nsnp} SNPs, missing rate {miss:.4f}"
    )


if __name__ == "__main__":
    main(*sys.argv[1:4])
