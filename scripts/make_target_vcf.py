"""Write a Beagle target VCF for one chromosome, with REF/ALT taken from the 1000G reference.

Usage: make_target_vcf.py <chrom> <ref_sites.tsv> <out.vcf>
ref_sites.tsv: POS REF ALT for biallelic SNVs in the reference (bcftools query output).
Sites where the called alleles are not within {REF, ALT}, or where the reference has more than
one record at a position, are dropped and counted.
"""

import sys
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)


def main(chrom: str, sites_fp: str, out_fp: str) -> None:
    g = pd.read_csv(
        ROOT / "data" / "genome_23andme.txt",
        sep="\t",
        comment="#",
        names=["rsid", "chrom", "pos", "gt"],
        dtype={"chrom": str},
    )
    g = g[(g["chrom"] == chrom) & (g["gt"] != "--")].drop_duplicates("pos", keep=False)
    ref = pd.read_csv(sites_fp, sep="\t", names=["pos", "ref", "alt"])
    ref = ref.drop_duplicates("pos", keep=False)
    m = g.merge(ref, on="pos")
    n = Counter(target=len(g), in_ref=len(m))
    lines = []
    for r in m.itertuples(index=False):
        codes = {r.ref: "0", r.alt: "1"}
        if r.gt[0] not in codes or r.gt[1] not in codes:
            n["allele_mismatch"] += 1
            continue
        a, b = sorted((codes[r.gt[0]], codes[r.gt[1]]))
        lines.append(
            f"{chrom}\t{r.pos}\t{r.rsid}\t{r.ref}\t{r.alt}\t.\tPASS\t.\tGT\t{a}/{b}"
        )
        n["written"] += 1
    with open(out_fp, "w") as fh:
        fh.write(
            '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        fh.write(
            f"##contig=<ID={chrom}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSELF\n"
        )
        fh.write("\n".join(lines) + "\n")
    print(chrom, dict(n))


if __name__ == "__main__":
    main(*sys.argv[1:4])
