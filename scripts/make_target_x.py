"""Beagle target VCF for chrX (non-PAR). A male SELF is written haploid (GT 0 or 1), a female diploid.

Usage: make_target_x.py <ref_sites.tsv> <eur.afreq> <out.vcf>
Left out: heterozygous calls in a male (errors), sites with alleles outside 1000G REF/ALT, positions with
more than one reference record, and ALT calls at alleles < 1% in 1000G Europeans (chip artefacts; see
scripts/qc_rare_hom.py for the autosomal version).
"""

import sys
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
RARE = 0.01


def main(sites_fp: str, afreq_fp: str, out_fp: str) -> None:
    male = (ROOT / "data" / "genome.fam").read_text().split()[4] == "1"
    g = pd.read_csv(
        ROOT / "data" / "genome_23andme.txt",
        sep="\t",
        comment="#",
        names=["rsid", "chrom", "pos", "gt"],
        dtype={"chrom": str},
    )
    g = g[(g["chrom"] == "X") & (g["gt"] != "--")].drop_duplicates("pos", keep=False)
    ref = pd.read_csv(sites_fp, sep="\t", names=["pos", "ref", "alt"]).drop_duplicates(
        "pos", keep=False
    )
    f = pd.read_csv(afreq_fp, sep="\t")
    f = f[~f["ALT_FREQS"].astype(str).str.contains(",")]
    af = dict(zip(f["ID"], f["ALT_FREQS"].astype(float)))
    m = g.merge(ref, on="pos")
    n = Counter(chip_x=len(g), in_ref=len(m))
    lines = []
    for r in m.itertuples(index=False):
        codes = {r.ref: "0", r.alt: "1"}
        if r.gt[0] not in codes or r.gt[1] not in codes:
            n["allele_mismatch"] += 1
            continue
        a, b = sorted((codes[r.gt[0]], codes[r.gt[1]]))
        if male and a != b:
            n["het_in_male"] += 1
            continue
        freq = af.get(f"X:{r.pos}:{r.ref}:{r.alt}")
        if "1" in (a, b) and freq is not None and freq < RARE:
            n["rare_alt_artefact"] += 1
            continue
        gt = a if male else f"{a}/{b}"
        lines.append(f"X\t{r.pos}\t{r.rsid}\t{r.ref}\t{r.alt}\t.\tPASS\t.\tGT\t{gt}")
        n["written"] += 1
    with open(out_fp, "w") as fh:
        fh.write(
            '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##contig=<ID=X>\n'
        )
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSELF\n")
        fh.write("\n".join(lines) + "\n")
    print("X target:", dict(n), "haploid" if male else "diploid")


if __name__ == "__main__":
    main(*sys.argv[1:4])
