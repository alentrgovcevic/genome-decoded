"""Measured hemizygous chrX genotypes in the PharmCAT G6PD region, as a haploid GRCh37 VCF.

G6PD is on X, which is not imputed. For a male sample the chip's X calls are hemizygous, so they can go to
PharmCAT directly as haploid genotypes (GT 0 or 1). Heterozygous X calls in a male are genotyping errors
and are left out; so are sites where the called base is neither 1000G REF nor ALT.
REF/ALT come from the 1000G pvar (data/pgx/kg_chrX_g6pd.pvar, made by scripts/pharmcat.sh).

Usage: g6pd_input.py <regions_hg19.nochr.bed> <kg_chrX.pvar> <out.vcf>
"""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)


def main(bed: str, pvar: str, out: str) -> None:
    fam_sex = (ROOT / "data" / "genome.fam").read_text().split()[4]
    reg = pd.read_csv(
        bed,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chr", "start", "end"],
        dtype={"chr": str},
    )
    reg = reg[reg["chr"] == "X"]
    g = pd.read_csv(
        ROOT / "data" / "genome_23andme.txt",
        sep="\t",
        comment="#",
        names=["rsid", "chrom", "pos", "gt"],
        dtype=str,
    )
    g = g[g["chrom"] == "X"].astype({"pos": int})
    g = g[[any(s < p <= e for s, e in zip(reg["start"], reg["end"])) for p in g["pos"]]]
    pv = pd.read_csv(
        pvar,
        sep="\t",
        comment="#",
        header=None,
        usecols=[1, 3, 4],
        names=["pos", "ref", "alt"],
    )
    pv = pv[pv["ref"].str.len().eq(1) & pv["alt"].str.len().eq(1)].drop_duplicates(
        "pos", keep=False
    )
    m = g.merge(pv, on="pos")
    lines, skipped_het, skipped_allele = [], 0, 0
    for r in m.itertuples(index=False):
        if r.gt[0] != r.gt[1]:
            skipped_het += 1
            continue
        a = r.gt[0]
        code = {r.ref: "0", r.alt: "1"}.get(a)
        if code is None:
            skipped_allele += 1
            continue
        gt = code if fam_sex == "1" else f"{code}/{code}"
        lines.append(f"X\t{r.pos}\t{r.rsid}\t{r.ref}\t{r.alt}\t.\tPASS\t.\tGT\t{gt}")
    with open(out, "w") as fh:
        fh.write(
            '##fileformat=VCFv4.2\n##contig=<ID=X>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSELF\n")
        fh.write("\n".join(lines) + ("\n" if lines else ""))
    print(
        f"G6PD region: {len(g)} chip X sites, {len(m)} in 1000G, written {len(lines)}; "
        f"skipped heterozygous {skipped_het}, allele mismatch {skipped_allele}; haploid={'yes' if fam_sex == '1' else 'no'}"
    )


if __name__ == "__main__":
    main(*sys.argv[1:4])
