"""Rescue SNVs that CrossMap drops as Fail(REF==ALT) during GRCh37 -> GRCh38 liftover.

At these sites GRCh38 changed the reference base to what was the GRCh37 ALT allele (e.g. CYP2C19
rs3758581, which PharmCAT needs to call CYP2C19 at all). For each such SNV: lift the position,
check the GRCh38 base equals the old ALT, then swap REF/ALT and recode GT (0<->1), DS (2-DS) and GP
(reversed). Indels are left out and counted.

Usage: rescue_swapped.py <crossmap.unmap> <chain.gz> <GRCh38 fasta> <out.vcf>
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import pysam

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)


def main(unmap: str, chain: str, fasta: str, out: str) -> None:
    header, recs = [], []
    for line in open(unmap):
        if line.startswith("#"):
            header.append(line)
            continue
        f = line.rstrip("\n").split("\t")
        if f[-1] == "Fail(REF==ALT)" and len(f[3]) == 1 and len(f[4]) == 1:
            recs.append(f[:-1])
    n_indel = sum(1 for l in open(unmap) if not l.startswith("#")) - len(recs)
    with tempfile.TemporaryDirectory() as td:
        bed_in, bed_out = Path(td) / "in.bed", Path(td) / "out.bed"
        bed_in.write_text(
            "".join(
                f"{r[0]}\t{int(r[1]) - 1}\t{r[1]}\t{i}\n" for i, r in enumerate(recs)
            )
        )
        subprocess.run(
            [
                str(ROOT / "env" / "bin" / "CrossMap"),
                "bed",
                chain,
                str(bed_in),
                str(bed_out),
            ],
            check=True,
            capture_output=True,
        )
        lifted = {
            int(l.split("\t")[3]): (l.split("\t")[0], int(l.split("\t")[2]))
            for l in bed_out.read_text().splitlines()
        }
    fa = pysam.FastaFile(fasta)
    rescued, mismatch = [], 0
    for i, r in enumerate(recs):
        if i not in lifted:
            continue
        chrom, pos = lifted[i]
        base = fa.fetch(chrom, pos - 1, pos).upper()
        ref, alt = r[3], r[4]
        if base != alt:
            mismatch += 1
            continue
        fmt = r[8].split(":")
        vals = r[9].split(":")
        for k, key in enumerate(fmt):
            if key == "GT":
                vals[k] = "".join({"0": "1", "1": "0"}.get(c, c) for c in vals[k])
            elif key == "DS":
                vals[k] = f"{2 - float(vals[k]):g}"
            elif key == "GP":
                vals[k] = ",".join(reversed(vals[k].split(",")))
        rescued.append([chrom, str(pos), r[2], alt, ref, *r[5:9], ":".join(vals)])
    with open(out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.writelines(
            h
            for h in header
            if h.startswith("##") and not h.startswith(("##fileformat", "##contig"))
        )
        fh.write(
            header[-1]
            if header and header[-1].startswith("#CHROM")
            else "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSELF\n"
        )
        fh.writelines("\t".join(r) + "\n" for r in rescued)
    print(
        f"swapped-reference SNVs rescued: {len(rescued)}; base mismatch after lift: {mismatch}; "
        f"indels/other left out: {n_indel}"
    )


if __name__ == "__main__":
    main(*sys.argv[1:5])
