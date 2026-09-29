"""Measured chip calls at PharmCAT positions that the imputation never saw (site absent from 1000G).

GSA carries custom probes for rare pharmacogenomic alleles; because those sites aren't in the 1000G panel
they never reach data/imputed/, and PharmCAT then treats them as missing (more ambiguous calls, e.g. NAT2,
CYP2C9 *18/*68). Here each such chip site is lifted to GRCh38 and matched to PharmCAT's own position list:
  - both letters equal the GRCh38 reference -> passed to PharmCAT as 0/0 (confirms the allele is absent);
  - any other call -> logged, not passed: rare non-reference array calls are the artefact-prone kind
    (Weedon et al. 2021) and need confirmation.
Usage: pgx_nonpanel.py <regions_hg19.nochr.bed> <imputed_pgx_hg19.vcf.gz> <out_hg38.vcf> <out_log.tsv>
"""

import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from types import SimpleNamespace

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
CHAIN = ROOT / "reference" / "liftover" / "hg19ToHg38.over.chain.gz"
POS = ROOT / "env" / "bin" / "pharmcat_positions.vcf.bgz"


def main(bed: str, imputed: str, out_vcf: str, out_log: str) -> None:
    reg = pd.read_csv(
        bed,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["chr", "start", "end", "gene"],
        dtype={"chr": str},
    )
    g = pd.read_csv(
        ROOT / "data" / "genome_23andme.txt",
        sep="\t",
        comment="#",
        names=["rsid", "chrom", "pos", "call"],
        dtype=str,
    )
    g = g[g["call"] != "--"].astype({"pos": int})
    g = g.merge(reg, left_on="chrom", right_on="chr")
    g = g[(g["pos"] > g["start"]) & (g["pos"] <= g["end"])].drop_duplicates(
        ["chrom", "pos"]
    )
    seen = set()
    with gzip.open(imputed, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                c = line.split("\t", 3)
                seen.add((c[0], int(c[1])))
    g = g[[(c, p) not in seen for c, p in zip(g["chrom"], g["pos"])]].reset_index(
        drop=True
    )
    with tempfile.TemporaryDirectory() as td:
        bi, bo = Path(td) / "in.bed", Path(td) / "out.bed"
        bi.write_text(
            "".join(
                f"chr{c}\t{p - 1}\t{p}\t{i}\n"
                for i, (c, p) in enumerate(zip(g["chrom"], g["pos"]))
            )
        )
        subprocess.run(
            [
                str(ROOT / "env" / "bin" / "CrossMap"),
                "bed",
                str(CHAIN),
                str(bi),
                str(bo),
            ],
            check=True,
            capture_output=True,
        )
        lifted = {
            int(l.split("\t")[3]): (l.split("\t")[0], int(l.split("\t")[2]))
            for l in bo.read_text().splitlines()
        }
    # PharmCAT's own GRCh38 position list (1,202 records). pysam reports this bgzf file as truncated,
    # so it is read once through bcftools and looked up by position.
    q = subprocess.run(
        [
            str(ROOT / "env" / "bin" / "bcftools"),
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\n",
            str(POS),
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    table = {}
    for l in q.strip().splitlines():
        c, p, ref, alt = l.split("\t")
        table.setdefault((c, int(p)), []).append(
            SimpleNamespace(pos=int(p), ref=ref, alts=tuple(alt.split(",")))
        )
    fetch = lambda c, p: table.get((c, p), [])
    # A male's X is haploid: write X genotypes as "0", matching g6pd_input.py (mixed ploidy -> G6PD no call).
    haploid_x = (ROOT / "data" / "genome.fam").read_text().split()[4] == "1"
    ref_gt = lambda c: "0" if haploid_x and c == "chrX" else "0/0"
    recs, log = [], []
    for i, r in g.iterrows():
        if i not in lifted:
            log.append((r.gene, r.rsid, r.chrom, r.pos, r.call, "", "not lifted"))
            continue
        chrom, pos = lifted[i]
        hits = [v for v in fetch(chrom, pos) if v.pos == pos and len(v.ref) == 1]
        if not hits:
            continue  # not a PharmCAT-defining position
        v = hits[0]
        alts = [a for a in v.alts if len(a) == 1]
        if not alts:
            continue
        if r.call[0] == v.ref and r.call[1] == v.ref:
            recs.append(
                f"{chrom}\t{pos}\t{r.rsid}\t{v.ref}\t{alts[0]}\t.\tPASS\tNONPANEL\tGT\t{ref_gt(chrom)}"
            )
            log.append(
                (
                    r.gene,
                    r.rsid,
                    chrom,
                    pos,
                    r.call,
                    v.ref,
                    "reference confirmed, passed",
                )
            )
        else:
            log.append(
                (
                    r.gene,
                    r.rsid,
                    chrom,
                    pos,
                    r.call,
                    v.ref,
                    "non-reference, NOT passed (unconfirmed rare call)",
                )
            )
    with open(out_vcf, "w") as fh:
        fh.write(
            '##fileformat=VCFv4.2\n##INFO=<ID=NONPANEL,Number=0,Type=Flag,Description="Measured chip site absent from 1000G">\n'
        )
        fh.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSELF\n'
        )
        fh.writelines(r + "\n" for r in recs)
    pd.DataFrame(
        log, columns=["gene", "rsid", "chrom", "pos", "chip_gt", "grch38_ref", "action"]
    ).to_csv(out_log, sep="\t", index=False)
    n_nonref = sum(1 for x in log if x[-1].startswith("non-reference"))
    print(
        f"non-panel chip sites in PGx regions: {len(g)}; at PharmCAT positions passed as reference: {len(recs)}; "
        f"non-reference calls held back: {n_nonref}"
    )


if __name__ == "__main__":
    main(*sys.argv[1:5])
