"""Filter the PharmCAT input for uninformative, implausible or low-confidence genotypes.

1. Imputed records where >= 99% of 1000G haplotypes carry ALT: Beagle imputes ALT/ALT for everyone,
   so the call says nothing about this person. Dropped.
2. Measured records homozygous for an allele with < 1% frequency in 1000G. On consumer arrays these are
   far more often probe artefacts than real (CYP2C19 rs367543002/3 read T/T and C/C at 0.26% variants,
   which fits no CYP2C19 haplotype). Set to missing.
3. Imputed genotypes with max genotype probability < 0.99. The masked check found 99% accuracy at
   >= 0.99 and 86% at 0.9-0.99, too low for a drug-dosing input. Set to missing.
Every action is written to a provenance TSV for the report.

Usage: pgx_site_filter.py <kg.afreq> <in.vcf.gz> <out.vcf> <provenance.tsv>
kg.afreq: plink2 --freq output with IDs chr:pos:ref:alt.
"""

import gzip
import sys

import pandas as pd

MAX_ALT_FREQ = 0.99
RARE = 0.01
MIN_GP = 0.99


def main(afreq: str, vcf_in: str, vcf_out: str, prov_out: str) -> None:
    f = pd.read_csv(afreq, sep="\t")
    f = f[~f["ALT_FREQS"].astype(str).str.contains(",")]  # biallelic only
    af = dict(zip(f["ID"], f["ALT_FREQS"].astype(float)))
    n = {"kept": 0, "dropped_fixed_alt": 0, "masked_rare_hom": 0, "masked_low_gp": 0}
    prov = []
    with gzip.open(vcf_in, "rt") as fin, open(vcf_out, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, rid, ref, alt, info, fmt, sample = (
                fields[0],
                fields[1],
                fields[2],
                fields[3],
                fields[4],
                fields[7],
                fields[8],
                fields[9],
            )
            freq = af.get(f"{chrom}:{pos}:{ref}:{alt}")
            imputed = "IMP" in info.split(";")
            if freq is not None and imputed and freq >= MAX_ALT_FREQ:
                n["dropped_fixed_alt"] += 1
                prov.append(
                    (
                        chrom,
                        pos,
                        rid,
                        "dropped: imputed at a site fixed for ALT in 1000G",
                    )
                )
                continue
            keys, vals = fmt.split(":"), sample.split(":")
            gt = vals[0].replace("|", "/")
            action = None
            if (
                freq is not None
                and not imputed
                and ((gt == "1/1" and freq < RARE) or (gt == "0/0" and 1 - freq < RARE))
            ):
                action = f"masked: measured homozygote for a <1% allele (1000G ALT freq {freq:.4f})"
                n["masked_rare_hom"] += 1
            elif imputed and "GP" in keys:
                gp = max(map(float, vals[keys.index("GP")].split(",")))
                if gp < MIN_GP:
                    action = (
                        f"masked: imputed with genotype probability {gp:.2f} < {MIN_GP}"
                    )
                    n["masked_low_gp"] += 1
            if action:
                vals[0] = "./."
                fields[9] = ":".join(vals)
                line = "\t".join(fields) + "\n"
                prov.append((chrom, pos, rid, action))
            fout.write(line)
            n["kept"] += 1
    pd.DataFrame(prov, columns=["chrom", "pos_grch37", "rsid", "action"]).to_csv(
        prov_out, sep="\t", index=False
    )
    print("PGx site filter:", n)


if __name__ == "__main__":
    main(*sys.argv[1:5])
