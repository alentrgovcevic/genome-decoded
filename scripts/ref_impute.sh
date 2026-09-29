#!/bin/bash
# Symmetric reference for PGS and an imputation benchmark with truth.
# Per chromosome: the 503 unrelated 1000G Europeans are cut down to SELF's chip sites (phase removed) and
# imputed together with SELF against a panel of 1000G *minus* those 503, so SELF and the reference are
# processed identically (dosage vs dosage), and the 503 give truth for non-chip sites.
# Output: data/refimp/chrN.vcf.gz (SELF + 503 EUR). Log: data/refimp/refimp.log. Run from the run directory runs/<sample>/.
# Usage: bash scripts/ref_impute.sh [chrom ...]   (default 1..22)
set -euo pipefail

P2=env-x86/bin/plink2
BCF=env/bin/bcftools
B=reference/beagle
KG=reference/panel
O=data/refimp
mkdir -p "$O" "$B/ref_noeur"
cp reference/panel/eur_unrelated.keep "$O/eur503.txt"
CHROMS=${*:-$(seq 1 22)}

for c in $CHROMS; do
    [ -f "$O/chr$c.vcf.gz.csi" ] && { echo "chr$c done, skipping"; continue; }
    echo "== chr$c $(date +%H:%M:%S)"
    [ -f "$B/ref_noeur/chr$c.bref3" ] || { echo "missing $B/ref_noeur/chr$c.bref3: run scripts/build_bref3.sh noeur $c"; exit 1; }
    # SELF target (already harmonised to 1000G REF/ALT) and its positions.
    # Drop SELF's measured homozygous calls for alleles < 1% in 1000G Europeans (chip artefacts, ~1,400
    # genome-wide; scripts/qc_rare_hom.py) so they don't steer the imputation.
    awk -v c="$c" 'NR>1 && $2==c {print $2"\t"$3}' results/qc_rare_hom.tsv > "$O/rarehom.chr$c.tsv"
    # bcftools refuses an empty -T file; a chromosome with no rare homozygotes (common for clean input) needs none
    excl=(); [ -s "$O/rarehom.chr$c.tsv" ] && excl=(-T ^"$O/rarehom.chr$c.tsv")
    $BCF view data/imputed/target.chr$c.vcf ${excl[@]+"${excl[@]}"} -Oz -o "$O/self.chr$c.vcf.gz" \
        && $BCF index -f "$O/self.chr$c.vcf.gz"
    $BCF query -f '%CHROM\t%POS\n' "$O/self.chr$c.vcf.gz" > "$O/pos.chr$c.tsv"
    awk '{print $1"\t"$2"\t"$2}' "$O/pos.chr$c.tsv" > "$O/pos.chr$c.range"
    # EUR truth genotypes at exactly those sites, phase stripped (the chip gives unphased calls).
    $P2 --pfile "$KG/all_phase3" vzs --chr "$c" --keep "$O/eur503.txt" --extract range "$O/pos.chr$c.range" \
        --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all --export vcf bgz --out "$O/eur.chr$c.raw" --silent
    $BCF view "$O/eur.chr$c.raw.vcf.gz" | sed '/^#/!s/|/\//g' | $BCF view -Oz -o "$O/eur.chr$c.vcf.gz"
    $BCF index -f "$O/eur.chr$c.vcf.gz"
    $BCF merge -m none "$O/self.chr$c.vcf.gz" "$O/eur.chr$c.vcf.gz" | $BCF view -T "$O/pos.chr$c.tsv" -m2 -M2 -Oz \
        -o "$O/target.chr$c.vcf.gz"
    java -Xmx16g -jar env/share/beagle-5.5*/beagle.jar ref="$B/ref_noeur/chr$c.bref3" gt="$O/target.chr$c.vcf.gz" \
        map="$B/plink.chr$c.GRCh37.map" out="$O/chr$c" gp=true nthreads=10 seed=20260928 > "$O/beagle.chr$c.out" 2>&1 \
        || { echo "Beagle failed on chr$c; see $O/beagle.chr$c.out"; exit 1; }
    $BCF index -f "$O/chr$c.vcf.gz"
    rm -f "$O/eur.chr$c.raw.vcf.gz"
done
echo "== all done $(date +%H:%M:%S)"
