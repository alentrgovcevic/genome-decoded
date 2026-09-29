#!/bin/bash
# Local imputation of autosomes with Beagle 5.5 against 1000 Genomes phase 3 (GRCh37, phased).
# Per chromosome: export reference (MAC >= 2, biallelic) -> bref3 -> harmonised target -> Beagle gp=true.
# Output: data/imputed/chrN.vcf.gz (+ .csi). Log: data/imputed/impute.log. Run from the run directory runs/<sample>/.
# Usage: bash scripts/impute.sh [chrom ...]   (default 1..22)
set -euo pipefail

P2=env-x86/bin/plink2
BCF=env/bin/bcftools
JAVA=java
B=reference/beagle
KG=reference/panel
O=data/imputed
mkdir -p "$O" "$B/ref"
CHROMS=${*:-$(seq 1 22)}

for c in $CHROMS; do
    [ -f "$O/chr$c.vcf.gz" ] && { echo "chr$c done, skipping"; continue; }
    echo "== chr$c $(date +%H:%M:%S)"
    [ -f "$B/ref/chr$c.bref3" ] || { echo "missing $B/ref/chr$c.bref3: run scripts/build_bref3.sh full $c"; exit 1; }
    env/bin/python scripts/make_target_vcf.py "$c" "$B/ref/chr$c.sites.tsv" "$O/target.chr$c.vcf"
    $JAVA -Xmx16g -jar env/share/beagle-5.5*/beagle.jar ref="$B/ref/chr$c.bref3" gt="$O/target.chr$c.vcf" \
        map="$B/plink.chr$c.GRCh37.map" out="$O/chr$c" gp=true nthreads=8 > "$O/beagle.chr$c.out" 2>&1
    $BCF index -f "$O/chr$c.vcf.gz"
done
echo "== all done $(date +%H:%M:%S)"
