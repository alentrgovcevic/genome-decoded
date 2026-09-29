#!/bin/bash
# Reference-level Beagle panels, built once per chromosome (shared by every sample; the workflow runs this
# before any imputation so parallel samples never race on the same file).
# Usage: build_bref3.sh full|noeur <chrom>
#   full   reference/beagle/ref/chrN.bref3 (+ chrN.sites.tsv)  from reference/panel (1000G minus exclude list)
#          for X: non-PAR only (1000G codes PAR as PAR1/PAR2), plus chrX_nonpar pgen and chrX.vcf.gz for the X target
#   noeur  reference/beagle/ref_noeur/chrN.bref3               same panel minus the unrelated Europeans, for the
#          symmetric PGS reference (scripts/ref_impute.sh)
set -euo pipefail  # note: macOS bash 3.2 treats an empty array as unbound, hence ${filter[@]+...} below
P2=${P2:-env-x86/bin/plink2}
BCF=env/bin/bcftools
B=reference/beagle
KG=reference/panel
kind=$1 c=$2
jar="$B/bref3.27Feb25.75f.jar"

if [ "$kind" = full ] && [ "$c" = X ]; then
    O=$B/ref; mkdir -p "$O"
    [ -f "$O/chrX.bref3" ] && exit 0
    $P2 --pfile "$KG/all_phase3" vzs --chr X --max-alleles 2 --mac 2 --nonfounders \
        --rm-dup exclude-all --make-pgen --out "$O/chrX_nonpar" --silent
    $P2 --pfile "$O/chrX_nonpar" --chr X --export vcf bgz --out "$O/chrX.raw" --silent
    $BCF norm -d exact "$O/chrX.raw.vcf.gz" -Oz -o "$O/chrX.vcf.gz"
    $BCF index -f "$O/chrX.vcf.gz"
    $BCF query -i 'TYPE="snp"' -f '%POS\t%REF\t%ALT\n' "$O/chrX.vcf.gz" > "$O/chrX.sites.tsv"
    java -jar "$jar" "$O/chrX.vcf.gz" > "$O/chrX.bref3.tmp" && mv "$O/chrX.bref3.tmp" "$O/chrX.bref3"
    rm "$O/chrX.raw.vcf.gz"
    exit 0
fi

case $kind in
    full)  O=$B/ref;       filter=() ;;
    noeur) O=$B/ref_noeur; filter=(--remove "$KG/eur_unrelated.keep") ;;
    *) echo "usage: build_bref3.sh full|noeur <chrom>" >&2; exit 2 ;;
esac
mkdir -p "$O"
[ -f "$O/chr$c.bref3" ] && exit 0
$P2 --pfile "$KG/all_phase3" vzs --chr "$c" ${filter[@]+"${filter[@]}"} --max-alleles 2 --mac 2 --nonfounders --rm-dup exclude-all \
    --export vcf bgz --out "$O/chr$c.raw" --silent
# --rm-dup works on IDs; identical records with ID "." (chr8:89329148 GTT>G) slip through and make Beagle
# abort with "Duplicate marker", so also drop exact duplicates by position+alleles.
$BCF norm -d exact "$O/chr$c.raw.vcf.gz" -Oz -o "$O/chr$c.vcf.gz"
rm "$O/chr$c.raw.vcf.gz"
if [ "$kind" = full ]; then
    $BCF query -i 'TYPE="snp"' -f '%POS\t%REF\t%ALT\n' "$O/chr$c.vcf.gz" > "$O/chr$c.sites.tsv"
fi
java -jar "$jar" "$O/chr$c.vcf.gz" > "$O/chr$c.bref3.tmp" && mv "$O/chr$c.bref3.tmp" "$O/chr$c.bref3"
rm "$O/chr$c.vcf.gz"
