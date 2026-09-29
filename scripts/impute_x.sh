#!/bin/bash
# Impute the non-pseudo-autosomal X chromosome with Beagle 5.5 against 1000 Genomes phase 3.
# PAR regions are excluded: the plink2 1000G files already code them as PAR1/PAR2, separate from X.
# A male SELF is haploid on X: calls are written as single alleles (0/1); heterozygous X calls (genotyping
# errors in a male) and ALT calls at alleles < 1% in 1000G Europeans (chip artefacts, as on the autosomes)
# are left out. Output: data/imputed/chrX.vcf.gz. Run from the run directory runs/<sample>/.
set -euo pipefail

P2=env-x86/bin/plink2
BCF=env/bin/bcftools
B=reference/beagle
KG=reference/panel
O=data/imputed
mkdir -p "$O" "$B/ref"

[ -f "$B/ref/chrX.bref3" ] || { echo "missing $B/ref/chrX.bref3: run scripts/build_bref3.sh full X"; exit 1; }

# EUR ALT frequencies on X, for the artefact filter.
$P2 --pfile "$B/ref/chrX_nonpar" --chr X --keep reference/panel/eur_unrelated.keep --snps-only just-acgt --nonfounders \
    --set-all-var-ids '@:#:$r:$a' --freq --out "$O/chrX_eur" --silent
env/bin/python scripts/make_target_x.py "$B/ref/chrX.sites.tsv" "$O/chrX_eur.afreq" "$O/target.chrX.vcf"
java -Xmx16g -jar env/share/beagle-5.5*/beagle.jar ref="$B/ref/chrX.bref3" gt="$O/target.chrX.vcf" \
    map="$B/plink.chrX.GRCh37.map" out="$O/chrX" gp=true nthreads=10 seed=20260928 > "$O/beagle.chrX.out" 2>&1 \
    || { echo "Beagle failed on chrX; see $O/beagle.chrX.out"; exit 1; }
$BCF index -f "$O/chrX.vcf.gz"
echo "chrX done: $($BCF view -H "$O/chrX.vcf.gz" | wc -l | tr -d ' ') records"
