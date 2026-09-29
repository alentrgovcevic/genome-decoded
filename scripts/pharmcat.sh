#!/bin/bash
# Pharmacogenomic diplotypes with PharmCAT from the imputed genome.
#   PharmCAT gene regions (GRCh38) -> lifted to GRCh37 -> extracted from data/imputed/ -> chr-renamed
#   -> CrossMap to GRCh38 -> PharmCAT preprocessor -> PharmCAT with HLA-A/B outside calls from HIBAG.
# Needs: data/imputed/chr*.vcf.gz, results/hla/hla_calls.tsv, reference/liftover/*.chain.gz,
#        env/bin/reference.fna.bgz (PharmCAT's GRCh38 FASTA). Run from the run directory runs/<sample>/.
set -euo pipefail
export PATH="$PWD/env/bin:$PATH"   # the preprocessor looks up bcftools/bgzip on PATH

BCF=env/bin/bcftools
PY=env/bin/python
O=data/pgx
R=results/pharmcat
mkdir -p "$O" "$R"
FASTA=env/bin/reference.fna.bgz

# 1. Gene regions to GRCh37 (plain chromosome names, to match the imputed VCFs).
env/bin/CrossMap bed reference/liftover/hg38ToHg19.over.chain.gz env/bin/pharmcat_regions.bed "$O/regions_hg19.bed"
sed 's/^chr//' "$O/regions_hg19.bed" | sort -k1,1V -k2,2n > "$O/regions_hg19.nochr.bed"

# 2. Extract those regions from each imputed chromosome, keep IMP/DR2, concatenate.
: > "$O/parts.txt"
for c in $(cut -f1 "$O/regions_hg19.nochr.bed" | sort -u | grep -E '^[0-9]+$'); do
    [ -f "data/imputed/chr$c.vcf.gz.csi" ] || { echo "WARNING: chr$c not imputed yet, its genes are missing"; continue; }
    $BCF view -R "$O/regions_hg19.nochr.bed" "data/imputed/chr$c.vcf.gz" -Oz -o "$O/pgx_chr$c.vcf.gz"
    $BCF index -f "$O/pgx_chr$c.vcf.gz"
    echo "$O/pgx_chr$c.vcf.gz" >> "$O/parts.txt"
done
$BCF concat -f "$O/parts.txt" -a -Oz -o "$O/pgx_hg19_all.vcf.gz"

# Drop *imputed* sites where >= 99% of 1000G haplotypes carry ALT: Beagle then imputes ALT/ALT for
# everyone, which says nothing about this person and can contradict PharmCAT's definitions
# Also mask measured homozygous calls for alleles < 1% in 1000G (array artefacts; see pgx_site_filter.py).
awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' "$O/regions_hg19.nochr.bed" > "$O/regions_hg19.range"
env-x86/bin/plink2 --pfile reference/panel/all_phase3 vzs --extract range "$O/regions_hg19.range" --nonfounders \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 missing --freq --out "$O/kg_pgx" --silent
$PY scripts/pgx_site_filter.py "$O/kg_pgx.afreq" "$O/pgx_hg19_all.vcf.gz" "$O/pgx_hg19.vcf" "$R/site_filter.tsv" | tee "$O/fixed_alt.log"
# G6PD (chrX, not imputed): add the measured hemizygous X calls as haploid genotypes.
awk '$1=="X"{print $1"\t"$2"\t"$3}' "$O/regions_hg19.nochr.bed" > "$O/g6pd.range"
env-x86/bin/plink2 --pfile reference/panel/all_phase3 vzs --extract range "$O/g6pd.range" --make-just-pvar cols= \
    --out "$O/kg_chrX_g6pd" --silent
$PY scripts/g6pd_input.py "$O/regions_hg19.nochr.bed" "$O/kg_chrX_g6pd.pvar" "$O/pgx_chrX_hg19.vcf" | tee "$O/g6pd.log"
$BCF view "$O/pgx_hg19.vcf" -Oz -o "$O/pgx_auto_hg19.vcf.gz" && $BCF index -f "$O/pgx_auto_hg19.vcf.gz"
$BCF view "$O/pgx_chrX_hg19.vcf" -Oz -o "$O/pgx_chrX_hg19.vcf.gz" && $BCF index -f "$O/pgx_chrX_hg19.vcf.gz"
$BCF concat -a "$O/pgx_auto_hg19.vcf.gz" "$O/pgx_chrX_hg19.vcf.gz" -Oz -o "$O/pgx_hg19.vcf.gz"

# 3. chr-prefix, lift to GRCh38, sort.
for c in $(seq 1 22) X Y; do echo -e "$c\tchr$c"; done > "$O/chrmap.txt"
$BCF annotate --rename-chrs "$O/chrmap.txt" "$O/pgx_hg19.vcf.gz" -Oz -o "$O/pgx_hg19_chr.vcf.gz"
env/bin/CrossMap vcf reference/liftover/hg19ToHg38.over.chain.gz "$O/pgx_hg19_chr.vcf.gz" "$FASTA" "$O/pgx_hg38.vcf" > "$O/crossmap.log" 2>&1
$PY scripts/rescue_swapped.py "$O/pgx_hg38.vcf.unmap" reference/liftover/hg19ToHg38.over.chain.gz "$FASTA" \
    "$O/rescued.vcf" | tee "$O/rescue.log"
$PY scripts/pgx_nonpanel.py "$O/regions_hg19.nochr.bed" "$O/pgx_hg19_all.vcf.gz" "$O/nonpanel.vcf" "$R/nonpanel_calls.tsv" \
    | tee "$O/nonpanel.log"
# Merge lifted + rescued + non-panel records explicitly (one header, records from all three), then sort.
{ grep '^##' "$O/pgx_hg38.vcf"; grep '^##INFO=<ID=NONPANEL' "$O/nonpanel.vcf"; grep '^#CHROM' "$O/pgx_hg38.vcf"
  # One record per position, priority lifted > rescued > non-panel (duplicates crash the preprocessor).
  grep -hv '^#' "$O/pgx_hg38.vcf" "$O/rescued.vcf" "$O/nonpanel.vcf" | awk -F'\t' '!seen[$1":"$2]++'; } \
  | $BCF sort -Oz -o "$O/pgx_hg38.vcf.gz"
$BCF index -f "$O/pgx_hg38.vcf.gz"

# Beagle's phase is statistical (one sample, no parents), so hand PharmCAT unphased genotypes and let it
# enumerate diplotypes. With the imputed phase, CYP2C19 came back "no call".
zcat < "$O/pgx_hg38.vcf.gz" | awk 'BEGIN{OFS="\t"} /^#/{print;next} {split($10,f,":"); gsub(/\|/,"/",f[1]); $10=f[1]; for(i=2;i in f;i++) $10=$10":"f[i]; print}' \
    | $BCF view -Oz -o "$O/pgx_hg38.unphased.vcf.gz"
mv "$O/pgx_hg38.unphased.vcf.gz" "$O/pgx_hg38.vcf.gz"
$BCF index -f "$O/pgx_hg38.vcf.gz"

# 4. HLA-A / HLA-B outside calls from HIBAG (two-field alleles).
awk -F'\t' 'NR>1 && ($1=="A"||$1=="B") {print "HLA-"$1"\t*"$2"/*"$3}' results/hla/hla_calls.tsv > "$R/outside_calls.tsv"

# 5. PharmCAT preprocessor, then PharmCAT itself (matcher + phenotyper + reporter).
$PY env/bin/pharmcat_vcf_preprocessor -vcf "$O/pgx_hg38.vcf.gz" -o "$O" -bf pgx
java -jar env/bin/pharmcat.jar -vcf "$O/pgx.preprocessed.vcf.bgz" -po "$R/outside_calls.tsv" \
    -o "$R" -bf pharmcat -reporterJson -reporterHtml -reporterCallsOnlyTsv
# CYP2D6 research mode, kept out of the main report: SNP-only calling, blind to *5 deletions and duplications.
java -jar env/bin/pharmcat.jar -vcf "$O/pgx.preprocessed.vcf.bgz" -research cyp2d6 -g CYP2D6 \
    -o "$R/research" -bf cyp2d6 -reporterJson -reporterCallsOnlyTsv
ls -la "$R"
# Every diplotype consistent with the unphased genotypes (e.g. TPMT *1/*3A vs *3B/*3C), for the report.
java -jar env/bin/pharmcat.jar -matcher -vcf "$O/pgx.preprocessed.vcf.bgz" -ma -o "$R/all" -bf all
