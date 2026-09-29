#!/bin/bash
# West-Eurasian ancestry against present-day Human Origins samples from AADR v66 (Reich lab, GRCh37).
#   scripts/aadr_to_bed.py (selected populations only) -> autosomal SNVs, IDs chr:pos -> strand check vs SELF
#   -> merge -> LD prune on the reference -> projected PCA (all, and a regional zoom from config) -> results/aadr/.
# Needs data/anc/self.{bed,bim,fam} from scripts/ancestry.sh. Run from the run directory runs/<sample>/.
set -euo pipefail

P1=env/bin/plink
P2=env-x86/bin/plink2
A=reference/aadr
D=data/aadr
R=results/aadr
mkdir -p "$D" "$R"

# Present-day HO groups; Reich-lab outlier/QC-flagged labels (-o, -QCremove, -Dup, -Discovery, -oPCA,
# _WGA) are deliberately not listed; single-person groups are left out.
# Reference populations = every population listed under ancestry.regions in the config.
env/bin/python scripts/run_config.py ancestry.regions > "$D/pops.txt"

# convertf 8.0.0 cannot read the TGENO layout AADR v66 ships, so convert with our own reader.
[ -f "$D/ho.bed" ] || env/bin/python scripts/aadr_to_bed.py "$A/v66.p1_HO.aadr.patch.PUB" "$D/pops.txt" "$D/ho"

# Autosomal biallelic ACGT SNVs, IDs chr:pos; keep sites also in SELF; drop strand-ambiguous pairs.
$P2 --bfile "$D/ho" --autosome --snps-only just-acgt --max-alleles 2 --set-all-var-ids @:# --rm-dup exclude-all \
    --make-bed --out "$D/ho_ids" --silent
env/bin/python scripts/aadr_strand.py "$D/ho_ids.bim" data/anc/self.bim "$D/keep.ids" "$R/strand.json"
$P1 --bfile "$D/ho_ids" --extract "$D/keep.ids" --geno 0.05 --make-bed --out "$D/ho_k" --silent
$P1 --bfile data/anc/self --extract "$D/keep.ids" --make-bed --out "$D/self_k" --silent
if ! $P1 --bfile "$D/ho_k" --bmerge "$D/self_k" --make-bed --out "$D/merged" --silent; then
    $P1 --bfile "$D/ho_k" --exclude "$D/merged-merge.missnp" --make-bed --out "$D/ho_k2" --silent
    $P1 --bfile "$D/self_k" --exclude "$D/merged-merge.missnp" --make-bed --out "$D/self_k2" --silent
    $P1 --bfile "$D/ho_k2" --bmerge "$D/self_k2" --make-bed --out "$D/merged" --silent
fi

# LD prune on the reference, then founders-only fam (convertf writes no parents; be safe).
$P2 --bfile "$D/ho_k" --indep-pairwise 200 50 0.2 --out "$D/prune" --silent
$P2 --bfile "$D/merged" --extract "$D/prune.prune.in" --make-bed --out "$D/mp" --silent
awk 'BEGIN{OFS="\t"}{$3=0;$4=0;print}' "$D/mp.fam" > "$D/mp.fam.tmp" && mv "$D/mp.fam.tmp" "$D/mp.fam"

# Population labels (aadr_to_bed.py also puts them in the FID column).
awk '{print $1"\t"$3}' "$A/v66.p1_HO.aadr.patch.PUB.ind" > "$D/id_pop.tsv"
awk 'NR==FNR{p[$1]=$2;next} ($2 in p){print $1"\t"$2}' "$D/id_pop.tsv" "$D/mp.fam" > "$D/ref.keep"
zoom=$(env/bin/python scripts/run_config.py ancestry.zoom_populations | paste -sd"|" -)
awk -v re="^($zoom)$" 'NR==FNR{p[$1]=$2;next} ($2 in p) && p[$2] ~ re {print $1"\t"$2}' "$D/id_pop.tsv" "$D/mp.fam" > "$D/zoom.keep"

pca() {  # $1 = label, $2 = keep file for the reference
    $P2 --bfile "$D/mp" --keep "$2" --maf 0.01 --freq counts --pca 10 allele-wts --out "$R/pca_$1" --silent
    $P2 --bfile "$D/mp" --read-freq "$R/pca_$1.acount" --score "$R/pca_$1.eigenvec.allele" 2 6 header-read \
        no-mean-imputation variance-standardize --score-col-nums 7-16 --out "$R/proj_$1" --silent
}
pca weur "$D/ref.keep"
# Regional zoom: drop reference outliers (> 6 SD on any of PC1-4, smartpca-style), up to 5 rounds, then project.
cp "$D/zoom.keep" "$D/zoom_clean.keep"
for i in 1 2 3 4 5; do
    $P2 --bfile "$D/mp" --keep "$D/zoom_clean.keep" --maf 0.01 --pca 4 --out "$D/zoom_iter" --silent
    n=$(env/bin/python -c "
import pandas as pd
e = pd.read_csv('$D/zoom_iter.eigenvec', sep='\t')
pcs = ['PC1', 'PC2', 'PC3', 'PC4']
z = (e[pcs] - e[pcs].mean()) / e[pcs].std()
out = (z.abs() > 6).any(axis=1)
e.loc[~out, ['#FID', 'IID']].to_csv('$D/zoom_clean.keep', sep='\t', header=False, index=False)
e.loc[out, ['#FID', 'IID']].to_csv('$D/zoom_outliers_$i.txt', sep='\t', header=False, index=False)
print(int(out.sum()))")
    echo "# zoom outlier round $i: removed $n"
    [ "$n" = 0 ] && break
done
pca zoom "$D/zoom_clean.keep"
cat "$D"/zoom_outliers_*.txt > "$R/zoom_outliers.tsv"
cp "$D/id_pop.tsv" "$R/id_pop.tsv"
wc -l "$D/keep.ids" "$D/prune.prune.in" "$D/ref.keep" "$D/zoom.keep" | sed 's/^/# /'
grep SELF "$R/proj_weur.sscore" | cut -f1-6
