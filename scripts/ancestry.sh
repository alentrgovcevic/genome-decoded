#!/bin/bash
# Ancestry against 1000 Genomes phase 3 (GRCh37, plink2 resources page):
#   intersect by chr:pos -> merge -> LD prune -> PCA (global + EUR) projected with allele weights
#   -> supervised ADMIXTURE on superpopulations (AFR, EUR, EAS, SAS; AMR left out as itself admixed).
# Needs env-x86 (plink2, admixture), env (plink 1.9), reference/panel/ built (workflow rule `panel`). Run from the run directory.
set -euo pipefail

P1=env/bin/plink
P2=env-x86/bin/plink2
ADM=env-x86/bin/admixture
KG=reference/panel
D=data/anc
R=results/ancestry
mkdir -p "$D" "$R"

[ -f "$KG/all_phase3.pgen" ] || $P2 --zst-decompress "$KG/all_phase3.pgen.zst" > "$KG/all_phase3.pgen"

# 1. Self: autosomal biallelic ACGT SNVs, IDs chr:pos, duplicate positions dropped.
# No-call sites carry "0 0" alleles, which plink2 rejects, so drop them with plink 1.9 first.
$P1 --bfile data/genome --autosome --geno 0 --make-bed --out "$D/self0" --silent
$P2 --bfile "$D/self0" --snps-only just-acgt --set-all-var-ids @:# --rm-dup exclude-all \
    --make-bed --out "$D/self" --silent
awk '{print $1"\t"$4"\t"$4"\t"$2}' "$D/self.bim" > "$D/self.range"

# 2. Reference: same positions, unrelated, biallelic ACGT SNVs, MAF >= 1%, no strand-ambiguous pairs.
$P2 --pfile "$KG/all_phase3" vzs --extract bed1 "$D/self.range" --autosome --snps-only just-acgt \
    --max-alleles 2 --remove "$KG/deg2_phase3.king.cutoff.out.id" --set-all-var-ids @:# --rm-dup exclude-all \
    --maf 0.01 --make-bed --out "$D/kg_all" --silent
awk '($5$6=="AT"||$5$6=="TA"||$5$6=="CG"||$5$6=="GC"){print $2}' "$D/kg_all.bim" > "$D/ambiguous.txt"
$P2 --bfile "$D/kg_all" --exclude "$D/ambiguous.txt" --make-bed --out "$D/kg" --silent

# 3. Merge (plink 1.9). Allele conflicts land in .missnp; drop them from both and retry once.
if ! $P1 --bfile "$D/kg" --bmerge "$D/self" --make-bed --out "$D/merged" --silent; then
    $P1 --bfile "$D/kg" --exclude "$D/merged-merge.missnp" --make-bed --out "$D/kg_f" --silent
    $P1 --bfile "$D/self" --exclude "$D/merged-merge.missnp" --make-bed --out "$D/self_f" --silent
    $P1 --bfile "$D/kg_f" --bmerge "$D/self_f" --make-bed --out "$D/merged" --silent
fi
# Keep only sites present in both.
cut -f2 "$D/self.bim" | sort > "$D/self.ids"; cut -f2 "$D/kg.bim" | sort > "$D/kg.ids"
comm -12 "$D/self.ids" "$D/kg.ids" > "$D/shared.ids"

# 4. LD prune on the reference only.
$P2 --bfile "$D/kg" --extract "$D/shared.ids" --indep-pairwise 200 50 0.2 --out "$D/prune" --silent
$P2 --bfile "$D/merged" --extract "$D/prune.prune.in" --make-bed --out "$D/mp" --silent
# Relatives are already removed; clear 1000G parent IDs so every sample counts as a founder.
awk 'BEGIN{OFS="\t"}{$3=0;$4=0;print}' "$D/mp.fam" > "$D/mp.fam.tmp" && mv "$D/mp.fam.tmp" "$D/mp.fam"
awk '$1!="SELF"{print $1"\t"$2}' "$D/mp.fam" > "$D/ref.keep"

# Population labels from the psam (columns: #IID PAT MAT SEX SuperPop Population).
awk 'NR>1{print $1"\t"$5"\t"$6}' "$KG/all_phase3.psam" > "$D/pops.tsv"
awk 'NR==FNR{sp[$1]=$2;next} sp[$2]=="EUR"{print $1"\t"$2}' "$D/pops.tsv" "$D/mp.fam" > "$D/eur.keep"

# 5. PCA on reference, then project everyone (reference included) with the same weights so scales match.
pca() {  # $1 = label, $2 = keep file for the reference
    # --maf re-applied within the kept reference: sites monomorphic there break variance-standardize.
    $P2 --bfile "$D/mp" --keep "$2" --maf 0.01 --freq counts --pca 10 allele-wts --out "$R/pca_$1" --silent
    $P2 --bfile "$D/mp" --read-freq "$R/pca_$1.acount" --score "$R/pca_$1.eigenvec.allele" 2 6 header-read \
        no-mean-imputation variance-standardize --score-col-nums 7-16 --out "$R/proj_$1" --silent
}
pca global "$D/ref.keep"
pca eur "$D/eur.keep"

# 6. Supervised ADMIXTURE: 20 random unrelated people from each of the 5 populations per superpopulation
#    (AMR excluded as itself admixed; ACB and ASW excluded for the same reason) + SELF. Seeded, so reruns match.
#    (Before 2026-09-28 this took the first 100 by ID per superpopulation: AFR was 93% ACB.)
env/bin/python - "$D/pops.tsv" "$D/mp.fam" "$D/adm.keep" <<'PY'
import sys
import pandas as pd
pops = pd.read_csv(sys.argv[1], sep="\t", names=["IID", "super", "pop"])
fam = pd.read_csv(sys.argv[2], sep=r"\s+", header=None, usecols=[0, 1], names=["FID", "IID"])
m = fam.merge(pops, on="IID")
m = m[(m["super"] != "AMR") & ~m["pop"].isin(["ACB", "ASW"])]
pick = m.groupby("pop", group_keys=False).apply(lambda g: g.sample(min(20, len(g)), random_state=20260928))
out = pd.concat([pick[["FID", "IID"]], pd.DataFrame({"FID": ["SELF"], "IID": ["SELF"]})])
out.to_csv(sys.argv[3], sep="\t", header=False, index=False)
print("# ADMIXTURE references per superpopulation:", pick.groupby("super").size().to_dict())
PY
$P1 --bfile "$D/mp" --keep "$D/adm.keep" --make-bed --out "$D/adm" --silent
awk 'NR==FNR{sp[$1]=$2;next} {print ($2 in sp) ? sp[$2] : "-"}' "$D/pops.tsv" "$D/adm.fam" > "$D/adm.pop"
K=$(grep -v '^-$' "$D/adm.pop" | sort -u | wc -l | tr -d ' ')
(cd "$D" && ../../$ADM -j8 -s 20260928 -B20 --supervised adm.bed "$K" > admixture.log)
paste <(cut -d' ' -f2 "$D/adm.fam") "$D/adm.$K.Q_se" > "$R/admixture_Q_se.tsv"
# ADMIXTURE orders ancestral columns by first appearance in the .pop file.
awk '$1!="-" && !seen[$1]++{print $1}' "$D/adm.pop" > "$R/admixture_columns.txt"
paste <(cut -d' ' -f2 "$D/adm.fam") "$D/adm.$K.Q" > "$R/admixture_Q.tsv"
cp "$D/pops.tsv" "$R/pops.tsv"

wc -l "$D/self.bim" "$D/kg.bim" "$D/shared.ids" "$D/prune.prune.in" | sed 's/^/# /'
grep -m1 SELF "$R/proj_global.sscore"; grep SELF "$R/admixture_Q.tsv"; cat "$R/admixture_columns.txt"
