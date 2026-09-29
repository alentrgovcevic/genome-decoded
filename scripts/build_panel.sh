#!/bin/bash
# Reference-level: derive the 1000 Genomes panel that every step uses, from the downloaded files.
#   reference/1kg/all_phase3.{pgen(.zst),pvar.zst,psam} + deg2_phase3.king.cutoff.out.id  (download_refs.sh)
#   reference/panel/exclude.txt  (one IID per line; written by the workflow from config exclude_from_reference)
# Outputs in reference/panel/:
#   all_phase3.{pgen,pvar.zst,psam}   1000G minus excluded samples (e.g. the demo individual, so its own
#                                     genotypes can't leak into imputation, ancestry or score references)
#   deg2_phase3.king.cutoff.out.id    relatives list (copied)
#   eur_unrelated.keep                EUR in the panel minus relatives (the PGS / frequency reference)
#   eur_seq.{pgen,pvar,psam}          those Europeans' sequenced autosomal SNVs, IDs chr:pos:REF:ALT
# Run from the reference directory's parent (the run directory symlinks reference/). Idempotent.
set -euo pipefail
P2=${P2:-env-x86/bin/plink2}
K=reference/1kg
O=reference/panel
mkdir -p "$O"
touch "$O/exclude.txt"
[ -f "$K/all_phase3.pgen" ] || $P2 --zst-decompress "$K/all_phase3.pgen.zst" > "$K/all_phase3.pgen"

if [ ! -f "$O/all_phase3.pgen" ] || ! cmp -s "$O/exclude.txt" "$O/.built_with_exclude.txt"; then
    if [ -s "$O/exclude.txt" ]; then
        $P2 --pfile "$K/all_phase3" vzs --remove "$O/exclude.txt" --make-pgen vzs --out "$O/all_phase3" --silent
    else
        ln -sf "$(cd "$K" && pwd)/all_phase3.pgen" "$O/all_phase3.pgen"
        ln -sf "$(cd "$K" && pwd)/all_phase3.pvar.zst" "$O/all_phase3.pvar.zst"
        cp "$K/all_phase3.psam" "$O/all_phase3.psam"
    fi
    cp "$O/exclude.txt" "$O/.built_with_exclude.txt"
    rm -f "$O/eur_seq.pgen"
fi
cp "$K/deg2_phase3.king.cutoff.out.id" "$O/deg2_phase3.king.cutoff.out.id"
awk -F'\t' 'NR==FNR{rel[$1]=1;next} FNR>1 && $5=="EUR" && !($1 in rel){print $1}' \
    "$O/deg2_phase3.king.cutoff.out.id" "$O/all_phase3.psam" > "$O/eur_unrelated.keep"
if [ ! -f "$O/eur_seq.pgen" ]; then
    $P2 --pfile "$O/all_phase3" vzs --keep "$O/eur_unrelated.keep" --autosome --snps-only just-acgt --max-alleles 2 \
        --set-all-var-ids '@:#:$r:$a' --rm-dup exclude-all --make-pgen --out "$O/eur_seq" --silent
fi
echo "panel: $(($(wc -l < "$O/all_phase3.psam") - 1)) samples ($(wc -l < "$O/exclude.txt" | tr -d ' ') excluded); \
EUR unrelated: $(wc -l < "$O/eur_unrelated.keep" | tr -d ' ')"
