#!/bin/bash
# Download the public reference data into ./reference (or $REF_DIR). Idempotent: existing files are skipped.
# Measured sizes (2026-09-28): 1000 Genomes pgen.zst 2.4 GB + pvar.zst 0.6 GB (7 GB once decompressed),
# ClinVar 0.45 GB, AADR Human Origins 4.1 GB, Beagle maps/jar 0.1 GB, HIBAG model 18 MB, PharmCAT GRCh38 FASTA ~0.9 GB.
# Measured: ~9 GB of downloads grow to ~37 GB once the workflow builds its panels, plus ~28 GB per analysed person. Nothing personal is sent anywhere: every request
# here is a plain download of a public file.
set -euo pipefail
cd "$(dirname "$0")"
R=${REF_DIR:-reference}
mkdir -p "$R"/{1kg,aadr,beagle,liftover,hibag,pgs}
get() { [ -s "$2" ] && { echo "  have $(basename "$2")"; return; }; echo "  get  $(basename "$2")"; curl -fSL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"; }

echo "== 1000 Genomes phase 3 (GRCh37), plink2 resources page"
page=$(curl -fsSL https://www.cog-genomics.org/plink/2.0/resources)
link() { echo "$page" | grep -oE "href=\"[^\"]*$1[^\"]*\"" | head -1 | sed 's/^href="//; s/"$//; s/&amp;/\&/g'; }
[ -s "$R/1kg/all_phase3.pgen" ] && echo "  have all_phase3.pgen (decompressed)" \
    || get "$(link '/all_phase3.pgen.zst')" "$R/1kg/all_phase3.pgen.zst"   # the panel step decompresses it
get "$(link 'all_phase3_noannot.pvar.zst')" "$R/1kg/all_phase3.pvar.zst"   # noannot, renamed as the page instructs
get "$(link 'phase3_corrected.psam')" "$R/1kg/all_phase3.psam"
get "$(link 'deg2_phase3.king.cutoff.out.id')" "$R/1kg/deg2_phase3.king.cutoff.out.id"

echo "== ClinVar variant_summary (GRCh37 + GRCh38 rows)"
CV=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
if [ ! -s "$R/variant_summary.txt.gz" ]; then
    get "$CV" "$R/variant_summary.txt.gz"
    # the report cites ClinVar's release date, not the download date
    curl -fsSI "$CV" | awk -F': ' 'tolower($1)=="last-modified"{print $2}' | tr -d '\r' > "$R/variant_summary.last_modified"
else echo "  have variant_summary.txt.gz"; fi

echo "== AADR v66 present-day Human Origins (Harvard Dataverse doi:10.7910/DVN/FFIDCW)"
env/bin/python - "$R/aadr" <<'PY'
import hashlib, json, os, subprocess, sys
out = sys.argv[1]
# curl, not urllib: Dataverse answers 403 to Python's default User-Agent
meta = json.loads(subprocess.run(["curl", "-fsSL", "https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/FFIDCW"],
                                 check=True, capture_output=True).stdout)
files = {f["dataFile"]["filename"]: f["dataFile"] for f in meta["data"]["latestVersion"]["files"]}
for name in ["v66.p1_HO.aadr.patch.PUB.geno", "v66.p1_HO.aadr.patch.PUB.snp", "v66.p1_HO.aadr.patch.PUB.ind", "v66.p1_HO.aadr.PUB.anno"]:
    f, dest = files[name], f"{out}/{name}"
    if os.path.exists(dest) and os.path.getsize(dest) == f["filesize"]:
        print("  have", name); continue
    print("  get ", name, f"({f['filesize'] / 1e6:.0f} MB)")
    subprocess.run(["curl", "-fSL", "--retry", "3", "-o", dest, f"https://dataverse.harvard.edu/api/access/datafile/{f['id']}"], check=True)
    h = hashlib.md5()
    with open(dest, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    got = h.hexdigest()
    assert got == f.get("md5"), f"md5 mismatch for {name}: {got}"
print("  AADR licence:", (meta["data"]["latestVersion"].get("license") or {}).get("name", "see the dataset page"),
      "- cite Mallick et al. 2024, Sci Data 11:182")
PY

echo "== Beagle 5.5 bref3 tool + GRCh37 genetic maps"
get https://faculty.washington.edu/browning/beagle/bref3.27Feb25.75f.jar "$R/beagle/bref3.27Feb25.75f.jar"
get https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip "$R/beagle/plink.GRCh37.map.zip"
[ -s "$R/beagle/plink.chr1.GRCh37.map" ] || (cd "$R/beagle" && unzip -oq plink.GRCh37.map.zip)

echo "== UCSC liftover chains"
get https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz "$R/liftover/hg19ToHg38.over.chain.gz"
get https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz "$R/liftover/hg38ToHg19.over.chain.gz"

echo "== HIBAG pre-fit model (official bucket; Illumina GSA v2.0, European, 2-field)"
get https://hibag.s3.amazonaws.com/download/HLARES/InfiniumGlobal-European-HLA4-hg19.RData "$R/hibag/InfiniumGlobal-European-HLA4-hg19.RData"

echo "== PGS Catalog scoring files (GRCh37-harmonised) for every ID in config/config.yaml"
for id in $(env/bin/python -c "import yaml; c=yaml.safe_load(open('config/config.yaml'))['pgs']; print(' '.join(sorted({i for v in c.values() for i in v})))"); do
    [ -s "$R/pgs/$id.json" ] || curl -fsSL "https://www.pgscatalog.org/rest/score/$id" -o "$R/pgs/$id.json"
    url=$(env/bin/python -c "import json; print(json.load(open('$R/pgs/$id.json'))['ftp_harmonized_scoring_files']['GRCh37']['positions'])")
    get "$url" "$R/pgs/${id}_hmPOS_GRCh37.txt.gz"
done

echo "== PharmCAT GRCh38 reference FASTA (PharmCAT's own download, into env/bin)"
[ -s env/bin/reference.fna.bgz ] && echo "  have reference.fna.bgz" \
    || env/bin/python -c "from pcat import utilities; utilities.download_reference_fasta_and_index('env/bin')"

echo "== done. Next: env/bin/python demo/make_demo.py && env/bin/snakemake -s workflow/Snakefile --cores 8"
