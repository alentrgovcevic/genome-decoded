"""Build the public demo input: a 1000 Genomes individual's genotypes at the GSA chip's sites, written as a
MyHeritage-format raw export, so the demo exercises exactly the path a real user's file takes.

Usage (from the repo root, after download_refs.sh):
    env/bin/python demo/make_demo.py [IID]          default IID = NA12878 (CEU; GIAB HG001; female)
Reads the ORIGINAL 1000G files (reference/1kg/), not the panel: the panel excludes this person by design
(config exclude_from_reference), so the demo never sees its own genotypes in any reference.
Writes demo/<IID>_GSA_myheritage.csv. Sites the person has no call for (or indels, multi-allelic sites) are "--".
"""

import gzip
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
P2 = str(REPO / "env-x86" / "bin" / "plink2")
BCF = str(REPO / "env" / "bin" / "bcftools")


def main(iid: str = "NA12878") -> None:
    work = REPO / "demo" / f".{iid}"
    work.mkdir(parents=True, exist_ok=True)
    sites = [
        l.rstrip("\n").split("\t")
        for l in gzip.open(REPO / "demo" / "gsa_sites.tsv.gz", "rt")
        if not l.startswith(("#", "rsid"))
    ]
    rng = work / "sites.range"
    rng.write_text("".join(f"{c}\t{p}\t{p}\n" for _, c, p in sites if c not in ("MT",)))
    (work / "keep.txt").write_text(iid + "\n")
    subprocess.run(
        [
            P2,
            "--pfile",
            str(REPO / "reference" / "1kg" / "all_phase3"),
            "vzs",
            "--keep",
            str(work / "keep.txt"),
            "--extract",
            "range",
            str(rng),
            "--snps-only",
            "just-acgt",
            "--max-alleles",
            "2",
            "--export",
            "vcf",
            "bgz",
            "--out",
            str(work / "sites"),
            "--silent",
        ],
        check=True,
    )
    q = subprocess.run(
        [
            BCF,
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
            str(work / "sites.vcf.gz"),
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    calls = {}
    for line in q.splitlines():
        c, p, ref, alt, gt = line.split("\t")
        key = (c, p)
        if key in calls:  # two records at one position: ambiguous, leave as no-call
            calls[key] = "--"
            continue
        al = [ref, alt]
        parts = gt.replace("|", "/").split("/")
        if "." in parts:
            calls[key] = "--"
        else:
            letters = "".join(al[int(x)] for x in parts)
            calls[key] = (
                letters * 2 if len(letters) == 1 else letters
            )  # haploid calls written as two letters
    out = REPO / "demo" / f"{iid}_GSA_myheritage.csv"
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    with open(out, "w") as fh:
        fh.write(
            f"##fileformat=MyHeritage\n##format=MHv1.0\n##chip=GSA\n##timestamp={stamp}\n##reference=build37\n#\n"
            f"# genome-decoded DEMO: 1000 Genomes {iid} genotypes at GSA sites, written in MyHeritage layout.\n"
            f"# Public data (1000 Genomes Project, open consent). Not a MyHeritage product.\n"
        )
        fh.write("RSID,CHROMOSOME,POSITION,RESULT\n")
        for rsid, c, p in sites:
            fh.write(f'"{rsid}","{c}","{p}","{calls.get((c, p), "--")}"\n')
    called = sum(1 for _, c, p in sites if calls.get((c, p), "--") != "--")
    print(
        f"{out.name}: {len(sites):,} sites, {called:,} with a genotype ({called / len(sites):.1%})"
    )


if __name__ == "__main__":
    main(*sys.argv[1:2])
