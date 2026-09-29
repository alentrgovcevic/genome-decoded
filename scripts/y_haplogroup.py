"""Call the Y haplogroup with 23andMe's yhaplo (ISOGG 2016 tree bundled with it).

Writes data/y.genos.txt (yhaplo sample-major format) and runs yhaplo into results/yhaplo/.
Het and no-call Y genotypes become "." (unobserved).
"""

import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()  # the run directory runs/<sample>/ (see workflow/Snakefile)
GENOME = ROOT / "data" / "genome_23andme.txt"
OUT = ROOT / "results" / "yhaplo"


def main() -> None:
    from run_config import sex

    if sex() != "male":
        OUT.mkdir(parents=True, exist_ok=True)
        for f in [
            *OUT.glob("*.txt"),
            *OUT.glob("*.nwk"),
        ]:  # yhaplo output from a run that took this sample as male
            f.unlink()
        (OUT / "SKIPPED").write_text(f"no Y chromosome: inferred sex {sex()}\n")
        print(f"Y haplogroup skipped (inferred sex: {sex()})")
        return
    g = pd.read_csv(
        GENOME,
        sep="\t",
        comment="#",
        names=["rsid", "chrom", "pos", "gt"],
        dtype={"chrom": str},
    )
    y = g[g["chrom"] == "Y"].sort_values("pos").drop_duplicates("pos")
    calls = [gt[0] if gt[0] == gt[1] and gt[0] in "ACGT" else "." for gt in y["gt"]]
    genos = ROOT / "data" / "y.genos.txt"
    genos.write_text(
        "ID\t" + "\t".join(map(str, y["pos"])) + "\nSELF\t" + "\t".join(calls) + "\n"
    )
    print(f"{len(y)} Y positions, {sum(c != '.' for c in calls)} usable calls")
    OUT.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(ROOT / "env" / "bin" / "yhaplo"),
            "-i",
            str(genos),
            "-o",
            str(OUT),
            "-hp",
            "-ds",
            "-as",
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
