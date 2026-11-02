"""Per-run settings written by the workflow into runs/<sample>/run_config.json.

Python:  from run_config import cfg; cfg()["pgs"]
Shell:   env/bin/python scripts/run_config.py <key.subkey>   prints a list one item per line, or a scalar
"""

import json
import subprocess
import sys
from pathlib import Path


def pipeline_commit(repo: Path | None = None) -> str:
    """genome-decoded's git commit (short), '+uncommitted' if tracked files differ from it; 'unknown' without git."""
    repo = repo or Path(__file__).resolve().parent.parent
    git = ["git", "-C", str(repo)]
    try:
        sha = subprocess.run(
            git + ["rev-parse", "--short=12", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            git + ["status", "--porcelain", "--untracked-files=no"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return sha + ("+uncommitted" if dirty else "")


def cfg() -> dict:
    fp = Path.cwd() / "run_config.json"
    if not fp.exists():
        raise SystemExit(
            f"{fp} missing: run steps through the workflow (it writes run_config.json)"
        )
    return json.loads(fp.read_text())


def get(key: str):
    v = cfg()
    for part in key.split("."):
        v = v[part]
    return v


def sex() -> str:
    """'male' / 'female' from plink's inferred sex in data/genome.fam (1/2), 'unknown' otherwise."""
    fam = Path.cwd() / "data" / "genome.fam"
    code = fam.read_text().split()[4] if fam.exists() else "0"
    return {"1": "male", "2": "female"}.get(code, "unknown")


if __name__ == "__main__":
    v = get(sys.argv[1])
    if isinstance(
        v, dict
    ):  # dict of lists (e.g. ancestry.regions): all items, flattened
        v = [x for xs in v.values() for x in xs]
    print("\n".join(map(str, v)) if isinstance(v, list) else v)
