import json
import subprocess
import sys
from pathlib import Path

import run_config

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_config.py"
CFG = {
    "pgs": {"both": ["PGS1"], "male": ["PGSm"], "female": ["PGSf"]},
    "ancestry": {"regions": {"North": ["A", "B"], "South": ["C"]}},
    "threads": 4,
}


def setup_run(tmp_path, monkeypatch, sex_code):
    (tmp_path / "run_config.json").write_text(json.dumps(CFG))
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "genome.fam").write_text(f"SELF SELF 0 0 {sex_code} -9\n")
    monkeypatch.chdir(tmp_path)


def test_get_and_sex(tmp_path, monkeypatch):
    setup_run(tmp_path, monkeypatch, 2)
    assert run_config.get("pgs.female") == ["PGSf"]
    assert run_config.sex() == "female"


def test_pgs_scores_follow_sex(tmp_path, monkeypatch):
    import pgs

    setup_run(tmp_path, monkeypatch, 1)
    assert pgs.scores() == ["PGS1", "PGSm"]
    (tmp_path / "data" / "genome.fam").write_text("SELF SELF 0 0 0 -9\n")
    assert pgs.scores() == ["PGS1"]  # unknown sex: no sex-specific score


def test_pipeline_commit(tmp_path):
    sha = subprocess.run(
        ["git", "-C", str(SCRIPT.parents[1]), "rev-parse", "--short=12", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    assert run_config.pipeline_commit().split("+")[0] == sha
    assert (
        run_config.pipeline_commit(tmp_path) == "unknown"
    )  # not a git checkout (e.g. a zip download)


def test_cli_flattens(tmp_path, monkeypatch):
    setup_run(tmp_path, monkeypatch, 1)
    out = subprocess.run(
        [sys.executable, str(SCRIPT), "ancestry.regions"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    assert out.split() == ["A", "B", "C"]
    assert (
        subprocess.run(
            [sys.executable, str(SCRIPT), "threads"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        == "4"
    )
