"""Tests for platform_study_record.py on SYNTHETIC inputs: study identity, hashes, and endpoint failure.

Run: uv run --no-project --with numpy --with scipy --with pytest pytest validation/test_platform_study_record.py
"""

import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import platform_study_record as R

METRICS = {"cax_127": 1.5, "r80_mm": 100.25, "r80_crossings": 1}


@pytest.fixture
def run_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A synthetic run directory (cwd) holding every file the record hashes, plus outdir and retain dirs."""
    files = {"cfg.txt": b"Num_Threads 4\n", "E200_S150.txt": b"plan\n", "cube.mhd": b"mhd\n", "cube.raw": b"raw",
             "BDL/BDL_default_DN_RangeShifter.txt": b"bdl\n", "Scanners/default/HU_Density_Conversion.txt": b"d\n",
             "Scanners/default/HU_Material_Conversion.txt": b"m\n", "Materials/A/x.dat": b"x\n",
             "out/Dose.raw": b"\0\0\0\0", "out/Dose.mhd": b"dose\n", "MCsquare": b"binary"}
    for rel, data in files.items():
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_bytes(data)
    (tmp_path / "retain").mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(R, "load", lambda outdir: object())
    monkeypatch.setattr(R, "endpoints", lambda d, side: dict(METRICS))
    monkeypatch.delenv("STUDY_ID", raising=False)
    return tmp_path


def write_record() -> dict:
    R.main("retain", "out", "A-up", "960001", "MCsquare", "icl 2021", "c0ffee")
    return json.loads(Path("retain/record.json").read_text())


def test_default_study_is_pe1_and_existing_fields_are_unchanged(run_dir):
    rec = write_record()
    assert rec["study"] == "pe1"
    assert rec["platform"] == "A-up" and rec["seed"] == 960001 and rec["commit"] == "c0ffee"
    assert rec["metrics"] == METRICS
    assert rec["sha256"]["config"] == hashlib.sha256(b"Num_Threads 4\n").hexdigest()


def test_study_id_env_is_recorded(run_dir, monkeypatch):
    monkeypatch.setenv("STUDY_ID", "apples-a")
    assert write_record()["study"] == "apples-a"


@pytest.mark.parametrize("bad", ["", " ", "a b", "../x", "-x", "a/b"])
def test_malformed_study_id_is_refused_before_a_record_is_written(run_dir, monkeypatch, bad):
    monkeypatch.setenv("STUDY_ID", bad)
    with pytest.raises(SystemExit):
        R.main("retain", "out", "A-up", "960001", "MCsquare", "icl", "c0ffee")
    assert not Path("retain/record.json").exists()


def test_cfg_and_metrics_hashes_are_recorded_and_track_the_content(run_dir, monkeypatch):
    rec = write_record()
    assert rec["cfg_sha256"] == hashlib.sha256(b"Num_Threads 4\n").hexdigest()
    assert rec["metrics_sha256"] == hashlib.sha256(json.dumps(METRICS, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    assert rec["endpoint_status"] == "ok"
    Path("retain/record.json").unlink()
    monkeypatch.setattr(R, "endpoints", lambda d, side: {**METRICS, "cax_127": 1.6})
    Path("cfg.txt").write_bytes(b"Num_Threads 3\n")
    rec2 = write_record()
    assert rec2["cfg_sha256"] != rec["cfg_sha256"]
    assert rec2["metrics_sha256"] != rec["metrics_sha256"]


def test_endpoint_failure_writes_a_record_that_is_not_ok(run_dir, monkeypatch, capsys):
    def boom(d, side):
        raise ValueError("bad dose grid")

    monkeypatch.setattr(R, "endpoints", boom)
    rec = write_record()
    assert rec["endpoint_status"] != "ok"
    assert rec["endpoint_status"] == "error" and "ValueError: bad dose grid" in rec["endpoint_error"]
    assert rec["metrics"] is None and rec["metrics_sha256"] is None
    assert rec["sha256"]["Dose.raw"] == hashlib.sha256(b"\0\0\0\0").hexdigest()  # provenance still recorded
    assert "STUDY_RESULT" in capsys.readouterr().out


def test_endpoint_load_failure_is_also_recorded(run_dir, monkeypatch):
    def no_load(outdir):
        raise FileNotFoundError(outdir)

    monkeypatch.setattr(R, "load", no_load)
    assert write_record()["endpoint_status"] == "error"


def test_existing_record_is_never_overwritten(run_dir):
    write_record()
    with pytest.raises(SystemExit):
        write_record()
