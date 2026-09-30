"""Behavioural tests for platform_study_analyse.py on SYNTHETIC records only (no study data is read).

Run: uv run --no-project --with numpy --with scipy --with pytest pytest validation/test_platform_study_analyse.py
"""

import copy
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import platform_study_analyse as A  # noqa: E402

SHARED = {"plan E200_S150.txt": "p", "cube.mhd": "m", "cube.raw": "r", "BDL": "b", "HU_Density": "d", "HU_Material": "t"}


def rec(plat: str, seed: int, shift: float = 0.0, commit: str | None = None) -> dict:
    """One synthetic record; values are platform-independent unless shifted, with small seed-dependent noise."""
    noise = ((seed * 7919) % 97 - 48) / 48 * 0.01
    metrics = {ep: 3.0 + shift + noise for ep in A.LATERAL}
    metrics |= {ep: 1.0e-3 * (1 + noise / 10) for ep in A.CAX}
    metrics |= {"r80_mm": 250.0 + noise, "r20_mm": 255.0 + noise, "r80_crossings": 1, "r20_crossings": 1}
    if commit is None:
        commit = A.RERUN_COMMIT if (plat, seed) in A.RERUN_SEEDS else A.STUDY_COMMIT
    return {"study": "pe1", "platform": plat, "seed": seed, "host": "h", "commit": commit,
            "sha256": dict(SHARED, binary=plat), "materials": {"combined_sha256": "M"}, "metrics": metrics}


def wave(w: int, plats=("windows", "linux", "macos")) -> list[dict]:
    return [rec(p, s) for p in plats for s in A.WAVE_SEEDS[w][p]]


def run(tmp_path: Path, recs: list[dict], look: int, verdicts: Path | None = None) -> subprocess.CompletedProcess:
    data = tmp_path / f"records_{look}.jsonl"
    data.write_text("".join(json.dumps(r) + "\n" for r in recs))
    v = verdicts or tmp_path / "look1.json"
    return subprocess.run([sys.executable, str(HERE / "platform_study_analyse.py"), str(data), str(look), str(v)],
                          capture_output=True, text=True)


def refused(p: subprocess.CompletedProcess, why: str) -> bool:
    return p.returncode != 0 and "REFUSING TO ANALYSE" in p.stderr and why in p.stderr


def test_complete_wave1_runs_and_writes_verdicts(tmp_path):
    p = run(tmp_path, wave(1), 1)
    assert p.returncode == 0, p.stderr
    v = json.loads((tmp_path / "look1.json").read_text())
    assert set(v) == {"windows", "macos"}
    assert "VERDICT windows vs linux (look 1)" in p.stdout


def test_look1_refuses_to_overwrite_its_verdicts(tmp_path):
    assert run(tmp_path, wave(1), 1).returncode == 0
    assert refused(run(tmp_path, wave(1), 1), "already analysed")


@pytest.mark.parametrize("look", [0, 3, -1])
def test_invalid_look(tmp_path, look):
    assert refused(run(tmp_path, wave(1), look), "look must be 1 or 2")


def test_no_data(tmp_path):
    assert refused(run(tmp_path, [], 1), "no records")


def test_missing_seed(tmp_path):
    assert refused(run(tmp_path, wave(1)[1:], 1), "missing")


def test_unexpected_seed(tmp_path):
    assert refused(run(tmp_path, wave(1) + [rec("linux", 2999)], 1), "unexpected")


def test_duplicate_seed(tmp_path):
    recs = wave(1)
    assert refused(run(tmp_path, recs + [copy.deepcopy(recs[0])], 1), "duplicate")


def test_missing_platform(tmp_path):
    assert refused(run(tmp_path, wave(1, ("windows", "linux")), 1), "missing")


def test_foreign_study(tmp_path):
    recs = wave(1)
    recs[5]["study"] = "other"
    assert refused(run(tmp_path, recs, 1), "study")


def test_foreign_commit(tmp_path):
    recs = wave(1)
    recs[3]["commit"] = "0" * 40
    assert refused(run(tmp_path, recs, 1), "commit")


def test_rerun_commit_allowed_only_for_the_recorded_seeds(tmp_path):
    recs = wave(1)
    assert any(r["commit"] == A.RERUN_COMMIT for r in recs)  # control: the two re-runs carry it and pass
    assert run(tmp_path, recs, 1).returncode == 0
    recs2 = wave(1)
    next(r for r in recs2 if (r["platform"], r["seed"]) == ("macos", 3001))["commit"] = A.RERUN_COMMIT
    other = tmp_path / "b"
    other.mkdir()
    assert refused(run(other, recs2, 1), "commit")


def test_input_differs_within_a_platform(tmp_path):
    recs = wave(1)
    next(r for r in recs if r["platform"] == "linux")["sha256"]["BDL"] = "other"
    assert refused(run(tmp_path, recs, 1), "differs between seeds")


def test_cube_raw_differs_across_platforms(tmp_path):
    recs = wave(1)
    for r in recs:
        if r["platform"] == "macos":
            r["sha256"]["cube.raw"] = "other"
    assert refused(run(tmp_path, recs, 1), "cube.raw differs between platforms")


def test_text_input_line_endings_across_platforms_are_reported_not_refused(tmp_path):
    recs = wave(1)
    for r in recs:
        if r["platform"] == "windows":
            r["sha256"]["BDL"] = "crlf"
    p = run(tmp_path, recs, 1)
    assert p.returncode == 0 and "note: text input BDL differs" in p.stdout


def test_nonfinite_endpoint_cannot_pass(tmp_path):
    recs = wave(1)
    next(r for r in recs if r["platform"] == "windows")["metrics"]["lateral_127_10"] = math.nan
    p = run(tmp_path, recs, 1)
    assert p.returncode == 0
    assert "lateral_127_10" in p.stdout and "CANNOT PASS: non-finite" in p.stdout
    assert "EQUIVALENT" not in json.loads((tmp_path / "look1.json").read_text())["windows"]


def test_look2_needs_look1_verdicts(tmp_path):
    assert refused(run(tmp_path, wave(1) + wave(2), 2), "needs look 1")


def test_look2_needs_frozen_wave2_commit(tmp_path):
    (tmp_path / "look1.json").write_text(json.dumps({"windows": "CONTINUE to wave 2", "macos": "CONTINUE to wave 2"}))
    assert refused(run(tmp_path, wave(1) + wave(2), 2), "WAVE2_COMMIT not frozen")


def test_stopped_comparison_is_frozen_at_look2(tmp_path, monkeypatch):
    """windows stopped at look 1; look 2 needs only linux+macos wave 2 and never recomputes windows."""
    (tmp_path / "look1.json").write_text(json.dumps({"windows": "EQUIVALENT", "macos": "CONTINUE to wave 2"}))
    monkeypatch.setattr(A, "WAVE2_COMMIT", "f" * 40)
    recs = wave(1) + [dict(r, commit="f" * 40) for r in wave(2, ("linux", "macos"))]
    data = tmp_path / "r.jsonl"
    data.write_text("".join(json.dumps(r) + "\n" for r in recs))
    with pytest.raises(SystemExit) as e:  # the in-process run exits 0 at the end only if nothing refused
        A.main(str(data), 2, str(tmp_path / "look1.json"))
        raise SystemExit(0)
    assert e.value.code == 0
    # windows wave-2 records present -> refused as unexpected
    recs_w = recs + [dict(r, commit="f" * 40) for r in wave(2, ("windows",))]
    data.write_text("".join(json.dumps(r) + "\n" for r in recs_w))
    with pytest.raises(SystemExit, match="unexpected"):
        A.main(str(data), 2, str(tmp_path / "look1.json"))


def test_no_comparison_continued_means_no_look2(tmp_path):
    (tmp_path / "look1.json").write_text(json.dumps({"windows": "EQUIVALENT", "macos": "NON-EQUIVALENT on cax_127"}))
    assert refused(run(tmp_path, wave(1), 2), "there is no look 2")
