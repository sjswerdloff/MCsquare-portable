"""report_run_times.py: population selection by the commit the analysis read each cell at."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import report_run_times as rrt

ACQ, RERUN = "a" * 40, "b" * 40


def _run(root: Path, name: str, arm: str, energy: int, seed: int, commit: str, wall: int) -> None:
    d = root / name
    d.mkdir(parents=True)
    rec = {"arm": arm, "host": "H", "cpu": "C", "threads": 3, "case": "E", "energy_mev": energy, "seed": seed,
           "simulated": 60_000_000, "wall_s": wall, "commit": commit}
    (d / "run.json").write_text(json.dumps(rec), encoding="utf-8")


def _doc(tmp_path: Path) -> Path:
    doc = {"acquisition_commit": ACQ, "rerun": {"cells": ["X 150 MeV"], "commit": RERUN},
           "runs": [{"arm": "X", "energy": 100, "seed": 1}, {"arm": "X", "energy": 150, "seed": 2}]}
    path = tmp_path / "doc.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path


def test_a_rerun_cell_is_taken_at_the_rerun_commit_and_the_interrupted_run_is_ignored(tmp_path: Path) -> None:
    root = tmp_path / "runs"
    _run(root, "orig100", "X", 100, 1, ACQ, 60)
    _run(root, "orig150", "X", 150, 2, ACQ, 999)  # the interrupted cell's run: not read
    _run(root, "rerun150", "X", 150, 2, RERUN, 120)
    out = rrt.collect([("e", root)], {"e": _doc(tmp_path)})
    assert out["n_runs"] == 2
    assert {g["energy_mev"]: g["median_s"] for g in out["groups"]} == {100: 10.0, 150: 20.0}


def test_a_population_run_missing_at_its_commit_is_refused(tmp_path: Path) -> None:
    root = tmp_path / "runs"
    _run(root, "orig100", "X", 100, 1, ACQ, 60)
    _run(root, "orig150", "X", 150, 2, ACQ, 120)  # only the acquisition's copy of a rerun cell
    with pytest.raises(rrt.RunTimesError, match="not found exactly once"):
        rrt.collect([("e", root)], {"e": _doc(tmp_path)})
