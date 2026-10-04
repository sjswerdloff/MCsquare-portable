"""Tests for the part E run lists (field_followup_runs.py) and the parameterised field endpoints.

Run: uv run --no-project --with numpy --with pytest pytest validation/tests/test_field_followup_runs.py
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import apples_seeds_check as sc
import field_followup_runs as fr
import platform_study_metrics as psm
from field_edge_analyse import SP, C, N

PLAN = Path(__file__).resolve().parents[2] / "docs" / "field_100_150_plan.md"


# ------------------------------------------------------------------------------------------------- run lists
def test_seeds_are_the_plans_table() -> None:
    want = {
        ("A-up", 100): range(960041, 960049), ("A-up", 150): range(960051, 960059),
        ("A-port", 100): range(961041, 961049), ("A-port", 150): range(961051, 961059),
        ("B-up", 100): range(962041, 962049), ("B-up", 150): range(962051, 962059),
        ("B-pgcc", 100): range(963041, 963049), ("B-pgcc", 150): range(963051, 963059),
        ("B-picc", 100): range(964041, 964049), ("B-picc", 150): range(964051, 964059),
    }
    for (arm, energy), seeds in want.items():
        assert fr.seeds(arm, energy, "full") == tuple(seeds)
        assert fr.seeds(arm, energy, "smoke") == (seeds[0] + 900,)
    text = PLAN.read_text(encoding="utf-8")
    for arm in fr.ARMS:  # the plan's table row for each arm states the same ranges
        b = fr.ARMS[arm].block
        assert f"| {arm} | {b}041–{b}048 | {b}051–{b}058 | {b}941, {b}951 |" in text


def test_design_constants_match_the_plan() -> None:
    assert fr.ENERGIES == (100, 150)
    assert fr.HISTORIES == {"full": 60_000_000, "smoke": 100_000}
    assert fr.RUNS_PER_CELL == 8 and fr.SIDE_MM == 150.0
    assert len(fr.all_runs()) == 90  # 80 full + 10 smoke
    assert [fr.plan_name(e) for e in fr.ENERGIES] == ["E100_S150", "E150_S150"]
    assert {a: (v.threads, v.host) for a, v in fr.ARMS.items()} == {
        "A-up": (4, "DESKTOP-SR5GKKA"), "A-port": (4, "DESKTOP-SR5GKKA"),
        "B-up": (3, "DESKTOP-5H86O9N"), "B-pgcc": (3, "DESKTOP-5H86O9N"), "B-picc": (3, "DESKTOP-5H86O9N"),
    }
    assert {a for a, v in fr.ARMS.items() if v.portable} == set(sc.PORTABLE_ARMS)
    text = PLAN.read_text(encoding="utf-8")
    for arm in fr.ARMS.values():
        assert re.fullmatch(r"[0-9a-f]{64}", arm.binary_sha256)
        assert f"`{arm.binary_sha256}`" in text
    assert len({a.binary_sha256 for a in fr.ARMS.values()}) == 5
    assert fr.CONTRASTS == (("A-port", "A-up", True), ("B-pgcc", "B-up", True), ("B-picc", "B-up", True),
                            ("B-pgcc", "B-picc", False))


@pytest.mark.parametrize(("arm", "energy", "mode"), [("A-port", 200, "full"), ("C-up", 100, "full"), ("A-up", 100, "pilot")])
def test_anything_outside_the_design_is_refused(arm: str, energy: int, mode: str) -> None:
    with pytest.raises(ValueError, match="not in the design"):
        fr.seeds(arm, energy, mode)


def test_the_committed_lists_pass_the_seed_screen() -> None:
    assert fr.check() == []


def test_a_seed_of_parts_a_and_b_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fr, "_FIRST", {100: 31, 150: 51})  # 96x031-038 are part A and B's field seeds
    errors = fr.check()
    assert any("A-port 100 MeV full seed 961031: repeats part A/B A-port" in e for e in errors)


def test_a_shared_generator_start_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    # 861041 + 1e5 * (c + 1) == 961041 + 1e5 * c on every thread: the same (initstate, stream) one batch call apart
    monkeypatch.setattr(sc, "EARLIER_PORTABLE", frozenset({861041}))
    errors = fr.check()
    assert any(e.startswith("A-port 100 MeV full seed 961041: generator start") and "earlier portable seed 861041" in e
               for e in errors)
    assert not any("A-up" in e for e in errors)  # the upstream arms use another generator


def test_a_seed_outside_its_block_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fr, "_FIRST", {100: 41, 150: 1051})
    assert any("outside the arm's six-digit block" in e for e in fr.check())


def test_command_line(capsys: pytest.CaptureFixture[str]) -> None:
    assert fr.main(["seeds", "B-picc", "150", "smoke"]) == 0
    assert capsys.readouterr().out.strip() == "964951"
    assert fr.main(["histories", "full"]) == 0
    assert capsys.readouterr().out.strip() == "60000000"
    assert fr.main(["binary", "A-up"]) == 0
    assert capsys.readouterr().out.strip() == fr.ARMS["A-up"].binary_sha256
    assert fr.main(["check"]) == 0
    assert capsys.readouterr().out.startswith("OK: 90 runs")
    assert fr.main(["seeds", "A-up", "200", "full"]) == 2
    assert fr.main(["binary", "nobody"]) == 2
    assert fr.main([]) == 2


# ----------------------------------------------------------------------------------------------- field endpoints
def synthetic(peak_row: int = 110, floor: float = 0.0) -> np.ndarray:
    """d[z, y, x]: a 150 mm field (1 inside, exp(-(r - 75) / 10) outside, per axis) times a depth curve with a peak."""
    c = (np.arange(N) - C + 0.5) * SP  # voxel centres, mm from the axis
    lat = np.where(np.abs(c) < 75, 1.0, np.exp(-(np.abs(c) - 75) / 10.0))
    iy = np.arange(N)
    depth = 1.0 + 3.0 * np.exp(-(((iy - peak_row) / 4.0) ** 2))
    depth[:peak_row - 8] = floor  # beyond the distal edge (smaller iy is deeper)
    return (lat[:, None, None] * depth[None, :, None] * lat[None, None, :]).astype(np.float32)


def test_row_of_depth_is_the_inverse_of_the_depth_formula() -> None:
    for depth in (39, 61, 79, 125, 127, 201):
        iy = psm.row_of_depth(depth)
        assert (N - iy - 0.5) * SP == depth
    assert psm.row_of_depth(127) == psm.DEPTH_IY[127] and psm.row_of_depth(201) == psm.DEPTH_IY[201]


@pytest.mark.parametrize("bad", [40, 0, -1, 39.0, True, 297, 299])
def test_row_of_depth_refuses_even_out_of_grid_and_non_integer_depths(bad: object) -> None:
    with pytest.raises(ValueError):
        psm.row_of_depth(bad)  # type: ignore[arg-type]


def test_field_endpoints_values_on_a_known_cube() -> None:
    d = synthetic()
    out = psm.field_endpoints(d, 150.0, (39, 61), (5, 10, 20, 30, 50, 70))
    for depth in (39, 61):
        for off in (5, 10, 20, 30, 50, 70):
            centre_mm = (int((75 + off) / SP) + 0.5) * SP  # the voxel the estimator reads
            assert out[f"lateral_{depth}_{off}"] == pytest.approx(100 * math.exp(-(centre_mm - 75) / 10.0), rel=1e-5)
        assert out[f"cax_{depth}"] == pytest.approx(float(d[C, psm.row_of_depth(depth), C]), rel=1e-6)
    assert out["cax_max"] == pytest.approx(4.0, rel=1e-6)
    assert out["r80_crossings"] == 1 and out["r20_crossings"] == 1
    assert out["r20_mm"] > out["r80_mm"] > (N - 110 - 0.5) * SP  # distal to the peak row
    assert len([k for k in out if k.startswith("lateral_")]) == 12
    assert not any(k.endswith("_invalid") for k in out)


def test_field_endpoints_agree_with_the_200_mev_estimator() -> None:
    rng = np.random.default_rng(7)
    d = (synthetic(peak_row=23, floor=0.0) * rng.uniform(0.9, 1.1, size=(N, N, N))).astype(np.float32)
    old = psm.endpoints(d, 150.0)
    new = psm.field_endpoints(d, 150.0, (127, 201), (5, 10, 20, 30))
    for k in [f"lateral_{dep}_{o}" for dep in (127, 201) for o in (5, 10, 20, 30)] + ["cax_127", "cax_201", "r80_mm", "r20_mm",
                                                                                     "r80_crossings", "r20_crossings"]:
        assert new[k] == old[k]


def test_a_strip_with_no_dose_gives_a_lateral_value_of_exactly_zero() -> None:
    d = synthetic()
    k = int((75 + 70) / SP)
    for a in (C + k, C - 1 - k):
        d[:, :, a] = 0.0
        d[a, :, :] = 0.0
    assert psm.field_endpoints(d, 150.0, (39,), (50, 70))["lateral_39_70"] == 0.0


def test_an_offset_beyond_the_grid_is_refused() -> None:
    with pytest.raises(ValueError, match="beyond the grid"):
        psm.field_endpoints(synthetic(), 150.0, (39,), (80,))


def test_non_positive_cax_dose_is_flagged_and_the_lateral_value_is_not_finite() -> None:
    d = synthetic()
    d[:, psm.row_of_depth(61) - 2:psm.row_of_depth(61) + 3, :] = 0.0
    out = psm.field_endpoints(d, 150.0, (39, 61), (5,))
    assert out["cax_61_invalid"] is True and "cax_39_invalid" not in out
    assert not math.isfinite(out["lateral_61_5"])


def test_wrong_cube_shape_is_refused() -> None:
    with pytest.raises(ValueError, match="expected"):
        psm.field_endpoints(np.zeros((N, N, N - 1), dtype=np.float32), 150.0, (39,), (5,))
