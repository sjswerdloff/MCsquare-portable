"""Behavioural tests for idd_r80.py on synthetic depth-only grids written the way MCsquare writes them."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "idd_r80.py"


def _write(tmp_path: Path, values_by_y: np.ndarray, spacing_y: float, offset_y: float, dims=None) -> Path:
    n = values_by_y.size
    dims = dims or (1, n, 1)
    (tmp_path / "Dose.raw").write_bytes(values_by_y.astype("<f4").tobytes())
    (tmp_path / "Dose.mhd").write_text(
        "ObjectType = Image\nNDims = 3\n"
        f"DimSize = {dims[0]} {dims[1]} {dims[2]}\n"
        f"ElementSpacing = 400.000000 {spacing_y:f} 400.000000\n"
        f"Offset = 0.000000 {offset_y:f} 0.000000\n"
        "ElementType = MET_FLOAT\nElementByteOrderMSB = False\nElementDataFile = Dose.raw\n"
    )
    return tmp_path / "Dose.mhd"


def _run(mhd: Path, entrance: float = 350.0) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), str(mhd), "--entrance-mm", str(entrance)],
                          capture_output=True, text=True, check=False)


def _record(r: subprocess.CompletedProcess) -> dict:
    assert r.returncode == 0, r.stderr
    return json.loads(r.stdout.split("IDD_R80 ", 1)[1])


def _profile(depth: np.ndarray, r80_true: float) -> np.ndarray:
    # rising to a FLAT top of 1.0 over [peak - 1.5, peak] (so the sampled max is exactly 1.0 on any grid of <= 1 mm),
    # then a linear falloff through 0.8 exactly at r80_true, zero beyond.
    peak = r80_true - 2.0
    rise = 0.5 + 0.5 * np.clip(depth / (peak - 1.5), 0.0, 1.0)
    d = np.where(depth <= peak, rise, 1.0 - 0.1 * (depth - peak))
    return np.clip(d, 0.0, None)


@pytest.mark.parametrize("spacing", [1.0, 0.1])
def test_r80_recovers_a_known_crossing_with_the_beam_entering_at_high_y(tmp_path, spacing):
    offset = 280.0
    n = int(round(70.0 / spacing))
    y_centre = offset + (np.arange(n) + 0.5) * spacing
    depth = 350.0 - y_centre  # index 0 is the DEEPEST bin: the file is stored in +y order
    rec = _record(_run(_write(tmp_path, _profile(depth, 40.8), spacing, offset)))
    assert rec["R80"] == pytest.approx(40.8, abs=1e-6)  # linear falloff: interpolation is exact
    assert rec["bin_mm"] == spacing and rec["R80_multiple_crossings"] is False


def test_depth_comes_from_the_header_offset_not_the_index(tmp_path):
    # ONE physical profile sampled on two grids whose origins differ by 1.3 mm: the same R80 only if depth is taken
    # from the header Offset (a reader using the index alone would report them 1.3 mm apart)
    spacing, n = 0.1, 650
    results = []
    for offset in (285.0, 283.7):
        d = tmp_path / f"o{offset}"
        d.mkdir()
        depth = 350.0 - (offset + (np.arange(n) + 0.5) * spacing)
        results.append(_record(_run(_write(d, _profile(depth, 40.8), spacing, offset)))["R80"])
    assert results[0] == pytest.approx(40.8, abs=1e-6) and results[1] == pytest.approx(40.8, abs=1e-6)


def test_non_depth_only_grid_is_refused(tmp_path):
    r = _run(_write(tmp_path, np.ones(4 * 70), 1.0, 280.0, dims=(4, 70, 1)))
    assert r.returncode != 0 and "depth-only" in r.stderr


def test_truncated_raw_is_refused(tmp_path):
    mhd = _write(tmp_path, np.ones(70), 1.0, 280.0)
    (tmp_path / "Dose.raw").write_bytes(b"\0" * (4 * 70 - 1))
    r = _run(mhd)
    assert r.returncode != 0 and "bytes" in r.stderr


def test_no_distal_crossing_reports_none(tmp_path):
    rec = _record(_run(_write(tmp_path, np.linspace(0.1, 1.0, 70)[::-1], 1.0, 280.0)))  # max at the deepest bin
    assert rec["R80"] is None


def test_grid_upstream_of_the_entrance_is_refused(tmp_path):
    r = _run(_write(tmp_path, np.ones(70), 1.0, 290.0))  # ends at y = 360 > 350
    assert r.returncode != 0 and "upstream" in r.stderr
