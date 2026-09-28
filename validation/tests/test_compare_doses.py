"""Tests for compare_doses.run_gamma axis handling."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compare_doses import build_axes, run_gamma  # noqa: E402


def _meta(dims: list[int]) -> dict:
    """Metadata as load_mhd returns it: dims and spacing in MHD (x, y, z) order."""
    return {"dims": dims, "spacing": [2.0, 2.0, 3.0], "offset": [0.0, 0.0, 0.0]}


def _dose(dims: list[int]) -> np.ndarray:
    """A smooth dose in numpy (z, y, x) order, as load_mhd reshapes it."""
    nx, ny, nz = dims
    z, y, x = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    return (100.0 * np.exp(-((x - nx / 2) ** 2 + (y - ny / 2) ** 2) / 20.0) * (1.0 + z / nz)).astype(np.float32)


@pytest.mark.parametrize("dims", [[12, 10, 6], [6, 10, 12], [8, 8, 5]])
def test_gamma_accepts_non_cubic_grid(dims: list[int]) -> None:
    """Every extent differing (or only z) must not trip pymedphys' axes/shape check."""
    dose = _dose(dims)
    axes = build_axes(_meta(dims))
    _, pass_rate = run_gamma(dose, dose, axes, axes)
    assert pass_rate == pytest.approx(100.0)


def test_gamma_detects_a_shift_along_z() -> None:
    """A dose scaled only along z must fail gamma where it exceeds 3%, so the z axis is really z."""
    dims = [8, 8, 6]
    ref = _dose(dims)
    evl = ref.copy()
    evl[3:, :, :] *= 1.2
    axes = build_axes(_meta(dims))
    _, pass_rate = run_gamma(ref, evl, axes, axes)
    assert pass_rate < 100.0
