"""Tests for compare_doses.run_gamma axis handling."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compare_doses import build_axes, run_gamma


def _meta(dims: list[int], spacing: list[float] | None = None) -> dict:
    """Metadata as load_mhd returns it: dims and spacing in MHD (x, y, z) order."""
    return {"dims": dims, "spacing": spacing or [2.0, 2.0, 3.0], "offset": [0.0, 0.0, 0.0]}


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


# Orientation on a CUBIC grid, where a wrong axis order cannot trip pymedphys' shape check and
# would instead pair z's spacing with x's array dimension, returning a wrong pass rate silently.
# Spacing is anisotropic: x = y = 2 mm, z = 5 mm. A one-voxel shift is therefore 5 mm along z
# (beyond the 3 mm DTA, so it must fail) and 2 mm along x (within DTA, so it must pass). With the
# spacings swapped, both verdicts flip, so the pair pins the orientation in both directions.
# (Design: cora-2f1e43dc, review of #7.)
_CUBE = [8, 8, 8]
_ANISO = [2.0, 2.0, 5.0]


def _shift_one_voxel(dose: np.ndarray, axis: int) -> np.ndarray:
    """Shift by one voxel along a numpy axis, repeating the first slice instead of wrapping."""
    out = np.roll(dose, 1, axis=axis)
    first = [slice(None)] * dose.ndim
    first[axis] = 0
    out[tuple(first)] = dose[tuple(first)]
    return out


def _ramp(axis: int) -> np.ndarray:
    """30..100 in steps of 10 (% of max) along one numpy axis, flat along the others.

    Each step is far beyond the 3% dose criterion, so a one-voxel shift can only pass on
    distance: the verdict depends on the spacing gamma assigns to that axis and nothing else.
    """
    shape = [1, 1, 1]
    shape[axis] = _CUBE[axis]
    profile = (30.0 + 10.0 * np.arange(_CUBE[axis])).reshape(shape)
    return np.broadcast_to(profile, tuple(_CUBE)).astype(np.float32).copy()


def test_gamma_fails_a_one_voxel_shift_along_z_of_5mm() -> None:
    ref = _ramp(axis=0)
    axes = build_axes(_meta(_CUBE, _ANISO))
    _, pass_rate = run_gamma(ref, _shift_one_voxel(ref, axis=0), axes, axes)
    assert pass_rate < 50.0  # only the repeated first slice matches exactly


def test_gamma_passes_a_one_voxel_shift_along_x_of_2mm() -> None:
    ref = _ramp(axis=2)
    axes = build_axes(_meta(_CUBE, _ANISO))
    _, pass_rate = run_gamma(ref, _shift_one_voxel(ref, axis=2), axes, axes)
    # 7 of 8 columns: the reference's 100% column has no counterpart once shifted (eval tops
    # out at 90%). With the axis order swapped, x gets 5 mm and only ~1 column in 8 passes.
    assert pass_rate > 80.0
