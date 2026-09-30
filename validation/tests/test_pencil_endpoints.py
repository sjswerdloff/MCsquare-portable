"""Tests for pencil_endpoints: layout mapping, slab selection, and each endpoint on known synthetic doses."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pencil_endpoints import (
    _integrated_gaussian,
    canonical_from_mcsquare,
    canonical_from_topas,
    endpoints,
    lateral_centres,
    r80,
    read_mhd,
    read_topas_bin,
    slab_indices,
)


def test_mcsquare_layout_on_asymmetric_phantom(tmp_path):
    """A marked voxel at MCsquare (x, y, z) lands at canonical (depth from the y = NY face, z, x)."""
    nx, ny, nz = 6, 10, 8
    zyx = np.zeros((nz, ny, nx), dtype="<f4")
    zyx[5, 9, 1] = 1.0  # z=5, y=9 (the entrance-face voxel), x=1
    zyx[2, 0, 4] = 2.0  # z=2, y=0 (the far face), x=4
    (tmp_path / "d.raw").write_bytes(zyx.tobytes())
    (tmp_path / "d.mhd").write_text(
        f"NDims = 3\nDimSize = {nx} {ny} {nz}\nElementSpacing = 1 1 1\nElementType = MET_FLOAT\n"
        "ElementByteOrderMSB = False\nElementDataFile = d.raw\n"
    )
    arr, dims, _ = read_mhd(tmp_path / "d.mhd")
    d = canonical_from_mcsquare(arr)
    assert dims == [nx, ny, nz]
    assert d.shape == (ny, nz, nx)
    assert d[0, 5, 1] == 1.0  # entrance voxel is depth bin 0
    assert d[ny - 1, 2, 4] == 2.0  # far face is the last depth bin
    assert d.sum() == 3.0


def test_read_mhd_rejects_size_mismatch(tmp_path):
    (tmp_path / "d.raw").write_bytes(np.zeros(5, dtype="<f4").tobytes())
    (tmp_path / "d.mhd").write_text(
        "DimSize = 2 2 2\nElementSpacing = 1 1 1\nElementType = MET_FLOAT\nElementDataFile = d.raw\n"
    )
    with pytest.raises(ValueError):
        read_mhd(tmp_path / "d.mhd")


def test_lateral_centres_put_axis_on_corner():
    c = lateral_centres(400)
    assert c[199] == -0.5 and c[200] == 0.5
    with pytest.raises(ValueError):
        lateral_centres(401)


@pytest.mark.parametrize("depth, expected", [(100, [97, 98, 99, 100, 101, 102]), (200, list(range(197, 203)))])
def test_slab_bins_centred_on_depth(depth, expected):
    idx = slab_indices(depth)
    assert list(idx) == expected
    assert (idx + 0.5).mean() == depth  # midpoint exactly d
    assert (idx[0] + 0.5, idx[-1] + 0.5) == (depth - 2.5, depth + 2.5)


def test_slab_rejects_fractional_depth():
    with pytest.raises(ValueError):
        slab_indices(100.5)


def test_r80_interpolates_first_distal_crossing():
    idd = np.array([1.0, 2.0, 5.0, 10.0, 9.0, 7.0, 3.0, 0.0])  # max 10 at bin 3, level 8
    depth, multi = r80(idd)
    # crossing between bins 4 (9, centre 4.5) and 5 (7, centre 5.5): 4.5 + (9 - 8)/(9 - 7)
    assert depth == pytest.approx(5.0)
    assert multi is False


def test_r80_absent_and_multiple():
    assert r80(np.array([1.0, 2.0, 3.0]))[0] is None  # no crossing distal to the max
    depth, multi = r80(np.array([10.0, 7.0, 9.0, 7.0]))
    assert depth == pytest.approx(0.5 + 2 / 3) and multi is True  # first crossing, 10 -> 7


def _pencil(n_depth=260, n_lat=80, sigma=(3.4, 3.4), mu=(0.0, 0.0), halo=0.0):
    """Canonical dose: a separable Gaussian core per depth, optionally with a flat halo over the plane."""
    ca = lateral_centres(n_lat)
    pa = _integrated_gaussian(ca, 1.0, mu[0], sigma[0])
    pb = _integrated_gaussian(ca, 1.0, mu[1], sigma[1])
    plane = np.outer(pa, pb) + halo / (n_lat * n_lat)
    idd = np.interp(np.arange(n_depth) + 0.5, [0, 200, 240, 250, 260], [1.0, 1.5, 4.0, 0.0, 0.0])
    return idd[:, None, None] * plane[None, :, :]


def test_sigma_and_centre_recovered_per_axis():
    e = endpoints(_pencil(sigma=(3.0, 4.0), mu=(0.3, -0.2)))
    assert e["sigma_100_a"] == pytest.approx(3.0, abs=1e-6)
    assert e["sigma_100_b"] == pytest.approx(4.0, abs=1e-6)
    assert e["sigma_100"] == pytest.approx(3.5, abs=1e-6)
    assert e["centre_200_a"] == pytest.approx(0.3, abs=1e-6)
    assert e["centre_200_b"] == pytest.approx(-0.2, abs=1e-6)


def test_ring_fractions_of_a_flat_plane_are_area_fractions():
    d = np.ones((210, 400, 400))
    e = endpoints(d)
    for lo, hi in ((20, 40), (40, 80)):
        c = lateral_centres(400)
        r = np.sqrt(c[:, None] ** 2 + c[None, :] ** 2)
        assert e[f"ring_100_{lo}_{hi}"] == pytest.approx(((r >= lo) & (r < hi)).sum() / 400**2)
    assert e["sigma_100"] is None  # a flat profile is not a Gaussian: failed, never a number


def test_rings_partition_the_plane_inside_200mm():
    d = _pencil(n_lat=400, halo=0.2)
    e = endpoints(d)
    total = sum(e[f"ring_200_{lo}_{hi}"] for lo, hi in ((0, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 200)))
    c = lateral_centres(400)
    inside = np.sqrt(c[:, None] ** 2 + c[None, :] ** 2) < 200
    slab = d[slab_indices(200)].sum(axis=0)
    assert total == pytest.approx(slab[inside].sum() / slab.sum())  # only the corners beyond 200 mm are unringed
    assert 0.9 < total < 1.0


def test_nonfinite_dose_is_refused():
    d = _pencil()
    d[5, 5, 5] = np.nan
    with pytest.raises(ValueError):
        endpoints(d)


def _write_topas(tmp_path, v, header_bins=None, report="Sum", scorer="Dose", extra_bytes=b"", drop_bytes=0):
    """Write v[ix, iy, kz] as TOPAS does: doubles, x fastest, plus a .binheader."""
    nx, ny, nz = header_bins or v.shape
    payload = np.asfortranarray(v, dtype="<f8").tobytes(order="F")
    payload = payload[: len(payload) - drop_bytes] + extra_bytes
    (tmp_path / "dose.bin").write_bytes(payload)
    (tmp_path / "dose.binheader").write_text(
        f"# TOPAS Version: 4.3\n# Results for scorer: {scorer}\n"
        f"# X in {nx} bins of 0.1 cm\n# Y in {ny} bins of 0.1 cm\n# Z in {nz} bins of 0.1 cm\n"
        f"# DoseToMedium ( Gy ) : {report}   \n"
        "# Binary file: dose.bin\n"
    )
    return tmp_path / "dose.bin"


def test_topas_layout_on_asymmetric_phantom(tmp_path):
    """A marked voxel at TOPAS (ix, iy, kz) lands at canonical (depth from the z = +HLZ face, ix, iy)."""
    nx, ny, nz = 4, 6, 10
    v = np.zeros((nx, ny, nz))
    v[1, 4, nz - 1] = 1.0  # entrance face (the beam comes in at +z)
    v[3, 0, 0] = 2.0  # far face
    arr, dims, spacing = read_topas_bin(_write_topas(tmp_path, v))
    assert dims == [nx, ny, nz] and spacing == pytest.approx([1.0, 1.0, 1.0])
    assert np.array_equal(arr, v)  # x-fastest byte order read back exactly
    d = canonical_from_topas(arr)
    assert d.shape == (nz, nx, ny)
    assert d[0, 1, 4] == 1.0  # entrance voxel is depth bin 0
    assert d[nz - 1, 3, 0] == 2.0  # far face is the last depth bin
    assert d.sum() == 3.0


def test_topas_reader_accepts_the_valid_control(tmp_path):
    v = np.arange(8.0).reshape(2, 2, 2)
    arr, _, _ = read_topas_bin(_write_topas(tmp_path, v))
    assert np.array_equal(arr, v)


@pytest.mark.parametrize("extra", range(1, 8))
def test_topas_reader_rejects_trailing_partial_double(tmp_path, extra):
    with pytest.raises(ValueError):
        read_topas_bin(_write_topas(tmp_path, np.zeros((2, 2, 2)), extra_bytes=b"\0" * extra))


@pytest.mark.parametrize("drop", [1, 7, 8])
def test_topas_reader_rejects_truncation(tmp_path, drop):
    with pytest.raises(ValueError):
        read_topas_bin(_write_topas(tmp_path, np.zeros((2, 2, 2)), drop_bytes=drop))


def test_topas_reader_rejects_header_demanding_more_voxels(tmp_path):
    with pytest.raises(ValueError):
        read_topas_bin(_write_topas(tmp_path, np.zeros((2, 2, 2)), header_bins=(2, 2, 3)))


@pytest.mark.parametrize(
    "report, scorer",
    [("Mean", "DoseSum"), ("Mean", "Dose"), ("Sum Mean", "Dose"), ("Mean Sum", "Dose"), ("Standard_Deviation", "Sum")],
)
def test_topas_reader_rejects_anything_but_a_single_sum_report(tmp_path, report, scorer):
    with pytest.raises(ValueError):
        read_topas_bin(_write_topas(tmp_path, np.zeros((2, 2, 2)), report=report, scorer=scorer))


def test_topas_reader_rejects_zero_bins(tmp_path):
    with pytest.raises(ValueError):
        read_topas_bin(_write_topas(tmp_path, np.zeros((2, 2, 2)), header_bins=(0, 2, 2)))
