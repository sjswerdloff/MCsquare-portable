"""Stage-1 pencil-beam endpoints (issue #32, validation/topas_design.md, "Endpoints").

Every endpoint is computed on one CANONICAL array D[k, a, b]:
  k  depth bin, centre (k + 0.5) mm below the entrance surface
  a, b  lateral bins, 1 mm, centres (i + 0.5 - N/2) mm, so the beam axis is the shared corner of the four
        central voxels (x = y = 0 between indices N/2 - 1 and N/2)
Loaders map each code's own layout onto it.
  MCsquare  established on an asymmetric probe (gantry 0: MCsquare's internal y is the beam axis, and the beam
            enters through the y = NY face, so depth bin k is MCsquare index j = NY - 1 - k).
  TOPAS     validation/topas/: binary doubles, x fastest (Fortran order), so the file is v[ix, iy, kz]; the beam
            enters at z = +HLZ travelling -z, so depth bin k is TOPAS index kz = NZ - 1 - k. Lateral bins are
            centred on the phantom, which puts the axis on the shared corner, as the canonical array requires.

Endpoints, per run:
  R80            whole-plane IDD, first downward crossing of 0.8 x max distal to the maximum, linear
                 interpolation between the voxel-centre depths bracketing it (the #31 rule). None if absent;
                 more than one downward crossing is flagged.
  sigma_{d}      slab = six bins with centres d - 2.5 ... d + 2.5 mm. Profiles projected onto each lateral
                 axis, fitted over |x| <= 10 mm with a voxel-integrated Gaussian (free amplitude, centre, sigma,
                 no background). sigma = mean of the two fits; failed if it does not converge or leaves [1, 20].
  ring_{d}_{lo}_{hi}  energy fraction in the slab: dose in voxels whose centre radius is in [lo, hi), divided
                 by the dose over the whole scored plane in the slab.
Diagnostics: sigma per axis, fitted centres (the centroid check), absolute IDD maximum, and the entrance slab
(d = 3 mm), whose sigma checks the sampled spot against the BDL on the real build.

Usage:
  python pencil_endpoints.py mcsquare <Dose.mhd> [--label TEXT]      prints one JSON line
  python pencil_endpoints.py topas <dose.bin> [--label TEXT]         reads <dose>.binheader beside it
"""

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf

SLAB_DEPTHS_MM = (3, 100, 200)  # 3 = entrance slab (bins 0-5): beam-model diagnostic, sigma ~ BDL spot size
FIT_HALF_WIDTH_MM = 10.0
SIGMA_VALID_MM = (1.0, 20.0)
RINGS_MM = ((0, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 200))
CONFIRMATORY_RINGS_MM = ((20, 40), (40, 80))


def read_mhd(path: Path) -> tuple[np.ndarray, list[int], list[float]]:
    """Read an MHD/raw float image; returns the array in numpy (z, y, x) order, DimSize and spacing (x, y, z)."""
    header = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            header[key.strip()] = value.strip()
    dims = [int(v) for v in header["DimSize"].split()]
    spacing = [float(v) for v in header["ElementSpacing"].split()]
    if header.get("ElementType") != "MET_FLOAT" or header.get("ElementByteOrderMSB", "False") != "False":
        raise ValueError(f"{path}: expected little-endian MET_FLOAT")
    raw = path.parent / header["ElementDataFile"]
    expected = 4 * dims[0] * dims[1] * dims[2]
    # Byte length, not element count: np.fromfile silently drops a trailing partial value (alden-ec2221c7, #38).
    if raw.stat().st_size != expected:
        raise ValueError(f"{raw}: {raw.stat().st_size} bytes, header says {dims} ({expected} bytes)")
    data = np.fromfile(raw, dtype="<f4")
    return data.reshape(dims[2], dims[1], dims[0]), dims, spacing


_TOPAS_AXIS = re.compile(r"^#\s*([XYZ]) in (\d+) bins? of ([0-9.eE+-]+) (cm|mm)\s*$")
# The quantity line, e.g. "# DoseToMedium ( Gy ) : Sum". The reports follow the colon after the unit.
_TOPAS_REPORT = re.compile(r"^#\s*\S+\s*\(\s*[^()]*\)\s*:\s*(.+?)\s*$")


def read_topas_bin(path: Path) -> tuple[np.ndarray, list[int], list[float]]:
    """Read a TOPAS binary Sum scorer; returns v[ix, iy, kz], DimSize (x, y, z) and spacing in mm (x, y, z).

    Dimensions and spacing come from the .binheader written beside the .bin, never assumed."""
    header = path.with_suffix(".binheader")
    axes: dict[str, tuple[int, float]] = {}
    reports: list[list[str]] = []
    for line in header.read_text().splitlines():
        m = _TOPAS_AXIS.match(line.strip())
        if m:
            if m.group(1) in axes:
                raise ValueError(f"{header}: axis {m.group(1)} declared twice")
            width = float(m.group(3)) * (10.0 if m.group(4) == "cm" else 1.0)
            axes[m.group(1)] = (int(m.group(2)), width)
            continue
        r = _TOPAS_REPORT.match(line.strip())
        if r:
            reports.append(r.group(1).split())
    if set(axes) != {"X", "Y", "Z"}:
        raise ValueError(f"{header}: expected X, Y and Z bin lines, found {sorted(axes)}")
    # Exactly one quantity line, reporting exactly one column, and that column is Sum. A multi-report file
    # interleaves columns, so reading it as one Sum array would be silently wrong.
    if reports != [["Sum"]]:
        raise ValueError(f"{header}: expected exactly one quantity line reporting only Sum, found {reports}")
    dims = [axes[a][0] for a in "XYZ"]
    spacing = [axes[a][1] for a in "XYZ"]
    if any(n <= 0 for n in dims) or not all(math.isfinite(w) and w > 0 for w in spacing):
        raise ValueError(f"{header}: non-positive or non-finite geometry {dims} {spacing}")
    expected_bytes = 8 * dims[0] * dims[1] * dims[2]
    actual_bytes = path.stat().st_size
    # Check BYTES, not values: np.fromfile silently drops a trailing partial double.
    if actual_bytes != expected_bytes:
        raise ValueError(f"{path}: {actual_bytes} bytes, header implies {expected_bytes}")
    data = np.fromfile(path, dtype="<f8")
    return data.reshape(dims, order="F"), dims, spacing


def canonical_from_topas(xyz: np.ndarray) -> np.ndarray:
    """TOPAS v[ix, iy, kz], beam entering at z = +HLZ and travelling to -z, to D[k, x, y]."""
    return np.ascontiguousarray(xyz[:, :, ::-1].transpose(2, 0, 1))


def canonical_from_mcsquare(zyx: np.ndarray) -> np.ndarray:
    """MCsquare (z, y, x), beam entering at y = NY and travelling to -y, to D[k, z, x]."""
    return np.ascontiguousarray(zyx[:, ::-1, :].transpose(1, 0, 2))


def lateral_centres(n: int) -> np.ndarray:
    """Voxel-centre coordinates (mm) of n 1 mm lateral bins with the axis on the central corner."""
    if n % 2:
        raise ValueError(f"lateral size {n} is odd: the axis cannot sit on the shared corner")
    return np.arange(n) + 0.5 - n / 2


def r80(idd: np.ndarray) -> tuple[float | None, bool]:
    """First downward crossing of 0.8 max distal to the maximum; returns (depth mm or None, multiple-crossing flag)."""
    depths = np.arange(idd.size) + 0.5
    imax = int(np.argmax(idd))
    level = 0.8 * idd[imax]
    crossings = [i for i in range(imax, idd.size - 1) if idd[i] >= level > idd[i + 1]]
    if not crossings:
        return None, False
    i = crossings[0]
    frac = (idd[i] - level) / (idd[i] - idd[i + 1])
    return float(depths[i] + frac * (depths[i + 1] - depths[i])), len(crossings) > 1


def slab_indices(depth_mm: float) -> np.ndarray:
    """The six depth bins whose centres are depth - 2.5 ... depth + 2.5 mm."""
    d = round(depth_mm)
    if d != depth_mm:
        raise ValueError("slab depths are whole millimetres")
    return np.arange(d - 3, d + 3)


def _integrated_gaussian(x: np.ndarray, amp: float, mu: float, sigma: float) -> np.ndarray:
    """Dose in 1 mm bins centred at x from a Gaussian of total amp: erf difference over each bin."""
    s = sigma * math.sqrt(2.0)
    return amp * 0.5 * (erf((x + 0.5 - mu) / s) - erf((x - 0.5 - mu) / s))


def fit_sigma(profile: np.ndarray, centres: np.ndarray) -> tuple[float | None, float | None]:
    """Voxel-integrated Gaussian fit over |x| <= 10 mm; returns (sigma, centre) or (None, None) if failed."""
    sel = np.abs(centres) <= FIT_HALF_WIDTH_MM
    x, y = centres[sel], profile[sel]
    if not np.all(np.isfinite(y)) or y.sum() <= 0:
        return None, None
    p0 = (float(y.sum()), float((x * y).sum() / y.sum()), 4.0)
    try:
        popt, _ = curve_fit(_integrated_gaussian, x, y, p0=p0, maxfev=10000)
    except (RuntimeError, ValueError):
        return None, None
    sigma, mu = abs(float(popt[2])), float(popt[1])
    if not (SIGMA_VALID_MM[0] <= sigma <= SIGMA_VALID_MM[1]) or not math.isfinite(mu):
        return None, None
    return sigma, mu


def endpoints(d: np.ndarray) -> dict:
    """All endpoints and diagnostics of one canonical dose array D[k, a, b]."""
    if not np.all(np.isfinite(d)):
        raise ValueError("dose contains non-finite values")
    ca, cb = lateral_centres(d.shape[1]), lateral_centres(d.shape[2])
    radius = np.sqrt(ca[:, None] ** 2 + cb[None, :] ** 2)
    idd = d.sum(axis=(1, 2))
    out: dict = {"idd_max": float(idd.max())}
    out["R80"], out["R80_multiple_crossings"] = r80(idd)
    for depth in SLAB_DEPTHS_MM:
        slab = d[slab_indices(depth)].sum(axis=0)
        sa, mua = fit_sigma(slab.sum(axis=1), ca)
        sb, mub = fit_sigma(slab.sum(axis=0), cb)
        out[f"sigma_{depth}_a"], out[f"sigma_{depth}_b"] = sa, sb
        out[f"centre_{depth}_a"], out[f"centre_{depth}_b"] = mua, mub
        out[f"sigma_{depth}"] = None if sa is None or sb is None else (sa + sb) / 2
        total = float(slab.sum())
        for lo, hi in RINGS_MM:
            ring = float(slab[(radius >= lo) & (radius < hi)].sum())
            out[f"ring_{depth}_{lo}_{hi}"] = ring / total if total > 0 else None
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("layout", choices=["mcsquare", "topas"])
    p.add_argument("dose", type=Path)
    p.add_argument("--label", default="")
    a = p.parse_args()
    if a.layout == "mcsquare":
        raw, dims, spacing = read_mhd(a.dose)
        canonical = canonical_from_mcsquare
    else:
        raw, dims, spacing = read_topas_bin(a.dose)
        canonical = canonical_from_topas
    if not np.allclose(spacing, [1.0, 1.0, 1.0]):
        raise SystemExit(f"expected 1 mm voxels, got {spacing}")
    record = {"label": a.label, "layout": a.layout, "dims": dims, **endpoints(canonical(raw))}
    print("ENDPOINTS " + json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
