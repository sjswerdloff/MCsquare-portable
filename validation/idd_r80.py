"""R80 from a depth-only MCsquare dose grid (an independent scoring grid of 1 x NY x 1 voxels), for the nuclear-off
diagnostic's fine-depth 70 MeV point (#32). Fixed before any fine-grid data exists.

Geometry (validation/mcsquare_pencil_case.py): the beam enters the CT at y = entrance_mm and travels -y, so the depth of
a scoring bin centre is entrance_mm - (Offset_y + (k + 0.5) * spacing_y), all in mm from the MetaImage header MCsquare
writes. R80 is the first downward crossing of 0.8 x max(IDD) distal to the maximum, by linear interpolation between the
two bin-centre depths bracketing it: the same rule as pencil_endpoints.r80, on whatever bin width the header gives.

Usage:
    python idd_r80.py <Dose.mhd> --entrance-mm 350 [--label L]   -> one JSON line prefixed "IDD_R80 "
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# MetaImage headers print %f (6 decimals); an edge within this of the entrance face counts as ON it, not past it.
EDGE_TOL_MM = 1e-4


def read_depth_idd(path: Path, entrance_mm: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Return (depth_mm of bin centres, ascending; IDD in the same order; bin width in mm).

    Raises ValueError on anything that is not a valid depth-only grid lying wholly downstream of the entrance face,
    and on any non-finite or negative dose: invalid input is never reported as a result."""
    if not np.isfinite(entrance_mm):
        raise ValueError(f"entrance {entrance_mm} mm is not finite")
    header = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            header[k.strip()] = v.strip()
    dims = [int(x) for x in header["DimSize"].split()]
    spacing = [float(x) for x in header["ElementSpacing"].split()]
    offset = [float(x) for x in header.get("Offset", "0 0 0").split()]
    if len(dims) != 3 or len(spacing) != 3 or len(offset) != 3:
        raise ValueError(f"{path}: expected 3 values each for DimSize/ElementSpacing/Offset, got {dims} {spacing} {offset}")
    if dims[1] < 2 or not all(np.isfinite(v) and v > 0 for v in spacing) or not all(np.isfinite(v) for v in offset):
        raise ValueError(f"{path}: invalid geometry: dims {dims}, spacing {spacing}, offset {offset}")
    if header.get("ElementType") != "MET_FLOAT" or header.get("ElementByteOrderMSB", "False") != "False":
        raise ValueError(f"{path}: expected little-endian MET_FLOAT")
    if dims[0] != 1 or dims[2] != 1:
        raise ValueError(f"{path}: expected a depth-only grid 1 x NY x 1, got {dims}")
    raw = path.parent / header["ElementDataFile"]
    if raw.stat().st_size != 4 * dims[1]:
        raise ValueError(f"{raw}: {raw.stat().st_size} bytes, header says {dims}")
    idd = np.fromfile(raw, dtype="<f4").astype(np.float64)
    if not np.all(np.isfinite(idd)) or np.any(idd < 0):
        raise ValueError(f"{raw}: dose contains non-finite or negative values")
    upper_edge = offset[1] + dims[1] * spacing[1]  # the grid's upstream edge, in the CT's y (mm)
    if upper_edge > entrance_mm + EDGE_TOL_MM:
        raise ValueError(f"scoring grid extends upstream of the entrance face (upper edge {upper_edge} mm > {entrance_mm} mm)")
    y_centre = offset[1] + (np.arange(dims[1]) + 0.5) * spacing[1]
    depth = entrance_mm - y_centre
    order = np.argsort(depth)
    depth, idd = depth[order], idd[order]
    if not np.all(np.diff(depth) > 0):
        raise ValueError(f"{path}: depth centres are not strictly ascending")
    return depth, idd, spacing[1]


def r80(depth: np.ndarray, idd: np.ndarray) -> tuple[float | None, bool]:
    """First downward 0.8 x max crossing distal to the maximum; (None, False) if there is none."""
    imax = int(np.argmax(idd))
    level = 0.8 * idd[imax]
    below = np.nonzero(idd[imax:] < level)[0]
    if below.size == 0:
        return None, False
    j = imax + int(below[0])
    x0, x1, y0, y1 = depth[j - 1], depth[j], idd[j - 1], idd[j]
    value = float(x0 + (level - y0) * (x1 - x0) / (y1 - y0))
    rest = idd[j:]
    multiple = bool(np.any((rest[:-1] < level) & (rest[1:] >= level)))
    return value, multiple


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mhd", type=Path)
    p.add_argument("--entrance-mm", type=float, required=True)
    p.add_argument("--label", default="")
    a = p.parse_args()
    try:
        depth, idd, width = read_depth_idd(a.mhd, a.entrance_mm)
    except ValueError as e:
        raise SystemExit(str(e)) from None
    if not idd.max() > 0:
        raise SystemExit(f"{a.mhd}: dose is zero everywhere")
    value, multiple = r80(depth, idd)
    imax = int(np.argmax(idd))
    record = {
        "label": a.label,
        "bin_mm": width,
        "n_bins": int(idd.size),
        "depth_range_mm": [float(depth[0]), float(depth[-1])],
        "peak_depth_mm": float(depth[imax]),
        "R80": value,
        "R80_multiple_crossings": multiple,
        "idd_sum_raw": float(idd.sum()),
    }
    print("IDD_R80 " + json.dumps(record, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
