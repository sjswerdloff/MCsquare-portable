"""Per-seed endpoints for the cross-platform study (validation/win_linux_equivalence_design.md).

Usage: python platform_study_metrics.py <output dir holding Dose.raw> <field side mm>

Prints one JSON object with the 13 endpoints, computed exactly as the design fixes them:
  lateral_<depth>_<out>  dose <out> mm outside the field edge at <depth> mm, % of the CAX dose there (8 values)
  cax_127, cax_201, cax_i23  absolute CAX dose at 127 mm, 201 mm and voxel index 23 (253 mm)
  r80_mm, r20_mm         distal depths where the CAX depth-dose falls to 80% / 20% of its maximum (None if absent)
plus r80_crossings / r20_crossings (number of downward crossings distal to the peak; >1 is flagged ambiguous).
"""

import json
import math
import sys

import numpy as np

from field_edge_analyse import C, N, OFFSETS_MM, SP, load, side_profile

DEPTH_IY = {127: N - 64, 201: N - 101}  # fixed voxel rows, the same as field_edge_profile.py


def distal_level(depth_mm: np.ndarray, dose: np.ndarray, frac: float) -> tuple[float | None, int]:
    """Depth where `dose` first falls through frac*max distal to its maximum, and the number of such crossings.

    depth_mm must be increasing. The crossing is the first adjacent pair (i, i+1) at or beyond the peak with
    dose[i] >= L > dose[i+1], L = frac * max(dose), linearly interpolated. Returns (None, 0) if there is none.
    """
    level = frac * float(np.max(dose))
    ipk = int(np.argmax(dose))
    first, count = None, 0
    for i in range(ipk, len(dose) - 1):
        if dose[i] >= level > dose[i + 1]:
            count += 1
            if first is None:
                first = float(depth_mm[i] + (dose[i] - level) / (dose[i] - dose[i + 1]) * (depth_mm[i + 1] - depth_mm[i]))
    return first, count


def endpoints(d: np.ndarray, side_mm: float) -> dict:
    half = side_mm / 2
    cax = d[C - 5:C + 5, :, C - 5:C + 5].mean(axis=(0, 2)).astype(np.float64)  # indexed by iy; beam along -y
    out: dict = {}
    for depth, iy in DEPTH_IY.items():
        prof = side_profile(d, iy)
        centre = prof[:5].mean()
        for off in OFFSETS_MM:
            out[f"lateral_{depth}_{off}"] = float(prof[int((half + off) / SP)] / centre * 100)
    out["cax_127"] = float(cax[DEPTH_IY[127]])
    out["cax_201"] = float(cax[DEPTH_IY[201]])
    out["cax_i23"] = float(cax[23])
    depth_mm = (N - np.arange(N) - 0.5) * SP  # depth of each iy
    order = np.argsort(depth_mm)  # increasing depth
    for name, frac in (("r80", 0.8), ("r20", 0.2)):
        r, n = distal_level(depth_mm[order], cax[order], frac)
        out[f"{name}_mm"] = r
        out[f"{name}_crossings"] = n
    for k in ("cax_127", "cax_201", "cax_i23"):
        if not (math.isfinite(out[k]) and out[k] > 0):
            out[f"{k}_invalid"] = True
    return out


def row_of_depth(depth_mm: int) -> int:
    """The voxel row iy whose centre is `depth_mm` below the entrance face. Rows are centred on odd millimetres
    (depth = (N - iy - 0.5) * SP with SP = 2 mm); the +-2 row band of side_profile must lie inside the grid."""
    if SP != 2.0 or isinstance(depth_mm, bool) or not isinstance(depth_mm, int) or depth_mm % 2 != 1:
        raise ValueError(f"depth must be an odd whole number of mm on the 2 mm grid, got {depth_mm!r}")
    iy = N - (depth_mm + 1) // 2
    if not 2 <= iy <= N - 3:
        raise ValueError(f"depth {depth_mm} mm (row {iy}) leaves no +-2 row band inside the {N}-row grid")
    return iy


def field_endpoints(d: np.ndarray, side_mm: float, depths_mm: tuple[int, ...], offsets_mm: tuple[int, ...]) -> dict:
    """The field-edge endpoints at chosen depths and offsets (docs/field_100_150_plan.md), estimators as endpoints().

      lateral_<depth>_<off>  dose <off> mm outside the field edge at <depth> mm, % of the CAX dose there
      cax_<depth>, cax_max   absolute CAX dose (central 10 x 10 voxels) at each depth, and the maximum of that curve
      r80_mm, r20_mm         as endpoints(), with r80_crossings / r20_crossings
    A CAX value that is not finite and positive gets <key>_invalid = True. A lateral value may be exactly 0.
    """
    if d.shape != (N, N, N):
        raise ValueError(f"dose is {d.shape}, expected {(N, N, N)}")
    half = side_mm / 2
    cax = d[C - 5:C + 5, :, C - 5:C + 5].mean(axis=(0, 2)).astype(np.float64)  # indexed by iy; beam along -y
    out: dict = {}
    cax_keys = []
    for depth in depths_mm:
        iy = row_of_depth(depth)
        prof = side_profile(d, iy)
        centre = prof[:5].mean()
        for off in offsets_mm:
            k = int((half + off) / SP)
            if not 0 <= k < len(prof):
                raise ValueError(f"offset {off} mm outside a {side_mm} mm field is beyond the grid (index {k} of {len(prof)})")
            out[f"lateral_{depth}_{off}"] = float(prof[k] / centre * 100)
        out[f"cax_{depth}"] = float(cax[iy])
        cax_keys.append(f"cax_{depth}")
    out["cax_max"] = float(cax.max())
    cax_keys.append("cax_max")
    depth_mm = (N - np.arange(N) - 0.5) * SP  # depth of each iy
    order = np.argsort(depth_mm)  # increasing depth
    for name, frac in (("r80", 0.8), ("r20", 0.2)):
        r, n = distal_level(depth_mm[order], cax[order], frac)
        out[f"{name}_mm"] = r
        out[f"{name}_crossings"] = n
    for k in cax_keys:
        if not (math.isfinite(out[k]) and out[k] > 0):
            out[f"{k}_invalid"] = True
    return out


if __name__ == "__main__":
    print(json.dumps(endpoints(load(sys.argv[1]), float(sys.argv[2])), sort_keys=True))
