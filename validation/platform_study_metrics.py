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


if __name__ == "__main__":
    print(json.dumps(endpoints(load(sys.argv[1]), float(sys.argv[2])), sort_keys=True))
