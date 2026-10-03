"""Lenovo (upstream + fix, icl) against Studio (Portable arm64) Lungman dose: where the beam lands, tumour dose."""
import json
import sys
from pathlib import Path

import numpy as np

L = Path(sys.argv[1])  # Lenovo run dir
S = Path(sys.argv[2])  # Studio package dir
M = Path(sys.argv[3])  # Lungman/mcsquare dir
shape = (426, 512, 512)  # z, y, x for DimSize 512 512 426
gy_per_mu = 1.602176e-16 * 114797310

ctype = {"MET_SHORT": np.int16, "MET_FLOAT": np.float32}[
    next(l.split("=")[1].strip() for l in (M / "CT.mhd").read_text().splitlines() if l.startswith("ElementType"))]
ct = np.fromfile(M / "CT.raw", dtype=ctype).reshape(shape)
mask = np.fromfile(M / "tumours_100HU_1_mask.raw", dtype=np.uint8).reshape(shape) > 0
dl = np.fromfile(L / "Outputs" / "Dose.raw", dtype=np.float32).reshape(shape) * gy_per_mu
ds = np.fromfile(S / "Outputs" / "Dose.raw", dtype=np.float32).reshape(shape) * gy_per_mu
body = ct > -900
spacing = np.array([0.7, 0.625, 0.625])  # z, y, x mm


def centroid(d: np.ndarray, sel: np.ndarray) -> np.ndarray:
    idx = np.argwhere(sel)
    w = d[sel]
    return (idx * w[:, None]).sum(0) / w.sum() * spacing


out = {}
for name, d in (("lenovo", dl), ("studio", ds)):
    t = d[mask]
    hi = body & (d > 0.5 * d[body].max())
    out[name] = {"tumour_mean_gy_per_mu": float(t.mean()), "tumour_min": float(t.min()), "tumour_max": float(t.max()),
                 "high_dose_centroid_zyx_mm": [round(float(v), 2) for v in centroid(d, hi)],
                 "high_dose_voxels": int(hi.sum())}
tc = np.argwhere(mask).mean(0) * spacing
out["tumour_centroid_zyx_mm"] = [round(float(v), 2) for v in tc]
out["centroid_shift_lenovo_minus_studio_mm"] = [round(a - b, 2) for a, b in zip(
    out["lenovo"]["high_dose_centroid_zyx_mm"], out["studio"]["high_dose_centroid_zyx_mm"], strict=True)]
out["tumour_mean_ratio_lenovo_over_studio"] = out["lenovo"]["tumour_mean_gy_per_mu"] / out["studio"]["tumour_mean_gy_per_mu"]
sel = body & (ds > 0.1 * ds[body].max())
diff = (dl - ds)[sel] / ds[body].max()
out["body_above_10pct_dose_diff_pct_of_studio_max"] = {"n": int(sel.sum()), "mean": float(diff.mean() * 100),
                                                        "p95_abs": float(np.percentile(np.abs(diff), 95) * 100)}
out["tumour_voxels_with_lenovo_dose_zero"] = int((dl[mask] == 0).sum())
print(json.dumps(out, indent=1))
