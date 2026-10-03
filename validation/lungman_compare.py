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

# Axis profile: the beam travels along (sin G, -cos G, 0) in array (x, y, z) terms (G = 135; Studio README), through
# the tumour centroid. Nearest-voxel samples every 0.5 mm from 80 mm upstream to 60 mm downstream of the centroid.
g = np.deg2rad(135.0)
direction_zyx = np.array([0.0, -np.cos(g), np.sin(g)])
steps = np.arange(-80.0, 60.0 + 1e-9, 0.5)
pts = (tc[None, :] + steps[:, None] * direction_zyx[None, :]) / spacing
idx = np.rint(pts).astype(int)
ok = np.all((idx >= 0) & (idx < np.array(shape)), axis=1)
pl, ps = dl[tuple(idx[ok].T)], ds[tuple(idx[ok].T)]
s_ok = steps[ok]


def distal_r80(prof: np.ndarray) -> float:
    i = int(np.argmax(prof))
    level = 0.8 * prof[i]
    for k in range(i, len(prof) - 1):
        if prof[k] >= level > prof[k + 1]:
            return float(s_ok[k] + (prof[k] - level) / (prof[k] - prof[k + 1]) * (s_ok[k + 1] - s_ok[k]))
    return float("nan")


out["axis_profile"] = {"step_mm": 0.5, "from_mm": float(s_ok[0]), "to_mm": float(s_ok[-1]),
                       "max_abs_diff_pct_of_studio_axis_max": float(np.abs(pl - ps).max() / ps.max() * 100),
                       "peak_mm": {"lenovo": float(s_ok[np.argmax(pl)]), "studio": float(s_ok[np.argmax(ps)])},
                       "distal_r80_mm_from_centroid": {"lenovo": distal_r80(pl), "studio": distal_r80(ps)}}

# Gamma 2%/2 mm, global, 10% cutoff, Studio as reference, on the box holding every body voxel above 10% of the body
# maximum; pass rate over body voxels (HU > -900) only (air voxels carry dose-per-gram noise; cora-2f1e43dc).
import pymedphys  # noqa: E402

sel10 = body & (ds > 0.1 * ds[body].max())
lo, hi = np.argwhere(sel10).min(0) - 5, np.argwhere(sel10).max(0) + 6
lo, hi = np.maximum(lo, 0), np.minimum(hi, np.array(shape))
box = tuple(slice(a, b) for a, b in zip(lo, hi, strict=True))
axes = tuple(np.arange(n_) * sp for n_, sp in zip(np.array(hi) - np.array(lo), spacing, strict=True))
gam = pymedphys.gamma(axes, ds[box], axes, dl[box], dose_percent_threshold=2, distance_mm_threshold=2,
                      lower_percent_dose_cutoff=10, global_normalisation=float(ds[body].max()), max_gamma=2,
                      random_subset=None)
evaluated = np.isfinite(gam) & body[box]
out["gamma_2pct_2mm_body"] = {"points": int(evaluated.sum()), "pass_rate_pct": float((gam[evaluated] <= 1).mean() * 100),
                              "box_zyx": [[int(a), int(b)] for a, b in zip(lo, hi, strict=True)]}
print(json.dumps(out, indent=1))
