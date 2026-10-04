"""Tabulate the #16 field-edge matrix: dose outside the field edge, before vs after the fix.

Usage: analyse.py <run dir> <case>... For each case E<energy>_S<side>, reads out_<case>_before_A, out_<case>_after_A and
out_<case>_before_B (same code as before_A, different seed: the noise floor). Every number is a
percentage of the central-axis dose at the same depth:
  out_mm     distance outside the geometric field edge
  before     dose there with main's code
  fix-bef    change the fix makes
  noise      seed B minus seed A with main's code
  noise_fix  seed B minus seed A with the fix, when out_<case>_after_B exists
The CAX rows give the relative change in central-axis dose (after/before - 1) at that depth.
"""

import os
import sys

import numpy as np

N = int(os.environ.get("FE_N", "150"))  # must match the cube the runs used
SP, C = 2.0, N // 2
OFFSETS_MM = (5, 10, 20, 30)


def load(path: str) -> np.ndarray:
    return np.fromfile(f"{path}/Dose.raw", dtype=np.float32).reshape(N, N, N)  # z, y, x; beam along -y


def side_profile(d: np.ndarray, iy: int) -> np.ndarray:
    """Lateral profile vs distance from the central axis: 4-side mean over a +-2 voxel depth band."""
    band = d[:, iy - 2:iy + 3, :].mean(axis=1)
    px = band[C - 10:C + 10, :].mean(axis=0)
    pz = band[:, C - 10:C + 10].mean(axis=1)
    halves = [px[C:], px[C - 1::-1], pz[C:], pz[C - 1::-1]]
    return np.mean(halves, axis=0)


def main(root: str, cases: list[str]) -> None:
    print(f"{'case':10} {'depth':>12} {'mm':>5} {'out_mm':>6} {'before%':>8} {'fix-bef':>8} {'noise':>7} {'noise_fix':>9}")
    for case in cases:
        side = float(case.split("_S")[1])
        half = side / 2
        a, f, b = (load(f"{root}/out_{case}_{v}") for v in ("before_A", "after_A", "before_B"))
        fb = load(f"{root}/out_{case}_after_B") if os.path.exists(f"{root}/out_{case}_after_B/Dose.raw") else None
        cax = lambda d: d[C - 5:C + 5, :, C - 5:C + 5].mean(axis=(0, 2))  # noqa: E731
        ca, cf, cb = cax(a), cax(f), cax(b)
        cfb = cax(fb) if fb is not None else None
        ipk = int(np.argmax(ca))
        depth = lambda iy: (N - iy - 0.5) * SP  # noqa: E731
        picks = {"mid-range": (N + ipk) // 2, "80% peak": N - int(0.8 * (N - ipk)), "peak": ipk}
        for label, iy in picks.items():
            print(f"{case:10} {'CAX ' + label:>12} {depth(iy):5.0f} {'':>6} {'':>8} "
                  f"{(cf[iy] / ca[iy] - 1) * 100:8.3f} {(cb[iy] / ca[iy] - 1) * 100:7.3f} "
                  + (f"{(cfb[iy] / cf[iy] - 1) * 100:9.3f}" if cfb is not None else ""))
            if label == "peak":
                continue
            pa, pf, pb = side_profile(a, iy), side_profile(f, iy), side_profile(b, iy)
            pfb = side_profile(fb, iy) if fb is not None else None
            dc = pa[:5].mean()
            for off in OFFSETS_MM:
                k = int((half + off) / SP)
                print(f"{case:10} {label:>12} {depth(iy):5.0f} {off:6d} {pa[k] / dc * 100:8.3f} "
                      f"{(pf[k] - pa[k]) / dc * 100:8.3f} {(pb[k] - pa[k]) / dc * 100:7.3f} "
                      + (f"{(pfb[k] - pf[k]) / dc * 100:9.3f}" if pfb is not None else ""))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
