"""One platform's field-edge numbers, for comparing two platforms without moving dose files between them.

Usage: python field_edge_profile.py <root> <case> <version>   e.g.  run E200_S150 before

Reads <root>/out_<case>_<version>_A and _B (two seeds). For each seed, and for their mean, prints the
central-axis dose (absolute, as MCsquare wrote it) at mid-range, 80%-of-peak and peak depth, and the
lateral dose 5/10/20/30 mm outside the field edge as % of the central dose at the same depth. The depths
are picked from seed A's own central-axis peak, exactly as field_edge_analyse.py picks them.
"""

import sys

import numpy as np

from field_edge_analyse import C, N, OFFSETS_MM, SP, load, side_profile


def main(root: str, case: str, version: str) -> None:
    """Print the per-seed and mean profile table for one platform."""
    half = float(case.split("_S")[1]) / 2
    doses = {s: load(f"{root}/out_{case}_{version}_{s}") for s in ("A", "B")}
    cax = {s: d[C - 5:C + 5, :, C - 5:C + 5].mean(axis=(0, 2)) for s, d in doses.items()}
    ipk = int(np.argmax(cax["A"]))
    picks = {"mid-range": (N + ipk) // 2, "80% peak": N - int(0.8 * (N - ipk)), "peak": ipk}
    print(f"{'depth':>12} {'mm':>5} {'out_mm':>6} {'seed A':>12} {'seed B':>12} {'mean':>12}")
    for label, iy in picks.items():
        depth = (N - iy - 0.5) * SP
        va, vb = (float(cax[s][iy]) for s in ("A", "B"))
        print(f"{'CAX ' + label:>12} {depth:5.0f} {'':>6} {va:12.6e} {vb:12.6e} {(va + vb) / 2:12.6e}")
        if label == "peak":
            continue
        prof = {s: side_profile(d, iy) for s, d in doses.items()}
        for off in OFFSETS_MM:
            k = int((half + off) / SP)
            pa, pb = (prof[s][k] / prof[s][:5].mean() * 100 for s in ("A", "B"))
            print(f"{label:>12} {depth:5.0f} {off:6d} {pa:12.4f} {pb:12.4f} {(pa + pb) / 2:12.4f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
