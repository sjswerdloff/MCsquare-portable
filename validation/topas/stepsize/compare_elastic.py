"""Elastic-off TOPAS against the halo-addendum control, 200 MeV (see PREDICTION_elastic.md).

Estimators as in validation/topas_halo_compare.py: ring fractions as a geometric-mean ratio (difference of mean
logs) with a 95% Welch interval; R80 and sigma as a difference of run-level means with a 95% Welch interval.
Diagnostic only.

Usage: compare_elastic.py <control.jsonl> <elastic_off.jsonl>
"""

import json
import math
import sys

from scipy import stats


def load(path: str) -> list[dict]:
    return [json.loads(line.split(" ", 1)[1]) for line in open(path) if line.strip()]


def welch(a: list[float], b: list[float]) -> tuple[float, float, float]:
    """Mean of b minus mean of a, and its 95% Welch interval."""
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    se2 = va / len(a) + vb / len(b)
    df = se2**2 / ((va / len(a)) ** 2 / (len(a) - 1) + (vb / len(b)) ** 2 / (len(b) - 1))
    t = stats.t.ppf(0.975, df)
    d = mb - ma
    return d, d - t * math.sqrt(se2), d + t * math.sqrt(se2)


ctrl, off = load(sys.argv[1]), load(sys.argv[2])
print(f"control runs {len(ctrl)}, elastic-off runs {len(off)}")
keys = sorted(k for k in ctrl[0] if k == "R80" or k.startswith(("sigma_", "ring_")) and not k.endswith(("_a", "_b")))
for k in keys:
    a = [r[k] for r in ctrl if r.get(k) is not None]
    b = [r[k] for r in off if r.get(k) is not None]
    if len(a) < 2 or len(b) < 2:
        print(f"{k:22s} not computed (runs {len(a)}, {len(b)})")
        continue
    if k.startswith("ring_"):
        if min(a + b) <= 0:
            print(f"{k:22s} a run scored zero; no ratio")
            continue
        d, lo, hi = welch([math.log(x) for x in a], [math.log(x) for x in b])
        print(f"{k:22s} off/control {math.exp(d):.3f}  95% [{math.exp(lo):.3f}, {math.exp(hi):.3f}]")
    else:
        d, lo, hi = welch(a, b)
        print(f"{k:22s} off - control {d:+.4f} mm  95% [{lo:+.4f}, {hi:+.4f}]")
