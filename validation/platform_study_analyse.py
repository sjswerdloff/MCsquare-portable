"""Pre-registered analysis for the cross-platform study (validation/win_linux_equivalence_design.md).

Usage: python platform_study_analyse.py <records.jsonl> <look: 1 or 2>

<records.jsonl> holds one STUDY_RESULT record per line (the JSON after "STUDY_RESULT "), for every seed of every
wave so far. All records are pooled; none is dropped. Needs numpy and scipy.

Per comparison (windows vs linux, macos vs linux), per endpoint, Welch two-sample on per-seed values
(ln dose for CAX), and:
  TOST       95% CI of W-L (alpha = 0.025 one-sided at each look) inside the margin      -> endpoint passes
  NONEQ      simultaneous CI at 1 - 0.025/26 two-sided entirely outside the margin      -> confirmatory non-equivalence
  FLAG       plain 95% CI entirely outside the margin (exploratory only)
  Holm       two-sided Welch p-values, Holm across all 26 endpoint-comparisons at 0.025 (secondary)
Verdict per comparison: EQUIVALENT if all 13 pass; NON-EQUIVALENT if any NONEQ; otherwise CONTINUE (look 1)
or INCONCLUSIVE (look 2). An endpoint with any missing R, or any invalid CAX dose, cannot pass.
"""

import json
import math
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

ALPHA_LOOK = 0.025
N_FAMILY = 26
LATERAL = [f"lateral_{d}_{o}" for d in (127, 201) for o in (5, 10, 20, 30)]
CAX = ["cax_127", "cax_201", "cax_i23"]
RANGE = ["r80_mm", "r20_mm"]
ENDPOINTS = LATERAL + CAX + RANGE


def margin(ep: str) -> tuple[float, float]:
    """(lower, upper) bounds for W - L on the analysis scale (ln ratio for CAX)."""
    if ep in CAX:
        return math.log(1 - 0.005), math.log(1 + 0.005)
    if ep in RANGE:
        return -0.5, 0.5
    d = 0.50 if ep.endswith("_5") else 0.20
    return -d, d


def values(recs: list[dict], ep: str) -> tuple[np.ndarray, str | None]:
    """Per-seed values on the analysis scale, or a reason the endpoint cannot be analysed."""
    out = []
    for r in recs:
        v = r["metrics"].get(ep)
        if v is None:
            return np.array([]), "missing (no crossing) in seed %d" % r["seed"]
        if ep in CAX:
            if not (math.isfinite(v) and v > 0):
                return np.array([]), "invalid dose in seed %d" % r["seed"]
            v = math.log(v)
        out.append(float(v))
    return np.array(out), None


def welch(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Mean difference a-b, its standard error and Welch-Satterthwaite df."""
    va, vb, na, nb = a.var(ddof=1), b.var(ddof=1), len(a), len(b)
    se2 = va / na + vb / nb
    df = se2 ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)) if se2 > 0 else float("inf")
    return float(a.mean() - b.mean()), math.sqrt(se2), df


def main(path: str, look: int) -> None:
    recs = [json.loads(line) for line in open(path) if line.strip()]
    by = defaultdict(list)
    for r in recs:
        by[r["platform"]].append(r)
    seen = set()
    for r in recs:
        key = (r["platform"], r["seed"])
        if key in seen:
            raise SystemExit(f"duplicate record for {key}: refusing to analyse")
        seen.add(key)
    print(f"look {look}; seeds per platform: " + ", ".join(f"{p} {len(v)}" for p, v in sorted(by.items())))
    for p, v in sorted(by.items()):
        hosts = defaultdict(int)
        for r in v:
            hosts[r["host"]] += 1
        flags = {ep: sum(1 for r in v if r["metrics"].get(ep.replace("_mm", "_crossings"), 1) > 1) for ep in RANGE}
        print(f"  {p}: hosts {dict(hosts)}; seeds with ambiguous crossings {flags}")
    rows, pvals = [], []
    for plat in ("windows", "macos"):
        if plat not in by or "linux" not in by:
            continue
        for ep in ENDPOINTS:
            a, why_a = values(by[plat], ep)
            b, why_b = values(by["linux"], ep)
            lo_m, hi_m = margin(ep)
            if why_a or why_b:
                rows.append((plat, ep, None, why_a or why_b, lo_m, hi_m))
                continue
            d, se, df = welch(a, b)
            t95 = stats.t.ppf(1 - ALPHA_LOOK, df)
            tsim = stats.t.ppf(1 - ALPHA_LOOK / N_FAMILY / 2, df)
            p = 2 * stats.t.sf(abs(d) / se, df) if se > 0 else (0.0 if d else 1.0)
            ci = (d - t95 * se, d + t95 * se)
            sim = (d - tsim * se, d + tsim * se)
            rows.append((plat, ep, dict(d=d, se=se, df=df, ci=ci, sim=sim, p=p), None, lo_m, hi_m))
            pvals.append((p, len(rows) - 1))
    # Holm across all analysable endpoint-comparisons, at the per-look alpha.
    holm = set()
    for rank, (p, i) in enumerate(sorted(pvals)):
        if p <= ALPHA_LOOK / (len(pvals) - rank):
            holm.add(i)
        else:
            break
    verdict = {}
    print(f"\n{'comparison':9} {'endpoint':16} {'W-L':>10} {'95% CI (TOST)':>24} {'sim CI':>24} {'margin':>18} result")
    for i, (plat, ep, s, why, lo_m, hi_m) in enumerate(rows):
        scale = (lambda x: 100 * (math.exp(x) - 1)) if ep in CAX else (lambda x: x)
        unit = "%" if ep in CAX else ("mm" if ep in RANGE else "pts")
        mtxt = f"[{scale(lo_m):+.3f},{scale(hi_m):+.3f}]{unit}"
        v = verdict.setdefault(plat, {"pass": 0, "noneq": [], "fail": []})
        if s is None:
            v["fail"].append(ep)
            print(f"{plat:9} {ep:16} {'':>10} {'':>24} {'':>24} {mtxt:>18} CANNOT PASS: {why}")
            continue
        passes = lo_m < s["ci"][0] and s["ci"][1] < hi_m
        noneq = s["sim"][0] > hi_m or s["sim"][1] < lo_m
        flag = s["ci"][0] > hi_m or s["ci"][1] < lo_m
        v["pass"] += passes
        if noneq:
            v["noneq"].append(ep)
        tags = ["PASS" if passes else "not shown"] + (["NONEQ"] if noneq else []) + (["flag"] if flag and not noneq else []) \
            + (["Holm-detected"] if i in holm else [])
        print(f"{plat:9} {ep:16} {scale(s['d']):+10.4f} [{scale(s['ci'][0]):+10.4f},{scale(s['ci'][1]):+10.4f}] "
              f"[{scale(s['sim'][0]):+10.4f},{scale(s['sim'][1]):+10.4f}] {mtxt:>18} {' '.join(tags)}")
    print()
    for plat, v in verdict.items():
        if v["pass"] == len(ENDPOINTS):
            res = "EQUIVALENT"
        elif v["noneq"]:
            res = "NON-EQUIVALENT on " + ", ".join(v["noneq"])
        else:
            res = "CONTINUE to wave 2" if look == 1 else "INCONCLUSIVE"
        print(f"VERDICT {plat} vs linux (look {look}): {res}  ({v['pass']}/{len(ENDPOINTS)} endpoints pass)")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
