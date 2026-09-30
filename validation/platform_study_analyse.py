"""Pre-registered analysis for the cross-platform study (validation/win_linux_equivalence_design.md).

Usage: python platform_study_analyse.py <records.jsonl> <look: 1 or 2> <look-1 verdicts.json>

<records.jsonl> holds one STUDY_RESULT record per line (the JSON after "STUDY_RESULT "), for every seed of every
wave so far. All records are pooled; none is dropped. Needs numpy and scipy.

VALIDATION, before any inference (fails closed; added 2026-10-01 on alden-ec2221c7's review 6821, before any
STUDY_RESULT was read): look must be 1 or 2; study id must be pe1; every record's (platform, seed) must be in the
expected set for that look and every expected seed present exactly once; each record's commit must be the frozen
study commit, or the re-run commit for exactly the two recorded interrupted seeds; within each platform every shared
input hash must be identical; cube.raw (written in binary mode) must be identical across ALL platforms; text inputs
may differ in bytes across platforms only (git line-ending translation on Windows), which is reported, since their
content is fixed by the commit. A non-finite endpoint value makes that endpoint unable to pass.

STOP STATE: look 1 writes its verdicts to <look-1 verdicts.json> (refusing to overwrite). Look 2 reads it; a
comparison that stopped at look 1 is reported as frozen and is not recomputed, and does not enter look 2's Holm
family. Look 2 requires the wave-2 commit to be frozen in WAVE2_COMMIT first.

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
import os
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

STUDY = "pe1"
STUDY_COMMIT = "7d07db0145671924428936d5a7fa896d8c3f8e22"
RERUN_COMMIT = "600c66e0fe19adaaee49e51fa023edf89c47f299"  # adds only the re-run workflow to STUDY_COMMIT
RERUN_SEEDS = {("windows", 1008), ("linux", 2003)}
WAVE2_COMMIT = None  # frozen here, before wave 2 starts, if wave 2 is run
WAVE_SEEDS = {1: {"windows": range(1001, 1017), "linux": range(2001, 2017), "macos": range(3001, 3017)},
              2: {"windows": range(1017, 1033), "linux": range(2017, 2033), "macos": range(3017, 3033)}}
COMPARED = ("windows", "macos")
TEXT_INPUTS = ["plan E200_S150.txt", "cube.mhd", "BDL", "HU_Density", "HU_Material"]


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
        if not math.isfinite(v):
            return np.array([]), "non-finite value in seed %d" % r["seed"]
        if ep in CAX:
            if not v > 0:
                return np.array([]), "invalid dose in seed %d" % r["seed"]
            v = math.log(v)
        out.append(float(v))
    return np.array(out), None


def refuse(why: str) -> None:
    raise SystemExit("REFUSING TO ANALYSE: " + why)


def expected(look: int, continuing: set[str]) -> set[tuple[str, int]]:
    """Every (platform, seed) the data must hold at this look."""
    exp = {(p, s) for p, r in WAVE_SEEDS[1].items() for s in r}
    if look == 2:
        for p in continuing | {"linux"}:
            exp |= {(p, s) for s in WAVE_SEEDS[2][p]}
    return exp


def validate(recs: list[dict], look: int, continuing: set[str]) -> None:
    """Fail closed before inference on anything but the complete, expected, same-provenance dataset."""
    if not recs:
        refuse("no records")
    for r in recs:
        if r.get("study") != STUDY:
            refuse(f"record from study {r.get('study')!r}, expected {STUDY!r}")
    keys = [(r.get("platform"), r.get("seed")) for r in recs]
    dup = {k for k in keys if keys.count(k) > 1}
    if dup:
        refuse(f"duplicate records for {sorted(dup)}")
    got, exp = set(keys), expected(look, continuing)
    if got != exp:
        refuse(f"seed set is not the expected one for look {look}: missing {sorted(exp - got)}, unexpected {sorted(got - exp)}")
    for r, key in zip(recs, keys):
        wave = 1 if any(key[1] in rng for rng in WAVE_SEEDS[1].values()) else 2
        ok = {STUDY_COMMIT} | ({RERUN_COMMIT} if key in RERUN_SEEDS else set())
        if wave == 2:
            ok = {WAVE2_COMMIT} if WAVE2_COMMIT else set()
        if r.get("commit") not in ok:
            refuse(f"{key} built at commit {r.get('commit')!r}, allowed {sorted(ok) or 'none (WAVE2_COMMIT not frozen)'}")
    for plat in {k[0] for k in keys}:
        mine = [r for r in recs if r["platform"] == plat]
        for name in TEXT_INPUTS + ["cube.raw"]:
            if len({r["sha256"][name] for r in mine}) != 1:
                refuse(f"{plat}: input {name} differs between seeds")
        if len({r["materials"]["combined_sha256"] for r in mine}) != 1:
            refuse(f"{plat}: Materials tree differs between seeds")
    if len({r["sha256"]["cube.raw"] for r in recs}) != 1:
        refuse("cube.raw differs between platforms")
    for name in TEXT_INPUTS:
        by_plat = {r["platform"]: r["sha256"][name] for r in recs}
        if len(set(by_plat.values())) != 1:
            print(f"note: text input {name} differs in bytes across platforms (content fixed by commit): {by_plat}")


def welch(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Mean difference a-b, its standard error and Welch-Satterthwaite df."""
    va, vb, na, nb = a.var(ddof=1), b.var(ddof=1), len(a), len(b)
    se2 = va / na + vb / nb
    df = se2 ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)) if se2 > 0 else float("inf")
    return float(a.mean() - b.mean()), math.sqrt(se2), df


def main(path: str, look: int, verdict_path: str) -> None:
    if look not in (1, 2):
        refuse(f"look must be 1 or 2, got {look}")
    frozen = {}
    if look == 1:
        if os.path.exists(verdict_path):
            refuse(f"{verdict_path} exists: look 1 was already analysed")
        continuing = set(COMPARED)
    else:
        if not os.path.exists(verdict_path):
            refuse(f"look 2 needs look 1's verdicts at {verdict_path}")
        frozen = json.load(open(verdict_path))
        continuing = {p for p in COMPARED if frozen.get(p, "").startswith("CONTINUE")}
        if not continuing:
            refuse("no comparison continued at look 1; there is no look 2")
    recs = [json.loads(line) for line in open(path) if line.strip()]
    validate(recs, look, continuing)
    by = defaultdict(list)
    for r in recs:
        by[r["platform"]].append(r)
    print(f"look {look}; seeds per platform: " + ", ".join(f"{p} {len(v)}" for p, v in sorted(by.items())))
    for p, v in sorted(by.items()):
        hosts = defaultdict(int)
        for r in v:
            hosts[r["host"]] += 1
        flags = {ep: sum(1 for r in v if r["metrics"].get(ep.replace("_mm", "_crossings"), 1) > 1) for ep in RANGE}
        print(f"  {p}: hosts {dict(hosts)}; seeds with ambiguous crossings {flags}")
    rows, pvals = [], []
    for plat in COMPARED:
        if plat not in continuing:
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
    for plat in COMPARED:
        if plat not in continuing:
            print(f"VERDICT {plat} vs linux: frozen at look 1: {frozen[plat]}")
    out = {}
    for plat, v in verdict.items():
        if v["pass"] == len(ENDPOINTS):
            res = "EQUIVALENT"
        elif v["noneq"]:
            res = "NON-EQUIVALENT on " + ", ".join(v["noneq"])
        else:
            res = "CONTINUE to wave 2" if look == 1 else "INCONCLUSIVE"
        print(f"VERDICT {plat} vs linux (look {look}): {res}  ({v['pass']}/{len(ENDPOINTS)} endpoints pass)")
        out[plat] = res
    if look == 1:
        with open(verdict_path, "x") as f:
            json.dump(out, f, indent=1, sort_keys=True)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        refuse("usage: platform_study_analyse.py <records.jsonl> <look> <look-1 verdicts.json>")
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
