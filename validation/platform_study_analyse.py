"""Pre-registered analysis for the cross-platform study (validation/win_linux_equivalence_design.md).

Usage: python platform_study_analyse.py <records.jsonl> <look: 1 or 2> <look-1 verdicts.json> [--manifest <p>] [--wave2 <p>]

<records.jsonl> holds one STUDY_RESULT record per line (the JSON after "STUDY_RESULT "), for every seed of every
wave so far. All records are pooled; none is dropped. Needs numpy and scipy.

VALIDATION, before any inference (fails closed; added 2026-10-01 on alden-ec2221c7's review 6821, before any
STUDY_RESULT was read): look must be 1 or 2; study id must be pe1; every record's (platform, seed) must be in the
expected set for that look and every expected seed present exactly once; each record's commit must be the frozen
study commit, or the re-run commit for exactly the two recorded interrupted seeds.

INPUT IDENTITY is checked against a FROZEN EXPECTED MANIFEST (platform_study_manifest.json, next to this file, produced
by platform_study_manifest.py from the study commit; --manifest overrides the path; the run refuses if it is absent or
records a different study commit). Per record: text inputs (plan, cube.mhd, BDL, HU_Density, HU_Material) and the
Materials tree digest must equal the manifest's LF variant on linux and macos, and its LF or CRLF variant on windows
(CRLF is allowed on windows only); cube.raw must equal the manifest's value exactly; the config hash must equal the
manifest's hash for that (platform, seed) exactly. All records of a platform must share one binary sha256, which is
printed. Anything else refuses, naming the record and field. There is no "differs but content is fixed by commit" path.
A non-finite endpoint value makes that endpoint unable to pass.

STOP STATE AND BINDING: look 1 writes <look-1 verdicts.json> (refusing to overwrite) holding "verdicts" (per comparison)
plus sha256 of this analyzer file, sha256 of the manifest file, and a fingerprint of the wave-1 records used (sha256 of
the canonical JSON, sort_keys, one record per line, records sorted by (platform, seed)). Look 2 recomputes all three
from what it is given and refuses on any mismatch, so the analyzer, the manifest or any wave-1 record cannot change
between the looks. The wave-2 commit lives OUTSIDE that binding, in platform_study_wave2.json ({"wave2_commit": null}
until wave 2 is frozen; --wave2 overrides the path); it is read at validation time and wave-2 records are refused
while it is null or the file is absent. A comparison that stopped at
look 1 is reported as frozen and is not recomputed, and does not enter look 2's Holm family.

Per comparison (windows vs linux, macos vs linux), per endpoint, Welch two-sample on per-seed values
(ln dose for CAX), and:
  TOST       95% CI of W-L (alpha = 0.025 one-sided at each look) inside the margin      -> endpoint passes
  NONEQ      simultaneous CI at 1 - 0.025/26 two-sided entirely outside the margin      -> confirmatory non-equivalence
  FLAG       plain 95% CI entirely outside the margin (exploratory only)
  Holm       two-sided Welch p-values, Holm across all 26 endpoint-comparisons at 0.025 (secondary)
Verdict per comparison: EQUIVALENT if all 13 pass; NON-EQUIVALENT if any NONEQ; otherwise CONTINUE (look 1)
or INCONCLUSIVE (look 2). An endpoint with any missing R, or any invalid CAX dose, cannot pass.
"""

import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

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
WAVE_SEEDS = {1: {"windows": range(1001, 1017), "linux": range(2001, 2017), "macos": range(3001, 3017)},
              2: {"windows": range(1017, 1033), "linux": range(2017, 2033), "macos": range(3017, 3033)}}
COMPARED = ("windows", "macos")
MANIFEST_PATH = Path(__file__).with_name("platform_study_manifest.json")
WAVE2_PATH = Path(__file__).with_name("platform_study_wave2.json")  # {"wave2_commit": null | "<40 hex>"}; not bound at look 1
CRLF_PLATFORMS = {"windows"}  # the only platform whose checkout/text-mode writes may legitimately give CRLF
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


def load_wave2_commit(path: Path | str) -> str | None:
    """The frozen wave-2 commit, or None while it is null or the file is absent."""
    path = Path(path)
    if not path.is_file():
        return None
    try:
        c = json.loads(path.read_text())["wave2_commit"]
    except (ValueError, KeyError, TypeError) as e:
        refuse(f"{path} is malformed ({type(e).__name__}: {e})")
    if c is not None and not (isinstance(c, str) and len(c) == 40 and set(c) <= set("0123456789abcdef")):
        refuse(f"{path}: wave2_commit must be null or a 40-hex commit id, got {c!r}")
    return c


def validate(recs: list[dict], look: int, continuing: set[str], manifest: dict, wave2_commit: str | None = None) -> None:
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
        refuse(f"seed set is not the expected one for look {look}: missing {sorted(exp - got)}, "
               f"unexpected {sorted(got - exp)}")
    for r, key in zip(recs, keys):
        wave = 1 if any(key[1] in rng for rng in WAVE_SEEDS[1].values()) else 2
        ok = {STUDY_COMMIT} | ({RERUN_COMMIT} if key in RERUN_SEEDS else set())
        if wave == 2:
            ok = {wave2_commit} if wave2_commit else set()
        if r.get("commit") not in ok:
            refuse(f"{key} built at commit {r.get('commit')!r}, allowed {sorted(ok) or 'none (wave2_commit not frozen)'}")
    check_inputs(recs, manifest)
    check_binaries(recs)


def sha_of(path: Path | str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_manifest(path: Path | str) -> dict:
    """The frozen expected-input manifest; refuses if absent, malformed, or for a different study commit."""
    path = Path(path)
    if not path.is_file():
        refuse(f"manifest {path} is absent")
    try:
        m = json.loads(path.read_text())
        if m["study_commit"] != STUDY_COMMIT:
            refuse(f"manifest is for study commit {m['study_commit']!r}, expected {STUDY_COMMIT!r}")
        for name in TEXT_INPUTS:
            m["text"][name]["LF"], m["text"][name]["CRLF"]
        m["cube.raw"], m["config"], m["materials"]["LF"], m["materials"]["CRLF"]
    except (ValueError, KeyError, TypeError) as e:
        refuse(f"manifest {path} is malformed ({type(e).__name__}: {e})")
    return m


def allowed_variants(plat: str, variants: dict[str, set[str]]) -> set[str]:
    """LF variant(s) for every platform; CRLF variant(s) added on windows only."""
    return set(variants["LF"]) | (set(variants["CRLF"]) if plat in CRLF_PLATFORMS else set())


def check_inputs(recs: list[dict], manifest: dict) -> None:
    """Every record's input hashes must be one of the manifest's approved values for its platform."""
    for r in recs:
        plat, seed = r["platform"], r["seed"]
        who = f"record ({plat}, {seed})"
        sh = r.get("sha256", {})
        for name in TEXT_INPUTS:
            ok = allowed_variants(plat, {v: {manifest["text"][name][v]} for v in ("LF", "CRLF")})
            if sh.get(name) not in ok:
                refuse(f"{who}: field sha256[{name!r}] = {sh.get(name)!r} is not an approved {plat} variant "
                       "of the frozen input")
        if sh.get("cube.raw") != manifest["cube.raw"]:
            refuse(f"{who}: field sha256['cube.raw'] = {sh.get('cube.raw')!r} is not the frozen cube.raw")
        mats = {v: set(manifest["materials"][v].values()) for v in ("LF", "CRLF")}
        got = r.get("materials", {}).get("combined_sha256")
        if got not in allowed_variants(plat, mats):
            refuse(f"{who}: field materials.combined_sha256 = {got!r} is not an approved {plat} Materials digest")
        want = manifest["config"].get(plat, {}).get(str(seed))
        if want is None:
            refuse(f"{who}: the manifest holds no config hash for ({plat}, {seed})")
        if sh.get("config") != want:
            refuse(f"{who}: field sha256['config'] = {sh.get('config')!r} is not the expected cfg.txt hash "
                   f"for ({plat}, {seed})")


def check_binaries(recs: list[dict]) -> None:
    """One binary sha256 per platform; the approved hashes are printed by main."""
    for plat in sorted({r["platform"] for r in recs}):
        hashes = {r.get("sha256", {}).get("binary") for r in recs if r["platform"] == plat}
        if len(hashes) != 1 or None in hashes:
            refuse(f"{plat}: records do not share one binary sha256 "
                   f"({len(hashes)} distinct values: {sorted(map(str, hashes))})")


def fingerprint(recs: list[dict]) -> str:
    """sha256 of the canonical JSON of these records: sort_keys, one record per line, sorted by (platform, seed)."""
    lines = (json.dumps(r, sort_keys=True) for r in sorted(recs, key=lambda r: (r["platform"], r["seed"])))
    return hashlib.sha256("".join(line + "\n" for line in lines).encode()).hexdigest()


def wave1_records(recs: list[dict]) -> list[dict]:
    """The wave-1 subset (by seed) of the supplied records."""
    return [r for r in recs if r["seed"] in {s for rng in WAVE_SEEDS[1].values() for s in rng}]


def check_binding(state: dict, recs: list[dict], manifest_path: Path | str) -> None:
    """Look 2: the analyzer, the manifest and the wave-1 records must be exactly what look 1 used."""
    for key in ("analyzer_sha256", "manifest_sha256", "wave1_fingerprint", "verdicts"):
        if key not in state:
            refuse(f"look-1 state lacks {key!r}")
    if state["analyzer_sha256"] != sha_of(__file__):
        refuse("the analyzer file changed since look 1")
    if state["manifest_sha256"] != sha_of(manifest_path):
        refuse("the manifest file changed since look 1")
    if state["wave1_fingerprint"] != fingerprint(wave1_records(recs)):
        refuse("the wave-1 records differ from the ones look 1 analysed")


def welch(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Mean difference a-b, its standard error and Welch-Satterthwaite df."""
    va, vb, na, nb = a.var(ddof=1), b.var(ddof=1), len(a), len(b)
    se2 = va / na + vb / nb
    df = se2 ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)) if se2 > 0 else float("inf")
    return float(a.mean() - b.mean()), math.sqrt(se2), df


def main(path: str, look: int, verdict_path: str, manifest_path: Path | str | None = None,
         wave2_path: Path | str | None = None) -> None:
    if look not in (1, 2):
        refuse(f"look must be 1 or 2, got {look}")
    manifest_path = manifest_path or MANIFEST_PATH
    manifest = load_manifest(manifest_path)
    frozen, state = {}, {}
    if look == 1:
        if os.path.exists(verdict_path):
            refuse(f"{verdict_path} exists: look 1 was already analysed")
        continuing = set(COMPARED)
    else:
        if not os.path.exists(verdict_path):
            refuse(f"look 2 needs look 1's verdicts at {verdict_path}")
        state = json.load(open(verdict_path))
        frozen = state.get("verdicts", {}) if isinstance(state, dict) else {}
        continuing = {p for p in COMPARED if frozen.get(p, "").startswith("CONTINUE")}
        if not continuing:
            refuse("no comparison continued at look 1; there is no look 2")
    recs = [json.loads(line) for line in open(path) if line.strip()]
    validate(recs, look, continuing, manifest, load_wave2_commit(wave2_path or WAVE2_PATH))
    if look == 2:
        check_binding(state, recs, manifest_path)
    by = defaultdict(list)
    for r in recs:
        by[r["platform"]].append(r)
    print(f"look {look}; seeds per platform: " + ", ".join(f"{p} {len(v)}" for p, v in sorted(by.items())))
    print("approved binary sha256 per platform: "
          + ", ".join(f"{p} {v[0]['sha256']['binary']}" for p, v in sorted(by.items())))
    print(f"analyzer sha256 {sha_of(__file__)}; manifest sha256 {sha_of(manifest_path)}")
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
            json.dump({"verdicts": out, "analyzer_sha256": sha_of(__file__), "manifest_sha256": sha_of(manifest_path),
                       "wave1_fingerprint": fingerprint(wave1_records(recs))}, f, indent=1, sort_keys=True)


if __name__ == "__main__":
    args, opts = sys.argv[1:], {}
    while len(args) >= 2 and args[-2] in ("--manifest", "--wave2"):
        opts[args[-2]], args = args[-1], args[:-2]
    if len(args) != 3:
        refuse("usage: platform_study_analyse.py <records.jsonl> <look> <look-1 verdicts.json> [--manifest <p>] [--wave2 <p>]")
    main(args[0], int(args[1]), args[2], opts.get("--manifest"), opts.get("--wave2"))
