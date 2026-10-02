"""TOPAS opt0 against Portable MCsquare in the lateral halo: the analysis of validation/topas_halo_addendum.md.

DESCRIPTIVE. Estimation only: no margin, no pass or fail, no attribution of cause.

Arms
  TOPAS     opt0, 1e7 histories, the frozen seeds in TOPAS_SEEDS: 8 runs at each of 100, 150 and 200 MeV, made by
            validation/topas/make_run.sh and run_topas.sh under --topas-root.
  MCsquare  the EXISTING same-host Portable runs at acquisition commit 2f9dab40 (the collected, hash-verified tree
            under --mcsquare-root). Their endpoints were read before this design was fixed, and they ran on other
            hosts (Lenovo, HP). All three Portable arms (A-port, B-pgcc, B-picc) are compared with TOPAS separately:
            none is selected or pooled, because every one of them has been read. They share one TOPAS arm, so the
            three comparisons are not independent.

Order of work, so that no endpoint is read from a partial set:
  1. every frozen TOPAS run must have exactly one run directory whose provenance.txt says COMPLETE, with run.txt and
     stage1_base.txt as recorded there and the base identical to this commit's validation/topas/stage1_base.txt.
     If any run is missing, incomplete, duplicated or unexpected, the analysis REFUSES here: no dose file is opened.
  2. each dose.bin is hashed against its provenance, then its endpoints are computed by pencil_endpoints.py with
     the slab depths of its energy.
  3. the MCsquare tree is verified exactly as apples_analyse.py verifies it (frozen population, collection manifest,
     dataset fingerprint), and its endpoint records are read.

Estimates, MCsquare arm against TOPAS, with pointwise 90% and 95% Welch intervals (not simultaneous):
  R80, sigma     difference of run-level means, MCsquare - TOPAS, in mm.
  ring fraction  geometric-mean ratio MCsquare / TOPAS, from the difference of mean logs.
  A ring in which any run of either code scored exactly zero has no log, so no ratio is given: the row carries, per
  code, the number of non-zero runs and the arithmetic mean of the per-run fractions. No pseudocount is added.
  An endpoint that is invalid in any run (failed fit, absent or multiple R80 crossing) is reported as not computed.
The 200 MeV rows carry TD's reference bands for display; the bands are not applied at 100 or 150 MeV.

Usage:
  python topas_halo_compare.py --topas-root DIR --status                 completeness only; opens no dose file
  python topas_halo_compare.py --topas-root DIR --mcsquare-root DIR --json OUT.json --md OUT.md
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import statistics as st
import sys
from dataclasses import dataclass
from pathlib import Path

import apples_analyse as aa
import numpy as np
import pencil_endpoints as pe

ENERGIES = (100, 150, 200)
# A block of its own: TD reserves 910001+ (opt0) and 911001+ (opt4) for the confirmatory arms.
TOPAS_SEEDS = {100: tuple(range(912001, 912009)), 150: tuple(range(912011, 912019)), 200: tuple(range(912021, 912029))}
TOPAS_HISTORIES = 10_000_000
TOPAS_EM = "opt0"
PORTABLE_ARMS = ("A-port", "B-pgcc", "B-picc")
PARTS = ("A", "B")
LEVELS = (0.90, 0.95)
BASE = Path(__file__).resolve().parent / "topas" / "stage1_base.txt"
# The same-host dataset this design reuses: apples_analyse.py's fingerprint of the collected tree at 2f9dab40.
MCSQUARE_FINGERPRINT = "0d5be975b294e0b77860e0ddbce5807caa526a0ca56481df5056effff40033ba"
LABEL = (
    "DESCRIPTIVE. The MCsquare endpoints are the same-host Portable runs at 2f9dab40: they were read before this "
    "design was fixed, and they ran on other hosts than TOPAS. No margin, no pass or fail, no attribution of cause. "
    "Intervals are pointwise. The three Portable arms share one TOPAS arm, so their rows are not independent."
)
# TD, "Reference bands, 200 MeV, arm - TOPAS opt0". Display only.
_BANDS_200 = {"R80": (-0.3, 0.3), "sigma_100": (-0.1, 0.1), "sigma_200": (-0.15, 0.15)}
_RING_BANDS_200 = {(20, 40): (0.90, 1.10), (40, 80): (0.90, 1.10), (80, 200): (0.75, 1.25)}

_RUN_DIR = re.compile(r"^E(\d+)_(opt[04])_seed(\d+)_n(\d+)_th(\d+)_a(\d+)$")
_DOSE_LINE = re.compile(r"^(\d+) sha256: ([0-9a-f]{64})$")
InputError = aa.InputError


@dataclass
class TopasRun:
    energy: int
    seed: int
    path: Path
    threads: int
    provenance: dict[str, str]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_provenance(text: str) -> dict[str, str]:
    """provenance.txt as run_topas.sh writes it: one `key: value` per line. A repeated key is kept as the LAST value
    except `verdict`, where any repeat makes the verdict unusable."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        key, sep, value = line.partition(": ")
        if not sep:
            continue
        if key == "verdict" and "verdict" in out:
            out["verdict"] = "AMBIGUOUS (more than one verdict line)"
            continue
        out[key] = value.strip()
    return out


def _run_txt_problems(run: TopasRun) -> list[str]:
    """run.txt must request exactly what the directory name says (make_run.sh writes both from the same arguments)."""
    text = (run.path / "run.txt").read_text(encoding="utf-8")
    want = (
        f"d:So/Beam/BeamEnergy = {run.energy} MeV",
        f"i:Ts/Seed = {run.seed}",
        f"i:So/Beam/NumberOfHistoriesInRun = {TOPAS_HISTORIES}",
        f"i:Ts/NumberOfThreads = {run.threads}",
    )
    lines = {ln.strip() for ln in text.splitlines()}
    bad = [f"run.txt lacks the line {w!r}" for w in want if w not in lines]
    if '"g4em-standard_opt0"' not in text or "g4em-standard_opt4" in text:
        bad.append("run.txt does not request g4em-standard_opt0 alone")
    return bad


def _input_problems(run: TopasRun, base_sha: str) -> list[str]:
    bad: list[str] = []
    for name in ("run.txt", "stage1_base.txt", "dose.bin", "dose.binheader"):
        if not (run.path / name).is_file():
            bad.append(f"{name} is missing")
    if bad:
        return bad
    bad += _run_txt_problems(run)
    for name in ("run.txt", "stage1_base.txt"):
        if sha256_file(run.path / name) != run.provenance.get(f"{name} sha256"):
            bad.append(f"{name} is not the file recorded in provenance.txt")
    if sha256_file(run.path / "stage1_base.txt") != base_sha:
        bad.append("stage1_base.txt differs from this commit's validation/topas/stage1_base.txt")
    if not _DOSE_LINE.match(run.provenance.get("dose.bin bytes", "")):
        bad.append("provenance.txt has no well-formed 'dose.bin bytes: N sha256: H' line")
    return bad


def discover(root: Path) -> tuple[dict[tuple[int, int], TopasRun], dict[str, list[str]]]:
    """The one COMPLETE run directory of every frozen (energy, seed), without opening any dose file.

    Raises InputError listing EVERY problem if the set is not exactly the frozen one. Returns the runs and the notes
    (attempts that did not complete; entries that are not run directories)."""
    if not root.is_dir():
        msg = f"{root} is not a directory"
        raise InputError(msg)
    base_sha = sha256_file(BASE)
    candidates: dict[tuple[int, int], list[TopasRun]] = {(e, s): [] for e in ENERGIES for s in TOPAS_SEEDS[e]}
    problems: list[str] = []
    notes: dict[str, list[str]] = {"other_attempts": [], "ignored": []}
    for entry in sorted(root.iterdir()):
        m = _RUN_DIR.match(entry.name)
        if m is None or not entry.is_dir():
            notes["ignored"].append(entry.name)
            continue
        energy, em, seed, histories, threads = int(m[1]), m[2], int(m[3]), int(m[4]), int(m[5])
        if em != TOPAS_EM or histories != TOPAS_HISTORIES or (energy, seed) not in candidates:
            problems.append(f"{entry.name}: not a run of the frozen list")
            continue
        prov_path = entry / "provenance.txt"
        prov = parse_provenance(prov_path.read_text(encoding="utf-8")) if prov_path.is_file() else {}
        verdict = prov.get("verdict")
        if verdict == "COMPLETE":
            candidates[(energy, seed)].append(TopasRun(energy, seed, entry, threads, prov))
        else:
            notes["other_attempts"].append(f"{entry.name}: {verdict or 'no verdict (not started, running or interrupted)'}")
    runs: dict[tuple[int, int], TopasRun] = {}
    for (energy, seed), found in candidates.items():
        if len(found) != 1:
            problems.append(f"E{energy} seed {seed}: {len(found)} COMPLETE run directories, need exactly 1")
            continue
        bad = _input_problems(found[0], base_sha)
        problems += [f"{found[0].path.name}: {b}" for b in bad]
        runs[(energy, seed)] = found[0]
    binaries = sorted({r.provenance.get("topas_bin sha256", "absent") for r in runs.values()})
    if len(binaries) > 1:
        problems.append(f"the runs used {len(binaries)} different TOPAS executables: {binaries}")
    if problems:
        msg = (f"TOPAS run set under {root} is not the frozen set of {len(candidates)} complete runs; "
               "no dose file was opened:\n  " + "\n  ".join(problems + notes["other_attempts"]))
        raise InputError(msg)
    return runs, notes


def topas_record(run: TopasRun) -> dict[str, object]:
    """Endpoints of one verified TOPAS run. The dose file must hash as its provenance recorded."""
    m = _DOSE_LINE.match(run.provenance["dose.bin bytes"])
    dose = run.path / "dose.bin"
    if m is None or sha256_file(dose) != m[2]:
        msg = f"{run.path.name}: dose.bin is not the file recorded in provenance.txt"
        raise InputError(msg)
    raw, _dims, spacing = pe.read_topas_bin(dose)
    if not np.allclose(spacing, [1.0, 1.0, 1.0]):
        msg = f"{run.path.name}: expected 1 mm voxels, got {spacing}"
        raise InputError(msg)
    return pe.endpoints(pe.canonical_from_topas(raw), aa.SLAB_DEPTHS[run.energy])


def _endpoints_line(path: Path, ledger: aa.Ledger) -> dict[str, object]:
    lines = [ln for ln in ledger.read_text(path).splitlines() if ln.startswith("ENDPOINTS ")]
    if len(lines) != 1:
        msg = f"{path}: {len(lines)} ENDPOINTS lines, expected exactly 1"
        raise InputError(msg)
    doc = json.loads(lines[0][len("ENDPOINTS "):])
    if not isinstance(doc, dict):
        msg = f"{path}: the ENDPOINTS line is not a JSON object"
        raise InputError(msg)
    return doc


def load_mcsquare(
    root: Path, *, commit: str = aa.ACQUISITION_COMMIT, expect_fingerprint: str | None = MCSQUARE_FINGERPRINT
) -> tuple[dict[tuple[str, int], list[dict[str, object]]], dict[str, object]]:
    """Endpoint records of the three Portable arms per energy, from the collected tree, verified as apples_analyse.py
    verifies it. Refuses a tree without a collection manifest, a Portable cell that is not 8 usable runs, and (unless
    `expect_fingerprint` is None, tests only) any dataset other than the one this design names."""
    frozen = aa.expected_population(commit)
    ledger = aa.Ledger(root)
    if not root.is_dir():
        msg = f"{root} is not a directory"
        raise InputError(msg)
    manifest = aa.read_manifest(root, ledger)
    if manifest is None:
        msg = f"{root}: no collection manifest; the MCsquare arm must be the collected, hash-verified tree"
        raise InputError(msg)
    aa.check_manifest_identity(manifest, commit, PARTS)
    listed = aa.listed_not_established(manifest, frozen)
    data, _issues = aa.load_root(root, PARTS, frozen, ledger, listed, attested=True)
    aa.verify_manifest(root, ledger, manifest, listed)
    out: dict[tuple[str, int], list[dict[str, object]]] = {}
    for arm in PORTABLE_ARMS:
        for energy in ENERGIES:
            runs = sorted((r for r in data[(arm, "P")] if r.energy == energy), key=lambda r: r.seed)
            bad = [f"{r.label}: {'; '.join(r.problems)}" for r in runs if not r.usable]
            if len(runs) != aa.RUNS_PER_CELL or bad:
                msg = f"{arm} at {energy} MeV: {len(runs)} runs, {len(bad)} unusable ({'; '.join(bad)}); need 8 usable"
                raise InputError(msg)
            out[(arm, energy)] = [_endpoints_line(root / arm / "P" / f"s{r.seed}" / "endpoints.json", ledger) for r in runs]
    fingerprint = ledger.fingerprint()
    if expect_fingerprint is not None and fingerprint["sha256"] != expect_fingerprint:
        msg = f"{root}: dataset fingerprint {fingerprint['sha256']} is not the reused dataset {expect_fingerprint}"
        raise InputError(msg)
    return out, fingerprint


def _finite(v: object) -> float | None:
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        return None
    return float(v)


def r80_value(rec: dict[str, object]) -> float | None:
    """R80, or None unless the multiple-crossing flag is present and exactly False (as apples_analyse.py)."""
    return _finite(rec.get("R80")) if rec.get("R80_multiple_crossings") is False else None


def sigma_value(rec: dict[str, object], depth: int) -> float | None:
    s = _finite(rec.get(f"sigma_{depth}"))
    return s if s is not None and aa.SIGMA_BOUNDS_MM[0] <= s <= aa.SIGMA_BOUNDS_MM[1] else None


def ring_value(rec: dict[str, object], depth: int, lo: int, hi: int) -> float | None:
    """A ring fraction in [0, 1]; zero is a valid value here (it is what the zero rows report)."""
    v = _finite(rec.get(f"ring_{depth}_{lo}_{hi}"))
    return v if v is not None and 0.0 <= v <= aa.RING_MAX else None


def _invalid(mc: list[float | None], tp: list[float | None]) -> dict[str, object] | None:
    n_mc, n_tp = sum(v is None for v in mc), sum(v is None for v in tp)
    if n_mc or n_tp:
        return {"status": "not computed",
                "reason": f"invalid in {n_mc} of {len(mc)} MCsquare runs and {n_tp} of {len(tp)} TOPAS runs"}
    return None


def difference_row(mc: list[float | None], tp: list[float | None]) -> dict[str, object]:
    """MCsquare - TOPAS difference of run-level means with pointwise Welch intervals."""
    bad = _invalid(mc, tp)
    if bad is not None:
        return bad
    d, se, df = aa.welch([v for v in mc if v is not None], [v for v in tp if v is not None])
    row: dict[str, object] = {"status": "difference", "estimate": d, "se": se, "df": None if math.isinf(df) else df,
                              "n_mcsquare": len(mc), "n_topas": len(tp)}
    for level in LEVELS:
        row[f"ci{round(level * 100)}"] = list(aa.interval(d, se, df, level))
    return row


def ring_row(mc: list[float | None], tp: list[float | None]) -> dict[str, object]:
    """Geometric-mean ratio MCsquare / TOPAS, or, when any run of either code is zero, the counts and means only."""
    bad = _invalid(mc, tp)
    if bad is not None:
        return bad
    a, b = [v for v in mc if v is not None], [v for v in tp if v is not None]
    if min(a) <= 0.0 or min(b) <= 0.0:
        return {"status": "no ratio (zero runs)",
                "mcsquare": {"n_nonzero": sum(v > 0 for v in a), "n": len(a), "mean_fraction": st.mean(a)},
                "topas": {"n_nonzero": sum(v > 0 for v in b), "n": len(b), "mean_fraction": st.mean(b)}}
    d, se, df = aa.welch([math.log(v) for v in a], [math.log(v) for v in b])
    row: dict[str, object] = {"status": "ratio", "estimate": math.exp(d), "se_log": se,
                              "df": None if math.isinf(df) else df, "n_mcsquare": len(a), "n_topas": len(b)}
    for level in LEVELS:
        row[f"ci{round(level * 100)}"] = [math.exp(x) for x in aa.interval(d, se, df, level)]
    return row


def _row(energy: int, arm: str, endpoint: str, row: dict[str, object], band: tuple[float, float] | None) -> dict[str, object]:
    return {"energy": energy, "arm": arm, "endpoint": endpoint, **row,
            "td_reference_band": list(band) if band is not None and energy == 200 else None}


def compare(
    mcsquare: dict[tuple[str, int], list[dict[str, object]]], topas: dict[int, list[dict[str, object]]]
) -> list[dict[str, object]]:
    """Every row: energy x Portable arm x (R80, then per depth sigma and the five rings)."""
    rows: list[dict[str, object]] = []
    for energy in ENERGIES:
        tp = topas[energy]
        for arm in PORTABLE_ARMS:
            mc = mcsquare[(arm, energy)]
            r80 = difference_row([r80_value(r) for r in mc], [r80_value(r) for r in tp])
            rows.append(_row(energy, arm, "R80", r80, _BANDS_200["R80"]))
            for depth in aa.SLAB_DEPTHS[energy]:
                sigma = difference_row([sigma_value(r, depth) for r in mc], [sigma_value(r, depth) for r in tp])
                rows.append(_row(energy, arm, f"sigma_{depth}", sigma, _BANDS_200.get(f"sigma_{depth}")))
                for lo, hi in aa.RINGS:
                    ring = ring_row([ring_value(r, depth, lo, hi) for r in mc], [ring_value(r, depth, lo, hi) for r in tp])
                    rows.append(_row(energy, arm, f"ring_{depth}_{lo}_{hi}", ring, _RING_BANDS_200.get((lo, hi))))
    return rows


def _num(x: object, digits: int) -> str:
    return f"{x:.{digits}f}" if isinstance(x, float) else "n/a"


def _pair(ci: object, digits: int) -> str:
    return f"[{ci[0]:.{digits}f}, {ci[1]:.{digits}f}]" if isinstance(ci, list) else "n/a"


def markdown(doc: dict[str, object]) -> list[str]:
    md = ["# TOPAS opt0 against Portable MCsquare: lateral halo (descriptive)", "", str(doc["label"]), ""]
    rows: list[dict[str, object]] = doc["rows"]  # type: ignore[assignment]
    for energy in ENERGIES:
        md += [f"## {energy} MeV", "", "| endpoint | arm | estimate | 90% | 95% | TD band (200 MeV only) |",
               "|---|---|---|---|---|---|"]
        zero: list[str] = []
        for r in (r for r in rows if r["energy"] == energy):
            digits = 3 if str(r["endpoint"]).startswith("ring") else 4
            band = _pair(r["td_reference_band"], 2) if r["td_reference_band"] else ""
            if r["status"] in ("difference", "ratio"):
                md.append(f"| {r['endpoint']} ({r['status']}) | {r['arm']} | {_num(r['estimate'], digits)} "
                          f"| {_pair(r['ci90'], digits)} | {_pair(r['ci95'], digits)} | {band} |")
            elif r["status"] == "not computed":
                md.append(f"| {r['endpoint']} | {r['arm']} | not computed: {r['reason']} | | | {band} |")
            else:
                m, t = r["mcsquare"], r["topas"]  # type: ignore[assignment]
                zero.append(f"| {r['endpoint']} | {r['arm']} | {m['n_nonzero']} of {m['n']} | {m['mean_fraction']:.3e} "
                            f"| {t['n_nonzero']} of {t['n']} | {t['mean_fraction']:.3e} |")
        if zero:
            header = "| endpoint | arm | MCsquare non-zero runs | MCsquare mean fraction | TOPAS non-zero runs | TOPAS mean fraction |"
            md += ["", "Rings with a zero run in either code (no ratio):", "", header, "|---|---|---|---|---|---|", *zero]
        md.append("")
    return md


def run(topas_root: Path, mcsquare_root: Path, *, expect_fingerprint: str | None = MCSQUARE_FINGERPRINT) -> dict[str, object]:
    runs, notes = discover(topas_root)  # refuses before any dose file is opened
    mcsquare, fingerprint = load_mcsquare(mcsquare_root, expect_fingerprint=expect_fingerprint)
    topas: dict[int, list[dict[str, object]]] = {e: [] for e in ENERGIES}
    run_docs: list[dict[str, object]] = []
    for (energy, seed), r in sorted(runs.items()):
        topas[energy].append(topas_record(r))
        m = _DOSE_LINE.match(r.provenance["dose.bin bytes"])
        run_docs.append({"dir": r.path.name, "energy": energy, "seed": seed, "threads": r.threads,
                         "host": r.provenance.get("host"), "wall_seconds": r.provenance.get("wall_seconds"),
                         "dose_sha256": m[2] if m else None})
    return {
        "label": LABEL,
        "mcsquare": {"acquisition_commit": aa.ACQUISITION_COMMIT, "dataset_fingerprint": fingerprint,
                     "arms": list(PORTABLE_ARMS), "runs_per_cell": aa.RUNS_PER_CELL},
        "topas": {"em": TOPAS_EM, "histories": TOPAS_HISTORIES,
                  "seeds": {str(e): list(TOPAS_SEEDS[e]) for e in ENERGIES}, "base_sha256": sha256_file(BASE),
                  "topas_bin_sha256": sorted({r.provenance.get("topas_bin sha256") for r in runs.values()}),
                  "runs": run_docs, **notes},
        "levels": list(LEVELS),
        "rows": compare(mcsquare, topas),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--topas-root", type=Path, required=True)
    p.add_argument("--mcsquare-root", type=Path)
    p.add_argument("--status", action="store_true", help="report completeness of the TOPAS set; opens no dose file")
    p.add_argument("--json", type=Path)
    p.add_argument("--md", type=Path)
    a = p.parse_args(argv)
    try:
        if a.status:
            runs, notes = discover(a.topas_root)
            print(f"COMPLETE: all {len(runs)} frozen TOPAS runs are present; no dose file was opened")
            for line in notes["other_attempts"]:
                print(f"other attempt: {line}")
            return 0
        if a.mcsquare_root is None or a.json is None or a.md is None:
            p.error("--mcsquare-root, --json and --md are required unless --status is given")
        doc = run(a.topas_root, a.mcsquare_root)
    except InputError as e:
        print(f"topas_halo_compare.py: REFUSED: {e}", file=sys.stderr)
        return 2
    a.json.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    a.md.write_text("\n".join(markdown(doc)) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
