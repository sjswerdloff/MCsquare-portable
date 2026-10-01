"""Pre-specified analysis for the same-host comparison programme, parts A and B (docs/apples_to_apples_plan.md,
"Pre-specified analysis" and "Endpoints").

Usage:
    uv run --no-project --with numpy --with scipy python validation/apples_analyse.py <root> [--json PATH] [--parts A,B]

Input layout, <root>/<arm>/<case>/s<seed>/ with arms A-up, A-port (part A) and B-up, B-pgcc, B-picc (part B), cases P
and F. Each seed directory holds run.json plus endpoints.json (case P: one line "ENDPOINTS " + JSON, as
pencil_endpoints.py prints it) or record.json (case F: as written by platform_study_record.py; the 13 endpoints are
its "metrics" object).

Contrasts, each arm minus reference (ratio for rings): A-port vs A-up; B-pgcc vs B-up and B-picc vs B-up (each its own
family); B-pgcc vs B-picc is DESCRIPTIVE only (estimates and intervals, no margin, no equivalence claim).

Per endpoint: estimate, Welch SE, Welch-Satterthwaite df, 90% and 95% intervals (same estimator conventions as
validation/report_tables.py: Welch on run-level values, log scale for ring fractions, difference of mean logs).
TOST at one-sided alpha 0.05 is the 90% interval strictly inside the margin; p_TOST is the larger of the two one-sided
t p-values. Outcomes: equivalent / not_equivalent / inconclusive / not_established. An endpoint is not_established if
any run of either arm is unusable for it (missing record, non-finite value, zero or negative ring fraction, failed
sigma fit, R80 with multiple crossings, transport_status not "ok") or if either arm has fewer than 2 runs. Nothing is
excluded and nothing is imputed.

Claims, kept separate: (1) the JOINT claim of a confirmatory contrast holds only if every endpoint is equivalent at
unadjusted alpha (intersection-union); (2) INDIVIDUAL claims use Holm on p_TOST over the whole family (39 P + 13 F = 52),
with not_established endpoints kept in the family as p = 1 (never rejected). Holm is never applied to the joint claim.

Identity. Arm, case, mode and commit come from run.json and the directory layout <root>/<arm>/<case>/s<seed> only; the
"study" field of record.json is never read, and the report header says so. Design completeness: the plan has 8 runs per
arm, case and energy (case F: 8 at 200 MeV). A contrast whose arm or reference has a cell with a different number of
runs present is reported PARTIAL with the counts, and its joint claim cannot be TRUE.

Exit status: 0 when the analysis ran, whatever the outcomes; 2 when the inputs are unusable (reason on stderr).
"""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics as st
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from platform_study_analyse import CAX as F_CAX
from platform_study_analyse import ENDPOINTS as F_ENDPOINTS
from platform_study_analyse import margin as f_margin
from scipy.stats import t as student_t

ALPHA = 0.05  # one-sided TOST alpha, and the Holm alpha
FAMILY_SIZE = 52  # 39 (P) + 13 (F)

PART_ARMS = {"A": ("A-up", "A-port"), "B": ("B-up", "B-pgcc", "B-picc")}
CASES = ("P", "F")
# (energy MeV) -> slab depths (mm), fixed in the plan
SLAB_DEPTHS = {100: (40, 60), 150: (80, 125), 200: (100, 200)}
RINGS = ((5, 10), (10, 20), (20, 40), (40, 80), (80, 200))
R80_MARGIN = 0.05
SIGMA_MARGIN = 0.02
RING_RATIO = (0.98, 1.02)
WIDE_RING_RATIO = (0.95, 1.05)  # ring 80-200 only
TRANSPORT_OK = "ok"
RUNS_PER_CELL = 8  # per arm x case x energy (case F: one cell, at 200 MeV)
IDENTITY_SOURCE = (
    "Run identity (arm, case, mode, commit) is taken from run.json and the directory layout "
    "<root>/<arm>/<case>/s<seed>; the `study` field of record.json is not read."
)

# (name, arm, reference, confirmatory?)
CONTRASTS_BY_PART = {
    "A": (("A-port vs A-up", "A-port", "A-up", True),),
    "B": (
        ("B-pgcc vs B-up", "B-pgcc", "B-up", True),
        ("B-picc vs B-up", "B-picc", "B-up", True),
        ("B-pgcc vs B-picc (descriptive)", "B-pgcc", "B-picc", False),
    ),
}

_SEED_DIR = re.compile(r"^s(\d+)$")


class InputError(Exception):
    """The inputs are unusable (missing arm, unexpected layout); the analysis cannot run."""


@dataclass(frozen=True)
class Spec:
    """One endpoint: where its value comes from, its scale and its margin on the analysis scale."""

    eid: str  # e.g. "P100/sigma_40", "F/lateral_127_5"
    case: str
    key: str  # key in the run's endpoint record
    log: bool
    lo: float  # margin on the analysis scale (log of the ratio when log)
    hi: float
    energy: int | None = None


@dataclass
class Run:
    """One seed directory. `usable` False means no endpoint of this run is established."""

    label: str
    seed: int
    energy: int | None
    usable: bool
    problems: list[str] = field(default_factory=list)
    mode: str | None = None  # run.json "mode", None if absent
    commit: str | None = None  # run.json "commit", None if absent
    raw: dict[str, float | None] = field(default_factory=dict)  # eid -> value, None = invalid for that endpoint


@dataclass
class Result:
    """One endpoint of one contrast."""

    eid: str
    spec: Spec
    outcome: str
    reason: str = ""
    estimate: float | None = None  # analysis scale
    se: float | None = None
    df: float | None = None
    ci90: tuple[float, float] | None = None
    ci95: tuple[float, float] | None = None
    p_tost: float | None = None
    holm_p: float | None = None
    holm_decision: str = ""


# ----------------------------------------------------------------------------------------------- specification


def p_specs() -> list[Spec]:
    """The 39 case-P endpoints (13 per energy)."""
    out: list[Spec] = []
    for energy, depths in SLAB_DEPTHS.items():
        out.append(Spec(f"P{energy}/R80", "P", "R80", False, -R80_MARGIN, R80_MARGIN, energy))
        for d in depths:
            out.append(Spec(f"P{energy}/sigma_{d}", "P", f"sigma_{d}", False, -SIGMA_MARGIN, SIGMA_MARGIN, energy))
        for d in depths:
            for lo, hi in RINGS:
                ratio = WIDE_RING_RATIO if (lo, hi) == (80, 200) else RING_RATIO
                key = f"ring_{d}_{lo}_{hi}"
                out.append(Spec(f"P{energy}/{key}", "P", key, True, math.log(ratio[0]), math.log(ratio[1]), energy))
    return out


def f_specs() -> list[Spec]:
    """The 13 case-F endpoints with the #31 margins, taken from platform_study_analyse.margin (not copied)."""
    out = []
    for ep in F_ENDPOINTS:
        lo, hi = f_margin(ep)
        out.append(Spec(f"F/{ep}", "F", ep, ep in F_CAX, lo, hi))
    return out


def all_specs() -> list[Spec]:
    specs = p_specs() + f_specs()
    if len(specs) != FAMILY_SIZE or len({s.eid for s in specs}) != FAMILY_SIZE:
        msg = f"endpoint family has {len(specs)} entries, expected {FAMILY_SIZE} distinct"
        raise AssertionError(msg)
    return specs


# ------------------------------------------------------------------------------------------------------ loading


def _num(v: object) -> float | None:
    """A finite real number, else None (booleans are not numbers here)."""
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return None
    return float(v) if math.isfinite(v) else None


def _read_json(path: Path) -> object:
    return json.loads(path.read_text())


def _p_values(rec: dict[str, object], energy: int) -> dict[str, float | None]:
    """Case-P endpoint values of one run for one energy; None where that endpoint is invalid for this run."""
    out: dict[str, float | None] = {}
    r80 = _num(rec.get("R80"))
    # fail closed: the multiple-crossing flag must be present and exactly False
    out[f"P{energy}/R80"] = r80 if rec.get("R80_multiple_crossings") is False else None
    for d in SLAB_DEPTHS[energy]:
        out[f"P{energy}/sigma_{d}"] = _num(rec.get(f"sigma_{d}"))
        for lo, hi in RINGS:
            v = _num(rec.get(f"ring_{d}_{lo}_{hi}"))
            out[f"P{energy}/ring_{d}_{lo}_{hi}"] = v if v is not None and v > 0 else None
    return out


def _f_values(metrics: dict[str, object]) -> dict[str, float | None]:
    """Case-F endpoint values of one run; None where invalid (missing, non-finite, bad dose, ambiguous R80)."""
    out: dict[str, float | None] = {}
    for ep in F_ENDPOINTS:
        v = _num(metrics.get(ep))
        if ep in F_CAX and (metrics.get(f"{ep}_invalid") is True or v is None or v <= 0):
            v = None
        if ep in ("r80_mm", "r20_mm"):  # #31 flags BOTH as ambiguous on more than one distal crossing
            crossings = metrics.get(ep.replace("_mm", "_crossings"))
            if isinstance(crossings, bool) or not isinstance(crossings, int) or crossings != 1:
                v = None
        out[f"F/{ep}"] = v
    return out


def _energy_of(run_json: dict[str, object] | None, rec: dict[str, object] | None) -> int | None:
    """Energy from run.json (energy_mev); failing that, the one energy whose slab depths the endpoint record lists."""
    if run_json is not None:
        e = run_json.get("energy_mev")
        if isinstance(e, int) and not isinstance(e, bool):
            return e
    if rec is not None:
        depths = rec.get("slab_depths_mm")
        if isinstance(depths, list):
            fits = [e for e, pair in SLAB_DEPTHS.items() if all(d in depths for d in pair)]
            if len(fits) == 1:
                return fits[0]
    return None


def _endpoint_record(path: Path, case: str) -> dict[str, object]:
    """The endpoint record of a run: case P's ENDPOINTS line, or case F's record["metrics"]."""
    if case == "P":
        lines = [ln for ln in path.read_text().splitlines() if ln.startswith("ENDPOINTS ")]
        if len(lines) != 1:
            msg = f"{len(lines)} ENDPOINTS lines, expected exactly 1"
            raise ValueError(msg)
        rec = json.loads(lines[0][len("ENDPOINTS ") :])
    else:
        rec = _read_json(path)
        if isinstance(rec, dict):
            rec = rec.get("metrics")
    if not isinstance(rec, dict):
        msg = "endpoint record is not a JSON object"
        raise TypeError(msg)
    return rec


def load_run(arm: str, case: str, run_dir: Path, seed: int) -> Run:
    """Load one seed directory. A damaged run stays in the population as unusable; only a mislabelled one is fatal."""
    label = f"{arm}/{case}/s{seed}"
    problems: list[str] = []
    run_json: dict[str, object] | None = None
    try:
        parsed = _read_json(run_dir / "run.json")
        if not isinstance(parsed, dict):
            msg = "run.json is not a JSON object"
            raise TypeError(msg)
        run_json = parsed
    except (OSError, ValueError, TypeError) as e:
        problems.append(f"run.json unreadable ({type(e).__name__})")
    if run_json is not None:
        for key, want in (("arm", arm), ("case", case)):
            if run_json.get(key) is not None and run_json.get(key) != want:
                msg = f"{label}: run.json says {key} {run_json.get(key)!r}, directory says {want!r}"
                raise InputError(msg)
        if run_json.get("seed") is not None and run_json.get("seed") != seed:
            msg = f"{label}: run.json says seed {run_json.get('seed')!r}, directory says {seed}"
            raise InputError(msg)
        status = run_json.get("transport_status")
        if status != TRANSPORT_OK:
            problems.append(f"transport_status {status!r}")
    rec: dict[str, object] | None = None
    try:
        rec = _endpoint_record(run_dir / ("endpoints.json" if case == "P" else "record.json"), case)
    except (OSError, ValueError, TypeError) as e:
        problems.append(f"endpoint record unusable ({type(e).__name__}: {e})")
    energy: int | None = None
    raw: dict[str, float | None] = {}
    if case == "P":
        energy = _energy_of(run_json, rec)
        if energy not in SLAB_DEPTHS:
            msg = f"{label}: cannot place the run at 100, 150 or 200 MeV (energy {energy!r})"
            raise InputError(msg)
        if rec is not None:
            raw = _p_values(rec, energy)
    elif rec is not None:
        raw = _f_values(rec)
    mode = run_json.get("mode") if run_json is not None else None
    commit = run_json.get("commit") if run_json is not None else None
    return Run(
        label, seed, energy, usable=not problems, problems=problems,
        mode=mode if isinstance(mode, str) else None, commit=commit if isinstance(commit, str) else None, raw=raw,
    )


def load_arm_case(root: Path, arm: str, case: str) -> list[Run]:
    case_dir = root / arm / case
    if not case_dir.is_dir():
        msg = f"missing {case_dir}"
        raise InputError(msg)
    runs: list[Run] = []
    for entry in sorted(case_dir.iterdir()):
        if entry.name.startswith("."):
            continue
        m = _SEED_DIR.match(entry.name)
        if m is None or not entry.is_dir():
            msg = f"unexpected entry {entry} (expected s<seed> directories)"
            raise InputError(msg)
        runs.append(load_run(arm, case, entry, int(m.group(1))))
    if not runs:
        msg = f"no s<seed> directories under {case_dir}"
        raise InputError(msg)
    return runs


def load_root(root: Path, parts: tuple[str, ...]) -> dict[tuple[str, str], list[Run]]:
    if not root.is_dir():
        msg = f"{root} is not a directory"
        raise InputError(msg)
    data: dict[tuple[str, str], list[Run]] = {}
    for part in parts:
        for arm in PART_ARMS[part]:
            if not (root / arm).is_dir():
                msg = f"missing arm directory {root / arm}"
                raise InputError(msg)
            for case in CASES:
                data[(arm, case)] = load_arm_case(root, arm, case)
    return data


# ----------------------------------------------------------------------------------------------- statistics


def welch(a: list[float], b: list[float]) -> tuple[float, float, float]:
    """Return (mean(a) - mean(b), Welch SE, Welch-Satterthwaite df). Needs at least 2 values per arm."""
    if len(a) < 2 or len(b) < 2:
        msg = "Welch needs at least 2 runs per arm"
        raise ValueError(msg)
    va, vb = st.variance(a) / len(a), st.variance(b) / len(b)
    se = math.sqrt(va + vb)
    if se == 0.0:
        return st.mean(a) - st.mean(b), 0.0, math.inf
    df = (va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1))
    return st.mean(a) - st.mean(b), se, df


def interval(d: float, se: float, df: float, level: float) -> tuple[float, float]:
    h = 0.0 if se == 0.0 else float(student_t.ppf(0.5 + level / 2, df)) * se
    return d - h, d + h


def tost_p(d: float, se: float, df: float, lo: float, hi: float) -> float:
    """max of the two one-sided t p-values for H0: diff <= lo and H0: diff >= hi."""
    if se == 0.0:
        return 0.0 if lo < d < hi else 1.0
    return max(float(student_t.sf((d - lo) / se, df)), float(student_t.sf((hi - d) / se, df)))


def classify(ci90: tuple[float, float], lo: float, hi: float) -> str:
    """equivalent: 90% interval inside the margin; not_equivalent: wholly outside; inconclusive: otherwise."""
    if lo < ci90[0] and ci90[1] < hi:
        return "equivalent"
    if ci90[0] >= hi or ci90[1] <= lo:
        return "not_equivalent"
    return "inconclusive"


def holm_adjust(pvalues: list[float | None]) -> list[float | None]:
    """Holm step-down adjusted p-values over the WHOLE family; a None p (not established) stays in the family as 1.

    Returns None for the entries that had no p-value, so they can never be reported as rejected.
    """
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: (1.0 if pvalues[i] is None else pvalues[i], i))  # type: ignore[type-var]
    adjusted: list[float | None] = [None] * m
    running = 0.0
    for rank, i in enumerate(order):
        p = pvalues[i]
        if p is None:
            continue  # sorted last: nothing after it can be rejected either
        running = max(running, min(1.0, (m - rank) * p))
        adjusted[i] = running
    return adjusted


def values_for(runs: list[Run], spec: Spec) -> tuple[list[float] | None, str]:
    """Analysis-scale values of one endpoint over the runs at its energy, or (None, reason)."""
    sel = [r for r in runs if spec.case == "F" or r.energy == spec.energy]
    out: list[float] = []
    bad: list[str] = []
    for r in sel:
        v = r.raw.get(spec.eid) if r.usable else None
        if v is None:
            bad.append(r.label)
        else:
            out.append(math.log(v) if spec.log else v)
    if bad:
        return None, f"unusable in {len(bad)} of {len(sel)} runs: " + ", ".join(bad[:3]) + (" ..." if len(bad) > 3 else "")
    if len(out) < 2:
        return None, f"{len(out)} run(s), need at least 2 per arm"
    return out, ""


def evaluate(spec: Spec, arm: list[Run], ref: list[Run], *, descriptive: bool) -> Result:
    va, why_a = values_for(arm, spec)
    vb, why_b = values_for(ref, spec)
    if va is None or vb is None:
        reason = "; ".join(w for w in (f"arm {why_a}" if why_a else "", f"reference {why_b}" if why_b else "") if w)
        return Result(spec.eid, spec, "not_established", reason)
    d, se, df = welch(va, vb)
    ci90, ci95 = interval(d, se, df, 0.90), interval(d, se, df, 0.95)
    if descriptive:
        return Result(spec.eid, spec, "descriptive", "", d, se, df, ci90, ci95)
    return Result(
        spec.eid, spec, classify(ci90, spec.lo, spec.hi), "", d, se, df, ci90, ci95, tost_p(d, se, df, spec.lo, spec.hi)
    )


def analyse_contrast(
    data: dict[tuple[str, str], list[Run]], arm: str, ref: str, specs: list[Spec], *, confirmatory: bool
) -> dict[str, object]:
    partial = partial_cells(data, (arm, ref))
    results = [evaluate(s, data[(arm, s.case)], data[(ref, s.case)], descriptive=not confirmatory) for s in specs]
    claims: dict[str, object] = {}
    if confirmatory:
        adjusted = holm_adjust([r.p_tost for r in results])
        for r, hp in zip(results, adjusted, strict=True):
            r.holm_p = hp
            if r.outcome == "not_established":
                r.holm_decision = "not_established"
            else:
                r.holm_decision = "equivalent" if hp is not None and hp < ALPHA else "not shown"
        claims = {
            "joint_claim_equivalent_on_all_endpoints": not partial and all(r.outcome == "equivalent" for r in results),
            "n_equivalent_unadjusted": sum(r.outcome == "equivalent" for r in results),
            "n_equivalent_holm": sum(r.holm_decision == "equivalent" for r in results),
            "n_not_established": sum(r.outcome == "not_established" for r in results),
        }
    return {"results": results, "claims": claims, "partial_cells": partial}


# ------------------------------------------------------------------------------------------ design completeness


def cell_counts(runs: list[Run], case: str) -> dict[str, int]:
    """Runs PRESENT (usable or not) per energy cell: 100/150/200 for case P, one "200" cell for case F."""
    keys = [str(e) for e in SLAB_DEPTHS] if case == "P" else ["200"]
    counts = dict.fromkeys(keys, 0)
    for r in runs:
        key = str(r.energy) if case == "P" else "200"
        counts[key] = counts.get(key, 0) + 1
    return counts


def partial_cells(data: dict[tuple[str, str], list[Run]], arms: tuple[str, ...]) -> list[str]:
    """Cells of the given arms (both cases) whose run count differs from RUNS_PER_CELL, as "arm case MeV: n of 8"."""
    out = []
    for arm in arms:
        for case in CASES:
            for energy, n in cell_counts(data[(arm, case)], case).items():
                if n != RUNS_PER_CELL:
                    out.append(f"{arm} {case} {energy} MeV: {n} of {RUNS_PER_CELL}")
    return out


# ----------------------------------------------------------------------------------------------- reporting


def _fmt(spec: Spec, x: float) -> str:
    return f"{math.exp(x):.4f}" if spec.log else f"{x:+.4f}"


def _margin_text(spec: Spec) -> str:
    if spec.log:
        return f"ratio [{math.exp(spec.lo):.3f}, {math.exp(spec.hi):.3f}]"
    return f"[{spec.lo:+.3f}, {spec.hi:+.3f}]"


def _json_num(x: float | None) -> float | None:
    return x if x is not None and math.isfinite(x) else None


def result_json(r: Result) -> dict[str, object]:
    """JSON-safe endpoint record. Estimates for log endpoints are on the log scale; `ratio_*` give the ratio scale.

    A non-finite value (df is infinite when both arms have zero variance) is written as null.
    """
    s = r.spec
    ratio = (lambda x: math.exp(x)) if s.log else (lambda x: x)
    return {
        "endpoint": r.eid,
        "scale": "log_ratio" if s.log else "difference",
        "outcome": r.outcome,
        "reason": r.reason,
        "estimate": _json_num(r.estimate),
        "se": _json_num(r.se),
        "df": _json_num(r.df),
        "ci90": None if r.ci90 is None else [_json_num(x) for x in r.ci90],
        "ci95": None if r.ci95 is None else [_json_num(x) for x in r.ci95],
        "estimate_reported_scale": None if r.estimate is None else _json_num(ratio(r.estimate)),
        "ci90_reported_scale": None if r.ci90 is None else [_json_num(ratio(x)) for x in r.ci90],
        "ci95_reported_scale": None if r.ci95 is None else [_json_num(ratio(x)) for x in r.ci95],
        "margin": None if r.outcome == "descriptive" else [_json_num(s.lo), _json_num(s.hi)],
        "p_tost": _json_num(r.p_tost),
        "holm_adjusted_p": _json_num(r.holm_p),
        "holm_decision": r.holm_decision,
    }


def markdown_table(name: str, confirmatory: bool, analysed: dict[str, object]) -> list[str]:
    results: list[Result] = analysed["results"]  # type: ignore[assignment]
    claims: dict[str, object] = analysed["claims"]  # type: ignore[assignment]
    partial: list[str] = analysed["partial_cells"]  # type: ignore[assignment]
    lines = [f"### {name}", ""]
    if partial:
        lines += [f"**PARTIAL** ({RUNS_PER_CELL} runs per arm, case and energy planned): " + "; ".join(partial) + ".", ""]
    if confirmatory:
        lines += [
            (f"Joint claim (equivalent on all {len(results)} endpoints, unadjusted alpha {ALPHA}): "
            f"**{'TRUE' if claims['joint_claim_equivalent_on_all_endpoints'] else 'FALSE'}**"
            f"{' (cannot be TRUE: PARTIAL)' if partial else ''}; "
            f"equivalent unadjusted {claims['n_equivalent_unadjusted']}, "
            f"equivalent after Holm {claims['n_equivalent_holm']}, not established {claims['n_not_established']}."),
            "",
        ]
    else:
        lines += ["Descriptive only: estimates and intervals; no margin, no equivalence claim.", ""]
    lines += [
        "Ring estimates and intervals are ratios; their SE is the SE of the log ratio.",
        "",
        "| endpoint | estimate | SE | df | 90% | 95% | margin | outcome | Holm-adjusted decision |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        s = r.spec
        if r.estimate is None or r.se is None or r.df is None or r.ci90 is None or r.ci95 is None:
            lines.append(f"| {r.eid} | n/a | n/a | n/a | n/a | n/a | {_margin_text(s)} | {r.outcome} ({r.reason}) | n/a |")
            continue
        if r.outcome == "descriptive":
            holm = "n/a (descriptive)"
        elif r.holm_p is not None:
            holm = f"{r.holm_decision} (p_adj {r.holm_p:.4f})"
        else:
            holm = r.holm_decision
        margin = "n/a" if r.outcome == "descriptive" else _margin_text(s)
        lines.append(
            f"| {r.eid} | {_fmt(s, r.estimate)} | {r.se:.4f} | {r.df:.1f} | "
            f"[{_fmt(s, r.ci90[0])}, {_fmt(s, r.ci90[1])}] | [{_fmt(s, r.ci95[0])}, {_fmt(s, r.ci95[1])}] | "
            f"{margin} | {r.outcome} | {holm} |"
        )
    return [*lines, ""]


def run_analysis(root: Path, parts: tuple[str, ...]) -> tuple[list[str], dict[str, object]]:
    """Return (markdown lines, JSON-able document). Raises InputError when the inputs are unusable."""
    data = load_root(root, parts)
    specs = all_specs()
    md = ["# Same-host comparison: parts " + ", ".join(parts), "", IDENTITY_SOURCE, ""]
    doc: dict[str, object] = {"parts": list(parts), "alpha": ALPHA, "family_size": FAMILY_SIZE,
                           "identity_source": IDENTITY_SOURCE, "runs_per_cell": RUNS_PER_CELL, "contrasts": []}
    md.append("Runs per arm and case (unusable runs stay in the population):")
    md.append("")
    runs_doc: dict[str, object] = {}
    for (arm, case), runs in sorted(data.items()):
        bad = [f"{r.label}: {'; '.join(r.problems)}" for r in runs if not r.usable]
        modes = sorted({r.mode or "unknown" for r in runs})
        commits = sorted({r.commit or "unknown" for r in runs})
        counts = cell_counts(runs, case)
        runs_doc[f"{arm}/{case}"] = {"n": len(runs), "unusable": bad, "modes": modes, "commits": commits, "per_energy": counts}
        md.append(
            f"- {arm}/{case}: {len(runs)} runs ({', '.join(f'{e} MeV: {n}' for e, n in counts.items())}), {len(bad)} unusable; "
            f"mode {'/'.join(modes)}, commit {'/'.join(c[:12] for c in commits)}" + "".join(f"\n  - {b}" for b in bad)
        )
    md.append("")
    doc["runs"] = runs_doc
    contrasts_doc: list[object] = []
    for part in parts:
        for name, arm, ref, confirmatory in CONTRASTS_BY_PART[part]:
            analysed = analyse_contrast(data, arm, ref, specs, confirmatory=confirmatory)
            md += markdown_table(name, confirmatory, analysed)
            results: list[Result] = analysed["results"]  # type: ignore[assignment]
            contrasts_doc.append(
                {
                    "name": name,
                    "arm": arm,
                    "reference": ref,
                    "kind": "confirmatory" if confirmatory else "descriptive",
                    "claims": analysed["claims"],
                    "partial": bool(analysed["partial_cells"]),
                    "partial_cells": analysed["partial_cells"],
                    "endpoints": [result_json(r) for r in results],
                }
            )
    doc["contrasts"] = contrasts_doc
    return md, doc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--json", type=Path, default=None, help="also write the full result here (allow_nan=False)")
    ap.add_argument("--parts", default="A,B", help="comma-separated parts to analyse (default A,B)")
    args = ap.parse_args(argv)
    try:
        parts = tuple(p.strip() for p in args.parts.split(","))
        if not parts or any(p not in PART_ARMS for p in parts) or len(set(parts)) != len(parts):
            msg = f"--parts must be distinct values from {sorted(PART_ARMS)}, got {args.parts!r}"
            raise InputError(msg)
        md, doc = run_analysis(args.root, parts)
        if args.json is not None:
            args.json.write_text(json.dumps(doc, indent=1, sort_keys=True, allow_nan=False) + "\n")
    except (InputError, OSError) as e:
        print(f"apples_analyse: unusable inputs or output: {e}", file=sys.stderr)
        return 2
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    sys.exit(main())
