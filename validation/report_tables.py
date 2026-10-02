"""Recompute the PLANNING tables of docs/validation_report_draft.md (the PRELIMINARY pencil-beam rows of sections 5.1, 5.2 and 5.4) from the committed endpoint records.

Inputs are the JSON-lines files in validation/report_data/, one per arm; each line carries one run's endpoints as written
by validation/pencil_endpoints.py, plus where the run lives and the sha256 of the record or script it came from.

Estimands (as in the merged design, #41): run-level mean difference for R80 and sigma (mm), geometric-mean ratio for
ring energy fractions (difference of mean logs). Uncertainty is the Welch standard error with Welch-Satterthwaite
degrees of freedom; pointwise 90% and 95% intervals are printed. Nothing is excluded: a missing or non-finite endpoint
in any run is reported as UNAVAILABLE for that row, never imputed.

Usage:
    uv run --no-project --with scipy python validation/report_tables.py [report_data_dir]
"""

from __future__ import annotations

import json
import math
import statistics as st
import sys
from pathlib import Path

from scipy.stats import t as student_t

DEPTHS = ("100", "200")
RINGS = ("5_10", "10_20", "20_40", "40_80", "80_200")
LINEAR = ("R80", "sigma_3", "sigma_100", "sigma_200")


def load(path: Path) -> list[dict]:
    """Return the endpoint dicts of one arm, one per run, in file order."""
    return [json.loads(line)["endpoints"] for line in path.read_text().splitlines() if line.strip()]


def welch(a: list[float], b: list[float]) -> tuple[float, float, float]:
    """Return (mean(a) - mean(b), Welch SE, Welch-Satterthwaite df). Needs at least 2 values per arm."""
    if len(a) < 2 or len(b) < 2:
        raise ValueError("Welch needs at least 2 runs per arm")
    va, vb = st.variance(a) / len(a), st.variance(b) / len(b)
    se = math.sqrt(va + vb)
    if se == 0.0:
        return st.mean(a) - st.mean(b), 0.0, float("inf")
    df = (va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1))
    return st.mean(a) - st.mean(b), se, df


def values(runs: list[dict], key: str, log: bool) -> list[float] | None:
    """The key from every run, logged if asked; None if any run lacks a finite positive value."""
    out = []
    for r in runs:
        v = r.get(key)
        if not isinstance(v, (int, float)) or not math.isfinite(v) or (log and v <= 0):
            return None
        out.append(math.log(v) if log else float(v))
    return out


def row(label: str, a: list[dict], b: list[dict], key: str, log: bool) -> str:
    """One table row: estimate, SE, df and 90/95% intervals (ratio scale for logged keys)."""
    va, vb = values(a, key, log), values(b, key, log)
    if va is None or vb is None:
        return f"| {label} | {key} | UNAVAILABLE (missing or non-finite in at least one run) |  |  |  |"
    d, se, df = welch(va, vb)
    ci = {}
    for level in (0.90, 0.95):
        h = 0.0 if se == 0.0 else student_t.ppf(0.5 + level / 2, df) * se
        ci[level] = (d - h, d + h)
    if log:
        fmt = lambda x: f"{math.exp(x):.4f}"
        est = f"ratio {fmt(d)} (SE of log {se:.4f})"
    else:
        fmt = lambda x: f"{x:+.4f}"
        est = f"{d:+.4f} mm (SE {se:.4f})"
    return (
        f"| {label} | {key} | {est} | {df:.1f} | [{fmt(ci[0.90][0])}, {fmt(ci[0.90][1])}] | "
        f"[{fmt(ci[0.95][0])}, {fmt(ci[0.95][1])}] |"
    )


def table(title: str, label: str, a: list[dict], b: list[dict]) -> list[str]:
    """All linear and ring rows for one contrast a - b (or a / b)."""
    out = [f"### {title}", "", "| contrast | endpoint | estimate | Welch df | 90% | 95% |", "|---|---|---|---|---|---|"]
    out += [row(label, a, b, k, log=False) for k in LINEAR]
    out += [row(label, a, b, f"ring_{d}_{r}", log=True) for d in DEPTHS for r in RINGS]
    return out + [""]


def main(argv: list[str]) -> int:
    root = Path(argv[0]) if argv else Path(__file__).parent / "report_data"
    lenovo = load(root / "pencil200_openmcsquare_fix_lenovo.jsonl")
    portable = load(root / "pencil200_portable_studio.jsonl")
    topas = load(root / "pencil200_topas_opt0_studio.jsonl")
    print(f"runs: openmcsquare+fix Lenovo {len(lenovo)}, portable Studio {len(portable)}, TOPAS opt0 Studio {len(topas)}\n")
    lines = table("Section 5: OpenMCsquare + fix (Lenovo) vs portable (Studio)", "upstream - portable", lenovo, portable)
    lines += table("Section 6: portable (Studio) vs TOPAS opt0 (Studio)", "portable - TOPAS", portable, topas)
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
