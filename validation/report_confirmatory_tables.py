"""Render the confirmatory same-host tables of docs/validation_report_draft.md from the committed analysis JSON.

The input is the document validation/apples_analyse.py wrote for the collection at the acquisition commit
(validation/report_data/apples_analysis_2f9dab40.json). Nothing is recomputed here: every number in a table is a
field of that document, formatted. The tables live in the report between marker comments

    <!-- BEGIN GENERATED: <name> -->  ...  <!-- END GENERATED: <name> -->

and `--write` replaces what is between them; `--check` exits 1 when the report differs from what this script renders
(the test suite runs the same comparison), so a table cannot be edited by hand and stay.

Direction is the analysis's: arm minus reference, i.e. Portable minus upstream (ratios Portable / upstream).

Usage:
    python validation/report_confirmatory_tables.py            # print the blocks
    python validation/report_confirmatory_tables.py --write    # splice them into the report
    python validation/report_confirmatory_tables.py --check    # exit 1 if the report is out of date
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ANALYSIS_JSON = REPO / "validation" / "report_data" / "apples_analysis_2f9dab40.json"
REPORT = REPO / "docs" / "validation_report_draft.md"

# report label -> the analysis's contrast name
CONFIRMATORY = {"A": "A-port vs A-up", "B1": "B-pgcc vs B-up", "B2": "B-picc vs B-up"}
DESCRIPTIVE = "B-pgcc vs B-picc (descriptive)"
SLABS = {100: (40, 60), 150: (80, 125), 200: (100, 200)}
RINGS = ("5_10", "10_20", "20_40", "40_80", "80_200")
OUTCOME_MARK = {"equivalent": "E", "inconclusive": "I", "not_equivalent": "NE"}
# case-F endpoint -> (row label, unit of the estimate, decimals)
F_ROWS = {
    "F/lateral_127_5": ("127 mm depth, 5 mm outside the edge", "points", 3),
    "F/lateral_127_10": ("127 mm depth, 10 mm outside", "points", 3),
    "F/lateral_127_20": ("127 mm depth, 20 mm outside", "points", 3),
    "F/lateral_127_30": ("127 mm depth, 30 mm outside", "points", 3),
    "F/lateral_201_5": ("201 mm depth, 5 mm outside the edge", "points", 3),
    "F/lateral_201_10": ("201 mm depth, 10 mm outside", "points", 3),
    "F/lateral_201_20": ("201 mm depth, 20 mm outside", "points", 3),
    "F/lateral_201_30": ("201 mm depth, 30 mm outside", "points", 3),
    "F/cax_127": ("central-axis dose at 127 mm", "ratio", 4),
    "F/cax_201": ("central-axis dose at 201 mm", "ratio", 4),
    "F/cax_i23": ("central-axis dose at 253 mm", "ratio", 4),
    "F/r80_mm": ("R80", "mm", 3),
    "F/r20_mm": ("R20", "mm", 3),
}
HOLM_NOTE = "†"
MINUS = "−"


class ReportError(Exception):
    """The analysis document or the report is not what this script renders from."""


def load(path: Path) -> dict[str, dict[str, object]]:
    """The analysis document's contrasts by name, each with its endpoint rows keyed by endpoint id."""
    doc = json.loads(path.read_text(encoding="utf-8"))
    out: dict[str, dict[str, object]] = {}
    for c in doc["contrasts"]:
        out[c["name"]] = {**c, "rows": {e["endpoint"]: e for e in c["endpoints"]}}
    missing = [n for n in (*CONFIRMATORY.values(), DESCRIPTIVE) if n not in out]
    if missing:
        msg = f"{path}: contrast(s) missing: {', '.join(missing)}"
        raise ReportError(msg)
    for name in CONFIRMATORY.values():
        if out[name]["partial"] or out[name]["kind"] != "confirmatory":
            msg = f"{path}: {name} is not a complete confirmatory contrast; the report's tables do not cover that"
            raise ReportError(msg)
    return out


def num(x: float, decimals: int, *, signed: bool) -> str:
    """A number with a typographic minus; differences carry an explicit sign, except a value that rounds to zero."""
    if round(x, decimals) == 0:
        return f"{0:.{decimals}f}"
    return f"{x:{'+' if signed else ''}.{decimals}f}".replace("-", MINUS)


def cell(row: dict[str, object], decimals: int = 4) -> str:
    """One endpoint of one contrast: estimate, 95% interval and the TOST outcome mark (none for a descriptive row).

    A not-established endpoint has no estimate. An endpoint equivalent unadjusted but not after Holm carries a dagger.
    """
    if row["outcome"] == "not_established":
        return "not established"
    est, ci = row["estimate_reported_scale"], row["ci95_reported_scale"]
    if est is None or ci is None:
        msg = f"{row['endpoint']}: outcome {row['outcome']} without an estimate"
        raise ReportError(msg)
    signed = row["scale"] == "difference"
    text = f"{num(est, decimals, signed=signed)} [{num(ci[0], decimals, signed=signed)}, {num(ci[1], decimals, signed=signed)}]"
    if row["outcome"] == "descriptive":
        return text
    mark = OUTCOME_MARK[str(row["outcome"])]
    if row["outcome"] == "equivalent" and row["holm_decision"] != "equivalent":
        mark += HOLM_NOTE
    return f"{text} {mark}"


def summary_table(contrasts: dict[str, dict[str, object]]) -> list[str]:
    out = [
        "| contrast | joint claim (all 52 equivalent) | equivalent, unadjusted / Holm | inconclusive | not established | not equivalent |",
        "|---|---|---|---|---|---|",
    ]
    for label, name in CONFIRMATORY.items():
        c = contrasts[name]
        rows, claims = c["rows"].values(), c["claims"]
        count = lambda outcome: sum(r["outcome"] == outcome for r in rows)  # noqa: E731
        joint = "holds" if claims["joint_claim_equivalent_on_all_endpoints"] else "not established"
        out.append(
            f"| {label} | {joint} | {claims['n_equivalent_unadjusted']} / {claims['n_equivalent_holm']} | "
            f"{count('inconclusive')} | {claims['n_not_established']} | {count('not_equivalent')} |"
        )
    return out


def r80_table(contrasts: dict[str, dict[str, object]]) -> list[str]:
    out = ["| energy | contrast | R80 difference (mm) |", "|---|---|---|"]
    for energy in SLABS:
        for i, (label, name) in enumerate(CONFIRMATORY.items()):
            first = f"{energy} MeV" if i == 0 else ""
            out.append(f"| {first} | {label} | {cell(contrasts[name]['rows'][f'P{energy}/R80'])} |")
    return out


def pencil_table(contrasts: dict[str, dict[str, object]]) -> list[str]:
    out = ["| energy, depth | contrast | σ difference (mm) | 5–10 mm | 10–20 mm | 20–40 mm | 40–80 mm | 80–200 mm |",
           "|---|---|---|---|---|---|---|---|"]
    for energy, depths in SLABS.items():
        for depth in depths:
            for i, (label, name) in enumerate(CONFIRMATORY.items()):
                rows = contrasts[name]["rows"]
                first = f"{energy} MeV, {depth} mm" if i == 0 else ""
                cells = [cell(rows[f"P{energy}/sigma_{depth}"])]
                cells += [cell(rows[f"P{energy}/ring_{depth}_{ring}"]) for ring in RINGS]
                out.append(f"| {first} | {label} | " + " | ".join(cells) + " |")
    return out


def field_table(contrasts: dict[str, dict[str, object]]) -> list[str]:
    out = ["| endpoint | scale | A | B1 | B2 |", "|---|---|---|---|---|"]
    for eid, (label, unit, decimals) in F_ROWS.items():
        cells = [cell(contrasts[name]["rows"][eid], decimals) for name in CONFIRMATORY.values()]
        out.append(f"| {label} | {unit} | " + " | ".join(cells) + " |")
    return out


def not_equivalent_table(contrasts: dict[str, dict[str, object]]) -> list[str]:
    """Every endpoint that is not equivalent (unadjusted or after Holm) in at least one confirmatory contrast, shown
    for all three."""
    order = list(next(iter(contrasts.values()))["rows"])
    keep = [
        e for e in order
        if any(contrasts[n]["rows"][e][k] != "equivalent" for n in CONFIRMATORY.values() for k in ("outcome", "holm_decision"))
    ]
    out = ["| endpoint | A | B1 | B2 |", "|---|---|---|---|"]
    out += [f"| {e} | " + " | ".join(cell(contrasts[n]["rows"][e]) for n in CONFIRMATORY.values()) + " |" for e in keep]
    return out


def compiler_table(contrasts: dict[str, dict[str, object]]) -> list[str]:
    """The descriptive contrast (Portable gcc against Portable icc on the HP): every endpoint, no margin, no claim."""
    out = ["| endpoint | Portable gcc − Portable icc, or ratio gcc / icc |", "|---|---|"]
    for eid, row in contrasts[DESCRIPTIVE]["rows"].items():
        decimals = F_ROWS[eid][2] if eid in F_ROWS else 4
        out.append(f"| {eid} | {cell(row, decimals)} |")
    return out


def render(contrasts: dict[str, dict[str, object]]) -> dict[str, list[str]]:
    """Every generated block of the report, by marker name."""
    return {
        "confirmatory-summary": summary_table(contrasts),
        "confirmatory-not-equivalent": not_equivalent_table(contrasts),
        "table-1-r80": r80_table(contrasts),
        "table-3-pencil": pencil_table(contrasts),
        "table-4-field": field_table(contrasts),
        "table-a1-compiler": compiler_table(contrasts),
    }


def splice(report: str, blocks: dict[str, list[str]]) -> str:
    """The report with each block's lines between its markers. Every block must have exactly one marker pair."""
    for name, lines in blocks.items():
        pattern = re.compile(
            rf"(<!-- BEGIN GENERATED: {re.escape(name)} -->\n).*?(<!-- END GENERATED: {re.escape(name)} -->)", re.DOTALL
        )
        if len(pattern.findall(report)) != 1:
            msg = f"the report must hold exactly one marker pair for {name!r}"
            raise ReportError(msg)
        report = pattern.sub(lambda m, body="\n".join(lines): f"{m.group(1)}{body}\n{m.group(2)}", report)
    return report


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--json", type=Path, default=ANALYSIS_JSON, help="the analysis document")
    ap.add_argument("--report", type=Path, default=REPORT)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help="splice the blocks into the report")
    mode.add_argument("--check", action="store_true", help="exit 1 if the report differs from the rendered blocks")
    args = ap.parse_args(argv)
    try:
        blocks = render(load(args.json))
        if not (args.write or args.check):
            for name, lines in blocks.items():
                print(f"## {name}\n\n" + "\n".join(lines) + "\n")
            return 0
        current = args.report.read_text(encoding="utf-8")
        wanted = splice(current, blocks)
    except (ReportError, OSError, KeyError, ValueError) as e:
        print(f"report_confirmatory_tables: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    if args.check:
        if wanted != current:
            print(f"{args.report} is not what {args.json.name} renders; run with --write", file=sys.stderr)
            return 1
        return 0
    args.report.write_text(wanted, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
