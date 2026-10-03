"""Render the confirmatory same-host tables of docs/validation_report_draft.md from the committed analysis JSON.

The input is the document validation/apples_analyse.py wrote for the collection at the acquisition commit
(validation/report_data/apples_analysis_2f9dab40.json), and, for the descriptive tables (1b, 2, 5), the document
validation/report_dose_descriptives.py wrote from the Dose files (validation/report_data/apples_descriptive_2f9dab40.json).
Table 6 (TOPAS against Portable, descriptive) is rendered from the document validation/topas_halo_compare.py wrote
(validation/report_data/halo_868d1910.json); its direction is MCsquare minus TOPAS (ratios MCsquare / TOPAS).
Nothing is recomputed here: every number in a table is a field of one of those documents, formatted. The tables live in the report between marker comments

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
DESCRIPTIVE_JSON = REPO / "validation" / "report_data" / "apples_descriptive_2f9dab40.json"
DESCRIPTIVE_SCHEMA = "apples_descriptive/1"
FRACTIONS_JSON = REPO / "validation" / "report_data" / "apples_ring_fractions_2f9dab40.json"
FRACTIONS_SCHEMA = "apples_ring_fractions/1"
TOPAS_HALO_JSON = REPO / "validation" / "report_data" / "halo_868d1910.json"
PART_E_JSON = REPO / "validation" / "report_data" / "field_e_analysis_26a93235.json"
REPORT = REPO / "docs" / "validation_report_draft.md"

# report label -> the analysis's contrast name
CONFIRMATORY = {"A": "A-port vs A-up", "B1": "B-pgcc vs B-up", "B2": "B-picc vs B-up"}
DESCRIPTIVE = "B-pgcc vs B-picc (descriptive)"
SLABS = {100: (40, 60), 150: (80, 125), 200: (100, 200)}
RINGS = ("5_10", "10_20", "20_40", "40_80", "80_200")
OUTCOME_MARK = {"equivalent": "E", "inconclusive": "I", "not_equivalent": "NE"}
TOPAS_ARMS = ("A-port", "B-pgcc", "B-picc")
TOPAS_NO_RATIO = "no ratio (zero runs)"
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
# descriptive tables: (metric key, column heading) of the range metrics
RANGE_COLUMNS = (("R90", "R90 difference (mm)"), ("R20", "R20 difference (mm)"), ("falloff", "distal fall-off, R20 − R80, difference (mm)"))
# pencil-beam gamma columns: (dimension, gamma spec key, heading); case-F columns: (gamma spec key, heading)
PENCIL_GAMMA_COLUMNS = (
    ("idd", "g2_2_c10", "IDD, 2%/2 mm"),
    ("idd", "g1_1_c10", "IDD, 1%/1 mm"),
    ("3d", "g2_2_c10", "3D, 2%/2 mm"),
    ("3d", "g1_1_c10", "3D, 1%/1 mm"),
)
FIELD_GAMMA_COLUMNS = (
    ("g2_2_c10", "2%/2 mm, 10% cutoff"),
    ("g1_1_c10", "1%/1 mm, 10% cutoff"),
    ("g2_2_c1", "low-dose: 2%/2 mm, 1% cutoff"),
    ("l2_2_c1", "local 2%/2 mm, 1% cutoff"),
)
# the upstream arm whose split halves are the noise control of each contrast
CONTROL_ARM = {"A": "A-up", "B1": "B-up", "B2": "B-up"}
PENCIL_PASS_DECIMALS = 3  # one failing point of the largest pencil-beam analysis (55 838 points) is 0.002 percentage points
PART_E_ENDPOINTS = 34  # 17 per energy at 100 and 150 MeV (docs/field_100_150_plan.md)
FIELD_PASS_DECIMALS = 4  # one failing point of the case-F analyses (about 7.7e5 to 1.0e6 points) is 0.0001 percentage points


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


def load_descriptive(path: Path) -> dict[str, object]:
    """The descriptive document, refused unless it is the expected schema and every binding control in it passed."""
    doc = json.loads(path.read_text(encoding="utf-8"))
    if doc.get("schema") != DESCRIPTIVE_SCHEMA or doc.get("descriptive_only") is not True:
        msg = f"{path}: not a {DESCRIPTIVE_SCHEMA} descriptive document"
        raise ReportError(msg)
    controls = doc["controls"]
    failed = [k for k in ("case_P_r80", "case_F_metrics", "dose_file_sha256") if controls.get(k, {}).get("passed") is not True]
    if failed:
        msg = f"{path}: binding control(s) did not pass: {', '.join(failed)}"
        raise ReportError(msg)
    return doc


def range_cell(entry: dict[str, object]) -> str:
    """One descriptive range-metric cell: estimate [95% interval] by the shared ``cell`` helper, or ``unavailable``.

    A metric is unavailable when some run has no single distal crossing; the descriptive document then carries no
    estimate, and none is invented.
    """
    if not entry["available"]:
        return "unavailable"
    ci = entry["ci95"]
    if ci is None:
        return f"{num(entry['estimate'], 4, signed=True)} [interval undefined: zero variance]"
    row = {
        "endpoint": "descriptive", "scale": "difference", "outcome": "descriptive", "holm_decision": "",
        "estimate_reported_scale": entry["estimate"], "ci95_reported_scale": ci,
    }  # fmt: skip
    return cell(row)


def pass_rate(entry: dict[str, object], decimals: int) -> str:
    """A pass rate in percent with its evaluated-point count; never a rate without the count.

    Raises:
        ReportError: If the rounded rate would read 100 while some evaluated point failed.
    """
    rate = entry["pass_rate_percent"]
    n = entry["n_evaluated"]
    if rate is None:
        return f"no evaluated points (n = {n})"
    text = f"{rate:.{decimals}f}"
    if float(text) >= 100 and entry["n_fail"] > 0:
        msg = f"a pass rate of {rate!r}% would print as 100 with {entry['n_fail']} failing point(s)"
        raise ReportError(msg)
    missing = entry["n_not_evaluated_in_region"]
    return f"{text} (n = {n}{f', {missing} not evaluated' if missing else ''})"


def gamma_cell(entry: dict[str, object], control: dict[str, object], decimals: int) -> str:
    """Pass rate of Portable against upstream, then the split-half noise control of the upstream arm."""
    return f"{pass_rate(entry, decimals)}; control {pass_rate(control, decimals)}"


def range_descriptive_table(desc: dict[str, object]) -> list[str]:
    """Table 1b: R90, R20 and fall-off differences (Portable − upstream) per energy and contrast; no outcome mark."""
    contrasts = desc["range_metrics"]["contrasts"]  # type: ignore[index]
    out = ["| energy | contrast | " + " | ".join(h for _, h in RANGE_COLUMNS) + " |", "|---|---|" + "---|" * len(RANGE_COLUMNS)]
    for energy in SLABS:
        for i, label in enumerate(CONFIRMATORY):
            first = f"{energy} MeV" if i == 0 else ""
            cells = [range_cell(contrasts[label][str(energy)][metric]) for metric, _ in RANGE_COLUMNS]
            out.append(f"| {first} | {label} | " + " | ".join(cells) + " |")
    return out


def gamma_pencil_table(desc: dict[str, object]) -> list[str]:
    """Table 2: pencil-beam gamma pass rates (IDD and 3D, two criteria) with their noise controls and point counts."""
    gamma = desc["gamma_pencil"]
    out = ["| energy | contrast | " + " | ".join(h for _, _, h in PENCIL_GAMMA_COLUMNS) + " |", "|---|---|" + "---|" * len(PENCIL_GAMMA_COLUMNS)]
    for energy in SLABS:
        for i, label in enumerate(CONFIRMATORY):
            first = f"{energy} MeV" if i == 0 else ""
            cells = []
            for dim, spec, _ in PENCIL_GAMMA_COLUMNS:
                block = gamma[str(energy)][dim]
                cells.append(
                    gamma_cell(block["contrasts"][label][spec], block["noise_control"][CONTROL_ARM[label]][spec], PENCIL_PASS_DECIMALS)
                )
            out.append(f"| {first} | {label} | " + " | ".join(cells) + " |")
    return out


def gamma_field_table(desc: dict[str, object]) -> list[str]:
    """Table 5: broad-field gamma pass rates (four analyses) with their noise controls and point counts."""
    gamma = desc["gamma_field"]
    out = ["| contrast | " + " | ".join(h for _, h in FIELD_GAMMA_COLUMNS) + " |", "|---|" + "---|" * len(FIELD_GAMMA_COLUMNS)]
    for label in CONFIRMATORY:
        cells = [
            gamma_cell(gamma["contrasts"][label][spec], gamma["noise_control"][CONTROL_ARM[label]][spec], FIELD_PASS_DECIMALS)
            for spec, _ in FIELD_GAMMA_COLUMNS
        ]
        out.append(f"| {label} | " + " | ".join(cells) + " |")
    return out


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


def load_fractions(path: Path) -> dict[str, float]:
    """The upstream arms' mean annulus energy fractions by endpoint id (validation/report_ring_fractions.py)."""
    doc = json.loads(path.read_text(encoding="utf-8"))
    if doc.get("schema") != FRACTIONS_SCHEMA:
        msg = f"{path}: schema {doc.get('schema')!r}, expected {FRACTIONS_SCHEMA!r}"
        raise ReportError(msg)
    return {str(k): float(v) for k, v in doc["mean_fraction"].items()}


def percent(x: float, decimals: int) -> str:
    return f"{100 * x:.{decimals}f}%"


def halo_bound_table(contrasts: dict[str, dict[str, object]], fractions: dict[str, float]) -> list[str]:
    """For each annulus endpoint not equivalent in some contrast: the share of the slab's energy it holds, the largest
    departure from 1 inside any of the three 95% intervals, and their product.

    The product is a plug-in illustration of scale, not a bound on dose: the fractions are normalised (they do not
    constrain the slab's absolute energy), the share's own uncertainty is not propagated, and an all-zero annulus has
    no interval (stated as such)."""
    out = [
        "| endpoint | share of the slab's energy (upstream mean) | largest \\|ratio − 1\\| within the three 95% intervals | product (illustrative) |",
        "|---|---|---|---|",
    ]
    order = list(next(iter(contrasts.values()))["rows"])
    for eid in order:
        rows = [contrasts[n]["rows"][eid] for n in CONFIRMATORY.values()]
        if "/ring_" not in eid or all(r["outcome"] == "equivalent" for r in rows):
            continue
        share = fractions[eid]
        bounds = [abs(b - 1) for r in rows if r["ci95_reported_scale"] is not None for b in r["ci95_reported_scale"]]
        if share == 0 or not bounds:
            if share != 0 or bounds:
                msg = f"{eid}: a zero share and an interval, or a non-zero share and no interval"
                raise ReportError(msg)
            out.append(f"| {eid} | no energy scored in any run | no interval | — |")
            continue
        worst = max(bounds)
        out.append(f"| {eid} | {percent(share, 4)} | {percent(worst, 1)} | {percent(share * worst, 4)} |")
    return out


def topas_halo_endpoints(energy: int) -> list[tuple[str, str]]:
    """The endpoints of one energy in table order, each with its label: R80, then sigma and the annuli per slab."""
    out = [("R80", "R80 (mm)")]
    for depth in SLABS[energy]:
        out.append((f"sigma_{depth}", f"σ at {depth} mm (mm)"))
        out += [(f"ring_{depth}_{r}", f"{r.replace('_', '–')} mm annulus at {depth} mm") for r in RINGS]
    return out


def load_topas_halo(path: Path) -> dict[tuple[int, str, str], dict[str, object]]:
    """The rows of the TOPAS halo document (validation/topas_halo_compare.py) by (energy, endpoint, arm).

    Refused unless the document holds exactly the rows the table prints, each once: a missing row would leave a
    cell that cannot be rendered, and a row outside the table would be a result the report does not show.
    """
    doc = json.loads(path.read_text(encoding="utf-8"))
    if tuple(doc["mcsquare"]["arms"]) != TOPAS_ARMS:
        msg = f"{path}: arms are {doc['mcsquare']['arms']}, not {list(TOPAS_ARMS)}"
        raise ReportError(msg)
    rows: dict[tuple[int, str, str], dict[str, object]] = {}
    for r in doc["rows"]:
        key = (int(r["energy"]), str(r["endpoint"]), str(r["arm"]))
        if key in rows:
            msg = f"{path}: {key} appears more than once"
            raise ReportError(msg)
        rows[key] = r
    wanted = {(e, ep, a) for e in SLABS for ep, _ in topas_halo_endpoints(e) for a in TOPAS_ARMS}
    if set(rows) != wanted:
        missing, extra = sorted(wanted - set(rows)), sorted(set(rows) - wanted)
        msg = f"{path}: rows are not the table's: missing {missing[:3]}, not in the table {extra[:3]}"
        raise ReportError(msg)
    return rows


def topas_halo_cell(row: dict[str, object]) -> str:
    """One row of the TOPAS halo document: estimate and 95% interval, or why there is none.

    A ring that is zero in some run has no ratio: the cell gives the non-zero run counts of both codes. A row whose
    precision was not estimated (zero sample variance in both codes) keeps its estimate and says so; it is never
    printed with an interval.
    """
    if row["status"] == TOPAS_NO_RATIO:
        m, t = row["mcsquare"], row["topas"]
        return f"no ratio: non-zero in {m['n_nonzero']} of {m['n']} MCsquare and {t['n_nonzero']} of {t['n']} TOPAS runs"
    if row["status"] not in ("difference", "ratio"):
        msg = f"{row['endpoint']}: unknown status {row['status']!r}"
        raise ReportError(msg)
    signed = row["status"] == "difference"
    decimals = 4 if signed else 3
    est = num(float(row["estimate"]), decimals, signed=signed)
    if row["uncertainty"] != "welch":
        if row.get("ci95") is not None:
            msg = f"{row['endpoint']}: an interval beside {row['uncertainty']!r}"
            raise ReportError(msg)
        return f"{est} ({row['uncertainty']})"
    lo, hi = row["ci95"]
    return f"{est} [{num(lo, decimals, signed=signed)}, {num(hi, decimals, signed=signed)}]"


def topas_halo_table(rows: dict[tuple[int, str, str], dict[str, object]]) -> list[str]:
    out = [f"| energy | endpoint | {' | '.join(TOPAS_ARMS)} | TD band |", "|---|---|" + "---|" * (len(TOPAS_ARMS) + 1)]
    for energy in SLABS:
        for i, (endpoint, label) in enumerate(topas_halo_endpoints(energy)):
            cells = [topas_halo_cell(rows[(energy, endpoint, arm)]) for arm in TOPAS_ARMS]
            bands = {json.dumps(rows[(energy, endpoint, arm)]["td_reference_band"]) for arm in TOPAS_ARMS}
            if len(bands) != 1:
                msg = f"{energy} MeV {endpoint}: the arms carry different TD bands"
                raise ReportError(msg)
            band = json.loads(bands.pop())
            signed = endpoint == "R80" or endpoint.startswith("sigma")
            shown = "" if band is None else f"[{num(band[0], 2, signed=signed)}, {num(band[1], 2, signed=signed)}]"
            first = f"{energy} MeV" if i == 0 else ""
            out.append(f"| {first} | {label} | {' | '.join(cells)} | {shown} |")
    return out


def load_part_e(path: Path) -> dict[str, dict[str, object]]:
    """The part E analysis document's three confirmatory contrasts by name, rows keyed by endpoint id.

    Refuses a document in which a confirmatory contrast is missing, partial, has its claims withheld or lacks rows.
    """
    doc = json.loads(path.read_text(encoding="utf-8"))
    out: dict[str, dict[str, object]] = {}
    for c in doc["contrasts"]:
        if c["confirmatory"]:
            out[f"{c['arm']} vs {c['reference']}"] = {**c, "rows": {r["endpoint"]: r for r in c["rows"]}}
    missing = [n for n in CONFIRMATORY.values() if n not in out]
    if missing:
        msg = f"{path}: confirmatory contrast(s) missing: {', '.join(missing)}"
        raise ReportError(msg)
    for name in CONFIRMATORY.values():
        c = out[name]
        if c["claims_withheld"] or c["partial_reasons"] or len(c["rows"]) != PART_E_ENDPOINTS:
            msg = f"{path}: {name} is partial or withheld; the report's part E tables do not cover that"
            raise ReportError(msg)
    if [list(out[n]["rows"]) for n in CONFIRMATORY.values()].count(list(out[CONFIRMATORY["A"]]["rows"])) != 3:
        msg = f"{path}: the confirmatory contrasts do not carry the same endpoints in the same order"
        raise ReportError(msg)
    return out


def part_e_cell(row: dict[str, object]) -> str:
    """One part E endpoint: estimate, 95% interval and the TOST outcome mark, as in Tables 3 and 4."""
    return cell({**row, "estimate_reported_scale": row.get("estimate"), "ci95_reported_scale": row.get("ci95")})


def part_e_summary_table(part_e: dict[str, dict[str, object]]) -> list[str]:
    out = [
        "| contrast | joint claim (all 34 equivalent) | secondary joint claim | equivalent, unadjusted | equivalent, Holm |",
        "|---|---|---|---|---|",
    ]
    for key, name in CONFIRMATORY.items():
        cl = part_e[name]["claims"]
        joint = "established" if cl["joint_claim_equivalent_on_all_endpoints"] else "not established"
        second = "established" if cl["secondary_joint_claim_holm_excluding_all_zero_rows"] else "not established"
        n = len(part_e[name]["rows"])
        out.append(f"| {key}: {name} | {joint} | {second} | {cl['n_equivalent_unadjusted']} of {n} | "
                   f"{cl['n_equivalent_holm']} of {n} |")
    return out


def part_e_endpoint_table(part_e: dict[str, dict[str, object]]) -> list[str]:
    out = ["| endpoint | scale | A | B1 | B2 |", "|---|---|---|---|---|"]
    for eid, row in part_e[CONFIRMATORY["A"]]["rows"].items():
        cells = [part_e_cell(part_e[name]["rows"][eid]) for name in CONFIRMATORY.values()]
        out.append(f"| {eid} | {row['scale']} | " + " | ".join(cells) + " |")
    return out


def render(
    contrasts: dict[str, dict[str, object]],
    descriptive: dict[str, object] | None = None,
    fractions: dict[str, float] | None = None,
    topas_halo: dict[tuple[int, str, str], dict[str, object]] | None = None,
    part_e: dict[str, dict[str, object]] | None = None,
) -> dict[str, list[str]]:
    """Every generated block of the report, by marker name.

    The descriptive blocks (Tables 1b, 2 and 5) are rendered when the descriptive document is given, the
    far-halo bound when the annulus fractions are, and Table 6 when the TOPAS halo document is.
    """
    blocks = {
        "confirmatory-summary": summary_table(contrasts),
        "confirmatory-not-equivalent": not_equivalent_table(contrasts),
        "table-1-r80": r80_table(contrasts),
        "table-3-pencil": pencil_table(contrasts),
        "table-4-field": field_table(contrasts),
        "table-a1-compiler": compiler_table(contrasts),
    }
    if descriptive is not None:
        blocks |= {
            "table-1b-range-descriptive": range_descriptive_table(descriptive),
            "table-2-gamma-pencil": gamma_pencil_table(descriptive),
            "table-5-gamma-field": gamma_field_table(descriptive),
        }
    if fractions is not None:
        blocks["halo-bound"] = halo_bound_table(contrasts, fractions)
    if topas_halo is not None:
        blocks["table-6-topas-halo"] = topas_halo_table(topas_halo)
    if part_e is not None:
        blocks["table-7-part-e-summary"] = part_e_summary_table(part_e)
        blocks["table-8-part-e-endpoints"] = part_e_endpoint_table(part_e)
    return blocks


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


def stale_blocks(report: str, blocks: dict[str, list[str]]) -> list[str]:
    """The names of the blocks whose text in the report is not what was rendered."""
    return [name for name, lines in blocks.items() if splice(report, {name: lines}) != report]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--json", type=Path, default=ANALYSIS_JSON, help="the analysis document")
    ap.add_argument("--descriptive-json", type=Path, default=DESCRIPTIVE_JSON, help="the descriptive document")
    ap.add_argument("--topas-halo-json", type=Path, default=TOPAS_HALO_JSON, help="the TOPAS halo document")
    ap.add_argument("--part-e-json", type=Path, default=PART_E_JSON, help="the part E analysis document")
    ap.add_argument("--report", type=Path, default=REPORT)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help="splice the blocks into the report")
    mode.add_argument("--check", action="store_true", help="exit 1 if the report differs from the rendered blocks")
    args = ap.parse_args(argv)
    try:
        blocks = render(
            load(args.json),
            load_descriptive(args.descriptive_json),
            load_fractions(FRACTIONS_JSON),
            load_topas_halo(args.topas_halo_json),
            load_part_e(args.part_e_json),
        )
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
            stale = ", ".join(stale_blocks(current, blocks))
            print(f"{args.report} is not what the committed documents render (blocks: {stale}); run with --write", file=sys.stderr)
            return 1
        return 0
    args.report.write_text(wanted, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
