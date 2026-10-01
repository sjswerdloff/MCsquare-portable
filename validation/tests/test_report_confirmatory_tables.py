"""Tests for validation/report_confirmatory_tables.py: cell formatting, splicing, and the report being up to date."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import report_confirmatory_tables as rt


def row(**over: object) -> dict[str, object]:
    base: dict[str, object] = {
        "endpoint": "P100/R80", "scale": "difference", "outcome": "equivalent", "holm_decision": "equivalent",
        "estimate_reported_scale": 0.0021, "ci95_reported_scale": [-0.0013, 0.0053],
    }
    return {**base, **over}


def test_difference_cell_is_signed_with_a_typographic_minus_and_the_outcome_mark() -> None:
    assert rt.cell(row()) == "+0.0021 [−0.0013, +0.0053] E"


def test_ratio_cell_is_unsigned() -> None:
    r = row(scale="log_ratio", estimate_reported_scale=1.035, ci95_reported_scale=[1.004, 1.067], outcome="inconclusive",
            holm_decision="not shown")
    assert rt.cell(r) == "1.0350 [1.0040, 1.0670] I"


def test_a_value_that_rounds_to_zero_has_no_sign() -> None:
    assert rt.cell(row(estimate_reported_scale=-0.00004)).startswith("0.0000 [")


def test_not_established_has_no_estimate() -> None:
    assert rt.cell(row(outcome="not_established", estimate_reported_scale=None, ci95_reported_scale=None)) == "not established"


def test_equivalent_unadjusted_but_not_after_holm_carries_the_dagger() -> None:
    assert rt.cell(row(holm_decision="not shown")).endswith(" E†")


def test_not_equivalent_is_marked_ne() -> None:
    assert rt.cell(row(outcome="not_equivalent", holm_decision="not shown")).endswith(" NE")


def test_descriptive_cell_has_no_mark() -> None:
    assert rt.cell(row(outcome="descriptive", holm_decision="")) == "+0.0021 [−0.0013, +0.0053]"


def test_an_outcome_with_no_estimate_that_is_not_not_established_is_refused() -> None:
    with pytest.raises(rt.ReportError):
        rt.cell(row(estimate_reported_scale=None))


def test_decimals_follow_the_argument() -> None:
    assert rt.cell(row(), 3) == "+0.002 [−0.001, +0.005] E"


def test_splice_replaces_only_what_is_between_the_markers_and_is_idempotent() -> None:
    text = "before\n<!-- BEGIN GENERATED: t -->\nold\n<!-- END GENERATED: t -->\nafter\n"
    once = rt.splice(text, {"t": ["| a |", "| b |"]})
    assert once == "before\n<!-- BEGIN GENERATED: t -->\n| a |\n| b |\n<!-- END GENERATED: t -->\nafter\n"
    assert rt.splice(once, {"t": ["| a |", "| b |"]}) == once


@pytest.mark.parametrize("text", ["no markers\n", "<!-- BEGIN GENERATED: t -->\n<!-- END GENERATED: t -->\n" * 2])
def test_splice_needs_exactly_one_marker_pair(text: str) -> None:
    with pytest.raises(rt.ReportError):
        rt.splice(text, {"t": ["x"]})


def test_load_refuses_a_partial_confirmatory_contrast(tmp_path: Path) -> None:
    doc = json.loads(rt.ANALYSIS_JSON.read_text(encoding="utf-8"))
    doc["contrasts"][0]["partial"] = True
    path = tmp_path / "a.json"
    path.write_text(json.dumps(doc))
    with pytest.raises(rt.ReportError):
        rt.load(path)


def test_load_refuses_a_missing_contrast(tmp_path: Path) -> None:
    doc = json.loads(rt.ANALYSIS_JSON.read_text(encoding="utf-8"))
    doc["contrasts"] = doc["contrasts"][:2]
    path = tmp_path / "a.json"
    path.write_text(json.dumps(doc))
    with pytest.raises(rt.ReportError):
        rt.load(path)


def test_every_endpoint_of_the_family_appears_in_a_table() -> None:
    contrasts = rt.load(rt.ANALYSIS_JSON)
    blocks = rt.render(contrasts)
    family = list(contrasts[rt.CONFIRMATORY["A"]]["rows"])
    assert len(family) == 52
    assert len(blocks["table-1-r80"]) - 2 == 9  # 3 energies x 3 contrasts
    assert len(blocks["table-3-pencil"]) - 2 == 18  # 6 slabs x 3 contrasts, 6 endpoints each = 108 cells
    assert len(blocks["table-4-field"]) - 2 == 13
    assert len(blocks["table-a1-compiler"]) - 2 == 52
    assert 3 + 6 * 6 + 13 == 52  # R80, (sigma + 5 rings) per slab, case F


def test_the_summary_counts_are_the_analysis_documents() -> None:
    lines = rt.summary_table(rt.load(rt.ANALYSIS_JSON))
    assert lines[2:] == [
        "| A | not established | 46 / 46 | 4 | 2 | 0 |",
        "| B1 | not established | 46 / 46 | 4 | 2 | 0 |",
        "| B2 | not established | 46 / 45 | 4 | 2 | 0 |",
    ]


def test_the_report_is_what_the_committed_analysis_document_renders() -> None:
    assert rt.main(["--check"]) == 0


def test_check_fails_when_a_table_cell_is_edited(tmp_path: Path) -> None:
    edited = tmp_path / "report.md"
    edited.write_text(rt.REPORT.read_text(encoding="utf-8").replace("| A | not established | 46 / 46 |", "| A | holds | 52 / 52 |"))
    assert rt.main(["--check", "--report", str(edited)]) == 1
