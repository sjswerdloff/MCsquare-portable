"""Tests for validation/report_confirmatory_tables.py: cell formatting, splicing, and the report being up to date."""

from __future__ import annotations

import json
import re
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


# ---------------------------------------------------------------------------------------------------------------------
# descriptive blocks: Table 1b (range), Table 2 (pencil-beam gamma), Table 5 (broad-field gamma)
# ---------------------------------------------------------------------------------------------------------------------


def gamma_entry(**over: object) -> dict[str, object]:
    base: dict[str, object] = {"pass_rate_percent": 99.5, "n_evaluated": 1000, "n_fail": 5, "n_not_evaluated_in_region": 0}
    return {**base, **over}


def test_range_cell_uses_the_shared_formatting_and_carries_no_outcome_mark() -> None:
    entry = {"available": True, "estimate": 0.0021, "ci95": [-0.0013, 0.0053]}
    assert rt.range_cell(entry) == "+0.0021 [−0.0013, +0.0053]"
    assert rt.range_cell({"available": True, "estimate": -0.00004, "ci95": [-0.001, 0.001]}).startswith("0.0000 [")


def test_an_unavailable_range_metric_is_printed_as_unavailable_and_never_as_a_number() -> None:
    assert rt.range_cell({"available": False, "reason": "x"}) == "unavailable"
    assert "undefined" in rt.range_cell({"available": True, "estimate": 0.5, "ci95": None})


def test_a_pass_rate_always_carries_its_point_count() -> None:
    assert rt.pass_rate(gamma_entry(), 3) == "99.500 (n = 1000)"
    assert rt.pass_rate(gamma_entry(n_not_evaluated_in_region=4), 3) == "99.500 (n = 1000, 4 not evaluated)"
    assert rt.pass_rate(gamma_entry(pass_rate_percent=None, n_evaluated=0), 3) == "no evaluated points (n = 0)"


def test_a_pass_rate_that_would_round_to_100_while_points_failed_is_refused() -> None:
    with pytest.raises(rt.ReportError):
        rt.pass_rate(gamma_entry(pass_rate_percent=99.99996, n_fail=3), 4)
    assert rt.pass_rate(gamma_entry(pass_rate_percent=100.0, n_fail=0), 3) == "100.000 (n = 1000)"
    assert rt.pass_rate(gamma_entry(pass_rate_percent=99.9991, n_fail=7), 4) == "99.9991 (n = 1000)"


def test_a_gamma_cell_shows_the_contrast_and_then_its_noise_control() -> None:
    cell_text = rt.gamma_cell(gamma_entry(), gamma_entry(pass_rate_percent=98.0, n_evaluated=990, n_fail=20), 3)
    assert cell_text == "99.500 (n = 1000); control 98.000 (n = 990)"


def committed_descriptive() -> dict[str, object]:
    return rt.load_descriptive(rt.DESCRIPTIVE_JSON)


def test_descriptive_load_refuses_a_failed_binding_control_and_another_schema(tmp_path: Path) -> None:
    doc = json.loads(rt.DESCRIPTIVE_JSON.read_text(encoding="utf-8"))
    for control in ("case_P_r80", "case_F_metrics", "dose_file_sha256"):
        broken = json.loads(json.dumps(doc))
        broken["controls"][control]["passed"] = False
        path = tmp_path / f"{control}.json"
        path.write_text(json.dumps(broken))
        with pytest.raises(rt.ReportError):
            rt.load_descriptive(path)
    wrong = {**doc, "schema": "something/else"}
    path = tmp_path / "schema.json"
    path.write_text(json.dumps(wrong))
    with pytest.raises(rt.ReportError):
        rt.load_descriptive(path)


def test_the_descriptive_tables_have_the_expected_shape() -> None:
    desc = committed_descriptive()
    rng, pencil, field = rt.range_descriptive_table(desc), rt.gamma_pencil_table(desc), rt.gamma_field_table(desc)
    assert len(rng) - 2 == 9  # 3 energies x 3 contrasts
    assert len(pencil) - 2 == 9
    assert len(field) - 2 == 3  # 3 contrasts
    assert rng[0].count("|") == 6  # energy, contrast and three metrics
    assert pencil[0].count("|") == 7  # energy, contrast and four analyses
    assert field[0].count("|") == 6  # contrast and four analyses
    assert all(line.count("|") == rng[0].count("|") for line in rng)
    assert all(line.count("|") == pencil[0].count("|") for line in pencil)
    assert all(line.count("|") == field[0].count("|") for line in field)
    assert all("unavailable" not in line for line in rng)  # every metric has a single distal crossing in all 120 runs


def test_every_gamma_cell_names_its_point_count_and_its_control() -> None:
    desc = committed_descriptive()
    for line in rt.gamma_pencil_table(desc)[2:] + rt.gamma_field_table(desc)[2:]:
        cells = [c for c in line.split("|")[1:-1]][2 if line.count("|") == 7 else 1 :]
        assert cells
        assert all("(n = " in c and "; control " in c and c.count("(n = ") == 2 for c in cells)


def test_a_table_cell_is_the_documents_number_not_a_recomputation() -> None:
    desc = committed_descriptive()
    entry = desc["gamma_field"]["contrasts"]["A"]["l2_2_c1"]  # type: ignore[index]
    control = desc["gamma_field"]["noise_control"]["A-up"]["l2_2_c1"]  # type: ignore[index]
    row = rt.gamma_field_table(desc)[2]
    assert f"{entry['pass_rate_percent']:.4f} (n = {entry['n_evaluated']}); control {control['pass_rate_percent']:.4f}" in row
    first = desc["range_metrics"]["contrasts"]["A"]["100"]["R90"]  # type: ignore[index]
    assert rt.num(first["estimate"], 4, signed=True) in rt.range_descriptive_table(desc)[2]


def test_the_noise_control_of_a_contrast_is_the_split_halves_of_its_own_upstream_arm() -> None:
    desc = committed_descriptive()
    assert rt.CONTROL_ARM == {"A": "A-up", "B1": "B-up", "B2": "B-up"}
    pencil = desc["gamma_pencil"]["100"]["3d"]  # type: ignore[index]
    assert set(pencil["noise_control"]) == {"A-up", "B-up"}
    # B1 and B2 share one control, so their control cells are identical
    rows = rt.gamma_pencil_table(desc)
    b1, b2 = rows[3].split(" | ")[2:], rows[4].split(" | ")[2:]
    assert [c.split("; control ")[1] for c in b1] == [c.split("; control ")[1] for c in b2]


def test_render_without_the_descriptive_document_renders_only_the_confirmatory_blocks() -> None:
    contrasts = rt.load(rt.ANALYSIS_JSON)
    assert "table-2-gamma-pencil" not in rt.render(contrasts)
    blocks = rt.render(contrasts, committed_descriptive())
    assert {"table-1b-range-descriptive", "table-2-gamma-pencil", "table-5-gamma-field"} <= set(blocks)


def test_the_report_holds_a_marker_pair_for_each_descriptive_block_and_no_placeholder_cells() -> None:
    report = rt.REPORT.read_text(encoding="utf-8")
    for name in ("table-1b-range-descriptive", "table-2-gamma-pencil", "table-5-gamma-field"):
        assert report.count(f"<!-- BEGIN GENERATED: {name} -->") == 1
        assert report.count(f"<!-- END GENERATED: {name} -->") == 1
    assert "| · |" not in report
    assert "· / ·" not in report


def test_every_figure_link_in_the_report_points_at_a_file_that_exists() -> None:
    report = rt.REPORT.read_text(encoding="utf-8")
    links = re.findall(r"!\[(Fig \d[^\]]*)\]\((figures/[^)]+\.png)\)", report)
    assert [name.split(".")[0] for name, _ in links] == ["Fig 1", "Fig 2", "Fig 2b", "Fig 3", "Fig 4"]
    for _, rel in links:
        assert (rt.REPORT.parent / rel).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_check_fails_when_a_descriptive_cell_is_edited(tmp_path: Path) -> None:
    text = rt.REPORT.read_text(encoding="utf-8")
    assert "| 100.000 (n = 79); control" in text
    edited = tmp_path / "report.md"
    edited.write_text(text.replace("| 100.000 (n = 79); control", "| 99.000 (n = 79); control", 1), encoding="utf-8")
    assert rt.main(["--check", "--report", str(edited)]) == 1


def test_check_refuses_to_render_from_a_descriptive_document_whose_controls_failed(tmp_path: Path) -> None:
    doc = json.loads(rt.DESCRIPTIVE_JSON.read_text(encoding="utf-8"))
    doc["controls"]["case_P_r80"]["passed"] = False
    path = tmp_path / "d.json"
    path.write_text(json.dumps(doc))
    assert rt.main(["--check", "--descriptive-json", str(path)]) == 2


def test_halo_bound_rows_are_the_annulus_endpoints_not_equivalent_somewhere() -> None:
    lines = rt.halo_bound_table(rt.load(rt.ANALYSIS_JSON), rt.load_fractions(rt.FRACTIONS_JSON))
    assert [ln.split(" | ")[0].lstrip("| ") for ln in lines[2:]] == [
        "P100/ring_40_40_80", "P100/ring_40_80_200", "P100/ring_60_40_80", "P100/ring_60_80_200",
        "P150/ring_80_80_200", "P150/ring_125_80_200",
    ]
    assert lines[3] == "| P100/ring_40_80_200 | no energy scored in any run | no interval | — |"


def test_halo_bound_product_is_share_times_the_widest_interval_bound() -> None:
    contrasts = rt.load(rt.ANALYSIS_JSON)
    eid = "P150/ring_80_80_200"
    for name in rt.CONFIRMATORY.values():
        contrasts[name]["rows"][eid]["ci95_reported_scale"] = [0.95, 1.02]
    contrasts[rt.CONFIRMATORY["B1"]]["rows"][eid]["ci95_reported_scale"] = [0.90, 1.04]
    fractions = {**rt.load_fractions(rt.FRACTIONS_JSON), eid: 0.002}
    row = next(ln for ln in rt.halo_bound_table(contrasts, fractions) if eid in ln)
    assert row == f"| {eid} | 0.2000% | 10.0% | 0.0200% |"


def test_halo_bound_refuses_a_zero_share_with_an_interval() -> None:
    fractions = {**rt.load_fractions(rt.FRACTIONS_JSON), "P150/ring_80_80_200": 0.0}
    with pytest.raises(rt.ReportError):
        rt.halo_bound_table(rt.load(rt.ANALYSIS_JSON), fractions)


def test_ring_fraction_extract_reads_only_seed_directories_of_the_upstream_arms(tmp_path: Path) -> None:
    import report_ring_fractions as rf

    for arm, value in (("A-up", 0.002), ("B-up", 0.004), ("A-port", 0.9)):
        run = tmp_path / arm / "P" / "s1"
        run.mkdir(parents=True)
        (run / "run.json").write_text(json.dumps({"energy_mev": 100}))
        (run / "endpoints.json").write_text("ENDPOINTS " + json.dumps({"ring_40_40_80": value, "R80": 77.0}) + "\n")
        (tmp_path / arm / "P" / ".provenance").mkdir()
    doc = rf.document(tmp_path)
    assert doc["mean_fraction"] == {"P100/ring_40_40_80": pytest.approx(0.003)}
    assert doc["n_runs_per_endpoint"] == 2


# ---------------------------------------------------------------------------------------------------------------------
# Table 6: TOPAS against Portable (descriptive), from the document validation/topas_halo_compare.py wrote
# ---------------------------------------------------------------------------------------------------------------------


def halo_row(**over: object) -> dict[str, object]:
    base = {"endpoint": "R80", "status": "difference", "estimate": -0.42478, "ci95": [-0.42617, -0.42340], "uncertainty": "welch"}
    return {**base, **over}


def test_a_topas_difference_cell_is_signed_to_four_decimals() -> None:
    assert rt.topas_halo_cell(halo_row()) == "−0.4248 [−0.4262, −0.4234]"


def test_a_topas_ratio_cell_is_unsigned_to_three_decimals() -> None:
    assert rt.topas_halo_cell(halo_row(endpoint="ring_60_40_80", status="ratio", estimate=0.4481, ci95=[0.4283, 0.4689])) == (
        "0.448 [0.428, 0.469]"
    )


def test_a_ring_with_a_zero_run_gives_the_run_counts_and_no_number() -> None:
    row = {
        "endpoint": "ring_40_80_200",
        "status": rt.TOPAS_NO_RATIO,
        "mcsquare": {"n": 8, "n_nonzero": 0, "mean_fraction": 0.0},
        "topas": {"n": 8, "n_nonzero": 3, "mean_fraction": 1e-9},
    }
    assert rt.topas_halo_cell(row) == "no ratio: non-zero in 0 of 8 MCsquare and 3 of 8 TOPAS runs"


def test_a_row_whose_precision_was_not_estimated_keeps_its_estimate_and_gets_no_interval() -> None:
    text = "not estimated: zero sample variance in both codes"
    assert rt.topas_halo_cell(halo_row(estimate=-0.5, ci95=None, uncertainty=text)) == f"−0.5000 ({text})"


def test_an_interval_beside_a_not_estimated_precision_is_refused() -> None:
    with pytest.raises(rt.ReportError, match="an interval beside"):
        rt.topas_halo_cell(halo_row(ci95=[-0.5, -0.5], uncertainty="not estimated: zero sample variance in both codes"))


def test_an_unknown_topas_row_status_is_refused() -> None:
    with pytest.raises(rt.ReportError, match="unknown status"):
        rt.topas_halo_cell(halo_row(status="equivalent"))


def test_table_6_has_one_row_per_endpoint_and_every_cell_is_the_documents_number() -> None:
    rows = rt.load_topas_halo(rt.TOPAS_HALO_JSON)
    assert len(rows) == 117  # 3 energies x (R80 + 2 slabs x (sigma + 5 rings)) x 3 arms
    lines = rt.topas_halo_table(rows)
    assert len(lines) - 2 == 39
    r = rows[(150, "ring_125_80_200", "B-picc")]
    wanted = f"{r['estimate']:.3f} [{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}]"
    line = next(x for x in lines if "80–200 mm annulus at 125 mm" in x)
    assert line.split(" | ")[4] == wanted
    assert sum("no ratio: non-zero in 0 of 8 MCsquare and 0 of 8 TOPAS runs" in x for x in lines) == 2


def test_td_bands_appear_only_where_the_document_carries_them() -> None:
    lines = rt.topas_halo_table(rt.load_topas_halo(rt.TOPAS_HALO_JSON))
    banded = [x for x in lines[2:] if not x.endswith("|  |")]
    assert len(banded) == 9
    assert all(x in lines[-13:] for x in banded)  # the 200 MeV rows
    assert lines[-13].endswith("| [−0.30, +0.30] |")
    assert lines[-1].endswith("| [0.75, 1.25] |")


@pytest.mark.parametrize("change", ["drop", "duplicate", "extra", "arms"])
def test_a_topas_document_that_is_not_exactly_the_tables_rows_is_refused(tmp_path: Path, change: str) -> None:
    doc = json.loads(rt.TOPAS_HALO_JSON.read_text(encoding="utf-8"))
    if change == "drop":
        doc["rows"].pop()
    elif change == "duplicate":
        doc["rows"].append(doc["rows"][0])
    elif change == "extra":
        doc["rows"].append({**doc["rows"][0], "endpoint": "sigma_3"})
    else:
        doc["mcsquare"]["arms"] = ["A-port", "B-pgcc"]
    path = tmp_path / "h.json"
    path.write_text(json.dumps(doc))
    with pytest.raises(rt.ReportError):
        rt.load_topas_halo(path)
    assert rt.main(["--check", "--topas-halo-json", str(path)]) == 2


def test_check_fails_when_a_table_6_cell_is_edited(tmp_path: Path) -> None:
    text = rt.REPORT.read_text(encoding="utf-8")
    assert text.count("| −0.4253 [−0.4260, −0.4247] |") == 1
    edited = tmp_path / "report.md"
    edited.write_text(text.replace("| −0.4253 [−0.4260, −0.4247] |", "| −0.0253 [−0.4260, −0.4247] |"), encoding="utf-8")
    assert rt.main(["--check", "--report", str(edited)]) == 1


def test_render_without_the_topas_document_has_no_table_6() -> None:
    assert "table-6-topas-halo" not in rt.render(rt.load(rt.ANALYSIS_JSON))


def test_check_names_the_block_that_is_out_of_date(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    text = rt.REPORT.read_text(encoding="utf-8")
    edited = tmp_path / "report.md"
    edited.write_text(text.replace("| −0.4253 [−0.4260, −0.4247] |", "| −0.0253 [−0.4260, −0.4247] |"), encoding="utf-8")
    assert rt.main(["--check", "--report", str(edited)]) == 1
    assert "(blocks: table-6-topas-halo)" in capsys.readouterr().err
    edited.write_text(text.replace("| A | not established | 46 / 46 |", "| A | holds | 52 / 52 |"), encoding="utf-8")
    assert rt.main(["--check", "--report", str(edited)]) == 1
    assert "(blocks: confirmatory-summary)" in capsys.readouterr().err


# ---------------------------------------------------------------------------------------------- part E (§5.5)
def test_part_e_tables_show_every_endpoint_of_every_confirmatory_contrast_with_the_documents_counts() -> None:
    part_e = rt.load_part_e(rt.PART_E_JSON)
    blocks = rt.render(rt.load(rt.ANALYSIS_JSON), part_e=part_e)
    rows = blocks["table-8-part-e-endpoints"][2:]
    assert len(rows) == rt.PART_E_ENDPOINTS == 34
    assert all(row.count(" E") >= 3 or "not established" in row or " I" in row or " NE" in row for row in rows)
    summary = blocks["table-7-part-e-summary"][2:]
    for line, name in zip(summary, rt.CONFIRMATORY.values(), strict=True):
        claims = part_e[name]["claims"]
        assert f"{claims['n_equivalent_unadjusted']} of 34" in line and f"{claims['n_equivalent_holm']} of 34" in line


@pytest.mark.parametrize("damage", ["withheld", "partial", "missing", "short", "renamed", "duplicate_endpoint",
                                    "duplicate_contrast"])
def test_part_e_load_refuses_an_incomplete_confirmatory_contrast(tmp_path: Path, damage: str) -> None:
    doc = json.loads(rt.PART_E_JSON.read_text(encoding="utf-8"))
    c = next(x for x in doc["contrasts"] if x["confirmatory"])
    if damage == "withheld":
        c["claims_withheld"] = True
    elif damage == "partial":
        c["partial_reasons"] = ["a run is missing"]
    elif damage == "missing":
        doc["contrasts"].remove(c)
    elif damage == "short":
        c["rows"] = c["rows"][:-1]
    elif damage == "renamed":  # same size, wrong identity, in every confirmatory contrast (alden-ec2221c7, #72)
        for x in doc["contrasts"]:
            if x["confirmatory"]:
                x["rows"][0]["endpoint"] = "NOT_A_STUDY_ENDPOINT"
    elif damage == "duplicate_endpoint":
        c["rows"][-1] = dict(c["rows"][0])
    else:
        doc["contrasts"].append(json.loads(json.dumps(c)))
    path = tmp_path / "e.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(rt.ReportError):
        rt.load_part_e(path)


# ---------------------------------------------------------------------------------------------- run times (Table 9)
def test_run_time_table_has_one_row_per_arm_on_one_host_each() -> None:
    lines = rt.run_time_table(rt.load_run_times(rt.RUN_TIMES_JSON))
    assert [ln.split(" | ")[0].lstrip("| ") for ln in lines[2:]] == ["A-up", "A-port", "B-up", "B-pgcc", "B-picc"]


def test_run_times_refuse_another_schema_and_a_missing_cell(tmp_path: Path) -> None:
    doc = json.loads(rt.RUN_TIMES_JSON.read_text(encoding="utf-8"))
    for damage in ("schema", "cell"):
        bad = json.loads(json.dumps(doc))
        if damage == "schema":
            bad["schema"] = "run_times/0"
        else:
            bad["groups"] = [g for g in bad["groups"] if not (g["arm"] == "B-picc" and g["case"] == "E")]
        path = tmp_path / f"{damage}.json"
        path.write_text(json.dumps(bad), encoding="utf-8")
        with pytest.raises(rt.ReportError):
            rt.load_run_times(path)
