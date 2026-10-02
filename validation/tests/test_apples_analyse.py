"""Tests for apples_analyse on SYNTHETIC trees only (no run output from any machine is read).

Run: uv run --no-project --with numpy --with scipy --with pytest pytest validation/tests/test_apples_analyse.py
"""

from __future__ import annotations

import json
import math
import random
import sys
from collections.abc import Callable
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import apples_analyse as aa

# isort: split
from apples_fixtures import (  # noqa: F401 - the autouse fixture must be in this module's namespace
    _REAL_EXPECTED,
    ACQ,
    ALL_ARMS,
    N_RUNS,
    SOURCES,
    V2_BINDING,
    V2_COMMIT,
    Mutate,
    RunMutate,
    _frozen_from_the_working_tree,
    attest,
    binary_of,
    build,
    cells,
    f_record,
    p_record,
    run_json_of,
)


def analyse(root: Path, parts: tuple[str, ...] = ("A",), *, attested: bool = True, commit: str = ACQ) -> dict:
    """The analysis of a flat tree. `attested` (default) first writes a collection manifest over the tree as it stands,
    standing in for apples_collect.py, unless one exists; without it every contrast is descriptive only."""
    if attested and not (root / aa.MANIFEST_NAME).exists():
        attest(root, parts, commit=commit)
    return aa.run_analysis(root, parts, commit=commit)[1]


def contrast(doc: dict, name: str) -> dict:
    return next(c for c in doc["contrasts"] if c["name"] == name)


def ep(doc: dict, name: str, eid: str) -> dict:
    return next(e for e in contrast(doc, name)["endpoints"] if e["endpoint"] == eid)


A_CONTRAST = "A-port vs A-up"


@pytest.fixture
def a_identical(tmp_path: Path) -> Path:
    return build(tmp_path, ("A-up", "A-port"), identical=True)


def assert_no_confirmatory_claims(c: dict) -> None:
    """A PARTIAL confirmatory contrast: no joint claim, no Holm decision, no TOST classification, labelled partial."""
    assert c["partial"] is True and c["partial_reasons"]
    assert c["claims_withheld"] is True
    assert c["claims"] == {}
    assert "joint_claim_equivalent_on_all_endpoints" not in c["claims"]
    for e in c["endpoints"]:
        assert e["outcome"] in ("descriptive", "not_established"), e
        assert e["holm_decision"] == "" and e["holm_adjusted_p"] is None and e["p_tost"] is None


# --------------------------------------------------------------------------- the family and its margins


def test_family_is_39_p_plus_13_f_endpoints() -> None:
    specs = aa.all_specs()
    assert len(specs) == 52
    assert sum(s.case == "P" for s in specs) == 39
    assert sum(s.case == "F" for s in specs) == 13
    assert sum(s.eid.startswith("P100/") for s in specs) == 13


def test_p_margins_are_as_pre_specified() -> None:
    specs = {s.eid: s for s in aa.p_specs()}
    for e, (d1, d2) in aa.SLAB_DEPTHS.items():
        assert (specs[f"P{e}/R80"].lo, specs[f"P{e}/R80"].hi) == (-0.05, 0.05)
        for d in (d1, d2):
            assert (specs[f"P{e}/sigma_{d}"].lo, specs[f"P{e}/sigma_{d}"].hi) == (-0.02, 0.02)
            for ring in ("5_10", "10_20", "20_40", "40_80"):
                s = specs[f"P{e}/ring_{d}_{ring}"]
                assert s.log
                assert (math.exp(s.lo), math.exp(s.hi)) == pytest.approx((0.98, 1.02))
            wide = specs[f"P{e}/ring_{d}_80_200"]
            assert (math.exp(wide.lo), math.exp(wide.hi)) == pytest.approx((0.95, 1.05))


def test_slab_depths_per_energy() -> None:
    assert aa.SLAB_DEPTHS == {100: (40, 60), 150: (80, 125), 200: (100, 200)}


def test_f_margins_are_the_31_margins() -> None:
    specs = {s.eid: s for s in aa.f_specs()}
    assert len(specs) == 13
    assert (specs["F/lateral_127_5"].lo, specs["F/lateral_127_5"].hi) == (-0.5, 0.5)
    assert (specs["F/lateral_201_5"].lo, specs["F/lateral_201_5"].hi) == (-0.5, 0.5)
    for off in (10, 20, 30):
        assert specs[f"F/lateral_127_{off}"].hi == 0.2
        assert specs[f"F/lateral_201_{off}"].hi == 0.2
    for k in ("cax_127", "cax_201", "cax_i23"):
        assert specs[f"F/{k}"].log
        assert (math.exp(specs[f"F/{k}"].lo), math.exp(specs[f"F/{k}"].hi)) == pytest.approx((0.995, 1.005))
    assert (specs["F/r80_mm"].lo, specs["F/r80_mm"].hi) == (-0.5, 0.5)
    assert (specs["F/r20_mm"].lo, specs["F/r20_mm"].hi) == (-0.5, 0.5)


# ------------------------------------------------------------------------------ estimator conventions


def test_welch_hand_checked() -> None:
    d, se, df = aa.welch([1.0, 2.0, 3.0], [2.0, 4.0, 6.0])
    assert d == pytest.approx(-2.0)
    assert se == pytest.approx(math.sqrt(5 / 3))
    assert df == pytest.approx(50 / 17)


def test_welch_needs_two_runs_per_arm() -> None:
    with pytest.raises(ValueError, match="at least 2"):
        aa.welch([1.0], [1.0, 2.0])


def test_zero_variance_both_arms_gives_degenerate_interval() -> None:
    d, se, df = aa.welch([1.0, 1.0], [1.0, 1.0])
    assert (d, se, df) == (0.0, 0.0, math.inf)
    assert aa.interval(d, se, df, 0.90) == (0.0, 0.0)


def test_equivalent_iff_tost_p_below_alpha() -> None:
    """The 90% interval inside the margin and p_TOST < 0.05 are one test, across a grid that avoids exact boundaries."""
    disagreements = 0
    for d in (-0.019, -0.012, -0.004, 0.0, 0.006, 0.011, 0.017, 0.0195, 0.03):
        for se in (0.0007, 0.003, 0.006, 0.009, 0.02):
            for df in (2.5, 5.0, 11.0, 80.0):
                ci90 = aa.interval(d, se, df, 0.90)
                inside = aa.classify(ci90, -0.02, 0.02) == "equivalent"
                if inside != (aa.tost_p(d, se, df, -0.02, 0.02) < aa.ALPHA):
                    disagreements += 1
    assert disagreements == 0


def test_tost_uses_90_not_95_interval() -> None:
    """t = 3 with df = 3: 90% interval (+-2.353 se) is inside the margin, the 95% one (+-3.182 se) is not."""
    se = 0.02 / 3.0
    ci90, ci95 = aa.interval(0.0, se, 3.0, 0.90), aa.interval(0.0, se, 3.0, 0.95)
    assert aa.classify(ci90, -0.02, 0.02) == "equivalent"
    assert aa.classify(ci95, -0.02, 0.02) == "inconclusive"
    assert aa.tost_p(0.0, se, 3.0, -0.02, 0.02) == pytest.approx(0.0288, abs=0.0005)


def test_p_tost_is_the_larger_one_sided_p() -> None:
    se, df = 0.005, 8.0
    p = aa.tost_p(0.012, se, df, -0.02, 0.02)  # closer to the upper bound
    from scipy.stats import t

    assert p == pytest.approx(float(t.sf((0.02 - 0.012) / se, df)))
    assert p > float(t.sf((0.012 + 0.02) / se, df))


# ------------------------------------------------------------------------------------ Holm


def test_holm_hand_checked() -> None:
    # sorted 0.005, 0.01, 0.03, 0.04 -> x4, x3, x2, x1 = 0.02, 0.03, 0.06, 0.04 -> running max 0.02, 0.03, 0.06, 0.06
    adj = aa.holm_adjust([0.03, 0.04, 0.01, 0.005])
    assert adj == pytest.approx([0.06, 0.06, 0.03, 0.02])


def test_holm_caps_at_one_and_is_monotone_in_p() -> None:
    # sorted 0.2, 0.4, 0.9 -> x3, x2, x1 = 0.6, 0.8, 0.9
    assert aa.holm_adjust([0.4, 0.9, 0.2]) == pytest.approx([0.8, 0.9, 0.6])
    capped = aa.holm_adjust([0.6, 0.7])  # 2 * 0.6 = 1.2 -> capped at 1
    assert capped == [1.0, 1.0]


def test_holm_keeps_not_established_in_the_family_and_never_rejects_it() -> None:
    # m = 3 including the None: 0.001*3, 0.02*2. Dropping the None (m = 2) would give 0.002 and 0.02.
    adj = aa.holm_adjust([0.001, None, 0.02])
    assert adj[0] == pytest.approx(0.003)
    assert adj[1] is None
    assert adj[2] == pytest.approx(0.04)


# --------------------------------------------------------------------------- end-to-end behaviour


def test_identical_arms_all_equivalent_and_joint_claim_true(a_identical: Path) -> None:
    doc = analyse(a_identical)
    c = contrast(doc, A_CONTRAST)
    assert len(c["endpoints"]) == 52
    assert {e["outcome"] for e in c["endpoints"]} == {"equivalent"}
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is True
    assert c["claims"]["n_equivalent_holm"] == 52
    assert {e["holm_decision"] for e in c["endpoints"]} == {"equivalent"}


def test_independent_draws_of_the_same_distribution_are_equivalent(tmp_path: Path) -> None:
    doc = analyse(build(tmp_path, ("A-up", "A-port")))
    c = contrast(doc, A_CONTRAST)
    assert all(e["df"] is not None and math.isfinite(e["df"]) for e in c["endpoints"])  # the Welch path ran
    assert {e["outcome"] for e in c["endpoints"]} == {"equivalent"}
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is True


def test_one_shifted_endpoint_is_not_equivalent_and_the_joint_claim_fails(tmp_path: Path) -> None:
    def shift(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 150:
            rec["R80"] += 0.2  # 4x the margin, noise sd is 0.005

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=shift))
    c = contrast(doc, A_CONTRAST)
    assert ep(doc, A_CONTRAST, "P150/R80")["outcome"] == "not_equivalent"
    assert [e["endpoint"] for e in c["endpoints"] if e["outcome"] != "equivalent"] == ["P150/R80"]
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is False
    assert ep(doc, A_CONTRAST, "P150/R80")["holm_decision"] == "not shown"


def test_wide_variance_endpoint_is_inconclusive(tmp_path: Path) -> None:
    def noisy(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if case == "P" and energy == 100:
            rec["sigma_40"] = 4.0 + random.Random(i * 13 + len(arm)).gauss(0.0, 0.05)  # sd 0.05 against a 0.02 margin

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=noisy))
    e = ep(doc, A_CONTRAST, "P100/sigma_40")
    assert e["outcome"] == "inconclusive"
    assert e["ci90"][0] < -0.02 or e["ci90"][1] > 0.02
    assert contrast(doc, A_CONTRAST)["claims"]["joint_claim_equivalent_on_all_endpoints"] is False


def test_ring_80_200_margin_is_5_percent_and_the_others_2_percent(tmp_path: Path) -> None:
    def ratio_3pc(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 200:
            rec["ring_100_80_200"] *= 1.03  # inside 5%, outside 2%
            rec["ring_100_40_80"] *= 1.03
            rec["ring_200_80_200"] *= 1.06  # outside 5%

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=ratio_3pc))
    assert ep(doc, A_CONTRAST, "P200/ring_100_80_200")["outcome"] == "equivalent"
    assert ep(doc, A_CONTRAST, "P200/ring_100_40_80")["outcome"] == "not_equivalent"
    assert ep(doc, A_CONTRAST, "P200/ring_200_80_200")["outcome"] == "not_equivalent"
    assert ep(doc, A_CONTRAST, "P200/ring_100_80_200")["estimate_reported_scale"] == pytest.approx(1.03, abs=0.002)


def test_joint_claim_is_not_touched_by_holm(tmp_path: Path) -> None:
    """Every endpoint at p_TOST of about 0.03: each equivalent at unadjusted 0.05, none rejected by Holm over 52."""

    def spread(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm != "A-port":
            return
        for s in aa.all_specs():
            if s.case != case or s.energy != energy:
                continue
            half = min(-s.lo, s.hi) / 2.5 * math.sqrt(N_RUNS - 1)  # alternating +-half: se = half / sqrt(n-1) = margin / 2.5
            sgn = -1.0 if i % 2 == 0 else 1.0
            rec[s.key] = rec[s.key] * math.exp(sgn * half) if s.log else rec[s.key] + sgn * half

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, mutate=spread))
    c = contrast(doc, A_CONTRAST)
    assert not c["partial"]
    assert all(0.01 < e["p_tost"] < 0.05 for e in c["endpoints"])
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is True
    assert c["claims"]["n_equivalent_unadjusted"] == 52
    assert c["claims"]["n_equivalent_holm"] == 0
    assert {e["holm_decision"] for e in c["endpoints"]} == {"not shown"}


# ---------------------------------------------------------------------------------- invalid values


def test_run_with_transport_status_not_ok_makes_its_endpoints_not_established_and_the_contrast_partial(
    tmp_path: Path,
) -> None:
    def fail(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 100 and i == 0:
            run["transport_status"] = "failed"

    doc = analyse(build(tmp_path, ("A-up", "A-port"), run_mutate=fail))
    c = contrast(doc, A_CONTRAST)
    bad = {e["endpoint"] for e in c["endpoints"] if e["outcome"] == "not_established"}
    assert bad == {s.eid for s in aa.p_specs() if s.energy == 100}  # that energy's 13, nothing else
    assert_no_confirmatory_claims(c)
    assert [r for r in c["partial_reasons"] if "run not usable: transport_status 'failed'" in r]
    assert doc["runs"]["A-port/P"]["n"] == 3 * N_RUNS  # the run stays in the population
    assert len(doc["runs"]["A-port/P"]["unusable"]) == 1


def test_f_run_not_ok_makes_all_13_f_endpoints_not_established(tmp_path: Path) -> None:
    def fail(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-up" and case == "F" and i == 2:
            run["transport_status"] = "timeout"

    doc = analyse(build(tmp_path, ("A-up", "A-port"), run_mutate=fail))
    bad = {e["endpoint"] for e in contrast(doc, A_CONTRAST)["endpoints"] if e["outcome"] == "not_established"}
    assert bad == {s.eid for s in aa.f_specs()}
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))


def test_missing_transport_status_is_not_ok(tmp_path: Path) -> None:
    def drop(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-up" and case == "P" and energy == 200 and i == 0:
            del run["transport_status"]

    doc = analyse(build(tmp_path, ("A-up", "A-port"), run_mutate=drop))
    assert ep(doc, A_CONTRAST, "P200/R80")["outcome"] == "not_established"
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))


def test_none_sigma_makes_only_that_endpoint_not_established(tmp_path: Path) -> None:
    def fit_failed(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 150 and i == 1:
            rec["sigma_125"] = None

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=fit_failed))
    c = contrast(doc, A_CONTRAST)
    assert [e["endpoint"] for e in c["endpoints"] if e["outcome"] == "not_established"] == ["P150/sigma_125"]
    assert ep(doc, A_CONTRAST, "P150/sigma_80")["outcome"] == "equivalent"
    assert ep(doc, A_CONTRAST, "P150/ring_125_40_80")["outcome"] == "equivalent"
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is False
    assert c["partial"] is False  # an invalid value in a usable run fails its endpoint, not the population


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_non_finite_value_is_not_established(tmp_path: Path, bad: float) -> None:
    def poison(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-up" and case == "F" and i == 0:
            rec["lateral_127_10"] = bad

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=poison))
    assert ep(doc, A_CONTRAST, "F/lateral_127_10")["outcome"] == "not_established"


@pytest.mark.parametrize("flag", [True, None, "no"])
def test_r80_multiple_crossings_or_unknown_flag_is_not_established(tmp_path: Path, flag: object) -> None:
    def crossing(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 100 and i == 0:
            rec["R80_multiple_crossings"] = flag

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=crossing))
    assert ep(doc, A_CONTRAST, "P100/R80")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "P100/sigma_40")["outcome"] == "equivalent"  # only R80 is affected


def test_missing_r80_is_not_established(tmp_path: Path) -> None:
    def gone(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 100 and i == 0:
            rec["R80"] = None

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=gone))
    assert ep(doc, A_CONTRAST, "P100/R80")["outcome"] == "not_established"


@pytest.mark.parametrize("v", [0.0, -0.01, None])
def test_zero_ring_fraction_is_not_established_and_does_not_crash(tmp_path: Path, v: float | None) -> None:
    def zero(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-up" and case == "P" and energy == 200 and i == 3:
            rec["ring_200_5_10"] = v

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=zero))
    assert ep(doc, A_CONTRAST, "P200/ring_200_5_10")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "P200/ring_200_10_20")["outcome"] == "equivalent"


def test_f_cax_invalid_flag_and_nonpositive_dose_not_established(tmp_path: Path) -> None:
    def bad(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-port" and case == "F" and i == 0:
            rec["cax_127_invalid"] = True
            rec["cax_201"] = 0.0

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=bad))
    assert ep(doc, A_CONTRAST, "F/cax_127")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "F/cax_201")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "F/cax_i23")["outcome"] == "equivalent"


def test_f_r80_ambiguous_or_absent_not_established(tmp_path: Path) -> None:
    def bad(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-up" and case == "F" and i == 0:
            rec["r80_crossings"] = 2
            rec["r20_mm"] = None

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=bad))
    assert ep(doc, A_CONTRAST, "F/r80_mm")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "F/r20_mm")["outcome"] == "not_established"


def test_f_r20_ambiguous_not_established_as_in_31(tmp_path: Path) -> None:
    def bad(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-up" and case == "F" and i == 0:
            rec["r20_crossings"] = 2

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=bad))
    assert ep(doc, A_CONTRAST, "F/r20_mm")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "F/r80_mm")["outcome"] == "equivalent"


def test_fewer_than_two_runs_in_an_arm_is_not_established_and_the_contrast_is_partial(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"), n=1)
    doc = analyse(root, attested=False)
    c = contrast(doc, A_CONTRAST)
    assert ep(doc, A_CONTRAST, "P150/R80")["outcome"] == "not_established"
    assert_no_confirmatory_claims(c)


def test_unreadable_run_json_keeps_the_run_unusable_and_the_contrast_partial(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    run = min((root / "A-port" / "P").iterdir())  # a 100 MeV run
    (run / "run.json").write_text("{not json")
    doc = analyse(root)
    assert doc["runs"]["A-port/P"]["n"] == 3 * N_RUNS
    assert len(doc["runs"]["A-port/P"]["unusable"]) == 1
    assert ep(doc, A_CONTRAST, "P100/R80")["outcome"] == "not_established"  # the run is placed by the frozen list
    assert ep(doc, A_CONTRAST, "P150/R80")["estimate"] is not None  # unaffected energies keep their estimates
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))


def test_missing_endpoint_file_is_an_unusable_run_not_a_crash(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    (min((root / "A-up" / "F").iterdir()) / "record.json").unlink()
    doc = analyse(root)
    assert {e["outcome"] for e in contrast(doc, A_CONTRAST)["endpoints"] if e["endpoint"].startswith("F/")} == {
        "not_established"
    }
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))


@pytest.mark.parametrize("n_lines", [0, 2])
def test_endpoints_file_without_exactly_one_endpoints_line_makes_the_contrast_partial(tmp_path: Path, n_lines: int) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    path = min((root / "A-port" / "P").iterdir()) / "endpoints.json"  # a 100 MeV run
    line = next(ln for ln in path.read_text().splitlines() if ln.startswith("ENDPOINTS "))
    path.write_text("".join(f"{line}\n" for _ in range(n_lines)) or "no endpoints here\n")
    doc = analyse(root)
    c = contrast(doc, A_CONTRAST)
    assert ep(doc, A_CONTRAST, "P100/R80")["outcome"] == "not_established"
    assert_no_confirmatory_claims(c)
    assert [r for r in c["partial_reasons"] if f"{n_lines} ENDPOINTS lines" in r]


# ------------------------------------------------------------------------------- part B and descriptive


def test_part_b_has_two_confirmatory_families_and_one_descriptive_contrast(tmp_path: Path) -> None:
    def shift(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "B-pgcc" and case == "P" and energy == 100:
            rec["sigma_40"] += 0.5  # far outside the margin

    doc = analyse(build(tmp_path, ("B-up", "B-pgcc", "B-picc"), mutate=shift), ("B",))
    names = [c["name"] for c in doc["contrasts"]]
    assert names == ["B-pgcc vs B-up", "B-picc vs B-up", "B-pgcc vs B-picc (descriptive)"]
    assert ep(doc, "B-pgcc vs B-up", "P100/sigma_40")["outcome"] == "not_equivalent"
    assert contrast(doc, "B-pgcc vs B-up")["claims"]["joint_claim_equivalent_on_all_endpoints"] is False
    assert contrast(doc, "B-picc vs B-up")["claims"]["joint_claim_equivalent_on_all_endpoints"] is True


def test_descriptive_contrast_makes_no_equivalence_claim(tmp_path: Path) -> None:
    def shift(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "B-pgcc" and case == "P" and energy == 100:
            rec["sigma_40"] += 0.5

    doc = analyse(build(tmp_path, ("B-up", "B-pgcc", "B-picc"), mutate=shift), ("B",))
    d = contrast(doc, "B-pgcc vs B-picc (descriptive)")
    assert d["kind"] == "descriptive"
    assert d["claims"] == {}
    assert {e["outcome"] for e in d["endpoints"]} == {"descriptive"}
    assert all(e["margin"] is None and e["p_tost"] is None and e["holm_adjusted_p"] is None for e in d["endpoints"])
    e = ep(doc, "B-pgcc vs B-picc (descriptive)", "P100/sigma_40")
    assert e["estimate"] == pytest.approx(0.5, abs=0.01)  # estimates and intervals are still reported
    assert e["ci90"][0] < e["estimate"] < e["ci90"][1]
    assert e["ci95"][0] < e["ci90"][0]


def test_contrast_direction_is_arm_minus_reference_and_ratio_for_rings(tmp_path: Path) -> None:
    def shift(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 200:
            rec["sigma_100"] += 0.1
            rec["ring_100_20_40"] *= 1.1

    doc = analyse(build(tmp_path, ("A-up", "A-port"), mutate=shift))
    assert ep(doc, A_CONTRAST, "P200/sigma_100")["estimate"] == pytest.approx(0.1, abs=0.01)
    assert ep(doc, A_CONTRAST, "P200/ring_100_20_40")["estimate_reported_scale"] == pytest.approx(1.1, abs=0.01)


# --------------------------------------------------------------------------------------- CLI and outputs


def test_cli_outputs_markdown_tables_and_strict_json(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    def constant(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if case == "P" and energy == 100:
            rec["R80"] = 77.5  # zero variance in both arms: df is infinite, which JSON must carry as null

    root = build(tmp_path / "in", ("A-up", "A-port"), identical=True, mutate=constant)
    attest(root)
    out = tmp_path / "result.json"
    rc = aa.main([str(root), "--json", str(out), "--parts", "A"])
    assert rc == 0
    md = capsys.readouterr().out
    assert md.count("### ") == 1
    assert "A-port vs A-up" in md
    assert "Joint claim" in md
    assert md.count("| P100/R80 |") == 1 and md.count("| F/cax_i23 |") == 1

    def refuse(token: str) -> None:
        raise AssertionError(token)

    doc = json.loads(out.read_text(), parse_constant=refuse)  # NaN/Infinity would be refused
    r80 = next(e for e in doc["contrasts"][0]["endpoints"] if e["endpoint"] == "P100/R80")
    assert r80["df"] is None and r80["se"] == 0.0 and r80["outcome"] == "equivalent"
    assert math.isfinite(next(e for e in doc["contrasts"][0]["endpoints"] if e["endpoint"] == "P150/R80")["df"])


def test_cli_markdown_has_one_table_per_contrast(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    root = build(tmp_path, ALL_ARMS, n=3)
    assert aa.main([str(root)]) == 0
    assert capsys.readouterr().out.count("### ") == 4


def test_exit_0_whatever_the_outcomes(tmp_path: Path) -> None:
    def shift(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-port" and case == "F":
            rec["r80_mm"] += 10.0

    root = build(tmp_path, ("A-up", "A-port"), mutate=shift)
    assert aa.main([str(root), "--parts", "A"]) == 0


def test_exit_2_on_missing_arm(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    root = build(tmp_path, ("A-up",))
    assert aa.main([str(root), "--parts", "A"]) == 2
    assert "A-port" in capsys.readouterr().err


def test_exit_2_on_missing_case_dir(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    for f in list((root / "A-port" / "F").rglob("*")):
        if f.is_file():
            f.unlink()
    for d in sorted((root / "A-port" / "F").iterdir()):
        d.rmdir()
    (root / "A-port" / "F").rmdir()
    assert aa.main([str(root), "--parts", "A"]) == 2
    assert "F" in capsys.readouterr().err


def test_exit_2_on_unexpected_layout(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    (root / "A-up" / "P" / "smoke").mkdir()
    assert aa.main([str(root), "--parts", "A"]) == 2
    assert "smoke" in capsys.readouterr().err


def test_exit_2_on_empty_case_dir(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    for d in (root / "A-port" / "F").iterdir():
        for f in d.iterdir():
            f.unlink()
        d.rmdir()
    assert aa.main([str(root), "--parts", "A"]) == 2


def test_exit_2_when_run_json_disagrees_with_its_directory(tmp_path: Path) -> None:
    def mislabel(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 100 and i == 0:
            run["arm"] = "A-up"

    root = build(tmp_path, ("A-up", "A-port"), run_mutate=mislabel)
    assert aa.main([str(root), "--parts", "A"]) == 2


def test_run_json_energy_that_contradicts_its_seeds_frozen_energy_is_a_provenance_failure(tmp_path: Path) -> None:
    def odd(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-up" and case == "P" and energy == 100 and i == 0:
            run["energy_mev"] = 70

    doc = analyse(build(tmp_path, ("A-up", "A-port"), run_mutate=odd))
    assert ep(doc, A_CONTRAST, "P100/R80")["outcome"] == "not_established"
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))


@pytest.mark.parametrize("parts", ["C", "A,A", ""])
def test_exit_2_on_bad_parts(tmp_path: Path, parts: str) -> None:
    assert aa.main([str(tmp_path), "--parts", parts]) == 2


def test_exit_2_when_root_is_not_a_directory(tmp_path: Path) -> None:
    assert aa.main([str(tmp_path / "nope")]) == 2


# ------------------------------------------------------------------ identity source and design completeness


def test_report_says_case_p_endpoint_values_are_taken_as_written(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    doc = analyse(root)
    assert doc["p_endpoints_not_verifiable"] == list(aa.P_ENDPOINTS_NOT_VERIFIABLE)
    assert aa.main([str(root), "--parts", "A"]) == 0
    assert "Case-P endpoints. Not verifiable: endpoint values: endpoints.json is written on the run host" in capsys.readouterr().out


def test_report_header_says_where_identity_comes_from(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    attest(root)
    assert aa.main([str(root), "--parts", "A"]) == 0
    out = capsys.readouterr().out
    assert "frozen workflow run lists at the acquisition commit" in out and "record binding" in out
    assert "Case-F record schema: platform_study_record@2f9dab40 (binding pe1-at-acquisition-2f9dab40)" in out
    assert "this schema carries no metrics hash" in out
    assert f"Acquisition commit (frozen run lists): {ACQ}" in out


def _set_study(root: Path, study: object) -> None:
    for rec in root.glob("*/F/s*/record.json"):
        body = json.loads(rec.read_text())
        body.pop("study", None)
        if study is not None:
            body["study"] = study
        rec.write_text(json.dumps(body))


def test_legacy_pe1_record_is_accepted_only_through_the_named_binding(tmp_path: Path) -> None:
    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True))
    c = contrast(doc, A_CONTRAST)
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is True and not c["partial"]
    assert doc["runs"]["A-up/F"]["record_bindings_used"] == ["pe1-at-acquisition-2f9dab40"]
    assert doc["runs"]["A-up/P"]["record_bindings_used"] == []
    assert doc["f_record_binding"]["schema"] == "platform_study_record@2f9dab40"
    assert any("no metrics hash" in x for x in doc["f_record_binding"]["not_verifiable"])


def test_the_fixture_records_are_the_pinned_writers_and_carry_none_of_the_later_fields(tmp_path: Path) -> None:
    """The F records the tests analyse come from running 2f9dab40's writer: exactly its ten fields, no endpoint status."""
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    for rec in root.glob("*/F/s*/record.json"):
        body = json.loads(rec.read_text())
        assert set(body) == aa.LEGACY_RECORD_FIELDS and set(body["sha256"]) == aa.F_RECORD_SHA256_KEYS
        assert not set(body) & aa.V2_ENDPOINT_FIELDS and body["study"] == "pe1"


@pytest.mark.parametrize("study", ["apples-a", "PE1", "garbage", ""])
def test_any_other_study_id_is_a_provenance_failure(tmp_path: Path, study: object) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    _set_study(root, study)
    c = contrast(analyse(root), A_CONTRAST)
    assert_no_confirmatory_claims(c)
    assert any("is not the study 'pe1' bound by pe1-at-acquisition-2f9dab40" in r for r in c["partial_reasons"])


def test_a_record_without_a_study_field_refuses(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    _set_study(root, None)
    with pytest.raises(aa.InputError, match="lacks field"):
        analyse(root)


def _legacy_doc(**over: object) -> dict:
    doc: dict = {"study": "pe1", "platform": "A-up", "seed": 1, "host": "h", "machine": "m", "commit": ACQ,
                 "compiler": "cc", "sha256": dict.fromkeys(aa.F_RECORD_SHA256_KEYS, "a" * 64),
                 "materials": {"files": 1, "combined_sha256": "c" * 64}, "metrics": {}}
    doc["sha256"]["binary"] = "b" * 64
    doc.update(over)
    return doc


def test_a_legacy_record_under_a_commit_without_its_binding_refuses() -> None:
    with pytest.raises(aa.InputError, match="no case-F record binding"):
        aa._reconcile_f(_legacy_doc(commit="f" * 40), {"binary_sha256": "b" * 64}, "A-up", 1, "f" * 40)
    problems, binding = aa._reconcile_f(_legacy_doc(), {"binary_sha256": "b" * 64}, "A-up", 1, ACQ)
    assert binding == "pe1-at-acquisition-2f9dab40" and problems == []


def test_a_legacy_record_under_a_commit_bound_to_the_v2_schema_refuses(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(aa, "F_RECORD_BINDINGS", (*aa.F_RECORD_BINDINGS, V2_BINDING))
    with pytest.raises(aa.InputError, match="bound to platform_study_record@0f5ef7c"):
        aa._reconcile_f(_legacy_doc(commit=V2_COMMIT), {"binary_sha256": "b" * 64}, "A-up", 1, V2_COMMIT)


def test_mode_and_commit_are_reported_from_run_json(tmp_path: Path) -> None:
    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True))
    assert doc["runs"]["A-port/P"]["modes"] == ["full"]
    assert doc["runs"]["A-port/P"]["commits"] == [ACQ]
    assert doc["runs"]["A-up/F"]["per_energy"] == {"200": N_RUNS}


def test_complete_design_is_not_partial_and_expected_count_is_8() -> None:
    assert aa.RUNS_PER_CELL == 8 == N_RUNS




# ------------------------------------------------------------------------------- the frozen population (item 1)


def _drop(root: Path, arm: str, case: str, seed: int | str) -> None:
    victim = root / arm / case / (seed if isinstance(seed, str) else f"s{seed}")
    for f in victim.iterdir():
        f.unlink()
    victim.rmdir()


def _clone(root: Path, arm: str, case: str, seed: int, new_name: str) -> Path:
    src = root / arm / case / f"s{seed}"
    dst = src.parent / new_name
    dst.mkdir()
    for f in src.iterdir():
        (dst / f.name).write_text(f.read_text())
    return dst


def test_frozen_population_is_five_arms_by_two_cases_with_planned_counts() -> None:
    frozen = aa.expected_population(ACQ)
    assert len(frozen.seeds) == 10
    assert sum(len(v) for v in frozen.seeds.values()) == 5 * (24 + 8)
    assert frozen.seeds[("A-port", "F")][961031] == 200 and frozen.seeds[("A-port", "P")][961011] == 150


def test_workflow_files_in_the_tree_are_those_of_the_acquisition_commit() -> None:
    """The tests read the working tree; the analysis reads `git show <acquisition commit>:<path>`. They must agree."""
    try:
        from_git = _REAL_EXPECTED(ACQ)
    except aa.InputError:
        pytest.skip("acquisition commit not in this clone")
    assert from_git.seeds == _REAL_EXPECTED(ACQ, sources=SOURCES).seeds


def test_abbreviated_or_malformed_acquisition_commit_is_refused() -> None:
    for bad in (ACQ[:12], "HEAD", "", ACQ.upper(), "--output=x"):
        with pytest.raises(aa.InputError, match="40-hex"):
            _REAL_EXPECTED(bad, sources=SOURCES)


def test_a_modified_run_list_is_refused_by_the_seed_checker() -> None:
    broken = {k: v.replace("100:961008 ", "", 1) for k, v in SOURCES.items()}
    with pytest.raises(aa.InputError, match="apples_seeds_check"):
        _REAL_EXPECTED(ACQ, sources=broken)


def test_complete_frozen_population_is_not_partial(a_identical: Path) -> None:
    c = contrast(analyse(a_identical), A_CONTRAST)
    assert c["partial"] is False and c["partial_reasons"] == [] and c["claims_withheld"] is False


def test_missing_seed_directory_makes_the_contrast_partial_and_withholds_every_claim(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    _drop(root, "A-port", "P", 961011)
    doc = analyse(root, attested=False)
    c = contrast(doc, A_CONTRAST)
    assert_no_confirmatory_claims(c)
    assert c["partial_reasons"] == [aa.NO_MANIFEST_REASON, "A-port P: 1 expected seed(s) missing: 961011"]
    assert doc["runs"]["A-port/P"]["n"] == 23 and doc["runs"]["A-port/P"]["expected"] == 24
    assert aa.main([str(root), "--parts", "A"]) == 0
    out = capsys.readouterr().out
    assert "**PARTIAL**" in out and "WITHHELD" in out and "Joint claim" not in out and "TRUE**" not in out
    assert "equivalent after Holm" not in out


def test_missing_f_seed_makes_the_contrast_partial(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    _drop(root, "A-up", "F", 960031)
    c = contrast(analyse(root, attested=False), A_CONTRAST)
    assert_no_confirmatory_claims(c)
    assert "960031" in c["partial_reasons"][1]


def test_missing_arm_or_case_directory_is_fatal(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    for f in (root / "A-port" / "F").rglob("*"):
        if f.is_file():
            f.unlink()
    for d in sorted((root / "A-port" / "F").iterdir()):
        d.rmdir()
    (root / "A-port" / "F").rmdir()
    assert aa.main([str(root), "--parts", "A"]) == 2
    assert "A-port/F" in capsys.readouterr().err.replace("\\", "/")


def test_same_count_substitution_with_an_unexpected_seed_is_partial(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    _drop(root, "A-port", "F", 961031)
    _clone(root, "A-port", "F", 961032, "s999999")  # still 8 directories
    assert len(list((root / "A-port" / "F").iterdir())) == N_RUNS
    doc = analyse(root, attested=False)
    c = contrast(doc, A_CONTRAST)
    assert_no_confirmatory_claims(c)
    text = " | ".join(c["partial_reasons"])
    assert "961031" in text and "s999999 (seed not in the frozen list)" in text
    assert doc["runs"]["A-port/F"]["n"] == 7  # the unexpected directory was not read as a run


def test_rejected_directory_is_never_read(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    stray = _clone(root, "A-port", "P", 961001, "s999999")
    (stray / "run.json").write_text("{definitely not json")
    (stray / "endpoints.json").unlink()
    assert "s999999" not in " ".join(aa.run_analysis(root, ("A",))[1]["dataset_fingerprint"].keys())
    aa_ledger = aa.Ledger(root)
    aa.load_arm_case(root, "A-port", "P", aa.expected_population(ACQ), aa_ledger)
    assert not any("s999999" in k for k in aa_ledger.files)


def test_leading_zero_alias_beside_the_real_seed_rejects_both_and_is_partial(tmp_path: Path) -> None:
    """The reviewer's probe: seven unique seeds plus s0960050, a copy of s960050, one intended run absent."""
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    _drop(root, "A-port", "F", 961038)
    _clone(root, "A-port", "F", 961031, "s0961031")
    assert len(list((root / "A-port" / "F").iterdir())) == N_RUNS
    doc = analyse(root, attested=False)
    c = contrast(doc, A_CONTRAST)
    assert_no_confirmatory_claims(c)
    text = " | ".join(c["partial_reasons"])
    assert "s0961031, s961031 (numeric identity 961031 repeated)" in text and "961038" in text
    assert doc["runs"]["A-port/F"]["n"] == 6  # neither copy of the repeated identity is used


def test_leading_zero_alias_alone_is_not_the_workflows_directory_name(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    (root / "A-port" / "F" / "s961031").rename(root / "A-port" / "F" / "s0961031")
    c = contrast(analyse(root, attested=False), A_CONTRAST)
    assert_no_confirmatory_claims(c)
    assert any("s0961031 (not written as s961031)" in r for r in c["partial_reasons"])
    assert any("961031" in r and "missing" in r for r in c["partial_reasons"])


def test_a_seed_expected_for_another_arm_is_unexpected_here(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True)
    _drop(root, "A-port", "F", 961031)
    foreign = root / "A-port" / "F" / "s960031"  # A-up's block, in A-port's case directory
    foreign.mkdir()
    for f in (root / "A-up" / "F" / "s960031").iterdir():
        (foreign / f.name).write_text(f.read_text())
    c = contrast(analyse(root, attested=False), A_CONTRAST)
    assert_no_confirmatory_claims(c)
    assert any("s960031 (seed not in the frozen list)" in r for r in c["partial_reasons"])


def test_only_the_contrasts_using_a_partial_arm_lose_their_claims(tmp_path: Path) -> None:
    root = build(tmp_path, ALL_ARMS, identical=True)
    _drop(root, "B-picc", "P", 964001)
    absent = {"arm": "B-picc", "case": "P", "seed": 964001, "energy": aa.expected_population(ACQ).seeds[("B-picc", "P")][964001],
              "outcome": "absent", "reason": "no directory: the job stopped", "files": []}
    attest(root, ("A", "B"), not_established=[absent])
    doc = analyse(root, ("A", "B"))
    assert contrast(doc, A_CONTRAST)["partial"] is False
    assert contrast(doc, A_CONTRAST)["claims"]["joint_claim_equivalent_on_all_endpoints"] is True
    assert contrast(doc, "B-pgcc vs B-up")["partial"] is False
    assert_no_confirmatory_claims(contrast(doc, "B-picc vs B-up"))
    d = contrast(doc, "B-pgcc vs B-picc (descriptive)")
    assert d["partial"] is True and d["claims"] == {} and d["claims_withheld"] is False


def test_extra_directory_that_is_not_a_seed_directory_stays_fatal(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    (root / "A-up" / "P" / "e100").mkdir()
    assert aa.main([str(root), "--parts", "A"]) == 2


# ------------------------------------------------------------------- identity and provenance (item 2)

SMOKE = {"mode": "smoke", "requested": 100_000, "simulated": 100_001}


def _only_first_of_each_cell(fn: Callable[[dict], None]) -> RunMutate:
    def go(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-port" and i == 0:
            fn(run)

    return go


@pytest.mark.parametrize(
    ("label", "edit"),
    [
        ("smoke mode at 1e5 histories", lambda r: r.update(SMOKE)),
        ("mode absent", lambda r: r.pop("mode")),
        ("mode full spelled FULL", lambda r: r.update(mode="FULL")),
        ("commit is another commit", lambda r: r.update(commit="0" * 40)),
        ("commit abbreviated", lambda r: r.update(commit=ACQ[:12])),
        ("commit absent", lambda r: r.pop("commit")),
        ("binary hash malformed", lambda r: r.update(binary_sha256="not-a-hash")),
        ("binary hash upper case", lambda r: r.update(binary_sha256=binary_of("A-port").upper())),
        ("binary hash absent", lambda r: r.pop("binary_sha256")),
        ("threads 99", lambda r: r.update(threads=99)),
        ("A arm run on B's 3 threads", lambda r: r.update(threads=3)),
        ("threads absent", lambda r: r.pop("threads")),
        ("simulated 1 against requested 1e7", lambda r: r.update(simulated=1)),
        ("simulated one below requested", lambda r: r.update(simulated=r["requested"] - 1)),
        ("simulated null", lambda r: r.update(simulated=None)),
        ("requested 1e5", lambda r: r.update(requested=100_000)),
        ("requested as a float", lambda r: r.update(requested=float(r["requested"]))),
        ("requested absent", lambda r: r.pop("requested")),
        ("arm absent", lambda r: r.pop("arm")),
        ("case absent", lambda r: r.pop("case")),
        ("seed absent", lambda r: r.pop("seed")),
        ("energy absent", lambda r: r.pop("energy_mev")),
        ("energy not the seed's energy", lambda r: r.update(energy_mev=150 if r["energy_mev"] == 100 else 100)),
        ("energy 70", lambda r: r.update(energy_mev=70)),
        ("energy as a string", lambda r: r.update(energy_mev=str(r["energy_mev"]))),
    ],
)
def test_run_json_probe_is_not_established_and_withholds_every_claim(
    tmp_path: Path, label: str, edit: Callable[[dict], None]
) -> None:
    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, run_mutate=_only_first_of_each_cell(edit)))
    c = contrast(doc, A_CONTRAST)
    assert_no_confirmatory_claims(c)
    assert doc["runs"]["A-port/P"]["unusable"], label
    # the affected cells' endpoints are not established, not silently estimated from the remaining runs
    assert sum(e["outcome"] == "not_established" for e in c["endpoints"]) >= 13


def test_the_reviewers_smoke_probe_one_observation_per_cell_is_refused(tmp_path: Path) -> None:
    """One smoke run per A-port cell, marked smoke at 1e5, and everything else full: the old code said joint TRUE."""
    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, run_mutate=_only_first_of_each_cell(lambda r: r.update(SMOKE))))
    assert len(doc["runs"]["A-port/P"]["unusable"]) == 3 and len(doc["runs"]["A-port/F"]["unusable"]) == 1
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))


def test_run_json_that_contradicts_its_directory_is_fatal(tmp_path: Path) -> None:
    for key, value in (("arm", "A-up"), ("case", "F"), ("seed", 1)):
        root = build(tmp_path / key, ("A-up", "A-port"), run_mutate=_only_first_of_each_cell(lambda r, k=key, v=value: r.update({k: v})))
        assert aa.main([str(root), "--parts", "A"]) == 2


def test_binary_hash_must_be_one_value_across_an_arm(tmp_path: Path) -> None:
    def other(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-up" and case == "F" and i == 0:
            run["binary_sha256"] = binary_of("something else")

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, run_mutate=other))
    c = contrast(doc, A_CONTRAST)
    assert_no_confirmatory_claims(c)
    assert any("binary_sha256 differs within arm A-up" in r for r in c["partial_reasons"])
    assert len(doc["runs"]["A-up/P"]["unusable"]) == 24 and len(doc["runs"]["A-up/F"]["unusable"]) == 8
    assert doc["runs"]["A-port/P"]["unusable"] == []


def test_an_arm_carrying_its_references_binary_is_partial(tmp_path: Path) -> None:
    def same(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-port":
            run["binary_sha256"] = binary_of("A-up")

    def same_f(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "A-port" and case == "F":
            rec["sha256"]["binary"] = binary_of("A-up")

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, run_mutate=same, doc_mutate=same_f))
    c = contrast(doc, A_CONTRAST)
    assert doc["runs"]["A-port/P"]["unusable"] == [] and doc["runs"]["A-port/F"]["unusable"] == []  # nothing else is wrong
    assert_no_confirmatory_claims(c)
    assert c["partial_reasons"] == [f"A-port and A-up carry the same binary_sha256 ({binary_of('A-up')})"]


def test_a_shared_binary_withholds_only_the_contrasts_between_those_two_arms(tmp_path: Path) -> None:
    def same(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "B-pgcc":
            run["binary_sha256"] = binary_of("B-picc")

    def same_f(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if arm == "B-pgcc" and case == "F":
            rec["sha256"]["binary"] = binary_of("B-picc")

    doc = analyse(build(tmp_path, ("B-up", "B-pgcc", "B-picc"), run_mutate=same, doc_mutate=same_f), ("B",))
    assert contrast(doc, "B-pgcc vs B-picc (descriptive)")["partial"] is True
    for name in ("B-pgcc vs B-up", "B-picc vs B-up"):
        assert contrast(doc, name)["partial"] is False and contrast(doc, name)["claims"]


def test_binary_hash_is_reported_per_arm(tmp_path: Path) -> None:
    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True))
    assert doc["runs"]["A-port/P"]["binary_sha256"] == [binary_of("A-port")]
    assert doc["runs"]["A-up/F"]["binary_sha256"] == [binary_of("A-up")]


@pytest.mark.parametrize(
    "edit",
    [
        lambda d: d.update(slab_depths_mm=[999]),
        lambda d: d.update(slab_depths_mm=[3, *d["slab_depths_mm"]]),
        lambda d: d.update(slab_depths_mm=list(reversed(d["slab_depths_mm"]))),
        lambda d: d.pop("slab_depths_mm"),
        lambda d: d.update(label="A-up_P_E100_N1e7_seed1"),
        lambda d: d.pop("label"),
        lambda d: d.update(layout="topas"),
    ],
)
def test_p_slab_depths_label_and_layout_must_match_the_frozen_energy(tmp_path: Path, edit: Callable[[dict], None]) -> None:
    def go(arm: str, case: str, energy: int | None, i: int, d: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 100 and i == 0:
            edit(d)

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, doc_mutate=go))
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))
    assert ep(doc, A_CONTRAST, "P100/R80")["outcome"] == "not_established"


def test_slab_depths_of_another_energy_do_not_move_the_run_to_that_energy(tmp_path: Path) -> None:
    """The old loader fell back to slab_depths_mm for placement; the frozen list now decides the energy."""

    def go(arm: str, case: str, energy: int | None, i: int, d: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 100 and i == 0:
            d["slab_depths_mm"] = list(aa.SLAB_DEPTHS[150])

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, doc_mutate=go))
    assert doc["runs"]["A-port/P"]["per_energy"] == {"100": 8, "150": 8, "200": 8}


@pytest.mark.parametrize(
    "edit",
    [
        lambda d: d.update(platform="A-up"),
        lambda d: d.update(seed=1),
        lambda d: d.update(commit="1" * 40),
        lambda d: d["sha256"].update(binary="0" * 64),
    ],
)
def test_f_record_top_level_is_reconciled_with_run_json(tmp_path: Path, edit: Callable[[dict], None]) -> None:
    def go(arm: str, case: str, energy: int | None, i: int, d: dict) -> None:
        if arm == "A-port" and case == "F" and i == 0:
            edit(d)

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, doc_mutate=go))
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))
    assert ep(doc, A_CONTRAST, "F/cax_127")["outcome"] == "not_established"


@pytest.mark.parametrize(
    ("label", "edit", "match"),
    [
        ("carries endpoint_status (mixed)", lambda d: d.update(endpoint_status="ok"), "mixes record schemas"),
        ("carries metrics_sha256 (mixed)", lambda d: d.update(metrics_sha256=aa.metrics_digest(d["metrics"])), "mixes"),
        ("carries cfg_sha256 (mixed)", lambda d: d.update(cfg_sha256=d["sha256"]["config"]), "mixes record schemas"),
        ("all three v2 fields at 2f9dab40", lambda d: d.update(cfg_sha256=d["sha256"]["config"], endpoint_status="ok",
                                                                metrics_sha256=aa.metrics_digest(d["metrics"])),
         "is bound to platform_study_record@2f9dab40"),
        ("metrics missing", lambda d: d.pop("metrics"), "lacks field"),
        ("metrics null", lambda d: d.update(metrics=None), "metrics is not an object"),
        ("a sha256 key missing", lambda d: d["sha256"].pop("binary"), "sha256 is not"),
        ("a sha256 key missing (HU)", lambda d: d["sha256"].pop("HU_Material"), "sha256 is not"),
        ("sha256 missing", lambda d: d.pop("sha256"), "lacks field"),
        ("materials malformed", lambda d: d.update(materials={"files": 1}), "materials"),
        ("an unexpected field", lambda d: d.update(note="x"), "does not write"),
    ],
)
def test_malformed_or_mixed_schema_f_record_refuses(
    tmp_path: Path, label: str, edit: Callable[[dict], None], match: str
) -> None:
    def go(arm: str, case: str, energy: int | None, i: int, d: dict) -> None:
        if arm == "A-port" and case == "F" and i == 0:
            edit(d)

    root = build(tmp_path, ("A-up", "A-port"), identical=True, doc_mutate=go)
    with pytest.raises(aa.InputError, match=match):
        analyse(root)
    assert aa.main([str(root), "--parts", "A"]) == 2, label


# ------------------------------------------------------------------- the v2 record schema (writers that emit it)


@pytest.fixture
def v2_bound(monkeypatch: pytest.MonkeyPatch) -> None:
    """A test-only binding of V2_COMMIT to the v2 schema (no acquisition is bound to it)."""
    monkeypatch.setattr(aa, "F_RECORD_BINDINGS", (*aa.F_RECORD_BINDINGS, V2_BINDING))


def test_v2_records_from_the_v2_writer_analyse_under_a_v2_binding(tmp_path: Path, v2_bound: None) -> None:
    root = build(tmp_path, ("A-up", "A-port"), identical=True, commit=V2_COMMIT, writer="v2")
    doc = analyse(root, commit=V2_COMMIT)
    c = contrast(doc, A_CONTRAST)
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is True
    assert doc["runs"]["A-up/F"]["record_bindings_used"] == [V2_BINDING.name]


@pytest.mark.parametrize(
    "edit",
    [
        lambda d: d["metrics"].update(cax_127=d["metrics"]["cax_127"] * 1.001),  # metrics no longer match their hash
        lambda d: d.update(metrics_sha256="0" * 64),
        lambda d: d.update(cfg_sha256="0" * 64),
    ],
)
def test_v2_record_is_strict(tmp_path: Path, v2_bound: None, edit: Callable[[dict], None]) -> None:
    def go(arm: str, case: str, energy: int | None, i: int, d: dict) -> None:
        if arm == "A-port" and case == "F" and i == 0:
            edit(d)

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, doc_mutate=go, commit=V2_COMMIT, writer="v2"),
                  commit=V2_COMMIT)
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))


@pytest.mark.parametrize("drop", sorted(aa.V2_ENDPOINT_FIELDS))
def test_v2_record_lacking_one_v2_field_is_mixed_and_refuses(tmp_path: Path, v2_bound: None, drop: str) -> None:
    def go(arm: str, case: str, energy: int | None, i: int, d: dict) -> None:
        if arm == "A-port" and case == "F" and i == 0:
            d.pop(drop)

    root = build(tmp_path, ("A-up", "A-port"), identical=True, doc_mutate=go, commit=V2_COMMIT, writer="v2")
    with pytest.raises(aa.InputError, match="mixes record schemas"):
        analyse(root, commit=V2_COMMIT)


def test_f_record_with_failed_endpoints_and_null_metrics_is_unusable(tmp_path: Path, v2_bound: None) -> None:
    def go(arm: str, case: str, energy: int | None, i: int, d: dict) -> None:
        if arm == "A-port" and case == "F" and i == 0:
            d.update(metrics=None, metrics_sha256=None, endpoint_status="error", endpoint_error="ValueError: x")

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, doc_mutate=go, commit=V2_COMMIT, writer="v2"),
                  commit=V2_COMMIT)
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))


def test_f_energy_70_probe(tmp_path: Path) -> None:
    def odd(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-port" and case == "F" and i == 0:
            run["energy_mev"] = 70

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, run_mutate=odd))
    assert_no_confirmatory_claims(contrast(doc, A_CONTRAST))
    assert ep(doc, A_CONTRAST, "F/r80_mm")["outcome"] == "not_established"


def test_part_b_threads_are_three_and_part_a_four(tmp_path: Path) -> None:
    doc = analyse(build(tmp_path, ALL_ARMS, identical=True), ("A", "B"))
    assert all(not c["partial"] for c in doc["contrasts"])
    assert {a: aa.ARM_THREADS[a] for a in ALL_ARMS} == {"A-up": 4, "A-port": 4, "B-up": 3, "B-pgcc": 3, "B-picc": 3}


# -------------------------------------------------------------------------------------- dataset fingerprint


def test_fingerprint_is_deterministic_and_covers_exactly_the_files_read(a_identical: Path) -> None:
    md1, doc1 = aa.run_analysis(a_identical, ("A",))
    _, doc2 = aa.run_analysis(a_identical, ("A",))
    fp = doc1["dataset_fingerprint"]
    assert fp == doc2["dataset_fingerprint"] and len(fp["sha256"]) == 64
    assert fp["n_files"] == 2 * (24 + 8) * 2  # two files per run directory, 64 runs
    assert any(fp["sha256"] in line for line in md1)
    (a_identical / "A-up" / "P" / ".provenance").mkdir()
    (a_identical / "A-up" / "P" / ".provenance" / "run_root.txt").write_text("not read by the analysis")
    assert aa.run_analysis(a_identical, ("A",))[1]["dataset_fingerprint"] == fp


def test_fingerprint_changes_with_any_byte_of_any_input(a_identical: Path) -> None:
    before = analyse(a_identical, attested=False)["dataset_fingerprint"]["sha256"]
    target = a_identical / "A-port" / "P" / "s961011" / "endpoints.json"
    original = target.read_text()
    target.write_text(original.replace("R80", "R81", 1) if "R80" in original else original + " ")
    changed = analyse(a_identical, attested=False)["dataset_fingerprint"]["sha256"]
    assert changed != before
    target.write_text(original)
    assert analyse(a_identical, attested=False)["dataset_fingerprint"]["sha256"] == before


def test_fingerprint_depends_on_which_path_holds_the_bytes(tmp_path: Path) -> None:
    root = build(tmp_path / "x", ("A-up", "A-port"), identical=True)
    before = analyse(root, attested=False)["dataset_fingerprint"]["sha256"]
    a, b = root / "A-port" / "P" / "s961001", root / "A-port" / "P" / "s961002"
    ta, tb = (a / "endpoints.json").read_text(), (b / "endpoints.json").read_text()
    (a / "endpoints.json").write_text(tb)
    (b / "endpoints.json").write_text(ta)
    assert analyse(root, attested=False)["dataset_fingerprint"]["sha256"] != before


# ------------------------------------------------------------------ domain bounds (item 3), unit and end to end


def _rec(energy: int = 100, **over: object) -> dict:
    rec = p_record(random.Random(1), energy, 1.0)
    rec.update(over)
    return rec


SIGMA_LO, SIGMA_HI = 1.0, 20.0


@pytest.mark.parametrize(
    ("value", "kept"),
    [
        (SIGMA_LO, True), (SIGMA_LO + 1e-9, True), (SIGMA_HI - 1e-9, True), (SIGMA_HI, True), (4.0, True),
        (math.nextafter(SIGMA_LO, 0.0), False), (SIGMA_LO - 1e-9, False), (math.nextafter(SIGMA_HI, 99.0), False),
        (SIGMA_HI + 1e-9, False), (0.0, False), (-0.0, False), (-1.0, False), (-4.0, False),
        (math.nan, False), (math.inf, False), (-math.inf, False), (None, False), ("4.0", False), (True, False),
    ],
)
def test_sigma_domain_is_the_producers_fit_contract(value: object, kept: bool) -> None:
    out = aa._p_values(_rec(sigma_40=value), 100)
    assert (out["P100/sigma_40"] is not None) is kept
    assert out["P100/sigma_60"] is not None  # only that endpoint


def test_sigma_bounds_are_the_producers() -> None:
    from pencil_endpoints import SIGMA_VALID_MM

    assert aa.SIGMA_BOUNDS_MM == SIGMA_VALID_MM == (SIGMA_LO, SIGMA_HI)


@pytest.mark.parametrize(
    ("value", "kept"),
    [
        (1.0, True), (1.0 - 1e-12, True), (0.5, True), (1e-300, True), (5e-324, True),
        (math.nextafter(1.0, 2.0), False), (1.0 + 1e-9, False), (1.01, False), (2.0, False),
        (0.0, False), (-0.0, False), (-1e-12, False), (-0.01, False),
        (math.nan, False), (math.inf, False), (-math.inf, False), (None, False), ("0.2", False), (False, False),
    ],
)
def test_ring_fraction_domain_is_0_to_1_and_zero_has_no_logarithm(value: object, kept: bool) -> None:
    out = aa._p_values(_rec(ring_40_20_40=value), 100)
    assert (out["P100/ring_40_20_40"] is not None) is kept
    assert out["P100/ring_40_40_80"] is not None


def test_out_of_domain_values_in_both_arms_are_not_established_end_to_end(tmp_path: Path) -> None:
    """The reviewer's probe: sigma=-1 mm and ring_5_10=1.01 in both arms used to give joint TRUE and 52 Holm decisions."""

    def bad(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if case == "P" and energy == 100:
            rec["sigma_40"], rec["ring_40_5_10"] = -1.0, 1.01

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, mutate=bad))
    c = contrast(doc, A_CONTRAST)
    assert ep(doc, A_CONTRAST, "P100/sigma_40")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "P100/ring_40_5_10")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "P100/sigma_60")["outcome"] == "equivalent"
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is False
    assert c["claims"]["n_not_established"] == 2
    assert ep(doc, A_CONTRAST, "P100/sigma_40")["holm_decision"] == "not_established"


def test_values_exactly_at_the_bounds_in_both_arms_are_established(tmp_path: Path) -> None:
    def edge(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if case == "P" and energy == 100:
            rec["sigma_40"], rec["sigma_60"], rec["ring_40_80_200"] = SIGMA_LO, SIGMA_HI, 1.0

    doc = analyse(build(tmp_path, ("A-up", "A-port"), identical=True, mutate=edge))
    for eid in ("P100/sigma_40", "P100/sigma_60", "P100/ring_40_80_200"):
        assert ep(doc, A_CONTRAST, eid)["outcome"] == "equivalent"


# ------------------------------------------------------------------------------ collection manifest (analysis side)


def _write_manifest(root: Path, *, skip: str | None = None, **override: object) -> None:
    manifest = attest(root, ("A",))
    manifest.update(override)
    manifest["files"] = [f for f in manifest["files"] if f["path"] != skip]
    (root / aa.MANIFEST_NAME).write_text(json.dumps(manifest))


def test_report_says_when_there_is_no_manifest(a_identical: Path) -> None:
    doc = analyse(a_identical, attested=False)
    assert doc["collection_manifest"]["present"] is False and "not attested" in doc["collection_manifest"]["note"]


def test_manifest_is_verified_when_present(a_identical: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write_manifest(a_identical)
    assert aa.main([str(a_identical), "--parts", "A"]) == 0
    assert "Collection manifest: verified, 128 file(s)" in capsys.readouterr().out
    doc = analyse(a_identical)
    assert doc["collection_manifest"]["present"] is True and doc["confirmatory_claims_possible"] is True


def test_file_changed_since_collection_is_refused(a_identical: Path) -> None:
    _write_manifest(a_identical)
    target = a_identical / "A-up" / "P" / "s960001" / "run.json"
    target.write_text(target.read_text().replace("HOST", "H0ST"))
    assert aa.main([str(a_identical), "--parts", "A"]) == 2


def test_file_read_but_not_in_the_manifest_is_refused(a_identical: Path) -> None:
    _write_manifest(a_identical, skip="A-up/P/s960001/run.json")
    assert aa.main([str(a_identical), "--parts", "A"]) == 2


def test_manifest_listing_a_missing_file_is_refused(a_identical: Path) -> None:
    _write_manifest(a_identical)
    (a_identical / "A-up" / "P" / "s960001" / "endpoints.json").unlink()
    assert aa.main([str(a_identical), "--parts", "A"]) == 2


def test_unreadable_manifest_is_refused(a_identical: Path) -> None:
    (a_identical / aa.MANIFEST_NAME).write_text("{")
    assert aa.main([str(a_identical), "--parts", "A"]) == 2


# ------------------------------------------- confirmatory claims need a verified collection (review 6935, item 2)


def test_flat_tree_without_a_manifest_makes_no_confirmatory_claim_anywhere(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The reviewer's probe: a complete, internally consistent flat tree with no manifest said joint equivalence TRUE."""
    root = build(tmp_path / "in", ALL_ARMS, identical=True)
    out = tmp_path / "result.json"
    assert aa.main([str(root), "--json", str(out), "--parts", "A,B"]) == 0
    md = capsys.readouterr().out
    doc = json.loads(out.read_text())
    assert doc["confirmatory_claims_possible"] is False and doc["claims_withheld_for_all"] == [aa.NO_MANIFEST_REASON]
    assert doc["collection_manifest"]["present"] is False
    for c in doc["contrasts"]:
        assert c["claims"] == {} and c["partial"] is True and aa.NO_MANIFEST_REASON in c["partial_reasons"]
        if c["kind"] == "confirmatory":
            assert_no_confirmatory_claims(c)
        for e in c["endpoints"]:
            assert e["outcome"] in ("descriptive", "not_established") and e["p_tost"] is None
            assert e["holm_decision"] == "" and e["holm_adjusted_p"] is None
    assert "Joint claim" not in md and "TRUE**" not in md and "equivalent after Holm" not in md
    assert "| equivalent |" not in md and "| not_equivalent |" not in md and "| inconclusive |" not in md
    assert "Confirmatory claims: WITHHELD for every contrast" in md
    # the estimates are still there: a preview, never confirmatory
    assert all(e["estimate"] is not None for c in doc["contrasts"] for e in c["endpoints"])


def test_the_same_tree_attested_is_confirmatory(tmp_path: Path) -> None:
    """Control for the test above: the only difference is the manifest."""
    root = build(tmp_path, ALL_ARMS, identical=True)
    doc = analyse(root, ("A", "B"))
    assert doc["confirmatory_claims_possible"] is True
    assert all(c["claims"]["joint_claim_equivalent_on_all_endpoints"] is True for c in doc["contrasts"] if c["kind"] == "confirmatory")


@pytest.mark.parametrize(
    ("label", "override", "match"),
    [
        ("another 40-hex acquisition commit", {"acquisition_commit": "1" * 40}, "is not the frozen commit in use"),
        ("an abbreviated acquisition commit", {"acquisition_commit": ACQ[:12]}, "is not the frozen commit in use"),
        ("no acquisition commit", {"acquisition_commit": None}, "is not the frozen commit in use"),
        ("parts B only", {"parts": ["B"]}, "do not cover the analysed part"),
        ("no parts", {"parts": []}, "do not cover the analysed part"),
        ("parts not a list", {"parts": "A"}, "is not a list of distinct parts"),
        ("schema 2", {"schema": 2}, "schema 2 is not supported"),
        ("schema 4", {"schema": 4}, "is not supported"),
        ("schema missing", {"schema": None}, "is not supported"),
        ("schema true", {"schema": True}, "is not supported"),
        ("another record binding", {"f_record_binding": {"name": "x", "schema": "platform_study_record@0f5ef7c"}},
         "f_record_binding"),
        ("no record binding", {"f_record_binding": None}, "f_record_binding"),
    ],
)
def test_manifest_that_does_not_attest_this_analysis_refuses_before_inference(
    a_identical: Path, label: str, override: dict, match: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write_manifest(a_identical, **override)
    read: list[str] = []
    real_load_run = aa.load_run
    monkeypatch.setattr(aa, "load_run", lambda *a, **k: read.append("run") or real_load_run(*a, **k))
    with pytest.raises(aa.InputError, match=match):
        aa.run_analysis(a_identical, ("A",))
    assert read == [], f"{label}: a run was read before the manifest was validated"
    assert aa.main([str(a_identical), "--parts", "A"]) == 2


def test_manifest_covering_more_parts_than_analysed_is_accepted(tmp_path: Path) -> None:
    root = build(tmp_path, ALL_ARMS, identical=True)
    attest(root, ("A", "B"))
    doc = analyse(root, ("A",))
    assert contrast(doc, A_CONTRAST)["claims"]["joint_claim_equivalent_on_all_endpoints"] is True
