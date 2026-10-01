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

Mutate = Callable[[str, str, int | None, int, dict], None]
RunMutate = Callable[[str, str, int | None, int, dict], None]

N_RUNS = 6
R80_BASE = {100: 77.5, 150: 158.0, 200: 259.0}
RING_FRAC = {(5, 10): 0.02, (10, 20): 0.05, (20, 40): 0.20, (40, 80): 0.40, (80, 200): 0.20}
ALL_ARMS = ("A-up", "A-port", "B-up", "B-pgcc", "B-picc")


def _noise(rng: random.Random, sd: float) -> float:
    return rng.gauss(0.0, sd)


def p_record(rng: random.Random, energy: int, scale: float) -> dict:
    rec: dict = {
        "label": "synthetic",
        "R80": R80_BASE[energy] + _noise(rng, 0.005 * scale),
        "R80_multiple_crossings": False,
        "slab_depths_mm": [3, *aa.SLAB_DEPTHS[energy]],
    }
    for d in aa.SLAB_DEPTHS[energy]:
        rec[f"sigma_{d}"] = 4.0 + _noise(rng, 0.001 * scale)
        for (lo, hi), frac in RING_FRAC.items():
            rec[f"ring_{d}_{lo}_{hi}"] = frac * math.exp(_noise(rng, 0.0005 * scale))
    return rec


def f_record(rng: random.Random, scale: float) -> dict:
    m: dict = {}
    for depth in (127, 201):
        for off in (5, 10, 20, 30):
            m[f"lateral_{depth}_{off}"] = 1.0 + _noise(rng, 0.01 * scale)
    for k in ("cax_127", "cax_201", "cax_i23"):
        m[k] = 1.0e-3 * math.exp(_noise(rng, 0.0005 * scale))
    m["r80_mm"], m["r20_mm"] = 250.0 + _noise(rng, 0.01 * scale), 260.0 + _noise(rng, 0.01 * scale)
    m["r80_crossings"] = m["r20_crossings"] = 1
    return m


def build(
    root: Path,
    arms: tuple[str, ...],
    *,
    n: int = N_RUNS,
    scale: float = 1.0,
    identical: bool = False,
    mutate: Mutate | None = None,
    run_mutate: RunMutate | None = None,
) -> Path:
    """Write a synthetic tree. `identical` gives every arm exactly the same values; mutate(arm, case, energy, i, rec)
    edits the endpoint record (case F: the metrics dict) and run_mutate edits run.json, before they are written."""
    for a_idx, arm in enumerate(arms):
        for case in aa.CASES:
            energies: list[int | None] = list(aa.SLAB_DEPTHS) if case == "P" else [None]
            for e_idx, energy in enumerate(energies):
                for i in range(n):
                    seed = 960000 + 1000 * a_idx + 100 * e_idx + (50 if case == "F" else 0) + i
                    rng = random.Random(7919 * (e_idx + 1) + 31 * i + (0 if identical else 104729 * (a_idx + 1)))
                    rec = p_record(rng, energy, scale) if energy is not None else f_record(rng, scale)
                    if mutate is not None:
                        mutate(arm, case, energy, i, rec)
                    run = {
                        "seed": seed,
                        "energy_mev": energy if energy is not None else 200,
                        "case": case,
                        "arm": arm,
                        "mode": "full",
                        "transport_status": "ok",
                    }
                    if run_mutate is not None:
                        run_mutate(arm, case, energy, i, run)
                    d = root / arm / case / f"s{seed}"
                    d.mkdir(parents=True)
                    (d / "run.json").write_text(json.dumps(run))
                    if case == "P":
                        (d / "endpoints.json").write_text("ENDPOINTS " + json.dumps(rec, sort_keys=True) + "\n")
                    else:
                        (d / "record.json").write_text(json.dumps({"study": "pe1", "metrics": rec}))
    return root


def analyse(root: Path, parts: tuple[str, ...] = ("A",)) -> dict:
    return aa.run_analysis(root, parts)[1]


def contrast(doc: dict, name: str) -> dict:
    return next(c for c in doc["contrasts"] if c["name"] == name)


def ep(doc: dict, name: str, eid: str) -> dict:
    return next(e for e in contrast(doc, name)["endpoints"] if e["endpoint"] == eid)


A_CONTRAST = "A-port vs A-up"


@pytest.fixture
def a_identical(tmp_path: Path) -> Path:
    return build(tmp_path, ("A-up", "A-port"), identical=True)


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
            half = min(-s.lo, s.hi) / math.sqrt(3.0)  # sample sd = 2 half / sqrt 3 -> se = margin / 3
            sgn = -1.0 if i % 2 == 0 else 1.0
            rec[s.key] = rec[s.key] * math.exp(sgn * half) if s.log else rec[s.key] + sgn * half

    doc = analyse(build(tmp_path, ("A-up", "A-port"), n=4, identical=True, mutate=spread))
    c = contrast(doc, A_CONTRAST)
    assert all(0.01 < e["p_tost"] < 0.05 for e in c["endpoints"])
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is True
    assert c["claims"]["n_equivalent_unadjusted"] == 52
    assert c["claims"]["n_equivalent_holm"] == 0
    assert {e["holm_decision"] for e in c["endpoints"]} == {"not shown"}


# ---------------------------------------------------------------------------------- invalid values


def test_run_with_transport_status_not_ok_makes_its_endpoints_not_established(tmp_path: Path) -> None:
    def fail(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-port" and case == "P" and energy == 100 and i == 0:
            run["transport_status"] = "failed"

    doc = analyse(build(tmp_path, ("A-up", "A-port"), run_mutate=fail))
    c = contrast(doc, A_CONTRAST)
    bad = {e["endpoint"] for e in c["endpoints"] if e["outcome"] == "not_established"}
    assert bad == {s.eid for s in aa.p_specs() if s.energy == 100}  # that energy's 13, nothing else
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is False
    assert c["claims"]["n_not_established"] == 13
    assert all(e["holm_decision"] == "not_established" for e in c["endpoints"] if e["endpoint"] in bad)
    assert doc["runs"]["A-port/P"]["n"] == 3 * N_RUNS  # the run stays in the population
    assert len(doc["runs"]["A-port/P"]["unusable"]) == 1


def test_f_run_not_ok_makes_all_13_f_endpoints_not_established(tmp_path: Path) -> None:
    def fail(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-up" and case == "F" and i == 2:
            run["transport_status"] = "timeout"

    doc = analyse(build(tmp_path, ("A-up", "A-port"), run_mutate=fail))
    bad = {e["endpoint"] for e in contrast(doc, A_CONTRAST)["endpoints"] if e["outcome"] == "not_established"}
    assert bad == {s.eid for s in aa.f_specs()}


def test_missing_transport_status_is_not_ok(tmp_path: Path) -> None:
    def drop(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-up" and case == "P" and energy == 200 and i == 0:
            del run["transport_status"]

    doc = analyse(build(tmp_path, ("A-up", "A-port"), run_mutate=drop))
    assert ep(doc, A_CONTRAST, "P200/R80")["outcome"] == "not_established"


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


def test_fewer_than_two_runs_in_an_arm_is_not_established(tmp_path: Path) -> None:
    def keep_one(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        pass

    root = build(tmp_path, ("A-up", "A-port"), mutate=keep_one)
    victims = sorted((root / "A-port" / "P").iterdir())[N_RUNS + 1 : 2 * N_RUNS]  # all but one 150 MeV run
    for v in victims:
        for f in v.iterdir():
            f.unlink()
        v.rmdir()
    doc = analyse(root)
    assert ep(doc, A_CONTRAST, "P150/R80")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "P100/R80")["outcome"] == "equivalent"


def test_unreadable_run_json_keeps_the_run_and_places_it_by_slab_depths(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    run = min((root / "A-port" / "P").iterdir())  # a 100 MeV run
    (run / "run.json").write_text("{not json")
    doc = analyse(root)
    assert doc["runs"]["A-port/P"]["n"] == 3 * N_RUNS
    assert len(doc["runs"]["A-port/P"]["unusable"]) == 1
    assert ep(doc, A_CONTRAST, "P100/R80")["outcome"] == "not_established"
    assert ep(doc, A_CONTRAST, "P150/R80")["outcome"] == "equivalent"


def test_missing_endpoint_file_is_an_unusable_run_not_a_crash(tmp_path: Path) -> None:
    root = build(tmp_path, ("A-up", "A-port"))
    (min((root / "A-up" / "F").iterdir()) / "record.json").unlink()
    doc = analyse(root)
    assert {e["outcome"] for e in contrast(doc, A_CONTRAST)["endpoints"] if e["endpoint"].startswith("F/")} == {
        "not_established"
    }


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


def test_exit_2_when_a_run_cannot_be_placed_at_an_energy(tmp_path: Path) -> None:
    def odd(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if arm == "A-up" and case == "P" and energy == 100 and i == 0:
            run["energy_mev"] = 70

    root = build(tmp_path, ("A-up", "A-port"), run_mutate=odd)
    assert aa.main([str(root), "--parts", "A"]) == 2


@pytest.mark.parametrize("parts", ["C", "A,A", ""])
def test_exit_2_on_bad_parts(tmp_path: Path, parts: str) -> None:
    assert aa.main([str(tmp_path), "--parts", parts]) == 2


def test_exit_2_when_root_is_not_a_directory(tmp_path: Path) -> None:
    assert aa.main([str(tmp_path / "nope")]) == 2
