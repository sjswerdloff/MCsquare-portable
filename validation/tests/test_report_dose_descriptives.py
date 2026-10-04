"""Tests for validation/report_dose_descriptives.py on SYNTHETIC arrays only: no acquisition data, no network."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Callable
from concurrent.futures import Executor, Future
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pencil_endpoints as pe
import platform_study_metrics as psm
import report_dose_descriptives as rd

SPEC_2_2 = rd.GammaSpec("g2_2_c10", 2.0, 2.0, 10.0, local=False)
SPEC_1_1 = rd.GammaSpec("g1_1_c10", 1.0, 1.0, 10.0, local=False)


def tent(n: int = 25, peak_depth: float = 10.0) -> np.ndarray:
    """Dose that rises and falls linearly through ``peak_depth`` mm: 1 - |depth - peak| / peak, depth bin k at k + 0.5 mm."""
    centres = np.arange(n) + 0.5
    return np.clip(1 - np.abs(centres - peak_depth) / peak_depth, 0, None)


class SerialExecutor(Executor):
    """Runs every job in the calling thread: spawned workers would not see the monkeypatched geometry, and numba's
    threaded kernels are not safe to call from several Python threads at once."""

    def __init__(self, max_workers: int = 1, initializer: Callable[[], None] | None = None) -> None:
        del max_workers, initializer

    def submit(self, fn: Callable[..., object], /, *args: object, **kwargs: object) -> Future[object]:  # type: ignore[override]
        future: Future[object] = Future()
        try:
            future.set_result(fn(*args, **kwargs))
        except BaseException as e:  # mirrors Executor semantics: the error is delivered by result()
            future.set_exception(e)
        return future


# ---------------------------------------------------------------------------------------------------------------------
# range metrics
# ---------------------------------------------------------------------------------------------------------------------


def test_range_levels_of_an_analytic_curve_have_the_known_answer() -> None:
    # distal of the maximum (0.95, at 9.5 and 10.5 mm) the curve is 1 - (d - 10) / 10, so level L is crossed at
    # d = 10 + 10 (1 - L) exactly, and linear interpolation between voxel centres reproduces it
    idd = tent()
    peak = idd.max()
    assert rd.range_at_level(idd, 0.9) == (pytest.approx(10 + 10 * (1 - 0.9 * peak), abs=1e-12), 1)
    assert rd.range_at_level(idd, 0.8) == (pytest.approx(10 + 10 * (1 - 0.8 * peak), abs=1e-12), 1)
    assert rd.range_at_level(idd, 0.2) == (pytest.approx(10 + 10 * (1 - 0.2 * peak), abs=1e-12), 1)


def test_the_levels_are_ordered_and_the_falloff_is_positive() -> None:
    m = rd.range_metrics(tent())
    assert m.r90 is not None
    assert m.r80 is not None
    assert m.r20 is not None
    assert m.r90 < m.r80 < m.r20
    falloff = m.value("falloff")
    assert falloff is not None
    assert falloff > 0
    assert falloff == pytest.approx(m.r20 - m.r80, abs=1e-12)
    assert falloff == pytest.approx(10 * (0.8 - 0.2) * tent().max(), abs=1e-12)


def test_a_sample_exactly_at_the_level_is_the_crossing_point() -> None:
    # pencil_endpoints.r80 takes idd[i] >= level > idd[i + 1]: a voxel centre exactly at the level is the interpolated depth
    idd = np.array([0.0, 1.0, 0.8, 0.5])
    assert rd.range_at_level(idd, 0.8) == (2.5, 1)
    assert pe.r80(idd) == (2.5, False)


def test_a_level_with_two_distal_crossings_is_counted_and_unusable() -> None:
    idd = np.array([0.0, 1.0, 0.5, 0.9, 0.0])  # 0.8 is crossed between bins 1-2 and again between 3-4
    depth, n = rd.range_at_level(idd, 0.8)
    assert n == 2
    assert depth == pytest.approx(1.5 + (1.0 - 0.8) / (1.0 - 0.5))  # the FIRST crossing is the one returned
    m = rd.range_metrics(idd)
    assert m.n80 == 2
    assert m.value("R80") is None  # never reported as if it were a single crossing
    assert m.value("falloff") is None


def test_a_level_never_reached_is_absent_not_imputed() -> None:
    idd = np.array([0.0, 1.0, 0.5, 0.4, 0.35])  # the 20% level (0.2) is never crossed distally
    assert rd.range_at_level(idd, 0.2) == (None, 0)
    m = rd.range_metrics(idd)
    assert m.r20 is None
    assert m.value("R20") is None
    assert m.value("falloff") is None
    assert m.value("R80") is not None  # an unrelated level is unaffected
    rising = np.array([0.0, 0.2, 0.5, 1.0])  # maximum in the last bin: nothing distal to it
    assert rd.range_at_level(rising, 0.8) == (None, 0)


def test_range_at_level_08_is_bit_identical_to_pencil_endpoints_r80() -> None:
    rng = np.random.default_rng(7)
    for _ in range(50):
        k = np.arange(350) + 0.5
        peak = rng.uniform(60, 300)
        idd = (
            np.where(k < peak, 0.3 + 0.7 * (k / peak) ** 3, np.clip(1 - (k - peak) / 5, 0, None)) * rng.uniform(1e6, 1e9)
        ).astype(np.float32)
        idd *= (1 + 1e-4 * rng.standard_normal(idd.size)).astype(np.float32)
        mine, n = rd.range_at_level(idd, 0.8)
        theirs, multiple = pe.r80(idd)
        assert mine == theirs
        assert (n > 1) == multiple


def test_cell_with_an_unusable_run_is_unavailable_in_both_directions() -> None:
    ok = [10.0, 10.1, 9.9, 10.2]
    assert rd.contrast_cell(ok, ok)["available"] is True
    for ev, ref in ((ok[:3] + [None], ok), (ok, [None, *ok[1:]])):
        cell = rd.contrast_cell(ev, ref)
        assert cell["available"] is False
        assert "estimate" not in cell


# ---------------------------------------------------------------------------------------------------------------------
# Welch statistics
# ---------------------------------------------------------------------------------------------------------------------


def test_welch_matches_scipy_and_the_direction_is_evaluated_minus_reference() -> None:
    rng = np.random.default_rng(3)
    ev, ref = rng.normal(5.002, 0.02, 8), rng.normal(5.0, 0.05, 8)
    w = rd.welch_difference(ev, ref)
    t = stats.ttest_ind(ev, ref, equal_var=False)
    ci = t.confidence_interval(0.95)
    assert w["estimate"] == pytest.approx(ev.mean() - ref.mean(), abs=1e-14)
    assert w["df"] == pytest.approx(t.df, rel=1e-12)
    assert w["ci95"][0] == pytest.approx(ci.low, abs=1e-12)
    assert w["ci95"][1] == pytest.approx(ci.high, abs=1e-12)
    assert w["se"] == pytest.approx((ev.var(ddof=1) / 8 + ref.var(ddof=1) / 8) ** 0.5, rel=1e-12)
    assert (w["n_evaluated"], w["n_reference"]) == (8, 8)


def test_welch_with_zero_variance_has_no_df_or_interval_and_refuses_too_few_runs() -> None:
    w = rd.welch_difference([1.0, 1.0, 1.0], [1.5, 1.5, 1.5])
    assert w["estimate"] == pytest.approx(-0.5)
    assert w["se"] == 0.0
    assert w["df"] is None
    assert w["ci95"] is None
    with pytest.raises(rd.DescriptiveError):
        rd.welch_difference([1.0], [1.0, 2.0])
    with pytest.raises(rd.DescriptiveError):
        rd.welch_difference([1.0, float("nan")], [1.0, 2.0])


def test_range_contrasts_use_the_per_run_falloff_and_the_reference_is_upstream() -> None:
    per_run: dict[tuple[str, int, int], rd.RangeMetrics] = {}
    rng = np.random.default_rng(11)
    for arm in rd.ARM_ORDER:
        for energy in rd.ENERGIES:
            for seed in range(1, 9):
                r80 = 100.0 + (0.01 if arm in ("apt", "bpg", "bpi") else 0.0) + rng.normal(0, 0.002)
                per_run[(arm, energy, seed)] = rd.RangeMetrics(r80 - 5, r80, r80 + 6 + rng.normal(0, 0.01), 1, 1, 1)
    out = rd.range_contrasts(per_run)
    cell = out["B1"]["150"]["falloff"]
    ev = [per_run[("bpg", 150, s)].value("falloff") for s in range(1, 9)]
    ref = [per_run[("bup", 150, s)].value("falloff") for s in range(1, 9)]
    assert cell["estimate"] == pytest.approx(np.mean(ev) - np.mean(ref), abs=1e-13)
    assert out["A"]["100"]["R80"]["estimate"] == pytest.approx(0.01, abs=0.01)
    assert set(out) == {"A", "B1", "B2"}
    assert set(out["A"]) == {"100", "150", "200"}
    assert set(out["A"]["100"]) == {"R90", "R80", "R20", "falloff"}


# ---------------------------------------------------------------------------------------------------------------------
# mean dose and the split-half partition
# ---------------------------------------------------------------------------------------------------------------------


def test_mean_dose_of_known_arrays_without_renormalisation() -> None:
    a = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
    b = 3 * a  # a run with three times the total must count three times, not be rescaled to the first run's total
    m = rd.mean_dose([a, b])
    assert m.dtype == np.float64
    np.testing.assert_array_equal(m, 2 * a.astype(np.float64))
    assert a.dtype == np.float32
    np.testing.assert_array_equal(a, [[1, 2], [3, 4]])  # inputs are not modified


def test_mean_dose_refuses_nothing_and_mismatched_shapes() -> None:
    with pytest.raises(rd.DescriptiveError):
        rd.mean_dose([])
    with pytest.raises(rd.DescriptiveError):
        rd.mean_dose([np.zeros((2, 2)), np.zeros((2, 3))])


def test_split_half_picks_the_four_lowest_and_the_four_highest_seeds() -> None:
    low, high = rd.split_half([964008, 964003, 964001, 964007, 964005, 964002, 964006, 964004])  # deliberately unordered
    assert low == [964001, 964002, 964003, 964004]
    assert high == [964005, 964006, 964007, 964008]


@pytest.mark.parametrize("seeds", [[1, 2, 3], [1], [], [1, 1, 2, 3]])
def test_split_half_refuses_an_odd_empty_or_repeated_seed_set(seeds: list[int]) -> None:
    with pytest.raises(rd.DescriptiveError):
        rd.split_half(seeds)


def test_accumulate_arm_builds_the_halves_from_the_seeds_not_the_load_order(tmp_path: Path) -> None:
    seeds = [50, 20, 80, 10, 70, 40, 60, 30]
    runs = [rd.Run("aup", "P", 100, s, tmp_path / f"s{s}") for s in seeds]  # presented out of seed order
    arm = rd.accumulate_arm(runs, lambda r: np.full((2, 2, 2), float(r.seed), dtype=np.float32), halves=True)
    assert arm.seeds == sorted(seeds)
    assert arm.low is not None
    assert arm.high is not None
    assert np.all(arm.low == np.mean([10, 20, 30, 40]))
    assert np.all(arm.high == np.mean([50, 60, 70, 80]))
    assert np.all(arm.mean == np.mean(seeds))
    no_halves = rd.accumulate_arm(runs, lambda _r: np.zeros((1, 1, 1)), halves=False)
    assert no_halves.low is None
    assert no_halves.high is None


# ---------------------------------------------------------------------------------------------------------------------
# binding controls
# ---------------------------------------------------------------------------------------------------------------------


def test_r80_control_accepts_agreement_and_returns_the_difference() -> None:
    assert rd.check_r80_control(recorded=77.5, computed=77.5, recorded_multiple=False, computed_crossings=1) == 0.0
    assert rd.check_r80_control(
        recorded=77.5, computed=77.5 + 5e-10, recorded_multiple=False, computed_crossings=1
    ) == pytest.approx(5e-10)
    assert rd.check_r80_control(recorded=None, computed=None, recorded_multiple=False, computed_crossings=0) == 0.0


@pytest.mark.parametrize(
    ("recorded", "computed", "multiple", "crossings"),
    [
        (77.5, 77.5 + 2e-9, False, 1),  # just over the 1e-9 mm tolerance
        (77.5, 77.6, False, 1),
        (77.5, None, False, 0),  # the recomputation lost the crossing
        (None, 77.5, False, 1),  # the record had none
        (77.5, 77.5, True, 1),  # the record flagged several crossings, the recomputation found one
        (77.5, 77.5, False, 2),  # the other way round
    ],
)
def test_r80_control_refuses_any_mismatch(
    recorded: float | None, computed: float | None, multiple: bool, crossings: int
) -> None:
    with pytest.raises(rd.ControlError):
        rd.check_r80_control(recorded=recorded, computed=computed, recorded_multiple=multiple, computed_crossings=crossings)


def field_metrics(**over: object) -> dict[str, object]:
    base: dict[str, object] = {k: 1.0 + i for i, k in enumerate(rd.FIELD_METRIC_KEYS)}
    base |= {"r80_crossings": 1, "r20_crossings": 1}
    return {**base, **over}


def test_field_control_accepts_a_match_and_reports_the_largest_relative_difference() -> None:
    rec = field_metrics()
    assert rd.check_field_control(rec, dict(rec)) == 0.0
    new = {**rec, "cax_127": float(rec["cax_127"]) * (1 + 5e-13)}  # type: ignore[arg-type]
    assert rd.check_field_control(rec, new) == pytest.approx(5e-13, rel=1e-3)


def test_field_control_has_exactly_the_13_endpoint_values() -> None:
    assert len(rd.FIELD_METRIC_KEYS) == 13
    assert len(set(rd.FIELD_METRIC_KEYS)) == 13


@pytest.mark.parametrize(
    "change",
    [
        {"cax_127": 2.0 * (1 + 1e-9)},  # relative difference above 1e-12
        {"lateral_127_5": None},  # a number became None
        {"r80_crossings": 2},  # a non-endpoint key must be identical
        {"cax_201": float("nan")},
    ],
)
def test_field_control_refuses_a_mismatch(change: dict[str, object]) -> None:
    rec = field_metrics(cax_127=2.0)
    with pytest.raises(rd.ControlError):
        rd.check_field_control(rec, {**rec, **change})


def test_field_control_requires_none_to_be_exact_and_the_key_sets_to_match() -> None:
    rec = field_metrics(r80_mm=None)
    assert rd.check_field_control(rec, dict(rec)) == 0.0
    with pytest.raises(rd.ControlError):
        rd.check_field_control(rec, {**rec, "r80_mm": 250.0})
    with pytest.raises(rd.ControlError):
        rd.check_field_control(rec, {k: v for k, v in rec.items() if k != "cax_127"})
    with pytest.raises(rd.ControlError):
        rd.check_field_control(
            {k: v for k, v in rec.items() if k != "cax_127"}, {k: v for k, v in rec.items() if k != "cax_127"}
        )


def test_the_case_f_endpoint_function_emits_every_key_the_control_compares() -> None:
    n = psm.N
    out = psm.endpoints(np.ones((n, n, n), dtype=np.float32), rd.F_SIDE_MM)
    assert set(rd.FIELD_METRIC_KEYS) <= set(out)
    assert {"r80_crossings", "r20_crossings"} <= set(out)


def stage1_result(**over: object) -> dict[str, object]:
    idd = tent(350, 100.0).astype(np.float32)
    rm = rd.range_metrics(idd)
    base = {
        "arm": "aup", "case": "P", "energy": 100, "seed": 1, "idd": idd, "range": rd.asdict(rm), "recorded_r80": rm.r80,
        "recorded_multiple": False, "pencil_endpoints_r80": pe.r80(idd)[0],
        "files": [{"path": "/x/Dose.raw", "bytes": 1, "sha256": "ab", "sha256_recorded": "ab"}],
    }  # fmt: skip
    return {**base, **over}


def test_apply_controls_reports_what_was_checked() -> None:
    summary = rd.apply_controls([stage1_result(), stage1_result(seed=2)])
    assert summary["case_P_r80"]["runs_checked"] == 2
    assert summary["case_P_r80"]["max_abs_diff_mm"] == 0.0
    assert summary["dose_file_sha256"]["files_checked"] == 2
    assert summary["case_F_metrics"]["runs_checked"] == 0


def test_apply_controls_refuses_a_wrong_recorded_r80_a_wrong_hash_and_a_range_function_that_drifts() -> None:
    good = stage1_result()
    with pytest.raises(rd.ControlError):
        rd.apply_controls([good, stage1_result(seed=2, recorded_r80=float(good["recorded_r80"]) + 1e-6)])  # type: ignore[arg-type]
    bad_hash = stage1_result(files=[{"path": "/x/Dose.raw", "bytes": 1, "sha256": "ab", "sha256_recorded": "cd"}])
    with pytest.raises(rd.ControlError):
        rd.apply_controls([bad_hash])
    with pytest.raises(rd.ControlError):
        rd.apply_controls([stage1_result(pencil_endpoints_r80=100.0)])


# ---------------------------------------------------------------------------------------------------------------------
# gamma
# ---------------------------------------------------------------------------------------------------------------------


def gaussian_profile(x: np.ndarray, centre: float = 0.0, sigma: float = 10.0) -> np.ndarray:
    return np.exp(-((x - centre) ** 2) / (2 * sigma**2))


def gamma_1d(ref: np.ndarray, ev: np.ndarray, x: np.ndarray, spec: rd.GammaSpec) -> tuple[np.ndarray, dict[str, object]]:
    norm = float(ref.max())
    g = rd.compute_gamma((x,), ref, ev, spec, norm)
    return g, rd.gamma_summary(g, ref, rd.cutoff_value(spec, norm))


def test_identical_distributions_pass_at_100_percent_with_zero_gamma() -> None:
    x = np.arange(-50, 50.5, 0.5)
    ref = gaussian_profile(x)
    g, s = gamma_1d(ref, ref.copy(), x, SPEC_2_2)
    assert s["pass_rate_percent"] == 100.0
    assert s["mean_gamma"] == 0.0
    assert s["n_evaluated"] == int((ref >= 0.1).sum()) > 0
    assert s["n_not_evaluated_in_region"] == 0
    assert np.isnan(g[ref < 0.1]).all()


def test_a_shift_beyond_the_distance_criterion_with_a_dose_change_beyond_the_dose_criterion_fails_where_expected() -> None:
    x = np.arange(-50, 50.5, 0.5)
    ref = gaussian_profile(x)
    ev = 1.05 * gaussian_profile(x, centre=3.0)  # 3 mm > 2 mm and 5% > 2%
    g, s = gamma_1d(ref, ev, x, SPEC_2_2)
    assert s["pass_rate_percent"] < 100.0  # type: ignore[operator]
    assert s["n_fail"] > 0
    steep = int(np.argmin(np.abs(x - 10.0)))  # slope 0.06 /mm: the 3 mm shift is a 0.18 dose difference
    assert g[steep] > 1.0
    # the same profile shifted by LESS than the criterion, with a dose change within it, passes everywhere
    _, ok = gamma_1d(ref, 1.01 * gaussian_profile(x, centre=1.0), x, SPEC_2_2)
    assert ok["pass_rate_percent"] == 100.0
    # tightening the criteria to 1%/1 mm turns that passing comparison into a failing one
    _, tight = gamma_1d(ref, 1.01 * gaussian_profile(x, centre=1.0), x, SPEC_1_1)
    assert tight["pass_rate_percent"] < 100.0  # type: ignore[operator]


def brute_force_gamma_1d(x: np.ndarray, ref: np.ndarray, ev: np.ndarray, spec: rd.GammaSpec, norm: float) -> np.ndarray:
    """Global gamma by dense search: the independent oracle for the pymedphys call."""
    fine = np.arange(x[0], x[-1] + 1e-9, 0.01)
    ev_fine = np.interp(fine, x, ev)
    out = np.full(x.size, np.nan)
    for i in np.flatnonzero(ref >= spec.cutoff_percent / 100 * norm):
        dose_term = (ev_fine - ref[i]) / (spec.dose_percent / 100 * norm)
        dist_term = (fine - x[i]) / spec.distance_mm
        out[i] = min(np.sqrt(dose_term**2 + dist_term**2).min(), rd.MAX_GAMMA)
    return out


def test_pymedphys_gamma_agrees_with_a_brute_force_search() -> None:
    x = np.arange(-50, 50.5, 0.5)
    ref = gaussian_profile(x)
    ev = 1.03 * gaussian_profile(x, centre=1.7, sigma=10.5)
    g, _ = gamma_1d(ref, ev, x, SPEC_2_2)
    oracle = brute_force_gamma_1d(x, ref, ev, SPEC_2_2, float(ref.max()))
    inside = np.isfinite(oracle)
    assert inside.sum() > 20
    assert np.isfinite(g[inside]).all()
    # pymedphys searches a 0.2 mm grid that is a subset of the oracle's 0.01 mm grid: never below the oracle, and above it by
    # at most the coarser sampling can cost (half a step in distance is 0.05 of the 2 mm criterion; a little more with the
    # dose term)
    assert (g[inside] >= oracle[inside] - 1e-9).all()
    assert (g[inside] - oracle[inside]).max() < 0.1
    far_from_one = np.abs(oracle[inside] - 1) > 0.1
    assert far_from_one.sum() > 10
    assert ((g[inside] <= 1) == (oracle[inside] <= 1))[far_from_one].all()


def test_the_low_dose_cutoff_excludes_points_and_changes_the_count() -> None:
    x = np.arange(-50, 50.5, 0.5)
    ref = gaussian_profile(x)
    counts = {}
    for cutoff in (1.0, 10.0, 50.0):
        spec = rd.GammaSpec(f"c{cutoff}", 2.0, 2.0, cutoff, local=False)
        g, s = gamma_1d(ref, ref.copy(), x, spec)
        counts[cutoff] = s["n_evaluated"]
        assert s["n_evaluated"] == int((ref >= cutoff / 100).sum())
        assert s["n_above_cutoff"] == s["n_evaluated"]
        assert np.isnan(g[ref < cutoff / 100]).all()
    assert counts[1.0] > counts[10.0] > counts[50.0] > 0


def test_local_normalisation_fails_a_low_dose_error_that_global_normalisation_passes() -> None:
    x = np.arange(0.0, 100.0, 1.0)
    ref = np.where(x < 50, 100.0, 5.0)
    ev = np.where(x < 50, 100.0, 5.2)  # +4% of the local dose, 0.2% of the maximum
    glob = rd.GammaSpec("g", 2.0, 2.0, 1.0, local=False)
    loc = rd.GammaSpec("l", 2.0, 2.0, 1.0, local=True)
    _, sg = gamma_1d(ref, ev, x, glob)
    _, sl = gamma_1d(ref, ev, x, loc)
    assert sg["pass_rate_percent"] == 100.0
    assert sl["pass_rate_percent"] == 50.0
    assert sl["n_evaluated"] == sg["n_evaluated"] == 100


def test_gamma_above_the_maximum_is_capped_at_it_and_counted() -> None:
    x = np.arange(-50, 50.5, 0.5)
    ref = gaussian_profile(x)
    g, s = gamma_1d(ref, 3.0 * ref, x, SPEC_2_2)  # a 200% dose error: the true gamma is far above 2
    assert np.nanmax(g) == rd.MAX_GAMMA
    assert s["n_gamma_capped"] > 0
    assert s["pass_rate_percent"] == 0.0
    assert s["mean_gamma"] <= rd.MAX_GAMMA


def test_depth_shift_moves_toward_larger_index_without_wrapping() -> None:
    a = np.arange(1.0, 6.0)
    np.testing.assert_array_equal(rd.depth_shifted(a, 2), [0, 0, 1, 2, 3])
    b = np.ones((4, 2, 2))
    shifted = rd.depth_shifted(b, 1)
    assert shifted[0].sum() == 0
    assert shifted[1:].sum() == 12
    for bad in (0, -1, 4, 9):
        with pytest.raises(rd.DescriptiveError):
            rd.depth_shifted(a if bad in (0, -1) else b[:4], bad)


def summary_with_rate(rate: float | None) -> dict[str, dict[str, object]]:
    return {"shift_depth_3mm": {"g2_2_c10": {"pass_rate_percent": rate}}}


def test_the_sensitivity_control_accepts_detected_perturbations_and_refuses_a_gamma_that_cannot_fail() -> None:
    ok = {"100": {"idd": summary_with_rate(91.0), "3d": summary_with_rate(97.5)}, "150": {"idd": summary_with_rate(80.0)}}
    assert rd.check_sensitivity(ok) == 3
    for bad in (100.0, None):
        with pytest.raises(rd.ControlError):
            rd.check_sensitivity({"100": {"idd": summary_with_rate(91.0), "3d": summary_with_rate(bad)}})
    with pytest.raises(rd.ControlError):
        rd.check_sensitivity({})


def test_a_real_gamma_call_detects_the_3_mm_shift_on_a_bragg_like_curve() -> None:
    x = np.arange(60) + 0.5
    idd = np.where(x < 40, 0.3 + 0.7 * (x / 40) ** 3, np.clip(1 - (x - 40) / 3, 0, None))
    for shift, expect_failures in ((3, True), (1, False)):
        shifted = rd.depth_shifted(idd, shift)
        _, s = gamma_1d(idd, shifted, x, SPEC_2_2)
        assert (s["pass_rate_percent"] < 100.0) is expect_failures  # type: ignore[operator]


def test_gamma_summary_counts_pass_fail_not_evaluated_and_capped_points() -> None:
    gamma = np.array([np.nan, 0.5, 1.0, 1.5, 2.0, np.inf])
    reference = np.array([1.0, 10, 10, 10, 10, 10])
    s = rd.gamma_summary(gamma, reference, cutoff=5.0)
    assert s["n_above_cutoff"] == 5
    assert s["n_evaluated"] == 4
    assert s["n_not_evaluated_in_region"] == 1  # the inf inside the region
    assert s["n_pass"] == 2  # 0.5 and exactly 1.0
    assert s["n_fail"] == 2
    assert s["pass_rate_percent"] == 50.0
    assert s["mean_gamma"] == pytest.approx(1.25)
    assert s["n_gamma_capped"] == 1
    empty = rd.gamma_summary(np.array([np.nan, np.nan]), np.array([1.0, 1.0]), cutoff=5.0)
    assert empty["n_evaluated"] == 0
    assert empty["pass_rate_percent"] is None


def test_gamma_summary_refuses_a_finite_gamma_below_the_cutoff() -> None:
    with pytest.raises(rd.GammaRegionError):
        rd.gamma_summary(np.array([0.5, 0.5]), np.array([1.0, 10.0]), cutoff=5.0)


def test_crop_slices_is_the_padded_bounding_box_clipped_to_the_array() -> None:
    a = np.zeros((20, 20))
    a[8:10, 12:13] = 1.0
    box = rd.crop_slices(a, 0.5, 3)
    assert box == (slice(5, 13), slice(9, 16))
    assert rd.crop_slices(a, 0.5, 100) == (slice(0, 20), slice(0, 20))
    assert rd.crop_slices(np.array([0.0, 1.0, 0.0]), 0.5, 0) == (slice(1, 2),)
    with pytest.raises(rd.GammaRegionError):
        rd.crop_slices(a, 2.0, 3)


def blob(shape: tuple[int, int, int], centre: tuple[float, float, float], sigma: float) -> np.ndarray:
    grids = np.meshgrid(*(np.arange(n, dtype=float) for n in shape), indexing="ij")
    r2 = sum((g - c) ** 2 for g, c in zip(grids, centre, strict=True))
    return np.exp(-r2 / (2 * sigma**2))


def test_cropping_cannot_change_a_3d_gamma_map() -> None:
    shape = (44, 44, 44)
    axes = tuple(np.arange(n, dtype=float) for n in shape)
    ref = blob(shape, (22, 22, 22), 2.0)
    ev = 1.03 * blob(shape, (22.7, 21.5, 22.2), 2.1)
    full = rd.compute_gamma(axes, ref, ev, SPEC_2_2, float(ref.max()))
    job = rd.prepare_gamma_job(axes, ref, ev, SPEC_2_2, spacing_mm=1.0, keep_map=True)
    assert job.reference.size < ref.size  # the test is only meaningful if it really crops
    res = rd.run_gamma_job(job)
    np.testing.assert_array_equal(res["map"], full)
    assert res["n_evaluated"] == int(np.isfinite(full).sum())
    uncropped = rd.run_gamma_job(rd.prepare_gamma_job(axes, ref, ev, SPEC_2_2, spacing_mm=1.0, crop=False, keep_map=True))
    check = rd.summarise_crop_check(res, uncropped)
    assert check["identical_gamma_maps"] is True
    assert check["identical_summary"] is True


def test_the_gamma_job_records_the_exact_parameters_and_returns_the_requested_plane() -> None:
    shape = (30, 30, 30)
    axes = tuple(np.arange(n, dtype=float) for n in shape)
    ref = blob(shape, (15, 15, 15), 2.0)
    res = rd.run_gamma_job(rd.prepare_gamma_job(axes, ref, ref.copy(), SPEC_1_1, spacing_mm=1.0, plane=(1, 15)))
    assert res["params"] == {
        "dose_percent_threshold": 1.0, "distance_mm_threshold": 1.0, "lower_percent_dose_cutoff": 10.0,
        "lower_dose_cutoff_absolute": pytest.approx(0.1), "local_gamma": False, "global_normalisation": pytest.approx(1.0),
        "interp_fraction": 10, "max_gamma": 2.0, "skip_once_passed": False, "random_subset": None, "interp_algo": "pymedphys",
        "ram_available_bytes": rd.GAMMA_RAM_BYTES,
    }  # fmt: skip
    assert res["crop"]["full_shape"] == [30, 30, 30]
    assert res["plane"].shape == (30, 30)  # (axis 0, axis 2) of the uncropped volume
    assert np.isfinite(res["plane"][15, 15])
    assert np.isnan(res["plane"][0, 0])  # below the cutoff, outside the crop
    assert res["pass_rate_percent"] == 100.0


def test_prepare_gamma_job_refuses_a_pad_shorter_than_the_search_and_mismatched_grids() -> None:
    axes = (np.arange(40, dtype=float),)
    ref = gaussian_profile(axes[0], 20.0, 3.0)
    wide = rd.GammaSpec("w", 2.0, 6.0, 10.0, local=False)  # 6 mm x max_gamma 2 = 12 mm search > the 10 mm pad
    with pytest.raises(rd.GammaRegionError):
        rd.prepare_gamma_job(axes, ref, ref, wide, spacing_mm=1.0)
    with pytest.raises(rd.GammaRegionError):
        rd.prepare_gamma_job(axes, ref, ref[:-1], SPEC_2_2, spacing_mm=1.0)
    with pytest.raises(rd.GammaRegionError):
        rd.prepare_gamma_job(axes, np.zeros(40), np.zeros(40), SPEC_2_2, spacing_mm=1.0)  # nothing reaches the cutoff


class RecordingExecutor(SerialExecutor):
    """A serial executor that also keeps every submitted job's first argument."""

    def __init__(self) -> None:
        super().__init__()
        self.jobs: list[object] = []

    def submit(self, fn: Callable[..., object], /, *args: object, **kwargs: object) -> Future[object]:  # type: ignore[override]
        self.jobs.append(args[0])
        return super().submit(fn, *args, **kwargs)


def test_contrasts_take_the_upstream_mean_as_reference_and_noise_controls_the_two_halves_of_the_upstream_arm() -> None:
    shape = (30, 30, 30)
    axes = tuple(np.arange(n, dtype=float) for n in shape)
    base = blob(shape, (15, 15, 15), 3.0)
    scale = {"aup": 1.0, "apt": 1.05, "bup": 2.0, "bpg": 2.2, "bpi": 2.4}  # every arm distinguishable by its level
    means = {
        arm: rd.ArmMean(
            mean=f * base,
            low=0.98 * f * base if arm in ("aup", "bup") else None,
            high=1.02 * f * base if arm in ("aup", "bup") else None,
            seeds=list(range(1, 9)),
        )
        for arm, f in scale.items()
    }
    pool = RecordingExecutor()
    futures: dict[tuple[str, ...], Future[dict[str, object]]] = {}
    rd.submit_gamma_set(pool, futures, ("P", "100", "3d"), means, axes, [SPEC_2_2], spacing_mm=1.0)  # type: ignore[arg-type]
    assert len(futures) == 3 + 2  # three contrasts and one control per upstream arm
    jobs = dict(zip(futures, pool.jobs, strict=True))
    for label, (ev, ref) in rd.CONTRASTS.items():
        job = jobs[("P", "100", "3d", "contrast", label, "g2_2_c10")]
        assert isinstance(job, rd.GammaJob)
        assert job.global_norm == pytest.approx(scale[ref])  # normalised to the UPSTREAM maximum
        assert job.evaluation.max() / job.reference.max() == pytest.approx(scale[ev] / scale[ref])
    for arm, label in (("aup", "A-up"), ("bup", "B-up")):
        job = jobs[("P", "100", "3d", "noise_control", label, "g2_2_c10")]
        assert isinstance(job, rd.GammaJob)
        assert job.global_norm == pytest.approx(0.98 * scale[arm])  # the LOW half is the control's reference
        assert job.evaluation.max() / job.reference.max() == pytest.approx(1.02 / 0.98)
    assert not any(k[3] == "noise_control" and k[4] in ("A-port", "B-pgcc", "B-picc") for k in futures)


# ---------------------------------------------------------------------------------------------------------------------
# profiles and axes
# ---------------------------------------------------------------------------------------------------------------------


def test_lateral_profile_averages_the_slab_depth_bins_and_the_two_central_voxels() -> None:
    nk, na, nb = 12, 6, 8
    k, a, b = np.meshgrid(np.arange(nk), np.arange(na), np.arange(nb), indexing="ij")
    d = (k * 10000 + a * 100 + b).astype(np.float64)
    slab = np.arange(4, 10)
    along_a = rd.lateral_profile(d, slab, 0)
    along_b = rd.lateral_profile(d, slab, 1)
    assert along_a.shape == (na,)
    assert along_b.shape == (nb,)
    expected_a = [np.mean([kk * 10000 + aa * 100 + bb for kk in slab for bb in (nb // 2 - 1, nb // 2)]) for aa in range(na)]
    expected_b = [np.mean([kk * 10000 + aa * 100 + bb for kk in slab for aa in (na // 2 - 1, na // 2)]) for bb in range(nb)]
    np.testing.assert_allclose(along_a, expected_a)
    np.testing.assert_allclose(along_b, expected_b)
    assert rd.centre_value(along_a) == pytest.approx((along_a[na // 2 - 1] + along_a[na // 2]) / 2)


def test_the_axis_must_sit_on_a_voxel_corner_and_the_axis_index_is_checked() -> None:
    with pytest.raises(rd.DescriptiveError):
        rd.central_pair(5)
    with pytest.raises(rd.DescriptiveError):
        rd.lateral_profile(np.zeros((4, 4, 4)), np.arange(2), 2)


def test_axes_are_voxel_centres_in_mm() -> None:
    k, a, b = rd.p_axes((350, 400, 400))
    assert k[0] == 0.5
    assert k[-1] == 349.5
    assert a[0] == -199.5
    assert a[199] == -0.5
    assert a[200] == 0.5
    assert b[-1] == 199.5
    z, y, x = rd.f_axes((150, 150, 150))
    assert z[0] == -149.0
    assert x[-1] == 149.0
    assert z[74] == -1.0
    assert z[75] == 1.0
    assert y.shape == (150,)


def test_field_profile_is_100_percent_for_a_uniform_dose_and_follows_the_endpoint_definition() -> None:
    n = psm.N
    iy = psm.DEPTH_IY[rd.F_PROFILE_DEPTH_MM]
    np.testing.assert_allclose(rd.field_profile_percent(np.full((n, n, n), 3.0), iy), 100.0)
    d = np.ones((n, n, n))
    d[:, :, : n // 2 - 30] = 0.5  # lower dose far to one side only
    prof = rd.field_profile_percent(d, iy)
    assert prof.shape == (n // 2,)
    assert prof[0] == pytest.approx(100.0)  # normalised to the central samples
    assert prof[-1] < 100.0  # the 4-side mean picks the low side up


# ---------------------------------------------------------------------------------------------------------------------
# discovery, JSON and the command line
# ---------------------------------------------------------------------------------------------------------------------


TINY_DIMS = (8, 40, 8)  # DimSize x, y, z
TINY_SLABS = {100: (10, 15), 150: (12, 16), 200: (14, 18)}


def tiny_dose(seed: int, portable: bool) -> np.ndarray:
    """A small MCsquare-layout (z, y, x) dose: a Bragg-like depth profile times a lateral Gaussian."""
    rng = np.random.default_rng(seed)
    nx, ny, nz = TINY_DIMS
    k = np.arange(ny)[::-1]  # depth bin of MCsquare y index j is ny - 1 - j
    shift = rng.uniform(-0.05, 0.05) + (0.02 if portable else 0.0)
    kk = k - shift
    depth = np.where(kk <= 28, 0.3 + 0.7 * (np.clip(kk, 0, None) / 28) ** 3, np.clip(1 - (kk - 28) / 4, 0, None))
    lateral = (
        np.exp(-((np.arange(nz) - (nz - 1) / 2) ** 2) / (2 * 1.5**2))[:, None]
        * np.exp(-((np.arange(nx) - (nx - 1) / 2) ** 2) / (2 * 1.5**2))[None, :]
    )
    dose = lateral[:, None, :] * depth[None, :, None] * 1e6
    return (dose * (1 + 0.002 * rng.standard_normal(dose.shape))).astype(np.float32)


def write_tiny_run(run_dir: Path, arm: str, energy: int, seed: int, *, r80_error: float = 0.0) -> None:
    out = run_dir / "out"
    out.mkdir(parents=True)
    arr = tiny_dose(seed, portable=arm in ("apt", "bpg", "bpi"))
    nx, ny, nz = TINY_DIMS
    (out / "Dose.raw").write_bytes(arr.astype("<f4").tobytes())
    (out / "Dose.mhd").write_text(
        "ObjectType = Image\nNDims = 3\n"
        f"DimSize = {nx} {ny} {nz}\nElementSpacing = 1.000000 1.000000 1.000000\nOffset = 0.000000 0.000000 0.000000\n"
        "ElementType = MET_FLOAT\nElementByteOrderMSB = False\nElementDataFile = Dose.raw\n",
        encoding="utf-8",
    )
    sums = {n: hashlib.sha256((out / n).read_bytes()).hexdigest() for n in ("Dose.raw", "Dose.mhd")}
    (out / "sha256.txt").write_text("".join(f"{h}  {n}\n" for n, h in sums.items()), encoding="utf-8")
    idd = pe.canonical_from_mcsquare(arr).sum(axis=(1, 2))
    r80, multiple = pe.r80(idd)
    assert r80 is not None
    record = {"R80": r80 + r80_error, "R80_multiple_crossings": multiple}
    (run_dir / "endpoints.json").write_text("ENDPOINTS " + json.dumps(record) + "\n", encoding="utf-8")
    meta = {
        "arm": rd.ARM_LABEL[arm],
        "case": "P",
        "seed": seed,
        "commit": rd.ACQUISITION_COMMIT,
        "transport_status": "ok",
        "energy_mev": energy,
    }
    (run_dir / "run.json").write_text(json.dumps(meta), encoding="utf-8")


def tiny_tree(root: Path, *, bad_run: tuple[str, int, int] | None = None) -> dict[str, Path]:
    roots = {"A": root / "raw_A" / rd.ACQUISITION_DIR, "B": root / "raw_B" / rd.ACQUISITION_DIR}
    for ai, arm in enumerate(rd.ARM_ORDER):
        for ei, energy in enumerate(rd.ENERGIES):
            for i in range(1, 9):
                seed = 1000 * (ai + 1) + 100 * (ei + 1) + i
                run_dir = roots[rd.ARM_PART[arm]] / arm / "p" / f"e{energy}" / f"s{seed}"
                error = 1e-6 if bad_run == (arm, energy, i) else 0.0
                write_tiny_run(run_dir, arm, energy, seed, r80_error=error)
    return roots


def namespace(roots: dict[str, Path], out: Path) -> argparse.Namespace:
    return argparse.Namespace(raw_a=roots["A"], raw_b=roots["B"], out=out, workers=2, cases="P", analysis_json=None)


@pytest.fixture
def tiny_geometry(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(rd, "P_DIMS_XYZ", TINY_DIMS)
    monkeypatch.setattr(rd, "SLABS", TINY_SLABS)


@pytest.mark.usefixtures("tiny_geometry")
def test_discovery_finds_the_eight_seeds_of_every_cell_and_checks_run_json(tmp_path: Path) -> None:
    roots = tiny_tree(tmp_path)
    runs = rd.discover_runs(roots, ["P"])
    assert len(runs) == 5 * 3 * 8
    cells = rd.runs_by_cell(runs)
    assert len(cells) == 15
    assert all(len(v) == 8 for v in cells.values())
    assert all([r.seed for r in v] == sorted(r.seed for r in v) for v in cells.values())
    victim = next(iter(roots["A"].glob("aup/p/e100/s*"))) / "run.json"
    meta = json.loads(victim.read_text())
    victim.write_text(json.dumps({**meta, "arm": "A-port"}))  # a run.json that says it belongs to another arm
    with pytest.raises(rd.DiscoveryError):
        rd.discover_runs(roots, ["P"])


@pytest.mark.usefixtures("tiny_geometry")
def test_discovery_refuses_a_missing_seed_and_a_missing_arm(tmp_path: Path) -> None:
    roots = tiny_tree(tmp_path)
    next(iter(roots["B"].glob("bpg/p/e150/s*"))).rename(roots["B"] / "bpg" / "p" / "e150" / "not_a_seed")
    with pytest.raises(rd.DiscoveryError):
        rd.discover_runs(roots, ["P"])
    with pytest.raises(rd.DiscoveryError):
        rd.discover_runs({"A": tmp_path / "nowhere", "B": roots["B"]}, ["P"])


@pytest.mark.usefixtures("tiny_geometry")
def test_a_raw_root_may_or_may_not_include_the_acquisition_directory(tmp_path: Path) -> None:
    roots = tiny_tree(tmp_path)
    parents = {k: v.parent for k, v in roots.items()}
    assert len(rd.discover_runs(parents, ["P"])) == len(rd.discover_runs(roots, ["P"])) == 120


@pytest.mark.usefixtures("tiny_geometry")
def test_a_failed_binding_control_exits_nonzero_and_writes_nothing(tmp_path: Path) -> None:
    roots = tiny_tree(tmp_path, bad_run=("bpi", 150, 3))
    out = tmp_path / "out"
    with pytest.raises(rd.ControlError, match="bpi|B-picc"):
        rd.run(namespace(roots, out), executor_factory=SerialExecutor)
    assert not out.exists()


def test_main_maps_a_control_failure_to_exit_3_and_other_errors_to_2(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    argv = ["--raw-a", str(tmp_path), "--raw-b", str(tmp_path), "--out", str(tmp_path / "o")]

    def control_fails(*_a: object, **_k: object) -> int:
        msg = "R80 differs"
        raise rd.ControlError(msg)

    def discovery_fails(*_a: object, **_k: object) -> int:
        msg = "no arm"
        raise rd.DiscoveryError(msg)

    monkeypatch.setattr(rd, "run", control_fails)
    assert rd.main(argv) == 3
    monkeypatch.setattr(rd, "run", discovery_fails)
    assert rd.main(argv) == 2


def test_the_worker_count_is_capped_at_eight(tmp_path: Path) -> None:
    for workers in (0, 9, 32):
        ns = namespace({"A": tmp_path, "B": tmp_path}, tmp_path / "o")
        ns.workers = workers
        assert rd.run(ns, executor_factory=SerialExecutor) == 2
    assert not (tmp_path / "o").exists()


@pytest.mark.usefixtures("tiny_geometry")
def test_the_whole_pipeline_on_a_small_synthetic_tree_writes_json_and_figures(tmp_path: Path) -> None:
    roots = tiny_tree(tmp_path)
    out = tmp_path / "out"
    assert rd.run(namespace(roots, out), executor_factory=SerialExecutor) == 0
    doc = json.loads((out / rd.OUTPUT_JSON).read_text())
    # controls ran on every run and are recorded
    assert doc["controls"]["case_P_r80"]["runs_checked"] == 120
    assert doc["controls"]["case_P_r80"]["max_abs_diff_mm"] == 0.0
    assert doc["controls"]["dose_file_sha256"]["files_checked"] == 240
    assert all(v["identical_gamma_maps"] for v in doc["controls"]["crop_invariance"].values())
    # range metrics: every contrast, energy and metric present, Welch fields consistent
    cell = doc["range_metrics"]["contrasts"]["B2"]["200"]["R90"]
    assert cell["available"] is True
    assert cell["n_evaluated"] == cell["n_reference"] == 8
    assert cell["ci95"][0] < cell["estimate"] < cell["ci95"][1]
    assert cell["estimate"] == pytest.approx(cell["mean_evaluated"] - cell["mean_reference"])
    # gamma: contrasts and a noise control, with counts, for every energy and dimension
    for energy in ("100", "150", "200"):
        for dim in ("idd", "3d"):
            block = doc["gamma_pencil"][energy][dim]
            assert set(block["contrasts"]) == {"A", "B1", "B2"}
            assert set(block["noise_control"]) == {"A-up", "B-up"}
            for spec in ("g2_2_c10", "g1_1_c10"):
                for entry in (
                    *(block["contrasts"][c][spec] for c in block["contrasts"]),
                    *(block["noise_control"][a][spec] for a in block["noise_control"]),
                ):
                    assert entry["n_evaluated"] > 0
                    assert 0.0 <= entry["pass_rate_percent"] <= 100.0
                    assert entry["params"]["max_gamma"] == 2.0
                    assert entry["params"]["interp_fraction"] == 10
    sens = doc["controls"]["gamma_sensitivity"]
    assert sens["cells_checked"] == 6  # 3 energies x (idd, 3d)
    for energy in ("100", "150", "200"):
        for dim in ("idd", "3d"):
            assert sens["results"][energy][dim]["shift_depth_3mm"]["g2_2_c10"]["pass_rate_percent"] < 100.0
            assert set(sens["results"][energy][dim]) == {"shift_depth_3mm", "scale_1.03"}
    assert "sqrt(2)" in doc["design"]["noise_control"]
    assert doc["provenance"]["acquisition_commit"] == rd.ACQUISITION_COMMIT
    assert len(doc["provenance"]["dose_files"]) == 240
    for key in ("python", "numpy", "scipy", "pymedphys", "matplotlib", "host"):
        assert doc["provenance"][key]
    assert doc["figures"] == ["fig1_idd.png", "fig2_pencil_lateral.png", "fig2b_pencil_lateral_y.png"]
    for name in doc["figures"]:
        assert (out / name).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_json_writer_rejects_nan_and_serialises_numpy(tmp_path: Path) -> None:
    path = tmp_path / "d.json"
    rd.write_json_atomic(path, {"a": np.float32(1.5), "b": np.arange(3), "c": np.int64(4)})
    assert json.loads(path.read_text()) == {"a": 1.5, "b": [0, 1, 2], "c": 4}
    with pytest.raises(ValueError, match="Out of range|nan|NaN"):
        rd.write_json_atomic(tmp_path / "bad.json", {"a": float("nan")})
    assert not (tmp_path / "bad.json").exists()


def test_the_analysis_cross_check_reports_the_largest_difference() -> None:
    per_run = {
        (arm, energy, s): rd.RangeMetrics(
            90 + 0.1 * s,
            100 + 0.1 * s + (0.01 if arm in ("apt", "bpg", "bpi") else 0) + 0.001 * (s % 3),
            106 + 0.1 * s,
            1,
            1,
            1,
        )
        for arm in rd.ARM_ORDER
        for energy in rd.ENERGIES
        for s in range(1, 9)
    }
    contrasts = rd.range_contrasts(per_run)
    names = {"A": "A-port vs A-up", "B1": "B-pgcc vs B-up", "B2": "B-picc vs B-up"}
    analysis = {
        "contrasts": [
            {
                "name": name,
                "endpoints": [
                    {
                        "endpoint": f"P{e}/R80",
                        **{k: contrasts[label][str(e)]["R80"][k] for k in ("estimate", "se", "df", "ci95")},
                    }
                    for e in rd.ENERGIES
                ],
            }
            for label, name in names.items()
        ]
    }
    exact = rd.compare_with_analysis(contrasts, analysis)
    assert exact["r80_rows_compared"] == 9
    assert all(v == 0.0 for v in exact["max_abs_difference"].values())
    analysis["contrasts"][1]["endpoints"][2]["estimate"] += 0.25
    assert rd.compare_with_analysis(contrasts, analysis)["max_abs_difference"]["estimate"] == pytest.approx(0.25)


# ---------------------------------------------------------------------------------------------------------------------
# figures (smoke: they run on small synthetic data and write a PNG)
# ---------------------------------------------------------------------------------------------------------------------


def is_png(path: Path) -> bool:
    return path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n" and path.stat().st_size > 2000


def test_fig1_writes_a_png(tmp_path: Path) -> None:
    depth = np.arange(60) + 0.5
    idd = {e: {rd.ARM_LABEL[a]: tent(60, 40.0) * (1 + 0.002 * i) for i, a in enumerate(rd.ARM_ORDER)} for e in rd.ENERGIES}
    rd.plot_fig1_idd(tmp_path / "f.png", depth, idd)
    assert is_png(tmp_path / "f.png")


def test_fig2_writes_a_png_for_each_lateral_axis_and_tolerates_zero_dose(tmp_path: Path) -> None:
    x = np.arange(40) - 19.5
    base = np.exp(-(x**2) / 20) + 1e-5
    base[:3] = 0.0  # a log axis must not choke on exact zeros
    lateral = {
        e: {
            d: {ax: {rd.ARM_LABEL[a]: base * (1 + 0.01 * i) for i, a in enumerate(rd.ARM_ORDER)} for ax in (0, 1)}
            for d in TINY_SLABS[e]
        }
        for e in rd.ENERGIES
    }
    for ax in (0, 1):
        rd.plot_fig2_lateral(tmp_path / f"f{ax}.png", x, lateral, ax)
        assert is_png(tmp_path / f"f{ax}.png")


def test_fig3_writes_a_png(tmp_path: Path) -> None:
    dist = (np.arange(75) + 0.5) * 2.0
    prof = {
        rd.ARM_LABEL[a]: np.where(dist < 75, 100.0, 7.0 * np.exp(-(dist - 75) / 20)) * (1 + 0.003 * i)
        for i, a in enumerate(rd.ARM_ORDER)
    }
    rd.plot_fig3_field_edge(tmp_path / "f.png", dist, prof)
    assert is_png(tmp_path / "f.png")


def test_fig4_writes_a_png_with_blank_points_below_the_cutoff(tmp_path: Path) -> None:
    axis = (np.arange(20) - 9.5) * 2.0
    plane = np.full((20, 20), np.nan)
    plane[5:15, 5:15] = np.linspace(0, 2, 100).reshape(10, 10)
    planes = {f"{c}|{k}": plane for c in rd.CONTRASTS for k in ("g1_1_c10", "l2_2_c1")}
    rd.plot_fig4_field_gamma(
        tmp_path / "f.png", axis, axis, planes, [("g1_1_c10", "1%/1 mm global"), ("l2_2_c1", "2%/2 mm local")]
    )
    assert is_png(tmp_path / "f.png")
