"""Descriptive dose analysis of the confirmatory same-host acquisition (commit 2f9dab40): figures, gamma, range metrics.

Everything here is DESCRIPTIVE: there are no margins and no pass/fail verdicts beyond the gamma index definition. The
confirmatory equivalence analysis is validation/apples_analyse.py; this script reads the same Dose files and reuses the
endpoint scripts the runs were made with (pencil_endpoints.py, platform_study_metrics.py, field_edge_analyse.py at
2f9dab40) for the dose reader, the canonical array and the case-F endpoints.

Direction of every difference: Portable minus upstream. The reference arm of a contrast is always the upstream arm.

    A  = A-port vs A-up    (Lenovo)        B1 = B-pgcc vs B-up    (HP)        B2 = B-picc vs B-up    (HP)

Stages
  1  Binding controls, before anything is written. Case P, 120 runs: R80 recomputed from the Dose file through the code
     path used for every other range metric must equal the recorded endpoints.json R80 to 1e-9 mm. Case F, 40 runs: all
     13 case-F endpoints recomputed by platform_study_metrics.endpoints must equal the recorded record.json metrics
     (relative difference <= 1e-12, None exact). Also: sha256 of every Dose.raw/Dose.mhd equals the run's sha256.txt.
     Any failure exits non-zero and writes nothing.
  2  Per-arm mean dose (voxel-wise mean of the 8 runs, runs NOT renormalised), range metrics, gamma, figure data.
  3  Figures and the JSON document.

Gamma is pymedphys.gamma (gamma_shell) with interp_fraction=10 and max_gamma=2 on the arm MEAN doses, axes in mm at voxel
centres. The noise control is the same gamma between two halves of the REFERENCE arm (its 4 lowest seeds against its 4
highest). Each half-mean has about sqrt(2) times the noise of an 8-run mean, so the control overstates the statistical
noise of the contrast it accompanies, by about that factor in the difference.

Usage:
    python report_dose_descriptives.py --raw-a <ROOT_A> --raw-b <ROOT_B> --out <DIR> [--workers N] [--cases P,F]
        [--analysis-json FILE]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import socket
import sys
import time
from collections.abc import Callable, Iterable, Mapping, Sequence
from concurrent.futures import Executor, Future, ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import field_edge_analyse
import numpy as np
import numpy.typing as npt
import pencil_endpoints as pe
import platform_study_metrics as psm
from scipy import stats

if TYPE_CHECKING:
    from types import ModuleType

    from matplotlib.axes import Axes

SCHEMA = "apples_descriptive/1"
ACQUISITION_COMMIT = "2f9dab404cea02f5072352f7ae5773dca27eb2a1"
ACQUISITION_DIR = ACQUISITION_COMMIT[:12]
OUTPUT_JSON = "apples_descriptive_2f9dab40.json"

ENERGIES = (100, 150, 200)
SEEDS_PER_CELL = 8
ARM_ORDER = ("aup", "apt", "bup", "bpg", "bpi")
ARM_LABEL = {"aup": "A-up", "apt": "A-port", "bup": "B-up", "bpg": "B-pgcc", "bpi": "B-picc"}
ARM_PART = {"aup": "A", "apt": "A", "bup": "B", "bpg": "B", "bpi": "B"}
UPSTREAM_OF_PART = {"A": "aup", "B": "bup"}
# contrast label -> (evaluated arm, reference arm); Portable minus upstream, reference = upstream
CONTRASTS = {"A": ("apt", "aup"), "B1": ("bpg", "bup"), "B2": ("bpi", "bup")}
SLABS = {100: (40, 60), 150: (80, 125), 200: (100, 200)}

# case P geometry as MCsquare writes it: DimSize (x, y, z), 1 mm voxels; canonical array D[k = depth, a = z, b = x]
P_DIMS_XYZ = (400, 350, 400)
P_SPACING_MM = 1.0
# case F: 300 mm cube, 2 mm voxels; array (z, y, x), beam along -y
F_DIMS_XYZ = (150, 150, 150)
F_SPACING_MM = 2.0
F_SIDE_MM = 150.0
F_FIELD_EDGE_MM = 75.0
F_PROFILE_DEPTH_MM = 127

LEVELS = {"R90": 0.9, "R80": 0.8, "R20": 0.2}
RANGE_METRICS = ("R90", "R80", "R20", "falloff")

R80_TOLERANCE_MM = 1e-9
FIELD_RTOL = 1e-12
FIELD_METRIC_KEYS = (
    "lateral_127_5",
    "lateral_127_10",
    "lateral_127_20",
    "lateral_127_30",
    "lateral_201_5",
    "lateral_201_10",
    "lateral_201_20",
    "lateral_201_30",
    "cax_127",
    "cax_201",
    "cax_i23",
    "r80_mm",
    "r20_mm",
)

INTERP_FRACTION = 10
SENSITIVITY_SHIFT_VOXELS = 3  # 3 mm at the 1 mm case-P voxel: beyond the 2 mm distance criterion
SENSITIVITY_SCALE = 1.03  # +3%: beyond the 2% dose criterion
MAX_GAMMA = 2.0
GAMMA_PAD_MM = 10.0
GAMMA_RAM_BYTES = 2 * 2**30
MAX_WORKERS = 8
CONFIDENCE = 0.95

SEED_DIR = re.compile(r"^s(\d+)$")

# categorical hues from the dataviz reference palette, ordered so that adjacent pairs separate under colour-vision
# deficiency; identity is also carried by line style (upstream solid, Portable broken), never by colour alone.
ARM_COLOUR = {"aup": "#2a78d6", "apt": "#eb6834", "bup": "#1baf7a", "bpg": "#4a3aa7", "bpi": "#e87ba4"}
ARM_STYLE = {"aup": "-", "apt": "--", "bup": "-", "bpg": "--", "bpi": ":"}
CONTRAST_ARM = {label: pair[0] for label, pair in CONTRASTS.items()}


class DescriptiveError(Exception):
    """Base class: the inputs or the analysis cannot be trusted."""


class DiscoveryError(DescriptiveError):
    """The raw run trees are not the 160 runs the acquisition describes."""


class ControlError(DescriptiveError):
    """A binding control failed: nothing may be written."""


class GammaRegionError(DescriptiveError):
    """A gamma region or result is not what the definition allows."""


# --------------------------------------------------------------------------------------------------------------------
# Gamma specifications
# --------------------------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class GammaSpec:
    """One gamma analysis: criteria, low-dose cutoff (percent of the reference maximum) and normalisation."""

    key: str
    dose_percent: float
    distance_mm: float
    cutoff_percent: float
    local: bool

    def label(self) -> str:
        """Human-readable criteria, e.g. ``2%/2 mm global, 10% cutoff``."""
        norm = "local" if self.local else "global"
        return f"{self.dose_percent:g}%/{self.distance_mm:g} mm {norm}, {self.cutoff_percent:g}% cutoff"


PENCIL_SPECS = (
    GammaSpec("g2_2_c10", 2.0, 2.0, 10.0, local=False),
    GammaSpec("g1_1_c10", 1.0, 1.0, 10.0, local=False),
)
FIELD_SPECS = (
    GammaSpec("g2_2_c10", 2.0, 2.0, 10.0, local=False),
    GammaSpec("g1_1_c10", 1.0, 1.0, 10.0, local=False),
    GammaSpec("g2_2_c1", 2.0, 2.0, 1.0, local=False),
    GammaSpec("l2_2_c1", 2.0, 2.0, 1.0, local=True),
)
FIELD_PLANE_SPECS = ("g1_1_c10", "l2_2_c1")  # the two analyses drawn as maps in Figure 4


# --------------------------------------------------------------------------------------------------------------------
# Range metrics (pure)
# --------------------------------------------------------------------------------------------------------------------


def integrated_depth_dose(canonical: np.ndarray) -> np.ndarray:
    """Dose summed over both lateral axes for each depth bin, exactly as pencil_endpoints.endpoints computes it."""
    return canonical.sum(axis=(1, 2))


def range_at_level(idd: np.ndarray, fraction: float) -> tuple[float | None, int]:
    """First distal crossing of ``fraction`` x maximum, linearly interpolated, and the number of distal crossings.

    This is pencil_endpoints.r80 with the level as a parameter: the same statements in the same order, so for
    ``fraction == 0.8`` the result is bit-identical to ``pencil_endpoints.r80``. Depth bin ``k`` has its centre at
    ``k + 0.5`` mm.

    Args:
        idd: Integrated depth dose, one value per 1 mm depth bin.
        fraction: Level as a fraction of the maximum (0.9 for R90, 0.8 for R80, 0.2 for R20).

    Returns:
        ``(depth_mm, n_crossings)``; ``(None, 0)`` when the level is never crossed distal to the maximum.
    """
    depths = np.arange(idd.size) + 0.5
    imax = int(np.argmax(idd))
    level = fraction * idd[imax]
    crossings = [i for i in range(imax, idd.size - 1) if idd[i] >= level > idd[i + 1]]
    if not crossings:
        return None, 0
    i = crossings[0]
    frac = (idd[i] - level) / (idd[i] - idd[i + 1])
    return float(depths[i] + frac * (depths[i + 1] - depths[i])), len(crossings)


@dataclass(frozen=True)
class RangeMetrics:
    """R90, R80 and R20 of one run's integrated depth dose, with the number of distal crossings of each level."""

    r90: float | None
    r80: float | None
    r20: float | None
    n90: int
    n80: int
    n20: int

    def value(self, name: str) -> float | None:
        """The metric if it is usable (one distal crossing), else None. The fall-off is R20 - R80 (positive)."""
        raw = {"R90": (self.r90, self.n90), "R80": (self.r80, self.n80), "R20": (self.r20, self.n20)}
        if name == "falloff":
            r20, r80 = self.value("R20"), self.value("R80")
            return None if r20 is None or r80 is None else r20 - r80
        depth, n = raw[name]
        return depth if depth is not None and n == 1 else None


def range_metrics(idd: np.ndarray) -> RangeMetrics:
    """R90, R80, R20 of one integrated depth dose by the first distal crossing."""
    (r90, n90), (r80, n80), (r20, n20) = (range_at_level(idd, LEVELS[k]) for k in ("R90", "R80", "R20"))
    return RangeMetrics(r90, r80, r20, n90, n80, n20)


def welch_difference(evaluated: Sequence[float], reference: Sequence[float], confidence: float = CONFIDENCE) -> dict[str, Any]:
    """Welch difference of means, evaluated minus reference, with its standard error, Satterthwaite df and interval.

    Args:
        evaluated: Per-run values of the evaluated (Portable) arm.
        reference: Per-run values of the reference (upstream) arm.
        confidence: Two-sided coverage of the interval.

    Returns:
        A dict with estimate, se, df, ci (lower, upper) and the two means and run counts. When both arms have zero
        variance the standard error is zero and the df undefined: ``df`` and ``ci`` are then None.

    Raises:
        DescriptiveError: If either arm has fewer than two runs or a non-finite value.
    """
    a, b = np.asarray(evaluated, dtype=np.float64), np.asarray(reference, dtype=np.float64)
    if a.size < 2 or b.size < 2 or not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        msg = f"Welch needs >= 2 finite values per arm, got {a.size} and {b.size}"
        raise DescriptiveError(msg)
    va, vb = a.var(ddof=1) / a.size, b.var(ddof=1) / b.size
    estimate = float(a.mean() - b.mean())
    out: dict[str, Any] = {
        "estimate": estimate,
        "se": float(math.sqrt(va + vb)),
        "mean_evaluated": float(a.mean()),
        "mean_reference": float(b.mean()),
        "n_evaluated": int(a.size),
        "n_reference": int(b.size),
        "df": None,
        "ci95": None,
    }
    if va + vb == 0.0:
        return out
    df = float((va + vb) ** 2 / (va**2 / (a.size - 1) + vb**2 / (b.size - 1)))
    half = float(stats.t.ppf(0.5 + confidence / 2, df)) * out["se"]
    out["df"], out["ci95"] = df, [estimate - half, estimate + half]
    return out


def contrast_cell(evaluated: Sequence[float | None], reference: Sequence[float | None]) -> dict[str, Any]:
    """One (contrast, energy, metric) cell: Welch statistics, or ``available: false`` if any run lacks the metric.

    A run lacks the metric when its level is absent or crossed more than once distal to the maximum. Nothing is
    imputed and no run is dropped.
    """
    bad_e = sum(v is None for v in evaluated)
    bad_r = sum(v is None for v in reference)
    if bad_e or bad_r:
        return {
            "available": False,
            "reason": "level absent or crossed more than once distal to the maximum in some runs",
            "n_unusable_evaluated": bad_e,
            "n_unusable_reference": bad_r,
        }
    clean_e = [float(v) for v in evaluated if v is not None]
    clean_r = [float(v) for v in reference if v is not None]
    return {"available": True, **welch_difference(clean_e, clean_r)}


# --------------------------------------------------------------------------------------------------------------------
# Binding controls (pure)
# --------------------------------------------------------------------------------------------------------------------


def check_r80_control(
    *,
    recorded: float | None,
    computed: float | None,
    recorded_multiple: bool,
    computed_crossings: int,
    tolerance: float = R80_TOLERANCE_MM,
) -> float:
    """Refuse a recomputed R80 that differs from the recorded one.

    Args:
        recorded: ``R80`` from the run's endpoints.json (None if the run recorded none).
        computed: R80 recomputed from the Dose file by the code path used for every other range metric.
        recorded_multiple: ``R80_multiple_crossings`` from endpoints.json.
        computed_crossings: Number of distal crossings counted by the recomputation.
        tolerance: Largest allowed absolute difference in mm.

    Returns:
        The absolute difference in mm (0.0 when both are None).

    Raises:
        ControlError: On a presence mismatch, a difference above ``tolerance``, or a multiple-crossing flag mismatch.
    """
    if (recorded is None) != (computed is None):
        msg = f"R80 presence differs: recorded {recorded!r}, recomputed {computed!r}"
        raise ControlError(msg)
    if bool(recorded_multiple) != (computed_crossings > 1):
        msg = (
            f"R80 multiple-crossing flag differs: recorded {recorded_multiple!r}, recomputed {computed_crossings} crossing(s)"
        )
        raise ControlError(msg)
    if recorded is None or computed is None:
        return 0.0
    diff = abs(float(computed) - float(recorded))
    if not diff <= tolerance:
        msg = f"R80 recomputed {computed!r} differs from recorded {recorded!r} by {diff:.3e} mm (tolerance {tolerance:g})"
        raise ControlError(msg)
    return diff


def check_field_control(recorded: Mapping[str, Any], computed: Mapping[str, Any], rtol: float = FIELD_RTOL) -> float:
    """Refuse recomputed case-F metrics that differ from the recorded ones.

    The 13 endpoint values must agree to a relative difference of ``rtol`` (None exact). Any other recorded key
    (crossing counts, validity flags) must be identical, and the two key sets must be equal.

    Returns:
        The largest relative difference over the 13 endpoint values.

    Raises:
        ControlError: On a key-set mismatch, a None/number mismatch, a relative difference above ``rtol``, or a
            non-endpoint key that differs.
    """
    if set(recorded) != set(computed):
        only_recorded, only_computed = sorted(set(recorded) - set(computed)), sorted(set(computed) - set(recorded))
        msg = f"metric keys differ: only recorded {only_recorded}, only recomputed {only_computed}"
        raise ControlError(msg)
    missing = [k for k in FIELD_METRIC_KEYS if k not in recorded]
    if missing:
        msg = f"recorded metrics lack {missing}"
        raise ControlError(msg)
    worst = 0.0
    for key in recorded:
        rec, new = recorded[key], computed[key]
        if key not in FIELD_METRIC_KEYS or rec is None or new is None or isinstance(rec, bool):
            if rec != new:
                msg = f"{key}: recorded {rec!r}, recomputed {new!r}"
                raise ControlError(msg)
            continue
        rel = abs(float(new) - float(rec)) / abs(float(rec)) if float(rec) != 0.0 else abs(float(new))
        if not rel <= rtol:
            msg = f"{key}: recorded {rec!r}, recomputed {new!r}, relative difference {rel:.3e} > {rtol:g}"
            raise ControlError(msg)
        worst = max(worst, rel)
    return worst


# --------------------------------------------------------------------------------------------------------------------
# Mean dose and the split-half partition (pure)
# --------------------------------------------------------------------------------------------------------------------


class RunningMean:
    """Voxel-wise mean of equally weighted arrays, accumulated in float64. Runs are NOT renormalised."""

    def __init__(self) -> None:
        self._sum: np.ndarray | None = None
        self.count = 0

    def add(self, array: np.ndarray) -> None:
        """Add one array; every array must have the same shape."""
        if self._sum is None:
            self._sum = np.array(array, dtype=np.float64)
        else:
            if array.shape != self._sum.shape:
                msg = f"shape {array.shape} differs from {self._sum.shape}"
                raise DescriptiveError(msg)
            np.add(self._sum, array, out=self._sum)
        self.count += 1

    def mean(self) -> np.ndarray:
        """The mean of the arrays added so far."""
        if self._sum is None:
            msg = "no arrays were added"
            raise DescriptiveError(msg)
        return self._sum / self.count


def mean_dose(arrays: Iterable[np.ndarray]) -> np.ndarray:
    """Voxel-wise mean over runs, float64, no per-run renormalisation."""
    acc = RunningMean()
    for a in arrays:
        acc.add(a)
    return acc.mean()


def split_half(seeds: Sequence[int]) -> tuple[list[int], list[int]]:
    """The lowest half and the highest half of a set of distinct seeds, each sorted ascending.

    Raises:
        DescriptiveError: If the seeds are not distinct or their count is odd or below two.
    """
    ordered = sorted(seeds)
    if len(set(ordered)) != len(ordered) or len(ordered) < 2 or len(ordered) % 2:
        msg = f"split-half needs an even number (>= 2) of distinct seeds, got {list(seeds)}"
        raise DescriptiveError(msg)
    half = len(ordered) // 2
    return ordered[:half], ordered[half:]


# --------------------------------------------------------------------------------------------------------------------
# Gamma (pymedphys.gamma), crop, summary
# --------------------------------------------------------------------------------------------------------------------


def depth_shifted(array: np.ndarray, voxels: int) -> np.ndarray:
    """The array moved ``voxels`` bins toward larger index along axis 0, the vacated bins zero (no wrap-around)."""
    if not 0 < voxels < array.shape[0]:
        msg = f"shift of {voxels} voxels does not fit an axis of {array.shape[0]}"
        raise DescriptiveError(msg)
    out = np.zeros_like(array)
    out[voxels:] = array[:-voxels]
    return out


def check_sensitivity(
    sensitivity: Mapping[str, Any], perturbation: str = "shift_depth_3mm", spec_key: str = "g2_2_c10"
) -> int:
    """Refuse a gamma instrument that cannot fail: a 3 mm depth shift must give a pass rate below 100% everywhere.

    Every gamma of the pencil beam came back at 100% on the real data, so the instrument has to be shown to return less
    than 100% on this class of input through the identical pipeline.

    Args:
        sensitivity: ``sensitivity[energy][dimension][perturbation][spec key]`` gamma summaries.
        perturbation: The deliberate perturbation that must be detected.
        spec_key: The analysis in which it must be detected.

    Returns:
        The number of (energy, dimension) cells checked.

    Raises:
        ControlError: If there is nothing to check, or any cell still reports a 100% pass rate.
    """
    checked = 0
    for energy, dims in sensitivity.items():
        for dim, perturbations in dims.items():
            rate = perturbations[perturbation][spec_key]["pass_rate_percent"]
            if rate is None or rate >= 100.0:
                msg = f"E{energy} {dim}: {perturbation} still passes at {rate!r}% ({spec_key}): the gamma cannot fail"
                raise ControlError(msg)
            checked += 1
    if not checked:
        msg = "no sensitivity cells to check"
        raise ControlError(msg)
    return checked


def cutoff_value(spec: GammaSpec, global_norm: float) -> float:
    """The absolute low-dose cutoff, computed as pymedphys computes it."""
    return spec.cutoff_percent / 100 * global_norm


def crop_slices(reference: np.ndarray, cutoff: float, pad_voxels: int) -> tuple[slice, ...]:
    """Bounding box of ``reference >= cutoff`` padded by ``pad_voxels`` on every side, clipped to the array.

    Raises:
        GammaRegionError: If no voxel reaches the cutoff.
    """
    region = reference >= cutoff
    if not region.any():
        msg = f"no reference voxel reaches the cutoff {cutoff:g}"
        raise GammaRegionError(msg)
    out = []
    for axis in range(reference.ndim):
        other = tuple(i for i in range(reference.ndim) if i != axis)
        idx = np.flatnonzero(region.any(axis=other)) if other else np.flatnonzero(region)
        out.append(slice(max(int(idx[0]) - pad_voxels, 0), min(int(idx[-1]) + 1 + pad_voxels, reference.shape[axis])))
    return tuple(out)


def compute_gamma(
    axes: Sequence[np.ndarray],
    reference: np.ndarray,
    evaluation: np.ndarray,
    spec: GammaSpec,
    global_norm: float,
    ram_bytes: int = GAMMA_RAM_BYTES,
) -> np.ndarray:
    """pymedphys.gamma on a shared grid with the module's fixed parameters.

    Points below the cutoff come back NaN, gamma above ``MAX_GAMMA`` is capped at it.
    """
    import pymedphys

    gamma = pymedphys.gamma(
        tuple(axes),
        reference,
        tuple(axes),
        evaluation,
        dose_percent_threshold=spec.dose_percent,
        distance_mm_threshold=spec.distance_mm,
        lower_percent_dose_cutoff=spec.cutoff_percent,
        interp_fraction=INTERP_FRACTION,
        max_gamma=MAX_GAMMA,
        local_gamma=spec.local,
        global_normalisation=global_norm,
        skip_once_passed=False,
        random_subset=None,
        ram_available=ram_bytes,
        interp_algo="pymedphys",
    )
    return np.asarray(gamma, dtype=np.float64)


def gamma_summary(gamma: np.ndarray, reference: np.ndarray, cutoff: float) -> dict[str, Any]:
    """Pass rate and counts of one gamma map.

    ``n_above_cutoff`` is the number of reference voxels the gamma was asked to evaluate (reference >= cutoff).
    ``n_evaluated`` is how many of them came back with a finite gamma; the rest are ``n_not_evaluated_in_region``.
    The pass rate is the percentage of evaluated points with gamma <= 1. The mean gamma is over evaluated points and
    is capped at ``MAX_GAMMA`` for the points whose true gamma exceeds it.

    Raises:
        GammaRegionError: If a finite gamma appears outside the region the cutoff defines.
    """
    region = reference >= cutoff
    finite = np.isfinite(gamma)
    if (finite & ~region).any():
        msg = "a finite gamma lies below the low-dose cutoff"
        raise GammaRegionError(msg)
    values = gamma[finite]
    n_eval = int(values.size)
    n_pass = int((values <= 1.0).sum())
    return {
        "n_above_cutoff": int(region.sum()),
        "n_evaluated": n_eval,
        "n_not_evaluated_in_region": int(region.sum()) - n_eval,
        "n_pass": n_pass,
        "n_fail": n_eval - n_pass,
        "pass_rate_percent": 100.0 * n_pass / n_eval if n_eval else None,
        "mean_gamma": float(values.mean()) if n_eval else None,
        "n_gamma_capped": int((values >= MAX_GAMMA).sum()),
    }


@dataclass
class GammaJob:
    """A gamma analysis with its (cropped) inputs, picklable for a worker process."""

    axes: tuple[np.ndarray, ...]
    reference: np.ndarray
    evaluation: np.ndarray
    spec: GammaSpec
    global_norm: float
    full_shape: tuple[int, ...]
    crop: tuple[slice, ...]
    plane: tuple[int, int] | None = None  # (axis, index) in uncropped coordinates: also return that plane of the map
    keep_map: bool = False  # also return the whole map in uncropped coordinates


def prepare_gamma_job(
    axes: Sequence[np.ndarray],
    reference: np.ndarray,
    evaluation: np.ndarray,
    spec: GammaSpec,
    *,
    spacing_mm: float,
    plane: tuple[int, int] | None = None,
    crop: bool = True,
    keep_map: bool = False,
) -> GammaJob:
    """Crop both volumes to the bounding box of the reference region padded by ``GAMMA_PAD_MM`` and bundle the job.

    The normalisation is the maximum of the UNCROPPED reference, and the pad is at least the largest search distance
    (``max_gamma`` x distance criterion), so the crop cannot change any gamma value. ``crop=False`` keeps the whole
    volume (used once per case to check that claim on real data).

    Raises:
        GammaRegionError: If the pad is shorter than the search distance, or the reference region is empty.
    """
    if MAX_GAMMA * spec.distance_mm > GAMMA_PAD_MM:
        msg = f"crop pad {GAMMA_PAD_MM} mm is shorter than the search distance {MAX_GAMMA * spec.distance_mm} mm"
        raise GammaRegionError(msg)
    if reference.shape != evaluation.shape or any(len(ax) != n for ax, n in zip(axes, reference.shape, strict=True)):
        msg = "reference, evaluation and axes do not share one grid"
        raise GammaRegionError(msg)
    global_norm = float(reference.max())
    if not (math.isfinite(global_norm) and global_norm > 0):
        msg = f"the reference maximum must be positive and finite, got {global_norm!r}"
        raise GammaRegionError(msg)
    box = (
        crop_slices(reference, cutoff_value(spec, global_norm), math.ceil(GAMMA_PAD_MM / spacing_mm))
        if crop
        else tuple(slice(0, n) for n in reference.shape)
    )
    return GammaJob(
        axes=tuple(np.ascontiguousarray(ax[s], dtype=np.float64) for ax, s in zip(axes, box, strict=True)),
        reference=np.ascontiguousarray(reference[box], dtype=np.float64),
        evaluation=np.ascontiguousarray(evaluation[box], dtype=np.float64),
        spec=spec,
        global_norm=global_norm,
        full_shape=tuple(reference.shape),
        crop=box,
        plane=plane,
        keep_map=keep_map,
    )


def run_gamma_job(job: GammaJob) -> dict[str, Any]:
    """Compute one gamma job (worker entry point): summary, exact parameters, crop, and the requested plane."""
    t0 = time.time()
    gamma = compute_gamma(job.axes, job.reference, job.evaluation, job.spec, job.global_norm)
    cutoff = cutoff_value(job.spec, job.global_norm)
    out: dict[str, Any] = {
        **gamma_summary(gamma, job.reference, cutoff),
        "criteria": job.spec.label(),
        "params": {
            "dose_percent_threshold": job.spec.dose_percent,
            "distance_mm_threshold": job.spec.distance_mm,
            "lower_percent_dose_cutoff": job.spec.cutoff_percent,
            "lower_dose_cutoff_absolute": cutoff,
            "local_gamma": job.spec.local,
            "global_normalisation": job.global_norm,
            "interp_fraction": INTERP_FRACTION,
            "max_gamma": MAX_GAMMA,
            "skip_once_passed": False,
            "random_subset": None,
            "interp_algo": "pymedphys",
            "ram_available_bytes": GAMMA_RAM_BYTES,
        },
        "crop": {
            "full_shape": list(job.full_shape),
            "slices": [[s.start, s.stop] for s in job.crop],
            "cropped_shape": list(job.reference.shape),
            "pad_mm": GAMMA_PAD_MM,
        },
        "elapsed_s": time.time() - t0,
    }
    if job.plane is not None or job.keep_map:
        full = np.full(job.full_shape, np.nan)
        full[job.crop] = gamma
        if job.plane is not None:
            out["plane"] = np.take(full, job.plane[1], axis=job.plane[0])
        if job.keep_map:
            out["map"] = full
    return out


# --------------------------------------------------------------------------------------------------------------------
# Profiles for the figures (pure)
# --------------------------------------------------------------------------------------------------------------------


def central_pair(n: int) -> slice:
    """The two voxels either side of the axis of an even-sized grid (indices n/2 - 1 and n/2)."""
    if n % 2:
        msg = f"lateral size {n} is odd: the axis is not on a voxel corner"
        raise DescriptiveError(msg)
    return slice(n // 2 - 1, n // 2 + 1)


def lateral_profile(mean: np.ndarray, slab: np.ndarray, axis: int) -> np.ndarray:
    """Dose along one lateral axis through the beam axis, averaged over the slab's depth bins.

    ``mean`` is the canonical D[k, a, b]. ``axis`` 0 gives the profile along the first lateral axis a (the two central b
    voxels are averaged), 1 along the second.
    """
    block = mean[slab]
    if axis == 0:
        return block[:, :, central_pair(block.shape[2])].mean(axis=(0, 2))
    if axis == 1:
        return block[:, central_pair(block.shape[1]), :].mean(axis=(0, 1))
    msg = f"lateral axis must be 0 or 1, got {axis}"
    raise DescriptiveError(msg)


def centre_value(profile: np.ndarray) -> float:
    """Mean of the two central samples of a profile on an even grid."""
    return float(profile[central_pair(profile.size)].mean())


# --------------------------------------------------------------------------------------------------------------------
# Run discovery and readers
# --------------------------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Run:
    """One run directory of the acquisition."""

    arm: str
    case: str
    energy: int
    seed: int
    run_dir: Path

    @property
    def out_dir(self) -> Path:
        """Directory holding Dose.mhd, Dose.raw and sha256.txt."""
        return self.run_dir / ("out" if self.case == "P" else "out_seed")

    @property
    def dose_mhd(self) -> Path:
        """The MHD header; the raw file sits beside it."""
        return self.out_dir / "Dose.mhd"


def arm_base(root: Path, arm: str) -> Path:
    """The arm's directory under a raw root, which may or may not include the acquisition directory."""
    for base in (root / arm, root / ACQUISITION_DIR / arm):
        if base.is_dir():
            return base
    msg = f"no directory for arm {arm} under {root}"
    raise DiscoveryError(msg)


def validate_identity(run: Run) -> None:
    """The run.json of a run must say it is this arm, case, seed and energy at the acquisition commit, transport ok."""
    meta = json.loads((run.run_dir / "run.json").read_text(encoding="utf-8"))
    expected = {
        "arm": ARM_LABEL[run.arm],
        "case": run.case,
        "seed": run.seed,
        "commit": ACQUISITION_COMMIT,
        "transport_status": "ok",
        "energy_mev": run.energy,
    }
    wrong = {k: (meta.get(k), v) for k, v in expected.items() if meta.get(k) != v}
    if wrong:
        msg = f"{run.run_dir}: run.json disagrees with its place in the tree (found, expected): {wrong}"
        raise DiscoveryError(msg)


def cell_runs(arm: str, case: str, energy: int, directory: Path) -> list[Run]:
    """The runs of one (arm, case, energy) cell: exactly ``SEEDS_PER_CELL`` seed directories, each identity-checked."""
    found = (
        sorted((int(m.group(1)), p) for p in directory.iterdir() if p.is_dir() and (m := SEED_DIR.match(p.name)))
        if directory.is_dir()
        else []
    )
    if len(found) != SEEDS_PER_CELL:
        msg = f"{directory}: {len(found)} seed directories, expected {SEEDS_PER_CELL}"
        raise DiscoveryError(msg)
    runs = [Run(arm, case, energy, seed, path) for seed, path in found]
    for run in runs:
        validate_identity(run)
    return runs


def discover_runs(roots: Mapping[str, Path], cases: Sequence[str]) -> list[Run]:
    """Every run of the requested cases, all five arms (40 per arm for case P at three energies, 8 for case F)."""
    runs: list[Run] = []
    for arm in ARM_ORDER:
        base = arm_base(roots[ARM_PART[arm]], arm)
        if "P" in cases:
            for energy in ENERGIES:
                runs += cell_runs(arm, "P", energy, base / "p" / f"e{energy}")
        if "F" in cases:
            runs += cell_runs(arm, "F", 200, base / "f")
    return runs


def load_p_canonical(mhd: Path, dims_xyz: Sequence[int] | None = None) -> np.ndarray:
    """A case-P Dose file on the canonical array D[k, a, b] (float32), geometry checked."""
    dims_xyz = P_DIMS_XYZ if dims_xyz is None else dims_xyz
    raw, dims, spacing = pe.read_mhd(mhd)
    if list(dims) != list(dims_xyz) or not np.allclose(spacing, [P_SPACING_MM] * 3):
        msg = f"{mhd}: dims {dims}, spacing {spacing}; expected {list(dims_xyz)} at {P_SPACING_MM} mm"
        raise DiscoveryError(msg)
    return pe.canonical_from_mcsquare(raw)


def load_f_array(mhd: Path) -> np.ndarray:
    """A case-F Dose file as (z, y, x) float32, geometry checked."""
    raw, dims, spacing = pe.read_mhd(mhd)
    if list(dims) != list(F_DIMS_XYZ) or not np.allclose(spacing, [F_SPACING_MM] * 3):
        msg = f"{mhd}: dims {dims}, spacing {spacing}; expected {list(F_DIMS_XYZ)} at {F_SPACING_MM} mm"
        raise DiscoveryError(msg)
    return raw


def sha256_file(path: Path) -> str:
    """Hex sha256 of a file, read in 16 MiB blocks."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def read_recorded_hashes(path: Path) -> dict[str, str]:
    """``<hash>  <name>`` lines of a run's sha256.txt."""
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[parts[1]] = parts[0]
    return out


def read_p_endpoints(path: Path) -> dict[str, Any]:
    """The ``ENDPOINTS {json}`` line of a case-P run."""
    text = path.read_text(encoding="utf-8").strip()
    if not text.startswith("ENDPOINTS "):
        msg = f"{path}: not an ENDPOINTS line"
        raise DiscoveryError(msg)
    return json.loads(text[len("ENDPOINTS ") :])


def stage1_run(run: Run) -> dict[str, Any]:
    """Worker for the binding controls: recompute from the Dose file and return both the new and the recorded values."""
    out_dir = run.out_dir
    raw_path = out_dir / "Dose.raw"
    hashes = read_recorded_hashes(out_dir / "sha256.txt")
    res: dict[str, Any] = {
        "arm": run.arm,
        "case": run.case,
        "energy": run.energy,
        "seed": run.seed,
        "files": [
            {"path": str(p), "bytes": p.stat().st_size, "sha256": sha256_file(p), "sha256_recorded": hashes.get(p.name)}
            for p in (raw_path, run.dose_mhd)
        ],
    }
    if run.case == "P":
        idd = integrated_depth_dose(load_p_canonical(run.dose_mhd))
        rm = range_metrics(idd)
        pe_r80, _ = pe.r80(idd)
        recorded = read_p_endpoints(run.run_dir / "endpoints.json")
        res |= {
            "idd": idd,
            "range": asdict(rm),
            "recorded_r80": recorded["R80"],
            "recorded_multiple": recorded["R80_multiple_crossings"],
            "pencil_endpoints_r80": pe_r80,
        }
    else:
        metrics = psm.endpoints(load_f_array(run.dose_mhd), F_SIDE_MM)
        recorded_metrics = json.loads((run.run_dir / "record.json").read_text(encoding="utf-8"))["metrics"]
        res |= {"metrics": metrics, "recorded_metrics": recorded_metrics}
    return res


def apply_controls(results: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Apply every binding control to the stage-1 results; raise ControlError listing all failures.

    Returns:
        The control summary for the JSON: runs checked and the largest differences seen.
    """
    failures: list[str] = []
    p_diffs: list[float] = []
    f_diffs: list[float] = []
    n_files = 0
    n_hash_bad = 0
    for r in results:
        tag = f"{ARM_LABEL[r['arm']]} {r['case']} E{r['energy']} s{r['seed']}"
        for f in r["files"]:
            n_files += 1
            if f["sha256"] != f["sha256_recorded"]:
                n_hash_bad += 1
                failures.append(f"{tag}: {Path(f['path']).name} sha256 {f['sha256']} != recorded {f['sha256_recorded']}")
        try:
            if r["case"] == "P":
                rm = r["range"]
                p_diffs.append(
                    check_r80_control(
                        recorded=r["recorded_r80"],
                        computed=rm["r80"],
                        recorded_multiple=r["recorded_multiple"],
                        computed_crossings=rm["n80"],
                    )
                )
                if r["pencil_endpoints_r80"] != rm["r80"]:
                    failures.append(
                        f"{tag}: range_at_level(0.8) {rm['r80']!r} != pencil_endpoints.r80 {r['pencil_endpoints_r80']!r}"
                    )
            else:
                f_diffs.append(check_field_control(r["recorded_metrics"], r["metrics"]))
        except ControlError as e:
            failures.append(f"{tag}: {e}")
    if failures:
        msg = f"{len(failures)} binding control failure(s); first: {failures[0]}"
        raise ControlError(msg)
    return {
        "case_P_r80": {
            "runs_checked": len(p_diffs),
            "tolerance_mm": R80_TOLERANCE_MM,
            "max_abs_diff_mm": max(p_diffs, default=None),
            "range_at_level_equals_pencil_endpoints_r80_runs": len(p_diffs),
            "passed": True,
        },
        "case_F_metrics": {
            "runs_checked": len(f_diffs),
            "metrics_compared_per_run": len(FIELD_METRIC_KEYS),
            "relative_tolerance": FIELD_RTOL,
            "max_relative_diff": max(f_diffs, default=None),
            "passed": True,
        },
        "dose_file_sha256": {"files_checked": n_files, "mismatches": n_hash_bad, "passed": True},
    }


# --------------------------------------------------------------------------------------------------------------------
# Analysis driver
# --------------------------------------------------------------------------------------------------------------------


def worker_init() -> None:
    """Worker start-up: one numba thread per worker process, so the worker count is the CPU budget."""
    os.environ["NUMBA_NUM_THREADS"] = "1"


def log(msg: str) -> None:
    """A timestamped progress line on stdout."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def runs_by_cell(runs: Iterable[Run]) -> dict[tuple[str, str, int], list[Run]]:
    """Runs grouped by (arm, case, energy), seed ascending."""
    out: dict[tuple[str, str, int], list[Run]] = {}
    for r in sorted(runs, key=lambda x: (x.arm, x.case, x.energy, x.seed)):
        out.setdefault((r.arm, r.case, r.energy), []).append(r)
    return out


@dataclass
class ArmMean:
    """The mean dose of one arm (and, for a reference arm, of its low and high seed halves)."""

    mean: np.ndarray
    low: np.ndarray | None
    high: np.ndarray | None
    seeds: list[int]


def accumulate_arm(runs: Sequence[Run], loader: Callable[[Run], np.ndarray], *, halves: bool) -> ArmMean:
    """Mean over the runs of one cell, streaming; with ``halves`` also the means of the 4 lowest and 4 highest seeds."""
    seeds = [r.seed for r in runs]
    low_seeds, _ = split_half(seeds) if halves else ([], [])
    total, low, high = RunningMean(), RunningMean(), RunningMean()
    for r in runs:
        d = loader(r)
        total.add(d)
        if halves:
            (low if r.seed in low_seeds else high).add(d)
    return ArmMean(total.mean(), low.mean() if halves else None, high.mean() if halves else None, sorted(seeds))


def p_axes(shape: Sequence[int]) -> tuple[np.ndarray, ...]:
    """Voxel-centre coordinates (mm) of a canonical case-P array: depth bin k at k + 0.5, lateral on the pe convention."""
    return (np.arange(shape[0]) + 0.5, pe.lateral_centres(shape[1]), pe.lateral_centres(shape[2]))


def f_axes(shape: Sequence[int]) -> tuple[np.ndarray, ...]:
    """Voxel-centre coordinates (mm) of a case-F array (z, y, x), 2 mm voxels, axis on the central corner."""
    return tuple((np.arange(n) + 0.5 - n / 2) * F_SPACING_MM for n in shape)


def range_contrasts(per_run: Mapping[tuple[str, int, int], RangeMetrics]) -> dict[str, Any]:
    """Welch statistics of R90, R80, R20 and the fall-off for every contrast and energy.

    Args:
        per_run: RangeMetrics by (arm, energy, seed).
    """
    out: dict[str, Any] = {}
    for label, (ev, ref) in CONTRASTS.items():
        out[label] = {}
        for energy in ENERGIES:
            cell: dict[str, Any] = {}
            for name in RANGE_METRICS:
                vals = {
                    arm: [
                        per_run[(arm, energy, s)].value(name)
                        for s in sorted(s for (a, e, s) in per_run if a == arm and e == energy)
                    ]
                    for arm in (ev, ref)
                }
                cell[name] = contrast_cell(vals[ev], vals[ref])
            out[label][str(energy)] = cell
    return out


def compare_with_analysis(contrasts: Mapping[str, Any], analysis: Mapping[str, Any]) -> dict[str, Any]:
    """Informational cross-check: our R80 Welch rows against the analysis document's ``P<E>/R80`` rows."""
    names = {"A": "A-port vs A-up", "B1": "B-pgcc vs B-up", "B2": "B-picc vs B-up"}
    by_name = {c["name"]: {e["endpoint"]: e for e in c["endpoints"]} for c in analysis["contrasts"]}
    worst = {"estimate": 0.0, "se": 0.0, "df": 0.0, "ci95": 0.0}
    n = 0
    for label, name in names.items():
        for energy in ENERGIES:
            mine = contrasts[label][str(energy)]["R80"]
            theirs = by_name[name][f"P{energy}/R80"]
            worst["estimate"] = max(worst["estimate"], abs(mine["estimate"] - theirs["estimate"]))
            worst["se"] = max(worst["se"], abs(mine["se"] - theirs["se"]))
            worst["df"] = max(worst["df"], abs(mine["df"] - theirs["df"]))
            worst["ci95"] = max(worst["ci95"], *(abs(a - b) for a, b in zip(mine["ci95"], theirs["ci95"], strict=True)))
            n += 1
    return {
        "r80_rows_compared": n,
        "max_abs_difference": worst,
        "against": "apples_analysis_2f9dab40.json (informational, not binding)",
    }


Futures = dict[tuple[str, ...], Future[dict[str, Any]]]


def p_loader(run: Run) -> np.ndarray:
    """Loader used for the case-P means: the canonical array."""
    return load_p_canonical(run.dose_mhd)


def f_loader(run: Run) -> np.ndarray:
    """Loader used for the case-F means: the (z, y, x) array."""
    return load_f_array(run.dose_mhd)


def field_profile_percent(mean: np.ndarray, iy: int) -> np.ndarray:
    """The case-F lateral profile of the endpoints (``side_profile``) in percent of the mean of its first 5 samples."""
    profile = field_edge_analyse.side_profile(mean, iy)
    return profile / profile[:5].mean() * 100


def submit_gamma_set(
    pool: Executor,
    futures: Futures,
    prefix: tuple[str, str, str],
    means: Mapping[str, ArmMean],
    axes: Sequence[np.ndarray],
    specs: Sequence[GammaSpec],
    *,
    spacing_mm: float,
    plane_row: int | None = None,
    reduce: Callable[[np.ndarray], np.ndarray] | None = None,
) -> None:
    """Submit the contrast gammas (Portable mean against upstream mean) and the split-half noise controls.

    Args:
        pool: Worker pool.
        futures: Filled in, keyed (case, energy, dimension, kind, contrast or reference arm, spec key).
        prefix: ``(case, energy, dimension)``.
        means: Arm means of one case and energy.
        axes: Voxel-centre axes of the arrays that are analysed (after ``reduce``).
        specs: The gamma analyses.
        spacing_mm: Voxel size, for the crop pad.
        plane_row: If given, also return the plane at this index of axis 1 for the specs in ``FIELD_PLANE_SPECS``.
        reduce: Applied to every mean before analysis (the lateral sum for the IDD).
    """
    prep = reduce or (lambda a: a)
    for label, (ev, ref) in CONTRASTS.items():
        for spec in specs:
            plane = (1, plane_row) if plane_row is not None and spec.key in FIELD_PLANE_SPECS else None
            job = prepare_gamma_job(
                axes, prep(means[ref].mean), prep(means[ev].mean), spec, spacing_mm=spacing_mm, plane=plane
            )
            futures[(*prefix, "contrast", label, spec.key)] = pool.submit(run_gamma_job, job)
    for arm in UPSTREAM_OF_PART.values():
        m = means[arm]
        if m.low is None or m.high is None:
            msg = f"{arm} has no split halves"
            raise DescriptiveError(msg)
        for spec in specs:
            job = prepare_gamma_job(axes, prep(m.low), prep(m.high), spec, spacing_mm=spacing_mm)
            futures[(*prefix, "noise_control", ARM_LABEL[arm], spec.key)] = pool.submit(run_gamma_job, job)


def submit_crop_check(
    pool: Executor,
    futures: Futures,
    prefix: tuple[str, str, str],
    means: Mapping[str, ArmMean],
    axes: Sequence[np.ndarray],
    spec: GammaSpec,
    *,
    spacing_mm: float,
) -> None:
    """Submit contrast A's analysis twice, cropped and on the whole volume, to compare the two maps on real data."""
    ref, ev = means["aup"].mean, means["apt"].mean
    for name, crop in (("cropped", True), ("full", False)):
        job = prepare_gamma_job(axes, ref, ev, spec, spacing_mm=spacing_mm, crop=crop, keep_map=True)
        futures[(*prefix, "crop_check", name, spec.key)] = pool.submit(run_gamma_job, job)


def submit_pencil_energy(
    energy: int,
    cells: Mapping[tuple[str, str, int], list[Run]],
    pool: Executor,
    futures: Futures,
    figure_data: dict[str, Any],
) -> None:
    """Mean doses of the five arms at one energy; store the figure data and submit every gamma job."""
    means: dict[str, ArmMean] = {}
    for arm in ARM_ORDER:
        t0 = time.time()
        means[arm] = accumulate_arm(cells[(arm, "P", energy)], p_loader, halves=arm in UPSTREAM_OF_PART.values())
        log(f"P E{energy} {ARM_LABEL[arm]}: mean of {len(means[arm].seeds)} runs in {time.time() - t0:.1f} s")
    axes = p_axes(means["aup"].mean.shape)
    figure_data["depth_mm"] = axes[0].tolist()
    figure_data["lateral_mm"] = axes[1].tolist()
    figure_data["idd"][str(energy)] = {ARM_LABEL[a]: integrated_depth_dose(m.mean).tolist() for a, m in means.items()}
    figure_data["lateral"][str(energy)] = {
        str(depth): {
            str(ax): {ARM_LABEL[a]: lateral_profile(m.mean, pe.slab_indices(depth), ax).tolist() for a, m in means.items()}
            for ax in (0, 1)
        }
        for depth in SLABS[energy]
    }
    submit_gamma_set(
        pool,
        futures,
        ("P", str(energy), "idd"),
        means,
        axes[:1],
        PENCIL_SPECS,
        spacing_mm=P_SPACING_MM,
        reduce=integrated_depth_dose,
    )
    submit_gamma_set(pool, futures, ("P", str(energy), "3d"), means, axes, PENCIL_SPECS, spacing_mm=P_SPACING_MM)
    if energy == ENERGIES[0]:
        submit_crop_check(pool, futures, ("P", str(energy), "3d"), means, axes, PENCIL_SPECS[0], spacing_mm=P_SPACING_MM)
    reference = means["aup"].mean
    perturbed = {
        f"shift_depth_{SENSITIVITY_SHIFT_VOXELS}mm": depth_shifted(reference, SENSITIVITY_SHIFT_VOXELS),
        f"scale_{SENSITIVITY_SCALE:g}": SENSITIVITY_SCALE * reference,
    }
    for name, evaluation in perturbed.items():
        for spec in PENCIL_SPECS:
            for dim, axs, ref_a, ev_a in (
                ("idd", axes[:1], integrated_depth_dose(reference), integrated_depth_dose(evaluation)),
                ("3d", axes, reference, evaluation),
            ):
                job = prepare_gamma_job(axs, ref_a, ev_a, spec, spacing_mm=P_SPACING_MM)
                futures[("P", str(energy), dim, "sensitivity", name, spec.key)] = pool.submit(run_gamma_job, job)


def submit_field(
    cells: Mapping[tuple[str, str, int], list[Run]], pool: Executor, futures: Futures, figure_data: dict[str, Any]
) -> None:
    """Mean doses of the five arms for case F; store the figure data and submit every gamma job."""
    means = {
        arm: accumulate_arm(cells[(arm, "F", 200)], f_loader, halves=arm in UPSTREAM_OF_PART.values()) for arm in ARM_ORDER
    }
    log("F: means of the five arms done")
    axes = f_axes(means["aup"].mean.shape)
    iy = psm.DEPTH_IY[F_PROFILE_DEPTH_MM]
    figure_data["field_profile"] = {
        "distance_from_axis_mm": ((np.arange(field_edge_analyse.C) + 0.5) * field_edge_analyse.SP).tolist(),
        "definition": "field_edge_analyse.side_profile at DEPTH_IY[127] (4-side mean, +-2 depth rows, +-10 voxels across), "
        "in percent of "
        "the mean of its first 5 samples: the profile behind the lateral_127_* endpoints of Table 4",
        "percent_of_central_axis": {ARM_LABEL[a]: field_profile_percent(m.mean, iy).tolist() for a, m in means.items()},
    }
    figure_data["field_plane_axes_mm"] = [axes[0].tolist(), axes[2].tolist()]
    submit_gamma_set(pool, futures, ("F", "200", "3d"), means, axes, FIELD_SPECS, spacing_mm=F_SPACING_MM, plane_row=iy)
    submit_crop_check(pool, futures, ("F", "200", "3d"), means, axes, FIELD_SPECS[0], spacing_mm=F_SPACING_MM)


def summarise_crop_check(cropped: Mapping[str, Any], full: Mapping[str, Any]) -> dict[str, Any]:
    """Compare a cropped and an uncropped gamma map of the same analysis: identical maps and identical counts."""
    same_map = bool(np.array_equal(cropped["map"], full["map"], equal_nan=True))
    keys = ("n_above_cutoff", "n_evaluated", "n_pass", "mean_gamma")
    same_counts = all(cropped[k] == full[k] for k in keys)
    both = np.isfinite(cropped["map"]) & np.isfinite(full["map"])
    return {
        "identical_gamma_maps": same_map,
        "identical_summary": same_counts,
        "max_abs_difference": float(np.abs(cropped["map"][both] - full["map"][both]).max()) if both.any() else None,
        "n_evaluated_cropped": cropped["n_evaluated"],
        "n_evaluated_full": full["n_evaluated"],
        "cropped_shape": cropped["crop"]["cropped_shape"],
        "full_shape": cropped["crop"]["full_shape"],
        "criteria": cropped["criteria"],
    }


@dataclass
class GammaResults:
    """Gamma results arranged for the JSON.

    ``pencil[energy][dim]["contrasts"|"noise_control"][label][spec]``; ``field["contrasts"|"noise_control"][label][spec]``;
    ``planes["<contrast>|<spec>"]`` the case-F maps; ``crop_checks["<case>_E<energy>_<dim>_<spec>"]``;
    ``sensitivity[energy][dim][perturbation][spec]`` the deliberately perturbed comparisons.
    """

    pencil: dict[str, Any]
    field: dict[str, Any]
    planes: dict[str, np.ndarray]
    crop_checks: dict[str, Any]
    sensitivity: dict[str, Any]


def nest_gamma_results(results: Mapping[tuple[str, ...], dict[str, Any]]) -> GammaResults:
    """Arrange the gamma results by case, kind and analysis."""
    pencil: dict[str, Any] = {}
    field: dict[str, Any] = {"contrasts": {}, "noise_control": {}}
    planes: dict[str, np.ndarray] = {}
    checks: dict[str, dict[str, Any]] = {}
    sensitivity: dict[str, Any] = {}
    for key, result in results.items():
        case, energy, dim, kind, who, spec = key
        res = {k: v for k, v in result.items() if k not in ("plane", "map")}
        if kind == "crop_check":
            checks.setdefault(f"{case}_E{energy}_{dim}_{spec}", {})[who] = result
        elif kind == "sensitivity":
            sensitivity.setdefault(energy, {}).setdefault(dim, {}).setdefault(who, {})[spec] = res
        else:
            slot = (
                pencil.setdefault(energy, {}).setdefault(dim, {"contrasts": {}, "noise_control": {}}) if case == "P" else field
            )
            slot["contrasts" if kind == "contrast" else "noise_control"].setdefault(who, {})[spec] = res
            if "plane" in result:
                planes[f"{who}|{spec}"] = result["plane"]
    crop_checks = {k: summarise_crop_check(v["cropped"], v["full"]) for k, v in checks.items()}
    return GammaResults(pencil, field, planes, crop_checks, sensitivity)


# --------------------------------------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------------------------------------


def _pyplot() -> ModuleType:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _style_axes(ax: Axes) -> None:
    ax.grid(visible=True, color="#d9d8d3", linewidth=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def plot_fig1_idd(path: Path, depth_mm: Sequence[float], idd: Mapping[int, Mapping[str, npt.ArrayLike]]) -> None:
    """Figure 1: mean IDD of the five arms (top) and Portable - upstream per contrast (bottom), one column per energy.

    Args:
        path: Output PNG.
        depth_mm: Depth-bin centres.
        idd: ``idd[energy][arm label]`` mean integrated depth dose.
    """
    plt = _pyplot()
    depth = np.asarray(depth_mm)
    fig, axes = plt.subplots(
        2, len(idd), figsize=(4.4 * len(idd), 7.2), sharex="col", gridspec_kw={"height_ratios": [3, 2]}, squeeze=False
    )
    for col, energy in enumerate(sorted(idd)):
        curves = {a: np.asarray(idd[energy][ARM_LABEL[a]]) for a in ARM_ORDER}
        top, bottom = axes[0, col], axes[1, col]
        norm = {part: curves[arm].max() for part, arm in UPSTREAM_OF_PART.items()}
        for arm in ARM_ORDER:
            top.plot(
                depth, curves[arm] / norm[ARM_PART[arm]], ARM_STYLE[arm], color=ARM_COLOUR[arm], lw=1.4, label=ARM_LABEL[arm]
            )
        for label, (ev, ref) in CONTRASTS.items():
            diff = (curves[ev] - curves[ref]) / norm[ARM_PART[ref]] * 100
            bottom.plot(
                depth, diff, ARM_STYLE[ev], color=ARM_COLOUR[ev], lw=1.2, label=f"{label}: {ARM_LABEL[ev]} - {ARM_LABEL[ref]}"
            )
        bottom.axhline(0, color="#52514e", lw=0.8)
        top.set_title(f"{energy} MeV")
        bottom.set_xlabel("depth (mm)")
        if col == 0:
            top.set_ylabel("IDD / maximum of the host's upstream arm")
            bottom.set_ylabel("Portable - upstream (% of upstream maximum)")
        for ax in (top, bottom):
            _style_axes(ax)
        top.legend(fontsize=8, frameon=False, loc="upper left")
        bottom.legend(fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_fig2_lateral(
    path: Path,
    lateral_mm: Sequence[float],
    lateral: Mapping[int, Mapping[int, Mapping[int, Mapping[str, Sequence[float]]]]],
    axis: int,
) -> None:
    """Figure 2 / 2b: lateral profiles through the beam axis at each energy's two depths, five arms, log axis.

    Args:
        path: Output PNG.
        lateral_mm: Lateral bin centres.
        lateral: ``lateral[energy][depth][axis][arm label]`` slab-averaged profile of the arm's mean dose.
        axis: 0 for the first lateral axis, 1 for the second.
    """
    plt = _pyplot()
    x = np.asarray(lateral_mm)
    energies = sorted(lateral)
    fig, axes = plt.subplots(2, len(energies), figsize=(4.6 * len(energies), 7.4), squeeze=False)
    for col, energy in enumerate(energies):
        for row, depth in enumerate(sorted(lateral[energy])):
            ax = axes[row, col]
            prof = {a: np.asarray(lateral[energy][depth][axis][ARM_LABEL[a]]) for a in ARM_ORDER}
            norm = {part: centre_value(prof[arm]) for part, arm in UPSTREAM_OF_PART.items()}
            for arm in ARM_ORDER:
                y = prof[arm] / norm[ARM_PART[arm]]
                ax.plot(x, np.where(y > 0, y, np.nan), ARM_STYLE[arm], color=ARM_COLOUR[arm], lw=1.1, label=ARM_LABEL[arm])
            ax.set_yscale("log")
            ax.set_title(f"{energy} MeV, depth {depth} mm")
            ax.set_xlabel(f"lateral position, axis {axis + 1} (mm)")
            if col == 0:
                ax.set_ylabel("dose / central value of the host's upstream arm")
            _style_axes(ax)
            if row == 0 and col == 0:
                ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_fig3_field_edge(
    path: Path, distance_mm: Sequence[float], profile: Mapping[str, Sequence[float]], edge_mm: float = F_FIELD_EDGE_MM
) -> None:
    """Figure 3: case-F lateral profile at 127 mm depth, five arms, log axis, with a linear inset 0-30 mm outside the edge.

    Args:
        path: Output PNG.
        distance_mm: Distance of each sample from the central axis.
        profile: ``profile[arm label]`` in percent of the central-axis dose.
        edge_mm: Nominal field edge (distance from the axis).
    """
    plt = _pyplot()
    x = np.asarray(distance_mm)
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    # the log-axis profile occupies the top and the lower right; the lower left is empty
    inset = ax.inset_axes((0.09, 0.12, 0.44, 0.42))
    for arm in ARM_ORDER:
        y = np.asarray(profile[ARM_LABEL[arm]])
        kw = {"color": ARM_COLOUR[arm], "ls": ARM_STYLE[arm], "lw": 1.3, "label": ARM_LABEL[arm]}
        ax.plot(x, np.where(y > 0, y, np.nan), **kw)
        sel = (x >= edge_mm) & (x <= edge_mm + 30)
        inset.plot(x[sel] - edge_mm, y[sel], **{**kw, "label": None})
    ax.axvline(edge_mm, color="#52514e", lw=0.9)
    ax.set_yscale("log")
    ax.set_xlabel("distance from the central axis (mm)")
    ax.set_ylabel("dose at 127 mm depth (% of central-axis dose)")
    ax.legend(fontsize=9, frameon=False, loc="upper right")
    inset.set_xlabel("outside the edge (mm)", fontsize=8)
    inset.set_ylabel("%", fontsize=8)
    inset.tick_params(labelsize=7)
    for a in (ax, inset):
        _style_axes(a)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_fig4_field_gamma(
    path: Path,
    x_mm: Sequence[float],
    z_mm: Sequence[float],
    planes: Mapping[str, np.ndarray],
    specs: Sequence[tuple[str, str]],
) -> None:
    """Figure 4: gamma maps on the 127 mm plane, one row per contrast, one column per analysis, scale 0-2.

    Args:
        path: Output PNG.
        x_mm: Coordinates of the plane's second (x) axis.
        z_mm: Coordinates of the plane's first (z) axis.
        planes: ``planes["<contrast>|<spec key>"]``, a (z, x) gamma map with NaN where not evaluated.
        specs: ``(spec key, panel label)`` for each column.
    """
    plt = _pyplot()
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("gamma", ["#2a78d6", "#cfcdc6", "#e34948"])
    cmap.set_bad("white")
    x, z = np.asarray(x_mm), np.asarray(z_mm)
    extent = (x[0] - (x[1] - x[0]) / 2, x[-1] + (x[1] - x[0]) / 2, z[-1] + (z[1] - z[0]) / 2, z[0] - (z[1] - z[0]) / 2)
    fig, axes = plt.subplots(
        len(CONTRASTS),
        len(specs),
        figsize=(4.8 * len(specs) + 1.2, 4.2 * len(CONTRASTS)),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    image = None
    for r, label in enumerate(CONTRASTS):
        for c, (key, panel) in enumerate(specs):
            ax = axes[r, c]
            image = ax.imshow(
                np.ma.masked_invalid(planes[f"{label}|{key}"]),
                cmap=cmap,
                vmin=0,
                vmax=2,
                extent=extent,
                origin="upper",
                interpolation="nearest",
            )
            ax.set_title(f"{label}: {panel}", fontsize=10)
            if r == len(CONTRASTS) - 1:
                ax.set_xlabel("x (mm)")
            if c == 0:
                ax.set_ylabel("z (mm)")
    fig.colorbar(image, ax=axes, label="gamma (capped at 2; blank: below the cutoff)", shrink=0.8)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------------------------------------
# Provenance and main
# --------------------------------------------------------------------------------------------------------------------


def provenance(dose_files: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Software, host and input provenance for the JSON."""
    import matplotlib
    import numba
    import pymedphys
    import scipy

    here = Path(__file__).resolve()
    scripts = {
        p.name: sha256_file(p) for p in (here, Path(pe.__file__), Path(psm.__file__), Path(field_edge_analyse.__file__))
    }
    return {
        "acquisition_commit": ACQUISITION_COMMIT,
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pymedphys": pymedphys.__version__,
        "matplotlib": matplotlib.__version__,
        "numba": numba.__version__,
        "numba_num_threads_in_workers": 1,
        "platform": platform.platform(),
        "host": socket.gethostname(),
        "script_sha256": scripts,
        "endpoint_scripts_note": "pencil_endpoints.py, platform_study_metrics.py and field_edge_analyse.py are the 2f9dab40 "
        "versions",
        "dose_files": list(dose_files),
    }


def json_default(o: object) -> float | int | list[Any]:
    """Serialise numpy scalars and arrays."""
    if isinstance(o, np.floating | np.integer):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    msg = f"not JSON serialisable: {type(o).__name__}"
    raise TypeError(msg)


def write_json_atomic(path: Path, doc: Mapping[str, Any]) -> None:
    """Write strict JSON (no NaN) via a temporary file, then rename."""
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(doc, indent=1, sort_keys=False, default=json_default, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def per_run_tables(
    results: Sequence[Mapping[str, Any]],
) -> tuple[dict[tuple[str, int, int], RangeMetrics], dict[str, Any], list[dict[str, Any]]]:
    """Range metrics by (arm, energy, seed), the same for the JSON, and the list of Dose files read."""
    per_run: dict[tuple[str, int, int], RangeMetrics] = {}
    per_run_json: dict[str, Any] = {}
    dose_files: list[dict[str, Any]] = []
    for r in results:
        dose_files += [
            {
                "arm": ARM_LABEL[r["arm"]],
                "case": r["case"],
                "energy": r["energy"],
                "seed": r["seed"],
                "path": f["path"],
                "bytes": f["bytes"],
                "sha256": f["sha256"],
            }
            for f in r["files"]
        ]
        if r["case"] != "P":
            continue
        rm = RangeMetrics(**r["range"])
        per_run[(r["arm"], r["energy"], r["seed"])] = rm
        per_run_json.setdefault(ARM_LABEL[r["arm"]], {}).setdefault(str(r["energy"]), {})[str(r["seed"])] = {
            **r["range"],
            **{name: rm.value(name) for name in RANGE_METRICS},
        }
    return per_run, per_run_json, dose_files


def design_block() -> dict[str, Any]:
    """The design and conventions block of the JSON."""
    return {
        "arms": {ARM_LABEL[a]: {"part": ARM_PART[a]} for a in ARM_ORDER},
        "contrasts": {k: {"evaluated": ARM_LABEL[e], "reference": ARM_LABEL[r]} for k, (e, r) in CONTRASTS.items()},
        "energies_mev": list(ENERGIES),
        "runs_per_cell": SEEDS_PER_CELL,
        "gamma_specs": {"pencil": {s.key: s.label() for s in PENCIL_SPECS}, "field": {s.key: s.label() for s in FIELD_SPECS}},
        "gamma_definition": (
            "pymedphys.gamma (gamma_shell) on the arm mean doses; reference = upstream mean, evaluation = Portable mean; "
            "axes in mm at "
            "voxel centres; interp_fraction=10, max_gamma=2; global normalisation = maximum of the reference; pass = gamma "
            "<= 1; pass "
            "rate over evaluated reference points (those >= the cutoff with a finite gamma); mean gamma is over those "
            "points with "
            "values above 2 capped at 2"
        ),
        "noise_control": (
            "The same gamma between two halves of the REFERENCE (upstream) arm: its 4 lowest seeds' mean (as reference, "
            "normalised to "
            "its own maximum, cutoff from its own maximum) against its 4 highest seeds' mean. Each half-mean has about "
            "sqrt(2) times the "
            "noise of an 8-run mean. One control per reference arm: A-up serves contrast A, B-up serves B1 and B2."
        ),
        "crop": (
            f"3D analyses are cropped to the bounding box of the reference >= its cutoff, padded by {GAMMA_PAD_MM:g} mm on "
            f"every side "
            f"(more than the {MAX_GAMMA * 2:g} mm largest search distance, so the crop cannot change a gamma value; checked "
            f"on real data "
            "in controls.crop_invariance). The normalisation is always the maximum of the uncropped reference."
        ),
        "mean_dose": "voxel-wise mean of the arm's 8 runs in float64; runs are not renormalised",
        "idd_gamma_input": "sum of the arm mean dose over both lateral axes",
        "figure_conventions": {
            "fig1": "IDD of each arm's mean dose, normalised to the maximum of its host's upstream arm; difference in % of "
            "that maximum",
            "fig2": "mean dose of the two central voxels across the other lateral axis, averaged over the six depth bins of "
            "slab_indices(depth); "
            "normalised to the mean of the two central samples of the host's upstream arm; signed lateral position; zero "
            "values are not drawn",
            "fig3": "see figure_data.field_profile.definition",
            "fig4": "gamma on the plane at DEPTH_IY[127] (z, x), colour scale 0-2, blank below the cutoff",
        },
    }


def make_figures(
    out_dir: Path, figure_data: Mapping[str, Any], planes: Mapping[str, np.ndarray], cases: Sequence[str]
) -> list[str]:
    """Write the figure PNGs for the requested cases and return their file names."""
    names: list[str] = []
    if "P" in cases:
        idd = {int(e): {lab: np.asarray(v) for lab, v in d.items()} for e, d in figure_data["idd"].items()}
        lateral = {
            int(e): {int(d): {int(ax): p for ax, p in axs.items()} for d, axs in ds.items()}
            for e, ds in figure_data["lateral"].items()
        }
        plot_fig1_idd(out_dir / "fig1_idd.png", figure_data["depth_mm"], idd)
        plot_fig2_lateral(out_dir / "fig2_pencil_lateral.png", figure_data["lateral_mm"], lateral, 0)
        plot_fig2_lateral(out_dir / "fig2b_pencil_lateral_y.png", figure_data["lateral_mm"], lateral, 1)
        names += ["fig1_idd.png", "fig2_pencil_lateral.png", "fig2b_pencil_lateral_y.png"]
    if "F" in cases:
        fp = figure_data["field_profile"]
        plot_fig3_field_edge(out_dir / "fig3_field_edge.png", fp["distance_from_axis_mm"], fp["percent_of_central_axis"])
        z_mm, x_mm = figure_data["field_plane_axes_mm"]
        plot_fig4_field_gamma(
            out_dir / "fig4_field_gamma.png",
            x_mm,
            z_mm,
            planes,
            [("g1_1_c10", "1%/1 mm global, 10% cutoff"), ("l2_2_c1", "2%/2 mm local, 1% cutoff")],
        )
        names += ["fig3_field_edge.png", "fig4_field_gamma.png"]
    return names


def run(args: argparse.Namespace, executor_factory: Callable[..., Executor] = ProcessPoolExecutor) -> int:
    """Run all stages; returns the exit code.

    Args:
        args: Parsed command line.
        executor_factory: Pool class for the workers (a thread pool in the tests, which cannot patch spawned processes).
    """
    if not 1 <= args.workers <= MAX_WORKERS:
        log(f"--workers must be 1..{MAX_WORKERS}")
        return 2
    cases = [c for c in args.cases.split(",") if c]
    out_dir = Path(args.out)
    t_start = time.time()
    runs = discover_runs({"A": Path(args.raw_a), "B": Path(args.raw_b)}, cases)
    log(f"discovered {len(runs)} runs ({sum(r.case == 'P' for r in runs)} case P, {sum(r.case == 'F' for r in runs)} case F)")

    futures: Futures = {}
    figure_data: dict[str, Any] = {"idd": {}, "lateral": {}}
    with executor_factory(max_workers=args.workers, initializer=worker_init) as pool:
        # Stage 1: the binding controls. Nothing is written, and no analysis starts, until every one has passed.
        results = list(pool.map(stage1_run, runs, chunksize=1))
        controls = apply_controls(results)
        log(f"binding controls passed: {json.dumps(controls)}")

        per_run, per_run_json, dose_files = per_run_tables(results)
        cells = runs_by_cell(runs)
        if "P" in cases:
            for energy in ENERGIES:
                submit_pencil_energy(energy, cells, pool, futures, figure_data)
        if "F" in cases:
            submit_field(cells, pool, futures, figure_data)
        log(f"{len(futures)} gamma jobs submitted; waiting")
        gamma_results: dict[tuple[str, ...], dict[str, Any]] = {}
        for key, fut in futures.items():
            res = fut.result()
            gamma_results[key] = res
            rate = res["pass_rate_percent"]
            log(
                f"gamma done {'/'.join(key)}: {'n/a' if rate is None else f'{rate:.3f}'}% of {res['n_evaluated']} "
                f"({res['elapsed_s']:.0f} s)"
            )

    nested = nest_gamma_results(gamma_results)
    controls["crop_invariance"] = nested.crop_checks
    bad = [k for k, v in nested.crop_checks.items() if not (v["identical_gamma_maps"] and v["identical_summary"])]
    if bad:
        msg = f"cropping changed the gamma map: {bad}"
        raise ControlError(msg)
    if nested.sensitivity:
        controls["gamma_sensitivity"] = {
            "definition": (
                f"The A-up mean dose against itself shifted {SENSITIVITY_SHIFT_VOXELS} mm in depth (zero-filled) and "
                f"against itself scaled by {SENSITIVITY_SCALE:g}, through the identical gamma path. The shift must give a "
                "pass rate below 100% in every energy and dimension (binding), so that the 100% pass rates of the real "
                "contrasts are not a dead instrument."
            ),
            "cells_checked": check_sensitivity(nested.sensitivity),
            "results": nested.sensitivity,
        }
    contrasts = range_contrasts(per_run) if per_run else {}
    doc: dict[str, Any] = {
        "schema": SCHEMA,
        "descriptive_only": True,
        "direction": "Portable minus upstream; the reference arm of every contrast is the upstream arm",
        "cases": cases,
        "controls": controls,
        "design": design_block(),
        "range_metrics": {
            "definition": (
                "R90/R80/R20: first distal crossing of 0.9/0.8/0.2 x the maximum of the whole-plane IDD with linear "
                "interpolation "
                "(pencil_endpoints.r80 with the level as a parameter); fall-off = R20 - R80 per run; Welch difference of "
                "the 8 run values "
                "per arm; a metric is unavailable in a cell if any of its 16 runs has no distal crossing or more than one; "
                "nothing is imputed"
            ),
            "contrasts": contrasts,
            "per_run": per_run_json,
        },
        "gamma_pencil": nested.pencil,
        "gamma_field": nested.field,
        "figure_data": figure_data,
        "provenance": provenance(dose_files),
        "wall_s": time.time() - t_start,
    }
    if args.analysis_json and contrasts:
        analysis = json.loads(Path(args.analysis_json).read_text(encoding="utf-8"))
        doc["cross_checks"] = {"r80_vs_analysis": compare_with_analysis(contrasts, analysis)}

    out_dir.mkdir(parents=True, exist_ok=True)
    doc["figures"] = make_figures(out_dir, figure_data, nested.planes, cases)
    write_json_atomic(out_dir / OUTPUT_JSON, doc)
    log(f"wrote {out_dir / OUTPUT_JSON} and {len(doc['figures'])} figures in {time.time() - t_start:.0f} s")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point. Exit 0 on success, 3 if a binding control failed (nothing written), 2 on other errors."""
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument(
        "--raw-a", required=True, help="raw run tree of part A (arm dirs aup, apt), with or without the 2f9dab404cea level"
    )
    ap.add_argument("--raw-b", required=True, help="raw run tree of part B (arm dirs bup, bpg, bpi)")
    ap.add_argument("--out", required=True, help="output directory (created only after the binding controls pass)")
    ap.add_argument("--workers", type=int, default=MAX_WORKERS, help=f"worker processes, 1..{MAX_WORKERS}")
    ap.add_argument("--cases", default="P,F", help="comma-separated cases to analyse (default P,F)")
    ap.add_argument("--analysis-json", default=None, help="analysis document for an informational R80 cross-check")
    args = ap.parse_args(argv)
    try:
        return run(args)
    except ControlError as e:
        print(f"report_dose_descriptives: BINDING CONTROL FAILED: {e}; nothing written", file=sys.stderr)
        return 3
    except (DescriptiveError, OSError, KeyError, ValueError) as e:
        print(f"report_dose_descriptives: {type(e).__name__}: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
