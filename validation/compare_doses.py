#!/usr/bin/env python3
"""Compare two MCsquare dose distributions.

Developer validation tool for verifying portability changes
(e.g., MKL→PCG, Cilk→OpenMP, Intel→ARM64) produce equivalent results.

Outputs:
  - Gamma analysis pass rate (pymedphys)
  - Depth-dose curve (PDD) comparison along central axis
  - Lateral profile comparison at specified depths
  - Dose difference statistics

Usage:
  python validation/compare_doses.py reference.mhd evaluation.mhd
  python validation/compare_doses.py ref.mhd eval.mhd --criteria 2 2 --output report/
  python validation/compare_doses.py ref.mhd eval.mhd --profile-depths 5.0 10.0 15.0
"""

import argparse
import sys
from pathlib import Path

import numpy as np


def load_mhd(filepath: str) -> tuple[np.ndarray, dict]:
    """Load an MHD file and return (dose_array, metadata).

    Reads the MHD text header and associated raw binary data.
    Returns the dose as a 3D numpy array and metadata dict with
    spacing, offset, and dimension info.
    """
    filepath = Path(filepath)
    metadata = {}

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            metadata[key.strip()] = value.strip()

    dims = [int(x) for x in metadata["DimSize"].split()]
    spacing = [float(x) for x in metadata["ElementSpacing"].split()]
    offset = [float(x) for x in metadata["Offset"].split()]

    type_map = {
        "MET_FLOAT": np.float32,
        "MET_DOUBLE": np.float64,
        "MET_SHORT": np.int16,
        "MET_USHORT": np.uint16,
    }
    dtype = type_map.get(metadata["ElementType"], np.float32)

    raw_path = filepath.parent / metadata["ElementDataFile"]
    data = np.fromfile(raw_path, dtype=dtype)
    dose = data.reshape(dims[::-1])  # MHD uses x,y,z but numpy is z,y,x

    meta = {
        "dims": dims,
        "spacing": spacing,
        "offset": offset,
        "filepath": str(filepath),
    }
    return dose, meta


def build_axes(meta: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build coordinate axes (x, y, z) in mm from metadata."""
    axes = []
    for i in range(3):
        ax = (
            np.arange(meta["dims"][i]) * meta["spacing"][i]
            + meta["offset"][i]
        )
        axes.append(ax)
    return tuple(axes)


def extract_pdd(dose: np.ndarray, axes: tuple, axis: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Extract depth-dose curve along central axis.

    Args:
        dose: 3D dose array (z, y, x)
        axes: (x, y, z) coordinate arrays
        axis: depth axis index in axes tuple (default 2 = z)

    Returns:
        (depths, doses) along the central axis
    """
    # Central indices for the non-depth axes
    nz, ny, nx = dose.shape
    cy, cx = ny // 2, nx // 2
    pdd = dose[:, cy, cx]
    depths = axes[axis]
    return depths, pdd


def extract_profile(
    dose: np.ndarray,
    axes: tuple,
    depth_mm: float,
    direction: str = "x",
) -> tuple[np.ndarray, np.ndarray]:
    """Extract lateral profile at a given depth.

    Args:
        dose: 3D dose array (z, y, x)
        axes: (x, y, z) coordinate arrays
        depth_mm: depth in mm at which to extract the profile
        direction: 'x' for crossline or 'y' for inline

    Returns:
        (positions, doses) along the profile direction
    """
    z_axis = axes[2]
    z_idx = int(np.argmin(np.abs(z_axis - depth_mm)))

    nz, ny, nx = dose.shape
    if direction == "x":
        cy = ny // 2
        profile = dose[z_idx, cy, :]
        positions = axes[0]
    else:
        cx = nx // 2
        profile = dose[z_idx, :, cx]
        positions = axes[1]

    return positions, profile


def run_gamma(
    dose_ref: np.ndarray,
    dose_eval: np.ndarray,
    axes_ref: tuple,
    axes_eval: tuple,
    dose_threshold: float = 3.0,
    distance_threshold: float = 3.0,
    lower_cutoff: float = 20.0,
) -> tuple[np.ndarray, float]:
    """Run gamma analysis using pymedphys.

    Returns:
        (gamma_array, pass_rate_percent)
    """
    import pymedphys

    gamma = pymedphys.gamma(
        axes_ref,
        dose_ref,
        axes_eval,
        dose_eval,
        dose_threshold,
        distance_threshold,
        lower_percent_dose_cutoff=lower_cutoff,
        interp_fraction=10,
        max_gamma=2.0,
        quiet=True,
    )

    valid = np.isfinite(gamma)
    if valid.sum() == 0:
        return gamma, 0.0

    pass_rate = 100.0 * np.sum(gamma[valid] <= 1.0) / valid.sum()
    return gamma, pass_rate


def dose_difference_stats(
    dose_ref: np.ndarray, dose_eval: np.ndarray
) -> dict:
    """Compute voxel-wise dose difference statistics."""
    # Only compare where dose is above noise floor (1% of max)
    threshold = 0.01 * max(dose_ref.max(), dose_eval.max())
    mask = (dose_ref > threshold) | (dose_eval > threshold)

    if mask.sum() == 0:
        return {"max_abs_diff": 0.0, "mean_abs_diff": 0.0, "rms_diff": 0.0}

    diff = dose_eval[mask] - dose_ref[mask]
    ref_max = dose_ref.max()

    return {
        "max_abs_diff": float(np.max(np.abs(diff))),
        "max_abs_diff_pct": float(100.0 * np.max(np.abs(diff)) / ref_max) if ref_max > 0 else 0.0,
        "mean_abs_diff": float(np.mean(np.abs(diff))),
        "mean_abs_diff_pct": float(100.0 * np.mean(np.abs(diff)) / ref_max) if ref_max > 0 else 0.0,
        "rms_diff": float(np.sqrt(np.mean(diff**2))),
        "rms_diff_pct": float(100.0 * np.sqrt(np.mean(diff**2)) / ref_max) if ref_max > 0 else 0.0,
        "voxels_compared": int(mask.sum()),
    }


def save_plots(
    axes_ref, dose_ref, axes_eval, dose_eval,
    gamma, profile_depths, output_dir,
):
    """Save comparison plots to output directory."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # PDD comparison
    depths_ref, pdd_ref = extract_pdd(dose_ref, axes_ref)
    depths_eval, pdd_eval = extract_pdd(dose_eval, axes_eval)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(depths_ref, pdd_ref, "b-", label="Reference", linewidth=1.5)
    ax.plot(depths_eval, pdd_eval, "r--", label="Evaluation", linewidth=1.5)
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Dose")
    ax.set_title("Depth-Dose Curve (Central Axis)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.savefig(output_dir / "pdd_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Lateral profiles at specified depths
    for depth in profile_depths:
        pos_ref, prof_ref = extract_profile(dose_ref, axes_ref, depth, "x")
        pos_eval, prof_eval = extract_profile(dose_eval, axes_eval, depth, "x")

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(pos_ref, prof_ref, "b-", label="Reference", linewidth=1.5)
        ax.plot(pos_eval, prof_eval, "r--", label="Evaluation", linewidth=1.5)
        ax.set_xlabel("Position (mm)")
        ax.set_ylabel("Dose")
        ax.set_title(f"Lateral Profile at depth = {depth:.1f} mm")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.savefig(
            output_dir / f"profile_depth_{depth:.0f}mm.png",
            dpi=150, bbox_inches="tight",
        )
        plt.close(fig)

    # Gamma histogram
    if gamma is not None:
        valid_gamma = gamma[np.isfinite(gamma)]
        if len(valid_gamma) > 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.hist(valid_gamma, bins=100, range=(0, 2), edgecolor="black", linewidth=0.5)
            ax.axvline(x=1.0, color="r", linestyle="--", label="gamma = 1.0")
            ax.set_xlabel("Gamma value")
            ax.set_ylabel("Voxel count")
            ax.set_title("Gamma Distribution")
            ax.legend()
            fig.savefig(output_dir / "gamma_histogram.png", dpi=150, bbox_inches="tight")
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Compare two MCsquare dose distributions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("reference", help="Reference dose MHD file")
    parser.add_argument("evaluation", help="Evaluation dose MHD file")
    parser.add_argument(
        "--criteria", nargs=2, type=float, default=[3.0, 3.0],
        metavar=("DOSE_PCT", "DIST_MM"),
        help="Gamma criteria: dose%% and distance mm (default: 3 3)",
    )
    parser.add_argument(
        "--lower-cutoff", type=float, default=20.0,
        help="Lower dose cutoff for gamma as %% of max (default: 20)",
    )
    parser.add_argument(
        "--profile-depths", nargs="+", type=float,
        help="Depths (mm) for lateral profiles (default: auto at 25%%, 50%%, 75%% of max depth)",
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Directory for plot output (default: no plots, text summary only)",
    )
    parser.add_argument(
        "--no-gamma", action="store_true",
        help="Skip gamma analysis (faster, no pymedphys dependency)",
    )

    args = parser.parse_args()

    # Load doses
    print(f"Loading reference: {args.reference}")
    dose_ref, meta_ref = load_mhd(args.reference)
    axes_ref = build_axes(meta_ref)

    print(f"Loading evaluation: {args.evaluation}")
    dose_eval, meta_eval = load_mhd(args.evaluation)
    axes_eval = build_axes(meta_eval)

    print(f"Reference:  {meta_ref['dims']} voxels, spacing {meta_ref['spacing']} mm")
    print(f"Evaluation: {meta_eval['dims']} voxels, spacing {meta_eval['spacing']} mm")
    print()

    # Dose difference
    print("=== Dose Difference Statistics ===")
    stats = dose_difference_stats(dose_ref, dose_eval)
    print(f"  Max absolute difference:  {stats['max_abs_diff']:.6f} ({stats['max_abs_diff_pct']:.3f}%)")
    print(f"  Mean absolute difference: {stats['mean_abs_diff']:.6f} ({stats['mean_abs_diff_pct']:.3f}%)")
    print(f"  RMS difference:           {stats['rms_diff']:.6f} ({stats['rms_diff_pct']:.3f}%)")
    print(f"  Voxels compared:          {stats['voxels_compared']}")
    print()

    # PDD summary
    print("=== Depth-Dose (Central Axis) ===")
    depths_ref, pdd_ref = extract_pdd(dose_ref, axes_ref)
    depths_eval, pdd_eval = extract_pdd(dose_eval, axes_eval)
    if pdd_ref.max() > 0:
        pdd_diff_pct = 100.0 * np.max(np.abs(pdd_ref - pdd_eval)) / pdd_ref.max()
        print(f"  Max PDD difference: {pdd_diff_pct:.3f}% of max dose")
    print()

    # Gamma analysis
    gamma = None
    if not args.no_gamma:
        print(f"=== Gamma Analysis ({args.criteria[0]}%/{args.criteria[1]}mm) ===")
        try:
            gamma, pass_rate = run_gamma(
                dose_ref, dose_eval, axes_ref, axes_eval,
                args.criteria[0], args.criteria[1], args.lower_cutoff,
            )
            print(f"  Pass rate: {pass_rate:.1f}%")
            valid = np.isfinite(gamma)
            if valid.sum() > 0:
                print(f"  Mean gamma: {np.mean(gamma[valid]):.3f}")
                print(f"  Max gamma:  {np.max(gamma[valid]):.3f}")
                print(f"  Voxels evaluated: {valid.sum()}")
        except ImportError:
            print("  pymedphys not installed - skipping gamma analysis")
            print("  Install with: pip install pymedphys")
        print()

    # Determine profile depths
    if args.profile_depths:
        profile_depths = args.profile_depths
    else:
        z_max = axes_ref[2][-1]
        profile_depths = [z_max * f for f in (0.25, 0.50, 0.75)]

    # Profiles summary
    print("=== Lateral Profiles ===")
    for depth in profile_depths:
        _, prof_ref = extract_profile(dose_ref, axes_ref, depth, "x")
        _, prof_eval = extract_profile(dose_eval, axes_eval, depth, "x")
        if prof_ref.max() > 0:
            diff_pct = 100.0 * np.max(np.abs(prof_ref - prof_eval)) / prof_ref.max()
            print(f"  Depth {depth:.1f} mm: max difference {diff_pct:.3f}% of profile max")
    print()

    # Plots
    if args.output:
        print(f"Saving plots to {args.output}/")
        try:
            save_plots(
                axes_ref, dose_ref, axes_eval, dose_eval,
                gamma, profile_depths, args.output,
            )
            print("  Plots saved.")
        except ImportError:
            print("  matplotlib not installed - skipping plots")

    # Summary verdict
    if gamma is not None and np.isfinite(gamma).sum() > 0:
        valid = np.isfinite(gamma)
        pass_rate = 100.0 * np.sum(gamma[valid] <= 1.0) / valid.sum()
        if pass_rate >= 95.0:
            print(f"RESULT: PASS (gamma pass rate {pass_rate:.1f}% >= 95%)")
        else:
            print(f"RESULT: FAIL (gamma pass rate {pass_rate:.1f}% < 95%)")
            sys.exit(1)


if __name__ == "__main__":
    main()
