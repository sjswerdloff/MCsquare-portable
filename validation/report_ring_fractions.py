"""Extract the mean annulus energy fractions of the upstream arms from the collected same-host endpoint records.

The analysis document holds ratios between arms, not the fractions themselves. The report's bound on the unresolved
far-halo endpoints needs to know how much energy each annulus holds, so this writes, for every case-P annulus
endpoint, the mean fraction over the upstream runs (A-up and B-up, 8 runs each per energy).

Usage:
    python validation/report_ring_fractions.py <collected dir> [--out validation/report_data/apples_ring_fractions_2f9dab40.json]
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from pathlib import Path

SCHEMA = "apples_ring_fractions/1"
UPSTREAM_ARMS = ("A-up", "B-up")
SEED_DIR = re.compile(r"^s\d+$")  # the collected tree also holds a .provenance directory per case
DEFAULT_OUT = Path(__file__).resolve().parent / "report_data" / "apples_ring_fractions_2f9dab40.json"


def ring_fractions(collected: Path) -> dict[str, list[float]]:
    """Every upstream run's annulus fractions, keyed by endpoint id (``P<energy>/ring_<depth>_<inner>_<outer>``)."""
    out: dict[str, list[float]] = {}
    for arm in UPSTREAM_ARMS:
        for run_dir in sorted(d for d in (collected / arm / "P").iterdir() if SEED_DIR.match(d.name)):
            energy = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))["energy_mev"]
            lines = [ln for ln in (run_dir / "endpoints.json").read_text(encoding="utf-8").splitlines() if ln.startswith("ENDPOINTS ")]
            if len(lines) != 1:
                msg = f"{run_dir}: {len(lines)} ENDPOINTS lines, expected 1"
                raise ValueError(msg)
            record = json.loads(lines[0][len("ENDPOINTS ") :])
            for key, value in record.items():
                if key.startswith("ring_"):
                    out.setdefault(f"P{energy}/{key}", []).append(float(value))
    return out


def document(collected: Path) -> dict[str, object]:
    """The extract: mean fraction and run count per endpoint."""
    per_run = ring_fractions(collected)
    counts = {len(v) for v in per_run.values()}
    if len(counts) != 1:
        msg = f"endpoints have differing run counts: {sorted(counts)}"
        raise ValueError(msg)
    return {
        "schema": SCHEMA,
        "arms": list(UPSTREAM_ARMS),
        "n_runs_per_endpoint": counts.pop(),
        "mean_fraction": {k: statistics.fmean(v) for k, v in sorted(per_run.items())},
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("collected", type=Path)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    try:
        doc = document(args.collected)
    except (OSError, ValueError, KeyError) as e:
        print(f"report_ring_fractions: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    args.out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
