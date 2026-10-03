"""Run times per 1e7 histories, by acquisition, arm and case, from the runs' own records (resource information only).

Reads every run.json under each --root (study label, then directory). For a study given --population (an analysis
document with a "runs" list of arm, energy and seed), only those runs are kept, and every one of them must be found.
Writes a document the report's Table 9 is rendered from (validation/report_confirmatory_tables.py).

    report_run_times.py --root apples=<collected> --root part_e=<tree> ... [--population part_e=<analysis.json>] --json <out>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import sys
from pathlib import Path

SCHEMA = "run_times/1"
PER = 1e7  # histories the table is normalised to


class RunTimesError(Exception):
    """The inputs cannot give a complete, unambiguous table."""


def collect(roots: list[tuple[str, Path]], populations: dict[str, Path]) -> dict[str, object]:
    """Every run's wall time, grouped by study, arm, host, threads, case and energy; refuses a missing population run."""
    runs: list[dict[str, object]] = []
    digest = hashlib.sha256()
    for study, root in roots:
        for path in sorted(root.rglob("run.json")):
            rec = json.loads(path.read_text(encoding="utf-8"))
            digest.update(hashlib.sha256(path.read_bytes()).hexdigest().encode() + b"\n")
            if not isinstance(rec.get("wall_s"), int | float) or not rec.get("simulated"):
                msg = f"{path}: no wall_s or simulated count"
                raise RunTimesError(msg)
            runs.append({"study": study, **{k: rec[k] for k in ("arm", "host", "cpu", "threads", "case", "energy_mev",
                                                                 "seed", "simulated", "wall_s", "commit")}})
    for study, doc_path in populations.items():
        doc = json.loads(doc_path.read_text(encoding="utf-8"))
        rerun = doc.get("rerun") or {}
        rerun_cells = set(rerun.get("cells", []))

        def commit_of(arm: str, energy: int, doc: dict = doc, rerun: dict = rerun, cells: set = rerun_cells) -> str:
            """The commit the analysis read this cell at: the rerun's for a rerun cell, else the acquisition's."""
            return rerun["commit"] if f"{arm} {energy} MeV" in cells else doc["acquisition_commit"]

        want = {(r["arm"], r["energy"], r["seed"]) for r in doc["runs"]}
        found: dict[tuple, list[dict[str, object]]] = {}
        for r in runs:
            key = (r["arm"], r["energy_mev"], r["seed"])
            if r["study"] == study and key in want and r["commit"] == commit_of(r["arm"], r["energy_mev"]):
                found.setdefault(key, []).append(r)
        bad = sorted(k for k in want if len(found.get(k, [])) != 1)
        if bad:
            msg = f"{study}: {len(bad)} population run(s) not found exactly once at their commit, e.g. {bad[:3]}"
            raise RunTimesError(msg)
        runs = [r for r in runs if r["study"] != study] + [found[k][0] for k in sorted(want)]
    groups: dict[tuple, list[float]] = {}
    for r in runs:
        key = (r["study"], r["arm"], r["host"], r["threads"], r["case"], r["energy_mev"])
        groups.setdefault(key, []).append(float(r["wall_s"]) * PER / float(r["simulated"]))
    return {
        "schema": SCHEMA,
        "per_histories": PER,
        "inputs_sha256": digest.hexdigest(),
        "n_runs": len(runs),
        "groups": [
            {"study": k[0], "arm": k[1], "host": k[2], "threads": k[3], "case": k[4], "energy_mev": k[5], "n": len(v),
             "median_s": statistics.median(v), "min_s": min(v), "max_s": max(v)}
            for k, v in sorted(groups.items())
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--root", action="append", required=True, help="STUDY=DIR; repeatable")
    ap.add_argument("--population", action="append", default=[], help="STUDY=analysis document; repeatable")
    ap.add_argument("--json", type=Path, required=True)
    a = ap.parse_args(argv)
    try:
        roots = [(s, Path(d)) for s, d in (x.split("=", 1) for x in a.root)]
        pops = {s: Path(d) for s, d in (x.split("=", 1) for x in a.population)}
        doc = collect(roots, pops)
    except (RunTimesError, OSError, KeyError, ValueError) as e:
        print(f"report_run_times: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    a.json.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
