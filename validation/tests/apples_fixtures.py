"""Shared synthetic-data helpers for the apples tests: trees of the FROZEN seeds, with full-run metadata. No run output
from any machine is read."""

from __future__ import annotations

import hashlib
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

N_RUNS = 8
R80_BASE = {100: 77.5, 150: 158.0, 200: 259.0}
RING_FRAC = {(5, 10): 0.02, (10, 20): 0.05, (20, 40): 0.20, (40, 80): 0.40, (80, 200): 0.20}
ALL_ARMS = ("A-up", "A-port", "B-up", "B-pgcc", "B-picc")
ACQ = aa.ACQUISITION_COMMIT
SOURCES = {Path(p).name: (aa.REPO / p).read_text(encoding="utf-8") for p in aa.WORKFLOW_PATHS}
_REAL_EXPECTED = aa.expected_population


@pytest.fixture(autouse=True)
def _frozen_from_the_working_tree(monkeypatch: pytest.MonkeyPatch) -> None:
    """The committed workflow files are unchanged since the acquisition commit (asserted below), so tests read them
    from the working tree instead of needing the acquisition commit in the local object store."""

    def frozen(commit: str = ACQ, repo: Path = aa.REPO, sources: dict[str, str] | None = None) -> aa.Frozen:
        return _REAL_EXPECTED(commit, repo, SOURCES if sources is None else sources)

    monkeypatch.setattr(aa, "expected_population", frozen)


def _noise(rng: random.Random, sd: float) -> float:
    return rng.gauss(0.0, sd)


def p_record(rng: random.Random, energy: int, scale: float, arm: str = "A-up", seed: int = 0) -> dict:
    rec: dict = {
        "label": aa.P_LABEL.format(arm=arm, energy=energy, seed=seed),
        "layout": "mcsquare",
        "R80": R80_BASE[energy] + _noise(rng, 0.005 * scale),
        "R80_multiple_crossings": False,
        "slab_depths_mm": list(aa.SLAB_DEPTHS[energy]),
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


def binary_of(arm: str) -> str:
    return hashlib.sha256(f"binary of {arm}".encode()).hexdigest()


def run_json_of(arm: str, case: str, energy: int, seed: int) -> dict:
    """A complete full-run run.json as the workflows write it."""
    return {
        "seed": seed, "requested": aa.REQUESTED[case], "simulated": aa.REQUESTED[case] + 17,
        "overshoot": 17, "threads": aa.ARM_THREADS[arm], "energy_mev": energy, "case": case, "arm": arm,
        "mode": "full", "host": "HOST", "cpu": "cpu", "binary_sha256": binary_of(arm), "commit": ACQ,
        "start_utc": "2026-10-01T00:00:00Z", "end_utc": "2026-10-01T01:00:00Z", "wall_s": 3600,
        "transport_status": "ok",
    }


def f_document(arm: str, seed: int, metrics: dict) -> dict:
    """record.json as platform_study_record.py writes it at the acquisition commit (study defaults to pe1)."""
    cfg = hashlib.sha256(f"cfg {seed}".encode()).hexdigest()
    return {
        "study": "pe1", "platform": arm, "seed": seed, "host": "host", "machine": "x86_64", "commit": ACQ,
        "compiler": "cc", "sha256": {"config": cfg, "binary": binary_of(arm), "Dose.raw": "0" * 64},
        "materials": {"files": 1, "combined_sha256": "0" * 64}, "cfg_sha256": cfg, "metrics": metrics,
        "metrics_sha256": aa.metrics_digest(metrics), "endpoint_status": "ok",
    }


def cells(arm: str, case: str) -> list[tuple[int | None, int, int, int]]:
    """(energy, index within the energy cell, seed, ordinal in the arm/case) for every frozen run, in seed order."""
    seeds = aa.expected_population(ACQ).seeds[(arm, case)]
    seen: dict[int, int] = {}
    out = []
    for n, seed in enumerate(sorted(seeds)):
        energy = seeds[seed]
        out.append((energy if case == "P" else None, seen.get(energy, 0), seed, n))
        seen[energy] = seen.get(energy, 0) + 1
    return out


def build(
    root: Path,
    arms: tuple[str, ...],
    *,
    n: int = N_RUNS,
    scale: float = 1.0,
    identical: bool = False,
    mutate: Mutate | None = None,
    run_mutate: RunMutate | None = None,
    doc_mutate: Mutate | None = None,
) -> Path:
    """Write a synthetic tree of the FROZEN seeds (first `n` of each cell). `identical` gives every arm exactly the same
    values; mutate(arm, case, energy, i, rec) edits the endpoint record (case F: the metrics dict), run_mutate edits
    run.json and doc_mutate the whole endpoint document (case F: record.json), before they are written."""
    for a_idx, arm in enumerate(arms):
        for case in aa.CASES:
            for energy, i, seed, _ordinal in cells(arm, case):
                if i >= n:
                    continue
                e_idx = list(aa.SLAB_DEPTHS).index(energy) if energy is not None else 0
                rng = random.Random(7919 * (e_idx + 1) + 31 * i + (0 if identical else 104729 * (a_idx + 1)))
                rec = p_record(rng, energy, scale, arm, seed) if energy is not None else f_record(rng, scale)
                if mutate is not None:
                    mutate(arm, case, energy, i, rec)
                run = run_json_of(arm, case, energy if energy is not None else 200, seed)
                if run_mutate is not None:
                    run_mutate(arm, case, energy, i, run)
                doc = rec if energy is not None else f_document(arm, seed, rec)
                if doc_mutate is not None:
                    doc_mutate(arm, case, energy, i, doc)
                d = root / arm / case / f"s{seed}"
                d.mkdir(parents=True)
                (d / "run.json").write_text(json.dumps(run))
                if case == "P":
                    (d / "endpoints.json").write_text("ENDPOINTS " + json.dumps(doc, sort_keys=True) + "\n")
                else:
                    (d / "record.json").write_text(json.dumps(doc))
    return root


# ----------------------------------------------------------------------------------- acquisition-shaped native trees

NATIVE_ARM = {"A-up": "aup", "A-port": "apt", "B-up": "bup", "B-pgcc": "bpg", "B-picc": "bpi"}
WINDOWS_ARMS = ("A-up", "A-port")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hash_lines(files: dict[str, bytes], *, windows: bool) -> bytes:
    """A hash list as the workflow writes it: Windows `<hash>  rel\\path` with CRLF, Linux sha256sum `<hash>  ./rel`."""
    lines = []
    for rel in sorted(files):
        name = rel.replace("/", "\\") if windows else "./" + rel
        lines.append(f"{sha(files[rel])}  {name}")
    return (("\r\n" if windows else "\n").join(lines) + ("\r\n" if windows else "\n")).encode()


def write_tree(root: Path, files: dict[str, bytes]) -> None:
    for rel, data in files.items():
        f = root / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_bytes(data)


def make_native(parent: Path, arms: tuple[str, ...] = tuple(NATIVE_ARM), *, commit: str = ACQ) -> dict[str, Path]:
    """Native workflow trees of the frozen full runs, one root per host (`<parent>/win/<sha12>`, `<parent>/lin/<sha12>`).

    Returns {"win": root, "lin": root}. The layout, file names and hash-list formats follow apples-a-windows.yml and
    apples-b-linux.yml: <arm>/<p|f>/{run_root.txt, snapshot/, snapshot_sha256.txt, ...}, P under e<E>/s<seed>, F under
    s<seed>, run.json/endpoints.json/record.json from the same generators the analysis tests use."""
    roots = {"win": parent / "win" / commit[:12], "lin": parent / "lin" / commit[:12]}
    for host in roots.values():
        host.mkdir(parents=True)
    frozen = aa.expected_population(commit)
    for a_idx, arm in enumerate(arms):
        windows = arm in WINDOWS_ARMS
        host = roots["win" if windows else "lin"]
        binary = binary_of(arm)
        snapshot = {
            "Materials/Water.txt": b"water\n", "Scanners/Water_Phantom/HU_Density_Conversion.txt": b"hu density\n",
            "Scanners/Water_Phantom/HU_Material_Conversion.txt": b"hu material\n",
            "Scanners/default/HU_Density_Conversion.txt": b"default density\n",
            "Scanners/default/HU_Material_Conversion.txt": b"default material\n",
            "BDL/BDL_default_DN_RangeShifter.txt": b"bdl\n", "validation/pencil_endpoints.py": b"# script\n",
            "build.txt": b"host HOST\ncompiler cc\n", "binary_sha256.txt": f"{binary}  MCsquare.exe\n".encode(),
        }
        for case in aa.CASES:
            case_dir = host / NATIVE_ARM[arm] / case.lower()
            case_dir.mkdir(parents=True)
            (case_dir / "run_root.txt").write_bytes(f"{arm} case {case} mode full, commit {commit}, run 12345\n".encode())
            write_tree(case_dir / "snapshot", snapshot)
            (case_dir / "snapshot_sha256.txt").write_bytes(hash_lines(snapshot, windows=windows))
            # the runs themselves leave bytecode in the snapshot AFTER the hash list is written
            write_tree(case_dir / "snapshot", {"validation/__pycache__/x.pyc": b"pyc"})
            seeds = frozen.seeds[(arm, case)]
            if case == "F":
                fcase = {"cube.mhd": b"mhd", "cube.raw": b"raw", "E200_S150.txt": b"plan"}
                write_tree(case_dir / "fcase", fcase)
                (case_dir / "fcase_sha256.txt").write_bytes(hash_lines(fcase, windows=windows))
            inputs: dict[int, dict[str, bytes]] = {}
            if case == "P":
                for energy in sorted(set(seeds.values())):
                    inputs[energy] = {
                        "water.mhd": b"ct", f"bdl_E{energy}.txt": b"bdl", f"plan_E{energy}.txt": b"plan",
                        **{k: v for k, v in snapshot.items() if k.startswith(("Materials/", "Scanners/Water_Phantom/"))},
                    }
                    write_tree(case_dir / f"e{energy}", inputs[energy])
                    (case_dir / f"e{energy}_inputs_sha256.txt").write_bytes(hash_lines(inputs[energy], windows=windows))
            for energy, i, seed, _n in cells(arm, case):
                e_idx = list(aa.SLAB_DEPTHS).index(energy) if energy is not None else 0
                rng = random.Random(7919 * (e_idx + 1) + 31 * i + 104729 * (a_idx + 1))
                run = run_json_of(arm, case, energy if energy is not None else 200, seed)
                n_text = "1e7" if case == "P" else "3e7"
                dose = {"Dose.raw": f"raw {arm} {seed}".encode(), "Dose.mhd": f"mhd {arm} {seed}".encode()}
                if case == "P":
                    d = case_dir / f"e{energy}" / f"s{seed}"
                    rec = p_record(rng, energy, 1.0, arm, seed)
                    cfg_name, out_name, rec_name = "config.txt", "out", "endpoints.json"
                    rec_text = "ENDPOINTS " + json.dumps(rec, sort_keys=True) + "\n"
                else:
                    d = case_dir / f"s{seed}"
                    cfg_name, out_name, rec_name = "cfg.txt", "out_seed", "record.json"
                extra = {}
                cfg = f"Num_Threads {run['threads']}\nNum_Primaries {n_text}\nRNG_Seed {seed}\nOutput_Directory x\n".encode()
                if case == "F":
                    document = f_document(arm, seed, f_record(rng, 1.0))
                    document["cfg_sha256"] = document["sha256"]["config"] = sha(cfg)
                    document["sha256"].update(
                        {"Dose.raw": sha(dose["Dose.raw"]), "Dose.mhd": sha(dose["Dose.mhd"]), "binary": binary,
                         "cube.mhd": sha(fcase["cube.mhd"]), "cube.raw": sha(fcase["cube.raw"]),
                         "plan E200_S150.txt": sha(fcase["E200_S150.txt"])}
                    )
                    rec_text = json.dumps(document)
                    extra = {
                        **{k: v for k, v in fcase.items()},
                        **{k: v for k, v in snapshot.items() if k.startswith(("Materials/", "Scanners/default/", "BDL/"))},
                    }
                write_tree(d, {"run.json": json.dumps(run).encode(), rec_name: rec_text.encode(), cfg_name: cfg,
                               "log.txt": f"Nbr primaries simulated: {run['simulated']}\n".encode(), **extra})
                write_tree(d / out_name, {**dose, "sha256.txt": hash_lines(dose, windows=windows)})
    return roots
