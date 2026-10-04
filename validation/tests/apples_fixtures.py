"""Shared synthetic-data helpers for the apples tests: trees of the FROZEN seeds, with full-run metadata. No run output
from any machine is read.

Case-F record.json files are not hand-written: they are produced by RUNNING a record writer's main() in a synthetic run
directory (cwd = run directory, as both workflows run it), with synthetic transport inputs and the endpoint computation
stubbed to return (or raise) synthetic metrics. The default writer is the PINNED one, `git show 2f9dab40:validation/
platform_study_record.py` (the acquisition snapshot's), imported from a temporary copy and checked against its sha256;
the "v2" writer is the working tree's platform_study_record.py (cfg_sha256, metrics_sha256, endpoint_status, 0f5ef7c on).
"""

from __future__ import annotations

import contextlib
import functools
import hashlib
import importlib.util
import io
import json
import math
import os
import random
import subprocess
import sys
import tarfile
import tempfile
from collections.abc import Callable, Iterator
from pathlib import Path
from types import ModuleType

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


def run_json_of(arm: str, case: str, energy: int, seed: int, commit: str = ACQ) -> dict:
    """A complete full-run run.json as the workflows write it."""
    return {
        "seed": seed, "requested": aa.REQUESTED[case], "simulated": aa.REQUESTED[case] + 17,
        "overshoot": 17, "threads": aa.ARM_THREADS[arm], "energy_mev": energy, "case": case, "arm": arm,
        "mode": "full", "host": "HOST", "cpu": "cpu", "binary_sha256": binary_of(arm), "commit": commit,
        "start_utc": "2026-10-01T00:00:00Z", "end_utc": "2026-10-01T01:00:00Z", "wall_s": 3600,
        "transport_status": "ok",
    }


# ------------------------------------------------------------------------------------ the case-F record writers

WRITER_PATH = "validation/platform_study_record.py"
PINNED_WRITER_SHA256 = "c86444fa54e7f219488a4d95bed1b81230128b3b5f9bb07a92610cb24c916a0a"  # git show ACQ:WRITER_PATH
V2_COMMIT = "0f5ef7c76eb70e3e3089988f0fb17a964f86573f"  # the writer adds cfg_sha256/metrics_sha256/endpoint_status here
V2_BINDING = aa.RecordBinding("test-only-v2-at-0f5ef7c", V2_COMMIT, "pe1", aa.RECORD_SCHEMA_V2,
                              "tests only: no acquisition is bound to the v2 schema")


class WriterUnavailableError(RuntimeError):
    """The pinned writer cannot be read from the acquisition commit, or is not the expected bytes."""


@functools.cache
def pinned_writer() -> ModuleType:
    """platform_study_record.py as committed at the acquisition commit, imported from a temporary copy of its
    `git show` source (never from the working tree). Refuses unless the bytes hash as pinned."""
    proc = subprocess.run(["git", "-C", str(aa.REPO), "show", f"{ACQ}:{WRITER_PATH}"], capture_output=True, check=False)
    if proc.returncode != 0:
        msg = f"cannot read {WRITER_PATH} at {ACQ}: {proc.stderr.decode(errors='replace').strip()}"
        raise WriterUnavailableError(msg)
    if sha(proc.stdout) != PINNED_WRITER_SHA256:
        msg = f"{WRITER_PATH} at {ACQ} hashes {sha(proc.stdout)}, pinned {PINNED_WRITER_SHA256}"
        raise WriterUnavailableError(msg)
    copy = Path(tempfile.mkdtemp(prefix="writer_2f9dab40_")) / "platform_study_record_2f9dab40.py"
    copy.write_bytes(proc.stdout)
    spec = importlib.util.spec_from_file_location("platform_study_record_2f9dab40", copy)
    if spec is None or spec.loader is None:
        msg = f"cannot import {copy}"
        raise WriterUnavailableError(msg)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def current_writer() -> ModuleType:
    import platform_study_record

    return platform_study_record


WRITERS: dict[str, Callable[[], ModuleType]] = {"2f9dab40": pinned_writer, "v2": current_writer}


@contextlib.contextmanager
def _writer_context(writer: ModuleType, run_dir: Path, metrics: dict | BaseException) -> Iterator[None]:
    """cwd = the run directory, STUDY_ID unset (the workflows never set it), the endpoint computation stubbed."""

    def endpoints(_dose: object, _side: float) -> dict:
        if isinstance(metrics, BaseException):
            raise metrics
        return dict(metrics)

    saved = (writer.load, writer.endpoints, os.getcwd(), os.environ.pop("STUDY_ID", None))
    writer.load = lambda outdir: outdir
    writer.endpoints = endpoints
    os.chdir(run_dir)
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            yield
    finally:
        writer.load, writer.endpoints = saved[0], saved[1]
        os.chdir(saved[2])
        if saved[3] is not None:
            os.environ["STUDY_ID"] = saved[3]


def run_writer(
    run_dir: Path, arm: str, seed: int, binary: Path, metrics: dict | BaseException, *, commit: str = ACQ,
    writer: str = "2f9dab40",
) -> dict | None:
    """Run the writer's main() as the workflows do: `(cd <run dir> && platform_study_record.py <run dir> out_seed <arm>
    <seed> <binary> <compiler> <commit>)`. Returns record.json as written, or None when the writer wrote none (the
    2f9dab40 writer raises on an endpoint failure, before it opens record.json)."""
    module = WRITERS[writer]()
    with _writer_context(module, run_dir, metrics):
        try:
            module.main(str(run_dir), "out_seed", arm, str(seed), str(binary), "cc 1.0\n", commit)
        except Exception:  # the pinned writer's only failure shape: no record
            if writer != "2f9dab40":
                raise
    path = run_dir / "record.json"
    return json.loads(path.read_text()) if path.is_file() else None


def write_binary(path: Path, arm: str) -> Path:
    """A synthetic MCsquare whose sha256 is binary_of(arm)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_bytes(f"binary of {arm}".encode())
    return path


F_CASE = {"cube.mhd": b"mhd", "cube.raw": b"raw", "E200_S150.txt": b"plan"}
F_DATA = {"Materials/Water.txt": b"water\n", "Materials/Air/air.txt": b"air\n",
          "Scanners/default/HU_Density_Conversion.txt": b"default density\n",
          "Scanners/default/HU_Material_Conversion.txt": b"default material\n",
          "BDL/BDL_default_DN_RangeShifter.txt": b"bdl\n"}


def f_cfg(arm: str, seed: int) -> bytes:
    return f"Num_Threads {aa.ARM_THREADS[arm]}\nNum_Primaries 3e7\nRNG_Seed {seed}\nOutput_Directory out_seed\n".encode()


def synthetic_f_record(
    work: Path, arm: str, seed: int, metrics: dict | BaseException, *, commit: str = ACQ, writer: str = "2f9dab40"
) -> dict | None:
    """A case-F record from the writer, run in a synthetic run directory <work>/<arm>/s<seed> (transport inputs and
    Dose files are synthetic; the binary hashes as binary_of(arm))."""
    d = work / arm / f"s{seed}"
    write_tree(d, {**F_CASE, **F_DATA, "cfg.txt": f_cfg(arm, seed), "out_seed/Dose.raw": f"raw {arm} {seed}".encode(),
                   "out_seed/Dose.mhd": f"mhd {arm} {seed}".encode()})
    return run_writer(d, arm, seed, write_binary(work / arm / "MCsquare", arm), metrics, commit=commit, writer=writer)


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
    commit: str = ACQ,
    writer: str = "2f9dab40",
) -> Path:
    """Write a synthetic tree of the FROZEN seeds (first `n` of each cell). `identical` gives every arm exactly the same
    values; mutate(arm, case, energy, i, rec) edits the endpoint record (case F: the metrics dict, BEFORE the writer
    runs), run_mutate edits run.json and doc_mutate the whole endpoint document (case F: the writer's record.json),
    before they are written. Case-F records come from `writer` (synthetic_f_record), run under <root>/.writer/."""
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
                run = run_json_of(arm, case, energy if energy is not None else 200, seed, commit)
                if run_mutate is not None:
                    run_mutate(arm, case, energy, i, run)
                if energy is None:
                    written = synthetic_f_record(root / ".writer", arm, seed, rec, commit=commit, writer=writer)
                    if written is None:
                        msg = f"the {writer} writer wrote no record for {arm} s{seed}"
                        raise WriterUnavailableError(msg)
                    doc = written
                else:
                    doc = rec
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


def write_git_archive(path: Path, files: dict[str, bytes], commit: str) -> None:
    """A tar shaped like `git archive` output: a pax global header carrying the commit, directory entries, then files."""
    with tarfile.open(path, "w", format=tarfile.PAX_FORMAT, pax_headers={"comment": commit}) as tf:
        dirs = sorted({str(Path(rel).parent) for rel in files if "/" in rel})
        for d in dirs:
            info = tarfile.TarInfo(d + "/")
            info.type = tarfile.DIRTYPE
            tf.addfile(info)
        for rel, data in sorted(files.items()):
            info = tarfile.TarInfo(rel)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))


def make_native(
    parent: Path, arms: tuple[str, ...] = tuple(NATIVE_ARM), *, commit: str = ACQ, writer: str = "2f9dab40"
) -> dict[str, Path]:
    """Native workflow trees of the frozen full runs, one root per host (`<parent>/win/<sha12>`, `<parent>/lin/<sha12>`).

    Returns {"win": root, "lin": root}. The layout, file names and hash-list formats follow apples-a-windows.yml and
    apples-b-linux.yml: <arm>/<p|f>/{run_root.txt, snapshot/, snapshot_sha256.txt, ...}, P under e<E>/s<seed>, F under
    s<seed>, run.json/endpoints.json from the same generators the analysis tests use, and each record.json written by
    RUNNING `writer` in its run directory (see run_writer), with the binary at <parent>/bin/<arm>/MCsquare."""
    roots = {"win": parent / "win" / commit[:12], "lin": parent / "lin" / commit[:12]}
    for host in roots.values():
        host.mkdir(parents=True)
    frozen = aa.expected_population(commit)
    for a_idx, arm in enumerate(arms):
        windows = arm in WINDOWS_ARMS
        host = roots["win" if windows else "lin"]
        binary = binary_of(arm)
        binary_file = write_binary(parent / "bin" / arm / "MCsquare", arm)
        snapshot = {
            **F_DATA, "Scanners/Water_Phantom/HU_Density_Conversion.txt": b"hu density\n",
            "Scanners/Water_Phantom/HU_Material_Conversion.txt": b"hu material\n", "validation/pencil_endpoints.py": b"# script\n",
            "build.txt": b"host HOST\ncompiler cc\n", "binary_sha256.txt": f"{binary}  MCsquare.exe\n".encode(),
        }
        for case in aa.CASES:
            case_dir = host / NATIVE_ARM[arm] / case.lower()
            case_dir.mkdir(parents=True)
            (case_dir / "run_root.txt").write_bytes(f"{arm} case {case} mode full, commit {commit}, run 12345\n".encode())
            write_tree(case_dir / "snapshot", snapshot)
            (case_dir / "snapshot_sha256.txt").write_bytes(hash_lines(snapshot, windows=windows))
            if windows:  # apples-a-windows.yml extracts the snapshot from this git archive and leaves it in place
                write_git_archive(case_dir / "inputs.tar", {k: v for k, v in snapshot.items()
                                                            if k not in ("build.txt", "binary_sha256.txt")}, commit)
            # the runs themselves leave bytecode in the snapshot AFTER the hash list is written
            write_tree(case_dir / "snapshot", {"validation/__pycache__/x.pyc": b"pyc"})
            seeds = frozen.seeds[(arm, case)]
            if case == "F":
                fcase = dict(F_CASE)
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
                run = run_json_of(arm, case, energy if energy is not None else 200, seed, commit)
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
                    extra = {
                        **{k: v for k, v in fcase.items()},
                        **{k: v for k, v in snapshot.items() if k.startswith(("Materials/", "Scanners/default/", "BDL/"))},
                    }
                else:
                    extra = {rec_name: rec_text.encode()}
                write_tree(d, {"run.json": json.dumps(run).encode(), cfg_name: cfg,
                               "log.txt": f"Nbr primaries simulated: {run['simulated']}\n".encode(), **extra})
                write_tree(d / out_name, {**dose, "sha256.txt": hash_lines(dose, windows=windows)})
                if case == "F" and run_writer(d, arm, seed, binary_file, f_record(rng, 1.0), commit=commit, writer=writer) is None:
                    msg = f"the {writer} writer wrote no record for {arm} s{seed}"
                    raise WriterUnavailableError(msg)
    return roots


# ------------------------------------------------------------- failure shapes, as the two workflows leave them on disk


def native_case_dir(native: dict[str, Path], arm: str, case: str) -> Path:
    return native["win" if arm in WINDOWS_ARMS else "lin"] / NATIVE_ARM[arm] / case.lower()


def native_seed_dir(native: dict[str, Path], arm: str, case: str, seed: int) -> Path:
    case_dir = native_case_dir(native, arm, case)
    if case == "F":
        return case_dir / f"s{seed}"
    return case_dir / f"e{aa.expected_population(ACQ).seeds[(arm, case)][seed]}" / f"s{seed}"


def stop_job_after(native: dict[str, Path], arm: str, case: str, seed: int | None) -> list[int]:
    """The job exited at `seed` (None: before its first run): remove every LATER seed directory of its frozen loop, and
    (case P) every energy directory the loop never reached, with its inputs list. Returns the removed (absent) seeds."""
    import shutil

    seeds = aa.expected_population(ACQ).seeds[(arm, case)]
    order = list(seeds)
    later = order if seed is None else order[order.index(seed) + 1 :]
    for s in later:
        shutil.rmtree(native_seed_dir(native, arm, case, s))
    if case == "P":
        reached = {seeds[s] for s in order if s not in later}
        for energy in sorted(set(seeds.values()) - reached):
            case_dir = native_case_dir(native, arm, case)
            shutil.rmtree(case_dir / f"e{energy}")
            (case_dir / f"e{energy}_inputs_sha256.txt").unlink()
    return later


def edit_run_json(d: Path, **changes: object) -> None:
    run = json.loads((d / "run.json").read_text())
    run.update(changes)
    (d / "run.json").write_text(json.dumps(run))


def transport_failure(native: dict[str, Path], arm: str, case: str, seed: int, status: str, *, simulated: int | None = None,
                      keep_dose: bool = False) -> Path:
    """write_run <status>: run.json with that status (simulated null unless the log line had been read), no
    out/sha256.txt and no endpoint record (both written only after ok), Dose possibly absent."""
    d = native_seed_dir(native, arm, case, seed)
    edit_run_json(d, transport_status=status, simulated=simulated, overshoot=None)
    out = d / ("out" if case == "P" else "out_seed")
    (out / "sha256.txt").unlink()
    (d / ("endpoints.json" if case == "P" else "record.json")).unlink()
    if not keep_dose:
        (out / "Dose.raw").unlink()
        (out / "Dose.mhd").unlink()
    (d / "log.txt").write_text("MCsquare started\nSegmentation fault\n")
    return d


def interrupted(native: dict[str, Path], arm: str, case: str, seed: int) -> Path:
    """Cancelled mid-simulation: config written, log partial, out empty, no run.json and no endpoint record."""
    d = native_seed_dir(native, arm, case, seed)
    (d / "run.json").unlink()
    (d / ("endpoints.json" if case == "P" else "record.json")).unlink()
    (d / "log.txt").write_text("MCsquare started\n")
    out = d / ("out" if case == "P" else "out_seed")
    for f in out.iterdir():
        f.unlink()
    return d


# ------------------------------------------------------------------------- attesting a flat tree as a collection


def attest(root: Path, parts: tuple[str, ...] = ("A",), *, commit: str = ACQ, **override: object) -> dict:
    """Write a schema-3 collection_manifest.json over a flat tree as it now stands (every file under the parts' arm
    directories), standing in for apples_collect.py in analysis tests; `override` replaces top-level manifest keys."""
    arms = [arm for part in parts for arm in aa.PART_ARMS[part]]
    files = [
        {"path": f.relative_to(root).as_posix(), "source": "synthetic", "sha256": sha(f.read_bytes())}
        for arm in arms for f in sorted((root / arm).rglob("*")) if f.is_file()
    ]
    binding = aa.record_binding(commit)
    manifest: dict = {
        "schema": aa.MANIFEST_SCHEMA, "acquisition_commit": commit, "parts": list(parts), "files": files,
        "not_established": [],
        "f_record_binding": None if binding is None else {"name": binding.name, "schema": binding.schema.name},
    }
    manifest.update(override)
    (root / aa.MANIFEST_NAME).write_text(json.dumps(manifest))
    return manifest
