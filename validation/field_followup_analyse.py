"""Same-host comparison, part E: the analysis of docs/field_100_150_plan.md (field edge at 100 and 150 MeV).

Portable MCsquare against upstream OpenMCsquare with the fix, on the same host: A-port vs A-up (Lenovo), B-pgcc vs B-up
and B-picc vs B-up (HP), confirmatory; B-pgcc vs B-picc, descriptive. 17 endpoints per energy, 34 per contrast.

Inputs: the run trees the two workflows packed, extracted. Each --root is the directory that holds the arm directories
(aup, apt, bup, bpg, bpi): <root>/<arm dir>/e<energy>/s<seed>/. An arm missing from every root is a refusal, so that
one host's endpoints are never read while the other host's tree is absent.

No dose value is read until the WHOLE frozen population (every run of both hosts) has been verified. If any run is
missing, unexpected, duplicated or unverified, the analysis refuses, unless --collection-closed states who declares
the acquisition closed and on what evidence; that statement is written into the output, and the contrasts the issues
touch are PARTIAL. A directory being present is not evidence that its host has finished.

Order of work:
  1. the run lists and the case generator this script imports must be the files of the acquisition commit.
  2. every run of the frozen lists is VERIFIED without reading a dose value: run.json against the frozen seed, energy,
     threads, histories, host, commit and pinned binary; cfg.txt line for line; the log's simulated count; the plan
     and CT against what field_edge_make_cases.py writes; every Materials, scanner and beam-model file against the
     committed blob, with no file added or missing; Dose.raw and Dose.mhd against the sha256 written beside them on
     the run host; Dose.mhd against the grid this analysis assumes. An entry that looks like a run but is not a run
     of the frozen lists, anywhere under a root, is a population issue. --status stops here.
  3. only then, every endpoint is COMPUTED HERE from the verified Dose.raw (platform_study_metrics.field_endpoints). No endpoint
     value written on a run host is used.
  4. per endpoint: difference of run-level means (lateral, range) or geometric-mean ratio (central-axis doses), arm -
     reference; Welch SE and df; TOST at one-sided alpha 0.05; equivalent / not equivalent / inconclusive.

Rules fixed in the plan:
  - a run that is missing, unexpected, duplicated (two runs with the same Dose.raw) or fails any check makes its
    contrasts PARTIAL: descriptive estimates, no joint claim, no Holm decision.
  - an invalid value in any run (central-axis dose not finite and positive, R80 or R20 without exactly one crossing,
    lateral value not finite) makes that endpoint NOT ESTABLISHED. A lateral value of exactly zero is valid.
  - zero sample variance in both arms makes the endpoint NOT ESTABLISHED (no interval, no p-value; it stays in the
    Holm family as a non-rejection and counts against the joint claim).
  - joint claim: all 34 equivalent, unadjusted. Individual claims: Holm over the 34.
  - secondary joint claim: every endpoint other than those where every run of both arms is exactly zero is
    equivalent AFTER Holm over all 34. The rows set aside are chosen by the data, so the unadjusted tests give no
    error control for this claim; Holm over the fixed family does, whichever rows are set aside.
  - descriptive only: each arm's mean; at 50 and 70 mm the ratio of arm means with a 95% interval (delta method on
    the log of each mean); the ratio of the three-row mean around each run's maximum row, beside the maximum.

Usage:
  python field_followup_analyse.py --root DIR [--root DIR] --commit SHA --status [--mode smoke]
  python field_followup_analyse.py --root DIR [--root DIR] --commit SHA --json OUT.json --md OUT.md
                                   [--collection-closed "who, and the evidence"]

  Acquisition amendment 1 (the two cells cut by the 3 h timeout, run again at a later commit): add
  --rerun-root DIR [--rerun-root DIR] --rerun-commit SHA to either form. A-port and B-pgcc at 150 MeV are then
  read from the rerun roots only, verified against SHA, and the interrupted cells are not read.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
import os
import re
import statistics as st
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field, replace
from pathlib import Path

import apples_analyse as aa
import field_edge_analyse as fea
import field_edge_make_cases as fm
import field_followup_runs as fr
import numpy as np
import platform_study_metrics as psm
from scipy.stats import t as student_t

DEPTHS_MM = {100: (39, 61), 150: (79, 125)}
OFFSETS_MM = (5, 10, 20, 30, 50, 70)
RATIO_OFFSETS_MM = (50, 70)  # the rows that also get a descriptive ratio of arm means
FAMILY_SIZE = 34
ALPHA = aa.ALPHA
LATERAL_MARGIN_5MM, LATERAL_MARGIN = 0.50, 0.20  # points of the central-axis dose
CAX_RATIO = (0.995, 1.005)
RANGE_MARGIN_MM = 0.5
GRID_N = 150  # voxels per side, 2 mm: the cube field_edge_make_cases.py writes and Dose.raw must have
REPO = Path(__file__).resolve().parents[1]
INPUT_PATHS = ("Materials", "Scanners/default", "BDL/BDL_default_DN_RangeShifter.txt")  # copied into every run directory
RUN_LIST_FILES = ("validation/field_followup_runs.py", "validation/field_edge_make_cases.py")  # must be the acquisition's
ANALYSIS_FILES = ("validation/field_followup_analyse.py", "validation/platform_study_metrics.py",
                  "validation/field_edge_analyse.py", "validation/apples_analyse.py")  # identity recorded, not required
LABEL = ("Same host, same thread count, and the same fix in both codes. Margins are those of the 200 MeV field (#31). "
         "One field (15 x 15 cm) in a uniform medium, no range shifter. Every endpoint was computed by the analysis from "
         "hash-verified dose files.")
_SEED_DIR = re.compile(r"^s(\d{6})$")
_RUN_NAME = re.compile(r"^s\d")  # named like a run directory, whatever follows
_RUN_FILES = ("run.json", "Dose.raw")  # a directory holding either, at any depth, holds run data
CELL_ENTRIES = ("snapshot", "snapshot_sha256.txt", "fcase", "fcase_sha256.txt", "run_root.txt", "inputs.tar")  # beside the runs
# Acquisition amendment 1 (docs/field_100_150_plan.md): these two cells were stopped at the runners' 3 h job timeout and
# were run again in full at a later commit. With --rerun-root/--rerun-commit they are read from the rerun and only from it.
RERUN_CELLS = (("A-port", 150), ("B-pgcc", 150))
# All a rerun commit may change relative to the acquisition commit: the two run requests and the plan.
RERUN_MAY_CHANGE = (".gitea/workflows/field-e-linux.yml", ".gitea/workflows/field-e-windows.yml", "docs/field_100_150_plan.md")
_COUNT_LABEL = "Nbr primaries simulated"
_COUNT_LINE = re.compile(r"Nbr primaries simulated: (\d+)(?: \(\d+ generated outside the geometry\))?")  # compute_simulation.c
_SHA_LINE = re.compile(r"^([0-9a-f]{64})  (\S+)$")
InputError = aa.InputError
Spec = aa.Spec


# ------------------------------------------------------------------------------------------------- specification
def margin(key: str) -> tuple[float, float, bool]:
    """(lower, upper, log scale) for an endpoint key: the #31 margins."""
    if key.startswith("lateral_"):
        m = LATERAL_MARGIN_5MM if int(key.rsplit("_", 1)[1]) == 5 else LATERAL_MARGIN
        return -m, m, False
    if key.startswith("cax_"):
        return math.log(CAX_RATIO[0]), math.log(CAX_RATIO[1]), True
    if key in ("r80_mm", "r20_mm"):
        return -RANGE_MARGIN_MM, RANGE_MARGIN_MM, False
    msg = f"no margin is defined for {key!r}"
    raise ValueError(msg)


def endpoint_keys(energy: int) -> list[str]:
    depths = DEPTHS_MM[energy]
    return ([f"lateral_{d}_{o}" for d in depths for o in OFFSETS_MM] + [f"cax_{d}" for d in depths]
            + ["cax_max", "r80_mm", "r20_mm"])


def specs() -> list[Spec]:
    """The 34 endpoints of one contrast, 17 per energy."""
    out = []
    for energy in fr.ENERGIES:
        for key in endpoint_keys(energy):
            lo, hi, log = margin(key)
            out.append(Spec(f"E{energy}/{key}", "E", key, log, lo, hi, energy))
    if len(out) != FAMILY_SIZE or len({s.eid for s in out}) != FAMILY_SIZE:
        msg = f"endpoint family has {len(out)} entries, expected {FAMILY_SIZE} distinct"
        raise AssertionError(msg)
    return out


# --------------------------------------------------------------------------------------- what every run is held to
def blob_id(data: bytes) -> str:
    """The git blob id of these bytes."""
    return hashlib.sha1(b"blob %d\x00" % len(data) + data).hexdigest()


def committed_blobs(commit: str, repo: Path) -> dict[str, str]:
    """path -> git blob id, at `commit`, of every input a run directory copies and every script of interest."""
    paths = [*INPUT_PATHS, *RUN_LIST_FILES, *ANALYSIS_FILES]
    proc = subprocess.run(["git", "-C", str(repo), "ls-tree", "-r", commit, "--", *paths],
                          capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        msg = f"cannot list the tree of {commit}: {proc.stderr.strip()}"
        raise InputError(msg)
    out: dict[str, str] = {}
    for line in proc.stdout.splitlines():
        meta, _, path = line.partition("\t")
        fields = meta.split()
        if len(fields) == 3 and fields[1] == "blob":
            out[path] = fields[2]
    return out


@dataclass
class Expected:
    """Read once: the committed inputs and the generated case."""

    commit: str
    mode: str
    blobs: dict[str, str]
    inputs: dict[str, str]  # run-directory path -> blob id, for the INPUT_PATHS files
    case_text: dict[int, dict[str, bytes]]  # energy -> {cube.mhd, E<energy>_S150.txt}, LF line ends
    cube_raw_sha256: str


def expected(commit: str, mode: str, repo: Path = REPO) -> Expected:
    if not aa.FULL_SHA.match(commit):
        msg = f"acquisition commit must be a full 40-hex sha, got {commit!r}"
        raise InputError(msg)
    if mode not in fr.HISTORIES:
        msg = f"mode must be one of {sorted(fr.HISTORIES)}, got {mode!r}"
        raise InputError(msg)
    if not fm.N == fea.N == psm.N == GRID_N or fea.SP != 2.0:
        msg = f"this analysis is for the {GRID_N}-voxel, 2 mm cube; FE_N gives {fm.N}"
        raise InputError(msg)
    blobs = committed_blobs(commit, repo)
    for rel, module in zip(RUN_LIST_FILES, (fr, fm), strict=True):
        if blobs.get(rel) != blob_id(Path(module.__file__).read_bytes()):
            msg = f"{rel} as this analysis imports it is not the file of the acquisition commit {commit}"
            raise InputError(msg)
    inputs = {p: b for p, b in blobs.items() if p.startswith(("Materials/", "Scanners/default/")) or p == INPUT_PATHS[2]}
    if not any(p.startswith("Materials/") for p in inputs) or INPUT_PATHS[2] not in inputs:
        msg = f"the tree of {commit} lacks Materials or the beam model"
        raise InputError(msg)
    case_text: dict[int, dict[str, bytes]] = {}
    with tempfile.TemporaryDirectory() as tmp, contextlib.chdir(tmp), contextlib.redirect_stdout(None):
        fm.main([float(e) for e in fr.ENERGIES], [fr.SIDE_MM])
        cube_raw = hashlib.sha256(Path("cube.raw").read_bytes()).hexdigest()
        for e in fr.ENERGIES:
            plan = f"{fr.plan_name(e)}.txt"
            case_text[e] = {"cube.mhd": Path("cube.mhd").read_bytes(), plan: Path(plan).read_bytes()}
    return Expected(commit, mode, blobs, inputs, case_text, cube_raw)


def rerun_commit_problems(acquisition: str, rerun: str, repo: Path = REPO) -> list[str]:
    """Why `rerun` cannot stand in for the acquisition commit, if it cannot: it must descend from it and change nothing
    but RERUN_MAY_CHANGE. Then the inputs, run lists, seeds and binaries the rerun ran are the acquisition's."""
    if not aa.FULL_SHA.match(rerun):
        return [f"rerun commit must be a full 40-hex sha, got {rerun!r}"]
    if rerun == acquisition:
        return ["the rerun commit is the acquisition commit"]
    anc = subprocess.run(["git", "-C", str(repo), "merge-base", "--is-ancestor", acquisition, rerun],
                         capture_output=True, text=True, check=False)
    if anc.returncode != 0:
        return [f"{rerun} does not descend from the acquisition commit {acquisition}" + (f" ({anc.stderr.strip()})" if anc.stderr.strip() else "")]
    diff = subprocess.run(["git", "-C", str(repo), "diff", "--name-only", acquisition, rerun],
                          capture_output=True, text=True, check=False)
    if diff.returncode != 0:
        return [f"cannot compare {acquisition} with {rerun}: {diff.stderr.strip()}"]
    extra = sorted(set(diff.stdout.split()) - set(RERUN_MAY_CHANGE))
    return [f"the rerun commit changes files other than the run requests and the plan: {extra[:5]}"] if extra else []


@dataclass
class Rerun:
    """Acquisition amendment 1: the roots holding the rerun of RERUN_CELLS, and what its runs are verified against."""

    roots: list[Path]
    want: Expected


def expected_cfg(arm: str, energy: int, seed: int, mode: str) -> list[str]:
    """cfg.txt as both workflows write it."""
    return [
        f"Num_Threads {fr.ARMS[arm].threads}", f"Num_Primaries {fr.HISTORIES[mode]}", f"RNG_Seed {seed}", "CT_File cube.mhd",
        "HU_Density_Conversion_File Scanners/default/HU_Density_Conversion.txt",
        "HU_Material_Conversion_File Scanners/default/HU_Material_Conversion.txt",
        "BDL_Machine_Parameter_File BDL/BDL_default_DN_RangeShifter.txt",
        f"BDL_Plan_File {fr.plan_name(energy)}.txt", "Output_Directory out_seed", "Dose_MHD_Output True",
    ]


# ------------------------------------------------------------------------------------------------- verification
@dataclass
class Run:
    arm: str
    energy: int
    seed: int
    path: Path
    problems: list[str] = field(default_factory=list)
    record: dict[str, object] = field(default_factory=dict)  # run.json, when it parsed
    dose_sha256: str | None = None
    values: dict[str, float | None] = field(default_factory=dict)  # endpoint key -> value; None = invalid
    peak3: float | None = None

    @property
    def usable(self) -> bool:
        return not self.problems

    @property
    def label(self) -> str:
        return f"{self.arm} {self.energy} MeV s{self.seed}"


def _lf(data: bytes) -> bytes:
    return data.replace(b"\r\n", b"\n")


def _text(data: bytes) -> str:
    """A text file from either host: PowerShell's redirection writes UTF-16 with a byte-order mark."""
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        return data.decode("utf-16")
    return data.decode("utf-8", errors="replace")


def _distinct_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """json.loads keeps the last of a repeated key; a record that repeats one is not resolved either way."""
    if len({k for k, _ in pairs}) != len(pairs):
        msg = "a key is repeated"
        raise ValueError(msg)
    return dict(pairs)


def _run_json_problems(rec: object, arm: str, energy: int, seed: int, want: Expected) -> list[str]:
    if not isinstance(rec, dict):
        return ["run.json is not a JSON object with distinct keys"]
    bad: list[str] = []
    for key, value in (
        ("seed", seed), ("requested", fr.HISTORIES[want.mode]), ("threads", fr.ARMS[arm].threads), ("energy_mev", energy),
        ("case", "E"), ("arm", arm), ("mode", want.mode), ("host", fr.ARMS[arm].host),
        ("binary_sha256", fr.ARMS[arm].binary_sha256), ("commit", want.commit), ("transport_status", "ok"),
    ):
        got = rec.get(key)
        if got != value or type(got) is not type(value):
            bad.append(f"run.json {key} is {got!r}, expected {value!r}")
    simulated = aa._int(rec.get("simulated"))
    if simulated is None or simulated < fr.HISTORIES[want.mode]:
        bad.append(f"run.json simulated is {rec.get('simulated')!r}, expected an integer of at least {fr.HISTORIES[want.mode]}")
    elif rec.get("overshoot") != simulated - fr.HISTORIES[want.mode]:
        bad.append(f"run.json overshoot is {rec.get('overshoot')!r}, expected {simulated - fr.HISTORIES[want.mode]}")
    return bad


def _mhd_problems(text: str) -> list[str]:
    """Dose.mhd must describe the grid this analysis assumes when it reshapes Dose.raw."""
    got: dict[str, str] = {}
    for line in text.splitlines():
        key, sep, value = line.partition("=")
        if sep:
            if key.strip() in got:
                return [f"Dose.mhd sets {key.strip()} more than once"]
            got[key.strip()] = value.strip()
    bad = [f"Dose.mhd {key} is {got.get(key)!r}, expected {value!r}" for key, value in (
        ("NDims", "3"), ("DimSize", f"{GRID_N} {GRID_N} {GRID_N}"), ("ElementType", "MET_FLOAT"),
        ("ElementByteOrderMSB", "False"), ("ElementDataFile", "Dose.raw")) if got.get(key) != value]
    try:
        spacing = [float(v) for v in got.get("ElementSpacing", "").split()]
    except ValueError:
        spacing = []
    if spacing != [fea.SP] * 3:
        bad.append(f"Dose.mhd ElementSpacing is {got.get('ElementSpacing')!r}, expected {fea.SP} mm three times")
    return bad


def simulated_count(log: str) -> tuple[int | None, str]:
    """(count, "") from the one line MCsquare prints when a simulation completes, or (None, why).

    Every line that carries the label is counted, and the whole line must be what compute_simulation.c prints: a
    numeric prefix of something else, or a second and different count, is not resolved in favour of either."""
    lines = [line.strip() for line in log.splitlines() if _COUNT_LABEL in line]
    if len(lines) != 1:
        return None, f"{len(lines)} '{_COUNT_LABEL}' lines, expected exactly 1"
    hit = _COUNT_LINE.fullmatch(lines[0])
    if hit is None:
        return None, f"a malformed '{_COUNT_LABEL}' line: {lines[0][:80]!r}"
    return int(hit[1]), ""


def _sha256_stream(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_run(run: Run, want: Expected, ledger: aa.Ledger) -> None:
    """Fill run.problems. Reads every file of the run except the dose values (Dose.raw is hashed, not interpreted)."""
    d, bad = run.path, run.problems
    needed = ("run.json", "cfg.txt", "log.txt", "cube.mhd", "cube.raw", f"{fr.plan_name(run.energy)}.txt",
              "out_seed/Dose.raw", "out_seed/Dose.mhd", "out_seed/sha256.txt")
    missing = [n for n in needed if not (d / n).is_file()]
    if missing:
        bad.append("missing: " + ", ".join(missing))
        return
    try:
        rec = json.loads(_text(ledger.read_bytes(d / "run.json")), object_pairs_hook=_distinct_keys)
    except ValueError:
        rec = None
    bad += _run_json_problems(rec, run.arm, run.energy, run.seed, want)
    run.record = rec if isinstance(rec, dict) else {}
    if _text(ledger.read_bytes(d / "cfg.txt")).splitlines() != expected_cfg(run.arm, run.energy, run.seed, want.mode):
        bad.append("cfg.txt is not the configuration the workflow writes for this run")
    log = _text(ledger.read_bytes(d / "log.txt"))
    count, why = simulated_count(log)
    if count is None:
        bad.append(f"log.txt has {why}")
    elif count != run.record.get("simulated"):
        bad.append(f"log.txt reports {count} primaries simulated, run.json {run.record.get('simulated')!r}")
    if "Unknown tag" in log:
        bad.append("log.txt reports a configuration tag MCsquare did not recognise")
    for name, text in want.case_text[run.energy].items():
        if _lf(ledger.read_bytes(d / name)) != text:
            bad.append(f"{name} is not what field_edge_make_cases.py writes")
    if hashlib.sha256(ledger.read_bytes(d / "cube.raw")).hexdigest() != want.cube_raw_sha256:
        bad.append("cube.raw is not what field_edge_make_cases.py writes")
    found = {p.relative_to(d).as_posix() for top in ("Materials", "Scanners", "BDL") for p in (d / top).rglob("*") if p.is_file()}
    if found != set(want.inputs):
        extra, lacking = sorted(found - set(want.inputs)), sorted(set(want.inputs) - found)
        bad.append(f"inputs differ from the committed set: {len(extra)} added {extra[:3]}, {len(lacking)} missing {lacking[:3]}")
    changed = sorted(p for p in found & set(want.inputs) if blob_id(ledger.read_bytes(d / p)) != want.inputs[p])
    if changed:
        bad.append(f"{len(changed)} input file(s) differ from the committed bytes: {changed[:3]}")
    lines = [m.group(2, 1) for line in _text(ledger.read_bytes(d / "out_seed/sha256.txt")).splitlines()
             if (m := _SHA_LINE.match(line.strip()))]
    if sorted(name for name, _ in lines) != ["Dose.mhd", "Dose.raw"]:  # each once: a repeated name is not resolved
        bad.append(f"out_seed/sha256.txt lists {sorted(name for name, _ in lines)}, expected Dose.mhd and Dose.raw once each")
        return
    recorded = dict(lines)
    mhd = ledger.read_bytes(d / "out_seed/Dose.mhd")
    if hashlib.sha256(mhd).hexdigest() != recorded["Dose.mhd"]:
        bad.append("Dose.mhd is not the file hashed on the run host")
    bad += _mhd_problems(_text(mhd))
    if (d / "out_seed/Dose.raw").stat().st_size != 4 * GRID_N**3:
        bad.append(f"Dose.raw is {(d / 'out_seed/Dose.raw').stat().st_size} bytes, expected {4 * GRID_N**3}")
    run.dose_sha256 = _sha256_stream(d / "out_seed/Dose.raw")
    ledger.files[(d / "out_seed/Dose.raw").relative_to(ledger.root).as_posix()] = run.dose_sha256
    if run.dose_sha256 != recorded["Dose.raw"]:
        bad.append("Dose.raw is not the file hashed on the run host")


def endpoint_values(ep: dict[str, object], energy: int) -> dict[str, float | None]:
    """The 17 endpoint values of one run, None where the plan calls the value invalid: a central-axis dose that is
    not finite and positive, a range without exactly one crossing, a lateral value that is not finite. A lateral
    value of exactly zero is valid."""
    out: dict[str, float | None] = {}
    for key in endpoint_keys(energy):
        v = aa._num(ep.get(key))
        if key.startswith("cax_") and (ep.get(f"{key}_invalid") is True or v is None or v <= 0):
            v = None
        if key in ("r80_mm", "r20_mm") and ep.get(key.replace("_mm", "_crossings")) != 1:
            v = None
        out[key] = v
    return out


def measure(run: Run) -> None:
    """Compute the endpoints of a verified run from its Dose.raw: the bytes analysed are the bytes that hash."""
    data = (run.path / "out_seed/Dose.raw").read_bytes()
    if hashlib.sha256(data).hexdigest() != run.dose_sha256:
        run.problems.append("Dose.raw changed between verification and reading")
        return
    dose = np.frombuffer(data, dtype="<f4").reshape(GRID_N, GRID_N, GRID_N)
    run.values = endpoint_values(psm.field_endpoints(dose, fr.SIDE_MM, DEPTHS_MM[run.energy], OFFSETS_MM), run.energy)
    run.peak3 = psm.peak_three_row_mean(dose)


Cells = dict[tuple[str, int], list[Run]]
Issues = dict[tuple[str, int], list[str]]


def _run_like(entry: Path) -> bool:
    """Looks like run data: named s<digit>..., or a directory with a run's files anywhere beneath it."""
    return bool(_RUN_NAME.match(entry.name)) or (entry.is_dir() and any(p.name in _RUN_FILES for p in entry.rglob("*")))


def _stray(entries: list[Path], where: str, notes: list[str]) -> list[str]:
    """Entries that have no place in the layout: those that look like run data are returned as population issues
    (a run that is not a run of the frozen lists must not be silently left out), the rest are noted as ignored."""
    issues = []
    for e in entries:
        if _run_like(e):
            issues.append(f"{where}: {e.name} looks like run data and is not a run of the frozen lists")
        else:
            notes.append(f"{where}: entry {e.name} is ignored (not run data)")
    return issues


def load(roots: list[Path], want: Expected, rerun: Rerun | None = None) -> tuple[Cells, Issues, aa.Ledger, list[str]]:
    """Every run of the frozen lists, VERIFIED. No dose value is read here. Refuses if an arm is in no root.

    With `rerun`, the cells of RERUN_CELLS come from the rerun roots, verified against the rerun commit, and the
    interrupted cells of the original roots are not read. Anything else in a rerun root is a population issue."""
    for root in [*roots, *(rerun.roots if rerun else [])]:
        if not root.is_dir():
            msg = f"{root} is not a directory"
            raise InputError(msg)
    base = Path(os.path.commonpath([str(r.resolve()) for r in [*roots, *(rerun.roots if rerun else [])]]))
    ledger = aa.Ledger(base)
    home: dict[str, Path] = {}
    for arm, spec in fr.ARMS.items():
        holders = [r.resolve() for r in roots if (r / spec.directory).is_dir()]
        if len(holders) != 1:
            msg = (f"arm {arm} (directory {spec.directory}) is in {len(holders)} of the roots, need exactly 1: both hosts' "
                   "trees must be present before any run is read")
            raise InputError(msg)
        home[arm] = holders[0] / spec.directory
    notes: list[str] = []
    arm_dirs = {a.directory for a in fr.ARMS.values()}
    cell_names = {f"e{e}" for e in fr.ENERGIES}
    outside: dict[str, list[str]] = {}  # arm -> issues from entries outside its cells; they touch both energies
    for root in roots:
        stray = _stray([e for e in sorted(root.iterdir()) if e.name not in arm_dirs], str(root), notes)
        for arm in (a for a in fr.ARMS if home[a].parent == root.resolve()):
            outside.setdefault(arm, []).extend(stray)
    for arm in fr.ARMS:
        inside = [e for e in sorted(home[arm].iterdir()) if not (e.name in cell_names and e.is_dir())]
        outside.setdefault(arm, []).extend(_stray(inside, f"{arm} ({home[arm].name})", notes))
    rerun_home: dict[str, Path] = {}
    rerun_outside: dict[str, list[str]] = {}
    if rerun is not None:
        rerun_arms = {a for a, _ in RERUN_CELLS}
        for arm in sorted(rerun_arms):
            holders = [r.resolve() for r in rerun.roots if (r / fr.ARMS[arm].directory).is_dir()]
            if len(holders) != 1:
                msg = f"rerun arm {arm} is in {len(holders)} of the rerun roots, need exactly 1"
                raise InputError(msg)
            rerun_home[arm] = holders[0] / fr.ARMS[arm].directory
        own = {fr.ARMS[a].directory for a in rerun_arms}
        for root in rerun.roots:
            stray = _stray([e for e in sorted(root.iterdir()) if e.name not in own], f"rerun {root}", notes)
            stray += [f"rerun {root}: {e.name} is not a rerun arm" for e in sorted(root.iterdir())
                      if e.name in arm_dirs - own and e.is_dir()]
            for arm in (a for a in rerun_arms if rerun_home[a].parent == root.resolve()):
                rerun_outside.setdefault(arm, []).extend(stray)
        for arm in sorted(rerun_arms):
            mine = {f"e{e}" for a, e in RERUN_CELLS if a == arm}
            inside = [e for e in sorted(rerun_home[arm].iterdir()) if not (e.name in mine and e.is_dir())]
            rerun_outside.setdefault(arm, []).extend(_stray(inside, f"rerun {arm} ({rerun_home[arm].name})", notes))
            rerun_outside[arm] += [f"rerun {arm}: {e.name} is not a rerun cell" for e in inside if e.name in cell_names]
    cells: Cells = {}
    issues: Issues = {}
    for arm in fr.ARMS:
        for energy in fr.ENERGIES:
            w = want
            if rerun is not None and (arm, energy) in RERUN_CELLS:
                w = rerun.want
                cell, runs, why = rerun_home[arm] / f"e{energy}", [], list(outside[arm]) + rerun_outside.get(arm, [])
                notes.append(f"{arm} {energy} MeV: read from the rerun at {w.commit[:12]} (acquisition amendment 1); "
                             f"the cell of {want.commit[:12]} is not read")
            else:
                cell, runs, why = home[arm] / f"e{energy}", [], list(outside[arm])
            seeds = fr.seeds(arm, energy, want.mode)
            entries = sorted(cell.iterdir()) if cell.is_dir() else []
            present = {int(m[1]) for e in entries if e.is_dir() and (m := _SEED_DIR.match(e.name))}
            why += [f"{arm} {energy} MeV: run s{s} is missing" for s in seeds if s not in present]
            why += [f"{arm} {energy} MeV: run s{s} is not in the frozen list" for s in sorted(present - set(seeds))]
            others = [e for e in entries if not (e.is_dir() and _SEED_DIR.match(e.name))]  # a known entry holding run data is stray too
            why += _stray([e for e in others if e.name not in CELL_ENTRIES or _run_like(e)], f"{arm} {energy} MeV", notes)
            for rel in RUN_LIST_FILES:
                snap = cell / "snapshot" / rel
                if not snap.is_file() or blob_id(ledger.read_bytes(snap)) != w.blobs[rel]:
                    why.append(f"{arm} {energy} MeV: the run host's snapshot of {rel} is not the acquisition commit's")
            for s in seeds:
                if s not in present:
                    continue
                run = Run(arm, energy, s, cell / f"s{s}")
                verify_run(run, w, ledger)
                why += [f"{run.label}: {p}" for p in run.problems]
                runs.append(run)
            cells[(arm, energy)], issues[(arm, energy)] = runs, why
    first: dict[str, Run] = {}  # a dose file that two runs share is a duplicated run, whatever their seeds say
    for run in (r for runs in cells.values() for r in runs if r.dose_sha256 is not None):
        other = first.setdefault(run.dose_sha256, run)  # type: ignore[arg-type]
        if other is not run:
            for r in (run, other):
                issues[(r.arm, r.energy)].append(f"{run.label} and {other.label} have the same Dose.raw")
    return cells, issues, ledger, notes


def measure_all(cells: Cells, issues: Issues) -> None:
    """Compute the endpoints of every verified run. Called only after load() has seen the whole population."""
    for (arm, energy), runs in cells.items():
        for run in (r for r in runs if r.usable):
            measure(run)
            issues[(arm, energy)] += [f"{run.label}: {p}" for p in run.problems]  # only what measuring found


# --------------------------------------------------------------------------------------------------- statistics
def _values(runs: list[Run], key: str) -> tuple[list[float] | None, str]:
    bad = [r.label for r in runs if not r.usable or r.values.get(key) is None]
    if bad:
        return None, f"invalid or unusable in {len(bad)} of {len(runs)} runs: " + ", ".join(bad[:3]) + (" ..." if len(bad) > 3 else "")
    if len(runs) < 2:
        return None, f"{len(runs)} run(s), need at least 2"
    return [float(r.values[key]) for r in runs], ""  # type: ignore[arg-type]


def _simulated(runs: list[Run]) -> int:
    return sum(aa._int(r.record.get("simulated")) or 0 for r in runs)


def evaluate(spec: Spec, arm: list[Run], ref: list[Run], *, descriptive: bool) -> dict[str, object]:
    """One endpoint of one contrast, as a JSON-ready row."""
    row: dict[str, object] = {"endpoint": spec.eid, "energy": spec.energy, "key": spec.key, "scale": "ratio" if spec.log else "difference",
                              "margin": [math.exp(spec.lo), math.exp(spec.hi)] if spec.log else [spec.lo, spec.hi],
                              "all_runs_zero": False}
    va, why_a = _values(arm, spec.key)
    vb, why_b = _values(ref, spec.key)
    if va is None or vb is None:
        reason = "; ".join(w for w in (f"arm: {why_a}" if why_a else "", f"reference: {why_b}" if why_b else "") if w)
        return {**row, "outcome": "not_established", "reason": reason, "p_tost": None}
    row["arm_mean"], row["reference_mean"] = st.mean(va), st.mean(vb)
    a, b = ([math.log(v) for v in va], [math.log(v) for v in vb]) if spec.log else (va, vb)
    d, se, df = aa.welch(a, b)
    back = math.exp if spec.log else float
    row["estimate"] = back(d)
    if se == 0.0:  # zero sample variance in BOTH arms: not zero variance (plan; #55 review 7098)
        zero = all(v == 0.0 for v in va + vb)
        return {**row, "outcome": "not_established", "p_tost": None, "all_runs_zero": zero,
                "reason": ("all runs zero in both arms" if zero else "zero sample variance in both arms")
                + f" (simulated histories: arm {_simulated(arm)}, reference {_simulated(ref)})"}
    ci90, ci95 = aa.interval(d, se, df, 0.90), aa.interval(d, se, df, 0.95)
    row.update(se=se, df=df, ci90=[back(v) for v in ci90], ci95=[back(v) for v in ci95], reason="")
    if descriptive:
        return {**row, "outcome": "descriptive", "p_tost": None}
    return {**row, "outcome": aa.classify(ci90, spec.lo, spec.hi), "p_tost": aa.tost_p(d, se, df, spec.lo, spec.hi)}


def ratio_of_means(a: list[float], b: list[float]) -> dict[str, object]:
    """Descriptive: mean(a) / mean(b) with a 95% interval from the delta method on the log of each mean
    (standard error of a log mean = SE of the mean / mean), Welch-Satterthwaite degrees of freedom."""
    ma, mb = st.mean(a), st.mean(b)
    if not (ma > 0 and mb > 0):
        return {"status": "not given", "reason": "a mean is not positive", "arm_mean": ma, "reference_mean": mb}
    ra, rb = st.variance(a) / len(a) / ma**2, st.variance(b) / len(b) / mb**2
    out: dict[str, object] = {"status": "ratio", "estimate": ma / mb, "arm_mean": ma, "reference_mean": mb}
    if ra + rb == 0.0:
        return {**out, "status": "ratio without interval", "reason": "zero sample variance in both arms"}
    se = math.sqrt(ra + rb)
    df = (ra + rb) ** 2 / (ra**2 / (len(a) - 1) + rb**2 / (len(b) - 1))
    h = float(student_t.ppf(0.975, df)) * se
    return {**out, "se_log": se, "df": df, "ci95": [ma / mb * math.exp(-h), ma / mb * math.exp(h)]}


def peak_ratio(arm: list[Run], ref: list[Run]) -> dict[str, object]:
    """Descriptive: geometric-mean ratio of the three-row mean around each run's maximum row, with a 95% interval."""
    if any(not r.usable or r.peak3 is None for r in arm + ref) or len(arm) < 2 or len(ref) < 2:
        return {"status": "not given", "reason": "unusable run, or a maximum at the edge of the grid"}
    d, se, df = aa.welch([math.log(r.peak3) for r in arm], [math.log(r.peak3) for r in ref])  # type: ignore[arg-type]
    if se == 0.0:
        return {"status": "ratio without interval", "estimate": math.exp(d), "reason": "zero sample variance in both arms"}
    lo, hi = aa.interval(d, se, df, 0.95)
    return {"status": "ratio", "estimate": math.exp(d), "ci95": [math.exp(lo), math.exp(hi)],
            "interval_within_margin": CAX_RATIO[0] < math.exp(lo) and math.exp(hi) < CAX_RATIO[1]}


def analyse_contrast(cells: Cells, issues: Issues, arm: str, ref: str, *, confirmatory: bool) -> dict[str, object]:
    partial = [why for name in (arm, ref) for e in fr.ENERGIES for why in issues[(name, e)]]
    withheld = confirmatory and bool(partial)
    rows = [evaluate(s, cells[(arm, s.energy)], cells[(ref, s.energy)], descriptive=not confirmatory or withheld)  # type: ignore[index]
            for s in specs()]
    claims: dict[str, object] = {}
    if confirmatory and not withheld:
        for row, hp in zip(rows, aa.holm_adjust([r["p_tost"] for r in rows]), strict=True):  # type: ignore[misc]
            row["holm_p"] = hp
            row["holm_decision"] = ("not_established" if row["outcome"] == "not_established"
                                    else "equivalent" if hp is not None and hp < ALPHA else "not shown")
        zero = [r["endpoint"] for r in rows if r["all_runs_zero"]]
        claims = {
            "joint_claim_equivalent_on_all_endpoints": all(r["outcome"] == "equivalent" for r in rows),
            "secondary_joint_claim_holm_excluding_all_zero_rows": all(r["holm_decision"] == "equivalent" for r in rows if not r["all_runs_zero"]),
            "all_zero_rows_excluded_from_the_secondary_claim": zero,
            "n_equivalent_unadjusted": sum(r["outcome"] == "equivalent" for r in rows),
            "n_equivalent_holm": sum(r["holm_decision"] == "equivalent" for r in rows),
            "n_not_established": sum(r["outcome"] == "not_established" for r in rows),
        }
    ratios, maxima = [], []
    by_id = {r["endpoint"]: r for r in rows}
    for energy in fr.ENERGIES:
        a, b = cells[(arm, energy)], cells[(ref, energy)]
        for depth in DEPTHS_MM[energy]:
            for off in RATIO_OFFSETS_MM:
                key = f"lateral_{depth}_{off}"
                va, _ = _values(a, key)
                vb, _ = _values(b, key)
                ratios.append({"endpoint": f"E{energy}/{key}", **(ratio_of_means(va, vb) if va is not None and vb is not None
                                                                   else {"status": "not given", "reason": "invalid or unusable run"})})
        top, r80 = by_id[f"E{energy}/cax_max"], by_id[f"E{energy}/r80_mm"]
        maxima.append({"energy": energy, "maximum_outcome": top["outcome"], "maximum_ratio": top.get("estimate"),
                       "maximum_ci90": top.get("ci90"), "r80_difference_mm": r80.get("estimate"), "r80_ci90": r80.get("ci90"),
                       "three_row_mean_ratio": peak_ratio(a, b)})
    return {"arm": arm, "reference": ref, "confirmatory": confirmatory, "partial_reasons": partial, "claims_withheld": withheld,
            "claims": claims, "rows": rows, "ratio_of_arm_means_at_50_and_70_mm": ratios, "maximum": maxima}


# ----------------------------------------------------------------------------------------------------- reporting
def _num(x: object, digits: int = 4) -> str:
    """A number for the table: `digits` decimals, or three significant figures where that many decimals would show
    fewer than three, so that a very small mean is not printed as 0.0000 (#60). Only an exact zero prints as 0."""
    if not isinstance(x, (int, float)) or isinstance(x, bool):
        return ""
    if x == 0:
        return "0"
    return f"{x:.{digits}f}" if abs(x) >= 10.0 ** (2 - digits) else f"{x:.3g}"


def _pair(ci: object, digits: int = 4) -> str:
    return f"[{_num(ci[0], digits)}, {_num(ci[1], digits)}]" if isinstance(ci, list) else ""


def markdown(doc: dict[str, object]) -> list[str]:
    md = [f"# Same-host comparison, part E: field edge at 100 and 150 MeV ({doc['acquisition_commit']})", "", str(doc["label"]), ""]
    col = doc["collection"]
    if col["population_issues_before_any_dose_was_read"]:  # type: ignore[index]
        md += [(f"**The population was not complete and verified ({col['population_issues_before_any_dose_was_read']} issue(s)). "  # type: ignore[index]
                f"The collection was declared closed: {col['declared_closed']}** The contrasts the issues touch are PARTIAL."), ""]  # type: ignore[index]
    for c in doc["contrasts"]:  # type: ignore[union-attr]
        kind = "confirmatory" if c["confirmatory"] else "descriptive"
        md += [f"## {c['arm']} against {c['reference']} ({kind})", ""]
        if c["claims_withheld"]:
            md += ["**PARTIAL: descriptive estimates only; no joint claim and no Holm decision.**", ""]
            md += [f"- {why}" for why in c["partial_reasons"]] + [""]
        elif c["confirmatory"]:
            k = c["claims"]
            md += [f"- Joint claim (all {FAMILY_SIZE} endpoints equivalent): **{'ESTABLISHED' if k['joint_claim_equivalent_on_all_endpoints'] else 'NOT ESTABLISHED'}**",
                   (f"- Secondary joint claim (every endpoint other than rows where every run of both arms is zero is equivalent "
                   f"after Holm over all {FAMILY_SIZE}): "
                   f"**{'ESTABLISHED' if k['secondary_joint_claim_holm_excluding_all_zero_rows'] else 'NOT ESTABLISHED'}**; "
                   f"rows excluded: {', '.join(k['all_zero_rows_excluded_from_the_secondary_claim']) or 'none'}"),
                   (f"- Equivalent: {k['n_equivalent_unadjusted']} of {FAMILY_SIZE} unadjusted, {k['n_equivalent_holm']} after Holm; "
                   f"not established: {k['n_not_established']}"), ""]
        md += ["| endpoint | arm mean | reference mean | estimate | 90% interval | margin | outcome | Holm p |", "|---|---|---|---|---|---|---|---|"]
        for r in c["rows"]:
            outcome = r["outcome"] + (f" ({r['reason']})" if r["reason"] else "")
            md.append(f"| {r['endpoint']} | {_num(r.get('arm_mean'))} | {_num(r.get('reference_mean'))} | {_num(r.get('estimate'))} | "
                      f"{_pair(r.get('ci90'))} | {_pair(r['margin'])} ({r['scale']}) | {outcome} | {_num(r.get('holm_p'))} |")
        md += ["", "Descriptive only. Ratio of arm means (arm / reference) at 50 and 70 mm outside the edge, 95% interval:", "",
               "| endpoint | arm mean | reference mean | ratio | 95% interval | note |", "|---|---|---|---|---|---|"]
        md += [f"| {r['endpoint']} | {_num(r.get('arm_mean'))} | {_num(r.get('reference_mean'))} | {_num(r.get('estimate'))} | "
               f"{_pair(r.get('ci95'))} | {r.get('reason', '')} |" for r in c["ratio_of_arm_means_at_50_and_70_mm"]]
        md += ["", ("Descriptive only. Beside each maximum: the R80 difference and the ratio of the three-row mean around each run's "
               "maximum row. No cause is attributed."), "",
               "| energy | maximum: outcome | maximum: ratio, 90% | R80 difference (mm), 90% | three-row mean: ratio, 95% |", "|---|---|---|---|---|"]
        for m in c["maximum"]:
            t = m["three_row_mean_ratio"]
            md.append(f"| {m['energy']} MeV | {m['maximum_outcome']} | {_num(m['maximum_ratio'], 5)} {_pair(m['maximum_ci90'], 5)} | "
                      f"{_num(m['r80_difference_mm'])} {_pair(m['r80_ci90'])} | {_num(t.get('estimate'), 5)} {_pair(t.get('ci95'), 5)} "
                      f"{t.get('reason', '')} |")
            if c["confirmatory"] and not c["claims_withheld"] and m["maximum_outcome"] in ("inconclusive", "not_equivalent") \
                    and t.get("interval_within_margin") is True:
                md += ["", (f"At {m['energy']} MeV the maximum is {m['maximum_outcome'].replace('_', ' ')} while the three-row ratio's 95% "
                       "interval lies within [0.995, 1.005]: compatible with a shift of the peak against the 2 mm rows and with a "
                       "difference in peak dose; these runs do not distinguish them. The row stands as it is."), ""]
        md.append("")
    return md


def _own_blobs(repo: Path, at_acquisition: dict[str, str]) -> dict[str, object]:
    out = {}
    for rel in ANALYSIS_FILES:
        own = blob_id((repo / rel).read_bytes()) if (repo / rel).is_file() else None
        out[rel] = {"blob": own, "at_acquisition": at_acquisition.get(rel), "unchanged": own == at_acquisition.get(rel)}
    return out


def make_rerun(commit: str, want: Expected, rerun_roots: list[Path] | None, rerun_commit: str | None,
               repo: Path = REPO) -> Rerun | None:
    """The rerun of acquisition amendment 1, or None. Refuses a rerun commit that is not the acquisition's inputs."""
    if (rerun_roots is None) != (rerun_commit is None):
        msg = "--rerun-root and --rerun-commit go together"
        raise InputError(msg)
    if rerun_roots is None or rerun_commit is None:
        return None
    bad = rerun_commit_problems(commit, rerun_commit, repo)
    if bad:
        raise InputError("; ".join(bad))
    return Rerun(list(rerun_roots), replace(want, commit=rerun_commit))


def run(roots: list[Path], commit: str, *, closed: str | None = None, repo: Path = REPO,
        rerun_roots: list[Path] | None = None, rerun_commit: str | None = None) -> dict[str, object]:
    """Verify the whole population, and only then read dose values.

    With any population issue the analysis refuses before a dose value is read, unless `closed` states who declares the
    acquisition closed and on what evidence. That statement is a declaration this script cannot check; it is recorded."""
    want = expected(commit, "full", repo)
    rerun = make_rerun(commit, want, rerun_roots, rerun_commit, repo)
    cells, issues, ledger, notes = load(roots, want, rerun)
    before = [why for cell in issues.values() for why in cell]
    if closed is not None and not closed.strip():
        msg = "--collection-closed needs a statement: who declares the acquisition closed, and the evidence"
        raise InputError(msg)
    if before and closed is None:
        msg = (f"the frozen population is not complete and verified ({len(before)} issue(s)); no dose value was read. "
               "If the acquisition is closed and this is what it left, say so with --collection-closed: "
               + "; ".join(before[:6]) + (" ..." if len(before) > 6 else ""))
        raise InputError(msg)
    measure_all(cells, issues)
    return {
        "label": LABEL, "plan": "docs/field_100_150_plan.md", "acquisition_commit": commit, "notes": notes,
        "rerun": None if rerun is None else {"amendment": "acquisition amendment 1", "commit": rerun.want.commit,
                                             "cells": [f"{a} {e} MeV" for a, e in RERUN_CELLS]},
        "collection": {"state": "closed by declaration, with population issues" if before else "complete and verified",
                       "population_issues_before_any_dose_was_read": len(before), "declared_closed": closed},
        "analysis_files": _own_blobs(repo, want.blobs), "dataset_fingerprint": ledger.fingerprint(),
        "histories_per_run": fr.HISTORIES["full"], "runs_per_cell": fr.RUNS_PER_CELL,
        "binaries": {a: s.binary_sha256 for a, s in fr.ARMS.items()},
        "population": {f"{a} {e} MeV": {"expected": fr.RUNS_PER_CELL, "verified": sum(r.usable for r in cells[(a, e)]),
                                        "issues": issues[(a, e)]} for a, e in cells},
        "runs": [{"arm": r.arm, "energy": r.energy, "seed": r.seed, "usable": r.usable, "simulated": r.record.get("simulated"),
                  "wall_s": r.record.get("wall_s"), "dose_sha256": r.dose_sha256} for runs in cells.values() for r in runs],
        "contrasts": [analyse_contrast(cells, issues, a, b, confirmatory=c) for a, b, c in fr.CONTRASTS],
    }


def status(roots: list[Path], commit: str, mode: str, *, repo: Path = REPO,
           rerun_roots: list[Path] | None = None, rerun_commit: str | None = None) -> tuple[bool, list[str]]:
    """Verification only: no dose value is read. Returns (every run of the frozen list verified, report lines)."""
    want = expected(commit, mode, repo)
    rerun = make_rerun(commit, want, rerun_roots, rerun_commit, repo)
    cells, issues, _ledger, notes = load(roots, want, rerun)
    n = len(fr.seeds(next(iter(fr.ARMS)), fr.ENERGIES[0], mode))
    lines = [f"{a} {e} MeV: {sum(r.usable for r in runs)} of {n} verified" for (a, e), runs in cells.items()]
    lines += [f"  {why}" for cell in issues.values() for why in cell] + [f"note: {x}" for x in notes]
    return not any(issues.values()), lines


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", type=Path, action="append", required=True, help="directory holding arm directories; repeatable")
    p.add_argument("--commit", required=True, help="the acquisition commit, full sha")
    p.add_argument("--mode", choices=sorted(fr.HISTORIES), default="full")
    p.add_argument("--status", action="store_true", help="verify the runs only; no dose value is read")
    p.add_argument("--rerun-root", type=Path, action="append",
                   help="directory holding the rerun of the cells cut by the timeout (acquisition amendment 1); repeatable")
    p.add_argument("--rerun-commit", help="the commit the rerun ran at, full sha (with --rerun-root)")
    p.add_argument("--json", type=Path)
    p.add_argument("--md", type=Path)
    p.add_argument("--collection-closed", metavar="STATEMENT",
                   help="needed to analyse a population with issues: who declares the acquisition closed and the evidence "
                        "(workflow run ids, archive hashes); written into the output")
    a = p.parse_args(argv)
    try:
        if a.status:
            ok, lines = status(a.root, a.commit, a.mode, rerun_roots=a.rerun_root, rerun_commit=a.rerun_commit)
            print("\n".join(lines))
            print(f"{'VERIFIED' if ok else 'NOT VERIFIED'}: mode {a.mode}; no dose value was read")
            return 0 if ok else 1
        if a.mode != "full":
            p.error("smoke runs are checked for completion only: use --status")
        if a.json is None or a.md is None:
            p.error("--json and --md are required unless --status is given")
        doc = run(a.root, a.commit, closed=a.collection_closed, rerun_roots=a.rerun_root, rerun_commit=a.rerun_commit)
    except InputError as e:
        print(f"field_followup_analyse.py: REFUSED: {e}", file=sys.stderr)
        return 2
    a.json.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    a.md.write_text("\n".join(markdown(doc)) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
