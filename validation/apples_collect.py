"""Collect the same-host comparison workflows' native run trees into the input layout of apples_analyse.py.

Usage:
    python apples_collect.py --source <ROOT> [--source <ROOT> ...] --out <DIR> [--parts A,B] [--commit <40-hex sha>]

Each <ROOT> is a workflow's own acquisition directory for the full runs, `.../<sha12>` with sha12 the first 12 characters
of the acquisition commit (default 2f9dab404cea02f5072352f7ae5773dca27eb2a1): C:\\mcsq-win\\ts\\<sha12> on the Lenovo
(arms aup, apt) and $HOME/fe-study/ap/<sha12> on the HP (arms bup, bpg, bpi), each copied somewhere readable. The native
layout, read from .gitea/workflows/apples-a-windows.yml and apples-b-linux.yml, is

    <ROOT>/<aup|apt|bup|bpg|bpi>/<p|f>/
        run_root.txt  snapshot/  snapshot_sha256.txt
        part A only:  inputs.tar (the git archive the snapshot is extracted from; its files must equal the snapshot's)
        case p:  e<100|150|200>/ (case files, then s<seed>/{run.json, endpoints.json, config.txt, log.txt, out/})
                 e<E>_inputs_sha256.txt
        case f:  fcase/  fcase_sha256.txt
                 s<seed>/{run.json, record.json, cfg.txt, log.txt, out_seed/, cube.*, E200_S150.txt, Materials, ...}

The analysis does not read this layout; this collector builds the one it does read, <DIR>/<arm>/<case>/s<seed>/, BY
COPYING. It never moves, edits or deletes anything under a source, and it writes only under <DIR>, which must not
exist. Nothing is written unless every integrity check passes; any integrity failure refuses with all the reasons found.

A run that FAILED is an acquisition outcome, not an integrity problem, and does not refuse. Every frozen (arm, case,
seed) is classified as one of:
    collected  run.json transport_status ok, Dose files hashing as out/sha256.txt says, and a usable endpoint record;
    failed     (interrupted) a seed directory with no run.json: the job was cancelled or the host died mid-run;
               (transport) run.json transport_status != ok: write_run wrote it and exited the job;
               (endpoint) transport ok, but case p endpoints.json missing, empty or not one ENDPOINTS line, or case f
               record.json missing, unparseable or endpoint_status != ok;
    absent     no seed directory, and every seed of the job's frozen loop after it has none either (the job exited).
A failed run's diagnostic files (run.json, the endpoint record, the config, the log and out/sha256.txt, whichever exist)
are copied to <DIR>/.not_established/<arm>/<case>/s<seed>/, never to the <arm>/<case>/s<seed>/ path the analysis reads
as a run. The manifest's `not_established` list gives arm, case, seed, frozen energy, outcome, reason and those files.
The classification is not a route around the checks: a failed run's run.json must still name its directory's arm, case,
seed and energy and the snapshot's binary, its config its seed, and it may hold nothing unexpected; a run whose
transport_status is ok must have its Dose files and out/sha256.txt (only the endpoint step may fail after ok); a non-ok
run must not have an out/sha256.txt or endpoint record (the workflows write both only after ok); a missing seed directory
before one that exists (a hole, not a stopped job) refuses.

What it checks before copying (all against the frozen run lists at the acquisition commit, apples_analyse.expected_
population): the source directory is named for the commit; only the five arms and two cases appear, nothing else at any
level of the tree (no smoke/ directory, no stray file, no unexpected seed or alias, no product of a collected run
missing, no expected seed missing except as an absent suffix);
run_root.txt names the arm, case, mode full and the commit; snapshot_sha256.txt lists exactly the snapshot (files under
__pycache__, which the runs themselves create after the hashes are taken, are tolerated and recorded);
snapshot/binary_sha256.txt equals every run's binary_sha256; the case-input hash lists match the files; each seed
directory sits under the energy the frozen list gives it; run.json agrees with the directory, with the config
(threads, primaries, seed) and with the log's primaries count (when run.json simulated is set); the Dose files hash as out/sha256.txt (and as record.json
for case f) says; record.json's config, cube/plan and binary hashes agree with the files and snapshot.

What it copies, per collected run: run.json, the endpoint record, the config, the log and out/sha256.txt (as out_sha256.txt);
per arm and case, under .provenance/: run_root.txt, the *_sha256.txt lists and the snapshot's build.txt and
binary_sha256.txt. The Dose files, case inputs and the rest of the snapshot are verified in place and not copied.
collection_manifest.json records the source path and sha256 of every copied file, the roots, the commit, and what was
verified or tolerated; apples_analyse.py verifies it again when it reads the tree.

Prints `collected N, failed N, absent N`. Exit status: 0 collected (whatever the outcomes); 2 refused (reasons on stderr).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import tarfile
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import apples_analyse as aa

ARM_DIRS = {"aup": "A-up", "apt": "A-port", "bup": "B-up", "bpg": "B-pgcc", "bpi": "B-picc"}
CASE_DIRS = {"p": "P", "f": "F"}
PART_OF_ARM = {arm: part for part, arms in aa.PART_ARMS.items() for arm in arms}
SEED_NAME = re.compile(r"^s(\d{6})$")
RUN_ROOT = re.compile(r"^(?P<arm>\S+) case (?P<case>[PF]) mode (?P<mode>\w+), commit (?P<commit>\S+), run (?P<run>\d+)$")
PRIMARIES = re.compile(r"Nbr primaries simulated: (\d+)")
TOLERATED_DIR = "__pycache__"
P_SEED_ENTRIES = {"run.json", "endpoints.json", "config.txt", "log.txt", "out"}
F_SEED_ENTRIES = {"run.json", "record.json", "cfg.txt", "log.txt", "out_seed", "cube.mhd", "cube.raw", "E200_S150.txt",
                  "Materials", "Scanners", "BDL"}
F_CASE_FILES = ("cube.mhd", "cube.raw", "E200_S150.txt")
DOSE_FILES = ("Dose.raw", "Dose.mhd")
WINDOWS_SNAPSHOT_TAR = "inputs.tar"
SNAPSHOT_ADDED_FILES = frozenset({"build.txt", "binary_sha256.txt"})


class CollectionError(Exception):
    """The native trees are not exactly the expected acquisition; nothing was written."""


@dataclass
class Collection:
    """What the checks found and what will be copied. `copies` is (source file, path relative to the output)."""

    commit: str
    frozen: aa.Frozen
    problems: list[str] = field(default_factory=list)
    copies: list[tuple[Path, str]] = field(default_factory=list)
    tolerated: list[str] = field(default_factory=list)
    # failed and absent runs: arm, case, seed, energy, outcome, reason and `copies` (source, path under .not_established/)
    not_established: list[dict[str, object]] = field(default_factory=list)
    verified: dict[str, int] = field(default_factory=lambda: {"runs": 0, "hash_lists": 0, "files_hashed": 0})

    def bad(self, where: object, what: str) -> None:
        self.problems.append(f"{where}: {what}")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def normalise(rel: str) -> str:
    """A path as the workflows list it (backslashes on Windows, a leading ./ from find) as plain forward slashes."""
    rel = rel.strip().replace("\\", "/")
    while rel.startswith("./"):
        rel = rel[2:]
    return rel


def parse_hash_list(text: str) -> dict[str, str]:
    """`<sha256>  <path>` lines (sha256sum, or the Windows Get-FileHash form); refuses a malformed line or a repeat."""
    out: dict[str, str] = {}
    for n, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        m = re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", line.rstrip("\r"))
        if m is None:
            msg = f"line {n} is not '<sha256>  <path>'"
            raise ValueError(msg)
        rel = normalise(m.group(2))
        if rel in out:
            msg = f"{rel} listed twice"
            raise ValueError(msg)
        out[rel] = m.group(1)
    if not out:
        msg = "empty hash list"
        raise ValueError(msg)
    return out


def tree_hashes(root: Path, skip_top: frozenset[str] = frozenset()) -> dict[str, str]:
    """sha256 of every file under root (relative posix paths), without the top-level entries in skip_top."""
    out = {}
    for f in sorted(root.rglob("*")):
        rel = f.relative_to(root).as_posix()
        if f.is_file() and rel.split("/")[0] not in skip_top:
            out[rel] = sha256_file(f)
    return out


def _in_pycache(rel: str) -> bool:
    return TOLERATED_DIR in rel.split("/")


def compare_to_list(
    col: Collection, where: object, listed: dict[str, str], actual: dict[str, str], *, allow_unlisted: bool = False
) -> None:
    """Every listed file present with its recorded hash; unlisted files only under __pycache__ (or when allowed)."""
    col.verified["hash_lists"] += 1
    col.verified["files_hashed"] += len(actual)
    for rel, want in sorted(listed.items()):
        if rel not in actual:
            col.bad(where, f"{rel} is listed but absent")
        elif actual[rel] != want:
            col.bad(where, f"{rel} hashes differently from its list")
    for rel in sorted(set(actual) - set(listed)):
        if _in_pycache(rel) or allow_unlisted:
            col.tolerated.append(f"{where}: unlisted {rel}")
        else:
            col.bad(where, f"{rel} exists but is not in the list")


def read_hash_list(col: Collection, path: Path) -> dict[str, str]:
    try:
        return parse_hash_list(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        col.bad(path, f"unusable hash list ({e})")
        return {}


def read_json_object(col: Collection, path: Path) -> dict[str, object]:
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        col.bad(path, f"unreadable JSON ({type(e).__name__})")
        return {}
    if not isinstance(parsed, dict):
        col.bad(path, "not a JSON object")
        return {}
    return parsed


def parse_config(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2:
            out[parts[0]] = parts[1].strip()
    return out


def check_entries(col: Collection, where: Path, allowed: set[str], required: set[str]) -> None:
    present = {p.name for p in where.iterdir()}
    for name in sorted(present - allowed):
        col.bad(where, f"unexpected entry {name!r}")
    for name in sorted(required - present):
        col.bad(where, f"expected product {name!r} is missing")


def parse_p_endpoints(path: Path) -> str | None:
    """Why a case-P endpoints.json is not a usable endpoint record (missing, empty, not one ENDPOINTS line holding a JSON
    object), or None when it is one. Only the shape is judged here; its content is the analysis's to reconcile."""
    if not path.is_file():
        return "endpoints.json is missing"
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        return "endpoints.json is empty"
    lines = [ln for ln in text.splitlines() if ln.startswith("ENDPOINTS ")]
    try:
        doc = json.loads(lines[0][len("ENDPOINTS ") :]) if len(lines) == 1 else None
    except ValueError:
        doc = None
    return None if isinstance(doc, dict) else "endpoints.json does not hold exactly one ENDPOINTS line with a JSON object"


def read_f_record(path: Path) -> tuple[dict[str, object] | None, str | None]:
    """(record, None) for a parsable case-F record.json, else (None, why it is not one)."""
    if not path.is_file():
        return None, "record.json is missing"
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        return None, f"record.json is unparseable ({type(e).__name__})"
    if not isinstance(parsed, dict):
        return None, "record.json is not a JSON object"
    return parsed, None


def check_seed_dir(
    col: Collection, arm: str, case: str, seed: int, energy: int, d: Path, snap: dict[str, str], binary: str | None,
    fcase: dict[str, str],
) -> tuple[str, str]:
    """Verify one run directory, classify it and queue its copies. Returns (outcome, reason): outcome `collected` or
    `failed` (an integrity problem is recorded in col.problems instead, and refuses the whole collection).

    The workflows' failure shapes (write_run and the endpoint step of both workflow files):
      interrupted  no run.json: the job was cancelled or the host died while the run was in progress;
      transport    run.json transport_status != ok: write_run wrote it and exited the job; out/sha256.txt and the
                   endpoint record are written only after `ok`, so neither can exist;
      endpoint     transport ok (Dose files and out/sha256.txt present), but P endpoints.json missing, empty or not
                   one ENDPOINTS line, or F record.json missing, unparseable or endpoint_status != ok.
    Checks applied to every shape: no unexpected entry; run.json (when present) names this arm, case, seed and energy
    and the snapshot's binary; the config (when present) names this seed, the arm's threads and the requested count.
    Skipped where the shape makes them impossible: Dose hashes when there is no out/sha256.txt (a transport failure
    never writes one; a partial Dose is listed as unverified), the log's primaries line when run.json simulated is
    null (write_run sets it only from that line), and record.json's hashes when there is no parsable record.
    """
    cfg_name = "config.txt" if case == "P" else "cfg.txt"
    out_name = "out" if case == "P" else "out_seed"
    rec_name = "endpoints.json" if case == "P" else "record.json"
    allowed = P_SEED_ENTRIES if case == "P" else F_SEED_ENTRIES
    check_entries(col, d, allowed, set())
    present = {p.name for p in d.iterdir()}
    out = d / out_name
    has_sums = (out / "sha256.txt").is_file()
    run: dict[str, object] = {}
    if "run.json" not in present:
        shape = "interrupted"
        reason = "interrupted: the run directory has no run.json"
        if rec_name in present:
            col.bad(d, f"{rec_name} exists without run.json")
    else:
        run = read_json_object(col, d / "run.json")
        for key, want in (("arm", arm), ("case", case), ("seed", seed), ("energy_mev", energy)):
            if run.get(key) != want:
                col.bad(d / "run.json", f"{key} is {run.get(key)!r}, the directory says {want!r}")
        if binary is not None and run.get("binary_sha256") != binary:
            col.bad(d / "run.json", "binary_sha256 differs from snapshot/binary_sha256.txt")
        status = run.get("transport_status")
        if not isinstance(status, str) or not status:
            col.bad(d / "run.json", f"transport_status {status!r} is not a status the workflows write")
        for name in sorted({"log.txt", cfg_name, out_name} - present):
            col.bad(d, f"expected product {name!r} is missing")
        if status == aa.TRANSPORT_OK:
            shape, reason = "ok", ""
            for name in (*DOSE_FILES, "sha256.txt"):
                if not (out / name).is_file():
                    col.bad(out, f"expected product {name!r} is missing (transport_status is ok)")
        else:
            shape, reason = "transport", f"transport failure: transport_status {status!r}"
            for unwritten in (rec_name, f"{out_name}/sha256.txt"):
                if (d / unwritten).exists():
                    col.bad(d / unwritten, f"exists, but transport_status is {status!r}: the workflows write it only after ok")
    if (d / cfg_name).is_file():
        cfg = parse_config(d / cfg_name)
        try:
            primaries = int(float(cfg.get("Num_Primaries", "")))
        except ValueError:
            primaries = None
        threads = run.get("threads") if run else aa.ARM_THREADS[arm]
        requested = run.get("requested") if run else aa.REQUESTED[case]
        for key, got, want in (("RNG_Seed", cfg.get("RNG_Seed"), str(seed)), ("Num_Threads", cfg.get("Num_Threads"), str(threads))):
            if got != want:
                col.bad(d / cfg_name, f"{key} is {got!r}, expected {want!r}")
        if primaries != requested:
            col.bad(d / cfg_name, f"Num_Primaries {cfg.get('Num_Primaries')!r} is not the requested {requested!r}")
    if (d / "log.txt").is_file() and run and (shape == "ok" or run.get("simulated") is not None):
        hit = PRIMARIES.search((d / "log.txt").read_text(encoding="utf-8", errors="replace"))
        if hit is None or int(hit.group(1)) != run.get("simulated"):
            col.bad(d / "log.txt", f"primaries line {hit.group(0) if hit else None!r} does not give run.json simulated")
    dose: dict[str, str] = {}
    if out.is_dir():
        have = {p.name for p in out.iterdir()}
        if has_sums:
            listed = read_hash_list(col, out / "sha256.txt")
            dose = {n: sha256_file(out / n) for n in DOSE_FILES if (out / n).is_file()}
            compare_to_list(col, out, listed, dose)
            if set(listed) != set(DOSE_FILES):
                col.bad(out / "sha256.txt", f"lists {sorted(listed)}, expected {sorted(DOSE_FILES)}")
            for extra in sorted(have - {*DOSE_FILES, "sha256.txt"}):
                col.tolerated.append(f"{out}: unlisted output {extra}")
        elif have:
            col.tolerated.append(f"{out}: not verified, no sha256.txt ({shape}): {', '.join(sorted(have))}")
    record: dict[str, object] | None = None
    if shape == "ok" and case == "P":
        why = parse_p_endpoints(d / rec_name)
        if why is not None:
            reason = f"endpoint failure: {why}"
    elif shape == "ok":
        record, why = read_f_record(d / rec_name)
        if why is not None:
            reason = f"endpoint failure: {why}"
        elif record is not None and record.get("endpoint_status") != "ok":
            error = record.get("endpoint_error")
            reason = f"endpoint failure: record.json endpoint_status {record.get('endpoint_status')!r}" + (
                f" ({error})" if error else ""
            )
    if case == "F":
        check_f_files(col, d, run, cfg_name, dose, fcase, snap, record)
    col.verified["runs"] += 1
    outcome = "failed" if reason else "collected"
    if outcome == "collected":
        base = f"{arm}/{case}/s{seed}"
        col.copies += [(d / "run.json", f"{base}/run.json"), (d / rec_name, f"{base}/{rec_name}"),
                       (d / cfg_name, f"{base}/{cfg_name}"), (d / "log.txt", f"{base}/log.txt")]
        if has_sums:
            col.copies.append((out / "sha256.txt", f"{base}/out_sha256.txt"))
        return outcome, ""
    base = f"{aa.NOT_ESTABLISHED_DIR}/{arm}/{case}/s{seed}"
    diagnostics = [(d / name, f"{base}/{name}") for name in ("run.json", rec_name, cfg_name, "log.txt") if (d / name).is_file()]
    if has_sums:
        diagnostics.append((out / "sha256.txt", f"{base}/out_sha256.txt"))
    col.not_established.append({"arm": arm, "case": case, "seed": seed, "energy": energy, "outcome": "failed",
                                "reason": reason, "copies": diagnostics})
    return outcome, reason


def check_f_files(
    col: Collection, d: Path, run: dict[str, object], cfg_name: str, dose: dict[str, str], fcase: dict[str, str],
    snap: dict[str, str], record: dict[str, object] | None,
) -> None:
    """Case f: the case files copied into the run directory against the fcase list and the snapshot (every shape: the
    workflow copies them before the run starts), and, when there is a parsable record.json, its hashes against them."""
    rec_name = d / "record.json"
    raw_hashes = record.get("sha256") if record is not None else None
    hashes: dict[str, object] = raw_hashes if isinstance(raw_hashes, dict) else {}
    if record is not None:
        if (d / cfg_name).is_file():
            cfg_hash = sha256_file(d / cfg_name)
            for key, got in (("cfg_sha256", record.get("cfg_sha256")), ("sha256.config", hashes.get("config"))):
                if got != cfg_hash:
                    col.bad(rec_name, f"{key} does not match cfg.txt")
        if hashes.get("binary") != run.get("binary_sha256"):
            col.bad(rec_name, "sha256.binary differs from run.json binary_sha256")
        for name in DOSE_FILES:
            if name in dose and hashes.get(name) != dose[name]:
                col.bad(rec_name, f"sha256.{name} differs from the Dose file")
    for key, name in (("cube.mhd", "cube.mhd"), ("cube.raw", "cube.raw"), ("plan E200_S150.txt", "E200_S150.txt")):
        if (d / name).is_file():
            got = sha256_file(d / name)
            if got != fcase.get(name) or (record is not None and hashes.get(key) != got):
                col.bad(d, f"{name} differs from fcase_sha256.txt or from record.json sha256[{key!r}]")
    for sub, prefix in (("Materials", "Materials/"), ("Scanners/default", "Scanners/default/")):
        if (d / sub).is_dir():
            want = {r: h for r, h in snap.items() if r.startswith(prefix)}
            compare_to_list(col, d / sub, {r[len(prefix) :]: h for r, h in want.items()}, tree_hashes(d / sub))
    bdl = "BDL/BDL_default_DN_RangeShifter.txt"
    if (d / bdl).is_file() and sha256_file(d / bdl) != snap.get(bdl):
        col.bad(d / bdl, "differs from the snapshot")
    for sub, names in (("Scanners", {"default"}), ("BDL", {"BDL_default_DN_RangeShifter.txt"})):
        if (d / sub).is_dir():
            check_entries(col, d / sub, names, names)


def check_snapshot_tar(col: Collection, tar_path: Path, snap: dict[str, str]) -> None:
    """Part A: inputs.tar (a `git archive` of HEAD) must hold exactly the snapshot's files, byte for byte.

    The workflow extracts it into snapshot/ and then adds build.txt and binary_sha256.txt, so those two are the only
    snapshot files the archive may lack. Directory entries and the archive's pax comment header carry no file content.
    """
    col.verified["hash_lists"] += 1
    members: dict[str, str] = {}
    try:
        with tarfile.open(tar_path) as tf:
            for member in tf.getmembers():
                if member.isdir() or member.type == tarfile.XGLTYPE:
                    continue
                if not member.isfile():
                    col.bad(tar_path, f"member {member.name!r} is not a regular file")
                    continue
                fh = tf.extractfile(member)
                if fh is None:
                    col.bad(tar_path, f"member {member.name!r} has no content")
                    continue
                members[normalise(member.name)] = hashlib.sha256(fh.read()).hexdigest()
    except (OSError, tarfile.TarError) as e:
        col.bad(tar_path, f"unreadable archive ({type(e).__name__})")
        return
    col.verified["files_hashed"] += len(members)
    expected = {rel: h for rel, h in snap.items() if rel not in SNAPSHOT_ADDED_FILES and not _in_pycache(rel)}
    for rel in sorted(set(expected) | set(members)):
        if rel not in members:
            col.bad(tar_path, f"{rel} is in the snapshot but not in the archive")
        elif rel not in expected:
            col.bad(tar_path, f"{rel} is in the archive but not in the snapshot")
        elif members[rel] != expected[rel]:
            col.bad(tar_path, f"{rel} differs between the archive and the snapshot")


def check_case_dir(col: Collection, arm: str, case: str, case_dir: Path) -> None:
    """One native <arm>/<case> directory: its provenance files, snapshot, case inputs and every run."""
    expected = col.frozen.seeds[(arm, case)]
    energies = sorted(set(expected.values())) if case == "P" else []
    # The snapshot, the run root and (case F) the case are written before the seed loop starts, so they are required. A
    # case-P energy directory is created when the loop reaches its first seed: a job that stopped earlier leaves none,
    # so e<E> is optional, but e<E> and its inputs list exist together or not at all.
    required = {"run_root.txt", "snapshot", "snapshot_sha256.txt"}
    if case == "F":
        required |= {"fcase", "fcase_sha256.txt"}
    if PART_OF_ARM[arm] == "A":
        required |= {WINDOWS_SNAPSHOT_TAR}  # the Windows workflow builds its snapshot by extracting this git archive
    allowed = required | {f"e{e}" for e in energies} | {f"e{e}_inputs_sha256.txt" for e in energies}
    allowed |= {f"s{s}" for s in expected} if case == "F" else set()
    check_entries(col, case_dir, allowed, required)
    provenance: list[tuple[Path, str]] = []
    run_root = case_dir / "run_root.txt"
    m = RUN_ROOT.match(run_root.read_text(encoding="utf-8").strip()) if run_root.is_file() else None
    if m is None:
        col.bad(run_root, "missing or not '<arm> case <C> mode <m>, commit <sha>, run <id>'")
    else:
        for key, got, want in (("arm", m["arm"], arm), ("case", m["case"], case), ("mode", m["mode"], "full"),
                               ("commit", m["commit"], col.commit)):
            if got != want:
                col.bad(run_root, f"{key} is {got!r}, expected {want!r}")
    snap: dict[str, str] = {}
    if (case_dir / "snapshot").is_dir() and (case_dir / "snapshot_sha256.txt").is_file():
        snap = read_hash_list(col, case_dir / "snapshot_sha256.txt")
        compare_to_list(col, case_dir / "snapshot", snap, tree_hashes(case_dir / "snapshot"))
    if (case_dir / WINDOWS_SNAPSHOT_TAR).is_file():
        check_snapshot_tar(col, case_dir / WINDOWS_SNAPSHOT_TAR, snap)
    binary: str | None = None
    binary_file = case_dir / "snapshot" / "binary_sha256.txt"
    words = binary_file.read_text(encoding="utf-8").split() if binary_file.is_file() else []
    if words and aa.SHA256_HEX.match(words[0]):
        binary = words[0]
    else:
        col.bad(binary_file, "missing, or does not start with a sha256")
    fcase: dict[str, str] = {}
    if case == "F" and (case_dir / "fcase").is_dir() and (case_dir / "fcase_sha256.txt").is_file():
        fcase = read_hash_list(col, case_dir / "fcase_sha256.txt")
        compare_to_list(col, case_dir / "fcase", fcase, tree_hashes(case_dir / "fcase"), allow_unlisted=True)
        if set(fcase) != set(F_CASE_FILES):
            col.bad(case_dir / "fcase_sha256.txt", f"lists {sorted(fcase)}, expected {sorted(F_CASE_FILES)}")
    seen: dict[int, tuple[str, str]] = {}
    if case == "P":
        for energy in energies:
            edir, inputs = case_dir / f"e{energy}", case_dir / f"e{energy}_inputs_sha256.txt"
            if edir.is_dir() != inputs.is_file():
                col.bad(case_dir, f"e{energy} and e{energy}_inputs_sha256.txt must exist together (the workflow writes both)")
            if edir.is_dir():
                seen.update(check_energy_dir(col, arm, case, energy, edir, case_dir, snap, binary))
    else:
        for entry in sorted(case_dir.glob("s*")):
            sm = SEED_NAME.match(entry.name)
            if sm and entry.is_dir() and int(sm.group(1)) in expected:
                seed = int(sm.group(1))
                seen[seed] = check_seed_dir(col, arm, case, seed, expected[seed], entry, snap, binary, fcase)
    classify_absent(col, arm, case, case_dir, seen)
    for name in ("run_root.txt", "snapshot_sha256.txt"):
        provenance.append((case_dir / name, f"{arm}/{case}/{aa.PROVENANCE_DIR}/{name}"))
    for name in ("build.txt", "binary_sha256.txt"):
        provenance.append((case_dir / "snapshot" / name, f"{arm}/{case}/{aa.PROVENANCE_DIR}/snapshot_{name}"))
    for name in sorted(p.name for p in case_dir.glob("*_sha256.txt") if p.name != "snapshot_sha256.txt"):
        provenance.append((case_dir / name, f"{arm}/{case}/{aa.PROVENANCE_DIR}/{name}"))
    col.copies += [(src, dest) for src, dest in provenance if src.is_file()]


def check_energy_dir(
    col: Collection, arm: str, case: str, energy: int, edir: Path, case_dir: Path, snap: dict[str, str],
    binary: str | None,
) -> dict[int, tuple[str, str]]:
    """Case p, one e<E> directory: its inputs list, and the seed directories (each must belong to this energy).
    Returns {seed: (outcome, reason)} for the seed directories found."""
    expected = col.frozen.seeds[(arm, case)]
    seeds = {p.name for p in edir.iterdir() if p.is_dir() and SEED_NAME.match(p.name)}
    inputs = case_dir / f"e{energy}_inputs_sha256.txt"
    if inputs.is_file():  # a missing list is reported by check_case_dir
        compare_to_list(col, edir, read_hash_list(col, inputs), tree_hashes(edir, skip_top=frozenset(seeds)))
    for sub, prefix in (("Materials", "Materials/"), ("Scanners/Water_Phantom", "Scanners/Water_Phantom/")):
        if (edir / sub).is_dir():
            want = {r[len(prefix) :]: h for r, h in snap.items() if r.startswith(prefix)}
            compare_to_list(col, edir / sub, want, tree_hashes(edir / sub))
    found: dict[int, tuple[str, str]] = {}
    for name in sorted(seeds):
        seed = int(name[1:])
        if seed not in expected:
            col.bad(edir, f"seed directory {name} is not in the frozen list")
        elif expected[seed] != energy:
            col.bad(edir, f"seed {seed} is frozen at {expected[seed]} MeV but sits under e{energy}")
        else:
            found[seed] = check_seed_dir(col, arm, case, seed, energy, edir / name, snap, binary, {})
    for entry in sorted(edir.iterdir()):
        if entry.is_dir() and entry.name.startswith("s") and entry.name not in seeds:
            col.bad(edir, f"directory {entry.name} is not s<six digits>")
    return found


def classify_absent(col: Collection, arm: str, case: str, case_dir: Path, seen: dict[int, tuple[str, str]]) -> None:
    """Seeds with no directory. The job runs its frozen list in order (the list's order, which is the dict's insertion
    order) and exits on a transport or endpoint-step failure, so the seeds it never reached form a SUFFIX of that order:
    those are `absent`. A seed with no directory BEFORE a seed that has one cannot come from a stopped job and refuses."""
    order = list(col.frozen.seeds[(arm, case)])
    reached = [i for i, s in enumerate(order) if s in seen]
    last = reached[-1] if reached else -1
    for seed in order[: last + 1]:
        if seed not in seen:
            col.bad(case_dir, f"expected seed {seed} has no directory, but a later seed of the job's loop has one")
    if last == -1:
        why = "no directory: the job stopped before its first run"
    else:
        outcome, reason = seen[order[last]]
        why = f"no directory: the job's last run directory is s{order[last]} ({outcome}" + (f": {reason})" if reason else ")")
    for seed in order[last + 1 :]:
        col.not_established.append({"arm": arm, "case": case, "seed": seed, "energy": col.frozen.seeds[(arm, case)][seed],
                                    "outcome": "absent", "reason": why, "copies": []})


def survey(col: Collection, sources: list[Path], parts: tuple[str, ...]) -> None:
    """Walk every source root: names, arms, cases; then the checks of each case directory."""
    wanted = {arm for part in parts for arm in aa.PART_ARMS[part]}
    found: dict[str, Path] = {}
    for src in sources:
        if not src.is_dir():
            col.bad(src, "is not a directory")
            continue
        if src.name != col.commit[:12]:
            col.bad(src, f"is named {src.name!r}, the acquisition commit is {col.commit[:12]!r}")
        for entry in sorted(src.iterdir()):
            arm = ARM_DIRS.get(entry.name)
            if arm is None or not entry.is_dir():
                col.bad(src, f"unexpected entry {entry.name!r} (expected one of {sorted(ARM_DIRS)})")
            elif arm not in wanted:
                col.bad(src, f"arm directory {entry.name!r} ({arm}) is not in parts {','.join(parts)}")
            elif arm in found:
                col.bad(src, f"arm {arm} is also present in {found[arm]}")
            else:
                found[arm] = src
    for arm in sorted(wanted - set(found)):
        col.bad("sources", f"arm {arm} is not present in any source")
    for arm, src in sorted(found.items()):
        arm_dir = src / next(k for k, v in ARM_DIRS.items() if v == arm)
        check_entries(col, arm_dir, set(CASE_DIRS), set(CASE_DIRS))
        for low, case in CASE_DIRS.items():
            if (arm_dir / low).is_dir():
                check_case_dir(col, arm, case, arm_dir / low)


def write_copies(copies: list[tuple[Path, str]], out: Path) -> list[dict[str, str]]:
    """Copy every queued file to out (never overwriting), verifying each copy; return the manifest's file entries."""
    entries = []
    for src, rel in sorted(copies, key=lambda c: c[1]):
        dest = out / rel
        if dest.exists():
            msg = f"{dest} already exists: refusing to overwrite"
            raise CollectionError(msg)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        want = sha256_file(src)
        if sha256_file(dest) != want:
            msg = f"copy of {src} to {dest} does not hash the same"
            raise CollectionError(msg)
        entries.append({"path": rel, "source": str(src), "sha256": want})
    return entries


def collect(sources: list[Path], out: Path, parts: tuple[str, ...] = ("A", "B"), commit: str = aa.ACQUISITION_COMMIT) -> dict[str, object]:
    """Verify, then copy. Returns the manifest (also written to <out>/collection_manifest.json)."""
    try:
        frozen = aa.expected_population(commit)
    except aa.InputError as e:
        raise CollectionError(str(e)) from e
    if out.exists():
        msg = f"{out} exists: the output directory must not exist"
        raise CollectionError(msg)
    col = Collection(commit, frozen)
    survey(col, sources, parts)
    if col.problems:
        shown = col.problems[:25]
        more = f" (+{len(col.problems) - 25} more)" if len(col.problems) > 25 else ""
        msg = f"refusing, {len(col.problems)} problem(s):\n  " + "\n  ".join(shown) + more
        raise CollectionError(msg)
    entries = write_copies(col.copies, out)
    not_established = []
    for ne in sorted(col.not_established, key=lambda e: (e["arm"], e["case"], e["seed"])):
        files = write_copies(ne["copies"], out)  # type: ignore[arg-type]
        not_established.append({**{k: v for k, v in ne.items() if k != "copies"}, "files": files})
    counts = {k: sum(1 for ne in not_established if ne["outcome"] == k) for k in ("failed", "absent")}
    manifest: dict[str, object] = {
        "schema": 2,
        "acquisition_commit": commit,
        "parts": list(parts),
        "sources": [str(s) for s in sources],
        "files": entries,
        "verified": col.verified,
        "tolerated": sorted(col.tolerated),
        "not_established": not_established,
        "outcomes": {"collected": col.verified["runs"] - counts["failed"], **counts},
    }
    (out / aa.MANIFEST_NAME).write_text(json.dumps(manifest, indent=1, sort_keys=True, allow_nan=False) + "\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, action="append", required=True, help="a native <sha12> root (repeatable)")
    ap.add_argument("--out", type=Path, required=True, help="the analysis input directory to create (must not exist)")
    ap.add_argument("--parts", default="A,B", help="comma-separated parts to collect (default A,B)")
    ap.add_argument("--commit", default=aa.ACQUISITION_COMMIT, help="full acquisition commit sha")
    args = ap.parse_args(argv)
    parts = tuple(p.strip() for p in args.parts.split(","))
    if not parts or any(p not in aa.PART_ARMS for p in parts) or len(set(parts)) != len(parts):
        print(f"apples_collect: --parts must be distinct values from {sorted(aa.PART_ARMS)}", file=sys.stderr)
        return 2
    try:
        manifest = collect(args.source, args.out, parts, args.commit)
    except (CollectionError, OSError) as e:
        print(f"apples_collect: {e}", file=sys.stderr)
        return 2
    outcomes: dict[str, int] = manifest["outcomes"]  # type: ignore[assignment]
    print(f"collected {outcomes['collected']}, failed {outcomes['failed']}, absent {outcomes['absent']}")
    print(f"{len(manifest['files'])} file(s) from {len(args.source)} source(s) into {args.out}")  # type: ignore[arg-type]
    return 0


if __name__ == "__main__":
    sys.exit(main())
