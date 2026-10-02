"""Pre-specified analysis for the same-host comparison programme, parts A and B (docs/apples_to_apples_plan.md,
"Pre-specified analysis" and "Endpoints").

Usage:
    uv run --no-project --with numpy --with scipy python validation/apples_analyse.py <root> [--json PATH] [--parts A,B]

Input layout, <root>/<arm>/<case>/s<seed>/ (built from the workflows' native trees by apples_collect.py) with arms A-up, A-port (part A) and B-up, B-pgcc, B-picc (part B), cases P
and F. Each seed directory holds run.json plus endpoints.json (case P: one line "ENDPOINTS " + JSON, as
pencil_endpoints.py prints it) or record.json (case F: as written by platform_study_record.py; the 13 endpoints are
its "metrics" object).

Contrasts, each arm minus reference (ratio for rings): A-port vs A-up; B-pgcc vs B-up and B-picc vs B-up (each its own
family); B-pgcc vs B-picc is DESCRIPTIVE only (estimates and intervals, no margin, no equivalence claim).

Per endpoint: estimate, Welch SE, Welch-Satterthwaite df, 90% and 95% intervals (same estimator conventions as
validation/report_tables.py: Welch on run-level values, log scale for ring fractions, difference of mean logs).
TOST at one-sided alpha 0.05 is the 90% interval strictly inside the margin; p_TOST is the larger of the two one-sided
t p-values. Outcomes: equivalent / not_equivalent / inconclusive / not_established. An endpoint is not_established if
any run of either arm is unusable for it (missing record, non-finite value, zero or negative ring fraction, failed
sigma fit, R80 with multiple crossings, transport_status not "ok") or if either arm has fewer than 2 runs. Nothing is
excluded and nothing is imputed.

Claims, kept separate: (1) the JOINT claim of a confirmatory contrast holds only if every endpoint is equivalent at
unadjusted alpha (intersection-union); (2) INDIVIDUAL claims use Holm on p_TOST over the whole family (39 P + 13 F = 52),
with not_established endpoints kept in the family as p = 1 (never rejected). Holm is never applied to the joint claim.

Frozen population and identity (amended 2026-10-01 after review 6933, before any full-run endpoint was read by the
analysis). The expected runs are the `full` seeds of the workflow files at the ACQUISITION COMMIT (default
2f9dab404cea02f5072352f7ae5773dca27eb2a1, read with `git show <commit>:<path>` and parsed by apples_seeds_check.py), per
arm, case and energy. A seed directory is accepted only if it is exactly `s<seed>` for an expected seed of that arm and
case. Unexpected seeds, directories whose numeric identity repeats (`s960050` and `s0960050`; BOTH are rejected), and
non-canonical names are rejected and listed, never read; expected seeds with no accepted directory are listed as
missing. A run whose run.json/record does not carry the full-run provenance below is kept in the population as
unusable and recorded as a provenance failure. ANY of these makes the contrast PARTIAL, and a PARTIAL contrast emits NO
confirmatory claim: no joint claim, no Holm decision and no TOST classification, only labelled descriptive estimates.

Provenance required of every run (else the run is unusable and the contrast PARTIAL): run.json has arm, case and seed
equal to the directory (a contradiction is fatal), mode "full", commit equal to the acquisition commit, requested 1e7
(P) or 3e7 (F), simulated >= requested, threads 4 (A arms) or 3 (B arms), energy_mev equal to the energy of that seed in
the workflow (F: 200), transport_status "ok" and a sha256 binary_sha256 that is the same across the whole arm. The
endpoint record is reconciled with run.json: case P needs layout "mcsquare", the label
`<arm>_P_E<energy>_N1e7_seed<seed>` and slab_depths_mm exactly the energy's frozen depths (100: 40,60; 150: 80,125;
200: 100,200); case F is read through the RECORD BINDING of the acquisition commit (F_RECORD_BINDINGS), which names
the producer version's record schema and study. 2f9dab40 (the full runs' snapshot) is bound to the schema its
platform_study_record.py writes: exactly study, platform, seed, host, machine, commit, compiler, sha256 (ten keys),
materials and metrics, with NO cfg_sha256, metrics_sha256 or endpoint_status (that writer has no failure record: an
endpoint failure leaves no record.json). Such a record needs study "pe1", platform = arm, seed, commit and sha256.binary
= run.json binary_sha256; its endpoints are its `metrics`, which no hash covers (stated in the report). A record of the
later schema (cfg_sha256, metrics_sha256, endpoint_status, from 0f5ef7c) is read only under a commit bound to that schema
and also needs endpoint_status "ok", metrics_sha256 = digest of metrics and cfg_sha256 = sha256.config. A record that
mixes the schemas, has the other schema than its commit's binding, lacks or adds a field, or whose commit has no binding
REFUSES (exit 2); a wrong study, platform, seed, commit or binary keeps the run as unusable (PARTIAL).
Domain bounds on case P: sigma_d in [1, 20] mm (the fit contract in pencil_endpoints.py) and ring fractions in (0, 1]
(0 is inside the energy-fraction domain but has no logarithm); otherwise that endpoint is not_established.

Every file read (run.json and the endpoint record of each accepted seed directory, collection_manifest.json when
present and the not_established entries' copied files) enters the DATASET FINGERPRINT, the sha256 of the sorted "<sha256>  <path relative to the root>" lines, written
to the report and the JSON.

Collection attestation. Confirmatory claims need a VERIFIED collection. collection_manifest.json (apples_collect.py,
schema 3) is validated BEFORE any run is read: a schema other than 3, an acquisition_commit other than the frozen
commit in use, parts that do not cover the analysed parts, or a case-F record binding other than this analysis's refuse;
then every file it lists must hash as recorded and every file read must be listed. WITHOUT a manifest the tree may
still be analysed as a preview, but every contrast is PARTIAL through the same withholding path as a population issue:
descriptive estimates only, no joint claim, no Holm decision, no TOST classification, stated in the report header and
in the JSON (confirmatory_claims_possible false).

Runs that were not established (apples_collect.py, manifest schema 2): the collector lists every frozen run it did not
collect in the manifest's `not_established` list, outcome `failed` (interrupted, transport failure, endpoint failure;
its diagnostic files copied under .not_established/<arm>/<case>/s<seed>/, which is never read as a run) or `absent` (the
job stopped before reaching it), with a reason. Each entry is re-verified (its copied files hash as recorded and enter
the dataset fingerprint; it is a frozen arm, case, seed at its frozen energy; no seed is listed twice; no listed seed is
also a collected run) and then makes its contrasts PARTIAL by the rule above: it is a population issue, so descriptive
estimates only, no joint claim, no Holm decision, no TOST classification. With a manifest present every frozen seed
must be accounted for exactly once, as a collected run or a not_established entry; a seed that is neither refuses.
Without a manifest a missing seed is a population issue as before. The report and the JSON list each not-established
run (arm, case, seed, energy, outcome, reason).

Design completeness: the plan has 8 runs per arm, case and energy (case F: 8 at 200 MeV); the frozen population
carries exactly that.

Exit status: 0 when the analysis ran, whatever the outcomes; 2 when the inputs are unusable (reason on stderr).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import statistics as st
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from scipy.stats import t as student_t

sys.path.insert(0, str(Path(__file__).resolve().parent))

import apples_seeds_check as sc
from pencil_endpoints import SIGMA_VALID_MM
from platform_study_analyse import CAX as F_CAX
from platform_study_analyse import ENDPOINTS as F_ENDPOINTS
from platform_study_analyse import margin as f_margin

ALPHA = 0.05  # one-sided TOST alpha, and the Holm alpha
FAMILY_SIZE = 52  # 39 (P) + 13 (F)

PART_ARMS = {"A": ("A-up", "A-port"), "B": ("B-up", "B-pgcc", "B-picc")}
CASES = ("P", "F")
# (energy MeV) -> slab depths (mm), fixed in the plan
SLAB_DEPTHS = {100: (40, 60), 150: (80, 125), 200: (100, 200)}
RINGS = ((5, 10), (10, 20), (20, 40), (40, 80), (80, 200))
R80_MARGIN = 0.05
SIGMA_MARGIN = 0.02
RING_RATIO = (0.98, 1.02)
WIDE_RING_RATIO = (0.95, 1.05)  # ring 80-200 only
TRANSPORT_OK = "ok"
RUNS_PER_CELL = 8  # per arm x case x energy (case F: one cell, at 200 MeV)

# Domain bounds on case-P endpoint values. sigma: the producer's fit contract (pencil_endpoints.SIGMA_VALID_MM, a fit
# outside it is recorded as failed). Ring fractions are energy fractions: 0 < f <= 1 (0 is a valid fraction but has no
# logarithm, so it is unusable on the analysis scale).
SIGMA_BOUNDS_MM = SIGMA_VALID_MM
RING_MAX = 1.0

# The acquisition: the full runs were requested at this commit; the frozen run lists are the workflow files there.
ACQUISITION_COMMIT = "2f9dab404cea02f5072352f7ae5773dca27eb2a1"
REPO = Path(__file__).resolve().parents[1]
WORKFLOW_PATHS = (".gitea/workflows/apples-a-windows.yml", ".gitea/workflows/apples-b-linux.yml")
FULL_SHA = re.compile(r"^[0-9a-f]{40}$")
SHA256_HEX = re.compile(r"^[0-9a-f]{64}$")
ARM_THREADS = {"A-up": 4, "A-port": 4, "B-up": 3, "B-pgcc": 3, "B-picc": 3}
REQUESTED = {"P": 10_000_000, "F": 30_000_000}
P_LABEL = "{arm}_P_E{energy}_N1e7_seed{seed}"  # the --label the workflows pass to pencil_endpoints.py
MANIFEST_NAME = "collection_manifest.json"
MANIFEST_SCHEMA = 3  # written by apples_collect.py; the only one this analysis reads
SUPPORTED_MANIFEST_SCHEMAS = (MANIFEST_SCHEMA,)
NO_MANIFEST_REASON = (
    "collection not verified: there is no collection manifest, so the collector's checks (config, snapshot, log, Dose "
    "and record-hash verification) are not attested for this tree; every contrast is descriptive only"
)
PROVENANCE_DIR = ".provenance"  # apples_collect.py puts run_root.txt and the *_sha256.txt files here; skipped (hidden)
NOT_ESTABLISHED_DIR = ".not_established"  # apples_collect.py copies a failed run's diagnostics here, never as a run
NE_OUTCOMES = ("failed", "absent")

IDENTITY_SOURCE = (
    "Run identity comes from the frozen workflow run lists at the acquisition commit and is validated against "
    "run.json and the endpoint record (arm, case, seed, mode, commit, requested/simulated, threads, energy, binary "
    "hash, slab depths, record hashes). A case-F record.json is read only through the record binding of the "
    "acquisition commit, which names the schema its writer emitted and the study it wrote (pe1 at 2f9dab40); a "
    "record of another schema, or mixing two, refuses; a record whose study is not the bound one is an unusable run."
)

# ------------------------------------------------------------------------- case-F record schemas, by producer version

# The fields platform_study_record.py writes, per version of the writer. At the acquisition commit 2f9dab40 the writer
# (`git show 2f9dab40:validation/platform_study_record.py`) emits exactly LEGACY_RECORD_FIELDS, computes the endpoints
# BEFORE opening record.json and has no failure path: an endpoint failure raises, no record.json is written and the job
# stops. From 0f5ef7c the writer adds V2_ENDPOINT_FIELDS (and endpoint_error on failure) and writes a record either way.
F_RECORD_SHA256_KEYS = frozenset(
    {"Dose.raw", "Dose.mhd", "config", "plan E200_S150.txt", "cube.mhd", "cube.raw", "BDL", "HU_Density", "HU_Material",
     "binary"}
)
LEGACY_RECORD_FIELDS = frozenset(
    {"study", "platform", "seed", "host", "machine", "commit", "compiler", "sha256", "materials", "metrics"}
)
V2_ENDPOINT_FIELDS = frozenset({"cfg_sha256", "metrics_sha256", "endpoint_status"})
MATERIALS_KEYS = frozenset({"files", "combined_sha256"})

_RECORD_VERIFIED_COMMON = (
    "study equals the binding's study; platform, seed and commit equal the run directory, run.json and the binding",
    "sha256.config equals the sha256 of the run's cfg.txt",
    "sha256.binary equals run.json binary_sha256 and snapshot/binary_sha256.txt",
    "sha256.Dose.raw and sha256.Dose.mhd equal the out_seed files, which hash as out_seed/sha256.txt lists",
    "sha256 of the plan and cube equal the run directory's files and fcase_sha256.txt",
    "sha256.BDL, HU_Density and HU_Material equal the run directory's files and the snapshot's",
    "materials recomputed from the run directory's Materials/ with the writer's tree_digest formula",
)


@dataclass(frozen=True)
class RecordSchema:
    """The top-level fields one version of platform_study_record.py writes, and what can be verified from them."""

    name: str
    fields: frozenset[str]  # required, exactly
    optional: frozenset[str]
    endpoint_status: bool  # the writer records endpoint failures (else a failure leaves no record.json)
    verified: tuple[str, ...]  # what apples_collect.py (files) and this analysis (run.json) verify
    not_verifiable: tuple[str, ...]


RECORD_SCHEMA_2F9DAB40 = RecordSchema(
    "platform_study_record@2f9dab40",
    LEGACY_RECORD_FIELDS,
    frozenset(),
    endpoint_status=False,
    verified=_RECORD_VERIFIED_COMMON,
    not_verifiable=(
        "metrics: this schema carries no metrics hash; the endpoint values are taken from `metrics` as written",
        "a missing record.json after transport ok is the endpoint-failure shape (the writer has no failure record)",
    ),
)
RECORD_SCHEMA_V2 = RecordSchema(
    "platform_study_record@0f5ef7c",
    LEGACY_RECORD_FIELDS | V2_ENDPOINT_FIELDS,
    frozenset({"endpoint_error"}),
    endpoint_status=True,
    verified=(
        *_RECORD_VERIFIED_COMMON,
        "cfg_sha256 equals sha256.config",
        "metrics_sha256 equals the digest of metrics",
        "endpoint_status is ok (else the run is an endpoint failure)",
    ),
    not_verifiable=(),
)


@dataclass(frozen=True)
class RecordBinding:
    """An explicit, named binding of a producer version: case-F records under `commit` must have exactly `schema` and
    carry `study`. Nothing else is relaxed."""

    name: str
    commit: str
    study: str
    schema: RecordSchema
    reason: str


F_RECORD_BINDINGS = (
    RecordBinding(
        "pe1-at-acquisition-2f9dab40",
        ACQUISITION_COMMIT,
        "pe1",
        RECORD_SCHEMA_2F9DAB40,
        "the full runs use the snapshot of 2f9dab40, whose platform_study_record.py writes the pe1 study id and none "
        "of cfg_sha256, metrics_sha256 or endpoint_status",
    ),
)


def record_binding(commit: str) -> RecordBinding | None:
    """The case-F record binding for an acquisition commit, or None (no record from that commit can be read)."""
    return next((b for b in F_RECORD_BINDINGS if b.commit == commit), None)

# (name, arm, reference, confirmatory?)
CONTRASTS_BY_PART = {
    "A": (("A-port vs A-up", "A-port", "A-up", True),),
    "B": (
        ("B-pgcc vs B-up", "B-pgcc", "B-up", True),
        ("B-picc vs B-up", "B-picc", "B-up", True),
        ("B-pgcc vs B-picc (descriptive)", "B-pgcc", "B-picc", False),
    ),
}

_SEED_DIR = re.compile(r"^s(\d+)$")


class InputError(Exception):
    """The inputs are unusable (missing arm, unexpected layout); the analysis cannot run."""


@dataclass(frozen=True)
class Spec:
    """One endpoint: where its value comes from, its scale and its margin on the analysis scale."""

    eid: str  # e.g. "P100/sigma_40", "F/lateral_127_5"
    case: str
    key: str  # key in the run's endpoint record
    log: bool
    lo: float  # margin on the analysis scale (log of the ratio when log)
    hi: float
    energy: int | None = None


@dataclass
class Run:
    """One seed directory. `usable` False means no endpoint of this run is established."""

    label: str
    seed: int
    energy: int | None
    usable: bool
    problems: list[str] = field(default_factory=list)
    provenance: list[str] = field(default_factory=list)  # the subset of problems that are identity/provenance failures
    mode: str | None = None  # run.json "mode", None if absent
    commit: str | None = None  # run.json "commit", None if absent
    binary: str | None = None  # run.json "binary_sha256" when it is a well-formed sha256
    host: str | None = None
    binding: str | None = None  # name of the RecordBinding that admitted this run's case-F record, if any
    raw: dict[str, float | None] = field(default_factory=dict)  # eid -> value, None = invalid for that endpoint


@dataclass
class Result:
    """One endpoint of one contrast."""

    eid: str
    spec: Spec
    outcome: str
    reason: str = ""
    estimate: float | None = None  # analysis scale
    se: float | None = None
    df: float | None = None
    ci90: tuple[float, float] | None = None
    ci95: tuple[float, float] | None = None
    p_tost: float | None = None
    holm_p: float | None = None
    holm_decision: str = ""


@dataclass
class Frozen:
    """The frozen population: (arm, case) -> {seed: energy}, from the full run lists at `commit`."""

    commit: str
    seeds: dict[tuple[str, str], dict[int, int]]


# ----------------------------------------------------------------------------------------------- specification


def p_specs() -> list[Spec]:
    """The 39 case-P endpoints (13 per energy)."""
    out: list[Spec] = []
    for energy, depths in SLAB_DEPTHS.items():
        out.append(Spec(f"P{energy}/R80", "P", "R80", False, -R80_MARGIN, R80_MARGIN, energy))
        for d in depths:
            out.append(Spec(f"P{energy}/sigma_{d}", "P", f"sigma_{d}", False, -SIGMA_MARGIN, SIGMA_MARGIN, energy))
        for d in depths:
            for lo, hi in RINGS:
                ratio = WIDE_RING_RATIO if (lo, hi) == (80, 200) else RING_RATIO
                key = f"ring_{d}_{lo}_{hi}"
                out.append(Spec(f"P{energy}/{key}", "P", key, True, math.log(ratio[0]), math.log(ratio[1]), energy))
    return out


def f_specs() -> list[Spec]:
    """The 13 case-F endpoints with the #31 margins, taken from platform_study_analyse.margin (not copied)."""
    out = []
    for ep in F_ENDPOINTS:
        lo, hi = f_margin(ep)
        out.append(Spec(f"F/{ep}", "F", ep, ep in F_CAX, lo, hi))
    return out


def all_specs() -> list[Spec]:
    specs = p_specs() + f_specs()
    if len(specs) != FAMILY_SIZE or len({s.eid for s in specs}) != FAMILY_SIZE:
        msg = f"endpoint family has {len(specs)} entries, expected {FAMILY_SIZE} distinct"
        raise AssertionError(msg)
    return specs


# ------------------------------------------------------------------------------------------------------ loading


def _num(v: object) -> float | None:
    """A finite real number, else None (booleans are not numbers here)."""
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return None
    return float(v) if math.isfinite(v) else None


def _int(v: object) -> int | None:
    """An integer (not a bool), else None. A float is not an integer here: 1e7 must be written 10000000."""
    return v if isinstance(v, int) and not isinstance(v, bool) else None


class Ledger:
    """Every input file the analysis reads, with the sha256 of the bytes it read (the dataset fingerprint)."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.files: dict[str, str] = {}

    def read_bytes(self, path: Path) -> bytes:
        data = path.read_bytes()
        self.files[path.relative_to(self.root).as_posix()] = hashlib.sha256(data).hexdigest()
        return data

    def read_text(self, path: Path) -> str:
        return self.read_bytes(path).decode("utf-8")

    def fingerprint(self) -> dict[str, object]:
        lines = [f"{h}  {p}\n" for p, h in sorted(self.files.items())]
        return {
            "sha256": hashlib.sha256("".join(lines).encode()).hexdigest(),
            "n_files": len(lines),
            "definition": "sha256 of the sorted '<sha256>  <path relative to root>' lines of every file read",
        }


def _read_json(path: Path, ledger: Ledger) -> object:
    return json.loads(ledger.read_text(path))


def expected_population(
    commit: str = ACQUISITION_COMMIT, repo: Path = REPO, sources: dict[str, str] | None = None
) -> Frozen:
    """The frozen full run lists. `sources` (workflow name -> text) replaces `git show <commit>:<path>` (tests only).

    Refuses an abbreviated or malformed commit, an unreadable workflow, any apples_seeds_check violation, and a list
    that is not 5 arms x 2 cases with the planned counts per energy.
    """
    if not FULL_SHA.match(commit):
        msg = f"acquisition commit must be a full 40-hex sha, got {commit!r}"
        raise InputError(msg)
    if sources is None:
        sources = {}
        for path in WORKFLOW_PATHS:
            proc = subprocess.run(
                ["git", "-C", str(repo), "show", f"{commit}:{path}"], capture_output=True, text=True, check=False
            )
            if proc.returncode != 0:
                msg = f"cannot read {path} at {commit}: {proc.stderr.strip()}"
                raise InputError(msg)
            sources[Path(path).name] = proc.stdout
    errors, runs, *_ = sc._check(sources, sc.ALL_ARMS)
    if errors:
        msg = f"the workflow run lists at {commit} fail apples_seeds_check: " + "; ".join(errors[:5])
        raise InputError(msg)
    seeds: dict[tuple[str, str], dict[int, int]] = {}
    for arm, case, mode, energy, seed, _where in runs:
        if mode == "full":
            seeds.setdefault((arm, case), {})[seed] = energy
    for arm in sc.ALL_ARMS:
        for case in CASES:
            per_energy = sorted(
                {e: sum(1 for v in seeds.get((arm, case), {}).values() if v == e) for e in set(seeds.get((arm, case), {}).values())}.items()
            )
            want = [(200, RUNS_PER_CELL)] if case == "F" else [(e, RUNS_PER_CELL) for e in sorted(SLAB_DEPTHS)]
            if per_energy != want:
                msg = f"{arm} {case}: frozen full list is {per_energy}, expected {want}"
                raise InputError(msg)
    return Frozen(commit, seeds)


def _p_values(rec: dict[str, object], energy: int) -> dict[str, float | None]:
    """Case-P endpoint values of one run for one energy; None where that endpoint is invalid for this run.

    Domain: sigma in SIGMA_BOUNDS_MM inclusive; ring fractions in (0, 1].
    """
    out: dict[str, float | None] = {}
    r80 = _num(rec.get("R80"))
    # fail closed: the multiple-crossing flag must be present and exactly False
    out[f"P{energy}/R80"] = r80 if rec.get("R80_multiple_crossings") is False else None
    for d in SLAB_DEPTHS[energy]:
        sigma = _num(rec.get(f"sigma_{d}"))
        in_bounds = sigma is not None and SIGMA_BOUNDS_MM[0] <= sigma <= SIGMA_BOUNDS_MM[1]
        out[f"P{energy}/sigma_{d}"] = sigma if in_bounds else None
        for lo, hi in RINGS:
            v = _num(rec.get(f"ring_{d}_{lo}_{hi}"))
            out[f"P{energy}/ring_{d}_{lo}_{hi}"] = v if v is not None and 0 < v <= RING_MAX else None
    return out


def _f_values(metrics: dict[str, object]) -> dict[str, float | None]:
    """Case-F endpoint values of one run; None where invalid (missing, non-finite, bad dose, ambiguous R80)."""
    out: dict[str, float | None] = {}
    for ep in F_ENDPOINTS:
        v = _num(metrics.get(ep))
        if ep in F_CAX and (metrics.get(f"{ep}_invalid") is True or v is None or v <= 0):
            v = None
        if ep in ("r80_mm", "r20_mm"):  # #31 flags BOTH as ambiguous on more than one distal crossing
            crossings = metrics.get(ep.replace("_mm", "_crossings"))
            if isinstance(crossings, bool) or not isinstance(crossings, int) or crossings != 1:
                v = None
        out[f"F/{ep}"] = v
    return out


def _endpoint_document(path: Path, case: str, ledger: Ledger) -> dict[str, object]:
    """The whole endpoint document of a run: case P's ENDPOINTS line, or case F's record.json (top level kept)."""
    text = ledger.read_text(path)
    if case == "P":
        lines = [ln for ln in text.splitlines() if ln.startswith("ENDPOINTS ")]
        if len(lines) != 1:
            msg = f"{len(lines)} ENDPOINTS lines, expected exactly 1"
            raise ValueError(msg)
        doc = json.loads(lines[0][len("ENDPOINTS ") :])
    else:
        doc = json.loads(text)
    if not isinstance(doc, dict):
        msg = "endpoint record is not a JSON object"
        raise TypeError(msg)
    return doc


def metrics_digest(metrics: object) -> str:
    """The digest platform_study_record.metrics_digest writes into record.json."""
    return hashlib.sha256(json.dumps(metrics, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _check_run_json(
    run_json: dict[str, object], arm: str, case: str, seed: int, energy: int, commit: str
) -> list[str]:
    """Provenance problems of run.json against the frozen expectation for this seed (full-run metadata)."""
    bad: list[str] = []

    def need(key: str, want: object) -> None:
        got = run_json.get(key)
        if got is None:
            bad.append(f"run.json has no {key}")
        elif got != want or type(got) is not type(want):
            bad.append(f"run.json {key} {got!r}, expected {want!r}")

    need("mode", "full")
    need("commit", commit)
    need("energy_mev", energy)
    need("threads", ARM_THREADS[arm])
    need("requested", REQUESTED[case])
    requested, simulated = _int(run_json.get("requested")), _int(run_json.get("simulated"))
    if simulated is None:
        bad.append(f"run.json simulated {run_json.get('simulated')!r} is not an integer")
    elif requested is not None and simulated < requested:
        bad.append(f"run.json simulated {simulated} < requested {requested}")
    binary = run_json.get("binary_sha256")
    if not isinstance(binary, str) or not SHA256_HEX.match(binary):
        bad.append(f"run.json binary_sha256 {binary!r} is not a lowercase sha256")
    for key in ("arm", "case", "seed"):
        if run_json.get(key) is None:
            bad.append(f"run.json has no {key}")
    return bad


def _reconcile_p(doc: dict[str, object], arm: str, seed: int, energy: int) -> list[str]:
    bad: list[str] = []
    if doc.get("layout") != "mcsquare":
        bad.append(f"endpoints layout {doc.get('layout')!r}, expected 'mcsquare'")
    want_label = P_LABEL.format(arm=arm, energy=energy, seed=seed)
    if doc.get("label") != want_label:
        bad.append(f"endpoints label {doc.get('label')!r}, expected {want_label!r}")
    if doc.get("slab_depths_mm") != list(SLAB_DEPTHS[energy]):
        bad.append(f"endpoints slab_depths_mm {doc.get('slab_depths_mm')!r}, expected {list(SLAB_DEPTHS[energy])}")
    return bad


def _is_sha256(v: object) -> bool:
    return isinstance(v, str) and SHA256_HEX.match(v) is not None


def f_record_shape(doc: dict[str, object], commit: str) -> tuple[RecordBinding | None, list[str]]:
    """The record binding of `commit` and the SHAPE problems of a case-F record.json under it (empty when the record
    has exactly the bound schema). A shape problem is an integrity failure (the record cannot come from the bound
    writer) and refuses; it is never an endpoint failure:
      - no binding for the commit;
      - a record carrying some but not all of cfg_sha256 / metrics_sha256 / endpoint_status (mixed schemas);
      - a record of one schema under a commit bound to the other (e.g. the 2f9dab40 schema at a later commit, or a
        record carrying endpoint_status at 2f9dab40);
      - a missing or unexpected top-level field; sha256 not exactly the writer's ten keys of sha256 hex; materials not
        {files, combined_sha256}; metrics not an object (2f9dab40, and a v2 record with endpoint_status ok); a v2
        record whose endpoint_status is not ok/error or whose error record is not the writer's error shape.
    """
    binding = record_binding(commit)
    if binding is None:
        return None, [f"no case-F record binding for commit {commit[:12]}"]
    keys = set(doc)
    v2 = keys & V2_ENDPOINT_FIELDS
    if v2 and v2 != V2_ENDPOINT_FIELDS:
        return binding, [f"mixes record schemas: carries {sorted(v2)} but not {sorted(V2_ENDPOINT_FIELDS - v2)}"]
    schema = RECORD_SCHEMA_V2 if v2 else RECORD_SCHEMA_2F9DAB40
    if schema is not binding.schema:
        why = f"record has the {schema.name} schema, but commit {commit[:12]} is bound to {binding.schema.name} ({binding.name})"
        return binding, [why]
    bad: list[str] = []
    missing, extra = sorted(schema.fields - keys), sorted(keys - schema.fields - schema.optional)
    if missing:
        bad.append(f"record lacks field(s) {missing} of {schema.name}")
    if extra:
        bad.append(f"record has field(s) {extra} that {schema.name} does not write")
    hashes = doc.get("sha256")
    if not isinstance(hashes, dict) or set(hashes) != F_RECORD_SHA256_KEYS or not all(map(_is_sha256, hashes.values())):
        got = sorted(hashes) if isinstance(hashes, dict) else hashes
        bad.append(f"record sha256 is not the writer's {len(F_RECORD_SHA256_KEYS)} sha256 values (keys {got!r})")
    materials = doc.get("materials")
    if (not isinstance(materials, dict) or set(materials) != MATERIALS_KEYS or _int(materials.get("files")) is None
            or not _is_sha256(materials.get("combined_sha256"))):
        bad.append(f"record materials {materials!r} is not {{files, combined_sha256}}")
    status = doc.get("endpoint_status")
    if schema.endpoint_status and status not in ("ok", "error"):
        bad.append(f"record endpoint_status {status!r} is not one the writer writes")
    elif schema.endpoint_status and status == "error":
        if doc.get("metrics") is not None or doc.get("metrics_sha256") is not None or not isinstance(doc.get("endpoint_error"), str):
            bad.append("record endpoint_status error without the writer's error shape (null metrics, endpoint_error)")
    elif "metrics" in keys and not isinstance(doc.get("metrics"), dict):
        bad.append("record metrics is not an object")
    if schema.endpoint_status and status == "ok" and "endpoint_error" in keys:
        bad.append("record endpoint_status ok with an endpoint_error")
    return binding, bad


def _reconcile_f(
    doc: dict[str, object], run_json: dict[str, object] | None, arm: str, seed: int, commit: str
) -> tuple[list[str], str | None]:
    """Case-F record.json against its binding, run.json and the frozen expectation; returns (provenance problems,
    binding name). A record whose shape is not the bound schema refuses (InputError); see f_record_shape.

    Both schemas: study is the binding's, platform/seed/commit are this run's, sha256.binary is run.json's. The 2f9dab40
    schema has no endpoint status and no metrics hash: its endpoints are `metrics` as written (the collector verifies
    the record's file hashes; nothing can verify the metrics). The v2 schema also needs endpoint_status ok,
    metrics_sha256 equal to the digest of metrics and cfg_sha256 equal to sha256.config.
    """
    binding, shape = f_record_shape(doc, commit)
    if shape or binding is None:
        msg = f"{arm}/F/s{seed} record.json: " + "; ".join(shape)
        raise InputError(msg)
    bad: list[str] = []
    if doc.get("study") != binding.study:
        bad.append(f"record study {doc.get('study')!r} is not the study {binding.study!r} bound by {binding.name}")
    for key, want in (("platform", arm), ("seed", seed), ("commit", commit)):
        if doc.get(key) != want or isinstance(doc.get(key), bool):
            bad.append(f"record {key} {doc.get(key)!r}, expected {want!r}")
    hashes: dict[str, object] = doc["sha256"]  # type: ignore[assignment]  # shape checked
    if run_json is not None and run_json.get("binary_sha256") != hashes["binary"]:
        bad.append("record sha256.binary differs from run.json binary_sha256")
    if binding.schema.endpoint_status:
        if doc["endpoint_status"] != "ok":
            bad.append(f"record endpoint_status {doc['endpoint_status']!r}, expected 'ok'")
        elif doc.get("metrics_sha256") != metrics_digest(doc["metrics"]):
            bad.append("record metrics_sha256 does not match its metrics")
        if doc.get("cfg_sha256") != hashes["config"]:
            bad.append("record cfg_sha256 does not equal sha256.config")
    return bad, binding.name


def load_run(
    arm: str, case: str, run_dir: Path, seed: int, energy: int, commit: str, ledger: Ledger
) -> Run:
    """Load one accepted seed directory. A damaged or unprovenanced run stays in the population as unusable;
    only a run.json that contradicts its own directory is fatal."""
    label = f"{arm}/{case}/s{seed}"
    problems: list[str] = []
    provenance: list[str] = []
    run_json: dict[str, object] | None = None
    try:
        parsed = _read_json(run_dir / "run.json", ledger)
        if not isinstance(parsed, dict):
            msg = "run.json is not a JSON object"
            raise TypeError(msg)
        run_json = parsed
    except (OSError, ValueError, TypeError) as e:
        provenance.append(f"run.json unreadable ({type(e).__name__})")
    if run_json is not None:
        for key, want in (("arm", arm), ("case", case), ("seed", seed)):
            if run_json.get(key) is not None and run_json.get(key) != want:
                msg = f"{label}: run.json says {key} {run_json.get(key)!r}, directory says {want!r}"
                raise InputError(msg)
        provenance += _check_run_json(run_json, arm, case, seed, energy, commit)
        status = run_json.get("transport_status")
        if status != TRANSPORT_OK:
            problems.append(f"transport_status {status!r}")
    rec: dict[str, object] | None = None
    binding: str | None = None
    try:
        doc = _endpoint_document(run_dir / ("endpoints.json" if case == "P" else "record.json"), case, ledger)
        if case == "P":
            provenance += _reconcile_p(doc, arm, seed, energy)
            rec = doc
        else:
            bad, binding = _reconcile_f(doc, run_json, arm, seed, commit)
            provenance += bad
            metrics = doc.get("metrics")
            rec = metrics if isinstance(metrics, dict) else None
    except (OSError, ValueError, TypeError) as e:
        problems.append(f"endpoint record unusable ({type(e).__name__}: {e})")
    raw: dict[str, float | None] = {}
    if rec is not None:
        raw = _p_values(rec, energy) if case == "P" else _f_values(rec)
    mode = run_json.get("mode") if run_json is not None else None
    run_commit = run_json.get("commit") if run_json is not None else None
    binary = run_json.get("binary_sha256") if run_json is not None else None
    host = run_json.get("host") if run_json is not None else None
    problems = [*provenance, *problems]
    return Run(
        label, seed, energy, usable=not problems, problems=problems, provenance=provenance,
        mode=mode if isinstance(mode, str) else None, commit=run_commit if isinstance(run_commit, str) else None,
        binary=binary if isinstance(binary, str) and SHA256_HEX.match(binary) else None,
        host=host if isinstance(host, str) else None, binding=binding, raw=raw,
    )


def load_arm_case(
    root: Path, arm: str, case: str, frozen: Frozen, ledger: Ledger,
    listed: dict[int, dict[str, object]] | None = None, *, attested: bool = False,
) -> tuple[list[Run], list[str]]:
    """Accepted runs of one arm and case, and the population issues (missing, rejected, not established, provenance
    failures).

    Only `s<seed>` directories whose seed is expected for this arm and case, written exactly as the workflow writes it
    and not repeated by numeric identity, are opened. Everything else that looks like a seed directory is rejected
    and listed; anything that does not look like one is fatal. `listed` holds the collection manifest's
    not_established entries for this arm and case (seed -> entry): each becomes a population issue (so the contrast is
    PARTIAL), and a seed both collected and listed is fatal. With a manifest (`attested`), every frozen seed must be a
    collected run or listed, else fatal: the collector lists every seed it did not collect. Without one, a missing
    seed is a population issue.
    """
    listed = listed or {}
    case_dir = root / arm / case
    if not case_dir.is_dir():
        msg = f"missing {case_dir}"
        raise InputError(msg)
    expected = frozen.seeds[(arm, case)]
    by_identity: dict[int, list[str]] = {}
    for entry in sorted(case_dir.iterdir()):
        if entry.name.startswith("."):
            continue
        m = _SEED_DIR.match(entry.name)
        if m is None or not entry.is_dir():
            msg = f"unexpected entry {entry} (expected s<seed> directories)"
            raise InputError(msg)
        by_identity.setdefault(int(m.group(1)), []).append(entry.name)
    if not by_identity and not listed:
        msg = f"no s<seed> directories under {case_dir}"
        raise InputError(msg)
    runs: list[Run] = []
    issues: list[str] = []
    rejected: list[str] = []
    for ident, names in sorted(by_identity.items()):
        if len(names) > 1:
            rejected.append(f"{', '.join(names)} (numeric identity {ident} repeated)")
        elif ident not in expected:
            rejected.append(f"{names[0]} (seed not in the frozen list)")
        elif names[0] != f"s{ident}":
            rejected.append(f"{names[0]} (not written as s{ident})")
        else:
            runs.append(load_run(arm, case, case_dir / names[0], ident, expected[ident], frozen.commit, ledger))
    accepted = {r.seed for r in runs}
    both = sorted(accepted & set(listed))
    if both:
        msg = f"{arm} {case}: seed(s) {', '.join(map(str, both))} both collected as runs and listed as not established"
        raise InputError(msg)
    missing = sorted(set(expected) - accepted - set(listed))
    if missing and attested:
        msg = (f"{arm} {case}: frozen seed(s) {', '.join(map(str, missing))} neither collected nor listed as not "
               "established in the collection manifest")
        raise InputError(msg)
    for seed, ne in sorted(listed.items()):
        issues.append(f"{arm}/{case}/s{seed} ({ne['energy']} MeV): not established, {ne['outcome']}: {ne['reason']}")
    if missing:
        issues.append(f"{arm} {case}: {len(missing)} expected seed(s) missing: {', '.join(map(str, missing))}")
    if rejected:
        issues.append(f"{arm} {case}: rejected, not read: " + "; ".join(rejected))
    for r in runs:
        if r.provenance:
            issues.append(f"{r.label}: provenance failure: " + "; ".join(r.provenance))
    return runs, issues


def check_arm_binary(data: dict[tuple[str, str], list[Run]], issues: dict[tuple[str, str], list[str]]) -> None:
    """One binary per arm: every run of an arm (both cases) must carry the same run.json binary_sha256.

    When they differ no run can be called the right one, so every run of that arm becomes a provenance failure.
    """
    for arm in {a for a, _c in data}:
        hashes = {r.binary for c in CASES for r in data.get((arm, c), []) if r.binary is not None}
        if len(hashes) > 1:
            note = f"binary_sha256 differs within arm {arm} ({len(hashes)} distinct values)"
            for c in CASES:
                for r in data.get((arm, c), []):
                    r.provenance.append(note)
                    r.problems.append(note)
                    r.usable = False
                issues[(arm, c)].append(f"{arm} {c}: {note}")


def load_root(
    root: Path, parts: tuple[str, ...], frozen: Frozen, ledger: Ledger,
    listed: dict[tuple[str, str], dict[int, dict[str, object]]] | None = None, *, attested: bool = False,
) -> tuple[dict[tuple[str, str], list[Run]], dict[tuple[str, str], list[str]]]:
    if not root.is_dir():
        msg = f"{root} is not a directory"
        raise InputError(msg)
    data: dict[tuple[str, str], list[Run]] = {}
    issues: dict[tuple[str, str], list[str]] = {}
    for part in parts:
        for arm in PART_ARMS[part]:
            if not (root / arm).is_dir():
                msg = f"missing arm directory {root / arm}"
                raise InputError(msg)
            for case in CASES:
                data[(arm, case)], issues[(arm, case)] = load_arm_case(
                    root, arm, case, frozen, ledger, (listed or {}).get((arm, case)), attested=attested
                )
    check_arm_binary(data, issues)
    return data, issues


def read_manifest(root: Path, ledger: Ledger) -> dict[str, object] | None:
    """collection_manifest.json (written by apples_collect.py) as a JSON object, or None when there is none."""
    path = root / MANIFEST_NAME
    if not path.is_file():
        return None
    try:
        manifest = json.loads(ledger.read_text(path))
    except (OSError, ValueError) as e:
        msg = f"{path} is unusable: {type(e).__name__}: {e}"
        raise InputError(msg) from e
    if not isinstance(manifest, dict):
        msg = f"{path} is not a JSON object"
        raise InputError(msg)
    return manifest


def _file_entries(entries: object, where: str) -> dict[str, str]:
    """A manifest file list ([{path, source, sha256}]) as {path: sha256}; refuses an entry without a path and a sha256,
    or a repeated path (`source` is a record of where the copy came from and is not checked here)."""
    if not isinstance(entries, list):
        msg = f"collection manifest: {where} is not a list"
        raise InputError(msg)
    out: dict[str, str] = {}
    for e in entries:
        if not isinstance(e, dict) or not isinstance(e.get("path"), str) or not isinstance(e.get("sha256"), str) \
                or not SHA256_HEX.match(e["sha256"]):
            msg = f"collection manifest: {where} has an entry without a path and a sha256: {e!r}"
            raise InputError(msg)
        if e["path"] in out:
            msg = f"collection manifest: {e['path']} is listed twice"
            raise InputError(msg)
        out[e["path"]] = e["sha256"]
    return out


def listed_not_established(
    manifest: dict[str, object] | None, frozen: Frozen
) -> dict[tuple[str, str], dict[int, dict[str, object]]]:
    """The manifest's not_established entries, (arm, case) -> {seed: entry}, each checked against the frozen population:
    a frozen (arm, case, seed) at its frozen energy, outcome failed or absent, a reason, files only under
    .not_established/<arm>/<case>/s<seed>/ (none for absent), and no seed listed twice. Raises InputError otherwise."""
    raw = [] if manifest is None else manifest.get("not_established", [])
    if not isinstance(raw, list):
        msg = "collection manifest: not_established is not a list"
        raise InputError(msg)
    out: dict[tuple[str, str], dict[int, dict[str, object]]] = {}
    for e in raw:
        if not isinstance(e, dict):
            msg = f"collection manifest: not_established entry {e!r} is not an object"
            raise InputError(msg)
        arm, case, seed, energy = e.get("arm"), e.get("case"), e.get("seed"), e.get("energy")
        expected = frozen.seeds.get((arm, case)) if isinstance(arm, str) and isinstance(case, str) else None
        if expected is None or _int(seed) is None or seed not in expected:
            msg = f"collection manifest: not_established {arm}/{case}/s{seed} is not in the frozen population"
            raise InputError(msg)
        if _int(energy) is None or energy != expected[seed]:
            msg = f"collection manifest: not_established {arm}/{case}/s{seed} energy {energy!r}, frozen {expected[seed]}"
            raise InputError(msg)
        if e.get("outcome") not in NE_OUTCOMES or not isinstance(e.get("reason"), str) or not e["reason"]:
            msg = f"collection manifest: not_established {arm}/{case}/s{seed} needs an outcome in {NE_OUTCOMES} and a reason"
            raise InputError(msg)
        files = _file_entries(e.get("files"), f"not_established {arm}/{case}/s{seed} files")
        prefix = f"{NOT_ESTABLISHED_DIR}/{arm}/{case}/s{seed}/"
        stray = [p for p in files if not p.startswith(prefix) or ".." in p.split("/")]
        if stray or (e["outcome"] == "absent" and files):
            msg = f"collection manifest: not_established {arm}/{case}/s{seed} lists file(s) outside {prefix} or for an absent run"
            raise InputError(msg)
        cell = out.setdefault((str(arm), str(case)), {})
        if seed in cell:
            msg = f"collection manifest: {arm}/{case}/s{seed} is listed as not established twice"
            raise InputError(msg)
        cell[seed] = {**e, "files_by_path": files}
    return out


def _ne_files(entry: dict[str, object]) -> dict[str, str]:
    """{path: sha256} of one not_established entry, as listed_not_established stored it."""
    files = entry["files_by_path"]
    return files if isinstance(files, dict) else {}


def check_manifest_identity(manifest: dict[str, object], commit: str, parts: tuple[str, ...]) -> None:
    """Refuse, BEFORE any run is read, a manifest that does not attest THIS analysis: an unsupported `schema`, an
    `acquisition_commit` other than the frozen commit in use, `parts` that do not cover the parts being analysed, or a
    case-F record binding other than the one this analysis applies to the commit."""
    schema = manifest.get("schema")
    if _int(schema) is None or schema not in SUPPORTED_MANIFEST_SCHEMAS:
        msg = f"collection manifest: schema {schema!r} is not supported (supported: {list(SUPPORTED_MANIFEST_SCHEMAS)})"
        raise InputError(msg)
    if manifest.get("acquisition_commit") != commit:
        msg = (f"collection manifest: acquisition_commit {manifest.get('acquisition_commit')!r} is not the frozen commit "
               f"in use, {commit}")
        raise InputError(msg)
    mparts = manifest.get("parts")
    if not isinstance(mparts, list) or not all(isinstance(p, str) and p in PART_ARMS for p in mparts) \
            or len(set(mparts)) != len(mparts):
        msg = f"collection manifest: parts {mparts!r} is not a list of distinct parts from {sorted(PART_ARMS)}"
        raise InputError(msg)
    uncovered = sorted(set(parts) - set(mparts))
    if uncovered:
        msg = f"collection manifest: parts {mparts} do not cover the analysed part(s) {uncovered}"
        raise InputError(msg)
    binding = record_binding(commit)
    declared = manifest.get("f_record_binding")
    want = None if binding is None else {"name": binding.name, "schema": binding.schema.name}
    got = {k: declared.get(k) for k in ("name", "schema")} if isinstance(declared, dict) else declared
    if binding is None or got != want:
        msg = f"collection manifest: f_record_binding {got!r} is not this analysis's binding for {commit[:12]}, {want!r}"
        raise InputError(msg)


def verify_manifest(
    root: Path, ledger: Ledger, manifest: dict[str, object] | None,
    listed: dict[tuple[str, str], dict[int, dict[str, object]]],
) -> dict[str, object]:
    """Check collection_manifest.json against the tree, when present.

    Every file it lists (the runs' and provenance `files`, and each not_established entry's diagnostics) must exist and
    hash as recorded, and every file the analysis read must be listed. The not_established files are read through the
    ledger, so the dataset fingerprint covers them. Returns the report entry; raises InputError on a mismatch.
    """
    if manifest is None:
        return {"present": False, "note": "no collection manifest: the provenance of this copy is not attested; "
                "confirmatory claims are withheld for every contrast"}
    files = _file_entries(manifest.get("files"), "files")
    ne_files = {p: h for cell in listed.values() for entry in cell.values() for p, h in _ne_files(entry).items()}
    for rel, want in sorted({**files, **ne_files}.items()):
        target = root / rel
        if target.is_file():
            data = ledger.read_bytes(target) if rel in ne_files else target.read_bytes()
            got: str | None = hashlib.sha256(data).hexdigest()
        else:
            got = None
        if got != want:
            msg = f"collection manifest: {rel} is {'missing' if got is None else 'changed'} since collection"
            raise InputError(msg)
    unlisted = sorted(set(ledger.files) - set(files) - set(ne_files) - {MANIFEST_NAME})
    if unlisted:
        msg = "collection manifest does not list file(s) the analysis read: " + ", ".join(unlisted[:3])
        raise InputError(msg)
    return {"present": True, "files_verified": len(files) + len(ne_files),
            "not_established_files_verified": len(ne_files), "acquisition_commit": manifest.get("acquisition_commit"),
            "schema": manifest.get("schema"), "parts": manifest.get("parts"),
            "f_record_binding": manifest.get("f_record_binding")}


# ----------------------------------------------------------------------------------------------- statistics


def welch(a: list[float], b: list[float]) -> tuple[float, float, float]:
    """Return (mean(a) - mean(b), Welch SE, Welch-Satterthwaite df). Needs at least 2 values per arm."""
    if len(a) < 2 or len(b) < 2:
        msg = "Welch needs at least 2 runs per arm"
        raise ValueError(msg)
    va, vb = st.variance(a) / len(a), st.variance(b) / len(b)
    se = math.sqrt(va + vb)
    if se == 0.0:
        return st.mean(a) - st.mean(b), 0.0, math.inf
    df = (va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1))
    return st.mean(a) - st.mean(b), se, df


def interval(d: float, se: float, df: float, level: float) -> tuple[float, float]:
    h = 0.0 if se == 0.0 else float(student_t.ppf(0.5 + level / 2, df)) * se
    return d - h, d + h


def tost_p(d: float, se: float, df: float, lo: float, hi: float) -> float:
    """max of the two one-sided t p-values for H0: diff <= lo and H0: diff >= hi."""
    if se == 0.0:
        return 0.0 if lo < d < hi else 1.0
    return max(float(student_t.sf((d - lo) / se, df)), float(student_t.sf((hi - d) / se, df)))


def classify(ci90: tuple[float, float], lo: float, hi: float) -> str:
    """equivalent: 90% interval inside the margin; not_equivalent: wholly outside; inconclusive: otherwise."""
    if lo < ci90[0] and ci90[1] < hi:
        return "equivalent"
    if ci90[0] >= hi or ci90[1] <= lo:
        return "not_equivalent"
    return "inconclusive"


def holm_adjust(pvalues: list[float | None]) -> list[float | None]:
    """Holm step-down adjusted p-values over the WHOLE family; a None p (not established) stays in the family as 1.

    Returns None for the entries that had no p-value, so they can never be reported as rejected.
    """
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: (1.0 if pvalues[i] is None else pvalues[i], i))  # type: ignore[type-var]
    adjusted: list[float | None] = [None] * m
    running = 0.0
    for rank, i in enumerate(order):
        p = pvalues[i]
        if p is None:
            continue  # sorted last: nothing after it can be rejected either
        running = max(running, min(1.0, (m - rank) * p))
        adjusted[i] = running
    return adjusted


def values_for(runs: list[Run], spec: Spec) -> tuple[list[float] | None, str]:
    """Analysis-scale values of one endpoint over the runs at its energy, or (None, reason)."""
    sel = [r for r in runs if spec.case == "F" or r.energy == spec.energy]
    out: list[float] = []
    bad: list[str] = []
    for r in sel:
        v = r.raw.get(spec.eid) if r.usable else None
        if v is None:
            bad.append(r.label)
        else:
            out.append(math.log(v) if spec.log else v)
    if bad:
        return None, f"unusable in {len(bad)} of {len(sel)} runs: " + ", ".join(bad[:3]) + (" ..." if len(bad) > 3 else "")
    if len(out) < 2:
        return None, f"{len(out)} run(s), need at least 2 per arm"
    return out, ""


def evaluate(spec: Spec, arm: list[Run], ref: list[Run], *, descriptive: bool) -> Result:
    va, why_a = values_for(arm, spec)
    vb, why_b = values_for(ref, spec)
    if va is None or vb is None:
        reason = "; ".join(w for w in (f"arm {why_a}" if why_a else "", f"reference {why_b}" if why_b else "") if w)
        return Result(spec.eid, spec, "not_established", reason)
    d, se, df = welch(va, vb)
    ci90, ci95 = interval(d, se, df, 0.90), interval(d, se, df, 0.95)
    if descriptive:
        return Result(spec.eid, spec, "descriptive", "", d, se, df, ci90, ci95)
    return Result(
        spec.eid, spec, classify(ci90, spec.lo, spec.hi), "", d, se, df, ci90, ci95, tost_p(d, se, df, spec.lo, spec.hi)
    )


def analyse_contrast(
    data: dict[tuple[str, str], list[Run]],
    arm: str,
    ref: str,
    specs: list[Spec],
    *,
    confirmatory: bool,
    issues: dict[tuple[str, str], list[str]],
    collection: list[str] | None = None,
) -> dict[str, object]:
    """One contrast. A confirmatory contrast whose population is not exactly the frozen one, or whose collection is not
    verified (`collection`: reasons that apply to every contrast, e.g. no collection manifest), is PARTIAL: it is
    analysed descriptively (estimates and intervals, nothing classified) and emits no joint claim and no Holm decision."""
    partial = [*(collection or []), *partial_reasons(issues, (arm, ref))]
    withheld = confirmatory and bool(partial)
    descriptive = not confirmatory or withheld
    results = [evaluate(s, data[(arm, s.case)], data[(ref, s.case)], descriptive=descriptive) for s in specs]
    claims: dict[str, object] = {}
    if confirmatory and not withheld:
        adjusted = holm_adjust([r.p_tost for r in results])
        for r, hp in zip(results, adjusted, strict=True):
            r.holm_p = hp
            if r.outcome == "not_established":
                r.holm_decision = "not_established"
            else:
                r.holm_decision = "equivalent" if hp is not None and hp < ALPHA else "not shown"
        claims = {
            "joint_claim_equivalent_on_all_endpoints": all(r.outcome == "equivalent" for r in results),
            "n_equivalent_unadjusted": sum(r.outcome == "equivalent" for r in results),
            "n_equivalent_holm": sum(r.holm_decision == "equivalent" for r in results),
            "n_not_established": sum(r.outcome == "not_established" for r in results),
        }
    return {"results": results, "claims": claims, "partial_reasons": partial, "claims_withheld": withheld}


# ------------------------------------------------------------------------------------------ design completeness


def cell_counts(runs: list[Run], case: str) -> dict[str, int]:
    """Runs PRESENT (usable or not) per energy cell: 100/150/200 for case P, one "200" cell for case F."""
    keys = [str(e) for e in SLAB_DEPTHS] if case == "P" else ["200"]
    counts = dict.fromkeys(keys, 0)
    for r in runs:
        key = str(r.energy) if case == "P" else "200"
        counts[key] = counts.get(key, 0) + 1
    return counts


def partial_reasons(issues: dict[tuple[str, str], list[str]], arms: tuple[str, ...]) -> list[str]:
    """Why the population of these arms (both cases) is not exactly the frozen one; empty when it is."""
    return [reason for arm in arms for case in CASES for reason in issues[(arm, case)]]


# ----------------------------------------------------------------------------------------------- reporting


def _fmt(spec: Spec, x: float) -> str:
    return f"{math.exp(x):.4f}" if spec.log else f"{x:+.4f}"


def _margin_text(spec: Spec) -> str:
    if spec.log:
        return f"ratio [{math.exp(spec.lo):.3f}, {math.exp(spec.hi):.3f}]"
    return f"[{spec.lo:+.3f}, {spec.hi:+.3f}]"


def _json_num(x: float | None) -> float | None:
    return x if x is not None and math.isfinite(x) else None


def result_json(r: Result) -> dict[str, object]:
    """JSON-safe endpoint record. Estimates for log endpoints are on the log scale; `ratio_*` give the ratio scale.

    A non-finite value (df is infinite when both arms have zero variance) is written as null.
    """
    s = r.spec
    ratio = (lambda x: math.exp(x)) if s.log else (lambda x: x)
    return {
        "endpoint": r.eid,
        "scale": "log_ratio" if s.log else "difference",
        "outcome": r.outcome,
        "reason": r.reason,
        "estimate": _json_num(r.estimate),
        "se": _json_num(r.se),
        "df": _json_num(r.df),
        "ci90": None if r.ci90 is None else [_json_num(x) for x in r.ci90],
        "ci95": None if r.ci95 is None else [_json_num(x) for x in r.ci95],
        "estimate_reported_scale": None if r.estimate is None else _json_num(ratio(r.estimate)),
        "ci90_reported_scale": None if r.ci90 is None else [_json_num(ratio(x)) for x in r.ci90],
        "ci95_reported_scale": None if r.ci95 is None else [_json_num(ratio(x)) for x in r.ci95],
        "margin": None if r.outcome == "descriptive" else [_json_num(s.lo), _json_num(s.hi)],
        "p_tost": _json_num(r.p_tost),
        "holm_adjusted_p": _json_num(r.holm_p),
        "holm_decision": r.holm_decision,
    }


def markdown_table(name: str, confirmatory: bool, analysed: dict[str, object]) -> list[str]:
    results: list[Result] = analysed["results"]  # type: ignore[assignment]
    claims: dict[str, object] = analysed["claims"]  # type: ignore[assignment]
    partial: list[str] = analysed["partial_reasons"]  # type: ignore[assignment]
    withheld = bool(analysed["claims_withheld"])
    lines = [f"### {name}", ""]
    if partial:
        lines += [
            "**PARTIAL**: the population is not exactly the frozen one (" + str(RUNS_PER_CELL) + " runs per arm, case "
            "and energy planned), or the collection is not verified." + (" Confirmatory claims WITHHELD: no joint claim, no Holm decision, no equivalence "
                                      "classification; the estimates below are descriptive only." if withheld else ""),
            "",
            *[f"- {reason}" for reason in partial],
            "",
        ]
    if confirmatory and not withheld:
        lines += [
            (f"Joint claim (equivalent on all {len(results)} endpoints, unadjusted alpha {ALPHA}): "
            f"**{'TRUE' if claims['joint_claim_equivalent_on_all_endpoints'] else 'FALSE'}**; "
            f"equivalent unadjusted {claims['n_equivalent_unadjusted']}, "
            f"equivalent after Holm {claims['n_equivalent_holm']}, not established {claims['n_not_established']}."),
            "",
        ]
    elif not confirmatory:
        lines += ["Descriptive only: estimates and intervals; no margin, no equivalence claim.", ""]
    lines += [
        "Ring estimates and intervals are ratios; their SE is the SE of the log ratio.",
        "",
        "| endpoint | estimate | SE | df | 90% | 95% | margin | outcome | Holm-adjusted decision |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        s = r.spec
        if r.estimate is None or r.se is None or r.df is None or r.ci90 is None or r.ci95 is None:
            lines.append(f"| {r.eid} | n/a | n/a | n/a | n/a | n/a | {_margin_text(s)} | {r.outcome} ({r.reason}) | n/a |")
            continue
        if r.outcome == "descriptive":
            holm = "n/a (descriptive)"
        elif r.holm_p is not None:
            holm = f"{r.holm_decision} (p_adj {r.holm_p:.4f})"
        else:
            holm = r.holm_decision
        margin = "n/a" if r.outcome == "descriptive" else _margin_text(s)
        lines.append(
            f"| {r.eid} | {_fmt(s, r.estimate)} | {r.se:.4f} | {r.df:.1f} | "
            f"[{_fmt(s, r.ci90[0])}, {_fmt(s, r.ci90[1])}] | [{_fmt(s, r.ci95[0])}, {_fmt(s, r.ci95[1])}] | "
            f"{margin} | {r.outcome} | {holm} |"
        )
    return [*lines, ""]


def run_analysis(
    root: Path, parts: tuple[str, ...], *, commit: str = ACQUISITION_COMMIT, sources: dict[str, str] | None = None
) -> tuple[list[str], dict[str, object]]:
    """Return (markdown lines, JSON-able document). Raises InputError when the inputs are unusable.

    `commit` is the acquisition commit whose workflow files freeze the population; `sources` (tests only) supplies
    those workflow texts instead of `git show`.
    """
    frozen = expected_population(commit, sources=sources)
    ledger = Ledger(root)
    if not root.is_dir():
        msg = f"{root} is not a directory"
        raise InputError(msg)
    manifest_doc = read_manifest(root, ledger)
    if manifest_doc is not None:
        check_manifest_identity(manifest_doc, commit, parts)  # before any run is read or any inference is made
    listed = listed_not_established(manifest_doc, frozen)
    data, issues = load_root(root, parts, frozen, ledger, listed, attested=manifest_doc is not None)
    manifest = verify_manifest(root, ledger, manifest_doc, listed)
    collection = [] if manifest["present"] else [NO_MANIFEST_REASON]
    fingerprint = ledger.fingerprint()
    specs = all_specs()
    binding = record_binding(commit)
    md = ["# Same-host comparison: parts " + ", ".join(parts), "", IDENTITY_SOURCE, ""]
    md.append(f"Acquisition commit (frozen run lists): {commit}")
    md.append(f"Dataset fingerprint: sha256 {fingerprint['sha256']} over {fingerprint['n_files']} input file(s)")
    md.append(
        f"Collection manifest: verified, {manifest['files_verified']} file(s)" if manifest["present"]
        else f"Collection manifest: {manifest['note']}"
    )
    if collection:
        md.append("Confirmatory claims: WITHHELD for every contrast (descriptive only): " + "; ".join(collection))
    if binding is not None:
        md.append(f"Case-F record schema: {binding.schema.name} (binding {binding.name}). Verified: "
                  + "; ".join(binding.schema.verified) + "."
                  + (" Not verifiable: " + "; ".join(binding.schema.not_verifiable) + "." if binding.schema.not_verifiable else ""))
    md.append("")
    doc: dict[str, object] = {"parts": list(parts), "alpha": ALPHA, "family_size": FAMILY_SIZE,
                           "identity_source": IDENTITY_SOURCE, "runs_per_cell": RUNS_PER_CELL,
                           "acquisition_commit": commit, "dataset_fingerprint": fingerprint,
                           "collection_manifest": manifest,
                           "confirmatory_claims_possible": not collection, "claims_withheld_for_all": collection,
                           "f_record_binding": None if binding is None else {
                               "name": binding.name, "commit": binding.commit, "study": binding.study,
                               "schema": binding.schema.name, "verified": list(binding.schema.verified),
                               "not_verifiable": list(binding.schema.not_verifiable)},
                           "contrasts": []}
    md.append("Runs per arm and case (unusable runs stay in the population):")
    md.append("")
    runs_doc: dict[str, object] = {}
    for (arm, case), runs in sorted(data.items()):
        bad = [f"{r.label}: {'; '.join(r.problems)}" for r in runs if not r.usable]
        modes = sorted({r.mode or "unknown" for r in runs})
        commits = sorted({r.commit or "unknown" for r in runs})
        counts = cell_counts(runs, case)
        bound = sorted({r.binding for r in runs if r.binding})
        runs_doc[f"{arm}/{case}"] = {
            "n": len(runs), "expected": len(frozen.seeds[(arm, case)]), "unusable": bad, "modes": modes,
            "commits": commits, "per_energy": counts, "issues": issues[(arm, case)],
            "binary_sha256": sorted({r.binary for r in runs if r.binary}),
            "hosts": sorted({r.host for r in runs if r.host}), "record_bindings_used": bound,
        }
        md.append(
            f"- {arm}/{case}: {len(runs)} of {len(frozen.seeds[(arm, case)])} expected runs "
            f"({', '.join(f'{e} MeV: {n}' for e, n in counts.items())}), {len(bad)} unusable; "
            f"mode {'/'.join(modes)}, commit {'/'.join(c[:12] for c in commits)}"
            + (f", record binding {'/'.join(bound)}" if bound else "")
            + "".join(f"\n  - {b}" for b in bad)
        )
    md.append("")
    analysed_arms = {arm for part in parts for arm in PART_ARMS[part]}
    ne_doc = [
        {"arm": arm, "case": case, "seed": seed, "energy_mev": e["energy"], "outcome": e["outcome"], "reason": e["reason"],
         "files": sorted(_ne_files(e))}
        for (arm, case), cell in sorted(listed.items()) if arm in analysed_arms for seed, e in sorted(cell.items())
    ]
    md.append(f"Runs not established (listed by the collector, never read as runs; their contrasts are PARTIAL): {len(ne_doc)}")
    md += [f"- {n['arm']}/{n['case']}/s{n['seed']} ({n['energy_mev']} MeV): {n['outcome']}: {n['reason']}" for n in ne_doc]
    md.append("")
    doc["runs"] = runs_doc
    doc["not_established"] = ne_doc
    contrasts_doc: list[object] = []
    for part in parts:
        for name, arm, ref, confirmatory in CONTRASTS_BY_PART[part]:
            analysed = analyse_contrast(data, arm, ref, specs, confirmatory=confirmatory, issues=issues,
                                        collection=collection)
            md += markdown_table(name, confirmatory, analysed)
            results: list[Result] = analysed["results"]  # type: ignore[assignment]
            contrasts_doc.append(
                {
                    "name": name,
                    "arm": arm,
                    "reference": ref,
                    "kind": "confirmatory" if confirmatory else "descriptive",
                    "claims": analysed["claims"],
                    "claims_withheld": analysed["claims_withheld"],
                    "partial": bool(analysed["partial_reasons"]),
                    "partial_reasons": analysed["partial_reasons"],
                    "endpoints": [result_json(r) for r in results],
                }
            )
    doc["contrasts"] = contrasts_doc
    return md, doc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--json", type=Path, default=None, help="also write the full result here (allow_nan=False)")
    ap.add_argument("--parts", default="A,B", help="comma-separated parts to analyse (default A,B)")
    ap.add_argument("--acquisition-commit", default=ACQUISITION_COMMIT, help="full sha whose workflow files freeze the runs")
    args = ap.parse_args(argv)
    try:
        parts = tuple(p.strip() for p in args.parts.split(","))
        if not parts or any(p not in PART_ARMS for p in parts) or len(set(parts)) != len(parts):
            msg = f"--parts must be distinct values from {sorted(PART_ARMS)}, got {args.parts!r}"
            raise InputError(msg)
        md, doc = run_analysis(args.root, parts, commit=args.acquisition_commit)
        if args.json is not None:
            args.json.write_text(json.dumps(doc, indent=1, sort_keys=True, allow_nan=False) + "\n")
    except (InputError, OSError) as e:
        print(f"apples_analyse: unusable inputs or output: {e}", file=sys.stderr)
        return 2
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    sys.exit(main())
