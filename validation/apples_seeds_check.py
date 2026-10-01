"""Seed and run-list check for the same-host comparison workflows (parts A and B).

This is SCREENING for identical generator starts and for a run list that differs from the plan. It is not a proof
that the runs are statistically independent.

Every matrix row must be a single-line flow mapping with the keys arm, case, full and smoke, e.g.
    - {arm: A-port, case: P, full: "100:961001 100:961002", smoke: "100:961901"}
Every whitespace-separated token in every full and smoke field must be ENERGY:SEED with ENERGY in 100/150/200 and a
six-digit SEED (case F: ENERGY 200 only). Nothing is skipped: a malformed token, an empty field, a missing field, a
duplicate row or an unknown arm is an error, and the number of tokens discovered in the fields must equal the number
consumed (and the number found by a loose scan of the whole file).

Expected run counts, per arm in A-up, A-port, B-up, B-pgcc, B-picc:
    case P  full 24 (8 per energy at 100, 150, 200)   smoke 3 (one per energy)
    case F  full  8 (all at 200)                      smoke 1 (200)

Generator model. src/compute_simulation.c, Simulation_loop, thread t calls
    pcg32_srandom_r(rng, RNG_Seed + 1e4 * t + 1e5 * Num_call, t)
(2f9dab40: lines 396-397). Num_call is a function-local static incremented once per call (line 373-374). Run_simulation
(compute_simulation.c:38-82) calls Simulation_loop once per batch whenever Compute_stat_uncertainty is true (default
True, data_config.c:86), which is the case for these runs: MIN_NUM_BATCH = 10 (define.h:63). Num_batch only grows, and
the batch counter only resets, in the branch taken when Stat_uncertainty != 0 (compute_simulation.c:55-83); its default
is 0.0 (data_config.c:87) and the workflows do not set it, so the count is exactly K = 10 calls, Num_call = 1..10. A
nonzero Stat_uncertainty would make K unbounded; this check does not model that. A PCG32 generator is the pair
(initstate, stream) = (RNG_Seed + 1e4 * t + 1e5 * c, t). Two runs reuse a generator if they share a pair for some thread
t and calls c1, c2 (possibly different calls: seeds 1e5 apart collide at c and c+1). The check computes those pairs over
t in [0, threads) and c in [1, K] for every portable-PCG run (arms A-port, B-pgcc, B-picc and the earlier portable
seeds) and refuses any shared pair, naming both runs, the thread and the calls. Seeds that differ by a multiple of 1e4
but not of 1e5 share an initstate only at DIFFERENT t, hence different streams, so they are accepted. The upstream arms
(A-up, B-up) use a different generator and are checked for duplicate seeds only. This is screening for identical
generator starts across the batch calls, not proof of independence.

Usage: python apples_seeds_check.py [workflow.yml ...]   (default: both apples workflows in this repository)
"""

import re
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_FILES = [
    REPO / ".gitea" / "workflows" / "apples-a-windows.yml",
    REPO / ".gitea" / "workflows" / "apples-b-linux.yml",
]
ALL_ARMS = frozenset({"A-up", "A-port", "B-up", "B-pgcc", "B-picc"})
PORTABLE_ARMS = frozenset({"A-port", "B-pgcc", "B-picc"})
CASES = ("P", "F")
ENERGIES = (100, 150, 200)

# Earlier portable-PCG seeds that no apples seed may repeat, by origin:
#   planning   900101-4, 900111-4, 900121-2   (topas-mcsquare-planning.yml, seeds 9001xx; 900121-2 is 100 MeV)
#   diagnostic 940001-44                       (topas-mcsquare-diagnostic.yml; 940001-24 ct-grid cells, 940031-44 later cells)
#   what-if    950001-4                        (what-if runs named in the plan's "Seeds" section; no workflow is committed)
EARLIER_PORTABLE = frozenset(
    list(range(900101, 900105))
    + list(range(900111, 900115))
    + list(range(900121, 900123))
    + list(range(940001, 940045))
    + list(range(950001, 950005))
)
EARLIER_THREADS = 24  # those workflows run THREADS: "24"
MAX_BATCH_CALLS = 10  # K: Simulation_loop calls per process = MIN_NUM_BATCH (define.h:63), Stat_uncertainty default 0.0
STATE_STRIDE = 10000

# Expected tokens per (case, mode) and energy.
EXPECTED: dict[tuple[str, str], dict[int, int]] = {
    ("P", "full"): {100: 8, 150: 8, 200: 8},
    ("P", "smoke"): {100: 1, 150: 1, 200: 1},
    ("F", "full"): {200: 8},
    ("F", "smoke"): {200: 1},
}

TOKEN = re.compile(r"^(100|150|200):(\d{6})$")
LOOSE_SPEC = re.compile(r"(?<![\w:])\d{2,3}:\d{5,}(?![\w:])")
FIELD = re.compile(r'(\w+):\s*("[^"]*"|[^,}"]*)')
THREADS_ENV = re.compile(r'^\s*THREADS:\s*"?(\d+)"?\s*$', re.MULTILINE)


def generator_pairs(seed: int, threads: int, batch_calls: int | None = None) -> dict[tuple[int, int], tuple[int, int]]:
    """The (initstate, stream) pairs passed to pcg32_srandom_r by one run, mapped to (thread, call).

    Threads are 0..threads-1 and calls 1..batch_calls (default MAX_BATCH_CALLS).
    """
    calls = MAX_BATCH_CALLS if batch_calls is None else batch_calls
    return {(seed + STATE_STRIDE * t + 10 * STATE_STRIDE * c, t): (t, c) for t in range(threads) for c in range(1, calls + 1)}


def parse(name: str, text: str) -> tuple[list[tuple[str, str, str, int, int, str]], list[str], int, int]:
    """Parse one workflow's matrix rows completely.

    Returns (runs, errors, discovered, consumed); a run is (arm, case, mode, energy, seed, where).
    """
    runs: list[tuple[str, str, str, int, int, str]] = []
    errors: list[str] = []
    discovered = consumed = 0
    for lineno, line in enumerate(text.splitlines(), 1):
        where = f"{name}:{lineno}"
        row = re.search(r"\{\s*arm:[^{}]*\}", line)
        if not row:
            if LOOSE_SPEC.search(line):
                errors.append(f"{where}: run spec(s) on a line that is not a matrix row naming an arm")
            elif re.search(r"\barm:", line) and re.search(r"\b(full|smoke):", line):
                errors.append(f"{where}: matrix row not parseable as a single-line mapping")
            continue
        fields = dict(FIELD.findall(row.group(0)))
        arm, case = fields.get("arm", "").strip(), fields.get("case", "").strip()
        if arm not in ALL_ARMS:
            errors.append(f"{where}: unknown arm {arm!r}")
        if case not in CASES:
            errors.append(f"{where}: unknown case {case!r}")
        for mode in ("full", "smoke"):
            if mode not in fields:
                errors.append(f"{where}: row has no {mode} field")
                continue
            raw = fields[mode].strip()
            if not (raw.startswith('"') and raw.endswith('"') and len(raw) >= 2):
                errors.append(f"{where}: {mode} field is not a quoted string")
                continue
            tokens = raw[1:-1].split()
            discovered += len(tokens)
            if not tokens:
                errors.append(f"{where}: {arm} {case} {mode} cell is empty")
            for tok in tokens:
                m = TOKEN.match(tok)
                if not m:
                    errors.append(f"{where}: malformed run token {tok!r} in {arm} {case} {mode}")
                    continue
                energy, seed = int(m.group(1)), int(m.group(2))
                if case == "F" and energy != 200:
                    errors.append(f"{where}: case F token {tok!r} is not at 200 MeV")
                    continue
                consumed += 1
                runs.append((arm, case, mode, energy, seed, where))
    if not runs and not errors:
        errors.append(f"{name}: no run specs found")
    return runs, errors, discovered, consumed


def count_errors(runs: list[tuple[str, str, str, int, int, str]], required_arms: frozenset[str]) -> list[str]:
    errors = []
    cells: Counter[tuple[str, str, str, int]] = Counter((a, c, m, e) for a, c, m, e, _s, _w in runs)
    present = {(a, c) for a, c, *_ in runs}
    for arm in sorted(required_arms):
        for case in CASES:
            if (arm, case) not in present:
                errors.append(f"arm {arm} case {case} has no runs")
    for arm, case in sorted(present):
        for mode in ("full", "smoke"):
            expected = EXPECTED[(case, mode)]
            for energy in sorted(set(expected) | {e for a, c, m, e in cells if (a, c, m) == (arm, case, mode)}):
                got, want = cells[(arm, case, mode, energy)], expected.get(energy, 0)
                if got != want:
                    errors.append(f"{arm} {case} {mode} at {energy} MeV: {got} runs, expected {want}")
    return errors


def check(sources: dict[str, str], required_arms: frozenset[str] = ALL_ARMS) -> list[str]:
    """Return the violations (empty when the run lists satisfy every rule)."""
    return _check(sources, required_arms)[0]


def _check(sources: dict[str, str], required_arms: frozenset[str]) -> tuple[list[str], list, int, int, int]:
    runs: list[tuple[str, str, str, int, int, str]] = []
    errors: list[str] = []
    discovered = consumed = loose = 0
    threads_of: dict[str, int] = {}
    for name, text in sources.items():
        r, e, d, c = parse(name, text)
        runs += r
        errors += e
        discovered += d
        consumed += c
        loose += sum(len(LOOSE_SPEC.findall(line)) for line in text.splitlines())
        m = THREADS_ENV.search(text)
        if m:
            for arm in {a for a, *_ in r}:
                threads_of[arm] = int(m.group(1))
        elif r:
            errors.append(f"{name}: no THREADS value found, so the generator model cannot be computed")
    if discovered != consumed:
        errors.append(f"discovered {discovered} run tokens but consumed {consumed}")
    if loose != discovered:
        errors.append(f"a loose scan finds {loose} run-spec-shaped tokens but the matrix rows hold {discovered}")
    errors += count_errors(runs, required_arms)

    seen: dict[int, str] = {}
    for arm, _case, _mode, _energy, seed, where in runs:
        if seed in seen:
            errors.append(f"seed {seed} repeats ({arm} at {where}; first at {seen[seed]})")
        else:
            seen[seed] = f"{arm} at {where}"

    owner: dict[tuple[int, int], tuple[str, int, int]] = {}
    for seed in sorted(EARLIER_PORTABLE):
        for pair, (t, c) in generator_pairs(seed, EARLIER_THREADS).items():
            owner.setdefault(pair, (f"earlier portable seed {seed}", t, c))
    for arm, _case, _mode, _energy, seed, where in runs:
        if arm not in PORTABLE_ARMS or arm not in threads_of:
            continue
        label = f"{arm} seed {seed} at {where}"
        reported: set[str] = set()
        for pair, (t, c) in sorted(generator_pairs(seed, threads_of[arm]).items()):
            if pair in owner:
                other, ot, oc = owner[pair]
                if other not in reported:
                    reported.add(other)
                    what = "repeats an earlier portable seed" if other.startswith("earlier") else "shares a generator start"
                    errors.append(
                        f"{label} {what}: (initstate, stream) {pair} at thread {t} call {c} "
                        f"also used by {other} at thread {ot} call {oc}"
                    )
            else:
                owner[pair] = (label, t, c)
    return errors, runs, discovered, consumed, loose


def main(argv: list[str]) -> int:
    files = [Path(a) for a in argv] or DEFAULT_FILES
    sources = {}
    for f in files:
        if not f.is_file():
            print(f"missing workflow file: {f}")
            return 1
        sources[f.name] = f.read_text(encoding="utf-8")
    errors, _runs, discovered, consumed, _loose = _check(sources, ALL_ARMS if not argv else frozenset())
    for e in errors:
        print(e)
    print(
        f"{'FAIL' if errors else 'OK'}: discovered {discovered} run tokens, consumed {consumed}, "
        f"in {len(sources)} file(s), {len(errors)} violation(s)"
    )
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
