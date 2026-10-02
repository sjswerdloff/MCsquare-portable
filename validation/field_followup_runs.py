"""The frozen run lists of the same-host comparison, part E (docs/field_100_150_plan.md): field edge at 100, 150 MeV.

This module is the ONE place the arms, seeds, histories, threads, hosts and binaries are written down. The two run
workflows ask it for their seeds, and the analysis imports it, so a run list cannot differ between them.

Usage (the workflows):
  python field_followup_runs.py seeds <arm> <energy> <full|smoke>     prints the seeds, space separated
  python field_followup_runs.py histories <full|smoke>               prints the number of histories, as an integer
  python field_followup_runs.py binary <arm>                         prints the pinned binary sha256
  python field_followup_runs.py threads <arm>                        prints the arm's thread count
  python field_followup_runs.py check                                seed screen; prints OK or the violations, exit 1

The seed screen is the one parts A and B used (apples_seeds_check.py): every seed of this part is unique, repeats no
seed of parts A and B and no earlier portable seed, and no two portable-PCG runs share a generator start
(initstate, stream) on any thread at any of the 10 batch calls. It is screening, not a proof of independence.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass

import apples_seeds_check as sc

ENERGIES = (100, 150)
SIDE_MM = 150.0
RUNS_PER_CELL = 8
HISTORIES = {"full": 60_000_000, "smoke": 100_000}
# Portable's source at the acquisition commit must equal this commit's: the binaries below were built from it.
SOURCE_COMMIT = "2f9dab404cea02f5072352f7ae5773dca27eb2a1"


@dataclass(frozen=True)
class Arm:
    block: int  # first three digits of every seed of the arm
    threads: int
    host: str
    portable: bool  # PCG generator (the seed screen models it); the upstream arms use another generator
    directory: str  # the arm's directory name in the run tree
    binary_sha256: str  # the binary every part A or B run of this arm recorded


ARMS = {
    "A-up": Arm(960, 4, "DESKTOP-SR5GKKA", False, "aup", "f3a28398a399224c1e20d6e54e945571f6ee17693f4380ee7ad9c0e3ef84b966"),
    "A-port": Arm(961, 4, "DESKTOP-SR5GKKA", True, "apt", "de93714ba856fc30605b1533247d5fa2be4e17c4c3e7cdc296bb5317c2aa60f4"),
    "B-up": Arm(962, 3, "DESKTOP-5H86O9N", False, "bup", "b9106df25839a0e295565c0845dfbaff3ac3c80efa7f839d3f764eb9c87d768b"),
    "B-pgcc": Arm(963, 3, "DESKTOP-5H86O9N", True, "bpg", "a21d406efac7fd67c8405df1bb86a86b3f419cb57a0ee82e12750e5824e16b7f"),
    "B-picc": Arm(964, 3, "DESKTOP-5H86O9N", True, "bpi", "16e5d89911eba190cb9560be90c3e4b395c55399c5da6236362a629892e45ac1"),
}
CONTRASTS = (("A-port", "A-up", True), ("B-pgcc", "B-up", True), ("B-picc", "B-up", True), ("B-pgcc", "B-picc", False))
# Last three digits: full runs 041-048 (100 MeV) and 051-058 (150 MeV); the smoke run of each is 941 and 951.
_FIRST = {100: 41, 150: 51}


def seeds(arm: str, energy: int, mode: str) -> tuple[int, ...]:
    """The seeds of one arm, energy and mode. Raises ValueError for anything outside the frozen design."""
    if arm not in ARMS or energy not in ENERGIES or mode not in HISTORIES:
        msg = f"not in the design: arm {arm!r}, energy {energy!r}, mode {mode!r}"
        raise ValueError(msg)
    base = ARMS[arm].block * 1000 + _FIRST[energy]
    return (base + 900,) if mode == "smoke" else tuple(base + i for i in range(RUNS_PER_CELL))


def plan_name(energy: int) -> str:
    return f"E{energy}_S{SIDE_MM:.0f}"


def all_runs() -> list[tuple[str, int, str, int]]:
    """Every (arm, energy, mode, seed) of the part, full and smoke."""
    return [(a, e, m, s) for a in ARMS for e in ENERGIES for m in HISTORIES for s in seeds(a, e, m)]


def check(apples_sources: dict[str, str] | None = None) -> list[str]:
    """Violations of the seed rules (empty when there are none). `apples_sources` replaces the two committed apples
    workflow files (tests only)."""
    if apples_sources is None:
        apples_sources = {f.name: f.read_text(encoding="utf-8") for f in sc.DEFAULT_FILES}
    errors, apples_runs, *_ = sc._check(apples_sources, sc.ALL_ARMS)
    errors = [f"parts A and B run lists: {e}" for e in errors]
    seen: dict[int, str] = {s: f"part A/B {arm} at {where}" for arm, _c, _m, _e, s, where in apples_runs}
    seen.update({s: f"earlier portable seed {s}" for s in sc.EARLIER_PORTABLE})
    owner: dict[tuple[int, int], str] = {}
    for s in sorted(sc.EARLIER_PORTABLE):
        for pair in sc.generator_pairs(s, sc.EARLIER_THREADS):
            owner.setdefault(pair, f"earlier portable seed {s}")
    for arm, _c, _m, _e, s, _where in apples_runs:
        if arm in sc.PORTABLE_ARMS:
            for pair in sc.generator_pairs(s, ARMS[arm].threads):
                owner.setdefault(pair, f"part A/B {arm} seed {s}")
    for arm, energy, mode, s in all_runs():
        label = f"{arm} {energy} MeV {mode} seed {s}"
        if not 100_000 <= s <= 999_999 or s // 1000 != ARMS[arm].block:
            errors.append(f"{label}: outside the arm's six-digit block {ARMS[arm].block}xxx")
        if s in seen:
            errors.append(f"{label}: repeats {seen[s]}")
        seen.setdefault(s, label)
        if ARMS[arm].portable:
            for pair, (t, c) in sc.generator_pairs(s, ARMS[arm].threads).items():
                if pair in owner and owner[pair] != label:
                    errors.append(f"{label}: generator start {pair} (thread {t}, call {c}) is also used by {owner[pair]}")
                owner.setdefault(pair, label)
    return errors


def main(argv: list[str]) -> int:
    try:
        if len(argv) == 4 and argv[0] == "seeds":
            print(" ".join(str(s) for s in seeds(argv[1], int(argv[2]), argv[3])))
        elif len(argv) == 2 and argv[0] == "histories":
            print(HISTORIES[argv[1]])
        elif len(argv) == 2 and argv[0] == "binary":
            print(ARMS[argv[1]].binary_sha256)
        elif len(argv) == 2 and argv[0] == "threads":
            print(ARMS[argv[1]].threads)
        elif argv == ["check"]:
            errors = check()
            for e in errors:
                print(e)
            print(f"{'FAIL' if errors else 'OK'}: {len(all_runs())} runs of part E screened, {len(errors)} violation(s)")
            return 1 if errors else 0
        else:
            print(__doc__, file=sys.stderr)
            return 2
    except (KeyError, ValueError) as e:
        print(f"field_followup_runs.py: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
