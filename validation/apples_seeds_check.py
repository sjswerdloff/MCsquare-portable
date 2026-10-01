"""Seed check for the same-host comparison workflows (docs/apples_to_apples_plan.md, parts A and B).

Reads the run lists out of the part A and part B workflow files and fails (exit 1) if

  * any seed appears twice, in either file or across both, in full or smoke mode, or repeats an earlier portable seed;
  * any two PORTABLE seeds differ by a nonzero multiple of 10000. Portable MCsquare seeds thread t with state
    RNG_Seed + 1e4 * t (report section 7), so such a pair can reuse a generator state. The portable arms are A-port,
    B-pgcc and B-picc, plus the earlier portable seeds below. Pairs of two EARLIER seeds are not checked: the
    9400xx and 9500xx blocks already collide and are never compared (plan, "Seeds").

Run specs are the tokens ENERGY:SEED on a matrix line that names its arm, e.g.
    - {arm: A-port, case: P, full: "100:961001 100:961002", smoke: "100:961901"}
A run spec on a line with no arm, a file with no run specs, or a missing arm is an error, so a format change cannot
make the check pass by finding nothing.

Usage: python apples_seeds_check.py [workflow.yml ...]   (default: both apples workflows in this repository)
"""

import itertools
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_FILES = [
    REPO / ".gitea" / "workflows" / "apples-a-windows.yml",
    REPO / ".gitea" / "workflows" / "apples-b-linux.yml",
]
ALL_ARMS = frozenset({"A-up", "A-port", "B-up", "B-pgcc", "B-picc"})
PORTABLE_ARMS = frozenset({"A-port", "B-pgcc", "B-picc"})
EARLIER_PORTABLE = frozenset(
    list(range(900101, 900105)) + list(range(940001, 940045)) + list(range(950001, 950005))
)
STATE_STRIDE = 10000

SPEC = re.compile(r"(?<![\w:])(\d{2,3}):(\d{6})(?![\w:])")
ARM = re.compile(r"\barm:\s*([A-Za-z]+-[A-Za-z]+)\b")


def parse(name: str, text: str) -> tuple[list[tuple[str, int, str]], list[str]]:
    """Return ([(arm, seed, where)], errors) for one workflow file's text."""
    runs, errors = [], []
    for lineno, line in enumerate(text.splitlines(), 1):
        specs = SPEC.findall(line)
        if not specs:
            continue
        m = ARM.search(line)
        if not m:
            errors.append(f"{name}:{lineno}: run spec(s) on a line that names no arm")
            continue
        arm = m.group(1)
        if arm not in ALL_ARMS:
            errors.append(f"{name}:{lineno}: unknown arm {arm!r}")
            continue
        runs += [(arm, int(seed), f"{name}:{lineno}") for _energy, seed in specs]
    if not runs:
        errors.append(f"{name}: no run specs found")
    return runs, errors


def check(sources: dict[str, str], required_arms: frozenset[str] = ALL_ARMS) -> list[str]:
    """Return a list of violations (empty when the run lists satisfy both rules)."""
    runs, errors = [], []
    for name, text in sources.items():
        r, e = parse(name, text)
        runs += r
        errors += e
    missing = required_arms - {arm for arm, _s, _w in runs}
    if missing:
        errors.append(f"arms with no runs: {', '.join(sorted(missing))}")

    seen: dict[int, str] = {}
    for arm, seed, where in runs:
        if seed in seen:
            errors.append(f"seed {seed} repeats ({arm} at {where}; first at {seen[seed]})")
        else:
            seen[seed] = f"{arm} at {where}"
        if arm in PORTABLE_ARMS and seed in EARLIER_PORTABLE:
            errors.append(f"seed {seed} ({arm} at {where}) repeats an earlier portable seed")

    new_portable = sorted({seed: arm for arm, seed, _w in runs if arm in PORTABLE_ARMS}.items())
    for (a, arm_a), (b, arm_b) in itertools.combinations(new_portable, 2):
        if (b - a) % STATE_STRIDE == 0:
            errors.append(f"portable seeds {a} ({arm_a}) and {b} ({arm_b}) differ by a multiple of {STATE_STRIDE}")
    for seed, arm in new_portable:
        for old in sorted(EARLIER_PORTABLE):
            if seed != old and (seed - old) % STATE_STRIDE == 0:
                errors.append(f"portable seed {seed} ({arm}) and earlier portable seed {old} differ by a multiple of {STATE_STRIDE}")
    return errors


def main(argv: list[str]) -> int:
    files = [Path(a) for a in argv] or DEFAULT_FILES
    sources = {}
    for f in files:
        if not f.is_file():
            print(f"missing workflow file: {f}")
            return 1
        sources[f.name] = f.read_text(encoding="utf-8")
    errors = check(sources, ALL_ARMS if not argv else frozenset())
    n = sum(len(SPEC.findall(t)) for t in sources.values())
    for e in errors:
        print(e)
    print(f"{'FAIL' if errors else 'OK'}: {n} run specs in {len(sources)} file(s), {len(errors)} violation(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
