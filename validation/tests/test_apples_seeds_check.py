"""Tests for apples_seeds_check: complete parsing, run counts, the generator model, and the committed workflows."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import apples_seeds_check as S
from apples_seeds_check import (
    ALL_ARMS,
    DEFAULT_FILES,
    check,
    generator_pairs,
    main,
    parse,
)

BASES = {"A-up": 960000, "A-port": 961000, "B-up": 962000, "B-pgcc": 963000, "B-picc": 964000}


def full_p(base: int) -> str:
    return " ".join(f"{e}:{base + o + i}" for e, o in ((100, 0), (150, 10), (200, 20)) for i in range(1, 9))


def smoke_p(base: int) -> str:
    return f"100:{base + 901} 150:{base + 911} 200:{base + 921}"


def full_f(base: int) -> str:
    return " ".join(f"200:{base + 30 + i}" for i in range(1, 9))


def row(arm: str, case: str, full: str, smoke: str) -> str:
    return f'          - {{arm: {arm}, case: {case}, full: "{full}", smoke: "{smoke}"}}\n'


def arm_rows(arm: str) -> str:
    b = BASES[arm]
    return row(arm, "P", full_p(b), smoke_p(b)) + row(arm, "F", full_f(b), f"200:{b + 931}")


def workflow(arms: tuple[str, ...], threads: int = 4) -> str:
    return f'env:\n  THREADS: "{threads}"\nmatrix:\n' + "".join(arm_rows(a) for a in arms)


A_FILE = workflow(("A-up", "A-port"), 4)
B_FILE = workflow(("B-up", "B-pgcc", "B-picc"), 3)
CLEAN = {"a.yml": A_FILE, "b.yml": B_FILE}


def mutate(file: str, old: str, new: str) -> dict[str, str]:
    text = CLEAN[file]
    assert old in text, old
    return {**CLEAN, file: text.replace(old, new, 1)}


def has(errors: list[str], *needles: str) -> bool:
    return any(all(n in e for n in needles) for e in errors)


def test_clean_lists_pass_and_discovered_equals_consumed():
    assert check(CLEAN) == []
    errors, runs, discovered, consumed, loose = S._check(CLEAN, ALL_ARMS)
    assert errors == []
    assert discovered == consumed == loose == len(runs) == 5 * (24 + 3 + 8 + 1)


def test_seven_digit_token_is_refused_not_dropped():
    errors = check(mutate("a.yml", "100:961001", "100:9610010"))
    assert has(errors, "malformed run token", "9610010")
    assert has(errors, "discovered", "consumed")


def test_empty_smoke_cell_is_refused():
    errors = check(mutate("a.yml", f'smoke: "200:{BASES["A-port"] + 931}"', 'smoke: ""'))
    assert has(errors, "A-port F smoke cell is empty")


def test_row_missing_a_field_is_refused():
    text = A_FILE.replace(f', smoke: "200:{BASES["A-up"] + 931}"', "", 1)
    assert has(check({**CLEAN, "a.yml": text}), "no smoke field")


def test_wrong_count_is_refused():
    errors = check(mutate("b.yml", " 100:963008", ""))
    assert has(errors, "B-pgcc P full at 100 MeV: 7 runs, expected 8")


def test_missing_whole_row_is_refused():
    text = "".join(line for line in A_FILE.splitlines(True) if "case: F" not in line or "A-port" not in line)
    assert has(check({**CLEAN, "a.yml": text}), "arm A-port case F has no runs")


def test_extra_run_at_an_unexpected_energy_is_refused():
    errors = check(mutate("a.yml", "100:960001", "250:960001"))
    assert has(errors, "malformed run token", "250:960001")


def test_case_f_must_be_200_mev():
    errors = check(mutate("a.yml", "200:960031", "100:960031"))
    assert has(errors, "case F token", "not at 200 MeV")


def test_duplicate_row_is_refused():
    text = A_FILE + row("A-up", "F", full_f(BASES["A-up"]), f"200:{BASES['A-up'] + 931}")
    errors = check({**CLEAN, "a.yml": text})
    assert has(errors, "A-up F full at 200 MeV: 16 runs, expected 8")


def test_unknown_arm_is_refused():
    assert has(check({**CLEAN, "a.yml": A_FILE + row("C-port", "P", "100:965001", "100:965901")}), "unknown arm")


def test_spec_on_a_line_that_is_not_a_row_fails_closed():
    errors = check({**CLEAN, "a.yml": A_FILE + '        smoke: "100:961905"\n'})
    assert has(errors, "not a matrix row")
    assert has(errors, "loose scan")


def test_no_specs_and_missing_arms_fail_closed():
    errors = check({"w.yml": "name: nothing here\n"})
    assert has(errors, "no run specs")
    assert has(errors, "has no runs")


def test_duplicate_seed_is_refused_within_and_across_files():
    assert has(check(mutate("a.yml", "100:960002", "100:960001")), "seed 960001 repeats")
    assert has(check(mutate("b.yml", "100:962002", "100:960001")), "seed 960001 repeats")


def test_repeating_an_earlier_portable_seed_is_refused():
    errors = check(mutate("b.yml", "100:963001", "100:940010"))
    assert has(errors, "repeats an earlier portable seed", "940010")


@pytest.mark.parametrize("seed", [900111, 900114, 900121, 900122, 950004])
def test_inventory_includes_the_planning_and_whatif_blocks(seed):
    assert seed in S.EARLIER_PORTABLE
    assert has(check(mutate("b.yml", "100:963001", f"100:{seed}")), "repeats an earlier portable seed", str(seed))


def test_earlier_seed_with_upstream_arm_is_not_a_generator_collision():
    """Upstream arms use a different generator; they are checked for duplicate seeds only."""
    errors = check(mutate("b.yml", "100:962001", "100:940010"))
    assert not has(errors, "earlier portable seed")


def test_seeds_1e4_apart_are_different_generators_and_accepted():
    """Formerly refused. 971001 and 961001 share an initstate at different t, so their streams differ."""
    text_errors = check(mutate("b.yml", "100:963001", "100:971001"))
    assert text_errors == []
    a, b = generator_pairs(961001, 4), generator_pairs(971001, 3)
    assert {i for i, _t in a} & {i for i, _t in b}  # the same initstate occurs...
    assert not a & b  # ...never with the same stream


def test_generator_pairs_follow_the_c_constructor():
    assert generator_pairs(961001, 4) == {(961001 + 100000 + 10000 * t, t) for t in range(4)}


def test_a_shared_pair_would_be_refused(monkeypatch):
    """The pair model is what decides: make two distinct seeds map to the same pair and the check must refuse."""
    real = S.generator_pairs
    monkeypatch.setattr(S, "generator_pairs", lambda seed, threads, num_call=1: real(961001 if seed == 963001 else seed, threads, num_call))
    assert has(check(CLEAN), "shares a generator start")


def test_committed_workflows_pass_and_carry_the_full_plan():
    for f in DEFAULT_FILES:
        assert f.is_file(), f
    sources = {f.name: f.read_text(encoding="utf-8") for f in DEFAULT_FILES}
    assert check(sources) == []
    runs = [r for name, text in sources.items() for r in parse(name, text)[0]]
    by_arm = {arm: sorted(s for a, _c, _m, _e, s, _w in runs if a == arm) for arm in ALL_ARMS}
    for arm, seeds in by_arm.items():
        base = BASES[arm]
        full = [base + o + i for o in (0, 10, 20, 30) for i in range(1, 9)]
        smoke = [base + 901, base + 911, base + 921, base + 931]
        assert seeds == sorted(full + smoke), arm


def test_main_reports_ok_on_committed_files(capsys):
    assert main([]) == 0
    out = capsys.readouterr().out.strip()
    assert "discovered 180 run tokens, consumed 180" in out
    assert out.endswith("0 violation(s)")
