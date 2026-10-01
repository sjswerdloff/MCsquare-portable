"""Tests for apples_seeds_check: both seed rules, fail-closed parsing, and the committed workflow files."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from apples_seeds_check import ALL_ARMS, DEFAULT_FILES, check, main, parse


def line(arm: str, full: str, smoke: str = "") -> str:
    return f'          - {{arm: {arm}, case: P, full: "{full}", smoke: "{smoke}"}}\n'


CLEAN = (
    line("A-up", "100:960001 150:960011", "100:960901")
    + line("A-port", "100:961001 150:961011", "100:961901")
    + line("B-up", "100:962001", "100:962901")
    + line("B-pgcc", "100:963001", "100:963901")
    + line("B-picc", "100:964001", "100:964901")
)


def test_clean_lists_pass():
    assert check({"w.yml": CLEAN}) == []


def test_repeat_within_an_arm_fails():
    text = CLEAN + line("A-up", "200:960001")
    errors = check({"w.yml": text})
    assert any("960001 repeats" in e for e in errors)


def test_repeat_across_files_and_modes_fails():
    a = line("A-up", "100:960001") + line("A-port", "100:961001")
    b = line("B-up", "100:962001") + line("B-pgcc", "100:963001") + line("B-picc", "100:964001", "100:960001")
    assert any("960001 repeats" in e for e in check({"a.yml": a, "b.yml": b}))


def test_portable_pair_a_multiple_of_1e4_apart_fails():
    text = CLEAN.replace("100:963001", "100:971001")  # B-pgcc 971001 - A-port 961001 = 1e4
    errors = check({"w.yml": text})
    assert any("961001" in e and "971001" in e and "multiple of 10000" in e for e in errors)


def test_portable_seed_against_earlier_portable_block_fails():
    text = CLEAN.replace("100:964001", "100:960002").replace("100:960001", "100:960003")
    # B-picc 960002 vs earlier 940002 / 950002: both differences are multiples of 1e4
    errors = check({"w.yml": text})
    assert any("960002" in e and "earlier portable seed 950002" in e for e in errors)


def test_repeating_an_earlier_portable_seed_fails():
    text = CLEAN.replace("100:963001", "100:940010")
    assert any("repeats an earlier portable seed" in e for e in check({"w.yml": text}))


def test_upstream_seed_a_multiple_of_1e4_from_portable_is_allowed():
    """A-up is not portable: 960001 - 950001 = 1e4 is no state collision."""
    assert 960001 in {s for _a, s, _w in parse("w.yml", CLEAN)[0]}
    assert check({"w.yml": CLEAN}) == []


def test_earlier_blocks_colliding_with_each_other_are_not_reported():
    """940001 and 950001 differ by 1e4 and are never compared (plan); only new seeds are checked."""
    errors = check({"w.yml": CLEAN})
    assert not any("940001" in e for e in errors)


def test_new_portable_pair_not_a_multiple_passes():
    text = CLEAN + line("B-pgcc", "100:963002")  # 963002 - 961001 = 2001
    assert check({"w.yml": text}) == []


def test_run_spec_on_a_line_without_an_arm_fails_closed():
    text = CLEAN + '        smoke: "100:961905"\n'
    assert any("names no arm" in e for e in check({"w.yml": text}))


def test_no_specs_and_missing_arms_fail_closed():
    errors = check({"w.yml": "name: nothing here\n"})
    assert any("no run specs" in e for e in errors)
    assert any("arms with no runs" in e for e in errors)


def test_unknown_arm_fails():
    assert any("unknown arm" in e for e in check({"w.yml": CLEAN + line("C-port", "100:965001")}))


@pytest.mark.parametrize(
    "token, found",
    [("100:961001", True), ("10:30:45", False), ("C:961001", False), ("100:9610011", False)],
)
def test_spec_token_shape(token, found):
    runs, _ = parse("w.yml", f"- {{arm: A-port, full: \"{token}\"}}")
    assert bool(runs) is found


def test_committed_workflows_pass_and_carry_the_full_plan():
    for f in DEFAULT_FILES:
        assert f.is_file(), f
    sources = {f.name: f.read_text(encoding="utf-8") for f in DEFAULT_FILES}
    assert check(sources) == []
    runs = [r for name, text in sources.items() for r in parse(name, text)[0]]
    by_arm = {arm: sorted(s for a, s, _w in runs if a == arm) for arm in ALL_ARMS}
    for arm, seeds in by_arm.items():
        base = {"A-up": 960000, "A-port": 961000, "B-up": 962000, "B-pgcc": 963000, "B-picc": 964000}[arm]
        full = [base + o + i for o in (0, 10, 20, 30) for i in range(1, 9)]
        smoke = [base + 901, base + 911, base + 921, base + 931]
        assert seeds == sorted(full + smoke), arm


def test_main_reports_ok_on_committed_files(capsys):
    assert main([]) == 0
    assert capsys.readouterr().out.strip().endswith("0 violation(s)")
