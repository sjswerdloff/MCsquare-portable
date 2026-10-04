"""Tests for report_tables.py: the Welch arithmetic, invalid-value handling, and reproduction of the committed tables."""

import json
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import report_tables as rt


def test_welch_matches_hand_computation():
    a, b = [1.0, 2.0, 3.0, 4.0], [0.0, 0.0, 1.0, 1.0]
    d, se, df = rt.welch(a, b)
    va, vb = (5 / 3) / 4, (1 / 3) / 4
    assert d == pytest.approx(2.5 - 0.5)
    assert se == pytest.approx(math.sqrt(va + vb))
    assert df == pytest.approx((va + vb) ** 2 / (va**2 / 3 + vb**2 / 3))


def test_welch_refuses_single_run():
    with pytest.raises(ValueError):
        rt.welch([1.0], [1.0, 2.0])


@pytest.mark.parametrize("bad", [None, float("nan"), float("inf"), 0.0, -1.0])
def test_invalid_ring_value_is_unavailable_not_imputed(bad):
    good = [{"ring_100_80_200": 0.01}, {"ring_100_80_200": 0.011}]
    broken = [{"ring_100_80_200": 0.01}, {"ring_100_80_200": bad} if bad is not None else {}]
    assert "UNAVAILABLE" in rt.row("x", good, broken, "ring_100_80_200", log=True)


def test_zero_is_valid_on_the_linear_scale():
    a = [{"R80": 0.0}, {"R80": 0.0}]
    assert "UNAVAILABLE" not in rt.row("x", a, a, "R80", log=False)


def test_committed_data_reproduces_published_values(capsys):
    rt.main([])
    out = capsys.readouterr().out
    # Values published in the report draft (sections 5 and 6) and in #32 c27679 / c27757.
    for expected in (
        "| upstream - portable | R80 | +0.0008 mm (SE 0.0038)",
        "| upstream - portable | ring_200_40_80 | ratio 1.0070",
        "| portable - TOPAS | R80 | -0.5118 mm (SE 0.0029)",
        "| portable - TOPAS | ring_200_40_80 | ratio 0.8887 (SE of log 0.0032)",
    ):
        assert expected in out


def test_every_committed_record_names_its_source():
    root = Path(rt.__file__).parent / "report_data"
    files = sorted(root.glob("*.jsonl"))
    assert len(files) == 3
    for f in files:
        for line in f.read_text().splitlines():
            rec = json.loads(line)
            assert rec["source"] and rec["seed"] and rec["endpoints"]["R80_multiple_crossings"] is False
