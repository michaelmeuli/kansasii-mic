import math

import pytest

from kansasii_mic.parsing import MicParseError, parse_mic


def test_plain_number():
    p = parse_mic("0.5")
    assert p.point_estimate == 0.5
    assert p.censored is None
    assert p.denominator is None
    assert math.isclose(p.log2, -1.0)


def test_left_censored_lt():
    p = parse_mic("<0.25")
    assert p.point_estimate == 0.25
    assert p.censored == "left"


def test_left_censored_lte():
    p = parse_mic("<=0.12")
    assert p.point_estimate == 0.12
    assert p.censored == "left"


def test_right_censored():
    p = parse_mic(">8")
    assert p.point_estimate == 8.0
    assert p.censored == "right"


def test_range_geometric_mean():
    p = parse_mic("4-8")
    assert p.is_range
    assert math.isclose(p.point_estimate, math.sqrt(4 * 8))
    assert p.low == 4.0
    assert p.high == 8.0


def test_combo_ratio():
    p = parse_mic("0.12/2.38")
    assert p.point_estimate == 0.12
    assert p.denominator == 2.38


def test_combo_censored():
    p = parse_mic(">8/152")
    assert p.point_estimate == 8.0
    assert p.denominator == 152.0
    assert p.censored == "right"


def test_combo_range():
    p = parse_mic("0.25/4.75-0.5/9.5")
    assert p.is_range
    assert math.isclose(p.point_estimate, math.sqrt(0.25 * 0.5))
    assert math.isclose(p.denominator, math.sqrt(4.75 * 9.5))


def test_known_typo_combo_range_is_repaired():
    p = parse_mic("0.12/2.38-0.25-4.75")
    assert p.repaired_from == "0.12/2.38-0.25-4.75"
    assert p.is_range
    assert math.isclose(p.point_estimate, math.sqrt(0.12 * 0.25))
    assert math.isclose(p.denominator, math.sqrt(2.38 * 4.75))


def test_empty_raises():
    with pytest.raises(MicParseError):
        parse_mic("")


def test_garbage_raises():
    with pytest.raises(MicParseError):
        parse_mic("not-a-number")
