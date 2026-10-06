import math

from kansasii_mic.breakpoints import clsi_normalized


def test_clsi_normalized_anchors() -> None:
    # Amikacin: S <= 16, R >= 64
    assert clsi_normalized("Amikacin", math.log2(16)) == 0
    assert clsi_normalized("Amikacin", math.log2(64)) == 1
    assert clsi_normalized("Amikacin", math.log2(32)) == 0.5


def test_clsi_normalized_sides() -> None:
    below = clsi_normalized("Amikacin", math.log2(1))
    above = clsi_normalized("Amikacin", math.log2(128))
    assert below is not None and below < 0
    assert above is not None and above > 1


def test_clsi_normalized_none_without_breakpoint() -> None:
    assert clsi_normalized("Isoniazid", 1.0) is None
