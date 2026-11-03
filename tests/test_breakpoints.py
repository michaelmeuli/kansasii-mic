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


def test_categorize_interval() -> None:
    from kansasii_mic.breakpoints import categorize_interval

    nan = float("nan")
    assert categorize_interval("Amikacin", 1, 4) == "S"
    assert categorize_interval("Amikacin", 4, 20) is None  # 20 > S limit 16, below R limit
    assert categorize_interval("Rifampicin", 1, 4) == "R"
    assert categorize_interval("Rifampicin", 1, nan) == "R"
    assert categorize_interval("Clarithromycin", 4, 16) is None
    assert categorize_interval("Isoniazid", 0.1, 1) == "S"  # non-CLSI threshold 8
    assert categorize_interval("Isoniazid", 1, nan) is None  # MIC > 1 may be on either side of 8
    assert categorize_interval("Isoniazid", 10, nan) == "R"
    assert categorize_interval("Clofazimin", 1, 4) is None
    assert categorize_interval("Clofazimin", 2, 4) == "R"
    assert categorize_interval("Streptomycin", 1, 4) is None


def test_non_clsi_thresholds() -> None:
    from kansasii_mic.breakpoints import categorize, has_cutoff

    assert categorize("Isoniazid", 8.0) == "S" and categorize("Isoniazid", 16.0) == "R"
    assert categorize("Clofazimin", 2.0) == "S" and categorize("Clofazimin", 4.0) == "R"
    assert has_cutoff("Ethambutol") and has_cutoff("Rifampicin") and not has_cutoff("Streptomycin")
