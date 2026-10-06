import pandas as pd

from kansasii_mic.plotting import _tie_spread


def test_tie_spread_is_deterministic_and_centred() -> None:
    values = pd.Series([0.0, 0.0, 0.0, 1.0])
    y = _tie_spread(values)
    assert list(y) == list(_tie_spread(values))
    assert list(y[:3]) == [-0.05, 0.0, 0.05]
    assert y.iloc[3] == 0.0


def test_tie_spread_caps_span_for_large_groups() -> None:
    y = _tie_spread(pd.Series([0.0] * 41))
    assert y.max() - y.min() <= 0.8 + 1e-9
