import sys
from pathlib import Path

from kansasii_mic.loading import SCREENING_MAP_TNR_COLUMNS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))  # vendored kansasii_paths.py
from kansasii_paths import TNR_COLUMNS  # noqa: E402


def test_tnr_columns_match_kansasii_paths() -> None:
    assert SCREENING_MAP_TNR_COLUMNS == TNR_COLUMNS
