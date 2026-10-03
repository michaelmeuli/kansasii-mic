import numpy as np
import pandas as pd

from kansasii_mic.outliers import per_antibiotic_outliers


def _df(antibiotic, log2_values):
    return pd.DataFrame(
        {
            "TNR": range(len(log2_values)),
            "NR": range(len(log2_values)),
            "PROBENNUMMER": [f"Mkan329-{i:03d}" for i in range(len(log2_values))],
            "antibiotic": antibiotic,
            "mhk_raw": "x",
            "point_estimate": [2.0**v for v in log2_values],
            "log2_mic": log2_values,
            "erg": "S",
        }
    )


def test_modified_z_defined_when_mad_is_zero():
    # >50% tied at the median -> MAD is 0, but the isolate at 2 must still get a score.
    out = per_antibiotic_outliers(_df("Clarithromycin", [-2.0] * 5 + [-1.0, 1.0]))
    assert out["modified_z"].notna().all()
    assert (out.loc[out["log2_mic"] == -2.0, "modified_z"] == 0).all()
    assert out.loc[out["log2_mic"] == 1.0, "modified_z"].iloc[0] > 0


def test_modified_z_nan_when_all_identical():
    out = per_antibiotic_outliers(_df("Rifabutin", [-2.0] * 4))
    assert out["modified_z"].isna().all()
