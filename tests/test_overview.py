import pandas as pd

from kansasii_mic.overview import mgit_overview, mhk_overview


def test_mgit_overview_counts_worst_call_per_isolate():
    df = pd.DataFrame(
        {
            "TNR": [1, 1, 2, 2, 3],
            "antibiotic": ["Amikacin"] * 5,
            "concentration_mg_l": [1.0, 4.0, 1.0, 4.0, 1.0],
            "int_erg": ["R", "S", "I", "S", "U"],
        }
    )
    row = mgit_overview(df).iloc[0]
    assert (row.nr_susceptible, row.nr_intermediate, row.nr_resistant) == (0, 1, 1)
    assert row.susceptible_concentration == "4"
    assert row.resistant_concentration == "1"


def test_mhk_overview_uses_clsi_breakpoints_and_skips_others():
    df = pd.DataFrame(
        {
            "antibiotic": ["Rifampicin"] * 3 + ["Ethambutol"],
            "point_estimate": [0.5, 2.0, 1.0, 4.0],
        }
    )
    out = mhk_overview(df)
    assert list(out["antibiotic"]) == ["Rifampicin"]
    assert (out.nr_susceptible[0], out.nr_intermediate[0], out.nr_resistant[0]) == (2, 0, 1)


def test_mgit_mic_counts_lowest_s_and_above_range():
    from kansasii_mic.overview import mgit_mic_counts

    df = pd.DataFrame(
        {
            "TNR": [1, 1, 2, 2, 3, 3],
            "antibiotic": ["A"] * 6,
            "concentration_mg_l": [1.0, 4.0] * 3,
            "int_erg": ["R", "S", "S", "S", "R", "R"],
        }
    )
    out = mgit_mic_counts(df).set_index("mic_label")["n"].to_dict()
    assert out == {"1": 1, "4": 1, ">4": 1}


def test_mhk_intermediate_label_matches_readme():
    from kansasii_mic.overview import _intermediate_label

    assert _intermediate_label(16, 64) == "32"
    assert _intermediate_label(1, 8) == "2-4"
    assert _intermediate_label(1, 2) == "--"
