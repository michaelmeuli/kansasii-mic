import pandas as pd

from kansasii_mic.overview import mgit_overview, mhk_overview


def _mgit(rows: list[tuple[int, str, float, str]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["TNR", "antibiotic", "concentration_mg_l", "erg"])


def test_mgit_overview_classifies_mic_with_clsi() -> None:
    df = _mgit(
        [(1, "Rifampicin", 1.0, "R"), (1, "Rifampicin", 4.0, "S"),  # MIC 4 -> R (>= 2)
         (2, "Rifampicin", 1.0, "S"), (2, "Rifampicin", 4.0, "S"),  # MIC 1 -> S
         (3, "Amikacin", 1.0, "R"), (3, "Amikacin", 4.0, "S"), (3, "Amikacin", 20.0, "S"),  # MIC 4 -> S
         (4, "Amikacin", 4.0, "R"), (4, "Amikacin", 20.0, "R"),  # MIC > 20: 40 < 64, unplaceable
         (5, "Isoniazid", 0.1, "R"), (5, "Isoniazid", 1.0, "S")]  # no CLSI breakpoint
    )
    out = mgit_overview(df).set_index("antibiotic")
    assert list(out.index) == ["Amikacin", "Rifampicin"]
    rif, ami = out.loc["Rifampicin"], out.loc["Amikacin"]
    assert (rif.nr_susceptible, rif.nr_intermediate, rif.nr_resistant) == (1, 0, 1)
    assert (ami.nr_susceptible, ami.nr_resistant, ami.nr_unclassifiable) == (1, 0, 1)
    assert rif.susceptible_concentration == "<= 1" and rif.resistant_concentration == ">= 2"


def test_mhk_overview_uses_clsi_breakpoints_and_skips_others() -> None:
    df = pd.DataFrame(
        {
            "antibiotic": ["Rifampicin"] * 3 + ["Ethambutol"],
            "point_estimate": [0.5, 2.0, 1.0, 4.0],
        }
    )
    out = mhk_overview(df)
    assert list(out["antibiotic"]) == ["Rifampicin"]
    assert (out.nr_susceptible[0], out.nr_intermediate[0], out.nr_resistant[0]) == (2, 0, 1)


def test_mgit_mic_counts_lowest_s_and_above_range() -> None:
    from kansasii_mic.overview import mgit_mic_counts

    df = pd.DataFrame(
        {
            "TNR": [1, 1, 2, 2, 3, 3],
            "antibiotic": ["A"] * 6,
            "concentration_mg_l": [1.0, 4.0] * 3,
            "erg": ["R", "S", "S", "S", "R", "R"],
        }
    )
    out = mgit_mic_counts(df).set_index("mic_label")["n"].to_dict()
    assert out == {"1": 1, "4": 1, ">4": 1}


def test_mhk_intermediate_label_matches_readme() -> None:
    from kansasii_mic.overview import _intermediate_label

    assert _intermediate_label(16, 64) == "32"
    assert _intermediate_label(1, 8) == "2-4"
    assert _intermediate_label(1, 2) == "--"


def test_mgit_ranking_counts_clsi_resistant_drugs() -> None:
    from kansasii_mic.mgit import tnr_mgit_resistance_ranking

    rows = [
        (1, "Rifampicin", 1.0, "R"), (1, "Rifampicin", 4.0, "S"),  # CLSI R
        (1, "Isoniazid", 0.1, "R"), (1, "Isoniazid", 1.0, "S"),  # no breakpoint
        (2, "Rifampicin", 1.0, "S"), (2, "Rifampicin", 4.0, "S"),
        (2, "Amikacin", 1.0, "R"), (2, "Amikacin", 4.0, "R"),  # R calls, MIC > 4: not CLSI R
    ]
    df = _mgit(rows)
    df["NR"] = df["TNR"]
    df["PROBENNUMMER"] = "x"
    out = tnr_mgit_resistance_ranking(df).set_index("TNR")
    assert out.loc[1, "n_R_clsi"] == 1 and out.loc[2, "n_R_clsi"] == 0
    assert out.loc[1, "rank_most_resistant"] == 1
