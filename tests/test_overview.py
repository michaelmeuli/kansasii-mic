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
         (5, "Isoniazid", 0.1, "R"), (5, "Isoniazid", 1.0, "S"),  # MIC 1 <= 8 -> S
         (6, "Isoniazid", 0.1, "R"), (6, "Isoniazid", 1.0, "R"),  # tested only to 1: unplaceable
         (7, "Bedaquilin", 1.0, "R"), (7, "Bedaquilin", 4.0, "S")]  # no cutoff at all
    )
    out = mgit_overview(df).set_index("antibiotic")
    assert list(out.index) == ["Rifampicin", "Amikacin", "Isoniazid"]
    inh = out.loc["Isoniazid"]
    assert (inh.nr_susceptible, inh.nr_resistant, inh.nr_unclassifiable) == (1, 0, 1)
    assert (inh.susceptible_concentration, inh.intermediate_concentration, inh.resistant_concentration) == ("<= 8", "--", "> 8")
    rif, ami = out.loc["Rifampicin"], out.loc["Amikacin"]
    assert (rif.nr_susceptible, rif.nr_intermediate, rif.nr_resistant) == (1, 0, 1)
    assert (ami.nr_susceptible, ami.nr_resistant, ami.nr_unclassifiable) == (1, 0, 1)
    assert rif.susceptible_concentration == "<= 1" and rif.resistant_concentration == ">= 2"


def test_sort_overview_sections_and_order(tmp_path) -> None:  # type: ignore[no-untyped-def]
    from kansasii_mic.overview import overview_section, sort_overview
    from kansasii_mic.plotting import plot_overview_table

    drugs = ["Streptomycin", "Amikacin", "Rifampicin", "Clofazimin", "Ethambutol", "Clarithromycin", "Ethionamid"]
    df = pd.DataFrame({"antibiotic": drugs, "x": range(len(drugs))})
    assert list(sort_overview(df)["antibiotic"]) == [
        "Clarithromycin", "Rifampicin", "Amikacin", "Ethambutol", "Streptomycin", "Ethionamid", "Clofazimin",
    ]
    assert [overview_section(d) for d in ("Rifampicin", "Amikacin", "Clofazimin")] == ["First Line", "Second Line", "No CLSI Breakpoints"]

    mic = pd.DataFrame({"antibiotic": ["Rifampicin", "Amikacin", "Ethambutol"], "point_estimate": [0.5, 8.0, 4.0]})
    out = plot_overview_table(mhk_overview(mic), "t", tmp_path / "o.png", "foot")
    assert out.exists()


def test_mhk_overview_uses_clsi_breakpoints_and_skips_others() -> None:
    df = pd.DataFrame(
        {
            "antibiotic": ["Rifampicin"] * 3 + ["Ethambutol"] * 2 + ["Bedaquilin"],
            "point_estimate": [0.5, 2.0, 1.0, 8.0, 16.0, 4.0],
        }
    )
    out = mhk_overview(df).set_index("antibiotic")
    assert list(out.index) == ["Rifampicin", "Ethambutol"]  # Bedaquilin: no cutoff; CLSI drugs first
    rif, emb = out.loc["Rifampicin"], out.loc["Ethambutol"]
    assert (rif.nr_susceptible, rif.nr_intermediate, rif.nr_resistant) == (2, 0, 1)
    assert (emb.nr_susceptible, emb.nr_resistant) == (1, 1)  # 8 == threshold is S, 16 > 8 is R


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
    assert out.loc[1, "n_R_mic"] == 1 and out.loc[2, "n_R_mic"] == 0
    assert out.loc[1, "rank_most_resistant"] == 1


def test_mgit_mic_counts_above_bin_uses_isolates_own_highest_tested() -> None:
    from kansasii_mic.overview import mgit_mic_counts

    df = pd.DataFrame(
        {
            "TNR": [1, 1, 2, 2, 2, 2],
            "antibiotic": ["A"] * 6,
            "concentration_mg_l": [0.1, 1.0, 0.1, 1.0, 3.0, 10.0],
            "erg": ["R", "R", "R", "S", "S", "S"],
        }
    )
    out = mgit_mic_counts(df)
    assert list(out["mic_label"]) == ["0.1", "1", ">1", "3", "10", ">10"]
    assert out.set_index("mic_label")["n"].to_dict() == {"0.1": 0, "1": 1, "3": 0, "10": 0, ">1": 1, ">10": 0}


def test_mgit_overview_ethambutol_and_clofazimin_thresholds() -> None:
    df = _mgit(
        [(1, "Ethambutol", 5.0, "R"), (1, "Ethambutol", 12.5, "S"),  # (5, 12.5] straddles 8
         (2, "Ethambutol", 5.0, "S"), (2, "Ethambutol", 12.5, "S"),
         (3, "Clofazimin", 0.25, "R"), (3, "Clofazimin", 1.0, "S"),  # MIC 1 <= 2
         (4, "Clofazimin", 1.0, "R"), (4, "Clofazimin", 2.0, "R"), (4, "Clofazimin", 4.0, "R")]  # MIC > 4 -> R
    )
    out = mgit_overview(df).set_index("antibiotic")
    emb, cfz = out.loc["Ethambutol"], out.loc["Clofazimin"]
    assert (emb.nr_susceptible, emb.nr_resistant, emb.nr_unclassifiable) == (1, 0, 1)
    assert (cfz.nr_susceptible, cfz.nr_resistant, cfz.nr_unclassifiable) == (1, 1, 0)
