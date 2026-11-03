from pathlib import Path

import pandas as pd

from kansasii_mic.loading import load_mic_long


def _write(tmp_path: Path) -> tuple[Path, Path]:
    link = tmp_path / "link.csv"
    link.write_text(
        "NR,PROBENNUMMER,LNR,LNR2,TNR,TNR_JC,NGS_JC,TNR_MM,NGS_MM,TNR3,TNR4,TNR5,TNR6,MHK,LABEL\n"
        "1,P1,1,,100,,,,,,,,,,\n"
        "2,P2,2,,200,201,,,,,,,,,\n"
        "3,P3,3,,300,,,,,,,,,301,\n"
    )
    mic = tmp_path / "mic.csv"
    mic.write_text(
        "TNR,ANTIBIOTIKA,MHK,INT_ERG,ERG\n"
        "100,Amikacin,2,S,S\n"
        "201,Amikacin,4,S,S\n"
        "301,Amikacin,8,S,S\n"
        "999,Amikacin,1,S,S\n"
    )
    return mic, link


def test_matches_all_tnr_columns_and_mhk(tmp_path: Path) -> None:
    mic, link = _write(tmp_path)
    df, _ = load_mic_long(mic, link)
    assert sorted(df["PROBENNUMMER"]) == ["P1", "P2", "P3"]
    assert sorted(df["TNR"]) == [100, 200, 300]


def test_mhk_made_under_other_tnr_is_still_matched(tmp_path: Path) -> None:
    mic, link = _write(tmp_path)
    link.write_text(
        "NR,PROBENNUMMER,LNR,LNR2,TNR,TNR_JC,NGS_JC,TNR_MM,NGS_MM,TNR3,TNR4,TNR5,TNR6,TNR_MHK,MHK,MHK_TNR_DIFFERS,LABEL\n"
        "1,P1,1,,100,,,,,,,,,100,,,\n"
        "3,P3,3,,300,,,,,,,,,350,301,yes,\n"
    )
    df, _ = load_mic_long(mic, link)
    assert sorted(df["PROBENNUMMER"]) == ["P1", "P3"]


def test_drop_excluded(tmp_path: Path) -> None:
    import pandas as pd

    from kansasii_mic.loading import drop_excluded

    res = tmp_path / "res.csv"
    res.write_text("PROBENNUMMER,excluded\nA,true\nB,false\nC,true\n")
    df = pd.DataFrame({"PROBENNUMMER": ["A", "B", "B"], "x": [1, 2, 3]})
    out, dropped = drop_excluded(df, res)
    assert list(out["PROBENNUMMER"]) == ["B", "B"]
    assert dropped == ["A"]
    out, dropped = drop_excluded(df, tmp_path / "missing.csv")
    assert len(out) == 3 and dropped == []


def test_drop_unclear_when_result_exists() -> None:
    from kansasii_mic.loading import drop_unclear_when_result_exists

    df = pd.DataFrame(
        {
            "NR": [1, 1, 1, 2, 3, 3, 4, 4, 5, 5],
            "antibiotic": ["Isoniazid"] * 10,
            "concentration_mg_l": [10.0, 10.0, 3.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0],
            "erg": ["S", "U", "U", "U", "U", "U", "R", "S", "S", "S"],
        }
    )
    out = drop_unclear_when_result_exists(df)
    assert list(zip(out["NR"], out["concentration_mg_l"], out["erg"])) == [
        (1, 10.0, "S"),
        (1, 3.0, "U"),  # only U for this key -> kept
        (2, 10.0, "U"),
        (3, 10.0, "U"),  # only U rows for this key: one kept
        (4, 10.0, "S"),  # R and S: the more resistant one is dropped
        (5, 10.0, "S"),  # identical duplicates collapse to one
    ]
