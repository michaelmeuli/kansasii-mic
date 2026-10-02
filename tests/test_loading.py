import pandas as pd

from kansasii_mic.loading import load_mic_long


def _write(tmp_path):
    link = tmp_path / "link.csv"
    link.write_text(
        "NR,PROBENNUMMER,LNR,LNR2,TNR,TNR_NGS,TNR3,TNR4,TNR5,TNR6,NGS,MHK,LABEL\n"
        "1,P1,1,,100,,,,,,,,\n"
        "2,P2,2,,200,201,,,,,,,\n"
        "3,P3,3,,300,,,,,,,301,\n"
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


def test_matches_all_tnr_columns_and_mhk(tmp_path):
    mic, link = _write(tmp_path)
    df, _ = load_mic_long(mic, link)
    assert sorted(df["PROBENNUMMER"]) == ["P1", "P2", "P3"]
    assert sorted(df["TNR"]) == [100, 200, 300]
