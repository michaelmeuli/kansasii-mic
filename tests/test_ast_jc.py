from pathlib import Path

import pytest

from kansasii_mic.loading import load_ast_jc_long, snap_ast_jc_value
from kansasii_mic.parsing import MicParseError


@pytest.mark.parametrize(
    ("abbr", "raw", "expected"),
    [
        ("MXF", "0.1", "0.12"),
        ("MXF", "0.125", "0.12"),
        ("RFB", "0.13", "0.12"),
        ("CLR", "0.3", "0.25"),
        ("AMI", "2.0", "2"),
        ("AMI", "32.0", "32"),
        ("STR", "<64", "<64"),
        ("ETH", "0.2", "0.15"),
        ("ETH", "0.3", "0.3"),
        ("ETH", "1.2", "1.2"),
        ("SXT", ">0.12/2.49", ">0.12/2.38"),
        ("SXT", ">0.12/2.38", ">0.12/2.38"),
        ("SXT", "0.25/4.76", "0.25/4.75"),
        ("SXT", "0.5/9.5", "0.5/9.5"),
        ("SXT", "1.0/19", "1/19"),
    ],
)
def test_snap(abbr: str, raw: str, expected: str) -> None:
    assert snap_ast_jc_value(abbr, raw) == expected


def test_snap_rejects_unknown_sxt_numerator() -> None:
    with pytest.raises(MicParseError):
        snap_ast_jc_value("SXT", "3/57")


def _write(tmp_path: Path) -> tuple[Path, Path]:
    link = tmp_path / "link.csv"
    link.write_text(
        "NR,PROBENNUMMER,LNR,LNR2,TNR,TNR_JC,NGS_JC,TNR_MM,NGS_MM,TNR3,TNR4,TNR5,TNR6,TNR_MHK,MHK,MHK_TNR_DIFFERS,MIC_JC,LABEL\n"
        "1,P1,1,,100,,,,,,,,,,,,,\n"
        "2,P2,2,,200,201,,,,,,,,,,,yes,\n"
    )
    ast = tmp_path / "ast_jc.csv"
    ast.write_text(
        "TNR,AMI,CLR,EMB,INH,MXF ,RIF,SXT,CIP,DOX,ETH,LZD,RFB,STR\n"
        "100,1,0.12,2,2,0.1,0.1,>0.12/2.40,1,1,0.2,1,0.13,2\n"
        "201,2.0,0.3,4.0,4,0.12,0.25,0.25/4.76,2,4,0.3,2,0.25,<64\n"
        "999,1,1,1,1,1,1,1/19,1,1,0.3,1,0.25,1\n"
    )
    return ast, link


def test_load_ast_jc_long(tmp_path: Path) -> None:
    ast, link = _write(tmp_path)
    df, failures = load_ast_jc_long(ast, link, have_panel=set())
    assert failures.empty
    assert sorted(df["PROBENNUMMER"].unique()) == ["P1", "P2"]  # TNR 999 is out of scope
    assert len(df) == 26
    assert set(df["source"]) == {"ast_jc"}
    p2 = df[df["PROBENNUMMER"] == "P2"].set_index("antibiotic")
    assert p2.loc["Sulfamethoxazole/Trimethoprim", "mhk_parsed"] == "0.25/4.75"
    assert p2.loc["Amikacin", "mhk_parsed"] == "2"
    assert p2.loc["Streptomycin", "mhk_parsed"] == "<64"


def test_isolates_with_panel_data_are_skipped(tmp_path: Path) -> None:
    ast, link = _write(tmp_path)
    df, _ = load_ast_jc_long(ast, link, have_panel={"P1"})
    assert sorted(df["PROBENNUMMER"].unique()) == ["P2"]
