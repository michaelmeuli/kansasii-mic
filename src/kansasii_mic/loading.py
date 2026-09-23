"""Load and normalize the mic.csv / screening_map.csv export."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .parsing import MicParseError, parse_mic

# Raw ANTIBIOTIKA label -> canonical antibiotic name. Anything not listed
# here is used as-is. Two normalization decisions (see plan/README):
#  - "Sulfamethox.-Trimethop." and "Trimethoprim/Sulfamethoxazol" are the
#    same drug reported under two different labels over time.
#  - The three Amikacin labels (plain / i.v. / inhalativ) are treated as
#    the same underlying broth-microdilution test; the original label is
#    kept in `raw_antibiotic` for provenance.
CANONICAL_ANTIBIOTIC = {
    "Sulfamethox.-Trimethop.": "Sulfamethoxazole/Trimethoprim",
    "Trimethoprim/Sulfamethoxazol": "Sulfamethoxazole/Trimethoprim",
    "Amikacin i.v.": "Amikacin",
    "Amikacin inhalativ": "Amikacin",
}


def _is_breakpoint_style_label(label: str) -> bool:
    """True for fixed-concentration MGIT breakpoint testing rows, e.g.

    "Amikacin 4 mg/l". These report a single categorical growth call
    (INT_ERG/ERG: S/I/R/K/U) at one fixed tested concentration and always
    carry an *empty* ``MHK`` column in this dataset -- unlike the
    broth-microdilution panel rows, which carry a numeric/censored MIC
    reading there. Processed separately by :func:`load_mgit_long`.
    """
    return "mg/l" in label


class MgitLabelError(ValueError):
    """Raised when an MGIT ANTIBIOTIKA label cannot be split into drug + concentration."""


_MGIT_LABEL_RE = re.compile(r"^(?P<antibiotic>.+?)\s+(?P<concentration>[0-9]+(?:\.[0-9]+)?)\s*mg/l$")


def _parse_mgit_label(raw_antibiotic: str) -> tuple[str, float]:
    """Split an MGIT ANTIBIOTIKA label into (antibiotic name, tested concentration mg/L)."""
    m = _MGIT_LABEL_RE.match(raw_antibiotic.strip())
    if not m:
        raise MgitLabelError(f"unrecognized MGIT label: {raw_antibiotic!r}")
    return m.group("antibiotic"), float(m.group("concentration"))


def _lookup_screening_map(
    tnr: int, smap_by_tnr: pd.DataFrame, smap_by_mhk: pd.DataFrame
) -> pd.Series | None:
    """Match a mic.csv TNR to its screening_map metadata row, or None if out of scope.

    screening_map.csv is the source of truth for which isolates are in
    scope: a row is kept only if its TNR matches screening_map's TNR
    column, or -- for a TNR recorded there as a superseded/duplicate id --
    screening_map's MHK column. On an MHK-column match, the row's TNR is
    replaced by screening_map's canonical TNR for that isolate (the
    returned Series' ``TNR`` field).

    Note: this ``MHK`` (screening_map's secondary sample id) is unrelated
    to mic.csv's own ``MHK`` column (the parsed lab MIC reading) -- the
    label is reused across the two source files for different things.
    """
    if tnr in smap_by_tnr.index:
        return smap_by_tnr.loc[tnr]
    if tnr in smap_by_mhk.index:
        return smap_by_mhk.loc[tnr]
    return None


def _match_screening_map_and_parse(
    subset: pd.DataFrame, smap_by_tnr: pd.DataFrame, smap_by_mhk: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Match MHK broth-microdilution ``subset`` rows to screening_map metadata and parse MHK values.

    Returns (parsed_long_df, parse_failures_df). parsed_long_df columns:
    NR, PROBENNUMMER, TNR, MHK, raw_antibiotic, antibiotic, mhk_raw,
    mhk_parsed, point_estimate, log2_mic, denominator, censored, is_range,
    repaired_from, int_erg, erg, label.
    """
    rows = []
    failures = []
    for row in subset.itertuples(index=False):
        meta = _lookup_screening_map(row.TNR, smap_by_tnr, smap_by_mhk)
        if meta is None:
            continue

        try:
            parsed = parse_mic(row.MHK)
        except MicParseError as exc:
            failures.append({"TNR": meta["TNR"], "ANTIBIOTIKA": row.ANTIBIOTIKA, "MHK": row.MHK, "error": str(exc)})
            continue
        canonical = CANONICAL_ANTIBIOTIC.get(row.ANTIBIOTIKA, row.ANTIBIOTIKA)
        rows.append(
            {
                "NR": meta["NR"],
                "PROBENNUMMER": meta["PROBENNUMMER"],
                "TNR": meta["TNR"],
                "MHK": meta["MHK"],
                "raw_antibiotic": row.ANTIBIOTIKA,
                "antibiotic": canonical,
                "mhk_raw": row.MHK,
                "mhk_parsed": parsed.text,
                "point_estimate": parsed.point_estimate,
                "log2_mic": parsed.log2,
                "denominator": parsed.denominator,
                "censored": parsed.censored,
                "is_range": parsed.is_range,
                "repaired_from": parsed.repaired_from,
                "int_erg": row.INT_ERG,
                "erg": row.ERG,
                "label": meta["LABEL"],
            }
        )

    long_df = pd.DataFrame(rows)
    failures_df = pd.DataFrame(failures)
    return long_df, failures_df


def _match_screening_map_mgit(
    subset: pd.DataFrame, smap_by_tnr: pd.DataFrame, smap_by_mhk: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Match MGIT breakpoint ``subset`` rows to screening_map metadata and parse their labels.

    Same TNR/MHK matching as :func:`_match_screening_map_and_parse`, but
    these rows carry no numeric MHK reading -- the result is the
    categorical INT_ERG/ERG call at the concentration named in the label.

    Returns (parsed_long_df, parse_failures_df). parsed_long_df columns:
    NR, PROBENNUMMER, TNR, MHK, raw_antibiotic, antibiotic,
    concentration_mg_l, int_erg, erg, label.
    """
    rows = []
    failures = []
    for row in subset.itertuples(index=False):
        meta = _lookup_screening_map(row.TNR, smap_by_tnr, smap_by_mhk)
        if meta is None:
            continue

        try:
            antibiotic, concentration = _parse_mgit_label(row.ANTIBIOTIKA)
        except MgitLabelError as exc:
            failures.append({"TNR": meta["TNR"], "ANTIBIOTIKA": row.ANTIBIOTIKA, "error": str(exc)})
            continue

        rows.append(
            {
                "NR": meta["NR"],
                "PROBENNUMMER": meta["PROBENNUMMER"],
                "TNR": meta["TNR"],
                "MHK": meta["MHK"],
                "raw_antibiotic": row.ANTIBIOTIKA,
                "antibiotic": antibiotic,
                "concentration_mg_l": concentration,
                "int_erg": row.INT_ERG,
                "erg": row.ERG,
                "label": meta["LABEL"],
            }
        )

    long_df = pd.DataFrame(rows)
    failures_df = pd.DataFrame(failures)
    return long_df, failures_df


def _load_smap_indexes(screening_map_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    smap = pd.read_csv(
        screening_map_csv,
        dtype={"NR": "Int64", "PROBENNUMMER": "string", "TNR": "Int64", "MHK": "Int64"},
    )
    smap_cols = ["NR", "PROBENNUMMER", "TNR", "MHK", "LABEL"]
    smap_by_tnr = (
        smap.dropna(subset=["TNR"]).drop_duplicates(subset="TNR", keep="first").set_index("TNR", drop=False)[smap_cols]
    )
    smap_by_mhk = (
        smap.dropna(subset=["MHK"]).drop_duplicates(subset="MHK", keep="first").set_index("MHK", drop=False)[smap_cols]
    )
    return smap_by_tnr, smap_by_mhk


def load_mic_long(mic_csv: Path, screening_map_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (parsed_long_df, parse_failures_df) for the MHK broth-microdilution panel rows.

    Rows with a fixed-concentration MGIT label (e.g. "Amikacin 4 mg/l") are
    excluded here -- use :func:`load_mgit_long` for those.
    """
    mic = pd.read_csv(mic_csv, dtype={"TNR": "Int64"})
    mic["MHK"] = mic["MHK"].fillna("").astype(str).str.strip()

    has_mhk = mic["MHK"] != ""
    breakpoint_style = mic["ANTIBIOTIKA"].apply(_is_breakpoint_style_label)
    subset = mic[has_mhk & ~breakpoint_style].copy()

    smap_by_tnr, smap_by_mhk = _load_smap_indexes(screening_map_csv)
    return _match_screening_map_and_parse(subset, smap_by_tnr, smap_by_mhk)


def load_mgit_long(mic_csv: Path, screening_map_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (parsed_long_df, parse_failures_df) for the MGIT breakpoint-style rows.

    Same screening_map TNR/MHK matching as :func:`load_mic_long`, but
    selecting rows whose ANTIBIOTIKA label is a fixed-concentration MGIT
    reading (e.g. "Amikacin 4 mg/l") and reading the categorical
    INT_ERG/ERG growth call rather than parsing an MHK dilution value
    (these rows' MHK column is always empty in this dataset).
    """
    mic = pd.read_csv(mic_csv, dtype={"TNR": "Int64"})

    breakpoint_style = mic["ANTIBIOTIKA"].apply(_is_breakpoint_style_label)
    subset = mic[breakpoint_style].copy()

    smap_by_tnr, smap_by_mhk = _load_smap_indexes(screening_map_csv)
    return _match_screening_map_mgit(subset, smap_by_tnr, smap_by_mhk)
