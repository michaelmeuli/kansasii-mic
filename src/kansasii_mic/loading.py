"""Load and normalize the mic.csv / screening_map.csv export."""

from __future__ import annotations

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
    """True for fixed-concentration breakpoint testing rows, e.g.

    "Amikacin 4 mg/l" -- these never carry an MHK value and are out of
    scope per the user's instruction to focus on MHK-valued data.
    """
    return "mg/l" in label


def load_mic_long(mic_csv: Path, screening_map_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (parsed_long_df, parse_failures_df).

    screening_map.csv is the source of truth for which isolates are in
    scope: a mic.csv row is kept only if its TNR matches screening_map's
    TNR column, or -- for a TNR recorded there as a superseded/duplicate
    id -- screening_map's MHK column. Rows matching neither are dropped.
    On an MHK-column match, the mic.csv TNR is replaced by screening_map's
    canonical TNR for that isolate.

    parsed_long_df columns: NR, PROBENNUMMER, TNR, MHK, raw_antibiotic,
    antibiotic, mhk_raw, point_estimate, log2_mic, denominator, censored,
    is_range, repaired_from, int_erg, erg, label.

    Note: this ``MHK`` (screening_map's secondary sample id) is unrelated
    to ``mhk_raw`` (the parsed lab MIC reading from mic.csv's own MHK
    column) -- the label is reused across the two source files for
    different things.
    """
    mic = pd.read_csv(mic_csv, dtype={"TNR": "Int64"})
    mic["MHK"] = mic["MHK"].fillna("").astype(str).str.strip()

    has_mhk = mic["MHK"] != ""
    subset = mic[has_mhk].copy()

    breakpoint_style = subset["ANTIBIOTIKA"].apply(_is_breakpoint_style_label)
    if breakpoint_style.any():
        # Sanity check on the assumption above: no such row is expected to
        # carry an actual MHK value in this dataset.
        subset = subset[~breakpoint_style]

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

    rows = []
    failures = []
    for row in subset.itertuples(index=False):
        if row.TNR in smap_by_tnr.index:
            meta = smap_by_tnr.loc[row.TNR]
        elif row.TNR in smap_by_mhk.index:
            meta = smap_by_mhk.loc[row.TNR]
        else:
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
