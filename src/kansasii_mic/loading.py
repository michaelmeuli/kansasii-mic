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

    parsed_long_df columns: TNR, raw_antibiotic, antibiotic, mhk_raw,
    point_estimate, log2_mic, denominator, censored, is_range,
    repaired_from, int_erg, erg, probennummer, label.
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

    smap = pd.read_csv(screening_map_csv, dtype={"TNR": "Int64"})
    smap_by_tnr = (
        smap.dropna(subset=["TNR"])
        .drop_duplicates(subset="TNR", keep="first")
        .set_index("TNR")[["PROBENNUMMER", "LABEL"]]
    )

    rows = []
    failures = []
    for row in subset.itertuples(index=False):
        try:
            parsed = parse_mic(row.MHK)
        except MicParseError as exc:
            failures.append({"TNR": row.TNR, "ANTIBIOTIKA": row.ANTIBIOTIKA, "MHK": row.MHK, "error": str(exc)})
            continue
        canonical = CANONICAL_ANTIBIOTIC.get(row.ANTIBIOTIKA, row.ANTIBIOTIKA)
        meta = smap_by_tnr.loc[row.TNR] if row.TNR in smap_by_tnr.index else None
        rows.append(
            {
                "TNR": row.TNR,
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
                "probennummer": meta["PROBENNUMMER"] if meta is not None else None,
                "label": meta["LABEL"] if meta is not None else None,
            }
        )

    long_df = pd.DataFrame(rows)
    failures_df = pd.DataFrame(failures)
    return long_df, failures_df
