"""Summaries for MGIT fixed-concentration breakpoint results.

Unlike the MHK broth-microdilution panel, MGIT rows report a single
categorical growth call (ERG: S/I/R/K/U) at one fixed tested concentration
per antibiotic, not a continuous MIC value. The Tukey-fence/robust-z
outlier statistics in :mod:`kansasii_mic.outliers` are built for continuous
log2 MIC values and don't have a clean meaning here, so this module builds
counts and rankings directly off the categorical calls instead.
"""

from __future__ import annotations

import pandas as pd

# K = no assessment criteria defined, U = result unclear/indeterminate.
ERG_CATEGORIES = ["S", "I", "R", "K", "U"]


def mgit_summary(mgit_long_df: pd.DataFrame) -> pd.DataFrame:
    """Per antibiotic x tested concentration: counts of each ERG category."""
    counts = (
        mgit_long_df.groupby(["antibiotic", "concentration_mg_l"])["erg"]
        .value_counts()
        .unstack("erg", fill_value=0)
    )
    for cat in ERG_CATEGORIES:
        if cat not in counts.columns:
            counts[cat] = 0
    counts = counts[ERG_CATEGORIES]
    counts["n_tested"] = counts[ERG_CATEGORIES].sum(axis=1)
    counts["pct_resistant"] = (counts["R"] / counts["n_tested"] * 100).round(1)
    return counts.reset_index().sort_values(["antibiotic", "concentration_mg_l"]).reset_index(drop=True)


def tnr_mgit_resistance_ranking(mgit_long_df: pd.DataFrame) -> pd.DataFrame:
    """Per-TNR count of each ERG category across all antibiotics/concentrations tested.

    Ranked by number of R (resistant) calls, most first.
    """

    def _agg(group: pd.DataFrame) -> pd.Series:
        counts = group["erg"].value_counts()
        return pd.Series(
            {
                "PROBENNUMMER": group["PROBENNUMMER"].iloc[0],
                "n_tests": len(group),
                **{f"n_{cat}": int(counts.get(cat, 0)) for cat in ERG_CATEGORIES},
            }
        )

    ranking = mgit_long_df.groupby("TNR").apply(_agg, include_groups=False).reset_index()
    ranking = ranking.sort_values("n_R", ascending=False).reset_index(drop=True)
    ranking["rank_most_resistant"] = ranking.index + 1
    return ranking


def filter_int_erg_i_or_r(mgit_long_df: pd.DataFrame) -> pd.DataFrame:
    """Keep only int_erg I/R rows, collapsing R calls to the highest concentration tested.

    A drug is typically tested at several concentrations per isolate; an R
    at a higher concentration is the more definitive result, so when the
    same isolate x antibiotic has R at more than one concentration, only
    the highest-concentration R row is kept. I rows are kept as-is (I is
    normally reported at a single concentration in this dataset).
    """
    subset = mgit_long_df[mgit_long_df["int_erg"].isin(["I", "R"])]
    is_r = subset["int_erg"] == "R"
    r_highest = (
        subset[is_r]
        .sort_values("concentration_mg_l")
        .drop_duplicates(subset=["TNR", "antibiotic"], keep="last")
    )
    result = pd.concat([subset[~is_r], r_highest])
    return result.sort_values(["TNR", "antibiotic", "concentration_mg_l"]).reset_index(drop=True)
