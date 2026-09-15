"""Detect data-quality artifacts before running the outlier analysis.

The dataset contains at least one TNR (``2021311033``) where two whole
blocks of rows -- each covering many unrelated drug classes -- report the
exact same MIC value for every antibiotic in the block (all ``4``, then
separately all ``1``). That is biologically implausible (real MICs vary
across drug classes) and looks like placeholder/template rows rather than
real results. This is detected generically as a run of consecutive rows,
within one TNR and in original file order, spanning many distinct
antibiotics that all share one point estimate -- not by hardcoding the
TNR value.
"""

from __future__ import annotations

import pandas as pd

MIN_FLAT_RUN_LENGTH = 8  # antibiotics; chosen well above what plausibly
# ties by chance (few dilution levels x few antibiotics), well below the
# ~13-17 distinct drugs seen in the actual artifact runs.


def find_flat_profile_runs(long_df: pd.DataFrame) -> pd.DataFrame:
    """Return the subset of rows belonging to a flat-profile artifact run.

    Includes an ``exclusion_reason`` column explaining the flagged run.
    """
    flagged_parts = []
    for tnr, group in long_df.groupby("TNR", sort=False):
        values = group["point_estimate"].tolist()
        idx = group.index.tolist()
        run_start = 0
        n = len(values)
        i = 1
        while i <= n:
            same_as_start = i < n and values[i] == values[run_start]
            if not same_as_start:
                run_len = i - run_start
                if run_len >= MIN_FLAT_RUN_LENGTH:
                    run_rows = group.loc[idx[run_start:i]].copy()
                    run_rows["exclusion_reason"] = (
                        f"flat-profile artifact: {run_len} distinct antibiotics all "
                        f"report identical MIC point estimate {values[run_start]!r} "
                        f"for TNR {tnr}"
                    )
                    flagged_parts.append(run_rows)
                run_start = i
            i += 1

    if not flagged_parts:
        return long_df.iloc[0:0].assign(exclusion_reason=pd.Series(dtype=str))
    return pd.concat(flagged_parts, axis=0)


def apply_qc(long_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (clean_df, excluded_df) after removing flat-profile artifacts."""
    excluded = find_flat_profile_runs(long_df)
    clean = long_df.drop(index=excluded.index)
    return clean, excluded
