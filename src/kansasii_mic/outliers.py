"""Per-antibiotic outlier detection and per-TNR resistance ranking."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .breakpoints import categorize

MAD_SCALE = 1.4826  # consistency constant so MAD approximates SD for normal data


def _tukey_fences(values: pd.Series) -> tuple[float, float, float, float]:
    q1, q3 = values.quantile([0.25, 0.75])
    iqr = q3 - q1
    return q1, q3, q1 - 1.5 * iqr, q3 + 1.5 * iqr


def per_antibiotic_outliers(clean_df: pd.DataFrame) -> pd.DataFrame:
    """Flag TNRs whose log2 MIC is a Tukey-fence outlier, per antibiotic.

    Also attaches the CLSI category (S/I/R/None) for cross-reference.
    """
    results = []
    for antibiotic, group in clean_df.groupby("antibiotic"):
        q1, q3, lower_fence, upper_fence = _tukey_fences(group["log2_mic"])
        median = group["log2_mic"].median()
        mad = (group["log2_mic"] - median).abs().median() * MAD_SCALE
        for row in group.itertuples(index=False):
            direction = None
            if row.log2_mic < lower_fence:
                direction = "more_susceptible"
            elif row.log2_mic > upper_fence:
                direction = "more_resistant"
            robust_z = (row.log2_mic - median) / mad if mad > 0 else np.nan
            results.append(
                {
                    "TNR": row.TNR,
                    "antibiotic": antibiotic,
                    "mhk_raw": row.mhk_raw,
                    "point_estimate": row.point_estimate,
                    "log2_mic": row.log2_mic,
                    "cohort_median_log2": median,
                    "robust_z": robust_z,
                    "tukey_lower_fence": lower_fence,
                    "tukey_upper_fence": upper_fence,
                    "outlier_direction": direction,
                    "clsi_category": categorize(antibiotic, row.point_estimate),
                    "lab_erg": row.erg,
                }
            )
    return pd.DataFrame(results)


def tnr_resistance_ranking(outliers_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-TNR score across all antibiotics that TNR was tested for.

    score = mean robust z-score across drugs (higher = more resistant
    overall), plus a count of CLSI-defined resistant (R) results.
    """
    def _agg(group: pd.DataFrame) -> pd.Series:
        return pd.Series(
            {
                "n_antibiotics_tested": len(group),
                "mean_robust_z": group["robust_z"].mean(skipna=True),
                "n_clsi_resistant": (group["clsi_category"] == "R").sum(),
                "n_clsi_susceptible": (group["clsi_category"] == "S").sum(),
                "n_outlier_more_resistant": (group["outlier_direction"] == "more_resistant").sum(),
                "n_outlier_more_susceptible": (group["outlier_direction"] == "more_susceptible").sum(),
            }
        )

    ranking = outliers_df.groupby("TNR").apply(_agg, include_groups=False).reset_index()
    ranking = ranking.sort_values("mean_robust_z", ascending=False).reset_index(drop=True)
    ranking["rank_most_resistant"] = ranking.index + 1
    return ranking
