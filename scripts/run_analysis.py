#!/usr/bin/env python3
"""CLI entrypoint: load MIC data, run QC + outlier analysis, write outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from kansasii_mic.loading import load_mgit_long, load_mic_long
from kansasii_mic.mgit import filter_int_erg_i_or_r, mgit_summary, tnr_mgit_resistance_ranking
from kansasii_mic.outliers import per_antibiotic_outliers, tnr_resistance_ranking
from kansasii_mic.plotting import (
    plot_antibiotic_distributions,
    plot_heatmap,
    plot_mgit_category_counts,
    plot_mgit_heatmap,
    plot_mgit_resistance_ranking,
    plot_resistance_ranking,
)
from kansasii_mic.qc import apply_qc


def _run_pipeline(long_df, failures_df, out_dir: Path, label: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = out_dir / "figures"

    if not failures_df.empty:
        failures_df.to_csv(out_dir / "parse_failures.csv", index=False)
        print(f"WARNING [{label}]: {len(failures_df)} MHK values could not be parsed; see parse_failures.csv")

    clean_df, excluded_df = apply_qc(long_df)
    clean_df.to_csv(out_dir / "mic_parsed.csv", index=False)
    excluded_df.to_csv(out_dir / "excluded_qc.csv", index=False)

    outliers_df = per_antibiotic_outliers(clean_df)
    outliers_df.to_csv(out_dir / "outliers_per_antibiotic.csv", index=False)

    ranking_df = tnr_resistance_ranking(outliers_df)
    ranking_df.to_csv(out_dir / "tnr_resistance_ranking.csv", index=False)

    dist_paths = plot_antibiotic_distributions(outliers_df, figures_dir)
    heatmap_path = plot_heatmap(outliers_df, ranking_df, figures_dir)
    ranking_path = plot_resistance_ranking(ranking_df, figures_dir)

    print(f"[{label}] Parsed rows (post-QC): {len(clean_df)}")
    print(f"[{label}] Unique TNRs analyzed: {clean_df['TNR'].nunique()}")
    print(f"[{label}] Excluded QC rows: {len(excluded_df)} (TNRs: {sorted(excluded_df['TNR'].unique().tolist())})")
    print(f"[{label}] Outlier flags: {(outliers_df['outlier_direction'].notna()).sum()}")
    print(f"[{label}] Figures written: {len(dist_paths)} per-antibiotic + {heatmap_path.name} + {ranking_path.name}")
    print(f"[{label}] Output directory: {out_dir}")


def _run_mgit_pipeline(mgit_long_df, mgit_failures_df, out_dir: Path) -> None:
    """MGIT breakpoint rows report a categorical S/I/R/K/U call, not a continuous
    MIC, so this mirrors _run_pipeline's shape (parsed table, summary tables,
    figures) without the continuous-value QC/outlier statistics.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = out_dir / "figures"

    if not mgit_failures_df.empty:
        mgit_failures_df.to_csv(out_dir / "parse_failures.csv", index=False)
        print(f"WARNING [MGIT]: {len(mgit_failures_df)} ANTIBIOTIKA labels could not be parsed; see parse_failures.csv")

    mgit_long_df.to_csv(out_dir / "mgit_parsed.csv", index=False)

    ir_highest_df = filter_int_erg_i_or_r(mgit_long_df)
    ir_highest_df.to_csv(out_dir / "mgit_int_erg_I_R_highest_concentration.csv", index=False)

    summary_df = mgit_summary(mgit_long_df)
    summary_df.to_csv(out_dir / "mgit_summary_per_antibiotic.csv", index=False)

    ranking_df = tnr_mgit_resistance_ranking(mgit_long_df)
    ranking_df.to_csv(out_dir / "tnr_mgit_resistance_ranking.csv", index=False)

    count_paths = plot_mgit_category_counts(summary_df, figures_dir)
    heatmap_path = plot_mgit_heatmap(mgit_long_df, ranking_df, figures_dir)
    ranking_path = plot_mgit_resistance_ranking(ranking_df, figures_dir)

    print(f"[MGIT] Parsed rows: {len(mgit_long_df)}")
    print(f"[MGIT] Unique TNRs analyzed: {mgit_long_df['TNR'].nunique()}")
    print(f"[MGIT] I/R rows (R collapsed to highest concentration): {len(ir_highest_df)}")
    print(f"[MGIT] Figures written: {len(count_paths)} per-antibiotic + {heatmap_path.name} + {ranking_path.name}")
    print(f"[MGIT] Output directory: {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mic-csv", type=Path, required=True)
    parser.add_argument("--screening-map", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    long_df, failures_df = load_mic_long(args.mic_csv, args.screening_map)
    _run_pipeline(long_df, failures_df, args.out_dir / "mhk", "MHK")

    mgit_long_df, mgit_failures_df = load_mgit_long(args.mic_csv, args.screening_map)
    _run_mgit_pipeline(mgit_long_df, mgit_failures_df, args.out_dir / "mgit")


if __name__ == "__main__":
    main()
