#!/usr/bin/env python3
"""CLI entrypoint: load MIC data, run QC + outlier analysis, write outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from kansasii_mic.loading import load_mic_long
from kansasii_mic.outliers import per_antibiotic_outliers, tnr_resistance_ranking
from kansasii_mic.plotting import plot_antibiotic_distributions, plot_heatmap, plot_resistance_ranking
from kansasii_mic.qc import apply_qc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mic-csv", type=Path, required=True)
    parser.add_argument("--screening-map", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = args.out_dir / "figures"

    long_df, failures_df = load_mic_long(args.mic_csv, args.screening_map)
    if not failures_df.empty:
        failures_df.to_csv(args.out_dir / "parse_failures.csv", index=False)
        print(f"WARNING: {len(failures_df)} MHK values could not be parsed; see parse_failures.csv")

    clean_df, excluded_df = apply_qc(long_df)
    clean_df.to_csv(args.out_dir / "mic_parsed.csv", index=False)
    excluded_df.to_csv(args.out_dir / "excluded_qc.csv", index=False)

    outliers_df = per_antibiotic_outliers(clean_df)
    outliers_df.to_csv(args.out_dir / "outliers_per_antibiotic.csv", index=False)

    ranking_df = tnr_resistance_ranking(outliers_df)
    ranking_df.to_csv(args.out_dir / "tnr_resistance_ranking.csv", index=False)

    dist_paths = plot_antibiotic_distributions(outliers_df, figures_dir)
    heatmap_path = plot_heatmap(outliers_df, ranking_df, figures_dir)
    ranking_path = plot_resistance_ranking(ranking_df, figures_dir)

    print(f"Parsed rows (MHK-valued, post-QC): {len(clean_df)}")
    print(f"Unique TNRs analyzed: {clean_df['TNR'].nunique()}")
    print(f"Excluded QC rows: {len(excluded_df)} (TNRs: {sorted(excluded_df['TNR'].unique().tolist())})")
    print(f"Outlier flags: {(outliers_df['outlier_direction'].notna()).sum()}")
    print(f"Figures written: {len(dist_paths)} per-antibiotic + {heatmap_path.name} + {ranking_path.name}")
    print(f"Output directory: {args.out_dir}")


if __name__ == "__main__":
    main()
