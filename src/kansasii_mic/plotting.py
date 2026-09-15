"""Static matplotlib figures summarizing the MIC outlier analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .breakpoints import CLSI_KANSASII_BREAKPOINTS

POINT_COLOR = "#3B6FA0"
OUTLIER_RESISTANT_COLOR = "#C0392B"
OUTLIER_SUSCEPTIBLE_COLOR = "#1E8449"
BREAKPOINT_S_COLOR = "#2E7D32"
BREAKPOINT_R_COLOR = "#B71C1C"
HEATMAP_CMAP = "YlOrRd"

plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "#444444",
        "axes.grid": True,
        "grid.color": "#DDDDDD",
        "grid.linewidth": 0.6,
        "font.size": 10,
    }
)


def _safe_filename(name: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def plot_antibiotic_distributions(outliers_df: pd.DataFrame, out_dir: Path) -> list[Path]:
    """One strip plot per antibiotic: log2 MIC across TNRs, CLSI lines, outliers labeled."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for antibiotic, group in outliers_df.groupby("antibiotic"):
        group = group.sort_values("log2_mic")
        fig, ax = plt.subplots(figsize=(7, 3.2))
        rng = np.random.default_rng(abs(hash(antibiotic)) % (2**32))
        jitter = rng.uniform(-0.15, 0.15, size=len(group))
        colors = [
            OUTLIER_RESISTANT_COLOR
            if d == "more_resistant"
            else OUTLIER_SUSCEPTIBLE_COLOR
            if d == "more_susceptible"
            else POINT_COLOR
            for d in group["outlier_direction"]
        ]
        ax.scatter(group["log2_mic"], jitter, c=colors, s=36, alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)

        bp = CLSI_KANSASII_BREAKPOINTS.get(antibiotic)
        if bp is not None:
            if bp.susceptible_max is not None:
                ax.axvline(np.log2(bp.susceptible_max), color=BREAKPOINT_S_COLOR, linestyle="--", linewidth=1.2, label=f"CLSI S ≤ {bp.susceptible_max}")
            if bp.resistant_min is not None:
                ax.axvline(np.log2(bp.resistant_min), color=BREAKPOINT_R_COLOR, linestyle="--", linewidth=1.2, label=f"CLSI R ≥ {bp.resistant_min}")

        for _, row in group[group["outlier_direction"].notna()].iterrows():
            y = jitter[group.index.get_loc(row.name)]
            ax.annotate(str(row["TNR"]), (row["log2_mic"], y), textcoords="offset points", xytext=(4, 4), fontsize=7, color="#333333")

        ax.set_yticks([])
        ax.set_ylim(-0.5, 0.5)
        xticks = sorted(group["log2_mic"].unique())
        ax.set_xticks(xticks)
        ax.set_xticklabels([_fmt_mic(v) for v in xticks], rotation=45, ha="right", fontsize=8)
        ax.set_xlabel("MIC (mg/L, log2 scale)")
        ax.set_title(f"{antibiotic} -- MIC distribution across isolates (n={len(group)})")
        if bp is not None:
            ax.legend(loc="upper left", fontsize=7, frameon=False)
        fig.tight_layout()
        path = out_dir / f"{_safe_filename(antibiotic)}_distribution.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths.append(path)
    return paths


def _fmt_mic(log2_value: float) -> str:
    value = 2**log2_value
    if value == int(value):
        return str(int(value))
    return f"{value:g}"


def plot_heatmap(outliers_df: pd.DataFrame, ranking_df: pd.DataFrame, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    pivot = outliers_df.pivot_table(index="TNR", columns="antibiotic", values="log2_mic", aggfunc="mean")
    tnr_order = [t for t in ranking_df["TNR"] if t in pivot.index]
    pivot = pivot.loc[tnr_order]

    fig, ax = plt.subplots(figsize=(max(6, 0.45 * pivot.shape[1]), max(5, 0.22 * pivot.shape[0])))
    masked = np.ma.masked_invalid(pivot.values)
    cmap = matplotlib.colormaps[HEATMAP_CMAP].copy()
    cmap.set_bad("#EEEEEE")
    im = ax.imshow(masked, aspect="auto", cmap=cmap)
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=6)
    ax.set_title("MIC (log2 mg/L) by isolate (TNR) x antibiotic\nrows sorted from most- to least-resistant overall")
    cbar = fig.colorbar(im, ax=ax, shrink=0.6)
    cbar.set_label("log2 MIC (mg/L)")
    fig.tight_layout()
    path = out_dir / "heatmap_tnr_antibiotic.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_resistance_ranking(ranking_df: pd.DataFrame, out_dir: Path, top_n: int = 15) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    ranked = ranking_df.dropna(subset=["mean_robust_z"]).sort_values("mean_robust_z", ascending=False)
    top = ranked.head(top_n)
    bottom = ranked.tail(top_n)
    combined = pd.concat([top, bottom]).drop_duplicates(subset="TNR").sort_values("mean_robust_z")

    fig, ax = plt.subplots(figsize=(7, max(4, 0.32 * len(combined))))
    colors = [OUTLIER_RESISTANT_COLOR if v > 0 else OUTLIER_SUSCEPTIBLE_COLOR for v in combined["mean_robust_z"]]
    ax.barh(combined["TNR"].astype(str), combined["mean_robust_z"], color=colors)
    ax.axvline(0, color="#444444", linewidth=0.8)
    ax.set_xlabel("Mean robust z-score across tested antibiotics\n(> 0 = more resistant than cohort, < 0 = more susceptible)")
    ax.set_title(f"Most resistant / most susceptible isolates (top {top_n} each)")
    fig.tight_layout()
    path = out_dir / "resistance_ranking.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
