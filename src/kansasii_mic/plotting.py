"""Static matplotlib figures summarizing the MIC outlier analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .breakpoints import CLSI_KANSASII_BREAKPOINTS
from .mgit import ERG_CATEGORIES

POINT_COLOR = "#3B6FA0"
OUTLIER_RESISTANT_COLOR = "#C0392B"
OUTLIER_SUSCEPTIBLE_COLOR = "#1E8449"
BREAKPOINT_S_COLOR = "#2E7D32"
BREAKPOINT_R_COLOR = "#B71C1C"
HEATMAP_CMAP = "YlOrRd"

# S/I/R/K/U category colors for the MGIT breakpoint figures.
ERG_COLORS = {
    "S": BREAKPOINT_S_COLOR,
    "I": "#F9A825",
    "R": BREAKPOINT_R_COLOR,
    "K": "#9E9E9E",
    "U": "#607D8B",
}

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
    """One strip plot per antibiotic: log2 MIC across isolates, CLSI lines, outliers labeled by PROBENNUMMER."""
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
            ax.annotate(str(row["PROBENNUMMER"]), (row["log2_mic"], y), textcoords="offset points", xytext=(4, 4), fontsize=7, color="#333333")

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
    pivot = outliers_df.pivot_table(index="PROBENNUMMER", columns="antibiotic", values="log2_mic", aggfunc="mean")
    isolate_order = [p for p in ranking_df["PROBENNUMMER"] if p in pivot.index]
    pivot = pivot.loc[isolate_order]

    fig, ax = plt.subplots(figsize=(max(6, 0.45 * pivot.shape[1]), max(5, 0.22 * pivot.shape[0])))
    masked = np.ma.masked_invalid(pivot.values)
    cmap = matplotlib.colormaps[HEATMAP_CMAP].copy()
    cmap.set_bad("#EEEEEE")
    im = ax.imshow(masked, aspect="auto", cmap=cmap)
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=6)
    ax.set_title("MIC (log2 mg/L) by isolate (PROBENNUMMER) x antibiotic\nrows sorted from most- to least-resistant overall")
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
    ax.barh(combined["PROBENNUMMER"].astype(str), combined["mean_robust_z"], color=colors)
    ax.axvline(0, color="#444444", linewidth=0.8)
    ax.set_xlabel("Mean robust z-score across tested antibiotics\n(> 0 = more resistant than cohort, < 0 = more susceptible)")
    ax.set_title(f"Most resistant / most susceptible isolates (top {top_n} each)")
    fig.tight_layout()
    path = out_dir / "resistance_ranking.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_mgit_category_counts(summary_df: pd.DataFrame, out_dir: Path) -> list[Path]:
    """One stacked bar chart per antibiotic: ERG category counts by tested concentration."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for antibiotic, group in summary_df.groupby("antibiotic"):
        group = group.sort_values("concentration_mg_l")
        fig, ax = plt.subplots(figsize=(6, 3.2))
        x = np.arange(len(group))
        bottom = np.zeros(len(group))
        for cat in ERG_CATEGORIES:
            values = group[cat].to_numpy(dtype=float)
            ax.bar(x, values, bottom=bottom, color=ERG_COLORS[cat], label=cat, width=0.6)
            bottom += values
        ax.set_xticks(x)
        ax.set_xticklabels([f"{c:g} mg/l" for c in group["concentration_mg_l"]])
        ax.set_ylabel("Isolates (n)")
        ax.set_title(f"{antibiotic} -- MGIT breakpoint results (n={int(group['n_tested'].sum())})")
        ax.legend(loc="upper right", fontsize=7, frameon=False, ncols=len(ERG_CATEGORIES))
        fig.tight_layout()
        path = out_dir / f"{_safe_filename(antibiotic)}_mgit_counts.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_mgit_heatmap(mgit_long_df: pd.DataFrame, ranking_df: pd.DataFrame, out_dir: Path) -> Path:
    """Isolate x antibiotic@concentration heatmap of ERG category (categorical, not a MIC scale)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    df = mgit_long_df.copy()
    df["column"] = df["antibiotic"] + " " + df["concentration_mg_l"].map(lambda v: f"{v:g}") + " mg/l"
    erg_code = {cat: i for i, cat in enumerate(ERG_CATEGORIES)}
    df["erg_code"] = df["erg"].map(erg_code)
    pivot = df.pivot_table(index="PROBENNUMMER", columns="column", values="erg_code", aggfunc="mean")
    isolate_order = [p for p in ranking_df["PROBENNUMMER"] if p in pivot.index]
    pivot = pivot.loc[isolate_order]

    fig, ax = plt.subplots(figsize=(max(6, 0.45 * pivot.shape[1]), max(5, 0.22 * pivot.shape[0])))
    cmap = matplotlib.colors.ListedColormap([ERG_COLORS[c] for c in ERG_CATEGORIES])
    cmap.set_bad("#EEEEEE")
    masked = np.ma.masked_invalid(pivot.values)
    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-0.5, vmax=len(ERG_CATEGORIES) - 0.5)
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=6)
    ax.set_title(
        "MGIT breakpoint result by isolate (PROBENNUMMER) x antibiotic@concentration\n"
        "rows sorted from most- to least-resistant overall"
    )
    cbar = fig.colorbar(im, ax=ax, ticks=range(len(ERG_CATEGORIES)), shrink=0.6)
    cbar.ax.set_yticklabels(ERG_CATEGORIES)
    cbar.set_label("ERG category")
    fig.tight_layout()
    path = out_dir / "heatmap_tnr_antibiotic_mgit.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_mgit_resistance_ranking(ranking_df: pd.DataFrame, out_dir: Path, top_n: int = 15) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    ranked = ranking_df.sort_values("n_R", ascending=False)
    top = ranked.head(top_n)
    bottom = ranked.tail(top_n)
    combined = pd.concat([top, bottom]).drop_duplicates(subset="TNR").sort_values("n_R")

    fig, ax = plt.subplots(figsize=(7, max(4, 0.32 * len(combined))))
    ax.barh(combined["PROBENNUMMER"].astype(str), combined["n_R"], color=BREAKPOINT_R_COLOR)
    ax.set_xlabel("Number of R (resistant) MGIT breakpoint calls")
    ax.set_title(f"Most MGIT-resistant isolates by R-call count (top {top_n})")
    fig.tight_layout()
    path = out_dir / "resistance_ranking_mgit.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
