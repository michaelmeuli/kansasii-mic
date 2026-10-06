"""Static matplotlib figures summarizing the MIC outlier analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

from .breakpoints import CLSI_KANSASII_BREAKPOINTS, clsi_normalized
from .mgit import ERG_CATEGORIES

POINT_COLOR = "#3B6FA0"
OUTLIER_RESISTANT_COLOR = "#C0392B"
OUTLIER_SUSCEPTIBLE_COLOR = "#1E8449"
BREAKPOINT_S_COLOR = "#2E7D32"
BREAKPOINT_R_COLOR = "#B71C1C"
HEATMAP_CMAP = "YlOrRd"
# Iglewicz-Hoaglin cutoff: |modified_z| above this is labeled on the distribution plots.
LABEL_Z_THRESHOLD = 1.0

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
    """One strip plot per antibiotic: log2 MIC across isolates, CLSI lines, points with |modified_z| > LABEL_Z_THRESHOLD labeled by NR."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for antibiotic_key, group in outliers_df.groupby("antibiotic"):
        antibiotic = str(antibiotic_key)
        group = group.sort_values(["log2_mic", "NR"])
        fig, ax = plt.subplots(figsize=(7, 4.2))
        jitter = _tie_spread(group["log2_mic"]).to_numpy()
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

        label_x, label_y, texts = [], [], []
        for _, row in group[group["modified_z"].abs() > LABEL_Z_THRESHOLD].iterrows():
            y = float(jitter[group.index.get_loc(row.name)])
            label_x.append(row["log2_mic"])
            label_y.append(y)
            texts.append(ax.text(row["log2_mic"], y, str(int(row["NR"])), fontsize=7, color="#333333"))

        ax.set_yticks([])
        ax.set_ylim(-0.5, 0.5)
        xticks = sorted(group["log2_mic"].unique())
        ax.set_xticks(xticks)
        ax.set_xticklabels([_fmt_mic(v) for v in xticks], rotation=45, ha="right", fontsize=8)
        ax.set_xlabel("MIC (mg/L, log2 scale)")
        ax.set_title(f"{antibiotic} -- MIC distribution across isolates (n={len(group)})")
        if texts:
            # Needs final axis limits/ticks; spreads labels apart with thin leader lines.
            adjust_text(texts, x=label_x, y=label_y, ax=ax, arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5))
        if bp is not None:
            ax.legend(loc="upper left", fontsize=7, frameon=False)
        fig.tight_layout()
        path = out_dir / f"{_safe_filename(antibiotic)}_distribution.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths.append(path)
    return paths


def _tie_spread(log2_mic: pd.Series, step: float = 0.05, max_span: float = 0.8) -> pd.Series:
    """Deterministic y position: isolates sharing a MIC are evenly spaced around y=0, in row order."""
    by_value = log2_mic.groupby(log2_mic)
    rank = by_value.cumcount()
    n = by_value.transform("size")
    step_per_group = np.minimum(step, max_span / (n - 1).clip(lower=1))
    return (rank - (n - 1) / 2) * step_per_group


def _fmt_mic(log2_value: float) -> str:
    value = 2**log2_value
    if value == int(value):
        return str(int(value))
    return f"{value:g}"


def plot_heatmap(outliers_df: pd.DataFrame, ranking_df: pd.DataFrame, out_dir: Path) -> Path:
    """Isolate x antibiotic heatmap, normalized per drug to CLSI (S limit = 0, R limit = 1).

    Drugs without a CLSI breakpoint pair go in a second panel, colored by cohort modified_z.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    df = outliers_df.copy()
    df["clsi_norm"] = [clsi_normalized(a, v) for a, v in zip(df["antibiotic"], df["log2_mic"])]
    has_bp = df.groupby("antibiotic")["clsi_norm"].apply(lambda s: s.notna().any())
    bp_drugs = sorted(has_bp[has_bp].index)
    other_drugs = sorted(has_bp[~has_bp].index)

    order = list(ranking_df["NR"])
    panels = [
        (df[df["antibiotic"].isin(bp_drugs)], "clsi_norm", bp_drugs),
        (df[df["antibiotic"].isin(other_drugs)], "modified_z", other_drugs),
    ]
    panels = [p for p in panels if p[2]]

    n_rows = len(order)
    widths = [len(drugs) for _, _, drugs in panels]
    fig, axes = plt.subplots(
        1,
        len(panels),
        figsize=(max(6, 0.5 * sum(widths) + 3 * len(panels)), max(5, 0.22 * n_rows)),
        gridspec_kw={"width_ratios": widths},
        squeeze=False,
    )
    for ax, (sub, value_col, drugs) in zip(axes[0], panels):
        pivot = sub.pivot_table(index="NR", columns="antibiotic", values=value_col, aggfunc="mean")
        pivot = pivot.reindex(index=order, columns=drugs)
        if value_col == "clsi_norm":
            cmap = matplotlib.colormaps["RdYlGn_r"].copy()
            norm = matplotlib.colors.TwoSlopeNorm(vmin=-2, vcenter=0.5, vmax=3)
            label, ticks = "MIC vs CLSI breakpoints", [0, 1]
        else:
            cmap = matplotlib.colormaps["coolwarm"].copy()
            norm = matplotlib.colors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
            label, ticks = "modified z (cohort; no CLSI breakpoint)", [-3, 0, 3]
        cmap.set_bad("white")
        im = ax.imshow(np.ma.masked_invalid(pivot.values), aspect="auto", cmap=cmap, norm=norm)
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_yticklabels(pivot.index, fontsize=6)
        ax.grid(False)
        cbar = fig.colorbar(im, ax=ax, shrink=0.6, ticks=ticks, location="bottom", pad=0.25)
        if value_col == "clsi_norm":
            cbar.ax.set_xticklabels(["S limit", "R limit"])
        cbar.set_label(label)
    axes[0][0].set_ylabel("isolate (NR), most to least resistant overall")
    fig.suptitle("MIC by isolate x antibiotic, normalized to CLSI (S limit = 0, R limit = 1); white = not tested")
    fig.tight_layout()
    path = out_dir / "heatmap_tnr_antibiotic.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_resistance_ranking(ranking_df: pd.DataFrame, out_dir: Path, top_n: int = 15) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    ranked = ranking_df.dropna(subset=["mean_modified_z"]).sort_values("mean_modified_z", ascending=False)
    top = ranked.head(top_n)
    bottom = ranked.tail(top_n)
    combined = pd.concat([top, bottom]).drop_duplicates(subset="TNR").sort_values("mean_modified_z")

    fig, ax = plt.subplots(figsize=(7, max(4, 0.32 * len(combined))))
    colors = [OUTLIER_RESISTANT_COLOR if v > 0 else OUTLIER_SUSCEPTIBLE_COLOR for v in combined["mean_modified_z"]]
    ax.barh(combined["NR"].astype(int).astype(str), combined["mean_modified_z"], color=colors)
    ax.axvline(0, color="#444444", linewidth=0.8)
    ax.set_xlabel("Mean modified z-score across tested antibiotics\n(> 0 = more resistant than cohort, < 0 = more susceptible)")
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
    for antibiotic_key, group in summary_df.groupby("antibiotic"):
        antibiotic = str(antibiotic_key)
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
    pivot = df.pivot_table(index="NR", columns="column", values="erg_code", aggfunc="mean")
    # pivot sorts the string labels ("20 mg/l" < "4 mg/l"); order numerically instead
    column_order = df.drop_duplicates("column").sort_values(["antibiotic", "concentration_mg_l"])["column"]
    pivot = pivot[column_order]
    isolate_order = [p for p in ranking_df["NR"] if p in pivot.index]
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
        "MGIT breakpoint result by isolate (NR) x antibiotic@concentration\n"
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
    ax.barh(combined["NR"].astype(int).astype(str), combined["n_R"], color=BREAKPOINT_R_COLOR)
    ax.set_xlabel("Number of R (resistant) MGIT breakpoint calls")
    ax.set_title(f"Most MGIT-resistant isolates by R-call count (top {top_n})")
    fig.tight_layout()
    path = out_dir / "resistance_ranking_mgit.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


_OVERVIEW_HEADERS = ["Antibiotic", "Susceptible\nconc. (mg/L)", "Intermediate\nconc. (mg/L)", "Resistant\nconc. (mg/L)", "n\nsusceptible", "n\nintermediate", "n\nresistant"]


def plot_overview_table(overview_df: pd.DataFrame, title: str, out_path: Path, footnote: str = "") -> Path:
    """Render an overview table (see overview.OVERVIEW_COLUMNS) as an image."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig_h = 0.9 + 0.36 * (len(overview_df) + 1)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.04)
    table = ax.table(
        cellText=overview_df.astype(str).to_numpy().tolist(),
        colLabels=_OVERVIEW_HEADERS,
        cellLoc="center",
        loc="upper center",
        colWidths=[0.29, 0.13, 0.16, 0.13, 0.09, 0.11, 0.09],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)
    count_colors = {4: BREAKPOINT_S_COLOR, 5: ERG_COLORS["I"], 6: BREAKPOINT_R_COLOR}
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#CCCCCC")
        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor(count_colors.get(col, "#444444"))
            cell.set_height(cell.get_height() * 1.5)
        elif col == 0:
            cell.set_text_props(weight="bold", ha="left")
            cell._loc = "left"  # type: ignore[attr-defined]  # no public setter for cell alignment
        elif row % 2 == 0:
            cell.set_facecolor("#F5F5F5")
    ax.set_title(title, fontsize=12, weight="bold")
    if footnote:
        fig.text(0.01, 0.01, footnote, fontsize=7, color="#555555", ha="left", va="bottom")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_mgit_mic_counts(mic_counts_df: pd.DataFrame, out_dir: Path) -> list[Path]:
    """One bar chart per antibiotic: MIC (lowest S concentration) on x, isolate count on y."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for antibiotic_key, group in mic_counts_df.groupby("antibiotic"):
        antibiotic = str(antibiotic_key)
        fig, ax = plt.subplots(figsize=(6, 3.4))
        x = np.arange(len(group))
        colors = [BREAKPOINT_R_COLOR if np.isnan(v) else POINT_COLOR for v in group["mic_mg_l"]]
        bars = ax.bar(x, group["n"], color=colors, width=0.6)
        ax.bar_label(bars, fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(group["mic_label"])
        ax.set_xlabel("MIC (mg/L) = lowest concentration with no growth (S)")
        ax.set_ylabel("Isolates (n)")
        ax.set_title(f"{antibiotic} -- MGIT MIC distribution (n={int(group['n'].sum())})")
        fig.tight_layout()
        path = out_dir / f"{_safe_filename(antibiotic)}_mgit_mic_counts.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_mic_count_table(counts_df: pd.DataFrame, title: str, out_path: Path) -> Path:
    """Antibiotic x MIC table of isolate counts (long df: antibiotic, mic_mg_l, n[, mic_label]).

    Rows with a ">x" ``mic_label`` (MIC above the tested range) sort after x. Blank = no
    isolates recorded / concentration not tested; an explicit 0 means tested, none found.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = counts_df.copy()
    if "mic_label" not in df:
        df["mic_label"] = [f"{v:.2g}" for v in df["mic_mg_l"]]
    df["key"] = [float(l[1:]) + 0.5 * 1e-9 if l.startswith(">") else float(l) for l in df["mic_label"]]
    order = df.drop_duplicates("mic_label").sort_values("key")["mic_label"].tolist()
    pivot = df.pivot_table(index="antibiotic", columns="mic_label", values="n", aggfunc="sum").reindex(columns=order)
    cells = [["" if pd.isna(v) else str(int(v)) for v in row] for row in pivot.to_numpy()]
    fig, ax = plt.subplots(figsize=(max(8, 0.55 * pivot.shape[1] + 3), 0.9 + 0.36 * (len(pivot) + 1)))
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.04)
    table = ax.table(
        cellText=cells,
        rowLabels=list(pivot.index),
        colLabels=list(pivot.columns),
        cellLoc="center",
        loc="upper center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    values = pivot.to_numpy()
    vmax = np.nanmax(values)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#CCCCCC")
        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor("#444444")
        elif col == -1:
            cell.set_text_props(weight="bold")
        else:
            v = values[row - 1, col]
            if not pd.isna(v) and v > 0:
                cell.set_facecolor(matplotlib.colors.to_hex(matplotlib.cm.Blues(0.15 + 0.6 * v / vmax)))
    ax.set_title(title, fontsize=12, weight="bold")
    fig.text(0.01, 0.01, "Columns: MIC (mg/L). Cells: number of isolates; blank = not tested / none recorded.", fontsize=7, color="#555555")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path
