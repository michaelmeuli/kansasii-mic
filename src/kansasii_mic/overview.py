"""Per-antibiotic overview tables: S/I/R concentrations and isolate counts."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from .breakpoints import CLSI_KANSASII_BREAKPOINTS, categorize, categorize_interval

OVERVIEW_COLUMNS = [
    "antibiotic",
    "susceptible_concentration",
    "intermediate_concentration",
    "resistant_concentration",
    "nr_susceptible",
    "nr_intermediate",
    "nr_resistant",
]

# Worst call wins when an isolate has several calls for one antibiotic.
_SEVERITY = {"S": 0, "I": 1, "R": 2}


def _fmt(values: Iterable[float]) -> str:
    return ", ".join(f"{v:g}" for v in sorted(set(values))) or "-"


MGIT_OVERVIEW_COLUMNS = [*OVERVIEW_COLUMNS, "nr_unclassifiable"]

# erg is the final call (also used by the duplicate resolution and the summary).
_CALL_COL = "erg"


def mgit_isolate_mics(mgit_long_df: pd.DataFrame) -> pd.DataFrame:
    """One row per TNR x antibiotic with the MIC interval implied by the MGIT calls.

    MIC = lowest tested concentration called S (no growth), so the true MIC lies in
    ``(lower, upper]``: ``upper`` is that concentration (NaN if no S call, i.e. MIC
    above the tested range) and ``lower`` the highest tested concentration below it
    (NaN if none). K/U calls are ignored; isolates with no S/I/R call are skipped.
    """
    rows = []
    calls = mgit_long_df[mgit_long_df[_CALL_COL].isin(_SEVERITY)]
    for (tnr, antibiotic), iso in calls.groupby(["TNR", "antibiotic"]):
        conc = iso["concentration_mg_l"]
        s = iso.loc[iso[_CALL_COL] == "S", "concentration_mg_l"]
        upper = float(s.min()) if len(s) else float("nan")
        below = conc[conc < upper] if len(s) else conc
        lower = float(below.max()) if len(below) else float("nan")
        rows.append({"TNR": tnr, "antibiotic": antibiotic, "lower": lower, "upper": upper})
    return pd.DataFrame(rows, columns=["TNR", "antibiotic", "lower", "upper"])


def mgit_overview(mgit_long_df: pd.DataFrame) -> pd.DataFrame:
    """MGIT: CLSI breakpoint ranges (mg/L) and isolates per category, from the MGIT MIC.

    Each isolate x antibiotic is classified once from its MIC interval (see
    :func:`mgit_isolate_mics`) with :func:`categorize_interval`; MICs the tested
    concentrations cannot place relative to the CLSI limits are counted in
    ``nr_unclassifiable``. Only antibiotics with a CLSI *M. kansasii* breakpoint are included.
    """
    mics = mgit_isolate_mics(mgit_long_df)
    rows = []
    for antibiotic_key, group in mics.groupby("antibiotic"):
        antibiotic = str(antibiotic_key)
        bp = CLSI_KANSASII_BREAKPOINTS.get(antibiotic)
        if bp is None or bp.susceptible_max is None or bp.resistant_min is None:
            continue
        cats = pd.Series(
            [categorize_interval(antibiotic, lo, up) for lo, up in zip(group["lower"], group["upper"], strict=True)],
            dtype=object,
        ).value_counts()
        s, r = bp.susceptible_max, bp.resistant_min
        rows.append(
            {
                "antibiotic": antibiotic,
                "susceptible_concentration": f"<= {s:g}",
                "intermediate_concentration": _intermediate_label(s, r),
                "resistant_concentration": f">= {r:g}",
                "nr_susceptible": int(cats.get("S", 0)),
                "nr_intermediate": int(cats.get("I", 0)),
                "nr_resistant": int(cats.get("R", 0)),
                "nr_unclassifiable": len(group) - int(sum(cats.get(c, 0) for c in "SIR")),
            }
        )
    return pd.DataFrame(rows, columns=MGIT_OVERVIEW_COLUMNS)


def _intermediate_label(s: float, r: float) -> str:
    """Doubling dilutions strictly between S and R limits, e.g. "32", "2-4", or "--" if none."""
    lo, hi = s * 2, r / 2
    if lo > hi:
        return "--"
    return f"{lo:g}" if lo == hi else f"{lo:g}-{hi:g}"


def mhk_overview(mic_df: pd.DataFrame) -> pd.DataFrame:
    """MHK: CLSI breakpoint ranges (mg/L) and isolates per category.

    Only antibiotics with a CLSI *M. kansasii* breakpoint are included.
    """
    rows = []
    for antibiotic_key, group in mic_df.groupby("antibiotic"):
        antibiotic = str(antibiotic_key)
        bp = CLSI_KANSASII_BREAKPOINTS.get(antibiotic)
        if bp is None or bp.susceptible_max is None or bp.resistant_min is None:
            continue
        calls = [categorize(antibiotic, v) for v in group["point_estimate"]]
        cats = pd.Series(calls, dtype=object).value_counts()
        s, r = bp.susceptible_max, bp.resistant_min
        rows.append(
            {
                "antibiotic": antibiotic,
                "susceptible_concentration": f"<= {s:g}",
                "intermediate_concentration": _intermediate_label(s, r),
                "resistant_concentration": f">= {r:g}",
                "nr_susceptible": int(cats.get("S", 0)),
                "nr_intermediate": int(cats.get("I", 0)),
                "nr_resistant": int(cats.get("R", 0)),
            }
        )
    return pd.DataFrame(rows, columns=OVERVIEW_COLUMNS)


def mgit_mic_counts(mgit_long_df: pd.DataFrame) -> pd.DataFrame:
    """MGIT: isolates per MIC, taking MIC = lowest tested concentration called S (no growth).

    Isolates with no S call but at least one I/R call are counted at
    ``">" + highest tested concentration`` (MIC above the tested range).
    Isolates with only K/U calls are skipped. Tested concentrations with no
    isolates are kept with n = 0. Columns: antibiotic, mic_mg_l (NaN for the
    ">" bin), mic_label, n.
    """
    mics = mgit_isolate_mics(mgit_long_df)
    rows = []
    for antibiotic, group in mgit_long_df.groupby("antibiotic"):
        tested = sorted(group["concentration_mg_l"].unique())
        iso = mics[mics["antibiotic"] == antibiotic]
        counts = iso["upper"].value_counts()
        for c in tested:
            rows.append({"antibiotic": antibiotic, "mic_mg_l": c, "mic_label": f"{c:g}", "n": int(counts.get(c, 0))})
        rows.append(
            {"antibiotic": antibiotic, "mic_mg_l": float("nan"), "mic_label": f">{tested[-1]:g}", "n": int(iso["upper"].isna().sum())}
        )
    return pd.DataFrame(rows, columns=["antibiotic", "mic_mg_l", "mic_label", "n"])


def mhk_mic_counts(mic_df: pd.DataFrame) -> pd.DataFrame:
    """MHK: isolates per MIC point estimate per antibiotic (long format: antibiotic, mic_mg_l, n)."""
    counts = mic_df.groupby(["antibiotic", "point_estimate"]).size().rename("n").reset_index()
    return counts.rename(columns={"point_estimate": "mic_mg_l"})
