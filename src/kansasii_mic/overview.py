"""Per-antibiotic overview tables: S/I/R concentrations and isolate counts."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from .breakpoints import CLSI_KANSASII_BREAKPOINTS, NON_CLSI_THRESHOLDS_MG_L, categorize, categorize_interval, has_cutoff

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
    ``nr_unclassifiable``. Only antibiotics with a CLSI *M. kansasii* breakpoint or a non-CLSI threshold are included.
    """
    mics = mgit_isolate_mics(mgit_long_df)
    rows = []
    for antibiotic_key, group in mics.groupby("antibiotic"):
        antibiotic = str(antibiotic_key)
        if not has_cutoff(antibiotic):
            continue
        cats = pd.Series(
            [categorize_interval(antibiotic, lo, up) for lo, up in zip(group["lower"], group["upper"], strict=True)],
            dtype=object,
        ).value_counts()
        rows.append(
            {
                "antibiotic": antibiotic,
                **_cutoff_columns(antibiotic),
                "nr_susceptible": int(cats.get("S", 0)),
                "nr_intermediate": int(cats.get("I", 0)),
                "nr_resistant": int(cats.get("R", 0)),
                "nr_unclassifiable": len(group) - int(sum(cats.get(c, 0) for c in "SIR")),
            }
        )
    return pd.DataFrame(rows, columns=MGIT_OVERVIEW_COLUMNS)


def _cutoff_columns(antibiotic: str) -> dict[str, str]:
    """S / I / R range strings: CLSI limits, or ``<= t`` / ``--`` / ``> t`` for a non-CLSI threshold."""
    if (t := NON_CLSI_THRESHOLDS_MG_L.get(antibiotic)) is not None:
        return {
            "susceptible_concentration": f"<= {t:g}",
            "intermediate_concentration": "--",
            "resistant_concentration": f"> {t:g}",
        }
    bp = CLSI_KANSASII_BREAKPOINTS[antibiotic]
    assert bp.susceptible_max is not None and bp.resistant_min is not None
    s, r = bp.susceptible_max, bp.resistant_min
    return {
        "susceptible_concentration": f"<= {s:g}",
        "intermediate_concentration": _intermediate_label(s, r),
        "resistant_concentration": f">= {r:g}",
    }


def _intermediate_label(s: float, r: float) -> str:
    """Doubling dilutions strictly between S and R limits, e.g. "32", "2-4", or "--" if none."""
    lo, hi = s * 2, r / 2
    if lo > hi:
        return "--"
    return f"{lo:g}" if lo == hi else f"{lo:g}-{hi:g}"


def mhk_overview(mic_df: pd.DataFrame) -> pd.DataFrame:
    """MHK: CLSI breakpoint ranges (mg/L) and isolates per category.

    Only antibiotics with a CLSI *M. kansasii* breakpoint or a non-CLSI threshold are included.
    """
    rows = []
    for antibiotic_key, group in mic_df.groupby("antibiotic"):
        antibiotic = str(antibiotic_key)
        if not has_cutoff(antibiotic):
            continue
        calls = [categorize(antibiotic, v) for v in group["point_estimate"]]
        cats = pd.Series(calls, dtype=object).value_counts()
        rows.append(
            {
                "antibiotic": antibiotic,
                **_cutoff_columns(antibiotic),
                "nr_susceptible": int(cats.get("S", 0)),
                "nr_intermediate": int(cats.get("I", 0)),
                "nr_resistant": int(cats.get("R", 0)),
            }
        )
    return pd.DataFrame(rows, columns=OVERVIEW_COLUMNS)


def mgit_mic_counts(mgit_long_df: pd.DataFrame) -> pd.DataFrame:
    """MGIT: isolates per MIC, taking MIC = lowest tested concentration called S (no growth).

    Isolates with no S call but at least one I/R call are counted at
    ``">" + the highest concentration tested for that isolate`` (MIC above its tested
    range; isolates not tested up to the drug's top concentration get a lower ">" bin).
    The ``">" + drug maximum`` bin is always listed. Isolates with only K/U calls are
    skipped. Tested concentrations with no isolates are kept with n = 0, rows ordered by
    concentration with each ">c" bin right after c. Columns: antibiotic, mic_mg_l (NaN
    for ">" bins), mic_label, n.
    """
    mics = mgit_isolate_mics(mgit_long_df)
    rows = []
    for antibiotic, group in mgit_long_df.groupby("antibiotic"):
        tested = sorted(group["concentration_mg_l"].unique())
        iso = mics[mics["antibiotic"] == antibiotic]
        counts = iso["upper"].value_counts()
        above = iso.loc[iso["upper"].isna(), "lower"].value_counts()
        bins = {c: int(counts.get(c, 0)) for c in tested}
        above_bins = {c: int(above.get(c, 0)) for c in {tested[-1], *above.index}}
        for c, n in bins.items():
            rows.append((c, 0, {"antibiotic": antibiotic, "mic_mg_l": c, "mic_label": f"{c:g}", "n": n}))
        for c, n in above_bins.items():
            rows.append((c, 1, {"antibiotic": antibiotic, "mic_mg_l": float("nan"), "mic_label": f">{c:g}", "n": n}))
        # keep each antibiotic's rows contiguous, ordered by concentration (">c" right after c)
    rows.sort(key=lambda r: (r[2]["antibiotic"], r[0], r[1]))
    return pd.DataFrame([r[2] for r in rows], columns=["antibiotic", "mic_mg_l", "mic_label", "n"])


def mhk_mic_counts(mic_df: pd.DataFrame) -> pd.DataFrame:
    """MHK: isolates per MIC point estimate per antibiotic (long format: antibiotic, mic_mg_l, n)."""
    counts = mic_df.groupby(["antibiotic", "point_estimate"]).size().rename("n").reset_index()
    return counts.rename(columns={"point_estimate": "mic_mg_l"})
