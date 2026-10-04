"""Per-antibiotic overview tables: S/I/R concentrations and isolate counts."""

from __future__ import annotations

import pandas as pd

from .breakpoints import CLSI_KANSASII_BREAKPOINTS, categorize

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


def _fmt(values) -> str:
    return ", ".join(f"{v:g}" for v in sorted(set(values))) or "-"


def mgit_overview(mgit_long_df: pd.DataFrame) -> pd.DataFrame:
    """MGIT: concentrations (mg/L) at which S/I/R was called, and isolates by worst call.

    An isolate tested at several concentrations is counted once per antibiotic,
    as R if any concentration was R, else I if any was I, else S. K/U calls
    carry no S/I/R information and are ignored.
    """
    df = mgit_long_df[mgit_long_df["int_erg"].isin(_SEVERITY)]
    rows = []
    for antibiotic, group in df.groupby("antibiotic"):
        worst = group.assign(sev=group["int_erg"].map(_SEVERITY)).groupby("TNR")["sev"].max()
        conc = {c: group.loc[group["int_erg"] == c, "concentration_mg_l"] for c in _SEVERITY}
        rows.append(
            {
                "antibiotic": antibiotic,
                "susceptible_concentration": _fmt(conc["S"]),
                "intermediate_concentration": _fmt(conc["I"]),
                "resistant_concentration": _fmt(conc["R"]),
                "nr_susceptible": int((worst == 0).sum()),
                "nr_intermediate": int((worst == 1).sum()),
                "nr_resistant": int((worst == 2).sum()),
            }
        )
    return pd.DataFrame(rows, columns=OVERVIEW_COLUMNS)


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
    for antibiotic, group in mic_df.groupby("antibiotic"):
        bp = CLSI_KANSASII_BREAKPOINTS.get(antibiotic)
        if bp is None:
            continue
        cats = pd.Series([categorize(antibiotic, v) for v in group["point_estimate"]]).value_counts()
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
    rows = []
    for antibiotic, group in mgit_long_df.groupby("antibiotic"):
        tested = sorted(group["concentration_mg_l"].unique())
        counts = {c: 0 for c in tested}
        above = 0
        for _, iso in group[group["int_erg"].isin(_SEVERITY)].groupby("TNR"):
            s = iso.loc[iso["int_erg"] == "S", "concentration_mg_l"]
            if len(s):
                counts[s.min()] += 1
            else:
                above += 1
        for c in tested:
            rows.append({"antibiotic": antibiotic, "mic_mg_l": c, "mic_label": f"{c:g}", "n": counts[c]})
        rows.append({"antibiotic": antibiotic, "mic_mg_l": float("nan"), "mic_label": f">{tested[-1]:g}", "n": above})
    return pd.DataFrame(rows, columns=["antibiotic", "mic_mg_l", "mic_label", "n"])


def mhk_mic_counts(mic_df: pd.DataFrame) -> pd.DataFrame:
    """MHK: isolates per MIC point estimate per antibiotic (long format: antibiotic, mic_mg_l, n)."""
    counts = mic_df.groupby(["antibiotic", "point_estimate"]).size().rename("n").reset_index()
    return counts.rename(columns={"point_estimate": "mic_mg_l"})
