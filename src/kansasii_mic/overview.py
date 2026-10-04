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
                "intermediate_concentration": f"> {s:g} and < {r:g}",
                "resistant_concentration": f">= {r:g}",
                "nr_susceptible": int(cats.get("S", 0)),
                "nr_intermediate": int(cats.get("I", 0)),
                "nr_resistant": int(cats.get("R", 0)),
            }
        )
    return pd.DataFrame(rows, columns=OVERVIEW_COLUMNS)
