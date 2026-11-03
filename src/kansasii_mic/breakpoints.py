"""CLSI reference breakpoints for *Mycobacterium kansasii*.

*M. kansasii* ATCC 12478 is the CLSI-designated quality-control strain used
with the Sensititre SLOMYCO broth microdilution panels (this dataset's
source), primarily to validate rifampin MIC testing. ATCC 12478 has no MIC
results of its own in this dataset's ``screening_map_link.csv`` (it appears only
as an unlinked control-row placeholder), so these are literature/CLSI values
rather than an in-dataset control run.

Values below (S = susceptible, I = intermediate, R = resistant, in mg/L)
are species-specific CLSI breakpoints for *M. kansasii* as summarized in:

- CLSI M24 (Susceptibility Testing of Mycobacteria, Nocardia spp., and
  Other Aerobic Actinomycetes) / CLSI M62 (QC ranges supplement).
- van Ingen J. "Antimycobacterial Susceptibility Testing of Nontuberculous
  Mycobacteria" (review summarizing the CLSI M. kansasii breakpoint table),
  PMC6760954.

Rifampin is the primary/first-line breakpoint; the second-line agents are
only formally indicated by CLSI when an isolate is rifampin-resistant, but
are listed here for every isolate for descriptive/QC purposes. No CLSI
*M. kansasii*-specific breakpoint exists for ethambutol, isoniazid,
streptomycin, ethionamide, bedaquiline or clofazimine -- consistent with
those drugs being marked "K" (no assessment criteria) in this dataset's
own ``INT_ERG``/``ERG`` columns. For ethambutol, isoniazid and clofazimine a
single non-CLSI threshold is used instead (``NON_CLSI_THRESHOLDS_MG_L``, defined in the shared
``kansasii_cutoffs``, together with the CLSI table).

Values are expressed as (susceptible upper bound, resistant lower bound);
``None`` means no CLSI cutoff is defined for that side. Antibiotic keys
match the canonical names produced by :mod:`kansasii_mic.loading`.
"""

from __future__ import annotations

import math

# The cutoff data (Breakpoint, CLSI_KANSASII_BREAKPOINTS, NON_CLSI_THRESHOLDS_MG_L) is shared with the iTOL
# tree script: vendored copy of kansasii-docs/templates/kansasii_cutoffs.py, edit the template only.
from .kansasii_cutoffs import CLSI_KANSASII_BREAKPOINTS, NON_CLSI_THRESHOLDS_MG_L, Breakpoint

__all__ = ["CLSI_KANSASII_BREAKPOINTS", "NON_CLSI_THRESHOLDS_MG_L", "Breakpoint", "categorize",
           "categorize_interval", "clsi_normalized", "has_cutoff"]


def has_cutoff(antibiotic: str) -> bool:
    """True if the drug has a CLSI breakpoint pair or a non-CLSI threshold."""
    bp = CLSI_KANSASII_BREAKPOINTS.get(antibiotic)
    if bp is not None:
        return bp.susceptible_max is not None and bp.resistant_min is not None
    return antibiotic in NON_CLSI_THRESHOLDS_MG_L


def categorize(antibiotic: str, point_estimate: float) -> str | None:
    """Return "S", "I" or "R" per the CLSI table / non-CLSI threshold, or None if undefined."""
    threshold = NON_CLSI_THRESHOLDS_MG_L.get(antibiotic)
    if threshold is not None:
        return "R" if point_estimate > threshold else "S"
    bp = CLSI_KANSASII_BREAKPOINTS.get(antibiotic)
    if bp is None:
        return None
    if bp.susceptible_max is not None and point_estimate <= bp.susceptible_max:
        return "S"
    if bp.resistant_min is not None and point_estimate >= bp.resistant_min:
        return "R"
    return "I"


def categorize_interval(antibiotic: str, lower: float, upper: float) -> str | None:
    """Classify a censored MGIT MIC lying in ``(lower, upper]`` (NaN = open end).

    ``upper`` is the lowest concentration called S (NaN: no S call, MIC above the tested
    range); ``lower`` the highest tested concentration below it (NaN: none). The next
    dilution above ``lower`` (``2 * lower``) is taken as the smallest possible MIC.
    Returns "S" / "I" / "R", or None if the interval straddles a limit/threshold or the
    drug has no cutoff. Non-CLSI drugs: S if ``upper <= t``, R if ``lower >= t``.
    """
    has_lower, has_upper = not math.isnan(lower), not math.isnan(upper)
    threshold = NON_CLSI_THRESHOLDS_MG_L.get(antibiotic)
    if threshold is not None:
        if has_upper and upper <= threshold:
            return "S"
        return "R" if has_lower and lower >= threshold else None
    bp = CLSI_KANSASII_BREAKPOINTS.get(antibiotic)
    if bp is None or bp.susceptible_max is None or bp.resistant_min is None:
        return None
    has_lower, has_upper = not math.isnan(lower), not math.isnan(upper)
    if has_upper and upper <= bp.susceptible_max:
        return "S"
    if not has_lower:
        return None
    if 2 * lower >= bp.resistant_min:
        return "R"
    if has_upper and lower >= bp.susceptible_max and upper < bp.resistant_min:
        return "I"
    return None


def clsi_normalized(antibiotic: str, log2_mic: float) -> float | None:
    """Place a log2 MIC on a per-drug CLSI scale: 0 = S limit, 1 = R limit.

    Values below 0 are on the susceptible side, between 0 and 1 intermediate,
    above 1 resistant; each unit spans the drug's S-to-R gap in dilutions.
    None if the drug has no usable (S and R) breakpoint pair.
    """
    bp = CLSI_KANSASII_BREAKPOINTS.get(antibiotic)
    if bp is None or bp.susceptible_max is None or bp.resistant_min is None:
        return None
    lo, hi = math.log2(bp.susceptible_max), math.log2(bp.resistant_min)
    if hi <= lo:
        return None
    return (log2_mic - lo) / (hi - lo)
