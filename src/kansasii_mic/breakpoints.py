"""CLSI reference breakpoints for *Mycobacterium kansasii*.

*M. kansasii* ATCC 12478 is the CLSI-designated quality-control strain used
with the Sensititre SLOMYCO broth microdilution panels (this dataset's
source), primarily to validate rifampin MIC testing. ATCC 12478 has no MIC
results of its own in this dataset's ``screening_map.csv`` (it appears only
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
own ``INT_ERG``/``ERG`` columns.

Values are expressed as (susceptible upper bound, resistant lower bound);
``None`` means no CLSI cutoff is defined for that side. Antibiotic keys
match the canonical names produced by :mod:`kansasii_mic.loading`.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Breakpoint:
    susceptible_max: float | None  # MIC <= this -> S
    resistant_min: float | None  # MIC >= this -> R


CLSI_KANSASII_BREAKPOINTS: dict[str, Breakpoint] = {
    "Clarithromycin": Breakpoint(susceptible_max=8, resistant_min=32),
    "Rifampicin": Breakpoint(susceptible_max=1, resistant_min=2),
    "Amikacin": Breakpoint(susceptible_max=16, resistant_min=64),
    "Ciprofloxacin": Breakpoint(susceptible_max=1, resistant_min=4),
    "Doxycyclin": Breakpoint(susceptible_max=1, resistant_min=8),
    "Linezolid": Breakpoint(susceptible_max=8, resistant_min=32),
    "Minocyclin": Breakpoint(susceptible_max=1, resistant_min=8),
    "Moxifloxacin": Breakpoint(susceptible_max=1, resistant_min=4),
    "Rifabutin": Breakpoint(susceptible_max=2, resistant_min=4),
    # Expressed via the trimethoprim (numerator) component, matching how
    # parsing.parse_mic reports the point estimate for combo values.
    "Sulfamethoxazole/Trimethoprim": Breakpoint(susceptible_max=2, resistant_min=4),
}


def categorize(antibiotic: str, point_estimate: float) -> str | None:
    """Return "S", "I" or "R" per the CLSI table, or None if undefined."""
    bp = CLSI_KANSASII_BREAKPOINTS.get(antibiotic)
    if bp is None:
        return None
    if bp.susceptible_max is not None and point_estimate <= bp.susceptible_max:
        return "S"
    if bp.resistant_min is not None and point_estimate >= bp.resistant_min:
        return "R"
    return "I"
