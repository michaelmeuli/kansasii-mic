"""Resistance cutoffs for the M. kansasii complex, shared by every repo that calls S/I/R or resistant.

CANONICAL COPY: kansasii-docs/templates/kansasii_cutoffs.py. Repos that need cutoffs vendor this file unchanged
(next to their code; `repos/check_sync.sh` flags drift and `--fix` restores it). Edit the template, never a copy.
Data only: the classification logic (categorize, MIC string parsing) stays in the repos.

Drug keys are the canonical names of kansasii-mic (`kansasii_mic.loading`); `DRUG_CODES` gives the three-letter
codes of the screening_map_results.csv columns (`<CODE>_mic`, `<CODE>_clsi_category`) and the iTOL rings.

CLSI breakpoints are the species-specific *M. kansasii* values of CLSI M24 / M62, as summarised in van Ingen,
PMC6760954 (mg/L). Drugs without one use a single threshold: resistant if MIC > threshold (MIC == threshold is S),
no I category.
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
    # kansasii_mic.parsing.parse_mic reports the point estimate for combo values.
    "Sulfamethoxazole/Trimethoprim": Breakpoint(susceptible_max=2, resistant_min=4),
}

# No CLSI M. kansasii breakpoint: one threshold in mg/L, resistant if MIC > threshold.
NON_CLSI_THRESHOLDS_MG_L: dict[str, float] = {
    "Ethambutol": 8.0,
    "Isoniazid": 8.0,
    "Ethionamid": 0.6,
    "Streptomycin": 16.0,
    "Clofazimin": 2.0,
}

# canonical drug name -> three-letter code (column prefix in screening_map_results.csv, iTOL ring)
DRUG_CODES: dict[str, str] = {
    "Amikacin": "AMK", "Ciprofloxacin": "CIP", "Clarithromycin": "CLR",
    "Doxycyclin": "DOX", "Ethambutol": "EMB", "Isoniazid": "INH",
    "Linezolid": "LZD", "Moxifloxacin": "MXF", "Rifabutin": "RFB",
    "Rifampicin": "RIF", "Streptomycin": "STR",
    "Sulfamethoxazole/Trimethoprim": "SXT", "Ethionamid": "ETO",
    "Clofazimin": "CFZ", "Minocyclin": "MNO",
}


def non_clsi_thresholds_by_code() -> dict[str, float]:
    """`NON_CLSI_THRESHOLDS_MG_L` keyed by three-letter code (EMB, INH, ETO, STR, CFZ)."""
    return {DRUG_CODES[name]: t for name, t in NON_CLSI_THRESHOLDS_MG_L.items()}
