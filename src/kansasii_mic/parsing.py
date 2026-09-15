"""Parsing of raw ``MHK`` (MIC) strings from the lab export.

Observed value shapes in ``data/imm/mic.csv`` (MHK column, non-empty rows only):

- plain doubling-dilution number: ``"0.5"``, ``"10"``
- left-censored (below lowest tested concentration): ``"<0.25"``, ``"<=0.12"``
- right-censored (above highest tested concentration): ``">8"``, ``">16"``
- a range spanning two adjacent tested concentrations: ``"4-8"``, ``"0.25-0.5"``
- a fixed-ratio combination MIC (Sulfamethoxazole/Trimethoprim, 19:1 ratio):
  ``"0.12/2.38"``, ``">8/152"``
- a combo range: ``"0.25/4.75-0.5/9.5"``
- one known data-entry typo, a combo range missing a slash on the upper
  bound: ``"0.12/2.38-0.25-4.75"`` (should read ``"0.12/2.38-0.25/4.75"``,
  confirmed against the correctly-formatted equivalent elsewhere in the
  file). Repaired generically, not by hardcoding this one row.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

_NUM = r"[0-9]+(?:\.[0-9]+)?"
_CENSOR_RE = re.compile(rf"^(<=|<|>=|>)\s*({_NUM})(?:/({_NUM}))?$")
_PLAIN_RE = re.compile(rf"^({_NUM})(?:/({_NUM}))?$")
_RANGE_RE = re.compile(
    rf"^({_NUM})(?:/({_NUM}))?\s*-\s*({_NUM})(?:/({_NUM}))?$"
)
# Known typo shape: "a/b-c-d" where a slash was dropped before the last
# number of the upper bound of a combo range (e.g. "0.12/2.38-0.25-4.75").
_TYPO_COMBO_RANGE_RE = re.compile(
    rf"^({_NUM})/({_NUM})-({_NUM})-({_NUM})$"
)


class MicParseError(ValueError):
    """Raised when an MHK string cannot be parsed at all."""


@dataclass(frozen=True)
class ParsedMic:
    raw: str
    point_estimate: float  # numerator component, on the original mg/L scale
    denominator: float | None  # SMX/TMP ratio partner, if a combo value
    low: float
    high: float
    censored: str | None  # None, "left" (<...) or "right" (>...)
    is_range: bool
    repaired_from: str | None  # original raw string if a typo was repaired

    @property
    def log2(self) -> float:
        return math.log2(self.point_estimate)


def _repair_typo(raw: str) -> str:
    m = _TYPO_COMBO_RANGE_RE.match(raw)
    if m:
        a, b, c, d = m.groups()
        return f"{a}/{b}-{c}/{d}"
    return raw


def parse_mic(raw: str) -> ParsedMic:
    """Parse a single raw MHK string into a :class:`ParsedMic`.

    Point estimate conventions (documented assumptions, see README):
    - censored ``<X`` / ``<=X`` / ``>X`` / ``>=X`` -> boundary value ``X``
    - range ``X-Y`` -> geometric mean of the bounds (arithmetic mean on the
      log2 scale doubling dilutions live on)
    - combo ``X/Y`` -> the first (numerator) component, matching how CLSI
      expresses the Sulfamethoxazole/Trimethoprim breakpoint (e.g. ``2/38``)
    """
    text = raw.strip()
    if not text:
        raise MicParseError("empty MHK value")

    repaired_from = None
    repaired = _repair_typo(text)
    if repaired != text:
        repaired_from = text
        text = repaired

    m = _CENSOR_RE.match(text)
    if m:
        direction, num, denom = m.groups()
        value = float(num)
        censored = "left" if direction in ("<", "<=") else "right"
        return ParsedMic(
            raw=raw,
            point_estimate=value,
            denominator=float(denom) if denom else None,
            low=value,
            high=value,
            censored=censored,
            is_range=False,
            repaired_from=repaired_from,
        )

    m = _RANGE_RE.match(text)
    if m:
        low_num, low_denom, high_num, high_denom = m.groups()
        low = float(low_num)
        high = float(high_num)
        point = math.sqrt(low * high)
        denom = None
        if low_denom and high_denom:
            denom = math.sqrt(float(low_denom) * float(high_denom))
        return ParsedMic(
            raw=raw,
            point_estimate=point,
            denominator=denom,
            low=low,
            high=high,
            censored=None,
            is_range=True,
            repaired_from=repaired_from,
        )

    m = _PLAIN_RE.match(text)
    if m:
        num, denom = m.groups()
        value = float(num)
        return ParsedMic(
            raw=raw,
            point_estimate=value,
            denominator=float(denom) if denom else None,
            low=value,
            high=value,
            censored=None,
            is_range=False,
            repaired_from=repaired_from,
        )

    raise MicParseError(f"unrecognized MHK format: {raw!r}")
