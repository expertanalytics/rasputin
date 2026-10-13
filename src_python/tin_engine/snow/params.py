"""The snow model's parameters as frozen Pydantic models.

They refuse what the C++ kernel would refuse, earlier and by field name; the
kernel still checks, since it can be called without them. The defaults are
chosen values for Norwegian spring snow, explained in
`docs/increments/46-hbv-snow.md`, *5. Constants*.
"""

from __future__ import annotations

import math
from statistics import NormalDist
from typing import Literal, Self

from pydantic import BaseModel, ConfigDict, Field

_FROZEN = ConfigDict(frozen=True, allow_inf_nan=False)


class HbvSnowParameters(BaseModel):
    """The HBV snow routine's five parameters, under Shyft's names."""

    model_config = _FROZEN

    #: Precipitation falls as snow below this temperature, in °C.
    tx: float = 0.0
    #: Melt threshold, in °C.
    ts: float = 0.0
    #: Degree-day factor, in mm per °C per day.
    cx: float = Field(default=3.0, ge=0.0)
    #: Liquid water the pack holds, as a fraction of its solid part.
    lw: float = Field(default=0.1, ge=0.0)
    #: Refreeze coefficient: refreeze is cfr·cx·(ts - T) per day.
    cfr: float = Field(default=0.05, ge=0.0)


class SnowDistribution(BaseModel):
    """Relative snow depth in equal-area bins within a triangle, ascending, mean 1."""

    model_config = _FROZEN

    factors: tuple[float, ...]

    @classmethod
    def lognormal(cls, cv: float, bins: int) -> Self:
        """Bins of a log-normal snow depth of coefficient of variation `cv`.

        Bin k's factor is the distribution's quantile at the bin's midpoint
        probability (k + 1/2)/bins, divided by the factors' mean so they average
        to 1. A cv of 0 is uniform snow, one bin. Refuses bins outside 1 to 64
        and a negative or non-finite cv with a ValueError naming the argument.
        """
        if not 1 <= bins <= 64:
            raise ValueError(f"bins must be between 1 and 64, not {bins}")
        if not (math.isfinite(cv) and cv >= 0.0):
            raise ValueError(f"cv must be finite and not negative, not {cv}")
        if cv == 0.0:
            return cls(factors=(1.0,))
        sigma = math.sqrt(math.log1p(cv * cv))
        unit = NormalDist()
        raw = [math.exp(sigma * unit.inv_cdf((k + 0.5) / bins)) for k in range(bins)]
        mean = math.fsum(raw) / bins
        return cls(factors=tuple(v / mean for v in raw))


class SunCorrection(BaseModel):
    """How the sun corrects melt: not at all, through cx (r) or through ts (d)."""

    model_config = _FROZEN

    mode: Literal["none", "melt_factor", "threshold"]
    #: Melt-factor form: mm per °C per day per W/m² of radiation anomaly.
    r: float = 0.015
    #: Threshold form: °C per W/m² of radiation anomaly.
    d: float = 0.02
