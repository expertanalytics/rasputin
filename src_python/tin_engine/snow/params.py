"""The snow model's parameters as frozen Pydantic models.

They refuse what the C++ kernel would refuse, earlier and by field name; the
kernel still checks, since it can be called without them. The defaults are
chosen values for Norwegian spring snow, explained in
`docs/increments/46-hbv-snow.md`, *5. Constants*.
"""

from __future__ import annotations

import datetime as dt
import math
from statistics import NormalDist
from typing import Literal, Self

from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, field_validator

from tin_engine.radiation.request import RadiationRequest
from tin_engine.resources import Resources

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


class SyntheticSpring(BaseModel):
    """The synthetic weather's time axis, into summer: when it starts, how long,
    how finely, and how long the warming lasts."""

    model_config = _FROZEN

    #: Start of the run, timezone-aware; kept in UTC.
    start: dt.datetime = dt.datetime(2027, 4, 1, tzinfo=dt.UTC)
    days: int = Field(default=153, ge=1)
    #: Days over which the temperature at 300 m rises from T_START to T_END,
    #: counted from `start`; after them it holds at T_END.
    warming_days: int = Field(default=122, ge=1)
    step_minutes: int = Field(default=60, ge=1)
    #: Sun positions per step, at the midpoints of equal sub-intervals.
    samples_per_step: int = Field(default=4, ge=1)

    @field_validator("start")
    @classmethod
    def _utc(cls, start: dt.datetime) -> dt.datetime:
        if start.utcoffset() is None:
            raise ValueError(f"{start.isoformat()} has no time zone; give it an offset or Z")
        return start.astimezone(dt.UTC)


class FrameSpec(BaseModel):
    """How often the run writes a frame of the animation."""

    model_config = _FROZEN

    every_hours: int = Field(default=24, ge=1)


class SnowJob(BaseModel):
    """Everything one snow run needs besides the mesh. Refuses frames that do
    not fall on whole steps, a spring that is not whole frames, and more
    frames than the animation's four-digit numbers can count."""

    model_config = _FROZEN

    correction: SunCorrection
    resources: Resources
    hbv: HbvSnowParameters = HbvSnowParameters()
    distribution: SnowDistribution = Field(
        default_factory=lambda: SnowDistribution.lognormal(0.0, 1)
    )
    spring: SyntheticSpring = SyntheticSpring()
    frames: FrameSpec = Field(default=FrameSpec(), validate_default=True)
    #: Snow water equivalent on every triangle at the start, mm.
    initial_swe_mm: float = Field(default=600.0, ge=0.0)
    #: The horizon search's settings for cast shadows; None: an open sky.
    shadows: RadiationRequest | None = None

    @field_validator("frames")
    @classmethod
    def _whole_steps(cls, frames: FrameSpec, info: ValidationInfo) -> FrameSpec:
        spring = info.data.get("spring")
        if spring is not None and frames.every_hours * 60 % spring.step_minutes:
            raise ValueError(
                f"a frame every {frames.every_hours} h is not a whole number of "
                f"{spring.step_minutes}-minute steps"
            )
        if spring is not None and spring.days * 24 % frames.every_hours:
            raise ValueError(
                f"{spring.days} days is not a whole number of {frames.every_hours}-hour frames"
            )
        # Frames are numbered from 0 at the start; the case file's swe****.var
        # pattern has four digits, so the last frame's number is at most 9999.
        if spring is not None and (n := spring.days * 24 // frames.every_hours + 1) > 10_000:
            raise ValueError(
                f"{spring.days} days with a frame every {frames.every_hours} h is {n:,} frames, "
                "more than the 10,000 the animation can number"
            )
        return frames
