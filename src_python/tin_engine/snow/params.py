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
    """Relative snow depth in equal-area bins within a triangle: each bin's mean
    (`factors`, ascending, mean 1) and its depth at the bin's deep edge
    (`deep`), the depth running linearly across the bin from 2·factor - deep to
    deep. An empty `deep` is flat: uniform depth in each bin."""

    model_config = _FROZEN

    factors: tuple[float, ...]
    deep: tuple[float, ...] = ()

    @field_validator("deep")
    @classmethod
    def _edges(cls, deep: tuple[float, ...], info: ValidationInfo) -> tuple[float, ...]:
        factors = info.data.get("factors")
        if not deep or factors is None:
            return deep
        if len(deep) != len(factors):
            raise ValueError(f"deep must be empty or one per bin: {len(deep)} for {len(factors)}")
        if not all(f <= h <= 2.0 * f for f, h in zip(factors, deep, strict=True)):
            raise ValueError("a bin's deep edge must lie between its mean and twice its mean")
        return deep

    @classmethod
    def lognormal(cls, cv: float, bins: int) -> Self:
        """Bins of a log-normal snow depth of mean 1 and coefficient of variation `cv`.

        For each of `bins` equal-probability bands, the factor is the band's
        mean and the deep edge its mean plus √3 times its standard deviation:
        the straight line of the band's mean and variance. Both are divided by
        the factors' mean. A cv of 0 is uniform snow, one flat bin. Refuses bins
        outside 1 to 64, a negative or non-finite cv, and a cv so large that a
        bin's shallow edge would be negative, with a ValueError naming cv.
        """
        if not 1 <= bins <= 64:
            raise ValueError(f"bins must be between 1 and 64, not {bins}")
        if not (math.isfinite(cv) and cv >= 0.0):
            raise ValueError(f"cv must be finite and not negative, not {cv}")
        if cv == 0.0:
            return cls(factors=(1.0,))
        factors, deep = _bands(cv, bins)
        if not _shallow_edges_ok(factors, deep):
            low, high = 0.0, cv
            for _ in range(50):
                mid = (low + high) / 2
                low, high = (mid, high) if _shallow_edges_ok(*_bands(mid, bins)) else (low, mid)
            unit = "bin" if bins == 1 else "bins"
            raise ValueError(f"snow-cv {cv:g} is too large for {bins} {unit}: at most {low:.3f}")
        return cls(factors=factors, deep=deep)


def _shallow_edges_ok(factors: tuple[float, ...], deep: tuple[float, ...]) -> bool:
    return all(2.0 * f >= h for f, h in zip(factors, deep, strict=True))


def _bands(cv: float, bins: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Each band's mean and mean + √3 sd, from the log-normal's partial moments
    E[X^j; a < X < b] = exp(j·mu + j²·s²/2)·(Φ(zb - j·s) - Φ(za - j·s)),
    divided by the means' mean."""
    s = math.sqrt(math.log1p(cv * cv))
    mu = -s * s / 2
    unit = NormalDist()

    def mass(z: float, shift: float) -> float:
        return 0.0 if z == -math.inf else 1.0 if z == math.inf else unit.cdf(z - shift)

    means, edges = [], []
    for k in range(bins):
        za = unit.inv_cdf(k / bins) if k > 0 else -math.inf
        zb = unit.inv_cdf((k + 1) / bins) if k + 1 < bins else math.inf
        m1 = math.exp(mu + s * s / 2) * (mass(zb, s) - mass(za, s)) * bins
        m2 = math.exp(2 * mu + 2 * s * s) * (mass(zb, 2 * s) - mass(za, 2 * s)) * bins
        means.append(m1)
        edges.append(m1 + math.sqrt(3.0) * math.sqrt(max(0.0, m2 - m1 * m1)))
    mean = math.fsum(means) / bins
    return tuple(v / mean for v in means), tuple(v / mean for v in edges)


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
        default_factory=lambda: SnowDistribution.lognormal(0.4, 1)
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
