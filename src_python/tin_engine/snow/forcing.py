"""The synthetic spring's weather and the sun's correction, per step, in NumPy.

The correction is the direct beam on each triangle's slope less the beam on
an open horizontal plane, averaged over the step's samples, so a flat
triangle gets exactly 0. The beam is clear-sky direct normal irradiance,
I0 (1/R)² Psi^(1/sin e), with one atmospheric pressure for the whole domain.
The weather is chosen, not measured: dry, with a temperature that rises
linearly over the run, swings daily and falls with height.
"""

from __future__ import annotations

import datetime as dt

import numpy as np
import numpy.typing as npt

from tin_engine import _core
from tin_engine.snow.params import SyntheticSpring

F64 = npt.NDArray[np.float64]
#: Solar constant (W/m²) and clear-sky transmissivity, a chosen value.
I0, PSI = 1368.0, 0.75
#: At 300 m: -3 °C at the start, +14 °C at the end, ±4 °C daily peaking at
#: 14:00 mean solar time; 0.6 °C colder per 100 m. Chosen values.
T_START, T_END, T_DAILY, PEAK_HOUR, BASE_M, LAPSE = -3.0, 14.0, 4.0, 14.0, 300.0, 0.006
_EPOCH = dt.datetime(1970, 1, 1, tzinfo=dt.UTC)


def beam(elevation_deg: npt.ArrayLike, distance_au: npt.ArrayLike) -> F64:
    """Clear-sky direct normal irradiance, W/m²; 0 with the sun at or below the horizon."""
    s = np.sin(np.radians(np.asarray(elevation_deg, dtype=np.float64)))
    up = s > 0.0
    r = np.asarray(distance_au, dtype=np.float64)
    return np.where(up, I0 / r**2 * PSI ** (1.0 / np.where(up, s, 1.0)), 0.0)


def radiation_anomaly(normals: F64, sun: F64, beam_wm2: F64) -> F64:
    """Per triangle, the mean over samples of beam · (max(0, n·u) - max(0, u_up)),
    W/m²; `normals` `(n, 3)`, `sun` `(S, 3)` unit, east-north-up, `beam_wm2` `(S,)`."""
    lit = np.maximum(normals @ sun.T, 0.0) - np.maximum(sun[:, 2], 0.0)
    out: F64 = lit @ beam_wm2 / len(beam_wm2)
    return out


def sun_and_beam(
    times_utc_ms: npt.NDArray[np.int64], lat_deg: float, lon_deg: float, height_m: float
) -> tuple[F64, F64]:
    """The sun's unit direction (east, north, up) `(S, 3)` and the beam `(S,)` at each time."""
    p = _core.sun_positions(np.asarray(times_utc_ms, dtype=np.int64), lat_deg, lon_deg, height_m)
    e, a = np.radians(p.elevation), np.radians(p.azimuth)
    sun = np.column_stack([np.cos(e) * np.sin(a), np.cos(e) * np.cos(a), np.sin(e)])
    return sun, beam(p.elevation, p.distance)


def sample_times(spring: SyntheticSpring, k: int) -> npt.NDArray[np.int64]:
    """UTC milliseconds of the midpoints of step `k`'s equal sub-intervals."""
    step_ms, s = spring.step_minutes * 60_000, spring.samples_per_step
    start = (spring.start - _EPOCH) // dt.timedelta(milliseconds=1) + k * step_ms
    return start + (2 * np.arange(s, dtype=np.int64) + 1) * step_ms // (2 * s)


def temperature(spring: SyntheticSpring, k: int, heights: F64, lon_deg: float) -> F64:
    """Air temperature, °C, at step `k`'s midpoint at each height (m); mean
    solar time is UTC + longitude/15 hours."""
    hours = (k + 0.5) * spring.step_minutes / 60.0
    of_day = spring.start.hour + spring.start.minute / 60.0 + spring.start.second / 3600.0
    solar = (of_day + hours + lon_deg / 15.0) % 24.0
    base = (
        T_START
        + (T_END - T_START) * hours / (spring.days * 24)
        + T_DAILY * np.cos(2.0 * np.pi * (solar - PEAK_HOUR) / 24.0)
    )
    out: F64 = base - LAPSE * (np.asarray(heights, dtype=np.float64) - BASE_M)
    return out
