"""The snow run: the synthetic spring stepped on every triangle, frame by frame.

An async generator, so a GUI or a server can drive it: the steps between two
frames run in one `asyncio.to_thread` call, and the kernel releases the GIL.
The sun's correction is computed under every mode, `none` included, so its
mean can be written whatever the mode; the kernel receives it only under a
correction.
"""

from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from tin_engine import _core
from tin_engine.snow.forcing import radiation_anomaly, sample_times, sun_and_beam, temperature
from tin_engine.snow.geometry import TriangleGeometry
from tin_engine.snow.params import SnowJob

F32 = npt.NDArray[np.float32]


@dataclass(frozen=True)
class Frame:
    """The state at `time_hours` after the start, float32 copies per triangle:
    `swe` (mm), `sca` (0 to 1), `outflow` (mm since the start) and
    `radiation_anomaly_mean` (W/m², over the steps so far; 0 at the start)."""

    time_hours: float
    swe: F32
    sca: F32
    outflow: F32
    radiation_anomaly_mean: F32


async def run(
    job: SnowJob, geometry: TriangleGeometry, *, lat_deg: float, lon_deg: float, height_m: float
) -> AsyncIterator[Frame]:
    """Frames of `job` on `geometry`, the sun as seen from (`lat_deg`,
    `lon_deg`, `height_m`): the initial state, then one every
    `job.frames.every_hours` to the end of the spring."""
    n, spring = len(geometry.centroid_z), job.spring
    factors = np.array(job.distribution.factors)
    state = _core.SnowState.uniform(n, _core.SnowDistribution(factors), job.initial_swe_mm)
    mode = getattr(_core.SunCorrection, job.correction.mode)
    per_frame = job.frames.every_hours * 60 // spring.step_minutes
    outflow, di_sum, dry, none = np.zeros(n), np.zeros(n), np.zeros(n), np.empty(0)
    hbv = job.hbv.model_dump()

    def advance(first: int) -> None:
        for k in range(first, first + per_frame):
            sun, beam = sun_and_beam(sample_times(spring, k), lat_deg, lon_deg, height_m)
            di = radiation_anomaly(geometry.normals, sun, beam)
            state.step(
                **hbv, mode=mode, r=job.correction.r, d=job.correction.d,
                temperature=temperature(spring, k, geometry.centroid_z, lon_deg),
                precipitation=dry, radiation_anomaly=none if job.correction.mode == "none" else di,
                dt_seconds=spring.step_minutes * 60.0, threads=job.resources.threads,
            )  # fmt: skip
            outflow[:] += state.outflow
            di_sum[:] += di

    def frame(steps: int) -> Frame:
        hours = steps * spring.step_minutes / 60.0
        mean = di_sum / max(steps, 1)
        copies = (
            np.asarray(a, dtype=np.float32).copy() for a in (state.swe, state.sca, outflow, mean)
        )
        return Frame(hours, *copies)

    yield frame(0)
    for q in range(spring.days * 24 // job.frames.every_hours):
        await asyncio.to_thread(advance, q * per_frame)
        yield frame((q + 1) * per_frame)
