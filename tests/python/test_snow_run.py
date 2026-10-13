"""The snow run as an async generator of frames: `snow/run.py`, test R7.

`docs/increments/46-hbv-snow.md`, *4.1 The blueprint* (per step: the
radiation anomaly, then `_core.SnowState.step` with the run's threads; per
frame: swe and sca copied out), *4.8 The Python surface* (`async def run(job,
geometry) -> AsyncIterator[Frame]`, `Frame(time_hours, swe, sca)` with
float32 copies; each frame's batch of steps runs in `asyncio.to_thread`) and
section 7's R7: a 2-triangle mesh, one 30-degree south facet and one
30-degree north, through ten April days: the south facet's SWE falls faster
under both corrections and equally under `none`; mass (I1) holds over the run.

Interface:

- `run.run(job, geometry, *, radiation, lat_deg, lon_deg, height_m)`: the
  site of the sun's beam, passed by the caller (the CLI takes it from the
  mesh's centre), and `radiation`, a `RadiationGeometry` (one receiver per
  triangle, in the triangles' order) whose evaluator gives b and sin e per
  triangle and sample (`docs/increments/46-hbv-snow.md`, *10.5*): dI = the
  mean over the step's samples of beam (b - max(0, sin e)). With
  `geometry_cls.open_sky` (no terrain) only self-shading is left.
- `Frame` has `time_hours` (float, hours since the start), `swe`, `sca`,
  `outflow` (cumulative since the start, mm) and `radiation_anomaly_mean`
  (the mean dI over the steps so far, W/m^2, 0 at the start), each float32
  `(n,)`, a copy. The first frame is the initial state at hour 0; then one
  every `frames.every_hours`. dI is computed under every correction, `none`
  included, so its mean is written whatever the mode.
- The steps between two frames run in one `asyncio.to_thread` call, looked up
  as `asyncio.to_thread` when called; the evaluator is reached through
  `_core.evaluate_instant`, looked up when called, with the job's threads.
"""

from __future__ import annotations

import asyncio
import datetime as dt
import importlib
import math
from collections.abc import Callable
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_array_equal

from tin_engine.resources import Resources

X0, Y0, Z = 500000.0, 6800000.0, 300.0
LAT, LON = 61.334, 15.0
TAN30 = math.tan(math.radians(30.0))
SIDE = 20.0


def two_facets() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """A 30-degree south facet and a 30-degree north facet, both with centroid z = 300 m."""

    def facet(dx: float, sign: float) -> list[list[float]]:
        plan = [(0.0, 0.0), (SIDE, 0.0), (0.0, SIDE)]
        cy = SIDE / 3
        # sign +1: z rises northward, so the facet faces south.
        return [[X0 + dx + x, Y0 + y, Z + sign * TAN30 * (y - cy)] for x, y in plan]

    vertices = np.array(facet(0.0, +1.0) + facet(100.0, -1.0))
    return vertices, np.array([[0, 1, 2], [3, 4, 5]])


SOUTH, NORTH = 0, 1
Mesh = tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]


@pytest.fixture(scope="module")
def mods() -> dict[str, ModuleType]:
    return {
        name: importlib.import_module(f"tin_engine.snow.{name}")
        for name in ("run", "params", "geometry")
    }


@pytest.fixture(scope="module")
def facets() -> Mesh:
    return two_facets()


def job(mods: dict[str, ModuleType], mode: str, *, threads: int = 1, **kw: Any) -> Any:
    p = mods["params"]
    spring = p.SyntheticSpring(start=dt.datetime(2027, 4, 10, tzinfo=dt.UTC), days=10)
    return p.SnowJob(
        correction=p.SunCorrection(mode=mode),
        resources=Resources(threads=threads, memory_budget=1 << 30),
        spring=spring,
        **kw,
    )


def sky_of(mesh: Mesh, horizon: int | None = None) -> Any:
    """The open-sky radiation geometry of `mesh`, or one with `horizon`
    (0.01 degrees) in every direction."""
    from tin_engine.radiation.frame import receivers_from_mesh
    from tin_engine.radiation.request import RadiationRequest

    geometry_cls = importlib.import_module("tin_engine.radiation.geometry").RadiationGeometry
    g = geometry_cls.open_sky(receivers_from_mesh(*mesh, "EPSG:25833"), RadiationRequest())
    if horizon is None:
        return g
    return geometry_cls(
        lat=g.lat, lon=g.lon, normal=g.normal, horizon=np.full_like(g.horizon, horizon),
        reach=g.reach, n_directions=g.n_directions, request=g.request,
    )  # fmt: skip


async def frames_of(
    mods: dict[str, ModuleType], the_job: Any, mesh: Mesh, radiation: Any = None
) -> list[Any]:
    geometry = mods["geometry"].triangle_geometry(*mesh, 0.0)
    run = mods["run"].run(
        the_job, geometry, radiation=sky_of(mesh) if radiation is None else radiation,
        lat_deg=LAT, lon_deg=LON, height_m=Z,
    )  # fmt: skip
    return [frame async for frame in run]


MODES = ("none", "melt_factor", "threshold")


@pytest.fixture(scope="module")
def runs(mods: dict[str, ModuleType], facets: Any) -> dict[str, list[Any]]:
    return {mode: asyncio.run(frames_of(mods, job(mods, mode), facets)) for mode in MODES}


class TestFrames:
    def test_a_frame_at_the_start_and_one_a_day(self, runs: dict[str, list[Any]]) -> None:
        for frames in runs.values():
            assert [f.time_hours for f in frames] == [24.0 * d for d in range(11)]

    def test_the_first_frame_is_the_initial_state(self, runs: dict[str, list[Any]]) -> None:
        first = runs["melt_factor"][0]
        assert_array_equal(first.swe, [600.0, 600.0])
        assert_array_equal(first.sca, [1.0, 1.0])
        assert_array_equal(first.outflow, [0.0, 0.0])
        assert_array_equal(first.radiation_anomaly_mean, [0.0, 0.0])

    def test_frames_are_float32_copies(self, runs: dict[str, list[Any]]) -> None:
        frames = runs["none"]
        for f in frames:
            for values in (f.swe, f.sca, f.outflow, f.radiation_anomaly_mean):
                assert values.dtype == np.float32 and values.shape == (2,)
        assert not np.shares_memory(frames[1].swe, frames[2].swe)

    def test_the_snow_melts_and_none_is_gone(self, runs: dict[str, list[Any]]) -> None:
        for frames in runs.values():
            assert (frames[-1].swe < 600.0).all() and (frames[-1].swe > 0.0).all()
            assert (frames[-1].outflow > 0.0).all()


class TestR7:
    def test_under_none_both_facets_melt_alike(self, runs: dict[str, list[Any]]) -> None:
        for f in runs["none"]:
            assert f.swe[SOUTH] == f.swe[NORTH]
            assert f.sca[SOUTH] == f.sca[NORTH]

    @pytest.mark.parametrize("mode", ["melt_factor", "threshold"])
    def test_the_south_facet_melts_faster_under_a_correction(
        self, runs: dict[str, list[Any]], mode: str
    ) -> None:
        frames = runs[mode]
        assert all(f.swe[SOUTH] <= f.swe[NORTH] for f in frames)
        assert frames[-1].swe[SOUTH] < frames[-1].swe[NORTH] - 1.0

    @pytest.mark.parametrize("mode", ["melt_factor", "threshold"])
    def test_the_correction_moves_melt_from_north_to_south(
        self, runs: dict[str, list[Any]], mode: str
    ) -> None:
        plain, corrected = runs["none"][-1].swe, runs[mode][-1].swe
        assert corrected[SOUTH] < plain[SOUTH]
        assert corrected[NORTH] > plain[NORTH]

    @pytest.mark.parametrize("mode", MODES)
    def test_the_south_facet_has_the_sun_and_the_north_the_shade(
        self, runs: dict[str, list[Any]], mode: str
    ) -> None:
        mean = runs[mode][-1].radiation_anomaly_mean
        assert mean[SOUTH] > 20.0 and mean[NORTH] < -20.0

    @pytest.mark.parametrize("mode", MODES)
    def test_mass_holds_over_the_run(self, runs: dict[str, list[Any]], mode: str) -> None:
        # I1, dry spring: swe + outflow = 600 mm. Scale: 600 mm in float32 (step 6.1e-5 mm).
        for f in runs[mode]:
            residual = 600.0 - f.swe.astype(np.float64) - f.outflow.astype(np.float64)
            assert np.abs(residual).max() <= 1e-3


class TestWiring:
    def test_a_flat_triangle_melts_the_same_under_every_correction(
        self, mods: dict[str, ModuleType]
    ) -> None:
        flat = (
            np.array([[X0, Y0, Z], [X0 + SIDE, Y0, Z], [X0, Y0 + SIDE, Z]]),
            np.array([[0, 1, 2]]),
        )
        out = {m: asyncio.run(frames_of(mods, job(mods, m), flat)) for m in MODES}
        for a, b, c in zip(out["none"], out["melt_factor"], out["threshold"], strict=True):
            assert_array_equal(a.swe, b.swe)
            assert_array_equal(a.swe, c.swe)
        # Not equal by doing nothing: melt is held as liquid water until it passes
        # a tenth of the ice, so SWE first falls after about four days, then does.
        assert out["none"][-1].swe[0] < 600.0
        # Open sky, flat: b is max(0, sin e) bit for bit, so dI is exactly 0.
        assert all(f.radiation_anomaly_mean[0] == 0.0 for f in out["none"])

    def test_the_frames_do_not_depend_on_the_thread_count(
        self, mods: dict[str, ModuleType], facets: Any, runs: dict[str, list[Any]]
    ) -> None:
        many = asyncio.run(frames_of(mods, job(mods, "threshold", threads=4), facets))
        for a, b in zip(runs["threshold"], many, strict=True):
            assert_array_equal(a.swe, b.swe)
            assert_array_equal(a.sca, b.sca)

    def test_the_distribution_and_initial_swe_come_from_the_job(
        self, mods: dict[str, ModuleType], facets: Any
    ) -> None:
        p = mods["params"]
        one_bin = job(
            mods, "none", distribution=p.SnowDistribution.lognormal(0.0, 1), initial_swe_mm=5.0
        )
        frames = asyncio.run(frames_of(mods, one_bin, facets))
        assert_array_equal(frames[0].swe, [5.0, 5.0])
        # One bin: a triangle is all white or all bare, never in between.
        assert all(set(f.sca.tolist()) <= {0.0, 1.0} for f in frames)
        assert frames[-1].sca.tolist() == [0.0, 0.0]

    async def test_each_frames_steps_run_off_the_event_loop(
        self, mods: dict[str, ModuleType], facets: Any, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        calls: list[Callable[..., Any]] = []
        original = asyncio.to_thread

        async def counting(func: Callable[..., Any], /, *args: Any, **kwargs: Any) -> Any:
            calls.append(func)
            return await original(func, *args, **kwargs)

        monkeypatch.setattr(asyncio, "to_thread", counting)
        frames = await frames_of(mods, job(mods, "melt_factor"), facets)
        assert len(frames) == 11
        assert len(calls) >= 10


class TestShadows:
    """The sun correction from the evaluator (section 10.5)."""

    def test_the_anomaly_is_the_beam_less_the_horizontal_one(
        self, mods: dict[str, ModuleType], facets: Mesh
    ) -> None:
        from tin_engine.snow.forcing import sample_times, sun_and_beam

        evaluate_instant = importlib.import_module("tin_engine.radiation.evaluate").evaluate_instant
        beam_anomaly = importlib.import_module("tin_engine.snow.forcing").beam_anomaly

        the_job = job(mods, "melt_factor")
        frames = asyncio.run(frames_of(mods, the_job, facets))
        sky = sky_of(facets)
        want = np.zeros(2)
        for k in range(24):  # the first frame's 24 hourly steps
            times = sample_times(the_job.spring, k)
            _, beam = sun_and_beam(times, LAT, LON, Z)
            f = evaluate_instant(sky, times, threads=1)
            want += beam_anomaly(f.b, f.sin_e, beam)
        # Scale: dI up to about 300 W/m^2, stored as float32 (about 3e-5).
        assert frames[1].radiation_anomaly_mean == pytest.approx(want / 24, abs=1e-3)

    def test_a_horizon_above_the_sun_shades_both_facets_alike(
        self, mods: dict[str, ModuleType], facets: Mesh
    ) -> None:
        # 60 degrees all round: above the April sun at 61 N (at most about 40).
        frames = asyncio.run(
            frames_of(mods, job(mods, "melt_factor"), facets, sky_of(facets, 6000))
        )
        mean = frames[-1].radiation_anomaly_mean
        assert mean[SOUTH] < -20.0 and mean[NORTH] < -20.0
        # Both get -beam sin e; the facets are 100 m apart, so their suns differ a little.
        assert mean[SOUTH] == pytest.approx(mean[NORTH], rel=1e-3)
        open_sky = asyncio.run(frames_of(mods, job(mods, "melt_factor"), facets))
        assert (mean < open_sky[-1].radiation_anomaly_mean).all()
        assert frames[-1].swe[SOUTH] > open_sky[-1].swe[SOUTH]

    def test_the_thread_count_reaches_the_evaluator(
        self, mods: dict[str, ModuleType], facets: Mesh, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import tin_engine._core as core

        seen: list[int] = []
        real = core.evaluate_instant

        def spy(*args: Any, **kwargs: Any) -> Any:
            seen.append(kwargs.get("threads", -1))
            return real(*args, **kwargs)

        monkeypatch.setattr(core, "evaluate_instant", spy)
        asyncio.run(frames_of(mods, job(mods, "none", threads=3), facets))
        assert seen and set(seen) == {3}
