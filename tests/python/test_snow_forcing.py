"""The synthetic spring's forcing and the sun correction in NumPy: `snow/forcing.py` and R2.

`docs/increments/46-hbv-snow.md`, *4.3 The sun correction*, *4.5 The
synthetic case* and *10.5 How the lit factor enters the sun correction*:

    dI_i = (1/S) sum_s beam_s (b_is - max(0, sin e_s))

with b_is the beam factor on triangle i (here self-shading only, max(0,
n_i . u_s), as under an open sky), u_s the sun's direction in east-north-up
at the midpoint of sub-interval s of the step (S equal sub-intervals),
beam_s = I0 (1/R)^2 Psi_a^(1/sin e) with I0 = 1368 W/m^2, Psi_a = 0.75, and
0 with the sun at or below the horizon. The temperature at 300 m runs
linearly from -3 to +14 degrees C over the run, plus a +-4 degree cosine
peaking at 14:00 mean solar time, minus 0.6 degrees per 100 m above 300 m,
evaluated at each step's midpoint.

R2 (through `_core.sun_positions` and this module's `beam_anomaly`,
2027-04-15 at 61.334 N 15 E, 96 samples of 15 minutes): +68.9, -113.5,
-12.6, -12.1 W/m^2 within 0.5. Only `TestSunAndBeam` needs the sun binding
(PR 2). The combiner's own edges are L2 in `test_radiation_evaluate.py`.

Interface:

- `forcing.beam_anomaly(b (S, n), sin_e (S, n), beam (S,)) -> (n,)`.
- `forcing.beam(elevation_deg, distance_au) -> array`, elementwise.
- `forcing.sun_and_beam(times_utc_ms, lat_deg, lon_deg, height_m) -> (sun
  (S, 3), beam (S,))`, sun = (cos e sin A, cos e cos A, sin e) from
  `_core.sun_positions`.
- `forcing.sample_times(spring, k) -> int64 (S,)`: UTC milliseconds of the
  S sub-interval midpoints of step k.
- `forcing.temperature(spring, k, heights, lon_deg) -> (n,)`: at step k's
  midpoint, mean solar time = UTC + lon/15 hours.
- `params.SyntheticSpring` fields `start` (timezone-aware), `days`,
  `step_minutes`, `samples_per_step`; defaults 2027-04-01T00:00Z, 91, 60, 4.
"""

from __future__ import annotations

import datetime as dt
import importlib
import math
from types import ModuleType

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_allclose, assert_array_equal

F64 = npt.NDArray[np.float64]
LAT, LON = 61.334, 15.0  # the synthetic domain's centre
START_MS = int(dt.datetime(2027, 4, 1, tzinfo=dt.UTC).timestamp() * 1000)


@pytest.fixture(scope="module")
def forcing() -> ModuleType:
    return importlib.import_module("tin_engine.snow.forcing")


@pytest.fixture(scope="module")
def params() -> ModuleType:
    return importlib.import_module("tin_engine.snow.params")


def plane(slope_deg: float, aspect_deg: float) -> F64:
    s, a = math.radians(slope_deg), math.radians(aspect_deg)
    return np.array([math.sin(s) * math.sin(a), math.sin(s) * math.cos(a), math.cos(s)])


PLANES = np.array([plane(30, 180), plane(30, 0), plane(30, 90), plane(30, 270), plane(0, 0)])


class TestBeam:
    def test_the_sun_at_30_degrees_at_1_au(self, forcing: ModuleType) -> None:
        assert forcing.beam(np.array([30.0]), np.array([1.0]))[0] == pytest.approx(769.5, rel=1e-12)

    def test_the_zenith_sun(self, forcing: ModuleType) -> None:
        assert forcing.beam(np.array([90.0]), np.array([1.0]))[0] == pytest.approx(
            1026.0, rel=1e-12
        )

    def test_the_distance_scales_by_its_inverse_square(self, forcing: ModuleType) -> None:
        b = forcing.beam(np.array([30.0]), np.array([0.99]))[0]
        assert b == pytest.approx(769.5 / 0.99**2, rel=1e-12)

    @pytest.mark.parametrize("elevation", [0.0, -0.01, -30.0, -90.0])
    def test_no_beam_at_or_below_the_horizon(self, forcing: ModuleType, elevation: float) -> None:
        assert forcing.beam(np.array([elevation]), np.array([1.0]))[0] == 0.0

    def test_a_grazing_sun_gives_a_tiny_finite_beam(self, forcing: ModuleType) -> None:
        b = forcing.beam(np.array([0.01]), np.array([1.0]))[0]
        assert math.isfinite(b) and 0.0 <= b < 1e-100


class TestSampleTimes:
    def test_four_midpoints_of_the_first_hour(
        self, forcing: ModuleType, params: ModuleType
    ) -> None:
        times = forcing.sample_times(params.SyntheticSpring(), 0)
        assert times.dtype == np.int64
        assert_array_equal(times, START_MS + np.array([7.5, 22.5, 37.5, 52.5]) * 60_000)

    def test_step_k_is_k_steps_later(self, forcing: ModuleType, params: ModuleType) -> None:
        spring = params.SyntheticSpring()
        assert_array_equal(
            forcing.sample_times(spring, 2183), forcing.sample_times(spring, 0) + 2183 * 3_600_000
        )

    def test_a_15_minute_step_with_one_sample(
        self, forcing: ModuleType, params: ModuleType
    ) -> None:
        spring = params.SyntheticSpring(step_minutes=15, samples_per_step=1)
        assert_array_equal(forcing.sample_times(spring, 3), [START_MS + 52.5 * 60_000])


class TestTemperature:
    """The ramp, the daily cosine and the lapse, at the step's midpoint."""

    @staticmethod
    def expected(k: int, z: float, lon: float) -> float:
        hours = k + 0.5  # midpoint of hourly step k, hours since 2027-04-01T00Z
        solar = (hours + lon / 15.0) % 24.0
        at_300 = (
            -3.0 + 17.0 * hours / (91 * 24) + 4.0 * math.cos(2 * math.pi * (solar - 14.0) / 24.0)
        )
        return at_300 - 0.006 * (z - 300.0)

    @pytest.mark.parametrize("k", [0, 12, 13, 500, 2183])
    def test_it_follows_section_4_5(self, forcing: ModuleType, params: ModuleType, k: int) -> None:
        heights = np.array([300.0, 1300.0, 234.0, 1924.4])
        t = forcing.temperature(params.SyntheticSpring(), k, heights, LON)
        assert_allclose(t, [self.expected(k, z, LON) for z in heights], atol=1e-12)

    def test_the_peak_is_at_13_utc_at_15_east(
        self, forcing: ModuleType, params: ModuleType
    ) -> None:
        spring = params.SyntheticSpring()
        day = [forcing.temperature(spring, k, np.array([300.0]), LON)[0] for k in range(24)]
        # Midpoints 12:30 and 13:30 UTC straddle 13:00 UTC, 14:00 mean solar time.
        assert day[12] == pytest.approx(day[13], abs=0.01)
        assert max(day) in (day[12], day[13])

    def test_the_peak_moves_with_the_longitude(
        self, forcing: ModuleType, params: ModuleType
    ) -> None:
        spring = params.SyntheticSpring()
        bygdin = [forcing.temperature(spring, k, np.array([300.0]), 8.4852)[0] for k in range(24)]
        assert int(np.argmax(bygdin)) == 13  # 14:00 solar is 13:26 UTC at 8.49 E

    def test_1000_m_higher_is_6_degrees_colder(
        self, forcing: ModuleType, params: ModuleType
    ) -> None:
        t = forcing.temperature(params.SyntheticSpring(), 40, np.array([500.0, 1500.0]), LON)
        assert t[0] - t[1] == pytest.approx(6.0, abs=1e-12)


class TestSunAndBeam:
    """Through `_core.sun_positions` (PR 2)."""

    def test_directions_are_unit_and_the_night_has_no_beam(self, forcing: ModuleType) -> None:
        midnight = int(dt.datetime(2027, 4, 15, 23, 0, tzinfo=dt.UTC).timestamp() * 1000)
        noon = int(dt.datetime(2027, 4, 15, 11, 0, tzinfo=dt.UTC).timestamp() * 1000)
        sun, beam = forcing.sun_and_beam(np.array([midnight, noon]), LAT, LON, 0.0)
        assert sun.shape == (2, 3)
        assert_allclose(np.linalg.norm(sun, axis=1), 1.0, atol=1e-12)
        assert sun[0, 2] < 0.0 and beam[0] == 0.0
        assert sun[1, 2] > 0.5 and 600.0 < beam[1] < 900.0
        assert sun[1, 1] < 0.0  # the noon sun is to the south

    def test_r2_day_means_on_2027_04_15(self, forcing: ModuleType, params: ModuleType) -> None:
        spring = params.SyntheticSpring(start=dt.datetime(2027, 4, 15, tzinfo=dt.UTC), days=1)
        total = np.zeros(len(PLANES))
        for k in range(24):
            sun, beam = forcing.sun_and_beam(forcing.sample_times(spring, k), LAT, LON, 0.0)
            b = np.maximum(sun @ PLANES.T, 0.0) * (sun[:, 2:] > 0.0)
            total += forcing.beam_anomaly(b, np.repeat(sun[:, 2:], len(PLANES), axis=1), beam)
        mean = total / 24
        # NOAA's calculator in the probe, TR1 here: about 0.01 degrees apart.
        assert_allclose(mean[:4], [68.9, -113.5, -12.6, -12.1], atol=0.5)
        assert mean[4] == 0.0
