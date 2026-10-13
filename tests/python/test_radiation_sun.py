"""The sun's position through the binding, and UTC at the Python boundary.

`docs/increments/46-hbv-snow.md`, *4.4 The sun's position: the radiation
design's TR1, not a second copy*, and the radiation design's verification of
the sun (`docs/increments/terrain-radiation.md`, R2, I1, I2): NREL's test case
through the binding; 10^4 random (time, place, height, weather) cases against
pvlib's independent SPA within 1e-6 degrees (pvlib is a test-only dependency,
as Ola ruled); `sun_positions` against NOAA's solar calculator at the synthetic
domain's centre; and `radiation/time.py`, which turns timezone-aware times into
UTC milliseconds and refuses a naive one.

Interface assumed (names chosen here and listed in the handback):

- `_core.delta_t_seconds(time_utc_ms: int) -> float`: the default delta T.
- `_core.solar_position(time_utc_ms, lat_deg, lon_deg, height_m, delta_t_s,
  pressure_hpa, temperature_c)`: SPA point by point, returning an object with
  `e0_deg` (topocentric elevation, no refraction), `e_deg` (with refraction,
  guarded below -1.8 degrees) and `azimuth_deg` (true, from north, clockwise).
- `_core.sun_positions(times_utc_ms, lat, lon, height_m)`: an int64 array of
  UTC milliseconds in; `elevation` (no refraction) and `azimuth` in degrees,
  `distance` in AU out, float64 arrays of the input's length, with the default
  delta T.
- `tin_engine.radiation.time.to_utc_ms(times)`: an iterable of timezone-aware
  `datetime`s or ISO 8601 strings with an offset, or a `datetime64` array
  (taken as UTC), to a 1-D int64 array of milliseconds since the Unix epoch. A
  naive `datetime` or a string without an offset is a `ValueError` saying the
  time zone is missing.
"""

from __future__ import annotations

import datetime as dt
import importlib
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest

NREL_MS = 1066419030000  # 2003-10-17 12:30:30 at UTC-7
UTC = dt.UTC
MINUS_7 = dt.timezone(dt.timedelta(hours=-7))
I64 = npt.NDArray[np.int64]


@pytest.fixture(scope="module")
def core() -> Any:
    import tin_engine._core as core

    for name in ("delta_t_seconds", "solar_position", "sun_positions"):
        assert hasattr(core, name), f"tin_engine._core has no {name}"
    return core


@pytest.fixture(scope="module")
def rtime() -> ModuleType:
    return importlib.import_module("tin_engine.radiation.time")


@pytest.fixture(scope="module")
def pvlib_spa() -> ModuleType:
    try:
        return importlib.import_module("pvlib.spa")
    except ImportError:
        pytest.fail(
            "pvlib is the test-only second SPA implementation; "
            "install the dev extra (pip install -e .[dev])"
        )


def day_ms(day: str, start_h: int, hours: int, step_s: int) -> I64:
    start = np.datetime64(f"{day}T{start_h:02d}:00", "ms").astype(np.int64)
    return start + np.arange(0, hours * 3600 + 1, step_s, dtype=np.int64) * 1000


class TestNrel:
    """I1 through the binding: NREL's published test case."""

    def test_the_reference_path_gives_the_published_zenith_and_azimuth(self, core: Any) -> None:
        s = core.solar_position(NREL_MS, 39.742476, -105.1786, 1830.14, 67.0, 820.0, 11.0)
        assert 90.0 - s.e_deg == pytest.approx(50.11162, abs=5e-6)
        assert s.azimuth_deg == pytest.approx(194.34024, abs=5e-6)
        assert s.e0_deg == pytest.approx(39.872046, abs=5e-7)

    def test_the_sun_earth_distance_at_the_nrel_case(self, core: Any) -> None:
        # Published R = 0.9965422974 AU; delta T for 2003 is 64-65 s, against
        # the case's 67 s, which moves R by under 1e-10 AU.
        p = core.sun_positions(np.array([NREL_MS], dtype=np.int64), 39.742476, -105.1786, 1830.14)
        assert p.distance[0] == pytest.approx(0.9965422974, abs=5e-10)


class TestAgainstPvlib:
    """10^4 random cases against pvlib's SPA, the same delta T, within 1e-6 degrees."""

    N = 10_000

    @pytest.fixture(scope="class")
    def cases(self, core: Any, pvlib_spa: ModuleType) -> dict[str, Any]:
        rng = np.random.default_rng(46_2)
        n = self.N
        t = rng.integers(-2208988800000, 4102444800000, n)  # 1900 to 2100
        lat = rng.uniform(-89.0, 89.0, n)
        lon = rng.uniform(-180.0, 180.0, n)
        h = rng.uniform(0.0, 3000.0, n)
        p = rng.uniform(500.0, 1050.0, n)
        temp = rng.uniform(-40.0, 40.0, n)
        delta_t = np.array([core.delta_t_seconds(int(v)) for v in t])
        ours = [
            core.solar_position(int(t[i]), lat[i], lon[i], h[i], delta_t[i], p[i], temp[i])
            for i in range(n)
        ]
        _, _, e, e0, az, _ = pvlib_spa.solar_position(
            t / 1000.0, lat, lon, h, p, temp, delta_t, 0.5667
        )
        return {
            "e0": np.array([s.e0_deg for s in ours]),
            "e": np.array([s.e_deg for s in ours]),
            "az": np.array([s.azimuth_deg for s in ours]),
            "pv_e0": e0,
            "pv_e": e,
            "pv_az": az,
        }

    def test_the_elevation_without_refraction(self, cases: dict[str, Any]) -> None:
        gap = np.abs(cases["e0"] - cases["pv_e0"])
        assert gap.max() <= 1e-6, f"worst case {int(gap.argmax())}: {gap.max():.3g} degrees"

    def test_the_azimuth(self, cases: dict[str, Any]) -> None:
        # Near the zenith the azimuth is ill-conditioned; scaled by cos(e0) the
        # gap is the arc between the two sun directions.
        d = np.abs(cases["az"] - cases["pv_az"]) % 360.0
        d = np.minimum(d, 360.0 - d) * np.cos(np.radians(cases["e0"]))
        assert d.max() <= 1e-6, f"worst case {int(d.argmax())}: {d.max():.3g} degrees"

    def test_the_refracted_elevation_where_both_refract(self, cases: dict[str, Any]) -> None:
        # pvlib, as SPA, refracts only from -0.8333 degrees up; the guard below
        # that is the C++ suite's (A1-A2).
        both = cases["e0"] >= -0.8333
        assert both.sum() > 1000  # roughly half the cases are by day
        gap = np.abs(cases["e"][both] - cases["pv_e"][both])
        assert gap.max() <= 1e-6

    def test_below_the_guard_the_refraction_is_held(self, cases: dict[str, Any]) -> None:
        low = cases["e0"] < -1.8
        assert low.sum() > 1000
        refraction = cases["e"][low] - cases["e0"][low]
        assert np.all(refraction > 0.0)
        assert np.all(np.isfinite(refraction))


class TestSunPositions:
    """The vectorised call against NOAA's calculator (`46-probes/probe_sun.out`)."""

    @pytest.mark.parametrize(("day", "highest"), [("2027-04-15", 38.44), ("2027-06-21", 52.10)])
    def test_the_days_highest_sun_at_the_synthetic_domain(
        self, core: Any, day: str, highest: float
    ) -> None:
        # 61.334 N 15 E, one-minute steps from 09:00 to 13:00 UTC. NOAA's
        # formulas are good to about 0.01 degrees; the pvlib SPA gives 38.4373
        # and 52.1017.
        p = core.sun_positions(day_ms(day, 9, 4, 60), 61.334, 15.0, 0.0)
        top = int(np.argmax(p.elevation))
        assert p.elevation[top] == pytest.approx(highest, abs=0.02)
        assert p.azimuth[top] == pytest.approx(180.0, abs=0.2)

    def test_it_agrees_with_the_reference_path_at_the_default_delta_t(self, core: Any) -> None:
        times = day_ms("2027-04-15", 0, 24, 900)
        p = core.sun_positions(times, 61.334, 15.0, 1000.0)
        for i, t in enumerate(times.tolist()):
            ref = core.solar_position(
                t, 61.334, 15.0, 1000.0, core.delta_t_seconds(t), 1010.0, 10.0
            )
            assert p.elevation[i] == pytest.approx(ref.e0_deg, abs=1e-9)
            assert p.azimuth[i] == pytest.approx(ref.azimuth_deg, abs=1e-9)

    def test_the_arrays_match_the_input(self, core: Any) -> None:
        times = day_ms("2027-06-21", 0, 24, 3600)
        p = core.sun_positions(times, 61.334, 15.0, 0.0)
        for name in ("elevation", "azimuth", "distance"):
            a = np.asarray(getattr(p, name))
            assert a.dtype == np.float64, name
            assert a.shape == times.shape, name
        # Midnight sun is not reached at 61 N: the sun is below the horizon at
        # midnight and above it at noon.
        assert p.elevation[0] < 0.0 < p.elevation[12]

    def test_no_timestamps_give_empty_arrays(self, core: Any) -> None:
        p = core.sun_positions(np.array([], dtype=np.int64), 61.334, 15.0, 0.0)
        for name in ("elevation", "azimuth", "distance"):
            assert np.asarray(getattr(p, name)).shape == (0,), name


class TestDeltaT:
    def test_2026_is_the_polynomial_not_the_legacy_defect(self, core: Any) -> None:
        # 2026-01-15T12:00Z: Espenak and Meeus give 75.1 s; the legacy's
        # parabola with /20 gave 3375 s.
        assert core.delta_t_seconds(1768478400000) == pytest.approx(75.10, abs=0.1)


class TestUtc:
    """I2: every timestamp is UTC; a naive one is refused."""

    def test_an_offset_datetime_becomes_utc_milliseconds(self, rtime: ModuleType) -> None:
        local = dt.datetime(2003, 10, 17, 12, 30, 30, tzinfo=MINUS_7)
        out = rtime.to_utc_ms([local])
        assert out.dtype == np.int64
        assert out.tolist() == [NREL_MS]

    def test_the_same_instant_in_any_zone_is_the_same_millisecond(self, rtime: ModuleType) -> None:
        oslo_summer = dt.timezone(dt.timedelta(hours=2))
        times = [
            dt.datetime(2027, 6, 21, 12, 0, tzinfo=UTC),
            dt.datetime(2027, 6, 21, 14, 0, tzinfo=oslo_summer),
            dt.datetime(2027, 6, 21, 5, 0, tzinfo=MINUS_7),
        ]
        out = rtime.to_utc_ms(times).tolist()
        assert out == [out[0]] * 3
        assert out[0] == int(np.datetime64("2027-06-21T12:00", "ms").astype(np.int64))

    def test_iso_strings_with_an_offset(self, rtime: ModuleType) -> None:
        out = rtime.to_utc_ms(["2003-10-17T12:30:30-07:00", "2003-10-17T19:30:30Z"])
        assert out.tolist() == [NREL_MS, NREL_MS]

    def test_milliseconds_are_kept_before_and_after_the_epoch(self, rtime: ModuleType) -> None:
        out = rtime.to_utc_ms(
            [
                dt.datetime(1969, 12, 31, 23, 59, 59, 999000, tzinfo=UTC),
                dt.datetime(1970, 1, 1, 0, 0, 0, 1000, tzinfo=UTC),
                dt.datetime(2003, 10, 17, 19, 30, 30, 123000, tzinfo=UTC),
            ]
        )
        assert out.tolist() == [-1, 1, NREL_MS + 123]

    def test_a_datetime64_array_is_utc_and_exact(self, rtime: ModuleType) -> None:
        a = np.array(["1969-12-31T23:59:59.999", "2003-10-17T19:30:30.123"], dtype="datetime64[ms]")
        assert rtime.to_utc_ms(a).tolist() == a.astype(np.int64).tolist()

    def test_no_times_give_an_empty_array(self, rtime: ModuleType) -> None:
        out = rtime.to_utc_ms([])
        assert out.dtype == np.int64
        assert out.shape == (0,)

    def test_a_naive_datetime_is_refused(self, rtime: ModuleType) -> None:
        with pytest.raises(ValueError, match="time zone"):
            rtime.to_utc_ms([dt.datetime(2027, 4, 15, 12, 0)])

    def test_a_naive_datetime_among_aware_ones_is_refused(self, rtime: ModuleType) -> None:
        with pytest.raises(ValueError, match="time zone"):
            rtime.to_utc_ms(
                [dt.datetime(2027, 4, 15, 12, 0, tzinfo=UTC), dt.datetime(2027, 4, 15, 13, 0)]
            )

    def test_a_string_without_an_offset_is_refused(self, rtime: ModuleType) -> None:
        with pytest.raises(ValueError, match="time zone"):
            rtime.to_utc_ms(["2027-04-15T12:00:00"])
