"""The instant evaluator through the binding and its synchronous Python call.

`docs/increments/46-hbv-snow.md`, *10.4 Which radiation steps, and how much of
each* (TR5's `instant` mode, its binding and a synchronous call, reading the
receiver set and horizon table as arrays; `RadiationGeometry` with TR2's and
TR3a's fields), *10.5 How the lit factor enters the sun correction* (the
combiner: dI = mean over samples of beam (b - max(0, sin e))) and *10.8 Red
suites for `@tester`* (L2, and I7 and I8 as they cross);
`docs/increments/terrain-radiation.md`, R5 and R10. The C++ suite
(`tests/cpp/unit/test_radiation_evaluate.cpp`) holds the lit test, I7 in full
and the sunrise; this file checks what crosses and the Python surface.

Interface assumed (names chosen here and listed in the handback):

- `_core.evaluate_instant(lat_deg, lon_deg, normal, horizon, times_utc_ms, *,
  limb_deg=0.0, pressure_hpa=1010.0, temperature_c=10.0, threads=0) ->
  InstantFactors`: `normal` (n, 3) float64, `horizon` (n, N) int16 in 0.01
  degrees, `times_utc_ms` int64; `InstantFactors` has `lit` uint8 and `cos_i`,
  `b`, `sin_e` float32, each (T, n). A `ValueError` naming the field or
  "receiver i"; a horizon that is not int16 is a `TypeError` naming int16.
  Releases the GIL.
- `radiation.geometry.RadiationGeometry`, frozen: `lat`, `lon` (n,),
  `normal` (n, 3), `horizon` (n, N) int16, `reach` (n, N) uint8,
  `n_directions`, `request` (a `RadiationRequest`).
  `RadiationGeometry.open_sky(frame, request)`: the frame's lat, lon and
  normal with a horizon of -9000 (-90 degrees, no terrain) and reach 0 in every
  direction.
- `radiation.evaluate.evaluate_instant(geometry, times_utc_ms, *,
  pressure_hpa=1010.0, temperature_c=10.0, threads=0) -> InstantFactors`, the
  limb from `geometry.request.limb` (0 for "centre", 0.2667 for "upper").
- `snow.forcing.beam_anomaly(b, sin_e, beam_wm2) -> (n,)`: b and sin_e (S, n),
  beam (S,): mean over s of beam_s (b_s - max(0, sin_e_s)), float64.
"""

from __future__ import annotations

import datetime as dt
import importlib
import math
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_array_equal

from gil_probe import ticks_during
from tin_engine.radiation.request import RadiationRequest

F64 = npt.NDArray[np.float64]
RELEASED_TICKS = 20  # as test_core_snow.py: a released call lets the ticker run
BYGDIN = (61.3694, 8.4852)
UPPER_LIMB = 0.2667
X0, Y0, Z = 500000.0, 6800000.0, 300.0


def day_ms(y: int, m: int, d: int) -> int:
    return int(dt.datetime(y, m, d, tzinfo=dt.UTC).timestamp() * 1000)


#: 96 samples of 15 minutes on 2027-04-15, UTC.
APRIL_15 = day_ms(2027, 4, 15) + 900_000 * np.arange(96, dtype=np.int64)


@pytest.fixture(scope="module")
def core() -> Any:
    import tin_engine._core as core

    assert hasattr(core, "evaluate_instant"), "tin_engine._core has no evaluate_instant"
    return core


def receivers(normals: list[list[float]], row: npt.ArrayLike) -> dict[str, Any]:
    n = len(normals)
    return {
        "lat_deg": np.full(n, BYGDIN[0]),
        "lon_deg": np.full(n, BYGDIN[1]),
        "normal": np.array(normals, dtype=np.float64),
        "horizon": np.tile(np.asarray(row, dtype=np.int16), (n, 1)),
    }


OPEN = np.full(64, -9000, dtype=np.int16)
UP, EAST, NORTH = [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]


def module(name: str) -> Any:
    """A module this increment adds, imported when a test runs (so a missing
    one fails that test, not the collection)."""
    return importlib.import_module(f"tin_engine.{name}")


class TestBinding:
    def test_shapes_and_types_time_major(self, core: Any) -> None:
        got = core.evaluate_instant(**receivers([UP, EAST, NORTH], OPEN), times_utc_ms=APRIL_15)
        assert got.lit.dtype == np.uint8
        for a in (got.cos_i, got.b, got.sin_e):
            assert a.dtype == np.float32
        for a in (got.lit, got.cos_i, got.b, got.sin_e):
            assert a.shape == (96, 3)

    def test_the_sun_matches_the_reference_path(self, core: Any) -> None:
        # I7 as it crosses: elevation from the up receiver, azimuth from the walls.
        got = core.evaluate_instant(**receivers([UP, EAST, NORTH], OPEN), times_utc_ms=APRIL_15)
        checked = 0
        for s, t in enumerate(APRIL_15):
            dt_s = core.delta_t_seconds(int(t))
            ref = core.solar_position(int(t), *BYGDIN, 0.0, dt_s, 1010.0, 10.0)
            if ref.e_deg <= -1.0:
                continue
            e = math.degrees(math.asin(float(got.sin_e[s, 0])))
            a = math.degrees(math.atan2(float(got.cos_i[s, 1]), float(got.cos_i[s, 2]))) % 360
            assert e == pytest.approx(ref.e_deg, abs=0.003)
            gap = abs(a - ref.azimuth_deg) % 360
            assert min(gap, 360 - gap) <= 0.003
            checked += 1
        assert checked > 40

    def test_pressure_zero_is_no_refraction(self, core: Any) -> None:
        got = core.evaluate_instant(
            **receivers([UP], OPEN), times_utc_ms=APRIL_15, pressure_hpa=0.0, temperature_c=10.0
        )
        for s, t in enumerate(APRIL_15):
            ref = core.solar_position(int(t), *BYGDIN, 0.0, core.delta_t_seconds(int(t)), 0.0, 10.0)
            assert ref.e_deg == ref.e0_deg
            e = math.degrees(math.asin(float(got.sin_e[s, 0])))
            assert e == pytest.approx(ref.e0_deg, abs=0.003)

    def test_the_limb_raises_the_sun(self, core: Any) -> None:
        # Sunrise behind a 0 horizon: the upper limb clears it before the centre.
        flat = receivers([EAST], np.zeros(64, dtype=np.int16))
        minutes = day_ms(2027, 4, 15) + 60_000 * np.arange(3 * 60, 5 * 60, dtype=np.int64)
        centre = core.evaluate_instant(**flat, times_utc_ms=minutes).lit[:, 0]
        upper = core.evaluate_instant(**flat, times_utc_ms=minutes, limb_deg=UPPER_LIMB).lit[:, 0]
        assert centre.any() and not centre.all()
        assert (upper >= centre).all()
        assert int(np.argmax(upper)) < int(np.argmax(centre))

    def test_the_thread_count_changes_nothing(self, core: Any) -> None:
        rng = np.random.default_rng(46)
        n = 500
        slope, aspect = rng.uniform(0, 1.3, n), rng.uniform(0, 2 * np.pi, n)
        args = {
            "lat_deg": rng.uniform(58, 71, n),
            "lon_deg": rng.uniform(4, 31, n),
            "normal": np.column_stack(
                [np.sin(slope) * np.sin(aspect), np.sin(slope) * np.cos(aspect), np.cos(slope)]
            ),
            "horizon": rng.integers(-500, 2500, (n, 64)).astype(np.int16),
            "times_utc_ms": APRIL_15,
        }
        one = core.evaluate_instant(**args, threads=1)
        eight = core.evaluate_instant(**args, threads=8)
        for name in ("lit", "cos_i", "b", "sin_e"):
            assert_array_equal(getattr(one, name), getattr(eight, name))
        assert one.lit.any() and not one.lit.all()

    def test_it_releases_the_gil(self, core: Any) -> None:
        n = 100_000
        args = receivers([UP], OPEN)
        args = {k: np.repeat(v, n, axis=0) for k, v in args.items()}

        def run() -> Any:
            return core.evaluate_instant(**args, times_utc_ms=APRIL_15, threads=1)

        _, ticks, elapsed = ticks_during(run)
        assert elapsed >= 0.05, f"took only {elapsed:.3f}s: too fast to see; enlarge the job"
        assert ticks >= RELEASED_TICKS, f"only {ticks} ticks in {elapsed:.3f}s: it holds the GIL"


class TestRefusals:
    def test_a_horizon_of_the_wrong_width(self, core: Any) -> None:
        args = receivers([UP, UP], OPEN)
        args["horizon"] = np.zeros((2, 63), dtype=np.int16)
        with pytest.raises(ValueError, match="horizon"):
            core.evaluate_instant(**args, times_utc_ms=APRIL_15)

    def test_a_normal_of_the_wrong_length(self, core: Any) -> None:
        args = receivers([UP, UP], OPEN)
        args["normal"] = np.zeros((1, 3))
        with pytest.raises(ValueError, match="normal"):
            core.evaluate_instant(**args, times_utc_ms=APRIL_15)

    def test_a_latitude_past_the_pole_limit(self, core: Any) -> None:
        args = receivers([UP, UP], OPEN)
        args["lat_deg"] = np.array([61.0, 89.95])
        with pytest.raises(ValueError, match="receiver 1"):
            core.evaluate_instant(**args, times_utc_ms=APRIL_15)

    def test_a_normal_that_is_not_finite(self, core: Any) -> None:
        args = receivers([UP, [0.0, math.nan, 1.0]], OPEN)
        with pytest.raises(ValueError, match="receiver 1"):
            core.evaluate_instant(**args, times_utc_ms=APRIL_15)

    def test_a_float_horizon(self, core: Any) -> None:
        args = receivers([UP], OPEN)
        args["horizon"] = args["horizon"].astype(np.float64)
        with pytest.raises(TypeError, match="int16"):
            core.evaluate_instant(**args, times_utc_ms=APRIL_15)


def facet_frame() -> Any:
    """The receiver frame of one 30-degree south facet and one flat triangle."""
    from tin_engine.radiation.frame import receivers_from_mesh

    t = math.tan(math.radians(30.0))
    pts = [[X0, Y0, Z], [X0 + 20, Y0, Z], [X0, Y0 + 20, Z + 20 * t]]
    pts += [[X0 + 100, Y0, Z], [X0 + 120, Y0, Z], [X0 + 100, Y0 + 20, Z]]
    return receivers_from_mesh(np.array(pts), np.array([[0, 1, 2], [3, 4, 5]]), "EPSG:25833")


class TestRadiationGeometry:
    def test_open_sky_carries_the_frame_and_no_terrain(self) -> None:
        geometry_cls = module("radiation.geometry").RadiationGeometry

        frame, request = facet_frame(), RadiationRequest()
        g = geometry_cls.open_sky(frame, request)
        assert_array_equal(g.lat, frame.lat)
        assert_array_equal(g.lon, frame.lon)
        assert_array_equal(g.normal, frame.normal)
        assert g.n_directions == 64 and g.request == request
        assert g.horizon.dtype == np.int16 and g.horizon.shape == (2, 64)
        assert (g.horizon == -9000).all()
        assert g.reach.dtype == np.uint8 and g.reach.shape == (2, 64)
        assert (g.reach == 0).all()

    def test_it_is_frozen(self) -> None:
        geometry_cls = module("radiation.geometry").RadiationGeometry

        g = geometry_cls.open_sky(facet_frame(), RadiationRequest())
        with pytest.raises((AttributeError, TypeError, ValueError)):
            g.n_directions = 32


class TestSynchronousCall:
    def test_it_is_the_binding_on_the_geometry(self, core: Any) -> None:
        evaluate_instant = module("radiation.evaluate").evaluate_instant
        geometry_cls = module("radiation.geometry").RadiationGeometry

        g = geometry_cls.open_sky(facet_frame(), RadiationRequest())
        got = evaluate_instant(g, APRIL_15, threads=1)
        want = core.evaluate_instant(
            g.lat, g.lon, g.normal, g.horizon, APRIL_15, limb_deg=0.0, threads=1
        )
        for name in ("lit", "cos_i", "b", "sin_e"):
            assert_array_equal(getattr(got, name), getattr(want, name))

    def test_the_upper_limb_comes_from_the_request(self, core: Any) -> None:
        evaluate_instant = module("radiation.evaluate").evaluate_instant
        geometry_cls = module("radiation.geometry").RadiationGeometry

        g = geometry_cls.open_sky(facet_frame(), RadiationRequest(limb="upper"))
        g = geometry_cls(
            lat=g.lat, lon=g.lon, normal=g.normal, horizon=np.zeros_like(g.horizon),
            reach=g.reach, n_directions=g.n_directions, request=g.request,
        )  # fmt: skip
        got = evaluate_instant(g, APRIL_15, pressure_hpa=900.0, temperature_c=-5.0)
        want = core.evaluate_instant(
            g.lat, g.lon, g.normal, g.horizon, APRIL_15,
            limb_deg=UPPER_LIMB, pressure_hpa=900.0, temperature_c=-5.0,
        )  # fmt: skip
        for name in ("lit", "cos_i", "b", "sin_e"):
            assert_array_equal(getattr(got, name), getattr(want, name))


class TestCombiner:
    """L2: dI = mean over samples of beam (b - max(0, sin e))."""

    def test_the_valley_floor_in_shadow_and_in_sun(self) -> None:
        beam_anomaly = module("snow.forcing").beam_anomaly

        s20, s40 = math.sin(math.radians(20.0)), math.sin(math.radians(40.0))
        shaded = beam_anomaly(np.array([[0.0]]), np.array([[s20]]), np.array([800.0]))
        assert shaded.shape == (1,)
        assert shaded[0] == pytest.approx(-273.6, abs=0.05)
        lit = beam_anomaly(np.array([[s40]]), np.array([[s40]]), np.array([800.0]))
        assert abs(lit[0]) <= 1e-9

    def test_it_is_the_mean_over_samples(self) -> None:
        beam_anomaly = module("snow.forcing").beam_anomaly

        b = np.array([[0.5, 0.0], [0.25, 0.1]])
        sin_e = np.array([[0.3, -0.2], [0.2, 0.1]])
        beam = np.array([800.0, 400.0])
        want = ((0.5 - 0.3) * 800 + (0.25 - 0.2) * 400) / 2, ((0.0 - 0.0) * 800 + 0.0) / 2
        assert beam_anomaly(b, sin_e, beam) == pytest.approx(want, abs=1e-12)

    def test_open_sky_flat_is_zero_at_every_sample(self, core: Any) -> None:
        # L2's last clause: --shadows off on a flat triangle gives dI = 0 within
        # 1e-9 W/m^2 (scale: beam up to 1400 W/m^2, b and sin_e equal bit for bit).
        beam_anomaly = module("snow.forcing").beam_anomaly

        got = core.evaluate_instant(**receivers([UP], OPEN), times_utc_ms=APRIL_15)
        assert_array_equal(got.b, np.maximum(got.sin_e, 0.0))
        for s in range(96):
            one = beam_anomaly(got.b[s : s + 1], got.sin_e[s : s + 1], np.array([1000.0]))
            assert abs(one[0]) <= 1e-9
