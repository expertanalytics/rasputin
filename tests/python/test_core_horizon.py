"""The horizon kernel through the binding: `_core.horizon_search`.

`docs/increments/46-hbv-snow.md`, *10.8 Red suites for `@tester`* ("PR 5,
TR3a"), and the radiation design (`docs/increments/terrain-radiation.md`, R4,
R10, I4, verification 3): the V-valley's table and reach, a far raster of
`None` equal to one never reached (the absent far raster ends the search at
D_near), the same table at 1 and 8 threads, refusals in plain words, float32
rasters only, and the GIL released during the search. The C++ suite
(`tests/cpp/unit/test_radiation_horizon.cpp`) holds the mutation fixtures and
the pruned-equals-exhaustive checks; this file checks what crosses.

Interface assumed (names chosen here and listed in the handback):

- `_core.horizon_search(x, y, z_obs, lat_deg, lon_deg, jacobian, normal,
  exclusion_radius, near, far, *, n_directions, curvature, d_near, d_far,
  threads=0) -> (horizon, reach)`: the receiver arrays as
  `_core.check_receivers` takes them; `near` a `_core.RasterView` over a
  float32 array, `far` one or `None` (then `d_far` is ignored and the search
  ends at `d_near`); `horizon` an int16 array (n, n_directions) in 0.01
  degrees, direction k at true azimuth k 360 / n_directions; `reach` a uint8
  array of the same shape, whole km, 255 when nothing was missing. Refused
  with a `ValueError` naming the field or "receiver i"; a float64 view is a
  `TypeError` naming float32. Releases the GIL.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest

from gil_probe import ticks_during

C_EARTH = (1 - 0.13) / (2 * 6371000.0)  # k = 0.13, R = 6371 km
RELEASED_TICKS = 20  # as test_core_snow.py: a released call lets the ticker run

F64 = npt.NDArray[np.float64]


@pytest.fixture(scope="module")
def core() -> Any:
    import tin_engine._core as core

    assert hasattr(core, "horizon_search"), "tin_engine._core has no horizon_search"
    return core


def receivers(xy: F64, z_obs: F64, r_i: F64, gamma_deg: F64 | None = None) -> dict[str, F64]:
    """Receiver arrays as `check_receivers` takes them; J turned by gamma."""
    n = len(z_obs)
    g = np.radians(np.zeros(n) if gamma_deg is None else gamma_deg)
    jacobian = np.stack([np.cos(g), np.sin(g), -np.sin(g), np.cos(g)], axis=1).reshape(n, 2, 2)
    return {
        "x": np.ascontiguousarray(xy[:, 0]),
        "y": np.ascontiguousarray(xy[:, 1]),
        "z_obs": z_obs,
        "lat_deg": np.full(n, 61.3694),
        "lon_deg": np.full(n, 8.4852),
        "jacobian": jacobian,
        "normal": np.tile([0.0, 0.0, 1.0], (n, 1)),
        "exclusion_radius": r_i,
    }


def one(z_obs: float, r_i: float) -> dict[str, F64]:
    return receivers(np.zeros((1, 2)), np.array([z_obs]), np.array([r_i]))


def centred(core: Any, half: float, d: float, z: npt.NDArray[np.float32]) -> Any:
    """A view of a square float32 grid centred on (0, 0)."""
    return core.raster_view(z, x_min=-half, y_max=half, delta_x=d, delta_y=d)


def grid(half: float, d: float) -> tuple[F64, F64]:
    n = round(2 * half / d) + 1
    c = -half + d * np.arange(n)
    return np.meshgrid(c, -c)  # x by column, y falling down the rows


@pytest.fixture(scope="module")
def valley(core: Any) -> Any:
    """The V-valley: axis on x = 0, 30 degree slopes, ridges at 1000 m."""
    x, _ = grid(1500.0, 10.0)
    z = (np.minimum(np.abs(x), 1000.0) * math.tan(math.radians(30))).astype(np.float32)
    return centred(core, 1500.0, 10.0, z)


@pytest.fixture(scope="module")
def rugged(core: Any) -> tuple[Any, Any]:
    """Random near and far rasters with NaN holes (seeded, so every run is the same)."""
    rng = np.random.default_rng(46)
    x, y = grid(2500.0, 10.0)
    near = 300.0 * np.sin(x / 170.0) * np.cos(y / 230.0) + rng.uniform(-5, 5, x.shape)
    near.flat[rng.integers(0, near.size, 60)] = np.nan
    xf, yf = grid(5000.0, 100.0)
    far = 1200.0 * np.exp(-((xf - 3000.0) ** 2 + (yf + 1000.0) ** 2) / 1.5e6) + rng.uniform(
        0, 50, xf.shape
    )
    far.flat[rng.integers(0, far.size, 20)] = np.nan
    return (
        centred(core, 2500.0, 10.0, near.astype(np.float32)),
        centred(core, 5000.0, 100.0, far.astype(np.float32)),
    )


def search(core: Any, rc: dict[str, F64], near: Any, far: Any, **kw: Any) -> tuple[Any, Any]:
    p = {"n_directions": 64, "curvature": 0.0, "d_near": 1200.0, "d_far": 1200.0, **kw}
    return core.horizon_search(**rc, near=near, far=far, **p)  # type: ignore[no-any-return]


def test_valley_table(core: Any, valley: Any) -> None:
    horizon, reach = search(core, one(0.0, 5.0), valley, None)
    assert horizon.dtype == np.int16 and reach.dtype == np.uint8
    assert horizon.shape == reach.shape == (1, 64)
    # tan h = |sin A| tan 30: h(0) = h(180) = 0, h(45) = 22.2077, h(90) = 30.
    assert horizon[0, [0, 8, 16, 32, 48]].tolist() == [0, 2221, 3000, 0, 3000]
    assert (reach == 255).all()


def test_valley_turned_by_gamma(core: Any, valley: Any) -> None:
    """True azimuth = grid azimuth + gamma: the valley's axis moves to true 15 degrees."""
    rc = receivers(np.zeros((1, 2)), np.array([0.0]), np.array([5.0]), np.array([15.0]))
    horizon, _ = search(core, rc, valley, None)
    assert horizon[0, 8] == 1610  # atan(|sin 30| tan 30) = 16.1021 degrees


def test_far_none_is_a_far_raster_never_reached(core: Any, rugged: tuple[Any, Any]) -> None:
    near, far = rugged
    rng = np.random.default_rng(7)
    rc = receivers(rng.uniform(-800, 800, (20, 2)), rng.uniform(0, 300, 20), np.full(20, 10.0))
    none = search(core, rc, near, None, curvature=C_EARTH, d_near=1000.0, d_far=60000.0)
    unreached = search(core, rc, near, far, curvature=C_EARTH, d_near=1000.0, d_far=1000.0)
    reached = search(core, rc, near, far, curvature=C_EARTH, d_near=1000.0, d_far=4500.0)
    np.testing.assert_array_equal(none[0], unreached[0])
    np.testing.assert_array_equal(none[1], unreached[1])
    assert not np.array_equal(reached[0], none[0])  # the far peak shows once searched


def test_threads_give_the_same_table(core: Any, rugged: tuple[Any, Any]) -> None:
    near, far = rugged
    rng = np.random.default_rng(8)
    n = 400
    rc = receivers(
        rng.uniform(-1500, 1500, (n, 2)),
        rng.uniform(-50, 600, n),
        rng.uniform(0, 25, n),
        rng.uniform(-15, 15, n),
    )
    kw = {"curvature": C_EARTH, "d_near": 1000.0, "d_far": 4500.0}
    one_thread = search(core, rc, near, far, threads=1, **kw)
    eight = search(core, rc, near, far, threads=8, **kw)
    np.testing.assert_array_equal(one_thread[0], eight[0])
    np.testing.assert_array_equal(one_thread[1], eight[1])


def test_receiver_off_the_raster(core: Any, valley: Any) -> None:
    rc = receivers(np.array([[9000.0, 9000.0]]), np.array([0.0]), np.array([10.0]))
    horizon, reach = search(core, rc, valley, None)
    assert (horizon == -9000).all() and (reach == 0).all()


def test_no_receivers(core: Any, valley: Any) -> None:
    rc = receivers(np.zeros((0, 2)), np.zeros(0), np.zeros(0))
    horizon, reach = search(core, rc, valley, None)
    assert horizon.shape == reach.shape == (0, 64)


@pytest.mark.parametrize(
    ("change", "says"),
    [
        ({"n_directions": 0}, "n_directions"),
        ({"curvature": -1e-9}, "curvature"),
        ({"curvature": math.nan}, "curvature"),
        ({"d_near": 0.0}, "d_near"),
        ({"d_near": math.inf, "d_far": math.inf}, "d_near"),
        ({"far": "valley", "d_far": 1000.0}, "d_far"),
        ({"far": "valley", "d_far": math.nan}, "d_far"),
    ],
)
def test_refusals_name_the_field(core: Any, valley: Any, change: dict[str, Any], says: str) -> None:
    kw = {**change}
    far = valley if kw.pop("far", None) == "valley" else None
    with pytest.raises(ValueError, match=says):
        search(core, one(0.0, 5.0), valley, far, **kw)


def test_refusal_names_the_receiver(core: Any, valley: Any) -> None:
    rc = receivers(np.zeros((3, 2)), np.array([0.0, 0.0, math.nan]), np.full(3, 5.0))
    with pytest.raises(ValueError, match="receiver 2 has a number that is not finite"):
        search(core, rc, valley, None)


def test_float64_raster_is_refused(core: Any) -> None:
    x, _ = grid(500.0, 10.0)
    view64 = centred(core, 500.0, 10.0, np.zeros_like(x))  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="float32"):
        search(core, one(0.0, 5.0), view64, None)


def test_search_releases_the_gil(core: Any) -> None:
    # A bowl rising faster than linearly with distance, so every sample beats
    # the last and nothing is pruned: 3000 receivers x 64 rays x 1000 samples
    # of 5 m, about 0.6 s at the radiation design's 3 ns a sample, a tenfold
    # margin over the 0.05 s floor.
    x, y = grid(5100.0, 10.0)
    bowl = ((x * x + y * y) / 10000.0).astype(np.float32)
    near = centred(core, 5100.0, 10.0, bowl)
    rng = np.random.default_rng(9)
    n = 3000
    rc = receivers(rng.uniform(-50, 50, (n, 2)), np.zeros(n), np.full(n, 5.0))

    def run() -> Any:
        return search(core, rc, near, None, d_near=5000.0, d_far=5000.0, threads=1)

    _, ticks, elapsed = ticks_during(run)
    assert elapsed >= 0.05, f"took only {elapsed:.3f}s: too fast to see; enlarge the job"
    assert ticks >= RELEASED_TICKS, f"only {ticks} ticks in {elapsed:.3f}s: it holds the GIL"
