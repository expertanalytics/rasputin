"""The near raster and the horizon search over it: `radiation/search.py`.

`docs/increments/46-hbv-snow.md`, *10.4* (TR3b's near raster only: `open_dem`
over the mesh box grown by D_near, `raster.to_core` for the view, the search
reach D_far equal to D_near with no far raster) and *10.8* (L4; and, marked
slow, L3 and I8 on the synthetic 2 m mesh). Every test here calls
`_core.horizon_search` (TR3a), so each also fails while that kernel is absent.

Interface assumed (names chosen here and listed in the handback):

- `search.near_raster(sources, box, request, resources=None) -> DemTile`:
  `open_dem(DemRequest(sources=tuple(sources), bounds=Bounds.of(box grown by
  request.d_near_m on every side)), resources).tile`; `box` is the mesh's
  `(x_min, y_min, x_max, y_max)`.
- `search.search_horizons(frame, near, request, *, threads=0) ->
  RadiationGeometry`: `_core.horizon_search` over the frame's receivers on
  `raster.to_core(near)`, far raster `None`, `n_directions` from the request,
  curvature (1 - k) / (2 R) with R the frame's `earth_radius_m`,
  `d_near = d_far = request.d_near_m`.
- E2-style figures and I8 run only with `RASPUTIN_SLOW=1`.
"""

from __future__ import annotations

import datetime as dt
import importlib
import os
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from cli_helpers import ran
from tin_engine.io.geotiff_write import encode_geotiff
from tin_engine.io.models import Bounds
from tin_engine.radiation.request import RadiationRequest
from tin_engine.snow.terrain import TerrainSpec, mountain

SPEC = TerrainSpec()
#: A receiver 3.2 km east of the synthetic raster's west edge (x0 494000).
RX, RY = 497200.0, 6800000.0
WEST, EAST = 48, 16  # of 64 directions: 270 and 90 degrees true


def module(name: str) -> Any:
    """A module this increment adds, imported when a test runs (so a missing
    one fails that test, not the collection)."""
    return importlib.import_module(f"tin_engine.{name}")


@pytest.fixture(scope="module")
def synth(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("synth") / "synth.tif"
    path.write_bytes(encode_geotiff(mountain(SPEC), SPEC.x0, SPEC.y0, SPEC.step, SPEC.epsg))
    return path


@pytest.fixture(scope="module")
def search() -> Any:
    import tin_engine._core as core

    assert hasattr(core, "horizon_search"), "tin_engine._core has no horizon_search"
    return module("radiation.search")


def small_triangle() -> tuple[Any, Any]:
    """One 10 m triangle with its centroid at (RX, RY), 1000 m up."""
    pts = np.array([[RX - 5, RY - 5, 1000.0], [RX + 5, RY - 5, 1000.0], [RX, RY + 10, 1000.0]])
    return pts, np.array([[0, 1, 2]])


def frame_of(points: Any, triangles: Any) -> Any:
    from tin_engine.radiation.frame import receivers_from_mesh

    return receivers_from_mesh(points, triangles, f"EPSG:{SPEC.epsg}")


class TestNearRaster:
    """L4: the mesh box grown by 20 km, from the mesh path's own `open_dem`."""

    def test_it_asks_open_dem_for_the_box_grown_by_d_near(
        self, search: Any, synth: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        asked: list[Any] = []
        real = search.open_dem

        def spy(request: Any, *args: Any, **kwargs: Any) -> Any:
            asked.append(request)
            return real(request, *args, **kwargs)

        monkeypatch.setattr(search, "open_dem", spy)
        box = (RX - 5, RY - 5, RX + 5, RY + 10)
        tile = search.near_raster((synth,), box, RadiationRequest())
        assert len(asked) == 1
        assert asked[0].sources == (synth,)
        assert asked[0].bounds == Bounds.of((RX - 20005, RY - 20005, RX + 20005, RY + 20010))
        # Coverage ends at the file's edge: the whole grid, no more.
        assert tile.array.shape == (1201, 1201)

    def test_a_smaller_d_near_asks_for_less(
        self, search: Any, synth: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        asked: list[Any] = []
        real = search.open_dem

        def spy(request: Any, *args: Any, **kwargs: Any) -> Any:
            asked.append(request)
            return real(request, *args, **kwargs)

        monkeypatch.setattr(search, "open_dem", spy)
        search.near_raster((synth,), (RX, RY, RX + 1, RY + 1), RadiationRequest(d_near_m=1000.0))
        assert asked[0].bounds == Bounds.of((RX - 1000, RY - 1000, RX + 1001, RY + 1001))


class TestSearch:
    @pytest.fixture(scope="class")
    def searched(self, search: Any, synth: Path) -> tuple[Any, Any, Any]:
        pts, tri = small_triangle()
        frame = frame_of(pts, tri)
        request = RadiationRequest()
        near = search.near_raster((synth,), (RX - 5, RY - 5, RX + 5, RY + 10), request)
        return frame, near, search.search_horizons(frame, near, request, threads=1)

    def test_reach_says_where_coverage_ends(self, searched: tuple[Any, Any, Any]) -> None:
        _, _, g = searched
        assert g.reach.dtype == np.uint8 and g.reach.shape == (1, 64)
        assert g.reach[0, WEST] == 3  # the west edge 3.2 km away, floored to whole km
        assert g.reach[0, EAST] == 8  # the east edge (506000) 8.8 km away

    def test_the_geometry_carries_the_frame_and_the_request(
        self, searched: tuple[Any, Any, Any]
    ) -> None:
        frame, _, g = searched
        assert_array_equal(g.lat, frame.lat)
        assert_array_equal(g.lon, frame.lon)
        assert_array_equal(g.normal, frame.normal)
        assert g.n_directions == 64 and g.request == RadiationRequest()
        assert g.horizon.dtype == np.int16 and g.horizon.shape == (1, 64)

    def test_it_is_the_kernel_with_no_far_raster_and_the_frames_curvature(
        self, searched: tuple[Any, Any, Any]
    ) -> None:
        import tin_engine._core as core
        from tin_engine.raster import to_core

        frame, near, g = searched
        c = (1 - 0.13) / (2 * frame.earth_radius_m)
        horizon, reach = core.horizon_search(
            frame.x, frame.y, frame.z_obs, frame.lat, frame.lon, frame.jacobian, frame.normal,
            frame.exclusion_radius, near=to_core(near), far=None, n_directions=64,
            curvature=c, d_near=20_000.0, d_far=20_000.0, threads=1,
        )  # fmt: skip
        assert_array_equal(g.horizon, horizon)
        assert_array_equal(g.reach, reach)

    def test_the_massif_stands_above_the_receiver(self, searched: tuple[Any, Any, Any]) -> None:
        # The receiver at 1000 m, 2.8 km west of the 1900 m summit: the east
        # horizon is high, the probe able to fail on a search that saw nothing.
        _, _, g = searched
        assert g.horizon[0, EAST] > 1000


def _slow() -> Any:
    return pytest.mark.skipif(
        os.environ.get("RASPUTIN_SLOW") != "1", reason="slow (minutes): RASPUTIN_SLOW=1"
    )


@_slow()
class TestOnTheTwoMetreMesh:
    """L3 (shadows only take away) and I8 (1 and 8 threads, bit for bit) on the
    synthetic 2 m mesh (294 770 triangles), 96 samples of 15 min on 2027-04-15."""

    @pytest.fixture(scope="class")
    def mesh2(self, synth: Path) -> Any:
        from tin_engine.io.vtk_read import read_mesh

        target = synth.parent / "mesh2.vtk"
        ran("mesh", "--dem", str(synth), "--tolerance", "2", "--out", str(target))
        return read_mesh(target.read_bytes())

    @pytest.fixture(scope="class")
    def geometries(self, search: Any, synth: Path, mesh2: Any) -> dict[str, Any]:
        geometry_cls = module("radiation.geometry").RadiationGeometry
        frame = frame_of(mesh2.vertices, mesh2.triangles)
        assert len(frame.x) == 294_770
        v = mesh2.vertices
        box = (v[:, 0].min(), v[:, 1].min(), v[:, 0].max(), v[:, 1].max())
        request = RadiationRequest()
        near = search.near_raster((synth,), box, request)
        return {
            "one": search.search_horizons(frame, near, request, threads=1),
            "eight": search.search_horizons(frame, near, request, threads=8),
            "open": geometry_cls.open_sky(frame, request),
            "centre": ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2, float(v[:, 2].mean())),
        }

    @staticmethod
    def times() -> Any:
        start = int(dt.datetime(2027, 4, 15, tzinfo=dt.UTC).timestamp() * 1000)
        return start + 900_000 * np.arange(96, dtype=np.int64)

    def test_i8_the_horizons_are_the_same_at_1_and_8_threads(
        self, geometries: dict[str, Any]
    ) -> None:
        assert_array_equal(geometries["one"].horizon, geometries["eight"].horizon)
        assert_array_equal(geometries["one"].reach, geometries["eight"].reach)

    def test_i8_and_l3(self, geometries: dict[str, Any]) -> None:
        from tin_engine.crs import site_of

        evaluate_instant = module("radiation.evaluate").evaluate_instant
        forcing = module("snow.forcing")
        beam_anomaly, sun_and_beam = forcing.beam_anomaly, forcing.sun_and_beam

        cx, cy, h = geometries["centre"]
        site = site_of(f"EPSG:{SPEC.epsg}", cx, cy)
        times = self.times()
        _, beam = sun_and_beam(times, site.lat_deg, site.lon_deg, h)
        n = len(geometries["one"].lat)
        shaded, open_sky = np.zeros(n), np.zeros(n)
        # Chunks of 8 samples: three (8, n) evaluations at a time, not (96, n).
        for chunk in np.array_split(np.arange(96), 12):
            one = evaluate_instant(geometries["one"], times[chunk], threads=1)
            eight = evaluate_instant(geometries["one"], times[chunk], threads=8)
            sky = evaluate_instant(geometries["open"], times[chunk], threads=8)
            for name in ("lit", "cos_i", "b", "sin_e"):
                assert_array_equal(getattr(one, name), getattr(eight, name))  # I8
            assert (one.b <= sky.b).all()  # L3, every triangle and sample
            assert_array_equal(one.cos_i, sky.cos_i)
            shaded += beam_anomaly(one.b, one.sin_e, beam[chunk]) * len(chunk)
            open_sky += beam_anomaly(sky.b, sky.sin_e, beam[chunk]) * len(chunk)
        shaded, open_sky = shaded / 96, open_sky / 96
        assert (shaded <= open_sky).all()  # L3, the day's mean
        assert (shaded < open_sky).any()  # and some triangle is shaded at all
        # Probe S1, on 3000 DEM nodes: a mean change of -1.1 W/m^2; the
        # unweighted triangle mean here, the band E2's.
        assert -2.5 <= float(np.mean(shaded - open_sky)) <= -0.5
