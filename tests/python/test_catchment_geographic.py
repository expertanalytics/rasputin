"""`catchment.delineate` on a DEM in degrees, flooding its own grid (increment 40).

`docs/increments/40-geographic-catchment.md`, sections 5.1 and 5.3 to 5.5 and
the red suites R3, R4, R7 and R8, with the window margin as a pair `(mx, my)`
through `_grow` and `_joined` (5.3). The DEMs are in memory
(`catchment_fixtures.MemoryRepository`); no file is read.

- R3: `catchment_fixtures.bowl()` on its projected 100 m grid and on 1" grids
  in EPSG:4326 at 0° and at 60° N: the same mask, node for node. The flood
  reads no coordinate.
- R4: on those 1" grids `nodes_area` is the sum over rows of the row's count
  times its cell's band area, by `geodesy_fixtures.band_area`'s closed form.
- R7: `gauge_fixtures.valley(dam=True)` and its twin in degrees
  (`geodesy_fixtures`, 10 m cells at the twin's centre latitude): the same
  chain, placed node and nodes; areas within 1e-3.
- R8: a reach in EPSG:4326 over the projected valley, and one in EPSG:25833
  over the twin: the same chain and nodes as the reach in the DEM's own CRS.

Interface from the design: `raster.to_core_lattice(tile)`;
`Catchment.nodes_area` and `Catchment.reduced_area` (m²); `_joined(m, out,
margin: tuple[float, float])`; `_grow(..., budget, base=(mx, my))` (base by
keyword); `gauge.place(..., cell_m=)`. Pinned here, where the design is
silent: `catchment_batch.seed_for(..., cell_m=)` by keyword.

Scales: degrees are compared to 1e-9 (about 0.1 mm) on values below 180;
metres of chain and offset to 1e-6 m on chains of 1.6 km; areas relative 1e-6
(R4's figure) or 1e-3 (R7's, the twin's 10 m cell holds only at its centre
latitude and UTM's scale is 0.9996).
"""

from __future__ import annotations

import math
import re
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest
from pyproj import Geod, Transformer
from shapely.geometry import box

import gauge_fixtures as gf
from catchment_fixtures import EPSG, MemoryRepository, bowl, lake_box, lat, tile_of
from geodesy_fixtures import (
    ARCSEC,
    TWIN_LAT0,
    TWIN_LON0,
    band_area,
    geo_tile,
    lattice_of,
    metres_per_degree,
    to_4326,
    twin_column_line,
    twin_tile,
    twin_xy,
)
from mosaic_fixtures import meta as projected_meta
from mosaic_fixtures import quadrants
from test_gauge import X, Y, lake, river_down
from tin_engine import _core
from tin_engine.catchment_batch import BatchRequest, seed_for
from tin_engine.gauge import Gauge
from tin_engine.io.models import Bounds
from tin_engine.raster import to_core
from tin_engine.resources import Resources

LON0 = 10.0


@pytest.fixture(scope="module")
def api() -> ModuleType:
    import tin_engine.catchment as catchment

    return catchment


@pytest.fixture(scope="module")
def gauge_api() -> ModuleType:
    import tin_engine.gauge as gauge

    return gauge


def split(tile: Any) -> MemoryRepository:
    rows, cols = tile.array.shape
    return MemoryRepository(quadrants(tile, row_cut=rows // 2, col_cut=cols // 2, overlap=1))


def on_whole(result: Any, whole: Any) -> np.ndarray:
    """The last window's mask placed on the whole grid `whole` (a meta)."""
    m = result.meta
    mask = np.asarray(result.mask) != 0
    r0 = round((whole.y_max - m.y_max) / whole.delta_y)
    c0 = round((m.x_min - whole.x_min) / whole.delta_x)
    out = np.zeros((whole.rows, whole.cols), dtype=np.uint8)
    out[r0 : r0 + mask.shape[0], c0 : c0 + mask.shape[1]] = mask
    return out


# ---------------------------------------------------------------- 5.1: a unit lattice


def test_to_core_lattice_takes_a_tile_in_degrees_and_floods_like_metres() -> None:
    """The hydrology reads the shape, the values and NoData only: the flood
    over the lattice view of a tile in degrees equals the flood over the
    projected view of the same array."""
    from tin_engine.raster import to_core_lattice

    z = bowl(rows=120, cols=60, rc=60, cc=30, ax=6.0, ay=20.0)
    geo = geo_tile(z, lon0=LON0, lat0=60.0, dx=ARCSEC, dy=ARCSEC)
    seed = np.zeros(z.shape, dtype=np.uint8)
    seed[60, 30] = 1
    here = _core.upstream(to_core_lattice(geo), seed)
    there = _core.upstream(to_core(tile_of(z)), seed)
    assert np.array_equal(np.asarray(here.mask), np.asarray(there.mask))
    assert here.nodes_in == there.nodes_in > 1


def test_to_core_lattice_needs_a_read_only_array() -> None:
    from tin_engine.io.models import DemTile
    from tin_engine.raster import to_core_lattice

    tile = tile_of(bowl(rows=20, cols=20, rc=10, cc=10, ax=3.0, ay=3.0))
    loose = DemTile.model_construct(meta=tile.meta, array=np.array(tile.array))
    assert loose.array.flags.writeable
    with pytest.raises(ValueError, match=r"read-only"):
        to_core_lattice(loose)


def test_to_core_still_refuses_degrees() -> None:
    """The mesh path still needs metres (5.1)."""
    geo = geo_tile(bowl(rows=20, cols=20, rc=10, cc=10), lon0=LON0, lat0=60.0, dx=ARCSEC, dy=ARCSEC)
    with pytest.raises(ValueError, match=r"(?i)geographic"):
        to_core(geo)


# ---------------------------------------------------------------- R3 and R4: the bowl in degrees


def geo_bowl(lat_c: float) -> Any:
    """The bowl on a 1" grid in EPSG:4326 whose centre row (200) is at `lat_c`."""
    return geo_tile(bowl(), lon0=LON0, lat0=lat_c + 200 * ARCSEC, dx=ARCSEC, dy=ARCSEC)


def geo_lake(tile: Any, rc: int = 200, cc: int = 100, half: float = 3.5) -> Any:
    m = tile.meta
    x0, y0 = m.node_xy(rc + half, cc - half)
    x1, y1 = m.node_xy(rc - half, cc + half)
    return box(x0, y0, x1, y1)


def bowl_in_degrees(api: ModuleType, lat_c: float) -> tuple[Any, Any]:
    tile = geo_bowl(lat_c)
    request = api.CatchmentRequest(
        seed=tile.meta.node_xy(200, 100),
        seed_crs="EPSG:4326",
        lakes=(geo_lake(tile),),
        lakes_crs="EPSG:4326",
    )
    return api.delineate(request, split(tile)), tile


@pytest.fixture(scope="module")
def projected_bowl(api: ModuleType) -> Any:
    request = api.CatchmentRequest(
        seed=lat(100, 200), seed_crs=EPSG, lakes=(lake_box(),), lakes_crs=EPSG
    )
    return api.delineate(request, split(tile_of(bowl())))


@pytest.mark.parametrize("lat_c", [0.0, 60.0], ids=["equator", "60N"])
def test_r3_the_same_array_on_two_grids_is_the_same_catchment(
    api: ModuleType, projected_bowl: Any, lat_c: float
) -> None:
    result, tile = bowl_in_degrees(api, lat_c)
    expected = on_whole(projected_bowl, tile_of(bowl()).meta)
    assert np.array_equal(on_whole(result, tile.meta), expected)
    assert result.nodes == projected_bowl.nodes == int(expected.sum())
    assert result.seed_nodes == projected_bowl.seed_nodes
    assert result.crs == "EPSG:4326"
    # The outline is traced on the native grid and comes out in degrees.
    x0, y0, x1, y1 = result.fine.bounds
    assert LON0 <= x0 < x1 <= LON0 + 200 * ARCSEC
    assert lat_c - 200 * ARCSEC <= y0 < y1 <= lat_c + 200 * ARCSEC


@pytest.mark.parametrize("lat_c", [0.0, 60.0], ids=["equator", "60N"])
def test_r4_nodes_area_is_the_row_weighted_geodesic_cell_area(
    api: ModuleType, lat_c: float
) -> None:
    """Each row's count times the band from half a spacing south of the row's
    latitude to half a spacing north, one spacing wide (WGS 84)."""
    result, tile = bowl_in_degrees(api, lat_c)
    whole = tile.meta
    counts = on_whole(result, whole).sum(axis=1)
    expected = 0.0
    for r in np.nonzero(counts)[0]:
        y = whole.y_max - r * whole.delta_y
        expected += counts[r] * band_area(y - ARCSEC / 2, y + ARCSEC / 2, ARCSEC)
    assert result.nodes_area == pytest.approx(expected, rel=1e-6)
    # A 1" cell is 950 m² at the equator and 480 m² at 60° N: not dx * dy.
    per_node = result.nodes_area / result.nodes
    assert per_node == pytest.approx(949.77 if lat_c == 0.0 else 479.69, rel=2e-3)


@pytest.mark.parametrize("lat_c", [0.0, 60.0], ids=["equator", "60N"])
def test_r4_outline_and_lake_areas_are_geodesic(api: ModuleType, lat_c: float) -> None:
    result, tile = bowl_in_degrees(api, lat_c)
    geod = Geod(ellps="WGS84")
    assert result.fine_area == pytest.approx(
        abs(geod.geometry_area_perimeter(result.fine)[0]), rel=1e-9
    )
    assert result.reduced_area == pytest.approx(
        abs(geod.geometry_area_perimeter(result.reduced)[0]), rel=1e-9
    )
    assert result.lake_area == pytest.approx(
        abs(geod.geometry_area_perimeter(geo_lake(tile))[0]), rel=1e-9
    )
    # The reduction keeps the area in its metre frame; on the ellipsoid the
    # frame's scale varies by 5e-4 over the outline, so to 1e-4 here.
    assert result.reduced_area == pytest.approx(result.fine_area, rel=1e-4)


def test_a_projected_catchment_keeps_its_areas_bit_for_bit(projected_bowl: Any) -> None:
    """Section 5.4: on a projected DEM each area is the float it is today."""
    r = projected_bowl
    assert r.nodes_area == r.nodes * r.meta.delta_x * r.meta.delta_y
    assert r.fine_area == r.fine.area
    assert r.reduced_area == r.reduced.area
    assert r.lake_area == lake_box().area


def test_the_default_tolerance_in_degrees_is_twice_the_longer_side_in_metres(
    api: ModuleType,
) -> None:
    """5.6: `2 * frame.cell_m(m)`, the frame at the fine outline's centre: at
    60° N the north-south side of a 1" cell, 30.95 m, so about 61.9 m."""
    result, _ = bowl_in_degrees(api, 60.0)
    assert result.tolerance == pytest.approx(2 * ARCSEC * metres_per_degree(60.0)[1], rel=1e-4)


# ---------------------------------------------------------------- 5.3: the margin pair


def test_joined_grows_x_by_mx_and_y_by_my(api: ModuleType) -> None:
    m = projected_meta(rows=10, cols=10, dx=10.0, dy=10.0)
    out = SimpleNamespace(row_min=2, row_max=5, col_min=3, col_max=7)
    x0, y0 = m.x_min, m.y_max
    got = api._joined(m, out, (100.0, 200.0))
    assert got == Bounds(x_min=x0 - 70.0, y_min=y0 - 250.0, x_max=x0 + 170.0, y_max=y0 + 180.0)


def spy(monkeypatch: pytest.MonkeyPatch, api: ModuleType, name: str, arg: int) -> list[Any]:
    """Record argument `arg` of every call of `api.<name>`, calling through."""
    seen: list[Any] = []
    real = getattr(api, name)

    def wrapper(*args: Any, **kwargs: Any) -> Any:
        seen.append(args[arg])
        return real(*args, **kwargs)

    monkeypatch.setattr(api, name, wrapper)
    return seen


def test_grow_starts_from_base_and_doubles_both_sides(
    api: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """From a window of a few cells round the seed, the bowl's catchment
    grows the window more than once: each growth test is `_joined(m, out,
    base)`, each next window `_joined` with the pair doubled, (2^k mx, 2^k my)."""
    tile = geo_bowl(60.0)
    m = tile.meta
    point = m.node_xy(200.0, 100.0)
    base = (2 * m.delta_x, 3 * m.delta_y)
    bounds = Bounds(
        x_min=point[0] - base[0], y_min=point[1] - base[1],
        x_max=point[0] + base[0], y_max=point[1] + base[1],
    )  # fmt: skip

    def flood(window: Any) -> Any:
        seed = api._seed_mask(window.meta, None, point)
        return api.upstream(window, seed), seed

    margins = spy(monkeypatch, api, "_joined", 2)
    budget = Resources.fallback().memory_budget
    repo = split(tile)
    api._grow(repo.footprints(), repo, bounds, flood, burnt=False, budget=budget, base=base)
    assert len(margins) >= 5 and len(margins) % 2 == 1, margins
    assert all(tuple(g) == base for g in margins[0::2]), margins
    doubled = [(2.0**k * base[0], 2.0**k * base[1]) for k in range(1, len(margins) // 2 + 1)]
    assert [tuple(g) for g in margins[1::2]] == doubled


def test_delineate_in_degrees_takes_the_margin_at_the_seed(
    api: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The first window is the seed's box grown by `units(2000)` at the seed;
    every growth test uses that pair."""
    tile = geo_bowl(60.0)
    x, y = tile.meta.node_xy(200.2, 100.3)
    kx, ky = metres_per_degree(y)
    mx, my = 2000.0 / kx, 2000.0 / ky
    plans = spy(monkeypatch, api, "_plan", 1)
    margins = spy(monkeypatch, api, "_joined", 2)
    api.delineate(api.CatchmentRequest(seed=(x, y), seed_crs="EPSG:4326"), split(tile))
    first = plans[0]
    assert (first.x_min, first.y_min, first.x_max, first.y_max) == pytest.approx(
        (x - mx, y - my, x + mx, y + my), rel=0.0, abs=1e-9
    )
    assert tuple(margins[0]) == pytest.approx((mx, my), rel=1e-12)
    assert mx == pytest.approx(2 * my, rel=1e-2)  # the pair is a pair, not one float


def test_delineate_on_a_projected_dem_keeps_todays_floats(
    api: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """R11 at the unit level: 2000.0, then 4000.0, 8000.0, ... on both axes."""
    x, y = lat(100.3, 200.2)
    plans = spy(monkeypatch, api, "_plan", 1)
    margins = spy(monkeypatch, api, "_joined", 2)
    api.delineate(api.CatchmentRequest(seed=(x, y), seed_crs=EPSG), split(tile_of(bowl())))
    first = plans[0]
    assert (first.x_min, first.y_min, first.x_max, first.y_max) == (
        x - 2000.0, y - 2000.0, x + 2000.0, y + 2000.0,
    )  # fmt: skip
    assert len(margins) >= 3, margins  # the bowl grows the window at least once
    assert all(tuple(g) == (2000.0, 2000.0) for g in margins[0::2]), margins
    assert [tuple(g) for g in margins[1::2]] == [
        (2000.0 * 2**k, 2000.0 * 2**k) for k in range(1, len(margins) // 2 + 1)
    ]


# ---------------------------------------------------------------- R7 and R8: the gauge in degrees

AT_ROW, ROW0, ROW1 = 150, 40, 200
CHAIN = [(r, gf.CC) for r in range(ROW0, ROW1 + 1)]


def proj_reach(gauge_api: ModuleType) -> Any:
    return gauge_api.Reach(
        line=gf.column_line(gf.CC, ROW0, ROW1),
        at=(AT_ROW - ROW0) * gf.D,
        uncertainty=30.0,
        corridor=5.0,
    )


def twin_reach(gauge_api: ModuleType) -> Any:
    """The same reach on the twin, the line moved by lattice position; `at`
    in metres, as the twin's cells are 10 m at the line's centre."""
    return gauge_api.Reach(
        line=twin_column_line(gf.CC, ROW0, ROW1),
        at=(AT_ROW - ROW0) * gf.D,
        uncertainty=30.0,
        corridor=5.0,
    )


def frame_length(a: tuple[float, float], b: tuple[float, float]) -> float:
    """Metres from `a` to `b` (degrees) in an affine frame at their midpoint."""
    kx, ky = metres_per_degree((a[1] + b[1]) / 2)
    return math.hypot((b[0] - a[0]) * kx, (b[1] - a[1]) * ky)


@pytest.fixture(scope="module")
def projected_gauge(api: ModuleType, gauge_api: ModuleType) -> Any:
    request = api.CatchmentRequest(
        seed=gf.lat(gf.CC + 2, AT_ROW), seed_crs=gf.EPSG, reach=proj_reach(gauge_api)
    )
    return api.delineate(request, split(gf.tile_of(gf.valley(dam=True))))


@pytest.fixture(scope="module")
def twin_gauge(api: ModuleType, gauge_api: ModuleType) -> Any:
    request = api.CatchmentRequest(
        seed=twin_xy(gf.CC + 2, AT_ROW), seed_crs="EPSG:4326", reach=twin_reach(gauge_api)
    )
    return api.delineate(request, split(twin_tile(gf.valley(dam=True))))


PROJ_WHOLE = gf.tile_of(gf.valley(dam=True)).meta
TWIN_WHOLE = twin_tile(gf.valley(dam=True)).meta


def assert_same_gauge(result: Any, whole: Any, expected: Any) -> None:
    """The chain, the placed node and the catchment of `expected` (the
    projected valley's run), node for node on the whole lattice."""
    g, e = result.gauge, expected.gauge
    assert g is not None and e is not None
    assert lattice_of(whole, g.chain) == lattice_of(PROJ_WHOLE, e.chain) == CHAIN
    assert lattice_of(whole, [g.node]) == lattice_of(PROJ_WHOLE, [e.node]) == [(AT_ROW, gf.CC)]
    assert result.nodes == expected.nodes
    assert np.array_equal(on_whole(result, whole), on_whole(expected, PROJ_WHOLE))
    assert g.lowered_nodes == e.lowered_nodes == len(gf.DAM_ROWS)


def test_r7_a_gauge_on_the_geographic_twin(twin_gauge: Any, projected_gauge: Any) -> None:
    assert_same_gauge(twin_gauge, TWIN_WHOLE, projected_gauge)
    g, e = twin_gauge.gauge, projected_gauge.gauge
    assert twin_gauge.crs == "EPSG:4326"
    assert g.node_offset_m == pytest.approx(0.0, abs=1e-6)
    assert g.downstream_checked == e.downstream_checked == "whole"
    assert twin_gauge.nodes_area == pytest.approx(projected_gauge.nodes_area, rel=1e-3)
    assert twin_gauge.fine_area == pytest.approx(projected_gauge.fine_area, rel=1e-3)
    assert twin_gauge.reduced_area == pytest.approx(projected_gauge.reduced_area, rel=1e-3)
    assert g.sensitivity.a0 == pytest.approx(e.sensitivity.a0, rel=1e-3)
    # The chain's own nodes are the reference: a0 is the placed node's area.
    assert g.sensitivity.a0 == pytest.approx(twin_gauge.nodes_area / 1e6, rel=1e-9)


def test_r7_the_gauged_window_pads_by_corridor_plus_margin_in_degrees(
    api: ModuleType, gauge_api: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """5.3: the pad is `frame.units(corridor + 2000)` and the base passed to
    `_grow` is `frame.units(2000)`, the burn frame at the centre of the
    line's box (5.5)."""
    reach = twin_reach(gauge_api)
    (x0, y0), (x1, y1) = reach.line[0], reach.line[-1]
    lo_x, lo_y, hi_x, hi_y = min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)
    kx, ky = metres_per_degree((lo_y + hi_y) / 2)
    pad = (5.0 + 2000.0) / kx, (5.0 + 2000.0) / ky
    plans = spy(monkeypatch, api, "_plan", 1)
    margins = spy(monkeypatch, api, "_joined", 2)
    request = api.CatchmentRequest(
        seed=twin_xy(gf.CC + 2, AT_ROW), seed_crs="EPSG:4326", reach=reach
    )
    api.delineate(request, split(twin_tile(gf.valley(dam=True))))
    first = plans[0]
    assert (first.x_min, first.y_min, first.x_max, first.y_max) == pytest.approx(
        (lo_x - pad[0], lo_y - pad[1], hi_x + pad[0], hi_y + pad[1]), rel=0.0, abs=1e-9
    )
    base = (2000.0 / kx, 2000.0 / ky)
    for g in margins:
        k = round(math.log2(g[0] / base[0]))
        assert tuple(g) == pytest.approx((2.0**k * base[0], 2.0**k * base[1]), rel=1e-9), margins


def test_r8_a_reach_in_degrees_over_the_projected_valley(
    api: ModuleType, gauge_api: ModuleType, projected_gauge: Any
) -> None:
    """The reach's line and `at` in EPSG:4326 (`at` in metres in its frame):
    moved into the DEM's CRS and `at` measured again there (5.5)."""
    a, b = to_4326(*gf.lat(gf.CC, ROW0)), to_4326(*gf.lat(gf.CC, ROW1))
    p = to_4326(*gf.lat(gf.CC, AT_ROW))
    reach = gauge_api.Reach(line=(a, b), at=frame_length(a, p), uncertainty=30.0, corridor=5.0)
    request = api.CatchmentRequest(
        seed=to_4326(*gf.lat(gf.CC + 2, AT_ROW)), seed_crs="EPSG:4326", reach=reach
    )
    result = api.delineate(request, split(gf.tile_of(gf.valley(dam=True))))
    assert result.crs == gf.EPSG
    assert_same_gauge(result, PROJ_WHOLE, projected_gauge)
    assert result.gauge.node_offset_m == pytest.approx(0.0, abs=0.05)


def test_r8_a_reach_in_metres_over_the_twin(
    api: ModuleType, gauge_api: ModuleType, projected_gauge: Any, twin_gauge: Any
) -> None:
    """The projected valley's reach, in EPSG:25833, over the twin in degrees."""
    request = api.CatchmentRequest(
        seed=gf.lat(gf.CC + 2, AT_ROW), seed_crs=gf.EPSG, reach=proj_reach(gauge_api)
    )
    result = api.delineate(request, split(twin_tile(gf.valley(dam=True))))
    assert result.crs == "EPSG:4326"
    assert_same_gauge(result, TWIN_WHOLE, projected_gauge)
    assert result.nodes_area == pytest.approx(twin_gauge.nodes_area, rel=1e-12)
    assert result.gauge.node_offset_m == pytest.approx(0.0, abs=1.0)


def test_a_reach_vertex_with_no_image_in_the_dem_crs_is_refused(
    api: ModuleType, gauge_api: ModuleType
) -> None:
    """5.7: in plain words, as a seed with no image is today. The reach is
    in an orthographic CRS; its second vertex lies off the projection's disc
    (radius about 6.4e6 m), so it has no longitude and latitude."""
    ortho = "+proj=ortho +lat_0=-60 +lon_0=-165 +ellps=WGS84 +units=m"
    reach = gauge_api.Reach(
        line=((0.0, 0.0), (1.0e7, 1.0e7)), at=10.0, uncertainty=30.0, corridor=5.0
    )
    request = api.CatchmentRequest(seed=(0.0, 0.0), seed_crs=ortho, reach=reach)
    with pytest.raises(api.CatchmentError, match=r"(?i)no image"):
        api.delineate(request, split(twin_tile(gf.valley(dam=True))))


# ---------------------------------------------------------------- seed_for: cells and frames


def segments_in(crs: str) -> list[Any]:
    """`test_gauge.river_down()`, a straight river south at x = X, in `crs`."""
    t = Transformer.from_crs(gf.EPSG, crs, always_xy=True)
    return [
        s.model_copy(
            update={"line": tuple(tuple(map(float, t.transform(x, y))) for x, y in s.line)}
        )
        for s in river_down()
    ]


@pytest.mark.parametrize("cell", [10.0, 29.8324])
def test_seed_for_passes_the_cell_on(cell: float) -> None:
    """`U`'s floor and the burn's half-width are three cells (5.3, 5.5)."""
    placement, _, made = seed_for(
        Gauge(x=X + 5.0, y=Y), river_down(), gf.EPSG, BatchRequest(), cell_m=cell
    )
    assert placement.reach.uncertainty == pytest.approx(3 * cell, rel=1e-12)
    assert made.reach.corridor == pytest.approx(3 * cell, rel=1e-12)


@pytest.mark.parametrize(("cell", "seeded"), [(10.0, True), (9.9, False)])
def test_seed_for_passes_the_lake_gap_as_three_cells(cell: float, seeded: bool) -> None:
    """The station 30 m from a lake its `P` lies in: seeded at 3 x 10 m, not
    at 3 x 9.9 m (the lake gap is `U`'s floor, 5.3)."""
    vatn = lake(X - 500.0, Y - 500.0, X + 5.0, Y + 500.0)
    _, seed, _ = seed_for(
        Gauge(x=X + 35.0, y=Y), river_down(kind="lake"), gf.EPSG, BatchRequest(), [vatn],
        cell_m=cell,
    )  # fmt: skip
    assert (seed is not None) is seeded


def test_seed_for_places_a_river_file_in_degrees_in_metres() -> None:
    """The river file in EPSG:4326: placed in a frame at the station, the
    reach and `P` moved back to degrees, distances in metres (5.5)."""
    lon, lat_ = to_4326(X + 20.0, Y)
    metres, _, _ = seed_for(
        Gauge(x=X + 20.0, y=Y), river_down(), gf.EPSG, BatchRequest(), cell_m=10.0
    )
    placed, _, made = seed_for(
        Gauge(x=lon, y=lat_), segments_in("EPSG:4326"), "EPSG:4326", BatchRequest(), cell_m=10.0
    )
    # UTM's scale is 0.9996 here: metres on the ground are 4e-4 longer.
    assert placed.distance_m == pytest.approx(metres.distance_m, rel=1e-3)
    assert placed.reach.uncertainty == 30.0
    assert placed.reach_up_m == pytest.approx(metres.reach_up_m, rel=1e-3)
    assert placed.position == pytest.approx(to_4326(*metres.position), abs=1e-7)
    xs, ys = zip(*made.reach.line, strict=True)
    assert all(abs(v) <= 180.0 for v in xs) and all(abs(v) <= 90.0 for v in ys)
    assert made.seed_crs == "EPSG:4326"


def test_seed_for_with_lakes_and_a_river_file_in_degrees_is_refused() -> None:
    """5.5: only the batch passes lakes, and it needs a DEM in metres; the
    refusal names the reason."""
    geo = segments_in("EPSG:4326")
    lon, lat_ = geo[2].line[0]
    vatn = lake(lon - 0.01, lat_ - 0.01, lon + 0.01, lat_ + 0.01)
    reason = re.escape("lakes need a river file in metres (a projected CRS)")
    with pytest.raises(ValueError, match=reason):
        seed_for(Gauge(x=lon, y=lat_), geo, "EPSG:4326", BatchRequest(), [vatn], cell_m=10.0)


# ---------------------------------------------------------------- the twin's own premises


def test_the_twins_cells_are_10_m_at_its_centre() -> None:
    """The fixture's premise, by the oracle: 10 m both ways at row 120."""
    lon, lat_ = twin_xy(gf.CC, gf.RIM_ROW)
    kx, ky = metres_per_degree(lat_)
    m = TWIN_WHOLE
    assert m.delta_x * kx == pytest.approx(10.0, rel=1e-12)
    assert m.delta_y * ky == pytest.approx(10.0, rel=1e-12)
    assert (m.x_min, m.y_max) == (TWIN_LON0, TWIN_LAT0)
    assert to_4326(*gf.lat(gf.CC, gf.RIM_ROW)) == pytest.approx((lon, lat_), abs=1e-12)
