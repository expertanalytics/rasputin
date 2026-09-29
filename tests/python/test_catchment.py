"""`catchment.delineate`: seed, window loop and fine outline (increment 22, PR 1).

`docs/increments/22-auto-catchment.md`, "Data flow", "The seed", "The window",
"The fine outline" and "The red suites" (PR 1, `test_catchment.py`). The DEM
is an in-memory repository (`catchment_fixtures.MemoryRepository`, as
`mosaic_fixtures` builds tiles); no file is read.

Interface assumed (the design names `CatchmentRequest` (frozen Pydantic),
`delineate(request, repository) -> Catchment` (frozen), `CatchmentError`,
`WINDOW_MARGIN_M`; the fields are chosen here and stated in the handback):

- `CatchmentRequest(seed=(x, y), seed_crs="EPSG:4326", lakes=None,
  lakes_crs=None)`. `lakes` is a tuple of shapely polygons or multipolygons
  in `lakes_crs` (the CLI reads them from `--lakes`; paths stop there).
- `Catchment`: `fine` (a shapely `Polygon` in the DEM's CRS, no interiors),
  `crs` (text pyproj reads as the DEM's CRS), `nodes` (in-nodes),
  `seed_nodes`, `fine_area` (m^2), `rings_dropped`, `holes_filled`,
  `windows` (one entry per flood), and the last window's `mask` (in-nodes
  non-zero) with its `meta` (`RasterMeta`).
- `CatchmentError` is a `ValueError`.

THE WINDOW'S REFERENCE is one `_core.upstream` over the whole raster with the
same seed mask: the loop must end equal to it, node for node.

THE FLAGS: the C++ suite pins `touches_edge` and `touches_nodata` only where
the design's words and any repair agree (see its header). The refusals here
are pinned by behaviour: a catchment cut by the data's edge, or by NoData, is
refused, whatever mechanism detects it.
"""

from __future__ import annotations

import asyncio
from typing import Any

import numpy as np
import pytest
import shapely
from pydantic import ValidationError
from pyproj import Transformer
from shapely.geometry import MultiPolygon, Point, Polygon, box

import tin_engine.mosaic as mosaic
from catchment_fixtures import (
    EPSG,
    MARGIN_CELLS,
    X0,
    Y0,
    D,
    MemoryRepository,
    bowl,
    filled,
    full_flood,
    lake_box,
    lat,
    placed,
    seeds_in,
    tile_of,
)
from mosaic_fixtures import quadrants
from test_outline import square_area
from tin_engine.crs import parse_crs


@pytest.fixture(scope="module")
def api() -> Any:
    import tin_engine.catchment as catchment

    return catchment


def repository_of(z: np.ndarray, nodata: float | None = None) -> MemoryRepository:
    rows, cols = z.shape
    return MemoryRepository(
        quadrants(tile_of(z, nodata), row_cut=rows // 2, col_cut=cols // 2, overlap=1)
    )


def to_4326(x: float, y: float) -> tuple[float, float]:
    lon, lat_ = Transformer.from_crs(EPSG, "EPSG:4326", always_xy=True).transform(x, y)
    return float(lon), float(lat_)


def moved(geometry: Polygon, crs: str) -> Polygon:
    t = Transformer.from_crs(EPSG, crs, always_xy=True)
    return shapely.transform(geometry, lambda xy: np.column_stack(t.transform(xy[:, 0], xy[:, 1])))


CENTRE = lat(100, 200)  # the bowl's centre node, inside the lake


def request(api: Any, **kwargs: Any) -> Any:
    kwargs.setdefault("seed", CENTRE)
    kwargs.setdefault("seed_crs", EPSG)
    if "lakes" not in kwargs:
        kwargs["lakes"], kwargs["lakes_crs"] = (lake_box(),), EPSG
    return api.CatchmentRequest(**kwargs)


def on_whole(result: Any, shape: tuple[int, int]) -> np.ndarray:
    m = result.meta
    assert np.asarray(result.mask).shape == (m.rows, m.cols)
    assert m.delta_x == D and m.delta_y == D
    return placed(np.asarray(result.mask) != 0, m.x_min, m.y_max, shape)


# ---------------------------------------------------------------------------
# The request and the error
# ---------------------------------------------------------------------------


def test_the_margin_is_2000_metres(api: Any) -> None:
    assert api.WINDOW_MARGIN_M == 2000
    assert MARGIN_CELLS * D == api.WINDOW_MARGIN_M


def test_catchment_error_is_a_value_error(api: Any) -> None:
    assert issubclass(api.CatchmentError, ValueError)


def test_the_request_is_frozen_and_seed_crs_defaults_to_wgs84(api: Any) -> None:
    req = api.CatchmentRequest(seed=(8.5425, 61.3512))
    assert req.seed_crs == "EPSG:4326"
    assert req.lakes is None
    with pytest.raises(ValidationError):
        req.seed = (0.0, 0.0)


# ---------------------------------------------------------------------------
# The seed
# ---------------------------------------------------------------------------


def test_a_lake_polygon_seeds_the_nodes_inside_it(api: Any) -> None:
    z = bowl()
    result = api.delineate(request(api), repository_of(z))
    seed = seeds_in(tile_of(z), lake_box())
    assert seed.sum() == 49  # the fixture's own premise: 7 x 7 nodes
    assert result.seed_nodes == 49
    whole_mask = on_whole(result, z.shape)
    assert np.all(whole_mask[seed == 1] == 1)


def test_a_lake_from_geojson_seeds_the_same_nodes(api: Any) -> None:
    doc = {"type": "Polygon", "coordinates": [list(lake_box().exterior.coords)]}
    lake = shapely.geometry.shape(doc)
    z = bowl()
    a = api.delineate(request(api, lakes=(lake,), lakes_crs=EPSG), repository_of(z))
    b = api.delineate(request(api), repository_of(z))
    assert a.seed_nodes == b.seed_nodes == 49
    assert np.array_equal(on_whole(a, z.shape), on_whole(b, z.shape))


def test_a_multipolygon_contributes_the_part_containing_the_point(api: Any) -> None:
    z = bowl()
    far = box(*lat(20, 380), *lat(30, 370))  # 100 nodes, far from the bowl
    lakes = (MultiPolygon([lake_box(), far]),)
    result = api.delineate(request(api, lakes=lakes, lakes_crs=EPSG), repository_of(z))
    assert result.seed_nodes == 49


def test_the_point_in_no_lake_is_refused_naming_the_point(api: Any) -> None:
    x, y = lat(60, 60)
    with pytest.raises(api.CatchmentError, match=f"{int(x)}"):
        api.delineate(request(api, seed=(x, y)), repository_of(bowl()))


def test_the_point_in_two_lakes_is_refused(api: Any) -> None:
    lakes = (lake_box(), lake_box().buffer(D))
    with pytest.raises(api.CatchmentError, match=f"{int(CENTRE[0])}"):
        api.delineate(request(api, lakes=lakes, lakes_crs=EPSG), repository_of(bowl()))


def test_without_lakes_the_seed_is_the_nearest_node(api: Any) -> None:
    z = bowl()
    x, y = CENTRE
    result = api.delineate(
        api.CatchmentRequest(seed=(x + 30.0, y - 20.0), seed_crs=EPSG), repository_of(z)
    )
    assert result.seed_nodes == 1
    seed = np.zeros(z.shape, dtype=np.uint8)
    seed[200, 100] = 1
    assert np.array_equal(on_whole(result, z.shape), full_flood(tile_of(z), seed))


def test_a_seed_in_wgs84_lands_on_the_same_node_as_in_the_dem_crs(api: Any) -> None:
    z = bowl()
    x, y = CENTRE[0] + 30.0, CENTRE[1] - 20.0
    utm = api.delineate(api.CatchmentRequest(seed=(x, y), seed_crs=EPSG), repository_of(z))
    wgs = api.delineate(api.CatchmentRequest(seed=to_4326(x, y)), repository_of(z))
    assert wgs.seed_nodes == 1
    assert np.array_equal(on_whole(wgs, z.shape), on_whole(utm, z.shape))


def test_lakes_and_seed_in_other_crss_are_moved_into_the_dem_crs(api: Any) -> None:
    z = bowl()
    here = api.delineate(request(api), repository_of(z))
    there = api.delineate(
        api.CatchmentRequest(
            seed=to_4326(*CENTRE), lakes=(moved(lake_box(), "EPSG:3035"),), lakes_crs="EPSG:3035"
        ),
        repository_of(z),
    )
    assert there.seed_nodes == 49
    assert np.array_equal(on_whole(there, z.shape), on_whole(here, z.shape))


# ---------------------------------------------------------------------------
# The window
# ---------------------------------------------------------------------------


def test_a_catchment_larger_than_the_first_window_grows_and_equals_one_flood(api: Any) -> None:
    z = bowl()
    tile = tile_of(z)
    reference = full_flood(tile, seeds_in(tile, lake_box()))
    rows = np.nonzero(reference.any(axis=1))[0]
    # The fixture's premise: the lake's box plus the margin does not hold it.
    assert rows.min() < 200 - 3 - MARGIN_CELLS

    result = api.delineate(request(api), repository_of(z))
    assert len(result.windows) >= 2
    assert np.array_equal(on_whole(result, z.shape), reference)
    assert result.nodes == int(reference.sum())
    assert parse_crs(result.crs) == parse_crs(EPSG)


def test_the_fine_outline_is_the_filled_outer_ring_round_the_seed(api: Any) -> None:
    z = bowl()
    tile = tile_of(z)
    reference = full_flood(tile, seeds_in(tile, lake_box()))
    result = api.delineate(request(api), repository_of(z))
    fine = result.fine
    assert isinstance(fine, Polygon)
    assert fine.is_valid
    assert len(fine.interiors) == 0
    assert fine.exterior.is_ccw
    assert fine.contains(Point(CENTRE))
    assert fine.area == pytest.approx(square_area(filled(reference)) * D * D, rel=1e-12)
    assert result.fine_area == pytest.approx(fine.area, rel=1e-12)
    r, c = np.nonzero(reference)
    x, y = X0 + D * c, Y0 - D * r
    assert shapely.contains_xy(fine, x, y).all()


@pytest.mark.parametrize(
    ("side", "kwargs"), [("north", {"rc": 40}), ("east", {"cc": 190, "ax": 15.0})]
)
def test_a_catchment_cut_by_the_data_edge_is_refused_naming_the_side(
    api: Any, side: str, kwargs: dict[str, Any]
) -> None:
    rc, cc = kwargs.get("rc", 200), kwargs.get("cc", 100)
    z = bowl(**kwargs)
    with pytest.raises(api.CatchmentError, match=rf"(?i)\b{side}\b"):
        api.delineate(
            request(api, seed=lat(cc, rc), lakes=(lake_box(rc, cc),), lakes_crs=EPSG),
            repository_of(z),
        )


def test_a_catchment_cut_by_nodata_is_refused(api: Any) -> None:
    z = bowl()
    z[150:153, 99:102] = np.nan  # inside the basin, upstream of the lake
    with pytest.raises(api.CatchmentError, match=r"(?i)no ?data"):
        api.delineate(request(api), repository_of(z))


def test_a_catchment_cut_by_the_sentinel_is_refused(api: Any) -> None:
    z = bowl()
    z[150:153, 99:102] = -32767.0
    with pytest.raises(api.CatchmentError, match=r"(?i)no ?data"):
        api.delineate(request(api), repository_of(z, nodata=-32767.0))


def test_the_memory_cap_refuses_before_the_flood(api: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    # The first window is the lake's box (7 cells) plus 2 x 20 cells: about 48
    # x 48 nodes. With physical memory 10 x that, 15a's cap (4 bytes a node
    # against half of it) passes and the flood's (4 + 2 bytes a node) refuses.
    first = 48 * 48
    fake = 10 * first
    monkeypatch.setattr(mosaic, "physical_memory", lambda: fake)
    monkeypatch.setattr(api, "physical_memory", lambda: fake, raising=False)
    with pytest.raises(api.CatchmentError, match=r"(?i)memory"):
        api.delineate(request(api), repository_of(bowl()))


def test_the_loop_stops_early_far_from_the_data_edge(api: Any) -> None:
    """The design's pin after 0e1c29a: grow-or-stop is decided on the base
    margin, and only the step doubles. The loop that checked the doubled
    margin grew on every step until it met the data's edge."""
    z = bowl()
    tile = tile_of(z)
    reference = full_flood(tile, seeds_in(tile, lake_box()))
    r, c = np.nonzero(reference)
    rows, cols = z.shape
    # The fixture's premises: more than four base margins of data beyond the
    # catchment on every side, and a catchment larger than the first window.
    assert min(r.min(), c.min(), rows - 1 - r.max(), cols - 1 - c.max()) > 4 * MARGIN_CELLS
    assert r.min() < 200 - 3 - MARGIN_CELLS

    result = api.delineate(request(api), repository_of(z))
    assert 2 <= len(result.windows) <= 3
    m = result.meta
    assert m.x_min > X0
    assert m.y_max < Y0
    assert m.x_min + (m.cols - 1) * D < X0 + (cols - 1) * D
    assert m.y_max - (m.rows - 1) * D > Y0 - (rows - 1) * D
    assert np.array_equal(on_whole(result, z.shape), reference)


def test_a_catchment_inside_the_first_window_takes_one_window(api: Any) -> None:
    """A peak with the lake on top: nothing drains into the lake, so the
    catchment is the lake's 49 nodes, inside the first window (the lake's box
    plus the base margin) with the margin to spare."""
    r, c = np.indices((100, 100)).astype(np.float64)
    z = (100.0 - np.hypot(r - 50, c - 50)).astype(np.float32)
    lake = lake_box(50, 50)
    result = api.delineate(
        api.CatchmentRequest(seed=lat(50, 50), seed_crs=EPSG, lakes=(lake,), lakes_crs=EPSG),
        repository_of(z),
    )
    assert result.nodes == result.seed_nodes == 49
    assert len(result.windows) == 1


# ---------------------------------------------------------------------------
# Holes filled, other pieces dropped
# ---------------------------------------------------------------------------


def plus_with_a_hole() -> tuple[np.ndarray, Polygon]:
    """A 70 x 70 terrain falling from the centre (35, 35) to every edge, and a
    lake whose nodes are the plus round the centre (not the centre) and one
    node (35, 45) joined to it by a corridor holding no node.

    The plus is at 100, the centre at 200, the rest 50 - distance: the
    centre's diagonal neighbours (48.6) are flooded before the plus, so they
    reach the centre first and it is out, a hole in a catchment of five
    nodes in two pieces."""
    r, c = np.indices((70, 70)).astype(np.float64)
    z = (50.0 - np.hypot(r - 35, c - 35)).astype(np.float32)
    for rr, cc in ((34, 35), (36, 35), (35, 34), (35, 36), (35, 45)):
        z[rr, cc] = 100.0
    z[35, 35] = 200.0

    def diamond(radius: float) -> Polygon:
        return Polygon(
            [lat(35 + radius, 35), lat(35, 35 - radius), lat(35 - radius, 35), lat(35, 35 + radius)]
        )

    annulus = diamond(1.4).difference(diamond(0.5))
    corridor = box(*lat(36, 35.4), *lat(45, 35.3))
    square = box(*lat(44.6, 35.4), *lat(45.4, 34.6))
    lake = shapely.union_all([annulus, corridor, square])
    assert isinstance(lake, Polygon) and len(lake.interiors) == 1
    return z, lake


def test_holes_are_filled_and_other_pieces_dropped_and_counted(api: Any) -> None:
    z, lake = plus_with_a_hole()
    tile = tile_of(z)
    seed = seeds_in(tile, lake)
    assert {(int(a), int(b)) for a, b in zip(*np.nonzero(seed), strict=True)} == {
        (34, 35), (36, 35), (35, 34), (35, 36), (35, 45)
    }  # fmt: skip
    result = api.delineate(
        api.CatchmentRequest(seed=lat(35, 34), seed_crs=EPSG, lakes=(lake,), lakes_crs=EPSG),
        MemoryRepository({"t.tif": tile}),
    )
    assert result.nodes == 5
    assert result.seed_nodes == 5
    assert result.rings_dropped == 1
    assert result.holes_filled == 1
    fine = result.fine
    assert len(fine.interiors) == 0
    assert fine.contains(Point(lat(35, 35)))  # the hole, filled
    assert not fine.contains(Point(lat(45, 35)))  # the other piece, dropped
    assert fine.area == pytest.approx(4.5 * D * D, rel=1e-12)


# ---------------------------------------------------------------------------
# Async callers
# ---------------------------------------------------------------------------


async def test_delineate_runs_in_a_worker_thread(api: Any) -> None:
    z = bowl()
    here = api.delineate(request(api), repository_of(z))
    there = await asyncio.to_thread(api.delineate, request(api), repository_of(z))
    assert np.array_equal(np.asarray(there.mask), np.asarray(here.mask))
    assert there.fine.equals_exact(here.fine, 0.0)
