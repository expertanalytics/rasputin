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
import hashlib
from typing import Any

import numpy as np
import pytest
import shapely
from pydantic import ValidationError
from pyproj import Transformer
from shapely.geometry import MultiPolygon, Point, Polygon, box

import gauge_fixtures as gf
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
from crs_fixtures import proj4_of
from mosaic_fixtures import quadrants
from test_outline import square_area
from tin_engine._core import accumulate, upstream
from tin_engine.crs import parse_crs
from tin_engine.raster import to_core as core_view


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


# ---------------------------------------------------------------------------
# Increment 29, PR 2: the gauge on the river (a `Reach` in the request)
#
# `docs/increments/29-nve-reference-catchments.md`, "The window loop and the
# catchment" and "The red suites", PR 2's `test_catchment.py`. The terrains are
# `gauge_fixtures`' 10 m ones; the reach is given directly, with
# `corridor=5.0` and the line on nodes, so the burnt chain is the line's own
# nodes and the hand burn (`gauge_fixtures.hand_burn`) is the reference.
#
# Pinned here (names the design leaves open are listed in the handback):
# `CatchmentRequest.reach` (default None) and `Catchment.gauge` (None without
# a reach), a `GaugeResult` with the design's fields; `gauge.sensitivity` is
# `sensitivity.Sensitivity` (`a0` in km2, `causes`, `well_posed`).
# ---------------------------------------------------------------------------


#: 22's results on the bowl, taken on the tree before PR 2 (529613a) over
#: the lattice-exact fields only (the mask, the fine ring, the counts, the
#: windows and the meta; not the reduced ring, whose C++ arithmetic may round
#: differently on another compiler): `reach=None` must leave them unchanged.
BOWL_LAKE_DIGEST = "0efdc506ba2c0d63409f13da1bd761388b34762d2cc10bd26e189755cd902c5f"
BOWL_POINT_DIGEST = "02ca6c38e46f445b5bf01996531616c130a0942766478a26d69b24346ac1740c"


def digest(result: Any) -> str:
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(np.asarray(result.mask, dtype=np.uint8)).tobytes())
    coords = np.asarray(result.fine.exterior.coords, dtype=np.float64)
    h.update(np.ascontiguousarray(coords).tobytes())
    fields = (
        result.nodes, result.seed_nodes, result.seed, result.rings_dropped,
        result.dropped_nodes, result.holes_filled, result.lake_area, result.tolerance,
        result.crs, [(w.bounds, w.rows, w.cols, w.grown) for w in result.windows], result.meta,
    )  # fmt: skip
    h.update(repr(fields).encode())
    return h.hexdigest()


@pytest.fixture(scope="module")
def gauge_api() -> Any:
    import tin_engine.gauge as gauge

    return gauge


def valley_repository(z: np.ndarray, nodata: float | None = None) -> MemoryRepository:
    rows, cols = z.shape
    tile = gf.tile_of(z, nodata)
    return MemoryRepository(quadrants(tile, row_cut=rows // 2, col_cut=cols // 2, overlap=1))


def on_lattice(result: Any, shape: tuple[int, int]) -> np.ndarray:
    """The last window's mask on the whole 10 m raster's lattice."""
    m = result.meta
    mask = np.asarray(result.mask) != 0
    r0 = round((gf.Y0 - m.y_max) / gf.D)
    c0 = round((m.x_min - gf.X0) / gf.D)
    out = np.zeros(shape, dtype=np.uint8)
    out[r0 : r0 + mask.shape[0], c0 : c0 + mask.shape[1]] = mask
    return out


def whole_flood(z: np.ndarray, node: tuple[int, int], nodata: float | None = None) -> Any:
    seed = np.zeros(z.shape, dtype=np.uint8)
    seed[node] = 1
    return upstream(core_view(gf.tile_of(z, nodata)), seed)


def gauge_request(
    api: Any,
    gauge_api: Any,
    *,
    col: int = gf.CC,
    row0: int = 40,
    row1: int = 200,
    at_row: int = 150,
) -> Any:
    reach = gauge_api.Reach(
        line=gf.column_line(col, row0, row1),
        at=(at_row - row0) * gf.D,
        uncertainty=30.0,
        corridor=5.0,
    )
    return api.CatchmentRequest(seed=gf.lat(col + 2, at_row), seed_crs=gf.EPSG, reach=reach)


def test_the_reach_defaults_to_none(api: Any) -> None:
    assert api.CatchmentRequest(seed=(8.5, 61.3)).reach is None


def test_a_reach_with_lakes_is_refused(api: Any, gauge_api: Any) -> None:
    reach = gauge_api.Reach(line=gf.column_line(gf.CC, 40, 200), at=1100.0, uncertainty=30.0)
    with pytest.raises(ValueError, match=r"(?i)lake"):
        api.CatchmentRequest(
            seed=CENTRE, seed_crs=EPSG, lakes=(lake_box(),), lakes_crs=EPSG, reach=reach
        )


def reach_refusal(api: Any, gauge_api: Any, seed_crs: str) -> str:
    """`delineate`'s refusal of a gauge request whose reach is in `seed_crs`,
    over the valley in EPSG:25833 (`gf.EPSG`)."""
    request = gauge_request(api, gauge_api).model_copy(update={"seed_crs": seed_crs})
    with pytest.raises(api.CatchmentError) as info:
        api.delineate(request, valley_repository(gf.valley(dam=True)))
    return str(info.value)


def test_a_reach_in_another_crs_is_refused(api: Any, gauge_api: Any) -> None:
    """A reach is burnt as given, so it must be in the DEM's CRS: EPSG:32633
    (UTM 33 on WGS 84, not ETRS89) is refused naming the DEM's."""
    message = reach_refusal(api, gauge_api, "EPSG:32633")
    assert message.startswith(f"the river reach must be in the DEM's CRS, {gf.EPSG}"), message


def test_a_reach_in_a_proj_string_naming_no_datum_is_refused_with_a_hint(
    api: Any, gauge_api: Any
) -> None:
    """Audit PR B, code review round 1, and Ola's D15 b: pyproj's PROJ string
    of EPSG:25833 names the GRS80 ellipsoid, not ETRS89, so it is not the
    DEM's CRS; the refusal names the code to write instead."""
    message = reach_refusal(api, gauge_api, proj4_of(25833))
    assert message.startswith(f"the river reach must be in the DEM's CRS, {gf.EPSG}"), message
    assert f"write {gf.EPSG}" in message, message


def test_no_reach_gives_22s_result_bit_for_bit(api: Any) -> None:
    z = bowl()
    lake = api.delineate(request(api, reach=None), repository_of(z))
    assert lake.gauge is None
    assert digest(lake) == BOWL_LAKE_DIGEST
    point = api.CatchmentRequest(seed=lat(100.3, 200.2), seed_crs=EPSG, reach=None)
    result = api.delineate(point, repository_of(z))
    assert result.gauge is None
    assert digest(result) == BOWL_POINT_DIGEST


def test_a_gauge_beside_the_river_drains_the_whole_valley_past_the_embankment(
    api: Any, gauge_api: Any
) -> None:
    """The catchment equals one flood, over the whole raster, from the
    hand-burnt placed node, and the valley above the embankment is in it."""
    z = gf.valley(dam=True)
    chain = [(r, gf.CC) for r in range(40, 201)]
    burnt = gf.hand_burn(z, chain)
    placed = (150, gf.CC)
    expected = whole_flood(burnt, placed)
    assert np.asarray(whole_flood(z, placed).mask)[50, gf.CC] == 0  # the premise
    result = api.delineate(gauge_request(api, gauge_api), valley_repository(z))
    mask = on_lattice(result, z.shape)
    assert np.array_equal(mask, np.asarray(expected.mask))
    assert mask[50, gf.CC] == 1
    assert result.nodes == expected.nodes_in
    g = result.gauge
    assert g is not None
    assert g.node == pytest.approx(gf.lat(gf.CC, 150))
    assert len(g.chain) == len(chain)
    assert np.allclose(np.asarray(g.chain), [gf.lat(c, r) for r, c in chain])
    assert g.chain_nodes == len(chain)
    assert g.node_offset_m == pytest.approx(0.0, abs=1e-6)
    assert g.lowered_nodes == len(gf.DAM_ROWS)
    assert g.direction_ok and g.end_closed
    assert g.end_extended_m == 0.0
    assert g.downstream_checked == "whole"
    count = np.asarray(accumulate(core_view(gf.tile_of(burnt))).count)
    assert g.sensitivity.a0 == pytest.approx(count[placed] * gf.CELL_KM2, rel=1e-12)
    assert count[placed] == result.nodes  # the exact oracle, through the real path
    assert g.sensitivity.well_posed
    assert tuple(g.sensitivity.causes) == ()


#: A NoData patch east of the floor, rows 153 to 157: a tributary from it joins
#: the river 20 m below the placed node (row 150), inside U, and never above it.
TRIBUTARY_NODATA = (slice(153, 158), slice(gf.CC + 3, gf.CC + 7))
#: The same patch in the valley above the placed node.
UPPER_NODATA = (slice(30, 34), slice(gf.CC + 3, gf.CC + 7))
NODATA = -32767.0


def test_nodata_reached_only_below_the_placed_node_is_not_a_refusal(
    api: Any, gauge_api: Any
) -> None:
    z = gf.valley(dam=True)
    z[TRIBUTARY_NODATA] = NODATA
    chain = [(r, gf.CC) for r in range(40, 201)]
    burnt = gf.hand_burn(z, chain)
    expected = whole_flood(burnt, (150, gf.CC), NODATA)
    assert not expected.touches_nodata  # the premises: the placed catchment is clear,
    assert whole_flood(burnt, (153, gf.CC), NODATA).touches_nodata  # D's is not
    result = api.delineate(gauge_request(api, gauge_api), valley_repository(z, NODATA))
    assert np.array_equal(on_lattice(result, z.shape), np.asarray(expected.mask))
    g = result.gauge
    assert g.downstream_checked in ("partly", "none")
    assert "downstream_unread" in g.sensitivity.causes
    assert g.sensitivity.well_posed is False


def test_nodata_in_the_placed_catchment_is_refused(api: Any, gauge_api: Any) -> None:
    z = gf.valley(dam=True)
    z[UPPER_NODATA] = NODATA
    with pytest.raises(api.CatchmentError, match=r"(?i)nodata"):
        api.delineate(gauge_request(api, gauge_api), valley_repository(z, NODATA))


def test_stage_b_grows_the_window_past_stage_a_and_equals_a_whole_raster_run(
    api: Any, gauge_api: Any
) -> None:
    """The placed node's catchment lies within about 700 m of the river, but
    D's (30 m below, past the tributary's join) reaches 5.4 km east: the final
    window must reach it, and the result is still the placed node's."""
    z = gf.two_basins()
    placed, below = (150, gf.BIG_CC), (153, gf.BIG_CC)
    chain = [(r, gf.BIG_CC) for r in range(120, 161)]
    assert np.array_equal(gf.hand_burn(z, chain), z)  # the premise: nothing to lower
    expected, d_side = whole_flood(z, placed), whole_flood(z, below)
    assert expected.col_max < gf.BIG_CC + 100 < d_side.col_max - 400  # the premise
    request = gauge_request(api, gauge_api, col=gf.BIG_CC, row0=120, row1=160)
    result = api.delineate(request, valley_repository(z))
    m = result.meta
    last_col = round((m.x_min + (m.cols - 1) * m.delta_x - gf.X0) / gf.D)
    assert last_col >= d_side.col_max
    assert np.array_equal(on_lattice(result, z.shape), np.asarray(expected.mask))
    g = result.gauge
    assert g.downstream_checked == "whole"
    count = np.asarray(accumulate(core_view(gf.tile_of(z))).count)
    assert g.sensitivity.area_down == pytest.approx(count[below] * gf.CELL_KM2, rel=1e-12)
    assert "swing" in g.sensitivity.causes


#: The jog of `valley_with_a_gap`: from this row on, down to the exit slot,
#: the valley floor runs two columns east, and the node between the two runs
#: has no data.
JOG_ROW = 170
JOG_GAP = (JOG_ROW, gf.CC + 1)


def valley_with_a_gap(gap: float = NODATA) -> np.ndarray:
    """`gf.valley(dam=True)` with its floor moved two columns east on rows
    `JOG_ROW` to the slot (the same 3 m walls), and `gap` at `JOG_GAP`."""
    z = gf.valley(dam=True)
    rows = np.arange(JOG_ROW, gf.SLOT_ROW)
    cols = np.arange(gf.CC - 5, gf.CC + 8)
    z[np.ix_(rows, cols)] = gf.floor_z(rows)[:, None] + 3.0 * np.abs(cols - (gf.CC + 2))[None, :]
    z[JOG_GAP] = gap
    return z


def no_data(z: np.ndarray, nodata: float | None) -> np.ndarray:
    """True where a node has no data: NaN, or the sentinel when there is one."""
    return np.isnan(z) | (False if nodata is None else z == nodata)


#: (the value at the gap, the DEM's `nodata`): the sentinel, and a NaN cell
#: with no sentinel, as a NaN-gapped float DEM arrives (PR 2's code review,
#: round 2: a NaN cell is NoData).
@pytest.mark.parametrize(
    ("gap", "nodata"), [(NODATA, NODATA), (np.nan, None)], ids=["sentinel", "nan-no-sentinel"]
)
def test_a_chain_through_nodata_below_the_placed_node_is_refused(
    api: Any, gauge_api: Any, gap: float, nodata: float | None
) -> None:
    """The line runs down the valley's column past the jog, with a 30 m
    corridor: rows `JOG_ROW - 1` and `JOG_ROW` choose the floor nodes
    (169, CC) and (170, CC + 2), and the straight join between them steps
    through `JOG_GAP`, NoData (test_burn.py's
    `test_a_chain_through_nodata_is_refused_naming_the_gap` pins that step on
    a small terrain). The gap is 200 m below the placed node (row 150) and
    past `D` (30 m below), so neither the seed nor any flood meets it first:
    the refusal is the burn's, passed on as a `CatchmentError`."""
    z = valley_with_a_gap(gap)
    # The premises: the floor nodes either side of the join, from the fixture.
    # A cross-section is the row's nodes within 30 m of column CC.
    section = slice(gf.CC - 3, gf.CC + 4)
    assert int(np.argmin(z[JOG_ROW - 1, section])) == 3  # (169, CC)
    row = z[JOG_ROW, section]
    with_data = np.where(no_data(row, nodata), np.inf, row)
    assert int(np.argmin(with_data)) == 5  # (170, CC + 2)
    assert no_data(z, nodata)[JOG_GAP]
    assert int(no_data(z, nodata).sum()) == 1  # the gap is the only node without data
    assert not whole_flood(z, (150, gf.CC), nodata).touches_nodata
    assert not whole_flood(z, (153, gf.CC), nodata).touches_nodata
    reach = gauge_api.Reach(
        line=gf.column_line(gf.CC, 40, 178), at=110 * gf.D, uncertainty=30.0, corridor=30.0
    )
    request = api.CatchmentRequest(seed=gf.lat(gf.CC + 2, 150), seed_crs=gf.EPSG, reach=reach)
    with pytest.raises(api.CatchmentError, match=r"(?i)gap \(NoData\) in the DEM") as raised:
        api.delineate(request, valley_repository(z, nodata))
    x, y = gf.lat(JOG_GAP[1], JOG_GAP[0])
    message = str(raised.value)
    assert f"{x:.0f}" in message and f"{y:.0f}" in message, message


# ---------------------------------------------------------------------------
# Only the burn's own refusals become a CatchmentError (code review, round 2)
# ---------------------------------------------------------------------------


def burn_raising(exc: Exception) -> Any:
    def raising(*_: Any) -> Any:
        raise exc

    return raising


def test_a_burn_refusal_is_passed_on_as_a_catchment_error_with_its_words(
    api: Any, gauge_api: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    import tin_engine.burn as burn

    refusal = burn.BurnRefusal("the river line crosses a gap (NoData) in the DEM at (1, 2)")
    monkeypatch.setattr(api, "burn_reach", burn_raising(refusal))
    with pytest.raises(api.CatchmentError) as raised:
        api.delineate(gauge_request(api, gauge_api), valley_repository(gf.valley(dam=True)))
    assert str(raised.value) == str(refusal)
    assert raised.value.__cause__ is refusal


def test_any_other_value_error_from_the_burn_is_not_a_catchment_error(
    api: Any, gauge_api: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Stage B's `except CatchmentError` would swallow an ordinary bug turned
    into one; a plain `ValueError` from inside `burn_reach` must reach the
    caller as itself."""
    bug = ValueError("an ordinary bug inside burn_reach")
    monkeypatch.setattr(api, "burn_reach", burn_raising(bug))
    with pytest.raises(ValueError, match="an ordinary bug inside burn_reach") as raised:
        api.delineate(gauge_request(api, gauge_api), valley_repository(gf.valley(dam=True)))
    assert raised.value is bug
    assert not isinstance(raised.value, api.CatchmentError)


# ---------------------------------------------------------------------------
# Increment 29, PR 4: the joined causes, and the mixed-grid refusal by type
#
# "The window loop and the catchment": PR 4, as the second reader, moves the
# join of the `direction` cause into `catchment.py`: one `causes` on
# `GaugeResult`, the sensitivity's own causes in their order, then
# `"direction"` when the line runs against the DEM's slope (the order
# `catchment --rivers` wrote them in PR 2). "The batch": the mosaic's
# mixed-lattice refusal reaches the caller of `delineate` as
# `MixedGridRefusal(CatchmentError)`, with the mosaic error's `tiles` and its
# words.
# ---------------------------------------------------------------------------


def test_a_well_posed_gauge_has_no_joined_cause(api: Any, gauge_api: Any) -> None:
    result = api.delineate(gauge_request(api, gauge_api), valley_repository(gf.valley(dam=True)))
    g = result.gauge
    assert g.direction_ok
    assert tuple(g.causes) == () == tuple(g.sensitivity.causes)


def test_a_line_against_the_slope_joins_direction_after_the_sensitivitys_causes(
    api: Any, gauge_api: Any
) -> None:
    """The valley's line digitised upstream (row 200 to row 40): measured on
    PR 2's code, the burn's direction check fails and the sensitivity raises
    `chain_not_draining` and `chain_end_open`; the joined list adds
    `direction` after them, and the sensitivity's own list stays without it."""
    reach = gauge_api.Reach(
        line=gf.column_line(gf.CC, 200, 40), at=500.0, uncertainty=30.0, corridor=5.0
    )
    request = api.CatchmentRequest(seed=gf.lat(gf.CC + 2, 150), seed_crs=gf.EPSG, reach=reach)
    g = api.delineate(request, valley_repository(gf.valley(dam=True))).gauge
    assert not g.direction_ok
    assert "direction" not in g.sensitivity.causes
    assert tuple(g.causes) == (*g.sensitivity.causes, "direction")


def test_a_window_on_two_grids_is_a_mixed_grid_refusal_naming_a_tile_of_each(
    api: Any, gauge_api: Any
) -> None:
    from batch_fixtures import mixed_tiles

    reach = gauge_api.Reach(
        line=(gf.lat(690, gf.JOIN_ROW), gf.lat(500, gf.JOIN_ROW)), at=900.0, uncertainty=30.0
    )
    request = api.CatchmentRequest(
        seed=gf.lat(600, gf.JOIN_ROW - 1.6), seed_crs=gf.EPSG, reach=reach
    )
    with pytest.raises(api.MixedGridRefusal) as raised:
        api.delineate(request, MemoryRepository(mixed_tiles()))
    exc = raised.value
    assert isinstance(exc, api.CatchmentError)
    assert set(exc.tiles) == {"a.tif", "b.tif"}
    assert str(exc).startswith("the request selects tiles on two lattices, ")
    assert isinstance(exc.__cause__, mosaic.MixedGridError)
    assert tuple(exc.tiles) == tuple(exc.__cause__.tiles)


def test_other_refusals_of_the_plan_are_not_mixed_grid(api: Any) -> None:
    """A seed whose box meets no tile is refused, but not as mixed grids."""
    far = api.CatchmentRequest(seed=lat(5000, 5000), seed_crs=EPSG)
    with pytest.raises(api.CatchmentError) as raised:
        api.delineate(far, repository_of(bowl()))
    assert not isinstance(raised.value, api.MixedGridRefusal)
