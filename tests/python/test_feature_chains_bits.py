"""Feature chains and their bits: increment 16b's invariant-critical Python suite.

`docs/increments/16b-terrain-polygons.md`, "Tests for @tester": I1 (bits from
the input), I2 (the clip keeps input vertices exact) and I5 (ring shape), with
the oracle built from the input features and the class map, on synthetic
partitions and on the committed CORINE extract, through the real `_engine`
(`build_pslg` -> `node` -> `triangulate`). The oracle is
`feature_fixtures.expected_masks`: for every noded constraint edge, the union
of the map masks of every feature whose **source** boundary, moved by pyproj
directly, lies along it within the snap. It never reads the clip's output, the
chains, or the noder's provenance.

**Mutation testing is required here** (README, "Cost constraints"). Mutants
killed against a scratch implementation that is not committed, each named at
the test that kills it (`MUTANT n`):

1. area clipping instead of line clipping (the domain boundary gets bits);
2. a doubled shared edge keeping only one side's bit (a Python-side merge
   that keeps the first chain's mask, "no union");
3. a wrong class-map entry (`corine`'s `52x` treated as land only; `511`
   missing);
4. read order not by primary key (the `ORDER BY` dropped);
5. a clipped ring closed by repeating its first index;
6. holes' rings dropped.

Pinned by this suite (with `test_feature_input.py`; see "Pinned by the red
suite (16b-1/2)"):

- `chains.start_chains(domain, features, vocabulary) -> StartChains`, with
  `features` an iterable of `feature_input.TerrainFeature`. `StartChains` has
  `vertices` (`(N, 2)` float64, the DEM's CRS) and `chains`, a sequence of
  `(indices, role, mask)` with `role` one of the strings `"outer"`, `"hole"`,
  `"breakline"`: the domain's exterior first, then its holes (mask 0), then
  every feature's lines in feature order.
- A closed feature line (an unclipped ring) is a breakline whose last index
  repeats its first; an open one's does not.
- `TerrainFeature` has `fid`, `mask` and `lines` (shapely `LineString`s in the
  DEM's CRS, clipped).

HOW THIS FILE GOES RED: `tin_engine.feature_input` and `tin_engine.chains` are
imported lazily (`feature_fixtures`), so each test fails on its own.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from shapely.geometry import LineString, Polygon

from feature_fixtures import (
    CROSSING,
    UTM33,
    Feat,
    Noded,
    all_vertices,
    at,
    corine_mask,
    domain_of,
    expected_masks,
    extract_features,
    mismatches,
    moved,
    on_boundary,
    open_one,
    run_engine,
    square,
    start,
    truths,
    within,
    write_geojson,
)
from gpkg_fixtures import EXTRACT, EXTRACT_TABLE, copy_reversed, needs_rtree
from test_cli_mesh_domain import quarter_circle
from tin_engine.domain import DomainPolygon
from tin_engine.features import DEFAULT_VOCABULARY

V = DEFAULT_VOCABULARY
LAAEA = "EPSG:3035"

# ----------------------------------------------------------- the partitions

#: 2 x 2 squares of 100 m, one property each, so every shared edge is two
#: different bits and a union is visible.
GRID = {
    "sw": (square(0, 0, 100, 100), "road"),
    "se": (square(100, 0, 200, 100), "river"),
    "ne": (square(100, 100, 200, 200), "wall"),
    "nw": (square(0, 100, 100, 200), "ditch"),
}


def grid_features() -> list[Feat]:
    return [Feat(name, poly, {"property": prop}) for name, (poly, prop) in GRID.items()]


def grid_truth() -> list[Any]:
    return truths((poly, V.mask(prop)) for poly, prop in GRID.values())


#: Contains the grid: every ring enters whole.
AROUND = domain_of(square(-50, -50, 250, 250))
#: Cuts the grid: its west side lies on the grid's west edge from y 0 to 200
#: (a feature edge along the boundary), its east side at x 150 slices the two
#: eastern squares, and a hole round the grid's centre cuts all four.
CUTTING = domain_of(
    Polygon(
        [at(0, -50), at(150, -50), at(150, 250), at(0, 250)],
        [[at(80, 80), at(80, 120), at(120, 120), at(120, 80)]],
    )
)

#: A square with a hole, and a second polygon filling the hole.
ISLAND_OUTER = Polygon(
    square(0, 0, 200, 200).exterior.coords, [square(60, 60, 140, 140).exterior.coords]
)
ISLAND_INNER = square(60, 60, 140, 140)

#: A lake, and a river ending on its shore at a T (mid-edge, not at a vertex).
LAKE = square(50, 50, 150, 150)
RIVER = LineString([at(0, 100.25), at(30, 97.5), at(50, 100.25)])


def mesh_of(tmp_path: Path, features: list[Feat], domain: DomainPolygon) -> tuple[Any, Any, Noded]:
    fs = open_one(write_geojson(tmp_path / "f.geojson", features), domain)
    started = start(domain, fs.features)
    return fs, started, run_engine(started)


def chain_roles(started: Any) -> list[str]:
    return [role for _, role, _ in started.chains]


def feature_chains(started: Any) -> list[tuple[list[int], str, int]]:
    return [(list(map(int, i)), r, int(m)) for i, r, m in started.chains if r == "breakline"]


# ------------------------------------------------------------------- I1


class TestI1BitsFromTheInput:
    """Every noded constraint edge's mask is exactly the union of the masks of
    the features whose linework lies along it; 0 where none does."""

    @pytest.mark.parametrize("domain", [AROUND, CUTTING], ids=["around", "cutting"])
    def test_the_grid(self, tmp_path: Path, domain: DomainPolygon) -> None:
        """MUTANT 1 (the cutting domain's x = 150 side and hole get bits under
        an area clip), MUTANT 2 (a shared edge keeps one side's bit)."""
        _, _, noded = mesh_of(tmp_path, grid_features(), domain)
        expected = expected_masks(noded, grid_truth())
        assert mismatches(noded, expected) == []

    def test_a_shared_edge_carries_both_sides(self, tmp_path: Path) -> None:
        """MUTANT 2, stated directly: the edge between `sw` (road) and `se`
        (river), x = 100 from y = 0 to 100, carries road | river."""
        _, _, noded = mesh_of(tmp_path, grid_features(), AROUND)
        shared = [
            m
            for ((ax, ay), (bx, by)), m in noded.masks.items()
            if ax == bx and abs(ax - at(100, 0)[0]) < 1e-3 and min(ay, by) < at(0, 100)[1] - 1
        ]
        assert shared, "no edge on the shared side was found"
        assert set(shared) == {V.mask("road", "river")}

    def test_no_boundary_edge_carries_a_bit_where_no_feature_edge_lies(
        self, tmp_path: Path
    ) -> None:
        """MUTANT 1, stated directly (the design's "direct assertion"): the
        cutting domain's east side (x = 150) and its hole cross features but
        lie along none, so every `Outer` and `Hole` edge there has mask 0."""
        _, _, noded = mesh_of(tmp_path, grid_features(), CUTTING)
        east, hole = at(150, 0)[0], shapely.box(*at(80, 80), *at(120, 120)).boundary
        crossing_only = [
            (e, noded.masks[e])
            for e in noded.boundary
            if (abs(e[0][0] - east) < 1e-3 and abs(e[1][0] - east) < 1e-3)
            or LineString(e).distance(hole) < 1e-3
        ]
        assert crossing_only, "found no boundary edge to check"
        assert [(e, m) for e, m in crossing_only if m != 0] == []

    def test_a_feature_edge_along_the_boundary_gives_it_its_bits(self, tmp_path: Path) -> None:
        """R6: a piece lying along the domain boundary is kept, and the noder's
        union gives that boundary edge the feature's bits."""
        _, _, noded = mesh_of(tmp_path, grid_features(), CUTTING)
        west = at(0, 0)[0]
        along = {
            e: noded.masks[e]
            for e in noded.boundary
            if abs(e[0][0] - west) < 1e-3 and abs(e[1][0] - west) < 1e-3
        }
        y = lambda e: (e[0][1] + e[1][1]) / 2  # noqa: E731
        assert {m for e, m in along.items() if at(0, 0)[1] < y(e) < at(0, 100)[1]} == {
            V.mask("road")
        }
        assert {m for e, m in along.items() if at(0, 100)[1] < y(e) < at(0, 200)[1]} == {
            V.mask("ditch")
        }
        assert {m for e, m in along.items() if y(e) < at(0, 0)[1] or y(e) > at(0, 200)[1]} == {0}

    def test_a_hole_ring_shared_with_the_polygon_filling_it(self, tmp_path: Path) -> None:
        """MUTANT 6: the outer polygon's hole ring carries its bit too, so the
        hole's edges are land_cover | water, not water alone."""
        features = [
            Feat("outer", ISLAND_OUTER, {"property": "land_cover"}),
            Feat("inner", ISLAND_INNER, {"property": "water"}),
        ]
        _, started, noded = mesh_of(tmp_path, features, AROUND)
        truth = truths([(ISLAND_OUTER, V.mask("land_cover")), (ISLAND_INNER, V.mask("water"))])
        assert mismatches(noded, expected_masks(noded, truth)) == []
        assert V.mask("land_cover", "water") in set(noded.masks.values())
        assert len(feature_chains(started)) == 3  # exterior, hole, the inner polygon

    def test_a_river_ending_on_the_lake_shore(self, tmp_path: Path) -> None:
        """A T-junction the noder resolves: the shore splits at the river's end
        and keeps `water` on both halves; the river's edges are `river` only."""
        features = [
            Feat("lake", LAKE, {"property": "water"}),
            Feat("river", RIVER, {"property": "river"}),
        ]
        _, _, noded = mesh_of(tmp_path, features, AROUND)
        truth = truths([(LAKE, V.mask("water")), (RIVER, V.mask("river"))])
        assert mismatches(noded, expected_masks(noded, truth)) == []
        end = shapely.Point(RIVER.coords[-1])
        at_end = [
            m
            for e, m in noded.masks.items()
            if min(shapely.Point(p).distance(end) for p in e) < 1e-3
        ]
        assert sorted(at_end) == sorted([V.mask("water"), V.mask("water"), V.mask("river")])


# ------------------------------------------------------------------- I2


class TestI2InputVerticesExact:
    def test_every_vertex_is_a_source_vertex_or_a_crossing(self, tmp_path: Path) -> None:
        """I2 on the cutting domain: bit for bit a source vertex, or within
        1e-6 m of the domain's boundary; never outside its closure."""
        fs, _, _ = mesh_of(tmp_path, grid_features(), CUTTING)
        xy = all_vertices([line for f in fs.features for line in f.lines])
        source = {p for poly, _ in GRID.values() for p in poly.exterior.coords}
        polygon = CUTTING.polygon
        for x, y in xy:
            if (x, y) not in source:
                assert on_boundary(polygon, np.array([[x, y]]))[0] <= CROSSING, (x, y)
        assert within(polygon, xy).max() <= CROSSING

    def test_reprojected_vertices_are_pyprojs_bit_for_bit(self, tmp_path: Path) -> None:
        """The grid given in EPSG:4326: every vertex that is not a crossing is
        pyproj's `always_xy` image of a source vertex, exactly."""
        utm = {name: (poly, prop) for name, (poly, prop) in GRID.items()}
        lonlat = [
            Feat(name, moved(poly, UTM33, "EPSG:4326"), {"property": prop})
            for name, (poly, prop) in utm.items()
        ]
        path = write_geojson(tmp_path / "f.geojson", lonlat, crs=None)
        fs = open_one(path, CUTTING)
        images = {
            p
            for f in lonlat
            for p in moved(f.geometry, "EPSG:4326", UTM33).exterior.coords  # type: ignore[arg-type,union-attr]
        }
        xy = all_vertices([line for f in fs.features for line in f.lines])
        off = [(x, y) for x, y in xy if (x, y) not in images]
        assert len(off) < len(xy)
        assert on_boundary(CUTTING.polygon, np.array(off)).max() <= CROSSING


# ------------------------------------------------------------------- I5


class TestI5RingShape:
    def test_domain_first_then_features_and_no_feature_is_outer_or_hole(
        self, tmp_path: Path
    ) -> None:
        _, started, _ = mesh_of(tmp_path, grid_features(), CUTTING)
        roles = chain_roles(started)
        assert roles[:2] == ["outer", "hole"]
        assert set(roles[2:]) == {"breakline"}
        assert [int(m) for _, _, m in started.chains[:2]] == [0, 0]

    def test_an_unclipped_ring_is_one_closed_breakline(self, tmp_path: Path) -> None:
        """Every grid ring inside `AROUND`: one breakline each, its last index
        its first, and its vertices the ring's."""
        fs, started, _ = mesh_of(tmp_path, grid_features(), AROUND)
        closed = feature_chains(started)
        assert len(closed) == 4
        xy = np.asarray(started.vertices)
        for (indices, _, mask), (poly, prop) in zip(closed, GRID.values(), strict=True):
            assert indices[0] == indices[-1] and len(set(indices)) == len(indices) - 1
            assert {tuple(xy[i]) for i in indices} == set(poly.exterior.coords)
            assert mask == V.mask(prop)
        assert [f.mask for f in fs.features] == [V.mask(p) for _, p in GRID.values()]

    def test_a_clipped_ring_is_open_breaklines(self, tmp_path: Path) -> None:
        """MUTANT 5: under `CUTTING` every ring is cut, so no feature chain may
        repeat its first index. A piece ends on the domain boundary or at a
        source vertex (where GEOS split it at the ring's own start)."""
        _, started, _ = mesh_of(tmp_path, grid_features(), CUTTING)
        xy = np.asarray(started.vertices)
        pieces = feature_chains(started)
        corners = {p for poly, _ in GRID.values() for p in poly.exterior.coords}
        assert pieces
        for indices, _, _ in pieces:
            assert indices[0] != indices[-1]
            for end in (tuple(xy[indices[0]]), tuple(xy[indices[-1]])):
                assert (
                    end in corners or on_boundary(CUTTING.polygon, np.array([end]))[0] <= CROSSING
                )

    def test_a_line_is_an_open_breakline(self, tmp_path: Path) -> None:
        _, started, _ = mesh_of(tmp_path, [Feat("r", RIVER, {"property": "river"})], AROUND)
        ((indices, _, mask),) = feature_chains(started)
        assert indices[0] != indices[-1] and len(indices) == 3
        assert mask == V.mask("river")

    def test_the_domain_half_is_the_domains_rings(self) -> None:
        """No features: the domain alone, as `cli._domain_chains` gave it
        before 16b (`chains.py` takes it over, R2)."""
        started = start(CUTTING, ())
        assert chain_roles(started) == ["outer", "hole"]
        xy = np.asarray(started.vertices)
        rings = [CUTTING.polygon.exterior, *CUTTING.polygon.interiors]
        for (indices, _, _), ring in zip(started.chains, rings, strict=True):
            assert [tuple(xy[int(i)]) for i in indices] == list(ring.coords)[:-1]


# ------------------------------------------------------------ order (R3, R7)


class TestReadOrder:
    @needs_rtree
    def test_the_extract_in_reverse_row_order_gives_the_same_chains(self, tmp_path: Path) -> None:
        """MUTANT 4: I3 at the chains. The extract rewritten with its rows and
        R-tree entries inserted in reverse primary-key order gives the same
        `StartChains`, bit for bit, because rows are read `ORDER BY` the key."""
        domain = domain_of(Polygon(quarter_circle()))
        reversed_path = copy_reversed(EXTRACT, tmp_path / "reversed.gpkg")
        a = open_one(EXTRACT, domain, "corine")
        b = open_one(reversed_path, domain, "corine")
        assert [f.fid for f in b.features] == [f.fid for f in a.features]
        assert [f.fid for f in a.features] == sorted(f.fid for f in a.features)
        sa, sb = start(domain, a.features), start(domain, b.features)
        assert np.array_equal(np.asarray(sa.vertices), np.asarray(sb.vertices))
        assert [(list(map(int, i)), r, int(m)) for i, r, m in sa.chains] == [
            (list(map(int, i)), r, int(m)) for i, r, m in sb.chains
        ]


# ------------------------------------------------------------ the extract


class TestCommittedExtract:
    """CORINE over the benchmark tile (EPSG:3035), the quarter circle
    (EPSG:25833), the `corine` map, through the real engine."""

    @pytest.fixture(scope="class")
    def run(self) -> tuple[Any, Any, Noded, DomainPolygon]:
        domain = domain_of(Polygon(quarter_circle()))
        fs = open_one(EXTRACT, domain, "corine")
        started = start(domain, fs.features)
        return fs, started, run_engine(started), domain

    def test_i1_every_edge_carries_exactly_its_features_bits(self, run: Any) -> None:
        """MUTANTS 1, 2 and 3 on real data. Coastal edges are sea (`523`) on
        one side and land on the other: land_cover | water."""
        _, _, noded, _ = run
        rows = extract_features(EXTRACT, EXTRACT_TABLE)
        truth = truths([(g, corine_mask(code)) for _, g, code in rows], src=LAAEA)
        assert mismatches(noded, expected_masks(noded, truth)) == []

    def test_the_bits_that_occur(self, run: Any) -> None:
        """Every feature edge carries `land_cover`; the coast and the rivers
        carry `water` as well; the domain boundary away from features, none."""
        _, _, noded, _ = run
        found = set(noded.masks.values())
        assert V.mask("land_cover") in found
        assert V.mask("land_cover", "water") in found
        assert found <= {0, V.mask("land_cover"), V.mask("land_cover", "water")}
        assert 0 in {noded.masks[e] for e in noded.boundary}

    def test_the_masks_are_the_corine_maps(self, run: Any) -> None:
        """MUTANT 3, per feature: each `TerrainFeature.mask` is R4's for its
        row's `Code_18`, read from the file independently."""
        fs, _, _, _ = run
        codes = {pk: code for pk, _, code in extract_features(EXTRACT, EXTRACT_TABLE)}
        assert [f.mask for f in fs.features] == [corine_mask(codes[f.fid]) for f in fs.features]
        assert {codes[f.fid] for f in fs.features} >= {"512", "523"}

    def test_i2_vertices_are_reprojected_source_vertices_or_crossings(self, run: Any) -> None:
        fs, _, _, domain = run
        rows = extract_features(EXTRACT, EXTRACT_TABLE)
        images = {
            p
            for _, g, _ in rows
            for line in shapely.get_parts(moved(g, LAAEA).boundary)
            for p in line.coords
        }
        xy = all_vertices([line for f in fs.features for line in f.lines])
        off = np.array([(x, y) for x, y in xy if (x, y) not in images])
        assert len(off) < len(xy) / 10
        assert on_boundary(domain.polygon, off).max() <= CROSSING
        assert within(domain.polygon, xy).max() <= CROSSING

    def test_i5_the_extracts_chains(self, run: Any) -> None:
        """A chain is closed exactly when it never reaches the domain boundary
        (an unclipped ring); every piece that does is open."""
        _, started, noded, domain = run
        roles = chain_roles(started)
        assert roles[0] == "outer" and set(roles[1:]) == {"breakline"}
        xy = np.asarray(started.vertices)
        for indices, _, _ in feature_chains(started):
            touches = on_boundary(domain.polygon, xy[indices]).min() <= CROSSING
            assert (indices[0] == indices[-1]) == (not touches)
        assert noded.status == "Ok"
