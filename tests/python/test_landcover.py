"""`tin_engine.landcover`: a land-cover class per triangle (increment 16c, R1).

`docs/increments/16c-landcover-labels.md`, "Tests for @tester": `regions` and
`label_triangles` on hand-made meshes, no `_core` call, no file. Every mesh
here is a rectilinear grid whose lines carry every input boundary, so the
constraint edges are exactly the grid edges lying on the input linework (the
domain's boundary, each polygon's rings, each road), found by distance and not
by the producer. Coordinates sit at UTM 33N magnitudes, off round numbers.

Each labelling fixture is checked three ways: against the centroid oracle of
`landcover_fixtures` (I2), for I1 (the spread), and against the codes the
design names for it. The fixtures are the design's CLI list as pure meshes:
two squares side by side, a forest with an empty hole, a lake in a holed and
in an unholed forest, a polygon clipped by the domain into two pieces, a road
across a polygon, a polygon covering the domain, boundaries 0.1 mm apart,
determinism, and overlaps (smallest area, ties to the smaller code; D2).

Pinned beyond the design's text (the design leaves them open):

- `label_triangles(vertices, triangles, edges, polygons=..., margin=...)`,
  `vertices` `(N, 3)` (the trimmed mesh's; only x and y are used), and
  `polygons` a sequence of `(geometry, code)` pairs.
- `CoverLabels.regions`, `.outside`, `.overlapped` and `.thin` are counts
  (ints), the numbers of R3's stderr line; `.codes` is the `(T,)` array.
- `regions(triangles, edges)` takes no vertex count.

Committed red at `196147e`: `tin_engine.landcover` did not exist yet, so
every test failed on `ModuleNotFoundError`. The module landed in `0487ed0`
and the suite has been green since.

Increment 30a (`docs/increments/30a-landcover-speed.md`, section 7) added the
classes from `TestBest` on: tests for the new helpers `_best` and `_member`,
and pins that were already green before the rewrite. Committed red at
`6cde0eb`: the helpers did not exist yet, so their tests failed and the pins
passed. The helpers landed in `3066d60` and the suite has been green since.
Pinned beyond that design's text: `_member` returns a bool array and
takes empty `values`; `_best` and `_member` take numpy arrays, positionally;
`_lookup(points, polygons)` keeps its signature and returns `(codes, hits)`;
P4 also holds with vertices present and no triangles.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pytest
import shapely
from numpy.testing import assert_array_equal
from shapely.geometry import LineString, MultiLineString, Polygon, box
from shapely.geometry.base import BaseGeometry

from landcover_fixtures import (
    MARGIN,
    edge_keys,
    interior_edge_count,
    landcover_oracle,
    spread_violations,
)
from tin_engine import landcover

X0, Y0 = 500_000.3, 6_600_000.7


def rect(x0: float, y0: float, x1: float, y1: float) -> Polygon:
    return box(X0 + x0, Y0 + y0, X0 + x1, Y0 + y1)


def line(*points: tuple[float, float]) -> LineString:
    return LineString([(X0 + x, Y0 + y) for x, y in points])


# ---------------------------------------------------------------- the meshes


@dataclass(frozen=True)
class Mesh:
    """A trimmed mesh as the writers see it: `(N, 3)` vertices, triangles, and
    the constraint edges."""

    vertices: np.ndarray
    triangles: np.ndarray
    edges: np.ndarray

    @property
    def centroids(self) -> np.ndarray:
        return self.vertices[self.triangles][:, :, :2].mean(axis=1)

    def inside(self, geometry: BaseGeometry) -> np.ndarray:
        """Triangles whose centroid is strictly inside `geometry`."""
        return np.asarray(shapely.contains_xy(geometry, *self.centroids.T))


def grid_mesh(
    xs: Sequence[float],
    ys: Sequence[float],
    lines: Sequence[BaseGeometry],
    domain: Polygon | None = None,
) -> Mesh:
    """A rectilinear grid, each cell split on its rising diagonal, keeping the
    cells whose centre is in `domain` (default: the whole grid). The
    constraint edges are the grid edges lying on `lines` or on the domain's
    boundary, within a micrometre at both ends and the midpoint."""
    gx, gy = np.meshgrid(np.asarray(xs) + X0, np.asarray(ys) + Y0, indexing="xy")
    vertices = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
    nx = len(xs)
    if domain is None:
        domain = box(X0 + xs[0], Y0 + ys[0], X0 + xs[-1], Y0 + ys[-1])
    triangles = []
    for j in range(len(ys) - 1):
        for i in range(nx - 1):
            v00, v10 = j * nx + i, j * nx + i + 1
            v01, v11 = v00 + nx, v10 + nx
            centre = vertices[[v00, v11], :2].mean(axis=0)
            if domain.contains(shapely.Point(*centre)):
                triangles += [(v00, v10, v11), (v00, v11, v01)]
    tri = np.array(triangles, dtype=np.int64)
    sides = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    _, first = np.unique(edge_keys(sides, len(vertices)), return_index=True)
    candidates = sides[np.sort(first)]
    linework = shapely.union_all([domain.boundary, *(_linework(g) for g in lines)])
    a, b = vertices[candidates[:, 0], :2], vertices[candidates[:, 1], :2]
    near = [
        np.asarray(shapely.distance(linework, shapely.points(p))) <= 1e-6
        for p in (a, b, (a + b) / 2)
    ]
    on = near[0] & near[1] & near[2]
    return Mesh(vertices=vertices, triangles=tri, edges=candidates[on])


def _linework(geometry: BaseGeometry) -> BaseGeometry:
    if isinstance(geometry, Polygon):
        return MultiLineString([r.coords for r in (geometry.exterior, *geometry.interiors)])
    return geometry


def label(mesh: Mesh, polygons: Sequence[tuple[BaseGeometry, int]], margin: float = MARGIN) -> Any:
    return landcover.label_triangles(
        mesh.vertices, mesh.triangles, mesh.edges, polygons=list(polygons), margin=margin
    )


def checked_label(
    mesh: Mesh, polygons: Sequence[tuple[BaseGeometry, int]], margin: float = MARGIN
) -> Any:
    """`label_triangles`, with I1 and I2 asserted on the way out."""
    labels = label(mesh, polygons, margin)
    codes = np.asarray(labels.codes)
    assert codes.shape == (len(mesh.triangles),)
    assert interior_edge_count(mesh.triangles) > 0
    assert spread_violations(mesh.triangles, mesh.edges, codes) == []
    oracle = landcover_oracle(mesh.vertices, mesh.triangles, polygons, margin)
    assert oracle.checked.any()
    assert oracle.mismatches(codes) == []
    return labels


STEPS = [float(v) for v in range(21)]
DOMAIN = rect(0, 0, 20, 20)


# ------------------------------------------------------------------- regions


def two_triangles() -> tuple[np.ndarray, np.ndarray]:
    return np.array([[0, 1, 2], [0, 2, 3]]), np.array([[0, 1], [1, 2], [2, 3], [3, 0]])


def strip(cells: int) -> np.ndarray:
    """`2 * cells` triangles in a row: bottom row 0..cells, top row after it."""
    top = cells + 1
    out = []
    for i in range(cells):
        out += [(i, i + 1, top + i + 1), (i, top + i + 1, top + i)]
    return np.array(out, dtype=np.int64)


class TestRegions:
    def test_an_unconstrained_shared_edge_joins_two_triangles(self) -> None:
        triangles, hull = two_triangles()
        assert_array_equal(landcover.regions(triangles, hull), [0, 0])

    def test_a_constrained_shared_edge_separates_them(self) -> None:
        triangles, hull = two_triangles()
        cut = np.vstack([hull, [[0, 2]]])
        assert_array_equal(landcover.regions(triangles, cut), [0, 1])

    def test_a_constraint_given_backwards_still_blocks(self) -> None:
        triangles, hull = two_triangles()
        cut = np.vstack([hull, [[2, 0]]])
        assert_array_equal(landcover.regions(triangles, cut), [0, 1])

    def test_no_constraints_at_all(self) -> None:
        triangles, _ = two_triangles()
        assert_array_equal(landcover.regions(triangles, np.zeros((0, 2), dtype=np.int64)), [0, 0])

    @pytest.mark.parametrize("order", ["forward", "reversed", "shuffled"])
    def test_a_strip_of_1000_triangles_is_one_component(self, order: str) -> None:
        """The union-find loops to a fixed point: a chain 1 000 long is one
        component whatever order its triangles are listed in."""
        triangles = strip(500)
        assert len(triangles) == 1000
        if order == "reversed":
            triangles = triangles[::-1]
        elif order == "shuffled":
            triangles = triangles[np.random.default_rng(16).permutation(len(triangles))]
        ids = np.asarray(landcover.regions(triangles, np.zeros((0, 2), dtype=np.int64)))
        assert ids.shape == (1000,)
        assert set(ids.tolist()) == {0}

    def test_ids_are_the_smallest_triangle_index_of_each_component(self) -> None:
        cut = line((10, 0), (10, 20))
        mesh = grid_mesh(STEPS, STEPS, [cut])
        order = np.random.default_rng(3).permutation(len(mesh.triangles))
        triangles = mesh.triangles[order]
        ids = np.asarray(landcover.regions(triangles, mesh.edges))
        west = np.asarray(shapely.contains_xy(rect(0, 0, 10, 20), *mesh.centroids[order].T))
        assert west.any() and (~west).any()
        assert set(ids[west].tolist()) == {int(np.flatnonzero(west).min())}
        assert set(ids[~west].tolist()) == {int(np.flatnonzero(~west).min())}

    def test_permuting_and_flipping_the_constraint_rows_changes_nothing(self) -> None:
        mesh = grid_mesh(STEPS, STEPS, [rect(4, 4, 16, 16), line((0, 10), (20, 10))])
        rng = np.random.default_rng(7)
        shuffled = mesh.edges[rng.permutation(len(mesh.edges))]
        flipped = shuffled[:, ::-1]
        first = np.asarray(landcover.regions(mesh.triangles, mesh.edges))
        assert len(set(first.tolist())) == 4
        assert_array_equal(landcover.regions(mesh.triangles, shuffled), first)
        assert_array_equal(landcover.regions(mesh.triangles, flipped), first)


# ---------------------------------------------------------- label_triangles


class TestLabelBasics:
    def test_codes_are_int32_one_per_triangle(self) -> None:
        mesh = grid_mesh(STEPS, STEPS, [rect(4, 4, 16, 16)])
        labels = label(mesh, [(rect(4, 4, 16, 16), 311)])
        codes = np.asarray(labels.codes)
        assert codes.dtype == np.int32
        assert codes.shape == (len(mesh.triangles),)

    def test_a_square_cut_by_a_constrained_diagonal(self) -> None:
        """Each side of the diagonal gets its own polygon's code."""
        diagonal = line((0, 0), (20, 20))
        lower = Polygon([(X0, Y0), (X0 + 20, Y0), (X0 + 20, Y0 + 20)])
        upper = Polygon([(X0, Y0), (X0 + 20, Y0 + 20), (X0, Y0 + 20)])
        mesh = grid_mesh(STEPS, STEPS, [diagonal])
        labels = checked_label(mesh, [(lower, 311), (upper, 512)])
        codes = np.asarray(labels.codes)
        assert set(codes[mesh.inside(lower)].tolist()) == {311}
        assert set(codes[mesh.inside(upper)].tolist()) == {512}
        assert (labels.regions, labels.outside, labels.overlapped) == (2, 0, 0)

    def test_no_polygons_is_all_zero_and_every_region_outside(self) -> None:
        mesh = grid_mesh(STEPS, STEPS, [rect(4, 4, 16, 16), line((0, 10), (20, 10))])
        labels = checked_label(mesh, [])
        assert set(np.asarray(labels.codes).tolist()) == {0}
        assert labels.regions == 4
        assert labels.outside == labels.regions
        assert labels.overlapped == 0

    def test_a_sliver_component_is_counted_thin_and_still_labelled(self) -> None:
        """A 1 mm sliver beside a big triangle, the shared edge constrained:
        two components, the sliver's largest inradius about 0.5 mm <= margin."""
        vertices = np.array(
            [
                [X0, Y0, 0.0],
                [X0 + 100, Y0, 0.0],
                [X0 + 100, Y0 + 100, 0.0],
                [X0 + 100.001, Y0 + 50, 0.0],
            ]
        )
        triangles = np.array([[0, 1, 2], [1, 3, 2]])
        edges = np.array([[0, 1], [1, 3], [3, 2], [2, 0], [1, 2]])
        cover = rect(-10, -10, 110, 110)
        labels = landcover.label_triangles(
            vertices, triangles, edges, polygons=[(cover, 311)], margin=MARGIN
        )
        assert_array_equal(labels.codes, [311, 311])
        assert (labels.regions, labels.thin) == (2, 1)
        finer = landcover.label_triangles(
            vertices, triangles, edges, polygons=[(cover, 311)], margin=1e-5
        )
        assert finer.thin == 0


# ------------------------------------------------------------ the fixtures


class TestFixtures:
    def test_two_squares_side_by_side(self) -> None:
        west, east = rect(2, 4, 10, 16), rect(10, 4, 18, 16)
        mesh = grid_mesh(STEPS, STEPS, [west, east])
        labels = checked_label(mesh, [(west, 311), (east, 512)])
        codes = np.asarray(labels.codes)
        assert set(codes[mesh.inside(west)].tolist()) == {311}
        assert set(codes[mesh.inside(east)].tolist()) == {512}
        assert set(codes[~mesh.inside(west.union(east))].tolist()) == {0}
        # The triangles touching the shared side carry both codes.
        seam = line((10, 5), (10, 15))
        touching = (
            np.asarray(
                shapely.distance(seam, shapely.polygons(mesh.vertices[mesh.triangles][:, :, :2]))
            )
            == 0
        )
        assert set(codes[touching].tolist()) == {311, 512}
        assert (labels.regions, labels.outside, labels.overlapped, labels.thin) == (3, 1, 0, 0)

    def test_a_forest_with_an_empty_hole(self) -> None:
        hole = rect(8, 8, 12, 12)
        forest = Polygon(rect(2, 2, 18, 18).exterior, [hole.exterior])
        mesh = grid_mesh(STEPS, STEPS, [forest])
        labels = checked_label(mesh, [(forest, 312)])
        codes = np.asarray(labels.codes)
        assert set(codes[mesh.inside(hole)].tolist()) == {0}
        assert set(codes[mesh.inside(forest)].tolist()) == {312}
        assert (labels.regions, labels.outside) == (3, 2)

    def test_a_lake_filling_the_hole_of_a_holed_forest(self) -> None:
        lake = rect(8, 8, 12, 12)
        forest = Polygon(rect(2, 2, 18, 18).exterior, [lake.exterior])
        mesh = grid_mesh(STEPS, STEPS, [forest, lake])
        labels = checked_label(mesh, [(forest, 312), (lake, 512)])
        codes = np.asarray(labels.codes)
        assert set(codes[mesh.inside(lake)].tolist()) == {512}
        assert set(codes[mesh.inside(forest)].tolist()) == {312}
        assert labels.overlapped == 0

    def test_a_lake_inside_an_unholed_forest_goes_to_the_smaller(self) -> None:
        """D2: the lake's component is in both polygons; the lake is smaller."""
        lake = rect(8, 8, 12, 12)
        forest = rect(2, 2, 18, 18)
        mesh = grid_mesh(STEPS, STEPS, [forest, lake])
        labels = checked_label(mesh, [(forest, 312), (lake, 512)])
        codes = np.asarray(labels.codes)
        assert set(codes[mesh.inside(lake)].tolist()) == {512}
        assert set(codes[mesh.inside(forest.difference(lake))].tolist()) == {312}
        assert labels.overlapped == 1

    def test_a_polygon_clipped_by_the_domain_into_two_pieces(self) -> None:
        """A band crossing a notch cut into the domain from the north: its two
        pieces are separate components, and both carry its code."""
        domain = rect(0, 0, 20, 20).difference(rect(8, 10, 12, 20))
        band = rect(4, 12, 16, 16)
        mesh = grid_mesh(STEPS, STEPS, [band], domain=domain)
        labels = checked_label(mesh, [(band, 324)])
        codes = np.asarray(labels.codes)
        west, east = mesh.inside(rect(4, 12, 8, 16)), mesh.inside(rect(12, 12, 16, 16))
        assert west.any() and east.any()
        assert set(codes[west | east].tolist()) == {324}
        ids = np.asarray(landcover.regions(mesh.triangles, mesh.edges))
        assert set(ids[west].tolist()).isdisjoint(ids[east].tolist())
        assert (labels.regions, labels.outside) == (3, 1)

    def test_a_road_across_a_polygon(self) -> None:
        """The road is a constraint and no polygon: both halves are the forest's."""
        forest = rect(4, 4, 16, 16)
        road = line((2, 10), (18, 10))
        mesh = grid_mesh(STEPS, STEPS, [forest, road])
        labels = checked_label(mesh, [(forest, 311)])
        codes = np.asarray(labels.codes)
        north, south = mesh.inside(rect(4, 10, 16, 16)), mesh.inside(rect(4, 4, 16, 10))
        assert set(codes[north | south].tolist()) == {311}
        ids = np.asarray(landcover.regions(mesh.triangles, mesh.edges))
        assert set(ids[north].tolist()).isdisjoint(ids[south].tolist())
        assert (labels.regions, labels.outside) == (3, 1)

    def test_a_polygon_covering_the_domain(self) -> None:
        cover = rect(-5, -5, 25, 25)
        mesh = grid_mesh(STEPS, STEPS, [])
        labels = checked_label(mesh, [(cover, 333)])
        assert set(np.asarray(labels.codes).tolist()) == {333}
        assert (labels.regions, labels.outside, labels.overlapped) == (1, 0, 0)

    def test_boundaries_a_tenth_of_a_millimetre_apart(self) -> None:
        """Two squares whose shared side is offset by 0.1 mm, less than the
        snap, with nothing merged (the pure form of the CLI fixture, where
        the noder may merge them): the sliver column between the two sides
        opens onto the outside at both ends, so it is part of the outside
        component; away from the seam each square has its own code (I2)."""
        xs = [*STEPS[:11], 10.0001, *STEPS[11:]]
        west, east = rect(2, 4, 10, 16), rect(10.0001, 4, 18, 16)
        mesh = grid_mesh(xs, STEPS, [west, east])
        labels = checked_label(mesh, [(west, 311), (east, 512)])
        codes = np.asarray(labels.codes)
        assert set(codes[mesh.inside(rect(2, 4, 9.99, 16))].tolist()) == {311}
        assert set(codes[mesh.inside(rect(10.01, 4, 18, 16))].tolist()) == {512}
        assert (labels.regions, labels.outside, labels.thin) == (3, 1, 0)
        gap = mesh.inside(rect(10, 4.5, 10.0001, 15.5))
        assert gap.any() and set(codes[gap].tolist()) == {0}


class TestOverlaps:
    """D2: the smallest area wins, ties to the smaller code, in any order."""

    @pytest.mark.parametrize("reverse", [False, True])
    def test_equal_areas_go_to_the_smaller_code(self, reverse: bool) -> None:
        square = rect(4, 4, 16, 16)
        mesh = grid_mesh(STEPS, STEPS, [square])
        polygons = [(square, 412), (rect(4, 4, 16, 16), 322)]
        labels = checked_label(mesh, polygons[::-1] if reverse else polygons)
        assert set(np.asarray(labels.codes)[mesh.inside(square)].tolist()) == {322}
        assert labels.overlapped == 1

    @pytest.mark.parametrize("reverse", [False, True])
    def test_nested_polygons_give_the_inner_code(self, reverse: bool) -> None:
        inner, outer = rect(6, 6, 14, 14), rect(-1, -1, 21, 21)
        mesh = grid_mesh(STEPS, STEPS, [inner])
        polygons = [(outer, 312), (inner, 512)]
        labels = checked_label(mesh, polygons[::-1] if reverse else polygons)
        codes = np.asarray(labels.codes)
        assert set(codes[mesh.inside(inner)].tolist()) == {512}
        assert set(codes[~mesh.inside(inner)].tolist()) == {312}
        assert (labels.regions, labels.overlapped, labels.outside) == (2, 1, 0)


class TestDeterminism:
    """R1: the codes are a function of the triangles, the constraint edges and
    the polygons as a set."""

    @pytest.fixture
    def case(self) -> tuple[Mesh, list[tuple[BaseGeometry, int]]]:
        lake = rect(8, 8, 12, 12)
        forest = rect(2, 2, 18, 18)
        bog = rect(12, 2, 18, 6)
        road = line((0, 15), (20, 15))
        mesh = grid_mesh(STEPS, STEPS, [forest, lake, bog, road])
        return mesh, [(forest, 312), (lake, 512), (bog, 412)]

    def test_the_same_call_twice_is_the_same_array(self, case: Any) -> None:
        mesh, polygons = case
        first, second = label(mesh, polygons), label(mesh, polygons)
        assert np.asarray(first.codes).tobytes() == np.asarray(second.codes).tobytes()

    def test_permuted_triangles_edges_and_polygons_give_the_same_labels(self, case: Any) -> None:
        mesh, polygons = case
        base = np.asarray(label(mesh, polygons).codes)
        rng = np.random.default_rng(29)
        order = rng.permutation(len(mesh.triangles))
        rows = rng.permutation(len(mesh.edges))
        shuffled = Mesh(mesh.vertices, mesh.triangles[order], mesh.edges[rows][:, ::-1])
        again = label(shuffled, polygons[::-1])
        assert_array_equal(np.asarray(again.codes), base[order])
        first = label(mesh, polygons)
        assert (again.regions, again.outside, again.overlapped, again.thin) == (
            first.regions,
            first.outside,
            first.overlapped,
            first.thin,
        )


# ------------------------------------------- increment 30a: the fast helpers
#
# `docs/increments/30a-landcover-speed.md`, section 7. The output does not
# change, so the red tests are the contracts of the two new helpers (R1-R3)
# and the pins (P1-P4) guard what the rewrite could break; the pins are green
# before the rewrite, which is their point.


def best(ids: np.ndarray, r: np.ndarray) -> np.ndarray:
    return np.asarray(landcover._best(ids, r))


def member(values: np.ndarray, table: np.ndarray) -> np.ndarray:
    return np.asarray(landcover._member(values, table))


def ids_of(groups: np.ndarray) -> np.ndarray:
    """`ids` as `regions` builds them: each triangle's group, named by the
    smallest triangle index in that group."""
    _, first, inverse = np.unique(groups, return_index=True, return_inverse=True)
    return first[inverse].astype(np.int64)


def best_by_sort(ids: np.ndarray, r: np.ndarray) -> np.ndarray:
    """The oracle: per id in increasing order, the largest `r`, ties to the
    lowest triangle index (the rule `np.lexsort` gave at `a7154ec`)."""
    out = []
    for root in np.unique(ids):
        members = np.flatnonzero(ids == root)
        top = r[members].max()
        out.append(members[r[members] == top].min())
    return np.array(out, dtype=np.int64)


class TestBest:
    """R1, R2: per component, the triangle of largest `r`, ties to the lowest
    index, one per component in increasing id."""

    def test_ties_go_to_the_lowest_index_whatever_the_listing(self) -> None:
        """Interleaved groups. Group 0 = {0, 3, 5, 8}: the maximum 5.0 is
        tied at 3 and 8, not adjacent and not first. Group 1 = {1, 4, 7}: every
        `r` is 0 (zero-area triangles). Groups 2 and 6 are single triangles.
        The answer for group 0 (3) is above group 1's (1): the order is by id,
        not by index."""
        ids = np.array([0, 1, 2, 0, 1, 0, 6, 1, 0], dtype=np.int64)
        r = np.array([1.0, 0.0, 0.25, 5.0, 0.0, 2.0, 0.5, 0.0, 5.0])
        chosen = np.asarray(best(ids, r))
        assert chosen.dtype == np.int64
        assert_array_equal(chosen, [3, 1, 2, 6])
        assert_array_equal(chosen, best_by_sort(ids, r))

    def test_a_tie_that_includes_the_root_goes_to_the_root(self) -> None:
        ids = np.array([0, 1, 0, 1, 0], dtype=np.int64)
        r = np.array([7.0, 1.0, 7.0, 3.0, 2.0])
        assert_array_equal(best(ids, r), [0, 3])

    def test_every_triangle_its_own_component(self) -> None:
        ids = np.arange(5, dtype=np.int64)
        r = np.array([0.0, 3.0, 0.0, 1.0, 2.0])
        assert_array_equal(best(ids, r), np.arange(5))

    def test_one_component_of_zero_area_triangles(self) -> None:
        ids = np.zeros(4, dtype=np.int64)
        assert_array_equal(best(ids, np.zeros(4)), [0])

    @pytest.mark.parametrize("seed", [30, 31, 32])
    def test_agrees_with_a_sort_on_a_random_case(self, seed: int) -> None:
        """3 000 triangles in 40 groups, `r` from five values: ties at a
        group's maximum are common, and the check below insists on it."""
        rng = np.random.default_rng(seed)
        ids = ids_of(rng.integers(0, 40, size=3000))
        r = rng.choice([0.0, 0.5, 1.0, 2.0, 3.0], size=3000)
        expected = best_by_sort(ids, r)
        tied = sum(
            int((r[ids == root] == r[ids == root].max()).sum() > 1) for root in np.unique(ids)
        )
        assert len(expected) == 40 and tied > 20
        chosen = np.asarray(best(ids, r))
        assert chosen.dtype == np.int64
        assert_array_equal(chosen, expected)


class TestMember:
    """R3: membership in a sorted, unique int64 table."""

    TABLE = np.array([10, 20, 30], dtype=np.int64)

    def test_values_around_the_table(self) -> None:
        values = np.array([5, 10, 15, 20, 25, 30, 35, 30, 5], dtype=np.int64)
        found = np.asarray(member(values, self.TABLE))
        assert found.dtype == np.bool_
        assert_array_equal(found, [False, True, False, True, False, True, False, True, False])

    def test_an_empty_table_finds_nothing(self) -> None:
        values = np.array([-1, 0, 7], dtype=np.int64)
        found = np.asarray(member(values, np.zeros(0, dtype=np.int64)))
        assert found.dtype == np.bool_
        assert_array_equal(found, [False, False, False])

    def test_no_values(self) -> None:
        found = np.asarray(member(np.zeros(0, dtype=np.int64), self.TABLE))
        assert found.shape == (0,)

    def test_a_table_of_one_entry(self) -> None:
        values = np.array([19, 20, 21], dtype=np.int64)
        assert_array_equal(member(values, np.array([20], dtype=np.int64)), [False, True, False])

    def test_keys_near_two_to_the_62(self) -> None:
        """Edge keys are `min * N + max`, about `N**2`: exact at 2**62, where a
        float64 round trip would merge neighbours."""
        big = 2**62
        table = np.array([big - 3, big, big + 5], dtype=np.int64)
        values = np.array([big - 4, big - 3, big - 1, big, big + 1, big + 5, big + 6])
        assert_array_equal(
            member(values.astype(np.int64), table),
            [False, True, False, True, False, True, False],
        )

    def test_agrees_with_isin_on_a_random_case(self) -> None:
        rng = np.random.default_rng(30)
        table = np.unique(rng.integers(0, 2_000, size=300)).astype(np.int64)
        values = rng.integers(-50, 2_050, size=5_000).astype(np.int64)
        assert_array_equal(member(values, table), np.isin(values, table))


class TestBoundaryPoints:
    """P1: the point-in-polygon test counts the boundary as inside, through
    `_lookup`. Every edge is axis-parallel and every point's coordinate is the
    same double as the edge's, so "on the edge" is exact at UTM magnitudes.

    `west` (100 m²) carries the larger code and `east` (200 m²) the smaller,
    so the area rule, not the code rule, decides the shared edge."""

    @pytest.fixture(params=[False, True], ids=["listed", "reversed"])
    def polygons(self, request: pytest.FixtureRequest) -> list[tuple[BaseGeometry, int]]:
        west, east, north = rect(0, 0, 10, 10), rect(10, 0, 30, 10), rect(0, 10, 10, 25)
        holed = Polygon(rect(40, 0, 60, 20).exterior, [rect(45, 5, 55, 15).exterior])
        out: list[tuple[BaseGeometry, int]] = [(west, 512), (east, 311), (north, 412), (holed, 324)]
        return out[::-1] if request.param else out

    @staticmethod
    def at(*points: tuple[float, float]) -> np.ndarray:
        return np.array([(X0 + x, Y0 + y) for x, y in points], dtype=np.float64)

    def test_shared_edge_counts_in_both_and_goes_to_the_smaller(self, polygons: Any) -> None:
        codes, hits = landcover._lookup(self.at((10, 5)), polygons)
        assert_array_equal(hits, [2])
        assert_array_equal(codes, [512])

    def test_corner_shared_by_three(self, polygons: Any) -> None:
        codes, hits = landcover._lookup(self.at((10, 10)), polygons)
        assert_array_equal(hits, [3])
        assert_array_equal(codes, [512])

    def test_hole_boundary_and_inside_the_hole(self, polygons: Any) -> None:
        """On the hole's side and at its corner: in the holed polygon. Inside
        the hole: in none."""
        codes, hits = landcover._lookup(self.at((45, 10), (45, 5), (50, 10)), polygons)
        assert_array_equal(hits, [1, 1, 0])
        assert_array_equal(codes, [324, 324, 0])


class TestCallerPolygons:
    """P2: `label_triangles` leaves the caller's polygons as they came:
    not prepared, the same WKB."""

    def test_unprepared_and_unchanged(self) -> None:
        west, east = rect(2, 4, 10, 16), rect(10, 4, 18, 16)
        pair = shapely.MultiPolygon([rect(0, 0, 2, 2), rect(18, 18, 20, 20)])
        polygons: list[tuple[BaseGeometry, int]] = [(west, 311), (east, 512), (pair, 412)]
        before = [shapely.to_wkb(p) for p, _ in polygons]
        mesh = grid_mesh(STEPS, STEPS, [west, east, pair])
        label(mesh, polygons)
        assert [bool(shapely.is_prepared(p)) for p, _ in polygons] == [False] * 3
        assert [shapely.to_wkb(p) for p, _ in polygons] == before


class TestDegenerateMeshes:
    def test_a_zero_area_component_is_labelled_and_counted_thin(self) -> None:
        """P3: three collinear vertices, all three sides constraints, beside an
        ordinary triangle in another polygon. Its `r` is 0 and its incentre is
        its middle vertex, inside the square coded 311."""
        vertices = np.array(
            [
                [X0 + 1, Y0 + 1, 0.0],
                [X0 + 2, Y0 + 1, 0.0],
                [X0 + 3, Y0 + 1, 0.0],
                [X0 + 20, Y0, 0.0],
                [X0 + 30, Y0, 0.0],
                [X0 + 20, Y0 + 10, 0.0],
            ]
        )
        triangles = np.array([[0, 1, 2], [3, 4, 5]])
        edges = np.array([[0, 1], [1, 2], [2, 0], [3, 4], [4, 5], [5, 3]])
        polygons = [(rect(0, 0, 10, 10), 311), (rect(15, -5, 35, 15), 512)]
        labels = landcover.label_triangles(
            vertices, triangles, edges, polygons=polygons, margin=MARGIN
        )
        assert_array_equal(labels.codes, [311, 512])
        assert (labels.regions, labels.outside, labels.overlapped, labels.thin) == (2, 0, 0, 1)

    @pytest.mark.parametrize("vertex_count", [0, 4])
    def test_no_triangles(self, vertex_count: int) -> None:
        """P4: empty codes (int32) and four zero counts; `regions` empty."""
        vertices = np.column_stack(
            [np.arange(vertex_count) + X0, np.full(vertex_count, Y0), np.zeros(vertex_count)]
        )
        triangles = np.zeros((0, 3), dtype=np.int64)
        edges = np.zeros((0, 2), dtype=np.int64)
        labels = landcover.label_triangles(
            vertices.reshape(-1, 3),
            triangles,
            edges,
            polygons=[(rect(0, 0, 10, 10), 311)],
            margin=MARGIN,
        )
        codes = np.asarray(labels.codes)
        assert codes.dtype == np.int32 and codes.shape == (0,)
        assert (labels.regions, labels.outside, labels.overlapped, labels.thin) == (0, 0, 0, 0)
        ids = np.asarray(landcover.regions(triangles, edges))
        assert ids.shape == (0,)
