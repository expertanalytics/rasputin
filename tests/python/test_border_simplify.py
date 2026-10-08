"""Land-cover borders simplified within a band: the binding and the adapter (increment 32).

``docs/increments/32-landcover-simplify.md``, sections 6 ("The Python
interface", "The C++ interface"), 7 (the guarantees) and 9, tests 12 and 13.
The invariant-critical suite is the C++ one,
``tests/cpp/unit/test_border_collapse.cpp``; this file checks that the
binding and the adapter carry its guarantees to shapely polygons, with
shapely as the oracle.

Interface, as section 6 fixes it, and PINNED HERE where it leaves the
Python spelling open (listed in the handback):

- ``_core.simplify_borders(points, ring_starts, band) -> BorderOutcome``:
  ``points`` an ``(N, 2)`` float64 array-like, ``ring_starts`` a 1-D array-like
  of ``R + 1`` unsigned integers; any other shape is a ``ValueError`` naming
  the shape. Releases the GIL.
- ``BorderOutcome.points`` an ``(M, 2)`` float64 array, ``.ring_starts`` a 1-D
  uint64 array, ``.status`` a ``_core.BorderStatus`` (``Ok``, ``InvalidBand``,
  ``BadRings``), ``.counts`` a ``_core.BorderCounts`` with ``junctions``,
  ``borders``, ``fixed_borders``, ``collinear``, ``collapses``,
  ``rejected_crossing`` and ``rejected_side``.
- ``tin_engine.border_simplify.simplify_borders(polygons, band_m) ->
  BorderResult``; ``BorderResult.polygons`` a tuple, same count and order as
  the input, each of the input's geometry type with its parts in the input's
  order and each part's holes in the input's order; ``.counts`` the binding's
  ``BorderCounts``. No polygons in, none out. Band 0: every polygon
  ``equals_exact`` its input, tolerance 0.
- A refusal is a ``ValueError`` whose message says "band" for
  ``InvalidBand`` and names a ring, a coordinate or finiteness for
  ``BadRings``.

Scale of the bounds: the hand-made coverage spans 1000 x 600 m at the origin
(areas to 1e-9 relative, the band to 1e-6 m); the CORINE extract is in
EPSG:3035 near (4.8e6, 5.4e6) m, where one ulp is 9.3e-10 m, and a part's
area is checked to 1e-9 relative or 1e-6 m2, whichever is larger (a few
collapses' rounding on its smallest parts).

Increment 32's fix (section 15.2, the clearance; 15.3, test 21), PINNED
HERE where it leaves the Python spelling open (listed in the handback):

- ``_core.simplify_borders(points, ring_starts, band, clearance=0.0)``, the
  clearance a keyword; ``_core.BorderStatus.InvalidClearance`` for a negative
  or non-finite one; ``BorderCounts`` gains ``rejected_clearance`` and
  ``skipped_placements``.
- ``border_simplify.simplify_borders(polygons, band_m, clearance_m=0.0)``; a
  refused clearance is a ``ValueError`` whose message says "clearance";
  ``clearance_m`` 0 gives the call without it, ``equals_exact`` tolerance 0.
- Test 21 also checks every output edge with a new vertex as an end against
  every output vertex it does not end (the other half of 15.2's added
  guarantee), with the same slack.

Scale of test 21's bound: 1 m less 1e-6 m. At EPSG:3035 coordinates near
(4.8e6, 5.4e6) m one ulp is 9.3e-10 m, so 1e-9 would leave about one ulp for
two distance computations (the kernel's and GEOS's) that may round
differently (15.3).

RED at the commit that added test 21 and the clearance tests (``2fba4349``):
``_core.simplify_borders`` took no clearance, so every call that passed one
was a ``TypeError``; the adapter had no ``clearance_m``; there was
no ``InvalidClearance`` and no ``rejected_clearance`` or
``skipped_placements``.

RED at the commit that added this file (``44f25968``): ``_core`` had no
``simplify_borders`` and ``tin_engine.border_simplify`` did not exist; every
test failed on the ``core`` or ``adapter`` fixture's assertion naming what was
missing.
"""

from __future__ import annotations

import importlib
from collections import defaultdict
from typing import Any

import numpy as np
import pytest
import shapely
from shapely.geometry import MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry

from gil_probe import ticks_during
from gpkg_fixtures import EXTRACT

COUNTS = (
    "junctions",
    "borders",
    "fixed_borders",
    "collinear",
    "collapses",
    "rejected_crossing",
    "rejected_side",
    "rejected_clearance",  # section 15.2
    "skipped_placements",
)
RELEASED_TICKS = 20  # as test_core_noding.py: a released call lets the ticker run


@pytest.fixture(scope="module")
def core() -> Any:
    import tin_engine._core as core

    assert hasattr(core, "simplify_borders"), "tin_engine._core has no simplify_borders"
    return core


@pytest.fixture(scope="module")
def adapter() -> Any:
    try:
        return importlib.import_module("tin_engine.border_simplify")
    except ModuleNotFoundError as exc:
        pytest.fail(f"tin_engine.border_simplify does not exist ({exc})")


# ---------------------------------------------------------------- helpers


def zigzag(
    a: tuple[float, float], b: tuple[float, float], n: int, amp: float
) -> list[tuple[float, float]]:
    """Inner vertices of a zig-zag from ``a`` to ``b`` in ``n`` steps, +-``amp`` across."""
    (ax, ay), (bx, by) = a, b
    length = float(np.hypot(bx - ax, by - ay))
    cx, cy = -(by - ay) / length, (bx - ax) / length
    out = []
    for k in range(1, n):
        s = amp if k % 2 else -amp
        out.append((ax + k / n * (bx - ax) + s * cx, ay + k / n * (by - ay) + s * cy))
    return out


def jagged_circle(cx: float, cy: float, r: float, n: int) -> list[tuple[float, float]]:
    """``n`` vertices, counter-clockwise, radii alternating ``r`` and ``0.9 r``."""
    t = 2 * np.pi * np.arange(n) / n
    rho = np.where(np.arange(n) % 2 == 0, r, 0.9 * r)
    return list(zip((cx + rho * np.cos(t)).tolist(), (cy + rho * np.sin(t)).tolist(), strict=True))


def flat(rings: list[list[tuple[float, float]]]) -> tuple[np.ndarray, np.ndarray]:
    points = np.array([p for r in rings for p in r], dtype=np.float64)
    starts = np.cumsum([0, *(len(r) for r in rings)]).astype(np.uint64)
    return points, starts


def rings_of(polygons: list[BaseGeometry]) -> list[list[tuple[float, float]]]:
    """Every ring of every part, open, in polygon, part, ring order."""
    return [
        [(x + 0.0, y + 0.0) for x, y in shapely.get_coordinates(r)[:-1].tolist()]
        for g in polygons
        for part in shapely.get_parts(g)
        for r in shapely.get_rings(part)
    ]


def junctions(polygons: list[BaseGeometry]) -> set[tuple[float, float]]:
    """Section 6's rule: a vertex with other than two distinct incident edges,
    or two whose sets of rings differ. Edges by exact coordinates."""
    users: dict[tuple[Any, Any], list[int]] = defaultdict(list)
    for k, ring in enumerate(rings_of(polygons)):
        for i, p in enumerate(ring):
            q = ring[(i + 1) % len(ring)]
            users[(min(p, q), max(p, q))].append(k)
    incident: dict[tuple[float, float], list[tuple[Any, Any]]] = defaultdict(list)
    for e in users:
        incident[e[0]].append(e)
        incident[e[1]].append(e)
    return {
        v
        for v, es in incident.items()
        if len(es) != 2 or sorted(users[es[0]]) != sorted(users[es[1]])
    }


def sharing_pairs(polygons: list[BaseGeometry]) -> set[tuple[int, int]]:
    """The pairs of polygons whose boundaries share a stretch of positive length."""
    return {
        (i, j)
        for i in range(len(polygons))
        for j in range(i + 1, len(polygons))
        if shapely.length(shapely.intersection(polygons[i].boundary, polygons[j].boundary)) > 0
    }


def shape_of(g: BaseGeometry) -> tuple[str, list[int]]:
    """Geometry type, and the number of holes of each part in order."""
    return g.geom_type, [len(p.interiors) for p in shapely.get_parts(g)]


def check_coverage(
    before: list[BaseGeometry],
    after: tuple[BaseGeometry, ...],
    band: float,
    area_abs: float = 0.0,
    densify: bool = True,
) -> None:
    """Test 12's checks. ``densify`` False measures the band on prepared
    geometries, sampled every 2 m and at every vertex, both ways (the
    CORINE extract is too large for GEOS's densified Hausdorff)."""
    assert len(after) == len(before)
    assert all(g.is_valid for g in after)
    assert bool(shapely.coverage_is_valid(np.array(after, dtype=object)))
    for a, b in zip(before, after, strict=True):
        assert shape_of(b) == shape_of(a)
        for pa, pb in zip(shapely.get_parts(a), shapely.get_parts(b), strict=True):
            assert abs(pb.area - pa.area) <= max(1e-9 * pa.area, area_abs), (pa.area, pb.area)
        if densify:
            # Every segment cut into 1 000: a lower estimate of the distance.
            assert shapely.hausdorff_distance(a.boundary, b.boundary, densify=0.001) <= band + 1e-6
        else:
            for src, dst in ((a.boundary, b.boundary), (b.boundary, a.boundary)):
                pts = shapely.points(shapely.get_coordinates(shapely.segmentize(src, 2.0)))
                shapely.prepare(dst)
                assert bool(shapely.dwithin(pts, dst, band + 1e-6).all())
    assert sharing_pairs(list(after)) == sharing_pairs(before)
    out = {
        (x + 0.0, y + 0.0)
        for x, y in shapely.get_coordinates(np.array(after, dtype=object)).tolist()
    }
    assert junctions(before) <= out


# ---------------------------------------------------------------- fixtures


def hand_made() -> list[BaseGeometry]:
    """Test 12's coverage of 1000 x 600 m, four classes:

    - A, a MultiPolygon: west of a zig-zag border near x = 400, and an
      island inside C;
    - B, south-east, below a zig-zag near y = 300, with a hole holding D;
    - C, north-east, with a hole holding A's island;
    - D, the island in B.

    The two zig-zags (8 m either side, 20 m steps) meet B and C's border at
    the junction (400, 300); the borders meet the outline at (400, 0),
    (400, 600) and (1000, 300). The islands are jagged circles of 32
    vertices."""
    j, s, n, e = (400.0, 300.0), (400.0, 0.0), (400.0, 600.0), (1000.0, 300.0)
    down, up, east = zigzag(j, s, 15, 8.0), zigzag(j, n, 15, 8.0), zigzag(j, e, 30, 8.0)
    west = [(0.0, 0.0), s, *down[::-1], j, *up, n, (0.0, 600.0)]
    a_island = jagged_circle(700.0, 450.0, 60.0, 32)
    d_island = jagged_circle(700.0, 150.0, 50.0, 32)
    south_east = [s, (1000.0, 0.0), e, *east[::-1], j, *down]
    north_east = [j, *east, e, (1000.0, 600.0), n, *up[::-1]]
    polygons = [
        MultiPolygon([Polygon(west), Polygon(a_island)]),
        Polygon(south_east, [d_island[::-1]]),
        Polygon(north_east, [a_island[::-1]]),
        Polygon(d_island),
    ]
    return polygons


def corine_coverage() -> list[BaseGeometry]:
    """The committed CORINE extract (EPSG:3035), repaired and merged by class
    as the land-cover stage's steps 2 and 3 do: 1 m, ``min_area``."""
    from tin_engine.feature_input import read_source

    rows = read_source(EXTRACT, None, "Code_18", lambda _crs: (-1e12, -1e12, 1e12, 1e12)).rows
    polys = np.array([g for _, g, _ in rows], dtype=object)
    codes = [str(v) for *_, v in rows]
    clean = shapely.coverage_clean(
        polys, snapping_distance=1.0, gap_width=1.0, merge_strategy="min_area"
    )
    groups: dict[str, list[int]] = defaultdict(list)
    for i, c in enumerate(codes):
        groups[c].append(i)
    merged = [shapely.coverage_union_all(clean[idx]) for idx in groups.values()]
    return [g for g in merged if not g.is_empty]


# ---------------------------------------------------------------- the binding


class TestBinding:
    def test_two_squares_through_the_binding(self, core: Any) -> None:
        up = zigzag((100.0, 0.0), (100.0, 100.0), 50, 3.0)
        left = [(0.0, 0.0), (100.0, 0.0), *up, (100.0, 100.0), (0.0, 100.0)]
        right = [(100.0, 0.0), (200.0, 0.0), (200.0, 100.0), (100.0, 100.0), *up[::-1]]
        points, starts = flat([left, right])
        out = core.simplify_borders(points, starts, 10.0)
        assert out.status == core.BorderStatus.Ok
        got = np.asarray(out.points)
        assert got.dtype == np.float64 and got.ndim == 2 and got.shape[1] == 2
        rs = np.asarray(out.ring_starts)
        assert rs.dtype == np.uint64 and rs.ndim == 1 and len(rs) == 3
        assert rs[0] == 0 and rs[-1] == len(got)
        assert len(got) < len(points)
        polys = [Polygon(got[rs[k] : rs[k + 1]]) for k in range(2)]
        for p, src in zip(polys, (left, right), strict=True):
            assert abs(p.area - Polygon(src).area) <= 1e-12 * Polygon(src).area
        assert isinstance(out.counts, core.BorderCounts)
        assert out.counts.collapses > 0
        for name in COUNTS:
            assert isinstance(getattr(out.counts, name), int), name
        again = core.simplify_borders(points, starts, 10.0)
        assert np.asarray(again.points).tobytes() == got.tobytes()

    def test_band_0_is_the_input(self, core: Any) -> None:
        points, starts = flat([[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (-0.0, 1.0)]])
        out = core.simplify_borders(points, starts, 0.0)
        assert out.status == core.BorderStatus.Ok
        assert np.asarray(out.points).tobytes() == points.tobytes()
        assert np.asarray(out.ring_starts).tolist() == starts.tolist()

    def test_the_statuses_cross_the_binding(self, core: Any) -> None:
        points, starts = flat([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]])
        status = core.BorderStatus
        # Section 15.2 adds InvalidClearance (this line said three, test 16).
        assert {s.name for s in status.__members__.values()} == {
            "Ok",
            "InvalidBand",
            "BadRings",
            "InvalidClearance",
        }
        for band in (-1.0, float("nan"), float("inf")):
            out = core.simplify_borders(points, starts, band)
            assert out.status == status.InvalidBand
            assert len(np.asarray(out.points)) == 0 and len(np.asarray(out.ring_starts)) == 0
        bad = points.copy()
        bad[1, 0] = float("nan")
        assert core.simplify_borders(bad, starts, 1.0).status == status.BadRings
        assert (
            core.simplify_borders(points, np.array([0, 3], dtype=np.uint64), 1.0).status
            == status.BadRings
        )
        assert (
            core.simplify_borders(points, np.array([0, 2, 4], dtype=np.uint64), 1.0).status
            == status.BadRings
        )

    @pytest.mark.parametrize("shape", [(4,), (4, 3), (2, 4, 2)])
    def test_points_of_another_shape_are_a_value_error(
        self, core: Any, shape: tuple[int, ...]
    ) -> None:
        with pytest.raises(ValueError, match=r"(?i)shape|\(N, ?2\)"):
            core.simplify_borders(np.zeros(shape), np.array([0, 4], dtype=np.uint64), 1.0)

    def test_starts_of_another_shape_are_a_value_error(self, core: Any) -> None:
        points, _ = flat([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]])
        with pytest.raises(ValueError, match=r"(?i)shape|1-D|one-dimensional"):
            core.simplify_borders(points, np.array([[0, 4]], dtype=np.uint64), 1.0)

    def test_it_releases_the_gil(self, core: Any) -> None:
        """A 40 x 40 grid of 100 m cells whose inner sides are zig-zags of 30
        steps: about 200 000 ring points."""
        cells, size, steps = 40, 100.0, 30

        def side(
            a: tuple[float, float], b: tuple[float, float], inner: bool
        ) -> list[tuple[float, float]]:
            return zigzag(a, b, steps, 2.0) if inner else []

        h = {
            (i, j): side((i * size, j * size), ((i + 1) * size, j * size), 0 < j < cells)
            for i in range(cells)
            for j in range(cells + 1)
        }
        v = {
            (i, j): side((i * size, j * size), (i * size, (j + 1) * size), 0 < i < cells)
            for i in range(cells + 1)
            for j in range(cells)
        }
        rings = []
        for i in range(cells):
            for j in range(cells):
                x0, y0, x1, y1 = i * size, j * size, (i + 1) * size, (j + 1) * size
                rings.append(
                    [
                        (x0, y0),
                        *h[i, j],
                        (x1, y0),
                        *v[i + 1, j],
                        (x1, y1),
                        *h[i, j + 1][::-1],
                        (x0, y1),
                        *v[i, j][::-1],
                    ]
                )
        points, starts = flat(rings)
        out, ticks, elapsed = ticks_during(lambda: core.simplify_borders(points, starts, 10.0))
        assert out.status == core.BorderStatus.Ok
        assert elapsed >= 0.05, f"took only {elapsed:.3f}s: too fast to see; enlarge the grid"
        assert ticks >= RELEASED_TICKS, f"only {ticks} ticks in {elapsed:.3f}s: it holds the GIL"


# ---------------------------------------------------------------- the adapter


class TestAdapter:
    BAND = 12.0

    def test_12_the_premise_a_valid_coverage(self) -> None:
        polygons = hand_made()
        assert all(p.is_valid for p in polygons)
        assert bool(shapely.coverage_is_valid(np.array(polygons, dtype=object)))
        assert (400.0, 300.0) in junctions(polygons)

    def test_12_a_hand_made_coverage(self, adapter: Any) -> None:
        polygons = hand_made()
        result = adapter.simplify_borders(polygons, self.BAND)
        assert isinstance(result.polygons, tuple)
        check_coverage(polygons, result.polygons, self.BAND)
        before = int(shapely.get_num_coordinates(np.array(polygons, dtype=object)).sum())
        after = int(shapely.get_num_coordinates(np.array(result.polygons, dtype=object)).sum())
        assert after < before
        assert result.counts.collapses > 0

    def test_12_the_corine_extract(self, adapter: Any) -> None:
        polygons = corine_coverage()
        result = adapter.simplify_borders(polygons, 50.0)
        check_coverage(polygons, result.polygons, 50.0, area_abs=1e-6, densify=False)
        assert result.counts.collapses > 0

    def test_band_0_gives_the_input(self, adapter: Any) -> None:
        polygons = hand_made()
        result = adapter.simplify_borders(polygons, 0.0)
        for a, b in zip(polygons, result.polygons, strict=True):
            assert shapely.equals_exact(a, b, tolerance=0.0)

    def test_no_polygons(self, adapter: Any) -> None:
        assert adapter.simplify_borders([], 50.0).polygons == ()

    def test_the_result_is_frozen(self, adapter: Any) -> None:
        result = adapter.simplify_borders(hand_made(), self.BAND)
        with pytest.raises((AttributeError, TypeError)):
            result.polygons = ()

    @pytest.mark.parametrize("band", [-1.0, float("nan"), float("inf")])
    def test_13_an_invalid_band_is_a_value_error(self, adapter: Any, band: float) -> None:
        with pytest.raises(ValueError, match=r"(?i)band"):
            adapter.simplify_borders(hand_made(), band)

    def test_13_a_non_finite_coordinate_is_a_value_error(self, adapter: Any) -> None:
        bad = Polygon([(0.0, 0.0), (10.0, 0.0), (float("nan"), 10.0), (0.0, 10.0)])
        with pytest.raises(ValueError, match=r"(?i)ring|coordinate|finite"):
            adapter.simplify_borders([bad], 1.0)


# ---------------------------------------------------------------- the clearance (section 15)


def new_vertex_clearance(
    before: list[BaseGeometry], after: tuple[BaseGeometry, ...]
) -> tuple[float, float, int]:
    """Test 21's distances, from input and output alone: the smallest
    distance from an output vertex the input did not have to an output edge
    it does not end, and from an output edge with such a vertex as an end to
    an output vertex it does not end; and the number of new vertices. Ends
    by exact coordinates, -0.0 and 0.0 alike."""
    old = {p for r in rings_of(before) for p in r}
    edges: set[tuple[tuple[float, float], tuple[float, float]]] = set()
    for ring in rings_of(list(after)):
        for i, p in enumerate(ring):
            q = ring[(i + 1) % len(ring)]
            edges.add((min(p, q), max(p, q)))
    vertices = sorted({v for e in edges for v in e})
    new = [v for v in vertices if v not in old]
    lines = list(edges)
    tree = shapely.STRtree([shapely.LineString(e) for e in lines])
    points = shapely.points(np.array(vertices, dtype=np.float64))
    vertex_side = edge_side = float("inf")
    # Candidates within 2 m only: test 21's floor is 1 m.
    hits = tree.query(points, predicate="dwithin", distance=2.0)
    new_set = set(new)
    for vi, ei in hits.T.tolist():
        v, (a, b) = vertices[vi], lines[ei]
        if v in (a, b):
            continue
        d = float(shapely.distance(points[vi], tree.geometries[ei]))
        if v in new_set:
            vertex_side = min(vertex_side, d)
        if a in new_set or b in new_set:
            edge_side = min(edge_side, d)
    return vertex_side, edge_side, len(new)


class TestClearance:
    """Section 15.2's clearance through the binding and the adapter."""

    def test_the_binding_takes_a_clearance(self, core: Any) -> None:
        points, starts = flat([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]])
        out = core.simplify_borders(points, starts, 1.0, clearance=1.0)
        assert out.status == core.BorderStatus.Ok
        assert out.counts.rejected_clearance == 0
        assert out.counts.skipped_placements == 0

    @pytest.mark.parametrize("clearance", [-1.0, float("nan"), float("inf")])
    def test_20_an_invalid_clearance_crosses_the_binding(self, core: Any, clearance: float) -> None:
        points, starts = flat([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]])
        out = core.simplify_borders(points, starts, 1.0, clearance=clearance)
        assert out.status == core.BorderStatus.InvalidClearance
        assert len(np.asarray(out.points)) == 0 and len(np.asarray(out.ring_starts)) == 0

    @pytest.mark.parametrize("clearance", [-1.0, float("nan"), float("inf")])
    def test_20_an_invalid_clearance_is_a_value_error(self, adapter: Any, clearance: float) -> None:
        with pytest.raises(ValueError, match=r"(?i)clearance"):
            adapter.simplify_borders(hand_made(), 12.0, clearance_m=clearance)

    def test_20_clearance_0_is_the_call_without_it(self, adapter: Any) -> None:
        polygons = hand_made()
        without = adapter.simplify_borders(polygons, 12.0)
        at_0 = adapter.simplify_borders(polygons, 12.0, clearance_m=0.0)
        for a, b in zip(without.polygons, at_0.polygons, strict=True):
            assert shapely.equals_exact(a, b, tolerance=0.0)
        assert at_0.counts.collapses == without.counts.collapses

    @pytest.mark.parametrize(
        ("coverage", "band", "area_abs", "densify"),
        [(hand_made, TestAdapter.BAND, 0.0, True), (corine_coverage, 50.0, 1e-6, False)],
        ids=["hand_made", "corine"],
    )
    def test_21_new_vertices_keep_the_clearance(
        self, adapter: Any, coverage: Any, band: float, area_abs: float, densify: bool
    ) -> None:
        polygons = coverage()
        result = adapter.simplify_borders(polygons, band, clearance_m=1.0)
        check_coverage(polygons, result.polygons, band, area_abs=area_abs, densify=densify)
        vertex_side, edge_side, made = new_vertex_clearance(polygons, result.polygons)
        assert made > 0  # the premise: the simplifier placed vertices
        assert vertex_side >= 1.0 - 1e-6, vertex_side
        assert edge_side >= 1.0 - 1e-6, edge_side
