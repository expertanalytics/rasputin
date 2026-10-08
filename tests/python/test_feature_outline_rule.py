"""Ola's outline rule: land-cover borders within D of the outline go onto it (20c-3).

``docs/increments/20c-soft-quality.md``, M6 (c), "Steps 2 to 4" and "Tests
``@tester`` writes red first", 20c-3: OR1 to OR8; OR9 and the T2 test are
added by "Rulings on 20c-3's timed check and gate runs". Invariant-critical
(mutation testing required: the rule rewrites input borders); the mutants
the design names are run after the green step, on the built code:

- the cut not in a fixed order: OR2;
- the vertex dropped at a rounded point: OR4
  (``test_or4_a_rounded_point_next_to_a_vertex_goes_onto_it``, added by the
  mutation round; ``TestTheCorner`` does not catch it, as its corner lies at
  a multiple of D); the outline's vertices dropped from a join:
  ``test_or5_a_shallow_notch_is_followed`` (added by the mutation round);
- the stretches on the outline kept: OR1, OR3;
- the rounding removed: OR4 (``TestRounding``);
- the inlet test removed: OR5 (``test_or5_an_inlet_is_joined_straight``).
- the cuts counted from the edge's own end, not along its line: OR8
  (``TestTheReadRegion``, added by ruling G4);
- ``_cut``'s end test removed (a multiple of D within D/2 of the edge's
  end kept): OR9 (``TestACutBesideAVertex``, added by ruling T1). The T2
  test (``TestOnlyWhatMovedIsRebuilt``) and the P1 test
  (``TestARingWithNoKeptVertex``, added by ruling P1 of "Rulings on the fix
  round") are not mutation-critical.

The rule's own function is called directly, on polygons and an outline in
EPSG:25833 at ``feature_fixtures``' UTM-shaped offset. The outline is a
``box(0, 0, 200, 200)`` unless a test says otherwise; the border under test
lies along its south side, 50 m or more from every other side.

Scale of the bounds: D = 5 m (ruling 6's default) on fixtures of 100 to
250 m, coordinates near 5e5 and 6.6e6 m. "On the outline" is within 1e-9 m
(a double there resolves about 1e-9 m); "near x" is within D / 2, the
rule's own rounding radius.

PINNED HERE, where the design leaves it open (listed for ``@architect``):

- ``feature_input.snap_to_outline(polygons, outline, distance)``:
  ``polygons`` a sequence of shapely ``Polygon`` (one coverage, in the
  computation CRS), ``outline`` the domain ``Polygon`` (every ring of it is
  outline), ``distance`` D in metres. It returns an object with ``lines``,
  one tuple of ``LineString`` per input polygon in input order (its
  linework after the rule, not yet clipped to the domain), and
  ``area_changed``, the land-cover area inside the outline that changed
  polygon, m² (the record's "land-cover area that changed class").
- D = 0 returns each polygon's rings, exterior first, as ``LineString``,
  coordinates identical (OR6).
- D at or above 100 m raises ``ValueError`` naming 100 (OR7); the CLI
  refusal is in ``test_cli_features_cleanup.py``.
- OR7's guard patches ``shapely.buffer`` and ``shapely.snap`` (the design:
  "no ``buffer``, no ``snap``"); both geometry methods go through them.

- T2: the polygon's parts are read back from ``polygons`` as shapely
  parts; the far part is the one holding its input's centroid, and
  ``equals_exact`` at 0 also pins its ring's start and orientation (it is
  the input part itself, not a union's rewrite of it).

RED at the commit that adds the P1 test (on ``7a45dabc``, the code as on
``8716f062``): ``test_p1_the_area_moved_is_a_alone`` (``area_changed``
1 000 000 m², not 20 m²).

RED at the commit that adds OR9 and the T2 test (on ``eb810a90``): OR9's
two cases (a 2.4 cm edge beside V) and ``test_t2_the_far_part_passes_through``
(the far part comes back through ``union_all``, its ring restarted).

RED at the commit that adds this file: ``feature_input`` has no
``snap_to_outline``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from itertools import pairwise
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import shapely
from shapely import affinity
from shapely.geometry import LineString, Point, Polygon
from shapely.geometry.base import BaseGeometry

import feature_fixtures as ff
from feature_fixtures import X0, Y0, at

D = 5.0  # ruling 6's default, metres
ON = 1e-9  # "on the outline", metres; see the docstring


@pytest.fixture(scope="module")
def fi() -> ModuleType:
    return ff.feature_input()


def poly(*points: tuple[float, float]) -> Polygon:
    return Polygon([at(x, y) for x, y in points])


def box(x0: float, y0: float, x1: float, y1: float) -> Polygon:
    return poly((x0, y0), (x1, y0), (x1, y1), (x0, y1))


OUTLINE = box(0, 0, 200, 200)


def rule(
    fi: ModuleType, polygons: Sequence[Polygon], outline: Polygon = OUTLINE, d: float = D
) -> Any:
    out = fi.snap_to_outline(list(polygons), outline, d)
    assert len(out.lines) == len(polygons)
    return out


def inside(lines: Sequence[LineString], outline: Polygon) -> BaseGeometry:
    """The linework clipped to the outline, as step 5 clips it."""
    return shapely.intersection(shapely.union_all(list(lines)), outline)


def segments(lines: Sequence[LineString]) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Every segment, its ends sorted, so direction does not count."""
    out = []
    for line in lines:
        xy = shapely.get_coordinates(line)
        for a, b in pairwise(xy):
            out.append(tuple(sorted((tuple(map(float, a)), tuple(map(float, b))))))
    return out  # type: ignore[return-value]


def on_outline(xy: np.ndarray, outline: Polygon) -> np.ndarray:
    return np.asarray(shapely.distance(outline.boundary, shapely.points(xy))) <= ON


def rel_x(xy: np.ndarray) -> np.ndarray:
    return xy[:, 0] - X0


# ---------------------------------------------------------------- OR1


class TestParallelBorder:
    """OR1: a border 3 m inside the south side, parallel to it for 100 m."""

    NEAR = poly((50, 60), (50, 3), (150, 3), (150, 60))
    FAR = poly((50, 60), (50, 6), (150, 6), (150, 60))

    def test_or1_it_goes_onto_the_outline(self, fi: ModuleType) -> None:
        (lines,) = rule(fi, [self.NEAR]).lines
        kept = inside(lines, OUTLINE)
        band = shapely.intersection(kept, box(0, 0, 200, D - 1e-6))
        xy = shapely.get_coordinates(band)
        # Nothing within D of the outline but where the two sides cross the band.
        assert len(xy) > 0
        assert np.all((np.abs(rel_x(xy) - 50) <= D / 2) | (np.abs(rel_x(xy) - 150) <= D / 2)), xy
        # Each side leaves the outline at a point on it.
        ends = np.array(
            [shapely.get_coordinates(p)[k] for p in shapely.get_parts(kept) for k in (0, -1)]
        )
        touching = ends[on_outline(ends, OUTLINE)]
        assert len(touching) == 2, touching
        for x in (50, 150):
            assert (np.abs(rel_x(touching) - x) <= D / 2).sum() == 1, (x, touching)
        # The parallel stretch is gone (OR1's mutant: the stretches on the outline kept).
        assert shapely.intersection(kept, box(55, 0.5, 145, 4)).is_empty

    def test_or1_six_metres_inside_is_untouched(self, fi: ModuleType) -> None:
        (lines,) = rule(fi, [self.FAR]).lines
        assert shapely.equals(shapely.union_all(list(lines)), self.FAR.exterior)
        got = {tuple(p) for line in lines for p in shapely.get_coordinates(line)}
        assert got == {tuple(p) for p in shapely.get_coordinates(self.FAR)}


# ---------------------------------------------------------------- OR2


class TestSharedBorder:
    """OR2: two polygons sharing a border that crosses the band slanted
    give the same linework along it from either side, with the edge given in
    either direction, on a straight and on a diagonal outline edge."""

    LEFT = poly((20, -50), (100, -50), (110, 60), (20, 60))
    RIGHT = poly((100, -50), (180, -50), (180, 60), (110, 60))
    SHARED = LineString([at(100, -50), at(110, 60)])

    @pytest.mark.parametrize("angle", [0.0, 30.0], ids=["straight", "diagonal"])
    @pytest.mark.parametrize("reverse", [False, True], ids=["same", "reversed"])
    def test_or2_one_border_from_either_side(
        self, fi: ModuleType, angle: float, reverse: bool
    ) -> None:
        def turn(g: BaseGeometry) -> Any:
            return affinity.rotate(g, angle, origin=(X0, Y0)) if angle else g

        left = Polygon(self.LEFT.exterior.coords[::-1]) if reverse else self.LEFT
        outline, shared = turn(OUTLINE), turn(self.SHARED)
        a, b = rule(fi, [turn(left), turn(self.RIGHT)], outline).lines

        def along(lines: Sequence[LineString]) -> list[Any]:
            mids = [Point((p[0] + q[0]) / 2, (p[1] + q[1]) / 2) for p, q in segments(lines)]
            near = np.asarray(shapely.distance(shared, mids)) <= D
            return sorted(s for s, n in zip(segments(lines), near, strict=True) if n)

        assert along(a), "no linework along the shared border"
        assert along(a) == along(b)
        # It reaches the outline (so the shared stretch is the one the rule moved).
        xy = np.array([p for s in along(a) for p in s])
        assert on_outline(xy, outline).any()


# ---------------------------------------------------------------- OR3


class TestThinStrip:
    """OR3: a strip 3 m wide between the outline and its neighbour. The
    outline has vertices at x = 50 and 150, where the strip's sides cross,
    so the rule's rounding puts those points exactly there."""

    OUTLINE = poly((0, 0), (50, 0), (150, 0), (200, 0), (200, 200), (0, 200))
    STRIP = box(50, -100, 150, 3)
    NEIGHBOUR = box(50, 3, 150, 100)

    def test_or3_the_strip_vanishes_and_its_neighbour_reaches_the_outline(
        self, fi: ModuleType
    ) -> None:
        out = rule(fi, [self.STRIP, self.NEIGHBOUR], self.OUTLINE)
        strip, neighbour = (inside(lines, self.OUTLINE) for lines in out.lines)
        assert shapely.length(strip) == pytest.approx(0.0, abs=1e-9)
        assert shapely.intersection(neighbour, box(55, 0.5, 145, 5)).is_empty
        for x in (50, 150):
            assert shapely.distance(neighbour, Point(at(x, 0))) <= ON
        changed = shapely.intersection(self.STRIP, self.OUTLINE).area  # 300 m²
        assert out.area_changed == pytest.approx(changed, abs=1e-6)


# ---------------------------------------------------------------- OR4


class TestRounding:
    """OR4: points the rule puts on the outline are an outline vertex or more
    than D/2 from every other such point. Seven borders cross the south side
    at right angles, some under D/2 apart, and the outline has a vertex at
    x = 61 among them."""

    XS = (50.0, 52.0, 53.6, 57.4, 59.0, 63.4, 64.0)
    OUTLINE = poly((0, 0), (61, 0), (200, 0), (200, 200), (0, 200))

    def test_or4_placed_points_are_apart(self, fi: ModuleType) -> None:
        edges = (30.0, *self.XS, 90.0)
        strips = [box(x0, -50, x1, 60) for x0, x1 in pairwise(edges)]
        out = rule(fi, strips, self.OUTLINE)
        xy = np.unique(
            np.concatenate(
                [shapely.get_coordinates(line) for lines in out.lines for line in lines]
            ),
            axis=0,
        )
        placed = xy[on_outline(xy, self.OUTLINE) & (np.abs(rel_x(xy) - 60) < 20)]
        assert len(placed) >= 2
        vertices = shapely.get_coordinates(self.OUTLINE)
        is_vertex = [bool((np.abs(vertices - p).max(axis=1) == 0).any()) for p in placed]
        assert any(is_vertex)  # the outline vertex at 61 was used
        for i in range(len(placed)):
            for j in range(i + 1, len(placed)):
                if is_vertex[i] and is_vertex[j]:
                    continue
                assert np.hypot(*(placed[i] - placed[j])) > D / 2, (placed[i], placed[j])

    def test_or4_a_rounded_point_next_to_a_vertex_goes_onto_it(self, fi: ModuleType) -> None:
        """The design's second rounding (step 4): a multiple of D within D/2 of
        an outline vertex goes to the vertex. A side crossing at x = 58.4 is
        2.6 m (over D/2) from the vertex at 61, so it rounds to the multiple
        60, which is 1 m from 61: it must land on (61, 0) itself. Without the
        second rounding (the prototype's defect) it lands on (60, 0), 1 m from
        the vertex the path then passes through (mutation round, 20c-3)."""
        (lines,) = rule(fi, [box(58.4, -50, 120, 60)], self.OUTLINE).lines
        kept = inside(lines, self.OUTLINE)
        ends = np.array(
            [shapely.get_coordinates(p)[k] for p in shapely.get_parts(kept) for k in (0, -1)]
        )
        west = ends[on_outline(ends, self.OUTLINE) & (rel_x(ends) < 90)]
        assert len(west) == 1, west
        assert tuple(west[0]) == at(61, 0), west  # the vertex's own coordinates


class TestTheCorner:
    """OR4: a border 2 m inside the outline across one of its corners, at
    V = (100, 0), where the outline bends by 11.3°. Its points round onto V,
    and V stays in the path, so the whole stretch lies on the outline and is
    dropped. V lies at a multiple of D along the outline, so this fixture
    does not separate the second rounding from the first; the prototype's
    defect is pinned by ``TestRounding``'s
    ``test_or4_a_rounded_point_next_to_a_vertex_goes_onto_it``."""

    OUTLINE = poly((0, 0), (100, 0), (200, 20), (200, 200), (0, 200))
    NEAR = poly((30, 60), (30, 2), (100, 2), (170, 16), (170, 60))

    def test_or4_the_corner_stays_in_the_path(self, fi: ModuleType) -> None:
        assert self.OUTLINE.exterior.distance(Point(at(170, 16))) < D  # the premise
        (lines,) = rule(fi, [self.NEAR], self.OUTLINE).lines
        kept = inside(lines, self.OUTLINE)
        assert shapely.intersection(kept, box(60, -1, 140, 12)).is_empty


# ---------------------------------------------------------------- OR5


class TestCrossings:
    def test_or5_a_crossing_keeps_one_point_near_where_it_was(self, fi: ModuleType) -> None:
        (lines,) = rule(fi, [box(52.3, -50, 130.7, 60)]).lines
        kept = inside(lines, OUTLINE)
        for x in (52.3, 130.7):
            window = box(x - 10, -1, x + 10, 1)
            touch = shapely.intersection(shapely.intersection(kept, OUTLINE.exterior), window)
            points = np.unique(shapely.get_coordinates(touch), axis=0)
            assert len(points) == 1, points
            assert abs(rel_x(points)[0] - x) <= D / 2

    def test_or5_an_inlet_is_joined_straight(self, fi: ModuleType) -> None:
        """The outline runs 30 m down and back round a 2 m wide inlet at
        x = 100..102 (62 m against 2 x 2 + 2 D = 14 m), so the two points the
        border's stretch beside it rounds onto, (100, 0) and (102, 0), are
        joined straight across the inlet's mouth, inside the domain."""
        outline = poly(
            (0, 0), (100, 0), (100, -30), (102, -30), (102, 0), (200, 0), (200, 200), (0, 200)
        )
        (lines,) = rule(fi, [poly((60, 60), (60, 3), (140, 3), (140, 60))], outline).lines
        kept = inside(lines, outline)
        mouth = shapely.intersection(kept, box(99.5, -1, 102.5, 1))
        assert shapely.distance(mouth, Point(at(101, 0))) <= ON
        assert shapely.length(mouth) == pytest.approx(2.0, abs=1e-9)
        assert shapely.intersection(kept, box(99, -30, 103, -1)).is_empty

    def test_or5_a_shallow_notch_is_followed(self, fi: ModuleType) -> None:
        """Not an inlet: the outline's way round a notch 2 m wide and 1 m deep
        at x = 100..102 is 4 m, under 2 x 2 + 2 D = 14 m, so the points the
        border rounds onto, (100, 0) and (102, 0), are joined along the outline
        through its vertices (100, -1) and (102, -1). The stretch lies on the
        outline and leaves no linework either way; the polygon after the rule
        shows it: it fills the notch, 2 m². A straight join (the outline's
        vertices on the way dropped, mutation round, 20c-3) leaves the notch
        out. Scale: areas in m² at coordinates near 5e5 m, 1e-6 m² bound."""
        outline = poly(
            (0, 0), (100, 0), (100, -1), (102, -1), (102, 0), (200, 0), (200, 200), (0, 200)
        )
        out = rule(fi, [poly((60, 60), (60, 3), (140, 3), (140, 60))], outline)
        notch = box(100, -1, 102, 0)
        assert shapely.intersection(out.polygons[0], notch).area == pytest.approx(2.0, abs=1e-6)


# ---------------------------------------------------------------- OR6, OR7


class TestOff:
    def test_or6_d_0_is_a_no_op(self, fi: ModuleType) -> None:
        holed = Polygon(
            [at(0, -50), at(80, -50), at(80, 60), at(0, 60)],
            [[at(20, 3), at(40, 3), at(40, 20), at(20, 20)]],
        )
        near = poly((50, 60), (50, 3), (150, 3), (150, 60))
        out = rule(fi, [holed, near], d=0.0)
        for polygon, lines in zip([holed, near], out.lines, strict=True):
            rings = [polygon.exterior, *polygon.interiors]
            assert len(lines) == len(rings)
            for line, ring in zip(lines, rings, strict=True):
                assert isinstance(line, LineString)
                assert np.array_equal(shapely.get_coordinates(line), shapely.get_coordinates(ring))
        assert out.area_changed == 0.0


def staircase(steps: int, step: float = 2.0) -> Polygon:
    """30d's case: a raster-traced outline, ``steps`` stairs up and to the
    right from (0, 0), closed by the north-west corner."""
    xy = [(0.0, 0.0)]
    for i in range(steps):
        xy += [((i + 1) * step, i * step), ((i + 1) * step, (i + 1) * step)]
    xy.append((0.0, steps * step))
    return poly(*xy)


class TestNoBuffer:
    """OR7: D is refused at 100 m (the read region's margin), and the rule
    never buffers (or snaps to) the outline, on a staircase of 5 000 steps."""

    @pytest.mark.parametrize("d", [100.0, 150.0])
    def test_or7_d_at_or_above_100_is_refused(self, fi: ModuleType, d: float) -> None:
        with pytest.raises(ValueError, match="100"):
            fi.snap_to_outline([box(50, 3, 150, 60)], OUTLINE, d)

    def test_or7_just_under_100_runs(self, fi: ModuleType) -> None:
        rule(fi, [box(50, 3, 150, 60)], d=99.0)

    def test_or7_no_buffer_on_a_staircase(
        self, fi: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        outline = staircase(5_000)
        assert shapely.get_num_coordinates(outline) > 10_000  # the premise: 5 000 steps

        def refused(*_: Any, **__: Any) -> Any:
            raise AssertionError("the outline rule called buffer or snap")

        monkeypatch.setattr(shapely, "buffer", refused)
        monkeypatch.setattr(shapely, "snap", refused)
        crossing = box(4_990, 4_900, 5_010, 5_100)  # straddles the stairs at (5 000, 5 000)
        (lines,) = rule(fi, [crossing], outline).lines
        assert sum(shapely.length(line) for line in lines) > 0


# ---------------------------------------------------------------- OR8


def local(*points: tuple[float, float]) -> Polygon:
    """A polygon in local coordinates, without ``at``'s UTM-shaped offset."""
    return Polygon(points)


class TestTheReadRegion:
    """OR8 (ruling G4): the read region does not move the cuts. One polygon
    whose long west edge runs from 3 m inside the square's west side, nearly
    parallel to it, south out across the read region's edge, is clipped by
    three regions whose margins do not differ by a multiple of D; the rule's
    linework inside the domain is then the same for all three. The mutant:
    the cuts counted from the edge's own end, not along its line (the
    173.3 m clip then differs by 0.62 m, per the design).

    Local coordinates, as ``@architect`` checked it: the square 0 to
    1000 m, coordinates up to about 1 200 m from the origin, D = 5 m. The
    1e-9 m bound is the design's ("to rounding, not bit for bit"); at this
    scale a double resolves about 2e-13 m, and on ``7334b82b`` the three
    agree exactly."""

    SQUARE = local((0, 0), (1000, 0), (1000, 1000), (0, 1000))
    SHAPE = local((3, -537.3), (53.1, -537.3), (53.1, 512.7), (3.4, 512.7))
    MARGINS = (100.0, 150.0, 173.3)

    def clipped(self, margin: float) -> BaseGeometry:
        m = margin
        region = local((-m, -m), (1000 + m, -m), (1000 + m, 1000 + m), (-m, 1000 + m))
        return shapely.intersection(self.SHAPE, region)

    def test_or8_the_premise_the_clips_differ(self) -> None:
        a, b, c = (self.clipped(m) for m in self.MARGINS)
        assert not shapely.equals(a, b) and not shapely.equals(b, c)
        assert shapely.distance(self.SHAPE.exterior, self.SQUARE.boundary) <= D

    def test_or8_the_same_lines_inside_the_domain(self, fi: ModuleType) -> None:
        got = []
        for margin in self.MARGINS:
            (lines,) = rule(fi, [self.clipped(margin)], self.SQUARE).lines
            got.append(shapely.normalize(inside(lines, self.SQUARE)))
        assert not got[0].is_empty
        for margin, other in zip(self.MARGINS[1:], got[1:], strict=True):
            assert shapely.equals_exact(got[0], other, tolerance=1e-9), (margin, other.wkt)


# ---------------------------------------------------------------- OR9


class TestACutBesideAVertex:
    """OR9 (ruling T1): a cut that falls a hair from an edge's own vertex
    leaves no linework edge shorter than D/2. A border edge runs from W, 1 m
    inside the square's south side, to V, 5.5 m inside it (beyond D), placed
    so that a multiple of D along the edge's line, counted as ``_cut``
    counts (from the foot of the origin, towards the lexicographically last
    end), falls 2.4 cm short of V. That cut is more than D from the outline,
    so it does not move, and step 3 keeps it beside the moved W: on
    ``4d52dec0`` the linework gets a 2.4 cm edge, as Numedalslagen's did
    (2.9 cm). The mutant: ``_cut``'s end test removed.

    Local coordinates, the square 0 to 1000 m, D = 5 m; coordinates up to
    about 100 m, where a double resolves about 1e-14 m. The bound D/2 is the
    ruling's: after the fix an end piece is D/2 to 3D/2 long."""

    SQUARE = local((0, 0), (1000, 0), (1000, 1000), (0, 1000))
    W, V = (100.37, 1.0), (102.28, 5.5)
    SHAPE = local(W, V, (102.28, 60.2), (60.37, 60.2), (58.1, 1.0))

    def test_or9_the_premise_a_cut_lies_beside_v(self) -> None:
        w, v = np.array(self.W), np.array(self.V)
        u = (v - w) / np.hypot(*(v - w))  # W is the lexicographically first end
        t = np.arange(np.floor(w @ u / D) + 1, np.ceil(v @ u / D)) * D
        cuts = w + (t - w @ u)[:, None] * u
        gaps = np.hypot(*(cuts - v).T)
        assert len(cuts) == 1 and 0 < gaps[0] < 0.05, gaps
        assert cuts[0][1] > D  # beyond D from the south side: it does not move

    @pytest.mark.parametrize("reverse", [False, True], ids=["same", "reversed"])
    def test_or9_no_edge_shorter_than_half_d(self, fi: ModuleType, reverse: bool) -> None:
        shape = Polygon(self.SHAPE.exterior.coords[::-1]) if reverse else self.SHAPE
        (lines,) = rule(fi, [shape], self.SQUARE).lines
        kept = inside(lines, self.SQUARE)
        short = [
            (a, b) for a, b in segments(list(shapely.get_parts(kept))) if math.dist(a, b) < D / 2
        ]
        assert not short, short


# ---------------------------------------------------------------- T2


class TestOnlyWhatMovedIsRebuilt:
    """Ruling T2: a polygon of two parts, one far from the outline and one
    with two separate stretches 3 m inside its south side. The far part
    comes back as it went in (``equals_exact`` at 0: passed through, not
    through a union), and ``area_changed`` equals the overlay measure the
    test computes itself: the input's and output's symmetric difference,
    intersected with the outline.

    Local coordinates, the square 0 to 1000 m, D = 5 m; areas up to about
    600 m² from coordinates up to 400 m, where an overlay's rounding is far
    under the 1e-6 m² bound (the ruling's)."""

    SQUARE = local((0, 0), (1000, 0), (1000, 1000), (0, 1000))
    FAR = local((300, 300), (400, 300), (400, 400), (300, 400))
    NEAR = local((20, 3), (40, 3), (40, 30), (80, 30), (80, 3), (100, 3), (100, 60), (20, 60))
    BOTH = shapely.MultiPolygon([NEAR, FAR])

    def test_t2_the_premise_two_stretches_near_one_far(self) -> None:
        band = shapely.intersection(self.NEAR.exterior, local((0, 0), (1000, 0), (1000, D), (0, D)))
        assert len(shapely.get_parts(shapely.line_merge(band))) == 2
        assert shapely.distance(self.FAR, self.SQUARE.boundary) > D

    def test_t2_the_far_part_passes_through(self, fi: ModuleType) -> None:
        out = rule(fi, [self.BOTH], self.SQUARE)
        (after,) = out.polygons
        far = [q for q in shapely.get_parts(after) if q.intersects(self.FAR.centroid)]
        assert len(far) == 1, after.wkt
        assert shapely.equals_exact(far[0], self.FAR, tolerance=0), far[0].wkt

    def test_t2_the_area_is_the_overlay_measure(self, fi: ModuleType) -> None:
        out = rule(fi, [self.BOTH], self.SQUARE)
        (after,) = out.polygons
        moved = shapely.intersection(shapely.symmetric_difference(self.BOTH, after), self.SQUARE)
        assert moved.area > 0
        assert out.area_changed == pytest.approx(moved.area, abs=1e-6)


# ---------------------------------------------------------------- P1


class TestARingWithNoKeptVertex:
    """Ruling P1 ("Rulings on the fix round"): a ring none of whose input
    vertices stays where it was adds the symmetric difference of its old and
    its new polygon to ``area_changed``, not the two polygons whole. A 1 km
    square outline; A, the 2 x 10 m rectangle (0,100)-(2,110) on its west
    side; B, the square less A, whose eight vertices all lie within D of the
    outline, so none is kept. A goes onto the outline and comes back empty;
    B comes back as the whole square; what changed polygon is A's 20 m².
    On ``8716f062`` ``area_changed`` is 1 000 000 m² (B's old and new
    polygon together cover the square).

    Local coordinates, the square 0 to 1000 m, D = 5 m; areas up to 1e6 m²
    from coordinates up to 1000 m, where a double resolves about 1e-13 m
    and an overlay of these axis-parallel rings is exact. The 1e-6 m² bound
    is the ruling's."""

    SQUARE = local((0, 0), (1000, 0), (1000, 1000), (0, 1000))
    A = local((0, 100), (2, 100), (2, 110), (0, 110))
    B = local((0, 0), (1000, 0), (1000, 1000), (0, 1000), (0, 110), (2, 110), (2, 100), (0, 100))

    def test_p1_the_premise_a_coverage_with_no_vertex_of_b_beyond_d(self) -> None:
        assert self.B.is_valid and self.B.area == pytest.approx(1e6 - 20, abs=1e-6)
        assert shapely.equals(shapely.union_all([self.A, self.B]), self.SQUARE)
        vertices = shapely.points(shapely.get_coordinates(self.B))
        assert (np.asarray(shapely.distance(self.SQUARE.boundary, vertices)) <= D).all()

    def test_p1_a_goes_empty_and_b_becomes_the_square(self, fi: ModuleType) -> None:
        a, b = rule(fi, [self.A, self.B], self.SQUARE).polygons
        assert a.is_empty, a.wkt
        assert shapely.equals(b, self.SQUARE), b.wkt

    def test_p1_the_area_moved_is_a_alone(self, fi: ModuleType) -> None:
        out = rule(fi, [self.A, self.B], self.SQUARE)
        assert out.area_changed == pytest.approx(20.0, abs=1e-6)
