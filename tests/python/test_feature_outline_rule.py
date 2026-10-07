"""Ola's outline rule: land-cover borders within D of the outline go onto it (20c-3).

``docs/increments/20c-soft-quality.md``, M6 (c), "Steps 2 to 4" and "Tests
``@tester`` writes red first", 20c-3: OR1 to OR8. Invariant-critical
(mutation testing required: the rule rewrites input borders); the mutants
the design names are run after the green step, on the built code:

- the cut not in a fixed order: OR2;
- the vertex dropped at a rounded point: OR4 (``TestTheCorner``);
- the stretches on the outline kept: OR1, OR3;
- the rounding removed: OR4 (``TestRounding``);
- the inlet test removed: OR5 (``test_or5_an_inlet_is_joined_straight``).
- the cuts counted from the edge's own end, not along its line: OR8
  (``TestTheReadRegion``, added by ruling G4).

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

RED at the commit that adds this file: ``feature_input`` has no
``snap_to_outline``.
"""

from __future__ import annotations

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


class TestTheCorner:
    """OR4: a border 2 m inside the outline across one of its corners, at
    V = (100, 0), where the outline bends by 11.3°. Its points round onto V,
    and V stays in the path, so the whole stretch lies on the outline and is
    dropped; losing V (round 3's prototype defect) leaves a chord off it."""

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
