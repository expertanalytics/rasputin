"""Exact land-cover class fractions per triangle from class polygons
(`tin_engine.polygon_fractions`): `subdivide` and `assign`.

`docs/increments/105-cell-table.md`, section 4.5 (fractions from polygons, the
subdivision, the int32 class numbers) and section 6, T6. The hand figures are
the design's: triangles `A = (0,0), (2,0), (0,2)` and `B = (2,0), (2,2),
(0,2)`, plan area 2 each, against class 311 = the strip `[0,1] x [0,2]` and
512 = the square `[1,2] x [0,1]`. Every hand figure is a binary fraction, so
the float32 fractions are compared exactly; the 1,000-gon's shares come from
the design's probe (`probe_cells.py`, shapely) and are compared at 1e-6, the
float32 rounding (6e-8) with room.

THE ORACLE where the figures are not by hand is GEOS's
`intersection(...).area` of the whole polygon with the whole triangle, never
the producer's pieces.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- `polygon_fractions.subdivide(polygon, most=64) -> Sequence[BaseGeometry]`:
  every piece has at most `most` coordinates (`shapely.get_num_coordinates`),
  none is empty, and a polygon at or under `most` is returned as the one piece.
- `polygon_fractions.assign(vertices, triangles, polygons, *, known, most=64)
  -> land_cover_fractions.Assigned`: `vertices` `(N, 2)` or `(N, 3)` in the
  polygons' CRS (z ignored), `triangles` `(T, 3)` indices, `polygons` an
  iterable of `(code, Polygon | MultiPolygon)` pairs, `known` the meshing
  rules' class codes. Uncovered area is class 0, a class like any other.
- A code not in `known` is a `ValueError`, "class <code> is not in the
  meshing rules" (the rules decide, not the legend); a non-finite
  polygon coordinate is a `ValueError` saying "finite"; a zero-area triangle
  is a `ValueError` (the shared `fractions_of` refusal, "no class area").
- No polygons: every triangle is class 0, fraction 1. No triangles: an empty
  table with one offset.
- `most` below 8 is a `ValueError` naming 8; a piece is not halved past
  `DEPTH` (64) levels; where the rectangle clip loses area the halves still
  tile the piece; an uncovered rest of at most 1e-9 of its triangle is no
  class 0 entry, one above it is.
"""

from __future__ import annotations

import importlib
from itertools import pairwise
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import shapely
from numpy.testing import assert_array_equal
from shapely.geometry import MultiPolygon, Point, Polygon, box

from fraction_fixtures import rows_of

#: A and B tile the square [0, 2]^2 along the line x + y = 2.
VERTICES = np.array([[0.0, 0.0, 5.0], [2.0, 0.0, 6.0], [2.0, 2.0, 7.0], [0.0, 2.0, 8.0]])
TRIANGLES = np.array([[0, 1, 3], [1, 2, 3]], dtype=np.int64)
STRIP = box(0.0, 0.0, 1.0, 2.0)
SQUARE = box(1.0, 0.0, 2.0, 1.0)
KNOWN = frozenset({311, 312, 512})
#: T6: A is 1.5 in the strip and 0.5 in the square; B is 0.5, 0.5 and 1.0 uncovered.
OFFSETS = np.array([0, 2, 5], dtype=np.uint32)
CLASSES = np.array([311, 512, 0, 311, 512], dtype=np.int32)
FRACTIONS = np.array([0.75, 0.25, 0.5, 0.25, 0.25], dtype=np.float32)
CODES = np.array([311, 0], dtype=np.int32)
#: The 1,000-gon: inscribed in [0, 2]^2, which A and B tile; 1,001 coordinates.
RING = Point(1.0, 1.0).buffer(1.0, quad_segs=250)
#: Its share of A and of B, half of it each (pi/4 = 0.785398 for the true circle).
RING_SHARE = 0.785393
#: A float32 fraction against a float64 figure: 6e-8 rounding, with room.
SHARE = 1e-6
#: The pieces' areas against the polygon's, relative; GEOS's clip at
#: coordinates below 3 (and below 7e6 in the far frame) loses less.
AREA = 1e-9
#: UTM-shaped, off every round number: section 3A's extreme scale.
FAR = np.array([500_000.3, 6_600_000.7])


@pytest.fixture(scope="module")
def pf() -> ModuleType:
    return importlib.import_module("tin_engine.polygon_fractions")


def assigned(pf: ModuleType, polygons: Any, **changes: Any) -> Any:
    kwargs: dict[str, Any] = {"known": KNOWN} | changes
    return pf.assign(VERTICES, TRIANGLES, polygons, **kwargs)


def table_of(result: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = result.table
    return t.offsets, t.classes, t.fractions


def shifted(geometry: Any, by: np.ndarray = FAR) -> Any:
    return shapely.transform(geometry, lambda xy: xy + by)


# ------------------------------------------------------------ T6, the hand figures


class TestHandFigures:
    def test_t6_the_strip_and_the_square_over_a_and_b(self, pf: ModuleType) -> None:
        offsets, classes, fractions = table_of(assigned(pf, [(311, STRIP), (512, SQUARE)]))
        assert_array_equal(offsets, OFFSETS)
        assert_array_equal(classes, CLASSES)
        assert_array_equal(fractions, FRACTIONS)

    def test_t6_the_code_is_the_largest_share_and_uncovered_is_class_0(
        self, pf: ModuleType
    ) -> None:
        result = assigned(pf, [(311, STRIP), (512, SQUARE)])
        assert_array_equal(result.codes, CODES)
        assert result.codes.dtype == np.int32

    def test_t6_the_class_numbers_are_int32(self, pf: ModuleType) -> None:
        """4.5: CORINE's 111 to 523 and AR50's 3031 to 3099 do not fit uint8."""
        _, classes, fractions = table_of(assigned(pf, [(311, STRIP), (512, SQUARE)]))
        assert (classes.dtype, fractions.dtype) == (np.int32, np.float32)

    def test_t6_the_input_order_does_not_change_the_table(self, pf: ModuleType) -> None:
        """Classes ascend within a row whatever order the polygons come in."""
        result = table_of(assigned(pf, [(512, SQUARE), (311, STRIP)]))
        for got, want in zip(result, (OFFSETS, CLASSES, FRACTIONS), strict=True):
            assert_array_equal(got, want)

    def test_t6_two_polygons_of_one_class_over_one_triangle_give_one_entry(
        self, pf: ModuleType
    ) -> None:
        """The strip in two halves, both 311: A meets both (1.0 and 0.5), and
        its row still has one 311 entry of 0.75."""
        halves = [(311, box(0.0, 0.0, 1.0, 1.0)), (311, box(0.0, 1.0, 1.0, 2.0))]
        offsets, classes, fractions = table_of(assigned(pf, [*halves, (512, SQUARE)]))
        assert_array_equal(offsets, OFFSETS)
        assert_array_equal(classes, CLASSES)
        assert_array_equal(fractions, FRACTIONS)

    def test_a_multipolygon_is_one_class(self, pf: ModuleType) -> None:
        """The strip and the box [1.5, 2] x [0, 0.5] as one 311 multipolygon:
        the box is 0.125 in A and 0.125 in B (the line x + y = 2 halves it),
        so A is 1.625 of 2 in 311 and B 0.625."""
        two = MultiPolygon([STRIP, box(1.5, 0.0, 2.0, 0.5)])
        offsets, classes, fractions = table_of(assigned(pf, [(311, two)]))
        assert_array_equal(offsets, [0, 2, 4])
        assert_array_equal(classes, [0, 311, 0, 311])
        want = np.array([0.1875, 0.8125, 0.6875, 0.3125], dtype=np.float32)
        assert_array_equal(fractions, want)

    def test_a_hole_is_uncovered(self, pf: ModuleType) -> None:
        """311 over all of [0, 2]^2 but a hole [0.25, 0.75]^2, inside A: A is
        1.75 in 311 and 0.25 uncovered, B all 311."""
        holed = Polygon(
            box(0.0, 0.0, 2.0, 2.0).exterior.coords,
            [box(0.25, 0.25, 0.75, 0.75).exterior.coords],
        )
        offsets, classes, fractions = table_of(assigned(pf, [(311, holed)]))
        assert_array_equal(offsets, [0, 2, 3])
        assert_array_equal(classes, [0, 311, 311])
        assert_array_equal(fractions, np.array([0.125, 0.875, 1.0], dtype=np.float32))

    def test_two_dimensional_vertices_give_the_same_table(self, pf: ModuleType) -> None:
        result = pf.assign(VERTICES[:, :2], TRIANGLES, [(311, STRIP), (512, SQUARE)], known=KNOWN)
        assert_array_equal(result.table.classes, CLASSES)
        assert_array_equal(result.table.fractions, FRACTIONS)

    def test_far_from_the_origin_the_shares_hold(self, pf: ModuleType) -> None:
        """Section 3A's extreme scales: the same figures 6,600 km from the
        origin, triangles of 2 m2. GEOS's clip there is good to about 1e-9 m,
        far under the float32 rounding of a fraction."""
        far = VERTICES.copy()
        far[:, :2] += FAR
        result = pf.assign(
            far, TRIANGLES, [(311, shifted(STRIP)), (512, shifted(SQUARE))], known=KNOWN
        )
        assert_array_equal(result.table.offsets, OFFSETS)
        assert_array_equal(result.table.classes, CLASSES)
        np.testing.assert_allclose(result.table.fractions, FRACTIONS, rtol=0, atol=SHARE)
        assert_array_equal(result.codes, CODES)


# ------------------------------------------------------------ edges of the input


class TestEdgesOfTheInput:
    def test_a_code_absent_from_the_rules_is_refused_naming_it(self, pf: ModuleType) -> None:
        with pytest.raises(ValueError, match=r"^class 9999 is not in the land-cover class file$"):
            assigned(pf, [(311, STRIP), (9999, SQUARE)])

    def test_no_polygons_leave_every_triangle_class_0(self, pf: ModuleType) -> None:
        offsets, classes, fractions = table_of(assigned(pf, []))
        assert_array_equal(offsets, [0, 1, 2])
        assert_array_equal(classes, [0, 0])
        assert_array_equal(fractions, np.ones(2, dtype=np.float32))

    def test_no_triangles_give_an_empty_table(self, pf: ModuleType) -> None:
        result = pf.assign(VERTICES, np.zeros((0, 3), dtype=np.int64), [(311, STRIP)], known=KNOWN)
        assert result.table.triangles == 0
        assert_array_equal(result.table.offsets, [0])
        assert result.table.classes.dtype == np.int32
        assert result.codes.shape == (0,)

    def test_a_polygon_far_from_every_triangle_adds_nothing(self, pf: ModuleType) -> None:
        _, classes, fractions = table_of(
            assigned(pf, [(311, STRIP), (512, SQUARE), (312, box(10.0, 10.0, 11.0, 11.0))])
        )
        assert_array_equal(classes, CLASSES)
        assert_array_equal(fractions, FRACTIONS)

    def test_a_polygon_touching_a_triangle_only_along_an_edge_adds_no_entry(
        self, pf: ModuleType
    ) -> None:
        """Section 3A's coincidence: [2, 3] x [0, 2] shares B's edge x = 2
        and A's corner (2, 0); zero area is no entry."""
        offsets, classes, _ = table_of(
            assigned(pf, [(311, STRIP), (512, SQUARE), (312, box(2.0, 0.0, 3.0, 2.0))])
        )
        assert_array_equal(offsets, OFFSETS)
        assert 312 not in classes.tolist()

    @pytest.mark.parametrize("bad", [np.nan, np.inf], ids=["nan", "inf"])
    def test_a_non_finite_polygon_coordinate_is_refused(self, pf: ModuleType, bad: float) -> None:
        broken = Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, bad), (0.0, 2.0)])
        with pytest.raises(ValueError, match="finite"):
            assigned(pf, [(311, broken)])

    def test_a_rounding_rest_is_no_class_0_entry(self, pf: ModuleType) -> None:
        """A 100 m triangle 6,600 km from the origin inside a 500 m 1,000-gon
        cut in pieces: the pieces' areas miss the triangle's by +5e-13 of it
        (checked here, so the case can fail), under the 1e-9 floor, so the row
        is 311 alone."""
        tri = np.array([[500_010.25, 6_600_020.5], [500_090.75, 6_600_035.25],
                        [500_040.5, 6_600_095.75]])  # fmt: skip
        big = Point(500_037.1, 6_600_041.9).buffer(500.0, quad_segs=250)
        shape = Polygon(tri)
        pieces = np.array(pf.subdivide(big), dtype=object)
        rest = shape.area - float(np.sum(shapely.area(shapely.intersection(shape, pieces))))
        assert 0.0 < rest <= 1e-9 * shape.area
        offsets, classes, fractions = table_of(
            pf.assign(tri, np.array([[0, 1, 2]]), [(311, big)], known=KNOWN)
        )
        assert_array_equal(offsets, [0, 1])
        assert_array_equal(classes, [311])
        assert_array_equal(fractions, [1.0])

    def test_an_uncovered_rest_above_the_floor_is_class_0(self, pf: ModuleType) -> None:
        """311 over [0, 2]^2 but a hole of 2e-4 x 2e-4 inside A: 2e-8 of A
        uncovered, twenty times the 1e-9 floor, is a class 0 entry."""
        holed = Polygon(
            box(0.0, 0.0, 2.0, 2.0).exterior.coords,
            [box(0.5, 0.5, 0.5002, 0.5002).exterior.coords],
        )
        offsets, classes, fractions = table_of(assigned(pf, [(311, holed)]))
        assert_array_equal(offsets, [0, 2, 3])
        assert_array_equal(classes, [0, 311, 311])
        np.testing.assert_allclose(fractions[0], 2e-8, rtol=1e-4)

    def test_a_zero_area_triangle_is_refused(self, pf: ModuleType) -> None:
        """Three collinear corners: no area to share out, as the raster route
        refuses."""
        line = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]])
        with pytest.raises(ValueError, match="no class area"):
            pf.assign(line, np.array([[0, 1, 2]]), [(311, STRIP)], known=KNOWN)


# ------------------------------------------------------------ the subdivision


def coordinate_counts(pieces: Any) -> list[int]:
    return [int(shapely.get_num_coordinates(p)) for p in pieces]


def area_sum(pieces: Any) -> float:
    return float(sum(p.area for p in pieces))


class TestSubdivide:
    def test_the_1000_gon_gives_32_pieces_of_at_most_58_coordinates(self, pf: ModuleType) -> None:
        assert shapely.get_num_coordinates(RING) == 1001
        pieces = list(pf.subdivide(RING, most=64))
        assert len(pieces) == 32
        assert max(coordinate_counts(pieces)) == 58
        assert not any(p.is_empty for p in pieces)

    def test_the_pieces_tile_the_polygon(self, pf: ModuleType) -> None:
        """The areas sum to the ring's, and no two pieces overlap: the union's
        area is the sum."""
        pieces = list(pf.subdivide(RING, most=64))
        assert abs(area_sum(pieces) / RING.area - 1.0) <= AREA
        assert abs(shapely.union_all(pieces).area / RING.area - 1.0) <= AREA
        assert shapely.symmetric_difference(shapely.union_all(pieces), RING).area <= AREA

    def test_the_default_is_64(self, pf: ModuleType) -> None:
        assert len(list(pf.subdivide(RING))) == 32

    @pytest.mark.parametrize("most", [8, 16, 64, 256])
    def test_every_piece_is_at_or_under_most(self, pf: ModuleType, most: int) -> None:
        pieces = list(pf.subdivide(RING, most=most))
        assert max(coordinate_counts(pieces)) <= most
        assert abs(area_sum(pieces) / RING.area - 1.0) <= AREA

    def test_a_polygon_under_most_comes_back_as_itself(self, pf: ModuleType) -> None:
        pieces = list(pf.subdivide(STRIP, most=64))
        assert len(pieces) == 1
        assert shapely.equals_exact(pieces[0], STRIP, tolerance=0.0)

    def test_a_polygon_of_exactly_most_coordinates_is_one_piece(self, pf: ModuleType) -> None:
        assert len(list(pf.subdivide(RING, most=1001))) == 1
        assert len(list(pf.subdivide(RING, most=1000))) > 1

    @pytest.mark.parametrize("most", [7, 5])
    def test_most_below_8_is_refused(self, pf: ModuleType, most: int) -> None:
        """A rectangle less a notch at one vertex is 7 vertices and the
        closing point; halving cannot always get under 8."""
        with pytest.raises(ValueError, match=f"at least 8, got {most}"):
            pf.subdivide(STRIP, most=most)

    def test_a_piece_is_not_halved_past_depth_levels(
        self, pf: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """With the backstop at 2 levels the 1,000-gon stops at 4 quarters,
        each over `most`, and they still tile it."""
        monkeypatch.setattr(pf, "DEPTH", 2)
        pieces = list(pf.subdivide(RING, most=8))
        assert len(pieces) == 4
        assert min(coordinate_counts(pieces)) > 8
        assert abs(area_sum(pieces) / RING.area - 1.0) <= AREA

    def test_vertices_an_ulp_apart_stop_at_64_levels(self, pf: ModuleType) -> None:
        """Section 3A's near-coincidence: a notch of 20 vertices one ulp apart
        near (0.3, 0.3) in [0, 1]^2. 64 halvings do not separate them, so one
        piece keeps more than `most` coordinates; the pieces still tile."""
        zigzag, x, y = [], 0.3, [0.3, float(np.nextafter(0.3, 1.0))]
        for k in range(20):
            zigzag.append((x, y[k % 2]))
            x = float(np.nextafter(x, 1.0))
        notched = Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), *zigzag[::-1], (0.0, 1.0)])
        assert notched.is_valid
        pieces = list(pf.subdivide(notched, most=8))
        assert max(coordinate_counts(pieces)) > 8
        assert abs(area_sum(pieces) / notched.area - 1.0) <= AREA

    def test_a_hole_touching_the_ring_on_the_cut_line_keeps_its_area(self, pf: ModuleType) -> None:
        """[0, 2] x [0, 4] less the hole (0, 2), (1, 1.5), (1, 2.5), which
        touches the outer ring at (0, 2), on the first cut y = 2: GEOS's
        rectangle clip gives the upper half 0.25 instead of 3.75; the pieces
        must still sum to 7.5."""
        holed = Polygon(box(0.0, 0.0, 2.0, 4.0).exterior.coords, [[(0, 2), (1, 1.5), (1, 2.5)]])
        assert holed.is_valid and shapely.get_num_coordinates(holed) == 9
        assert shapely.clip_by_rect(holed, 0.0, 2.0, 2.0, 4.0).area == 0.25
        pieces = list(pf.subdivide(holed, most=8))
        assert max(coordinate_counts(pieces)) <= 8
        assert abs(area_sum(pieces) / holed.area - 1.0) <= AREA
        assert abs(shapely.union_all(pieces).area / holed.area - 1.0) <= AREA

    def test_holes_count_and_are_kept(self, pf: ModuleType) -> None:
        """A square of 5 coordinates around a 1,000-gon hole: 1,006 in all, so
        it is cut, and the pieces' areas are the square's less the hole's."""
        holed = Polygon(box(-1.0, -1.0, 3.0, 3.0).exterior.coords, [RING.exterior.coords])
        assert shapely.get_num_coordinates(holed) == 1006
        pieces = list(pf.subdivide(holed, most=64))
        assert len(pieces) > 1 and max(coordinate_counts(pieces)) <= 64
        assert abs(area_sum(pieces) / holed.area - 1.0) <= AREA

    def test_a_long_collinear_edge_is_cut(self, pf: ModuleType) -> None:
        """Section 3A's collinear sequence: 1,000 vertices on the bottom edge
        of [0, 1]^2; every piece at most 64 coordinates, the area kept."""
        bottom = [(k / 1000.0, 0.0) for k in range(1001)]
        square = Polygon([*bottom, (1.0, 1.0), (0.0, 1.0)])
        pieces = list(pf.subdivide(square, most=64))
        assert max(coordinate_counts(pieces)) <= 64
        assert abs(area_sum(pieces) - 1.0) <= AREA

    def test_a_multipolygon_is_cut_part_by_part(self, pf: ModuleType) -> None:
        two = MultiPolygon([RING, shifted(RING, np.array([5.0, 0.0]))])
        pieces = list(pf.subdivide(two, most=64))
        assert max(coordinate_counts(pieces)) <= 64
        assert abs(area_sum(pieces) / two.area - 1.0) <= AREA

    def test_far_from_the_origin_the_pieces_tile_the_polygon(self, pf: ModuleType) -> None:
        far = shifted(RING)
        pieces = list(pf.subdivide(far, most=64))
        assert len(pieces) == 32
        assert abs(area_sum(pieces) / far.area - 1.0) <= AREA


class TestSubdividedShares:
    def test_the_1000_gons_share_is_the_same_over_the_pieces_and_the_whole(
        self, pf: ModuleType
    ) -> None:
        """T6: half the 1,000-gon in A and half in B, 0.785393 of each, cut
        in 32 pieces or not cut at all; GEOS's whole-polygon clip agrees."""
        oracle = [
            float(shapely.intersection(Polygon(VERTICES[t, :2]), RING).area) / 2.0
            for t in TRIANGLES
        ]
        np.testing.assert_allclose(oracle, [RING_SHARE, RING_SHARE], rtol=0, atol=SHARE)
        for most in (64, 2000):
            rows = rows_of(*table_of(assigned(pf, [(311, RING)], most=most)))
            for t, row in enumerate(rows):
                assert set(row) == {0, 311}, (most, t, row)
                assert abs(row[311] - oracle[t]) <= SHARE, (most, t, row)
                assert abs(row[0] - (1.0 - oracle[t])) <= SHARE, (most, t, row)

    def test_each_row_sums_to_one(self, pf: ModuleType) -> None:
        corner = box(1.9, 1.9, 2.0, 2.0)  # outside the 1,000-gon
        assert not corner.intersects(RING)
        offsets, _, fractions = table_of(assigned(pf, [(311, RING), (512, corner)]))
        for a, b in pairwise(offsets.astype(np.int64)):
            assert abs(float(fractions[a:b].astype(np.float64).sum()) - 1.0) <= SHARE
