"""The pure steps of `mesh --catchments`: select, repair, simplify, the band
from the budget, the ids.

`docs/increments/106-wednesday-cells.md`, section 4, tests 1 to 4 and 8's
pure half. Every figure was measured with shapely 2.2.0 (GEOS 3.14.1) and
reproduced for this suite at the design's head. All units in metres, in a
projected CRS; coordinates below 3,100 m, so every bound below holds with
some twelve digits to spare.

PINNED HERE, where the design names the steps but not their spelling
(listed for `@architect`), all in `tin_engine.sub_catchments`:

- `select(units: Sequence[BaseGeometry], domain: BaseGeometry) -> list[int]`:
  the indices, in order, of the units whose inner point lies in `domain`.
- `repair(units: Sequence[BaseGeometry]) -> list[BaseGeometry]`: the units
  after `coverage_clean` at 1 m, in order.
- `simplify(units: Sequence[BaseGeometry], band_m: float) ->
  list[BaseGeometry]`: the units, in order, after `simplify_borders` with
  the frame, clearance 1 m.
- `band_m(area_m2: float, triangles: int) -> float`.
- `check_ids(fids: Sequence[object], name: str) -> list[int]`: refuses a fid
  that is no integer from 1 to 2**31 - 1 with a `ValueError` saying
  `<name>: feature <fid> has no id Shyft's catchment_id can hold`.
- `cell_ids(units: Sequence[tuple[int, BaseGeometry]], centroids) ->
  tuple[np.ndarray, int]`: each centroid's unit id (the nearest unit's when
  it lies in none) and how many lay in none.
"""

from __future__ import annotations

import numpy as np
import pytest
import shapely
from numpy.testing import assert_allclose, assert_array_equal
from shapely.geometry import Polygon, box

from tin_engine import sub_catchments as sc

# ------------------------------------------------------------ 1. select


def test_select_by_the_inner_point_three_squares() -> None:
    squares = [box(0, 0, 1000, 1000), box(1000, 0, 2000, 1000), box(2000, 0, 3000, 1000)]
    assert sc.select(squares, box(0, 0, 2050, 1000)) == [0, 1]


#: Its centroid (1500, 454.35) lies in its notch, outside it; its inner point is (500, 650).
U_SHAPE = box(0, 0, 3000, 1000).difference(box(1000, 300, 2000, 1000))


def test_the_u_fixture_is_as_the_design_measured() -> None:
    assert U_SHAPE.area == 2_300_000.0
    assert not U_SHAPE.contains(U_SHAPE.centroid)
    assert U_SHAPE.point_on_surface().equals(shapely.Point(500, 650))


def test_a_domain_holding_the_centroid_but_not_the_inner_point_selects_nothing() -> None:
    """Kills a centroid rule and an intersects rule: this domain holds the
    centroid and crosses the unit."""
    domain = box(900, 0, 2100, 1000)
    assert domain.contains(U_SHAPE.centroid) and domain.intersects(U_SHAPE)
    assert sc.select([U_SHAPE], domain) == []


def test_a_domain_holding_the_inner_point_but_not_the_centroid_selects_it() -> None:
    domain = box(0, 0, 1000, 1000)
    assert not domain.contains(U_SHAPE.centroid)
    assert sc.select([U_SHAPE], domain) == [0]


# ------------------------------------------------------------ 2. repair


def test_repair_closes_a_slit_between_two_units() -> None:
    a = Polygon([(0, 0), (1000, 0), (1000, 500), (1000, 1000), (0, 1000)])
    b = Polygon([(1000, 0), (2000, 0), (2000, 1000), (1000, 1000), (1000.1, 500)])
    before = shapely.coverage_union_all([a, b])
    assert before.geom_type == "MultiPolygon" and before.area == pytest.approx(1_999_950, rel=1e-12)
    repaired = sc.repair([a, b])
    assert len(repaired) == 2
    after = shapely.coverage_union_all(repaired)
    assert after.geom_type == "Polygon"
    assert after.area == pytest.approx(2_000_000, rel=1e-12)
    assert repaired[1].area == pytest.approx(1_000_000, rel=1e-12)


# ------------------------------------------------------------ 3. ids


def test_ids_are_the_units_fids() -> None:
    assert sc.check_ids([7, 9], "units.gpkg") == [7, 9]
    assert sc.check_ids([1, 2**31 - 1], "units.gpkg") == [1, 2**31 - 1]
    units = [(7, box(0, 0, 1000, 1000)), (9, box(1000, 0, 2000, 1000))]
    ids, strays = sc.cell_ids(units, np.array([[500.0, 500.0], [1500.0, 300.0], [10.0, 990.0]]))
    assert_array_equal(ids, [7, 9, 7])
    assert strays == 0


def test_a_centroid_in_no_unit_takes_the_nearest_and_is_counted() -> None:
    units = [(7, box(0, 0, 1000, 1000)), (9, box(1000, 0, 2000, 1000))]
    ids, strays = sc.cell_ids(units, np.array([[2010.0, 500.0], [500.0, 500.0]]))
    assert_array_equal(ids, [9, 7])
    assert strays == 1


@pytest.mark.parametrize("fid", [0, 2**31, "x"], ids=["zero", "2^31", "text"])
def test_a_fid_shyft_cannot_hold_is_refused(fid: object) -> None:
    with pytest.raises(ValueError, match="has no id Shyft's catchment_id can hold") as caught:
        sc.check_ids([7, fid], "units.gpkg")
    assert str(caught.value).startswith(f"units.gpkg: feature {fid} "), caught.value


# ------------------------------------------------------------ 4. simplify


def zigzag_border(x: float) -> list[tuple[float, float]]:
    """x = `x` from y 0 to 1000 in 50 m steps, every other vertex 20 m east."""
    return [(x + 20.0 * (k % 2), 50.0 * k) for k in range(21)]


def zigzag_top(x0: float, x1: float) -> list[tuple[float, float]]:
    """y = 1000 from `x0` to `x1` in 50 m steps, every other vertex 20 m up."""
    return [(x0 + 50.0 * k, 1000.0 + 20.0 * (k % 2)) for k in range(int((x1 - x0) / 50) + 1)]


def three_zigzag_units() -> list[Polygon]:
    left = Polygon([(0, 0), *zigzag_border(1000), *zigzag_top(0, 1000)[::-1][1:]])
    middle = Polygon(
        [(1000, 0), (2000, 0), *zigzag_border(2000)[1:], *zigzag_top(1000, 2000)[::-1][1:],
         *zigzag_border(1000)[::-1][1:-1]]
    )  # fmt: skip
    right = Polygon(
        [(2000, 0), (3000, 0), *zigzag_top(2000, 3000)[::-1], *zigzag_border(2000)[::-1][1:-1]]
    )
    return [left, middle, right]


def vertices(polygons: list[Polygon]) -> int:
    return sum(len(p.exterior.coords) for p in polygons)


def test_the_zigzag_fixture_is_as_the_design_measured() -> None:
    units = three_zigzag_units()
    assert [u.area for u in units] == [1_020_000.0, 1_010_000.0, 1_000_000.0]
    assert bool(shapely.coverage_is_valid(units))
    assert vertices(units) == 148
    assert len(shapely.coverage_union_all(units).exterior.coords) == 66


def test_simplify_keeps_the_coverage_and_moves_the_outline() -> None:
    """The frame makes the union's outline a border: without it the outline keeps its 66."""
    units = three_zigzag_units()
    out = sc.simplify(units, 50.0)
    assert len(out) == 3
    assert_allclose([p.area for p in out], [u.area for u in units], rtol=1e-12, atol=0)
    assert bool(shapely.coverage_is_valid(out))
    union = shapely.coverage_union_all(out)
    assert union.geom_type == "Polygon" and len(union.interiors) == 0
    assert union.area == pytest.approx(3_030_000, rel=1e-12)
    assert vertices(out) == 110
    assert len(union.exterior.coords) == 28


# ------------------------------------------------------------ 8. the band, pure


def test_the_band_is_an_eighth_of_the_equilateral_edge_at_the_budget() -> None:
    """h = sqrt(4 A / (sqrt(3) N)): at 4 km2 and 400 triangles h is 151.967 m."""
    band = sc.band_m(4_000_000, 400)
    assert band == pytest.approx(18.995892141289815, rel=1e-9)
    assert sc.band_m(4_000_000, 40_000) == pytest.approx(band / 10, rel=1e-12)
