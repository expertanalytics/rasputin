"""`land_cover_map.coverage_gap(domain, covers, *, threads)`: the part of the
domain no land-cover or `--prefer` polygon covers.

`docs/increments/122-lakes-overrule.md`, section 3.6 (what a gap is: tiles,
parts, noise, place), section 6 (a domain no polygon reaches; a part on the
domain's border) and section 7, tests 13 to 16. Squares in metres, one
projected plane; `CoverageGap` is a `FeatureError`.

THE ORACLE: the gaps are boxes worked by hand, so their areas are exact up to
rounding (`M2` = 1e-9 m² on areas up to 2e6 m²: corners on whole or tenth
metres, below 1.2e4 or, shifted, 6.9e6, where float64 spacing is 1e-9 m).
A place is the centre of the part's largest inscribed circle, found by GEOS
to `GAP_NOISE / 10` = 0.001 m: so `PLACE` = 0.01 m, the design's bound.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- `Gap` is frozen and compares by value; `parts` is a tuple of
  `(area, x, y)` float tuples, `area` their sum; no gap is `Gap((), 0.0)`.
- `coverage_gap` takes `covers` as a sequence of shapely polygons (or
  multipolygons); an empty sequence leaves the whole domain one part.
- A strip of width 0.018 m (radius 0.009) is noise and one of 0.024 m
  (radius 0.012) a gap: the design's `GAP_NOISE` is a radius "at most"
  0.01 m, and GEOS's inscribed radius is found to within 0.001 m.
"""

from __future__ import annotations

from typing import Any

import pytest
import shapely
from shapely.geometry import Polygon, box
from shapely.geometry.base import BaseGeometry

#: On an area, m²: see the module docstring.
M2 = 1e-9
#: On a place, m: the inscribed circle's tolerance is 0.001 m.
PLACE = 0.01
HEATH, RIVER, LAKE = 322, 511, 512


@pytest.fixture(scope="module")
def lcm() -> Any:
    from tin_engine import land_cover_map

    return land_cover_map


def shifted(geometry: BaseGeometry, origin: tuple[float, float]) -> BaseGeometry:
    return shapely.affinity.translate(geometry, *origin)


ORIGINS = [(0.0, 0.0), (136000.0, 6830000.0)]
ORIGIN_IDS = ["origin", "utm-scale"]


def test_the_chosen_constants(lcm: Any) -> None:
    assert lcm.GAP_NOISE == 0.01
    assert lcm.GAP_TILE == 5000.0


def test_coverage_gap_is_a_feature_error(lcm: Any) -> None:
    from tin_engine.feature_input import FeatureError

    assert issubclass(lcm.CoverageGap, FeatureError)


# ------------------------------------------------------------ test 13


@pytest.mark.parametrize("origin", ORIGINS, ids=ORIGIN_IDS)
def test_13_one_part_on_the_domains_corner(lcm: Any, origin: tuple[float, float]) -> None:
    """Domain [0,100]²; covers [0,100] x [0,60] and [0,60] x [60,100]: the
    gap [60,100]², 1,600 m², inscribed radius 20 at (80, 80). It touches the
    domain's border and is still a gap (section 6)."""
    x0, y0 = origin
    covers = [shifted(box(0, 0, 100, 60), origin), shifted(box(0, 60, 60, 100), origin)]
    gap = lcm.coverage_gap(shifted(box(0, 0, 100, 100), origin), covers)
    assert isinstance(gap, lcm.Gap)
    assert len(gap.parts) == 1, gap
    ((area, x, y),) = gap.parts
    assert area == pytest.approx(1600.0, abs=M2)
    assert gap.area == pytest.approx(1600.0, abs=M2)
    assert (x - x0, y - y0) == pytest.approx((80.0, 80.0), abs=PLACE)


def test_13_parts_come_area_descending(lcm: Any) -> None:
    """Domain [0,100]²; covers all but [0,10]² (100 m²) and [50,80]² (900
    m²): the larger first; `area` the sum."""
    covers = [box(0, 0, 100, 100).difference(box(0, 0, 10, 10)).difference(box(50, 50, 80, 80))]
    gap = lcm.coverage_gap(box(0, 0, 100, 100), covers)
    assert [round(a, 6) for a, _, _ in gap.parts] == [900.0, 100.0]
    assert gap.area == pytest.approx(1000.0, abs=M2)
    assert gap.parts[0][1:] == pytest.approx((65.0, 65.0), abs=PLACE)
    assert gap.parts[1][1:] == pytest.approx((5.0, 5.0), abs=PLACE)


def test_a_domain_no_polygon_reaches_is_one_part_the_whole_domain(lcm: Any) -> None:
    """Section 6: a cover far off and none at all both leave the domain."""
    domain = Polygon([(0, 0), (40, 0), (40, 30), (0, 30)])
    for covers in ([box(500, 500, 600, 600)], []):
        gap = lcm.coverage_gap(domain, covers)
        assert len(gap.parts) == 1, (covers, gap)
        assert gap.area == pytest.approx(1200.0, abs=M2)
        assert gap.parts[0][1:] == pytest.approx((20.0, 15.0), abs=PLACE)


def test_a_covered_domain_has_no_gap(lcm: Any) -> None:
    gap = lcm.coverage_gap(box(0, 0, 100, 100), [box(-10, -10, 110, 110)])
    assert gap == lcm.Gap((), 0.0)


# ------------------------------------------------------------ test 14


def test_14_an_island_left_out_of_the_map_is_no_gap(lcm: Any) -> None:
    """107's case 5: heath on [0,100]² less [30,70]², lake [30,70]², the
    `--prefer` ring [30,70]² less [45,55]². The island [45,55]² is left out
    of the map (test 8), but the base's lake covers it: no part."""
    ring = box(30, 30, 70, 70).difference(box(45, 45, 55, 55))
    covers = [box(0, 0, 100, 100).difference(box(30, 30, 70, 70)), box(30, 30, 70, 70), ring]
    gap = lcm.coverage_gap(box(0, 0, 100, 100), covers)
    assert gap.parts == () and gap.area == 0.0


def test_14_test_7s_scene_is_no_gap(lcm: Any) -> None:
    """Heath over [0,100]² less the lake [20,80]² and the river [80,90] x
    [10,90]; NVE's lake [20,78] x [20,80]: every point is covered."""
    lake, river = box(20, 20, 80, 80), box(80, 10, 90, 90)
    heath = box(0, 0, 100, 100).difference(lake).difference(river)
    gap = lcm.coverage_gap(box(0, 0, 100, 100), [heath, lake, river, box(20, 20, 78, 80)])
    assert gap.parts == () and gap.area == 0.0


# ------------------------------------------------------------ test 15


@pytest.mark.parametrize("origin", ORIGINS, ids=ORIGIN_IDS)
def test_15_a_slit_of_a_millimetre_is_noise(lcm: Any, origin: tuple[float, float]) -> None:
    """Covers [0,50] x [0,100] and [50.001,100] x [0,100]: a slit of 0.1 m²,
    inscribed radius 0.0005 m, under `GAP_NOISE`: no part."""
    covers = [shifted(box(0, 0, 50, 100), origin), shifted(box(50.001, 0, 100, 100), origin)]
    gap = lcm.coverage_gap(shifted(box(0, 0, 100, 100), origin), covers)
    assert gap.parts == () and gap.area == 0.0


@pytest.mark.parametrize("origin", ORIGINS, ids=ORIGIN_IDS)
def test_15_a_slit_of_ten_centimetres_is_a_gap(lcm: Any, origin: tuple[float, float]) -> None:
    """[50.1,100] instead: 10 m², radius 0.05 m: one part, though thin."""
    covers = [shifted(box(0, 0, 50, 100), origin), shifted(box(50.1, 0, 100, 100), origin)]
    gap = lcm.coverage_gap(shifted(box(0, 0, 100, 100), origin), covers)
    assert len(gap.parts) == 1, gap
    assert gap.parts[0][0] == pytest.approx(10.0, abs=1e-6 if origin[0] else M2)
    assert gap.parts[0][1] - origin[0] == pytest.approx(50.05, abs=PLACE)


@pytest.mark.parametrize(("width", "parts"), [(0.018, 0), (0.024, 1)])
def test_15_noise_is_a_width_not_an_area(lcm: Any, width: float, parts: int) -> None:
    """A strip 1 km long: 18 m² at width 0.018 m (radius 0.009, noise), 24 m²
    at 0.024 m (radius 0.012, a gap). The area does not decide."""
    domain = box(0, 0, 1000, 100)
    covers = [box(0, 0, 1000, 50), box(0, 50 + width, 1000, 100)]
    gap = lcm.coverage_gap(domain, covers)
    assert len(gap.parts) == parts, gap


# ------------------------------------------------------------ test 16


def tiles_scene() -> tuple[Polygon, list[BaseGeometry]]:
    return box(0, 0, 12000, 1000), [box(0, 0, 4000, 1000), box(6000, 0, 12000, 1000)]


def test_16_a_part_cut_by_a_tile_border_is_one_part(lcm: Any) -> None:
    """Domain [0,12000] x [0,1000]; covers [0,4000] and [6000,12000] across:
    the gap [4000,6000] x [0,1000] is cut by the tile border x = 5000 and
    comes back one part, 2,000,000 m², radius 500 (its centre anywhere on
    y = 500 between x = 4500 and 5500)."""
    domain, covers = tiles_scene()
    gap = lcm.coverage_gap(domain, covers)
    assert len(gap.parts) == 1, gap
    ((area, x, y),) = gap.parts
    assert area == pytest.approx(2_000_000.0, abs=1e-6)
    assert y == pytest.approx(500.0, abs=PLACE)
    assert 4500.0 - PLACE <= x <= 5500.0 + PLACE


def test_16_threads_and_the_covers_order_change_nothing(lcm: Any) -> None:
    domain, covers = tiles_scene()
    one = lcm.coverage_gap(domain, covers, threads=1)
    assert lcm.coverage_gap(domain, covers, threads=4) == one
    assert lcm.coverage_gap(domain, covers[::-1], threads=1) == one


def test_16_threads_on_many_tiles_with_many_parts(lcm: Any) -> None:
    """A 30 km x 20 km domain (24 tiles) with a 10 m hole every 3 km: 70 parts
    of 100 m² each, equal in 1 and 8 threads."""
    domain = box(0, 0, 30000, 20000)
    holes = [
        box(x, y, x + 10, y + 10)
        for x in range(1000, 30000, 3000)
        for y in range(1000, 20000, 3000)
    ]
    covers = [domain.difference(shapely.union_all(holes))]
    one = lcm.coverage_gap(domain, covers, threads=1)
    assert len(one.parts) == len(holes) == 70
    assert one.area == pytest.approx(7000.0, abs=1e-6)
    assert lcm.coverage_gap(domain, covers, threads=8) == one
