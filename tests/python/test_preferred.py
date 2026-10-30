"""`tin_engine.preferred.replace`: a named class taken from preferred polygons,
the rest from the base land cover.

`docs/increments/107-prefer-classes.md`, section 2 (the rule) and section 5,
cases 1 to 5 (synthetic squares, metres, one CRS), on `replace` with `near`
the square [0,100]² (or the part a case names). `overlay`, which also
measured the areas inside a mesh's triangles, goes
(`docs/increments/111-one-land-cover-map.md`, section 8), and with it the
checks of `Overlaid.areas`; the class areas below are `replace`'s. Codes are
CORINE's: 322 heath, 312 conifer forest, 511 water courses, 512 water bodies
(the lake class), 335 glaciers.

THE ORACLE: class areas summed from the returned polygons by GEOS. Every box
here has corners on whole metres, so GEOS's areas are exact up to rounding;
the bound `M2` is 1e-6 m² at 100 m scale (areas up to 1e4 m²), the largest
input used.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- `Replaced.polygons` is a sequence of (code, geometry): kept, given,
  refilled.
- A piece that takes a class may come back as its own polygon or merged with
  a neighbour: only the summed area per code is checked.
- The returned polygons never overlap (their areas sum to their union's).
- The same replaced code twice in `chosen` unions the two sources.
- A preferred polygon that meets no geometry of `near` is not returned; one
  that meets one comes back whole (unclipped).
"""

from __future__ import annotations

import importlib
from collections.abc import Sequence
from types import ModuleType

import pytest
import shapely
from shapely.geometry import box
from shapely.geometry.base import BaseGeometry

#: On an area in m²: see the module docstring.
M2 = 1e-6

HEATH, CONIFER, RIVER, LAKE, GLACIER = 322, 312, 511, 512, 335


@pytest.fixture
def preferred() -> ModuleType:
    return importlib.import_module("tin_engine.preferred")


def square() -> list[BaseGeometry]:
    """`near`: the square [0,100]²."""
    return [box(0, 0, 100, 100)]


def per_code(polygons: Sequence[tuple[int, BaseGeometry]]) -> dict[int, float]:
    areas: dict[int, float] = {}
    for code, geometry in polygons:
        areas[code] = areas.get(code, 0.0) + float(geometry.area)
    return areas


def check(found: dict[int, float], expected: dict[int, float]) -> None:
    found = {c: a for c, a in found.items() if a > M2}
    assert found.keys() == expected.keys(), found
    for code, area in expected.items():
        assert found[code] == pytest.approx(area, abs=M2), (code, found)


def check_disjoint(polygons: Sequence[tuple[int, BaseGeometry]]) -> None:
    union = shapely.union_all([g for _, g in polygons]).area
    assert sum(g.area for _, g in polygons) == pytest.approx(union, abs=M2)


def heath_with(hole: BaseGeometry, code: int = LAKE) -> list[tuple[int, BaseGeometry]]:
    """Heath over the square less `hole`, and `hole` as `code`."""
    return [(HEATH, box(0, 0, 100, 100).difference(hole)), (code, hole)]


def test_case1_a_leftover_piece_takes_its_neighbours_class(preferred: ModuleType) -> None:
    """Lake [40,60]² in the base, preferred lake [30,50] x [40,60]: the
    preferred lake is the lake (400 m²), the leftover [50,60] x [40,60] borders
    heath on three sides and becomes heath; heath totals 9,600 m². 400 given
    replaces the 400 m² lake."""
    base = heath_with(box(40, 40, 60, 60))
    result = preferred.replace(base, [(LAKE, [box(30, 40, 50, 60)])], square())
    check(per_code(result.polygons), {LAKE: 400.0, HEATH: 9600.0})
    check_disjoint(result.polygons)


def test_case2_no_preferred_lake_the_longest_border_wins(preferred: ModuleType) -> None:
    """Heath left of x = 50, conifer right; lake [44,60] x [40,60] across the
    line, no preferred lake. The lake borders heath for 32 m and conifer for
    40 m, so all 320 m² become conifer."""
    lake = box(44, 40, 60, 60)
    base = [
        (HEATH, box(0, 0, 50, 100).difference(lake)),
        (CONIFER, box(50, 0, 100, 100).difference(lake)),
        (LAKE, lake),
    ]
    result = preferred.replace(base, [(LAKE, [])], square())
    check(per_code(result.polygons), {HEATH: 5000.0 - 120.0, CONIFER: 5000.0 - 200.0 + 320.0})
    check_disjoint(result.polygons)


def test_case2_a_tie_takes_the_lower_code(preferred: ModuleType) -> None:
    """Lake [40,60]² centred on the line: 40 m of border each side; conifer
    (312) is the lower code."""
    lake = box(40, 40, 60, 60)
    base = [
        (HEATH, box(0, 0, 50, 100).difference(lake)),
        (CONIFER, box(50, 0, 100, 100).difference(lake)),
        (LAKE, lake),
    ]
    result = preferred.replace(base, [(LAKE, [])], square())
    check(per_code(result.polygons), {HEATH: 4800.0, CONIFER: 5200.0})


def test_case2_each_connected_piece_takes_its_own_class(preferred: ModuleType) -> None:
    """Lake [30,70] x [40,60] across the line, preferred lake [45,55] x [40,60]
    splits the leftover in two: the west piece borders heath only, the east
    piece conifer only, 300 m² each."""
    lake = box(30, 40, 70, 60)
    base = [
        (HEATH, box(0, 0, 50, 100).difference(lake)),
        (CONIFER, box(50, 0, 100, 100).difference(lake)),
        (LAKE, lake),
    ]
    result = preferred.replace(base, [(LAKE, [box(45, 40, 55, 60)])], square())
    check(per_code(result.polygons), {HEATH: 4900.0, CONIFER: 4900.0, LAKE: 200.0})
    check_disjoint(result.polygons)


@pytest.mark.parametrize(
    ("order", "lake", "glacier"),
    [((LAKE, GLACIER), 400.0, 300.0), ((GLACIER, LAKE), 300.0, 400.0)],
    ids=["lake-first", "glacier-first"],
)
def test_case3_the_first_preferred_source_wins_an_overlap(
    preferred: ModuleType, order: tuple[int, int], lake: float, glacier: float
) -> None:
    """Preferred lake [0,20]², glacier [10,30]², over heath."""
    sources = {LAKE: [box(0, 0, 20, 20)], GLACIER: [box(10, 10, 30, 30)]}
    chosen = [(code, sources[code]) for code in order]
    result = preferred.replace([(HEATH, box(0, 0, 100, 100))], chosen, square())
    check(per_code(result.polygons), {LAKE: lake, GLACIER: glacier, HEATH: 9300.0})
    check_disjoint(result.polygons)


def test_case3_the_same_code_twice_unions_the_sources(preferred: ModuleType) -> None:
    """Two lake sources, [0,20]² and [10,30]²: one lake of 700 m²."""
    chosen = [(LAKE, [box(0, 0, 20, 20)]), (LAKE, [box(10, 10, 30, 30)])]
    result = preferred.replace([(HEATH, box(0, 0, 100, 100))], chosen, square())
    check(per_code(result.polygons), {LAKE: 700.0, HEATH: 9300.0})
    check_disjoint(result.polygons)


def test_case4_other_water_classes_stay_the_bases(preferred: ModuleType) -> None:
    """A river strip [0,100] x [45,55] (511) and no 512 in the base; preferred
    lake [40,60]²: lake 400, river 800, heath 8,800; nothing replaced."""
    river = box(0, 45, 100, 55)
    base = heath_with(river, RIVER)
    result = preferred.replace(base, [(LAKE, [box(40, 40, 60, 60)])], square())
    check(per_code(result.polygons), {LAKE: 400.0, RIVER: 800.0, HEATH: 8800.0})
    check_disjoint(result.polygons)


def test_case5_a_piece_in_a_preferred_hole_takes_the_nearest_kept_class(
    preferred: ModuleType,
) -> None:
    """Lake [40,60]² in the base; preferred lake the ring [30,70]² less
    [45,55]²: the 100 m² piece [45,55]², in the ring's hole (an island),
    borders only the preferred lake, so it takes the class of the nearest
    kept polygon (heath), and the polygons cover the square once (Bygdin's
    islands, 2026-10-29)."""
    ring = box(30, 30, 70, 70).difference(box(45, 45, 55, 55))
    base = heath_with(box(40, 40, 60, 60))
    result = preferred.replace(base, [(LAKE, [ring])], square())
    check(per_code(result.polygons), {LAKE: 1500.0, HEATH: 8500.0})
    check_disjoint(result.polygons)
    covered = shapely.union_all([g for _, g in result.polygons])
    assert covered.area == pytest.approx(10000.0, abs=M2)


def test_case5_an_island_takes_the_nearer_of_two_classes(preferred: ModuleType) -> None:
    """Heath west of x = 50, conifer east; lake [20,80]²; preferred lake
    [22,80] x [20,80] less the island [26,36] x [40,50]: the island is 6 m
    from heath across the lake and 44 m from conifer, so it is heath."""
    west, east, water = box(0, 0, 50, 100), box(50, 0, 100, 100), box(20, 20, 80, 80)
    base = [(HEATH, west.difference(water)), (CONIFER, east.difference(water)), (LAKE, water)]
    island = box(26, 40, 36, 50)
    result = preferred.replace(base, [(LAKE, [box(22, 20, 80, 80).difference(island)])], square())
    (got,) = [n for _, p, n in result.refilled if p.equals(island)]
    assert got == HEATH


def test_case5_a_piece_outside_every_preferred_polygon_is_left_out(
    preferred: ModuleType,
) -> None:
    """Lake [40,60]² alone in the square (no heath); preferred lake [30,50] x
    [40,60]: the piece [50,60] x [40,60] borders no kept polygon and lies in
    no hole, so no polygon covers it (class 0)."""
    base = [(LAKE, box(40, 40, 60, 60))]
    result = preferred.replace(base, [(LAKE, [box(30, 40, 50, 60)])], square())
    check(per_code(result.polygons), {LAKE: 400.0})
    covered = shapely.union_all([g for _, g in result.polygons])
    assert covered.intersection(box(50, 40, 60, 60)).area == pytest.approx(0.0, abs=M2)


def test_a_preferred_polygon_off_the_mesh_is_dropped_one_across_its_edge_kept_whole(
    preferred: ModuleType,
) -> None:
    """Lake [200,220]² meets not the square: not returned. Lake [90,110] x
    [40,60] crosses the square's east side: returned whole (400 m²)."""
    chosen = [(LAKE, [box(200, 200, 220, 220), box(90, 40, 110, 60)])]
    result = preferred.replace([(HEATH, box(0, 0, 100, 100))], chosen, square())
    check(per_code(result.polygons), {LAKE: 400.0, HEATH: 9800.0})


def test_replace_gives_kept_then_given_then_refilled(preferred: ModuleType) -> None:
    """Case 1's polygons in `Replaced.polygons`' order: the kept heath, the
    preferred lake, the refilled piece [50,60] x [40,60] as heath."""
    base = heath_with(box(40, 40, 60, 60))
    chosen = [(LAKE, [box(30, 40, 50, 60)])]
    replaced = preferred.replace(base, chosen, near=square())
    assert [c for c, _ in replaced.polygons] == [HEATH, LAKE, HEATH]
    lake = box(30, 40, 50, 60)
    want = [box(0, 0, 100, 100).difference(box(40, 40, 60, 60)).difference(lake), lake,
            box(50, 40, 60, 60)]  # fmt: skip
    for (_, got), expected in zip(replaced.polygons, want, strict=True):
        assert shapely.equals(got, expected), got.wkt
    check(per_code(replaced.polygons), {LAKE: 400.0, HEATH: 9600.0})
