"""`land_cover_map.build(base, preferred, classes, region)`: the one land-cover
map, the base files joined first given wins, cut by `--prefer`, the refilled
pieces joined to their class.

`docs/increments/111-one-land-cover-map.md`, *3. Design* (the map, several
`--land-cover` files) and *5. Tests*, test 1 (the scene's map; the join test
of `docs/increments/110-prefer-lines.md`, moved here) and test 10's overlap
rule, on the polygons alone. Areas by shapely; the scene's figures are the
design's, to 1 m² (`AREA`); the unit square's corners are whole metres at
most 100 m from the origin, so its areas are exact up to rounding (`EPS`).

PINNED HERE, where the design leaves it open (listed for `@architect`):

- ``base`` is one sequence of ``(code, polygon)`` pairs per `--land-cover`
  file, in the order given; ``preferred`` is ``(replaced code, polygons)``
  pairs, as `preferred.replace`'s ``chosen``; ``classes`` is a
  `land_cover_classes.load` result; ``region`` a polygon in the map's CRS.
  Each is passed by position.
- ``LandCoverMap.areas`` is a sequence of ``(code, polygon)`` pairs;
  ``LandCoverMap.overlap`` is ``(index of the file kept, area cut)`` pairs.
- Every polygon reaching the read box is repaired whole once, where the map
  is built (`feature_input._repaired`), whether or not it is in an
  overlapping pair: a bowtie, in a single file or in either file of a pair,
  comes out valid. ``LandCoverMap.repaired`` counts the invalid polygons
  reaching the read box.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import pytest
import shapely
from shapely.geometry import box
from shapely.geometry.base import BaseGeometry

from fraction_fixtures import dem_origin
from test_cli_land_cover_polygons import cover_polygons
from test_cli_mesh_prefer import preferred_lake

#: On an area of the scene, m²: the design's figures are to 1 m².
AREA = 1.0
#: On an area of the unit-square scene, m²: whole-metre corners, rounding only.
EPS = 1e-9
HEATH, LAKE, ARABLE = 322, 512, 211


@pytest.fixture(scope="module")
def lcm() -> Any:
    from tin_engine import land_cover_map

    return land_cover_map


@pytest.fixture(scope="module")
def corine() -> Any:
    from tin_engine import land_cover_classes

    return land_cover_classes.load("corine")


def everywhere() -> BaseGeometry:
    x0, y0 = dem_origin()
    return box(x0 - 1e4, y0 - 1e4, x0 + 1e4, y0 + 1e4)


def by_code(areas: Sequence[tuple[int, BaseGeometry]]) -> dict[int, float]:
    out: dict[int, float] = {}
    for code, polygon in areas:
        out[code] = out.get(code, 0.0) + polygon.area
    return out


def planted_box() -> BaseGeometry:
    x0, y0 = dem_origin()
    return box(x0 + 420.3, y0 - 1319.8, x0 + 620.4, y0 - 1119.1)


def test_1_the_scenes_map(lcm: Any, corine: Any) -> None:
    built = lcm.build([cover_polygons()], [(LAKE, [preferred_lake()])], corine, everywhere())
    found = by_code(built.areas)
    assert set(found) == {311, 312, 512}, found
    for code, area in ((311, 1_409_620.0), (312, 509_495.0), (512, 420_090.0)):
        assert found[code] == pytest.approx(area, abs=AREA), (code, found)
    (lake,) = [g for c, g in built.areas if c == LAKE]
    assert shapely.equals(lake, preferred_lake())


def test_1_no_prefer_passes_the_base_through(lcm: Any, corine: Any) -> None:
    built = lcm.build([cover_polygons()], [], corine, everywhere())
    assert by_code(built.areas) == pytest.approx(by_code(cover_polygons()), abs=EPS)


def test_1_a_refilled_piece_joins_its_class(lcm: Any, corine: Any) -> None:
    """110's unit square: heath over the square less [40,60]², water on
    [40,60]², the preferred lake [30,50] x [40,60]. The water left, [50,60] x
    [40,60], is refilled with heath and joined to it: one heath polygon, the
    square less the lake, with no old border inside."""
    square, old, preferred = box(0, 0, 100, 100), box(40, 40, 60, 60), box(30, 40, 50, 60)
    base = [(HEATH, square.difference(old)), (LAKE, old)]
    built = lcm.build([base], [(LAKE, [preferred])], corine, square)
    heath = [g for c, g in built.areas if c == HEATH]
    assert len(heath) == 1, heath
    assert heath[0].geom_type == "Polygon", heath[0].wkt
    assert heath[0].symmetric_difference(square.difference(preferred)).area <= EPS
    (lake,) = [g for c, g in built.areas if c == LAKE]
    assert shapely.equals(lake, preferred)


def test_10_a_later_files_polygon_wholly_inside_is_cut_away(lcm: Any, corine: Any) -> None:
    polygons = dict(cover_polygons())
    a = [(311, polygons[311]), (312, polygons[312])]
    b = [(512, polygons[512]), (ARABLE, planted_box())]
    found = by_code(lcm.build([a, b], [(LAKE, [preferred_lake()])], corine, everywhere()).areas)
    assert ARABLE not in found, found
    for code, area in ((311, 1_409_620.0), (312, 509_495.0), (512, 420_090.0)):
        assert found[code] == pytest.approx(area, abs=AREA), (code, found)


def test_10_the_box_first_given_cuts_the_forest(lcm: Any, corine: Any) -> None:
    polygons = dict(cover_polygons())
    a = [(311, polygons[311]), (312, polygons[312])]
    b = [(512, polygons[512]), (ARABLE, planted_box())]
    found = by_code(lcm.build([b, a], [(LAKE, [preferred_lake()])], corine, everywhere()).areas)
    expected = {311: 1_369_460.0, ARABLE: 40_160.0, 512: 420_090.0, 312: 509_495.0}
    assert set(found) == set(expected), found
    for code, area in expected.items():
        assert found[code] == pytest.approx(area, abs=AREA), (code, found)


def test_10_two_identical_polygons_the_first_keeps_it(lcm: Any, corine: Any) -> None:
    polygons = dict(cover_polygons())
    copy = [(ARABLE, shapely.Polygon(polygons[312].exterior, polygons[312].interiors))]
    found = by_code(lcm.build([cover_polygons(), copy], [], corine, everywhere()).areas)
    assert ARABLE not in found, found
    assert found == pytest.approx(by_code(cover_polygons()), abs=EPS)


def test_10_polygons_that_only_touch_are_both_kept(lcm: Any, corine: Any) -> None:
    """Pairs that meet on a border share no area: nothing is cut."""
    left, right = box(0, 0, 50, 100), box(50, 0, 100, 100)
    built = lcm.build([[(311, left)], [(ARABLE, right)]], [], corine, box(0, 0, 100, 100))
    assert by_code(built.areas) == pytest.approx({311: 5000.0, ARABLE: 5000.0}, abs=EPS)


def bowtie(x0: float, y0: float, x1: float, y1: float) -> BaseGeometry:
    """A ring that crosses itself once, at the box's centre: its two
    triangles, west and east, are half the box ``(x0, y0, x1, y1)``."""
    return shapely.Polygon([(x0, y0), (x1, y1), (x1, y0), (x0, y1)])


def test_10_a_bowtie_in_the_later_file_is_repaired_and_cut(lcm: Any, corine: Any) -> None:
    """Code review round 1: a self-crossing polygon in an overlapping pair
    made `_joined`'s intersection throw (GEOSException). Repaired as the
    map is built (`_repaired`), the bowtie is two triangles; the
    west one, 2,500 m², lies in the first file's square, so it is the overlap
    the first file kept, and the later file keeps the east one."""
    tie = bowtie(50, 0, 150, 100)
    assert not tie.is_valid and shapely.make_valid(tie).area == pytest.approx(5000.0, abs=EPS)
    base = [[(311, box(0, 0, 100, 100))], [(ARABLE, tie)]]
    built = lcm.build(base, [], corine, box(0, 0, 200, 100))
    assert [(k, pytest.approx(g.area, abs=EPS)) for k, g in built.overlap] == [(0, 2500.0)]
    (arable,) = [g for c, g in built.areas if c == ARABLE]
    assert arable.is_valid, arable.wkt
    east = shapely.Polygon([(100, 50), (150, 0), (150, 100)])
    assert arable.symmetric_difference(east).area <= EPS, arable.wkt
    (forest,) = [g for c, g in built.areas if c == 311]
    assert shapely.equals(forest, box(0, 0, 100, 100))
    assert built.repaired == 1


def test_10_a_bowtie_in_the_earlier_file_cuts_the_later_one(lcm: Any, corine: Any) -> None:
    """The same pair the other way round: the first file's bowtie keeps its
    west triangle's 2,500 m² of the later square, which keeps the rest. Code
    review round 2: the bowtie itself comes out repaired too, its two
    triangles, and is counted once."""
    tie = bowtie(50, 0, 150, 100)
    base = [[(311, tie)], [(ARABLE, box(0, 0, 100, 100))]]
    built = lcm.build(base, [], corine, box(0, 0, 200, 100))
    assert [(k, pytest.approx(g.area, abs=EPS)) for k, g in built.overlap] == [(0, 2500.0)]
    (arable,) = [g for c, g in built.areas if c == ARABLE]
    want = box(0, 0, 100, 100).difference(shapely.Polygon([(50, 0), (100, 50), (50, 100)]))
    assert arable.symmetric_difference(want).area <= EPS, arable.wkt
    (forest,) = [g for c, g in built.areas if c == 311]
    assert forest.is_valid, forest.wkt
    assert forest.symmetric_difference(shapely.make_valid(tie)).area <= EPS, forest.wkt
    assert built.repaired == 1


def test_10_a_bowtie_in_a_single_file_is_repaired(lcm: Any, corine: Any) -> None:
    """Code review round 2: a bowtie in one file, in no overlapping pair,
    stayed unrepaired, and `preferred.replace`'s cut threw GEOSException on
    it. Repaired once where the map is built, every polygon comes out valid
    (the bowtie as its two triangles, 5,000 m²), with and without
    ``preferred``, and the one repair is counted."""
    tie, region = bowtie(50, 0, 150, 100), box(0, 0, 200, 100)
    base = [[(311, box(0, 0, 50, 100)), (ARABLE, tie)]]
    for preferred in ([], [(LAKE, [box(120, 40, 140, 60)])]):
        built = lcm.build(base, preferred, corine, region)
        assert all(g.is_valid for _, g in built.areas), [g.wkt for _, g in built.areas]
        assert built.repaired == 1, preferred
    (arable,) = [g for c, g in lcm.build(base, [], corine, region).areas if c == ARABLE]
    assert arable.symmetric_difference(shapely.make_valid(tie)).area <= EPS, arable.wkt


def test_10_clean_polygons_count_no_repair(lcm: Any, corine: Any) -> None:
    assert lcm.build([cover_polygons()], [], corine, everywhere()).repaired == 0


def test_a_map_cut_by_prefer_covers_the_region_once(lcm: Any, corine: Any) -> None:
    """Bygdin, 2026-10-29: islands in a preferred lake, where the base had
    water, were left with no class. The base covers the region: heath less a
    lake [20,80]². The preferred lake [25,75]² has two islands as holes, one
    in the base's water, one across its shore ([10,30] x [60,70], partly
    heath). After the cut, the map's areas cover the region once: they sum
    to its area and to their union's."""
    region = box(0, 0, 100, 100)
    water = box(20, 20, 80, 80)
    base = [[(HEATH, box(-200, -200, 300, 300).difference(water)), (LAKE, water)]]
    islands = [box(40, 40, 50, 50), box(10, 60, 30, 70)]
    preferred = box(25, 25, 75, 75).difference(shapely.union_all(islands))
    built = lcm.build(base, [(LAKE, [preferred])], corine, region)
    inside = [g.intersection(region) for _, g in built.areas]
    assert sum(g.area for g in inside) == pytest.approx(region.area, abs=EPS)
    assert shapely.union_all(inside).area == pytest.approx(region.area, abs=EPS)
