"""Self-crossing land-cover and reservoir polygons are repaired when read.

`docs/increments/106-wednesday-cells.md`, section 2 (106a):
`feature_input.read_class_polygons` and `feature_input.read_polygons` repair
each polygon after moving it into the mesh's CRS. A valid polygon is kept as
it is; an invalid one is replaced by the polygonal parts of
`shapely.make_valid` (default method), and dropped when it has none. Both
also return how many they repaired or dropped. `mesh` prints one stderr line
when that count is above 0 and records `land_cover_repaired` and
`reservoirs_repaired` whenever polygons are read.

THE FIXTURES, in the design's figures (a projected CRS, metres):

- the bowtie `Polygon([(0,0),(4,4),(4,0),(0,4)])`, whose ring crosses itself
  at (2, 2); `make_valid` gives two triangles of area 4 each;
- the triangle `(0,0),(4,0),(0,4)`, plan area 8; the bowtie's west half lies
  inside it and its east half touches it only along the edge (4,0)-(2,2), so
  half of the triangle is under the repaired bowtie;
- the degenerate `Polygon([(0,0),(1,1),(2,2),(0,0)])`: no area; `make_valid`
  gives a multilinestring, so nothing polygonal is left.

Bounds: shares and areas here are GEOS intersections of polygons with
integer vertices at coordinates under 10, exact in binary but compared with
1e-12 absolute so a different GEOS overlay order is not a finding.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- Both readers return a pair, `(polygons, count)`, polygons first; the
  polygons keep today's element types (`(code, geometry)` pairs and
  geometries). Whether a repaired polygon comes back as one multipolygon or
  as one entry per part is left open: the tests look at the parts.
- The stderr sentence for one repaired polygon out of one read is exactly
  `land cover: 1 of 1 polygons crossed itself and was repaired` (the
  design's example with its figures replaced), where "of N" counts the
  polygons the reader read.
- `land_cover_repaired` is a count in the `--record` JSON (a number, not
  text).
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
import shapely.affinity
from shapely.errors import GEOSException
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry

from cli_helpers import invoke
from feature_fixtures import UTM33, Feat, rules, write_geojson
from fraction_fixtures import MESH_CRS, dem_origin, write_dem, write_domain
from tin_engine.cells import overlay
from tin_engine.feature_input import _repaired, read_class_polygons, read_polygons
from tin_engine.polygon_fractions import assign

#: On a share or an area: see the module docstring.
EXACT = 1e-12
CODE = 1
BOWTIE = Polygon([(0, 0), (4, 4), (4, 0), (0, 4)])
TRIANGLE = np.array([[0.0, 0.0], [4.0, 0.0], [0.0, 4.0]])
TRIANGLES = np.array([[0, 1, 2]])
DEGENERATE = Polygon([(0, 0), (1, 1), (2, 2), (0, 0)])
VALID = Polygon([(0, 0), (3, 0), (3, 1), (0, 1)])


def write_cover(path: Path, *polygons: BaseGeometry, code: int = CODE) -> Path:
    """A GeoJSON of `polygons` coded `code` under CORINE's attribute, in UTM 33."""
    feats = [Feat(k + 1, g, {"Code_18": str(code)}) for k, g in enumerate(polygons)]
    return write_geojson(path, feats, crs=UTM33)


def pair(out: Any) -> tuple[list[Any], int]:
    """A reader's result as `(polygons, count)`, or a failure saying it is not."""
    assert isinstance(out, tuple) and len(out) == 2, (
        f"the reader returns no repair count: {type(out).__name__} {out!r}"
    )
    polygons, count = out
    assert isinstance(count, int), count
    return list(polygons), count


def classes(path: Path) -> tuple[list[tuple[int, BaseGeometry]], int]:
    return pair(read_class_polygons(path, None, rules("corine"), UTM33, TRIANGLE))


def plain(path: Path) -> tuple[list[BaseGeometry], int]:
    return pair(read_polygons(path, None, UTM33, TRIANGLE))


def parts(geometries: Sequence[BaseGeometry]) -> list[BaseGeometry]:
    return [p for g in geometries for p in shapely.get_parts(g)]


# ------------------------------------------------------------ test 1: today


class TestTheBowtieBreaksTheOverlays:
    """Test 1: unrepaired, GEOS refuses the bowtie in both overlays. This is
    why the repair is needed, and it pins that the fixture is invalid."""

    def test_the_bowtie_is_invalid(self) -> None:
        assert shapely.is_valid_reason(BOWTIE) == "Self-intersection[2 2]"

    def test_assign_throws(self) -> None:
        with pytest.raises(GEOSException, match="TopologyException: side location conflict at 2 2"):
            assign(TRIANGLE, TRIANGLES, [(CODE, BOWTIE)], known=[CODE])

    def test_overlay_throws(self) -> None:
        with pytest.raises(GEOSException, match="TopologyException: side location conflict at 2 2"):
            overlay(TRIANGLE[TRIANGLES], [BOWTIE])


# ------------------------------------------------------------ tests 2 to 5: the readers


class TestReadClassPolygons:
    def test_2_the_bowtie_comes_back_as_two_valid_triangles_counted_once(
        self, tmp_path: Path
    ) -> None:
        found, count = classes(write_cover(tmp_path / "bowtie.geojson", BOWTIE))
        assert count == 1
        assert {code for code, _ in found} == {CODE}
        pieces = parts([g for _, g in found])
        assert len(pieces) == 2, [p.wkt for p in pieces]
        assert all(p.geom_type == "Polygon" and p.is_valid for p in pieces)
        assert all(len(p.exterior.coords) == 4 for p in pieces), "each a triangle"
        assert abs(sum(p.area for p in pieces) - 8.0) <= EXACT

    def test_2_assign_on_the_repaired_bowtie_gives_half_and_half(self, tmp_path: Path) -> None:
        found, _ = classes(write_cover(tmp_path / "bowtie.geojson", BOWTIE))
        table = assign(TRIANGLE, TRIANGLES, found, known=[CODE]).table
        assert table.triangles == 1
        row = dict(
            zip(
                np.asarray(table.classes).tolist(),
                np.asarray(table.fractions).tolist(),
                strict=True,
            )
        )
        assert set(row) == {0, CODE}, row
        assert abs(row[CODE] - 0.5) <= EXACT and abs(row[0] - 0.5) <= EXACT, row

    def test_4_a_valid_polygon_comes_back_unchanged_and_counts_0(self, tmp_path: Path) -> None:
        found, count = classes(write_cover(tmp_path / "valid.geojson", VALID))
        assert count == 0
        assert len(found) == 1 and found[0][0] == CODE
        assert found[0][1].equals_exact(VALID, tolerance=0.0), found[0][1].wkt

    def test_5_a_polygon_with_no_area_is_dropped_and_counts_1(self, tmp_path: Path) -> None:
        found, count = classes(write_cover(tmp_path / "flat.geojson", DEGENERATE, VALID))
        assert count == 1
        assert [g.wkt for _, g in found] == [VALID.wkt]


class TestReadPolygons:
    def test_3_overlay_on_the_repaired_bowtie_gives_a_share_of_half(self, tmp_path: Path) -> None:
        found, count = plain(write_cover(tmp_path / "bowtie.geojson", BOWTIE))
        assert count == 1
        pieces = parts(found)
        assert len(pieces) == 2 and all(p.is_valid for p in pieces), [p.wkt for p in pieces]
        shares = overlay(TRIANGLE[TRIANGLES], found)
        assert shares.shape == (1,)
        assert abs(shares[0] - 0.5) <= EXACT, shares

    def test_4_a_valid_polygon_comes_back_unchanged_and_counts_0(self, tmp_path: Path) -> None:
        found, count = plain(write_cover(tmp_path / "valid.geojson", VALID))
        assert count == 0
        assert len(found) == 1 and found[0].equals_exact(VALID, tolerance=0.0), found

    def test_5_a_polygon_with_no_area_is_dropped_and_counts_1(self, tmp_path: Path) -> None:
        found, count = plain(write_cover(tmp_path / "flat.geojson", DEGENERATE, VALID))
        assert count == 1
        assert [g.wkt for g in found] == [VALID.wkt]


# ------------------------------------------------------------ test 6: a mesh run


def scene_bowtie() -> Polygon:
    """The bowtie scaled into `fraction_fixtures`' scene (domain about
    x0+110 .. x0+2230, y0-2230 .. y0-110): crossing at (x0+1100.3,
    y0-1100.3), no vertex on a DEM node (step 60 m)."""
    x0, y0 = dem_origin()
    lo, hi = 500.3, 1700.3
    return Polygon([(x0 + lo, y0 - hi), (x0 + hi, y0 - lo), (x0 + hi, y0 - hi), (x0 + lo, y0 - lo)])


def test_6_a_mesh_run_with_a_bowtie_land_cover_completes_and_says_so(tmp_path: Path) -> None:
    """CORINE code 311, since the run's legend must know the class."""
    bowtie = scene_bowtie()
    assert not bowtie.is_valid, "the scene's bowtie must cross itself"
    cover = write_geojson(
        tmp_path / "bowtie.geojson", [Feat(1, bowtie, {"Code_18": "311"})], crs=MESH_CRS
    )
    record = tmp_path / "x.json"
    code, output = invoke(
        "mesh",
        *("--dem", str(write_dem(tmp_path / "dem.tif"))),
        *("--domain", str(write_domain(tmp_path / "domain.geojson"))),
        *("--tolerance", "3", "--out", str(tmp_path / "x.vtk"), "--record", str(record)),
        *("--land-cover", str(cover), "--land-cover-classes", "corine"),
    )
    assert code == 0, output
    assert "TopologyException" not in output, output
    assert "land cover: 1 of 1 polygons crossed itself and was repaired" in output, output
    assert json.loads(record.read_text(encoding="ascii"))["land_cover_repaired"] == 1


# ------------------------------------------------------------ a repair cut to the mesh's box

#: The mesh points' box is widened by 100 m before the cut (`feature_input.MARGIN`),
#: so these points put the box at [0, 1000]².
NEAR = np.array([[100.0, 100.0], [900.0, 900.0], [100.0, 900.0]])
BOX = shapely.box(0, 0, 1000, 1000)
#: `@reviewer`'s polygon on 106a's clip, scaled by 100: its ring crosses
#: itself at (1276, 824), outside the box, and cutting before the repair
#: swaps inside for outside there (area 658 519.5 in the box against the
#: whole repair's 341 480.5).
RING: list[tuple[float, float]] = [(17, 4), (11, 10), (15.24, 8.24), (6.76, 8.24), (5, 4)]
RING += [(6.76, -0.24), (11, -2), (15.24, -0.24)]
CROSSING_OUTSIDE = shapely.affinity.scale(Polygon(RING), 100, 100, origin=(0, 0))


def assert_same_in_box(got: BaseGeometry | None, polygon: BaseGeometry) -> None:
    """`got`, inside BOX, is `make_valid(polygon)` inside BOX: areas within
    1e-9 relative, symmetric difference under 1e-9 of that area. Areas are
    at most 10⁶ m² here, so 1e-9 is far above GEOS's overlay rounding."""
    want = shapely.make_valid(polygon).intersection(BOX)
    have = Polygon() if got is None else got.intersection(BOX)
    assert have.area == pytest.approx(want.area, rel=1e-9), (have.area, want.area)
    assert have.symmetric_difference(want).area <= 1e-9 * want.area, have.wkt


def test_a_crossing_outside_the_box_keeps_inside_and_outside(tmp_path: Path) -> None:
    """The reader's repair, inside the mesh's box, equals the whole polygon's
    repair cut to that box."""
    assert shapely.is_valid_reason(CROSSING_OUTSIDE) == "Self-intersection[1276 824]"
    found, count = pair(
        read_polygons(write_cover(tmp_path / "p.geojson", CROSSING_OUTSIDE), None, UTM33, NEAR)
    )
    assert count == 1
    assert_same_in_box(shapely.union_all(found) if found else None, CROSSING_OUTSIDE)


def circles_with_a_bowtie(n: int, seed: int) -> list[Polygon]:
    """`n` invalid 64-gons on random circles meeting BOX, each made to cross
    itself by swapping two neighbouring vertices that both lie outside BOX."""
    rng = np.random.default_rng(seed)
    out: list[Polygon] = []
    while len(out) < n:
        c, r = rng.uniform(-500, 1500, 2), rng.uniform(200, 1500)
        a = np.sort(rng.uniform(0, 2 * np.pi, 64))
        pts = np.c_[c[0] + r * np.cos(a), c[1] + r * np.sin(a)]
        inside = shapely.contains_xy(BOX, pts[:, 0], pts[:, 1])
        both_out = np.flatnonzero(~inside[:-1] & ~inside[1:])
        if not len(both_out):
            continue
        i = both_out[rng.integers(len(both_out))]
        pts[[i, i + 1]] = pts[[i + 1, i]]
        p = Polygon(pts)
        if not p.is_valid and p.intersects(BOX):
            out.append(p)
    return out


@pytest.mark.parametrize("polygon", circles_with_a_bowtie(300, seed=106), ids=lambda p: "")
def test_random_bowties_outside_the_box_keep_inside_and_outside(polygon: Polygon) -> None:
    """Seeded: 300 circles, one small self-crossing each (11 of them fail on
    the cut-then-repair of `be26db8f`)."""
    assert_same_in_box(_repaired(polygon, NEAR), polygon)
