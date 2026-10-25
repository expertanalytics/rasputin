"""The cell table: legend categories, `tin_engine.cells` and the two GeoPackage writers.

`docs/increments/105-cell-table.md`, section 3 (the columns), 4.1 (categories
in the legend files, the check only for the cells table) and section 6, T1
to T4. Every expected figure is hand arithmetic or recomputed here from the
formulas of section 3, never read from the producer.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- `cells.categorise(fractions, categories)` returns a (triangles, 5) float64
  array, columns in the blueprint's order: glacier, lake, forest, savanna,
  unspecified. A class absent from `categories` is unspecified.
- `cells.geometry(vertices, triangles)` returns `(midpoints (T, 3), area_m2,
  slope_deg, aspect_deg)`, degrees, aspect NaN where `n_x = n_y = 0`.
- `cells.cell_table(vertices, triangles, fractions, categories)` builds the
  `CellTable`; `geopackage_write.write_cells(conn, cells, srs_id)` and
  `write_fractions(conn, fractions)` write the two tables.
- The `cells` columns in section 3's order; `land_cover_fractions.triangle`
  is the triangle's `fid` (index + 1).
- `categories` holds only named classes: 0 and a class without `category`
  are absent from it.
"""

from __future__ import annotations

import math
import re
import sqlite3
import tomllib
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from pyproj import CRS
from pyproj.enums import WktVersion

from gpkg_fixtures import needs_rtree
from legend_fixtures import coded_map, corine_like, write
from tin_engine import cells, legends
from tin_engine.io import geopackage_write
from tin_engine.io.geopackage import layer_info, query_features
from tin_engine.io.models import ClassFractions, SrsEntry

GLACIER, LAKE, FOREST, SAVANNA, UNSPECIFIED = range(5)
#: Section 3: 44 classes at most of float32 rounding (6e-8 each) in one row.
FIVE_SUM = 3e-6
SHIPPED = Path(legends.__file__).parent / "data" / "legends"
CORINE = {311: "forest", 312: "forest", 313: "forest", 511: "lake", 512: "lake",
          521: "lake", 522: "lake", 335: "glacier", 523: "unspecified"}  # fmt: skip
WANT = {
    "corine": CORINE,
    "corine-water": CORINE,
    "clc18_kode": CORINE,
    "ar50": {30: "forest", 3031: "forest", 3032: "forest", 3033: "forest", 3099: "forest",
             81: "lake", 70: "glacier", 82: "unspecified"},
    "mapbiomas-c11": {3: "forest", 5: "forest", 6: "forest", 49: "forest", 9: "forest",
                      33: "lake", 31: "lake", 4: "savanna", 7: "savanna"},
}  # fmt: skip


def table(offsets: list[int], classes: list[int], fractions: list[float]) -> ClassFractions:
    return ClassFractions(
        offsets=np.array(offsets, np.uint32),
        classes=np.array(classes, np.int32),
        fractions=np.array(fractions, np.float32),
    )


def five(fractions: ClassFractions, legend: str) -> np.ndarray:
    got = cells.categorise(fractions, legends.load(legend).categories)
    assert got.dtype == np.float64 and got.shape == (fractions.triangles, 5), got
    return got


# ------------------------------------------------------------ T1 legends


@pytest.mark.parametrize("name", sorted(WANT))
def test_t1_shipped_categories(name: str) -> None:
    legend = legends.load(name)
    assert dict(legend.categories) == WANT[name]
    assert dict(cells.categories_of(legend)) == WANT[name]


@pytest.mark.parametrize("name", sorted(WANT))
def test_t1_water_codes_are_the_files_water_classes(name: str) -> None:
    data = tomllib.loads((SHIPPED / f"{name}.toml").read_text("utf-8"))
    water = {c["code"] for c in data["class"] if c.get("water")}
    assert legends.load(name).water_codes == frozenset(water)
    assert water <= set(WANT[name]), "every shipped water class is named"


def test_t1_class_0_is_unspecified() -> None:
    assert legends.load("corine").categories.get(0, "unspecified") == "unspecified"


def test_t1_uncategorised_water_loads_and_categories_of_refuses_it(tmp_path: Path) -> None:
    legend = legends.load(str(write(tmp_path / "w.toml", corine_like("w", "40"))))
    assert 512 in legend.water_codes and 512 not in legend.categories
    with pytest.raises(ValueError, match=r"\b512\b"):
        cells.categories_of(legend)


def test_t1_categories_of_names_every_uncategorised_water_class(tmp_path: Path) -> None:
    two = (
        '[[class]]\ncode = 2\nlabel = "Lake"\ncolour = "#2532e4"\nwater = true\n\n'
        '[[class]]\ncode = 7\nlabel = "Pond"\ncolour = "#2532e5"\nwater = true\n'
    )
    legend = legends.load(str(write(tmp_path / "two.toml", coded_map(classes=two))))
    with pytest.raises(ValueError) as caught:
        cells.categories_of(legend)
    assert re.search(r"\b2\b", str(caught.value)) and re.search(r"\b7\b", str(caught.value))


def test_t1_land_without_category_passes_the_check(tmp_path: Path) -> None:
    """Only water must be named: Forest (code 1) has no category and passes."""
    named = (
        '[[class]]\ncode = 1\nlabel = "Forest"\ncolour = "#1f8d49"\n\n'
        '[[class]]\ncode = 2\nlabel = "Lake"\ncolour = "#2532e4"\nwater = true\n'
        'category = "lake"\n'
    )
    legend = legends.load(str(write(tmp_path / "n.toml", coded_map(classes=named))))
    assert dict(cells.categories_of(legend)) == {2: "lake"}


def test_t1_an_unknown_category_is_refused_at_load(tmp_path: Path) -> None:
    sea = '[[class]]\ncode = 2\nlabel = "Sea"\ncolour = "#2532e4"\nwater = true\ncategory = "sea"\n'
    path = write(tmp_path / "sea.toml", coded_map(classes=sea))
    with pytest.raises(ValueError, match="category"):
        legends.load(str(path))


# ------------------------------------------------------------ T2 categorise


def test_t2_corine_rows() -> None:
    got = five(table([0, 2, 3, 5], [311, 512, 335, 0, 312], [0.25, 0.75, 1.0, 0.5, 0.5]), "corine")
    assert_array_equal(got[:, FOREST], [0.25, 0, 0.5])
    assert_array_equal(got[:, LAKE], [0.75, 0, 0])
    assert_array_equal(got[:, GLACIER], [0, 1, 0])
    assert_array_equal(got[:, SAVANNA], [0, 0, 0])
    assert_array_equal(got[:, UNSPECIFIED], [0, 0, 0.5])
    assert_array_equal(got.sum(1) - got[:, SAVANNA], [1, 1, 1])


def test_t2_sea_is_unspecified() -> None:
    got = five(table([0, 2], [311, 523], [0.25, 0.75]), "corine")
    assert_array_equal(got[0], [0, 0, 0.25, 0, 0.75])


def test_t2_savanna_sits_inside_unspecified() -> None:
    got = five(table([0, 2], [3, 4], [0.25, 0.75]), "mapbiomas-c11")
    assert_array_equal(got[0], [0, 0, 0.25, 0.75, 0.75])


def test_t2_unspecified_is_one_minus_the_sum_not_the_unnamed_classes() -> None:
    """A row short of 1 (forest 0.5 alone) and an empty row: Shyft's 1 - sum."""
    got = five(table([0, 1, 1], [311], [0.5]), "corine")
    assert_array_equal(got[:, UNSPECIFIED], [0.5, 1.0])


def test_t2_unspecified_is_clamped_at_0() -> None:
    fractions = np.array([0.4, 0.6], np.float32)
    assert fractions.astype(np.float64).sum() > 1.0, "the premise: float32 rounding overshoots"
    got = five(table([0, 2], [311, 512], fractions.tolist()), "corine")
    assert got[0, UNSPECIFIED] == 0.0
    assert abs(got[0].sum() - got[0, SAVANNA] - 1.0) <= FIVE_SUM


def test_t2_no_triangles() -> None:
    assert five(table([0], [], []), "corine").shape == (0, 5)


# ------------------------------------------------------------ T3 geometry

#: F1 of `test_io_geopackage_write.py`.
F1_VERTICES = np.array([[0, 0, 0], [1, 0, 1], [1, 1, 2], [0, 1, 3]], np.float64)
F1_TRIANGLES = np.array([[0, 1, 2], [0, 2, 3]], np.int64)


def one(corners: list[list[float]]) -> tuple[np.ndarray, ...]:
    return cells.geometry(np.array(corners, np.float64), np.array([[0, 1, 2]], np.int64))


def test_t3_f1() -> None:
    mid, area, slope, aspect = cells.geometry(F1_VERTICES, F1_TRIANGLES)
    assert_allclose(mid, [[2 / 3, 1 / 3, 1], [1 / 3, 2 / 3, 5 / 3]], rtol=1e-15)
    assert_allclose(area, [0.5, 0.5], rtol=1e-15)
    want_slope = [math.degrees(math.atan(math.sqrt(2))), math.degrees(math.atan(math.sqrt(10)))]
    assert_allclose(slope, want_slope, rtol=1e-12)
    assert_allclose(aspect, [225.0, math.degrees(math.atan2(1, -3))], rtol=1e-12)


@pytest.mark.parametrize(
    "corners",
    [[[0, 0, 0], [1, 0, 0], [0, 1, 1]], [[0, 0, 0], [0, 1, 1], [1, 0, 0]]],
    ids=["counter-clockwise", "clockwise"],
)
def test_t3_south_facing_either_winding(corners: list[list[float]]) -> None:
    _, area, slope, aspect = one(corners)
    assert_allclose((area[0], slope[0], aspect[0]), (0.5, 45.0, 180.0), rtol=1e-12)


def test_t3_flat_has_no_aspect() -> None:
    _, _, slope, aspect = one([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    assert slope[0] == 0.0 and np.isnan(aspect[0])


def test_t3_duplicate_vertex_has_no_area_and_no_aspect() -> None:
    _, area, slope, aspect = one([[0, 0, 0], [0, 0, 0], [1, 1, 1]])
    assert area[0] == 0.0 and slope[0] == 0.0 and np.isnan(aspect[0])


def test_t3_utm_scale_coordinates() -> None:
    """F1 moved by (500 km, 7000 km): slope and aspect unchanged to 1e-9."""
    near = cells.geometry(F1_VERTICES, F1_TRIANGLES)
    far = cells.geometry(F1_VERTICES + np.array([5e5, 7e6, 0]), F1_TRIANGLES)
    assert_allclose(far[0], near[0] + np.array([5e5, 7e6, 0]), rtol=1e-15)
    for a, b in zip(near[1:], far[1:], strict=True):
        assert_allclose(b, a, rtol=1e-9)


# ------------------------------------------------------------ T4 writer

#: F1 and a flat third triangle, whose aspect is NULL.
VERTICES = np.vstack([F1_VERTICES, [[2, 0, 5], [3, 0, 5], [2, 1, 5]]])
TRIANGLES = np.vstack([F1_TRIANGLES, [[4, 5, 6]]])
FRACTIONS = table([0, 2, 3, 4], [311, 512, 335, 0], [0.25, 0.75, 1.0, 1.0])
COLUMNS = [
    ("fid", "INTEGER"), ("geom", "POINT"), ("catchment_id", "INTEGER"), ("area_m2", "REAL"),
    ("z_m", "REAL"), ("slope_deg", "REAL"), ("aspect_deg", "REAL"), ("glacier", "REAL"),
    ("lake", "REAL"), ("reservoir", "REAL"), ("forest", "REAL"), ("unspecified", "REAL"),
    ("savanna", "REAL"),
]  # fmt: skip
UTM33 = 25833


@pytest.fixture
def conn() -> Iterator[sqlite3.Connection]:
    srs = SrsEntry(
        srs_id=UTM33, organization="EPSG", code=UTM33, name="ETRS89 / UTM zone 33N",
        definition=CRS.from_epsg(UTM33).to_wkt(version=WktVersion.WKT1_GDAL),
    )  # fmt: skip
    connection = sqlite3.connect(":memory:")
    geopackage_write.create_geopackage(connection, srs)
    categories = legends.load("corine").categories
    table_ = cells.cell_table(VERTICES, TRIANGLES, FRACTIONS, categories)
    geopackage_write.write_cells(connection, table_, UTM33)
    geopackage_write.write_fractions(connection, FRACTIONS)
    yield connection
    connection.close()


def test_t4_cells_columns_and_registration(conn: sqlite3.Connection) -> None:
    info = conn.execute("PRAGMA table_info(cells)").fetchall()
    assert [(r[1], r[2]) for r in info] == COLUMNS
    assert conn.execute(
        "SELECT geometry_type_name, srs_id, z FROM gpkg_geometry_columns WHERE table_name = 'cells'"
    ).fetchone() == ("POINT", UTM33, 1)
    kinds = dict(conn.execute("SELECT table_name, data_type FROM gpkg_contents").fetchall())
    assert (kinds["cells"], kinds["land_cover_fractions"]) == ("features", "attributes")


@needs_rtree
def test_t4_cells_rtree_has_one_point_box_per_triangle(conn: sqlite3.Connection) -> None:
    # SQLite's R-tree stores 32-bit floats; when the cast lands on the wrong
    # side it scales by 1 -+ 2**-23 and casts again, so each box contains its
    # midpoint within 1.5 * 2**-23 relative; 2**-22 bounds it, at any scale.
    rows = conn.execute("SELECT id, minx, maxx, miny, maxy FROM rtree_cells_geom ORDER BY id")
    boxes = np.array([r[1:] for r in rows])
    mid, *_ = cells.geometry(VERTICES, TRIANGLES)
    assert (boxes[:, [0, 2]] <= mid[:, :2]).all() and (mid[:, :2] <= boxes[:, [1, 3]]).all()
    assert_allclose(boxes, mid[:, [0, 0, 1, 1]], rtol=2.0**-22, atol=0)
    assert conn.execute("SELECT count(*) FROM rtree_cells_geom").fetchone() == (3,)


def test_t4_cells_rows(conn: sqlite3.Connection) -> None:
    rows = conn.execute(
        "SELECT fid, catchment_id, area_m2, z_m, slope_deg, aspect_deg, glacier, lake, "
        "reservoir, forest, unspecified, savanna FROM cells ORDER BY fid"
    ).fetchall()
    assert [r[:2] for r in rows] == [(1, 1), (2, 1), (3, 1)]
    mid, area, slope, aspect = cells.geometry(VERTICES, TRIANGLES)
    assert_allclose([r[2:5] for r in rows], np.column_stack([area, mid[:, 2], slope]), rtol=1e-15)
    assert rows[2][5] is None, "flat: aspect NULL"
    assert_allclose([r[5] for r in rows[:2]], aspect[:2], rtol=1e-15)
    assert [r[6:] for r in rows] == [
        (0.0, 0.75, 0.0, 0.25, 0.0, 0.0),
        (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
    ]


def test_t4_cells_read_back_as_points_in_fid_order(conn: sqlite3.Connection) -> None:
    layer = layer_info(conn, "cells")
    found = list(query_features(conn, layer, None, "z_m"))
    mid, *_ = cells.geometry(VERTICES, TRIANGLES)
    assert [f.fid for f in found] == [1, 2, 3]
    assert all(f.geometry is not None and f.geometry.geom_type == "Point" for f in found)
    assert_allclose([f.geometry.coords[0] for f in found if f.geometry], mid, rtol=1e-15)


def test_t4_fractions_long_form(conn: sqlite3.Connection) -> None:
    rows = conn.execute(
        "SELECT triangle, code, fraction FROM land_cover_fractions ORDER BY fid"
    ).fetchall()
    assert rows == [(1, 311, 0.25), (1, 512, 0.75), (2, 335, 1.0), (3, 0, 1.0)]
    indexed = [
        [c[2] for c in conn.execute(f"PRAGMA index_info('{i[1]}')")]
        for i in conn.execute("PRAGMA index_list(land_cover_fractions)")
    ]
    assert ["triangle"] in indexed, indexed
