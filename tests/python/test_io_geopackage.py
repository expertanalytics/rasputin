"""`tin_engine.io.geopackage` and `io.repository.open_geopackage` (increment 16b-1).

`docs/increments/16b-terrain-polygons.md` R2, R3 and "Tests for @tester"
(GeoPackage decoding, layer resolution). Every GeoPackage here is a real file
written to `tmp_path` by `gpkg_fixtures.write_gpkg`; every blob is built from
OGC 12-128's header layout by `gpkg_fixtures.blob`, never by the reader.

Pinned by this suite (see "Pinned by the red suite (16b-1/2)"):

- `io.repository.open_geopackage(path) -> sqlite3.Connection`, read-only: a
  write through it fails, and a missing file is not created.
- `io.geopackage.GeoPackageError`, a `ValueError`: every refusal below.
- `layer_info(conn, table) -> GpkgLayer`, `table=None` meaning "the only
  features table". `GpkgLayer` is frozen, with `table`, `column` (the
  geometry column), `pk` (the primary-key column), `srs_id`, `crs` (text
  `crs.parse_crs` reads) and `rtree` (the R-tree table's name, or `None`).
- `decode_geometry(blob) -> (srs_id, geometry)`; `geometry` is `None` (or
  empty) for a blob with the empty flag.
- `query_features(conn, layer, box, attribute)` yields objects with `fid`
  (the primary key), `geometry` and `value` (the `attribute` column), in
  ascending primary-key order whatever order the rows or the R-tree hold.
  `box` is `(minx, miny, maxx, maxy)` in the layer's CRS, closed.
- Refusals of a row name the layer's table and the row's primary key.

Committed red at `e99c8ea`: the modules under test did not exist yet, and
because they are imported inside fixtures each test failed on its own with
`ModuleNotFoundError` or `AttributeError` while collection was unaffected. They
landed in `5079da8` and the suite has been green since.
"""

from __future__ import annotations

import ast
import importlib
import sqlite3
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import shapely
from pyproj import CRS
from shapely.geometry import LineString, MultiPolygon, Polygon, box

from gpkg_fixtures import (
    HAS_RTREE,
    Layer,
    Row,
    Srs,
    blob,
    needs_rtree,
    write_gpkg,
)

SRC = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine"
SQUARE = box(0.0, 0.0, 10.0, 10.0)
LAEA = 3035


@pytest.fixture(scope="module")
def gpkg() -> ModuleType:
    return importlib.import_module("tin_engine.io.geopackage")


@pytest.fixture(scope="module")
def open_geopackage() -> Any:
    return importlib.import_module("tin_engine.io.repository").open_geopackage


def squares(*pks: int, spacing: float = 100.0) -> list[Row]:
    """One 10 m square per key, `spacing` apart in x, `Code_18` = `1<pk>`."""
    return [
        Row(pk, box(k * spacing, 0.0, k * spacing + 10.0, 10.0), {"Code_18": f"1{pk:02d}"})
        for k, pk in enumerate(pks)
    ]


def one_layer(tmp_path: Path, rows: list[Row], **kwargs: Any) -> Path:
    return write_gpkg(tmp_path / "f.gpkg", [Layer("clc", LAEA, rows, **kwargs)])


def read_all(
    gpkg: ModuleType,
    open_geopackage: Any,
    path: Path,
    table: str | None = None,
    bounds: tuple[float, float, float, float] | None = None,
) -> list[Any]:
    conn = open_geopackage(path)
    try:
        layer = gpkg.layer_info(conn, table)
        return list(gpkg.query_features(conn, layer, bounds, "Code_18"))
    finally:
        conn.close()


# ------------------------------------------------------------ open_geopackage


class TestOpenGeopackage:
    def test_it_returns_a_connection_that_reads(self, tmp_path: Path, open_geopackage: Any) -> None:
        path = one_layer(tmp_path, squares(1))
        conn = open_geopackage(path)
        try:
            assert isinstance(conn, sqlite3.Connection)
            assert conn.execute("PRAGMA application_id").fetchone()[0] == 0x47504B47
        finally:
            conn.close()

    def test_it_is_read_only(self, tmp_path: Path, open_geopackage: Any) -> None:
        """R3: read-only, which also keeps the SpatiaLite triggers from firing."""
        path = one_layer(tmp_path, squares(1))
        before = path.read_bytes()
        conn = open_geopackage(path)
        try:
            with pytest.raises(sqlite3.OperationalError, match="readonly"):
                conn.execute("DELETE FROM clc")
        finally:
            conn.close()
        assert path.read_bytes() == before

    def test_a_missing_file_is_refused_and_not_created(
        self, tmp_path: Path, open_geopackage: Any
    ) -> None:
        missing = tmp_path / "missing.gpkg"
        with pytest.raises((sqlite3.OperationalError, OSError, ValueError)):
            open_geopackage(missing).execute("SELECT 1").fetchone()
        assert not missing.exists()

    def test_a_path_with_uri_characters_opens_the_named_file(
        self, tmp_path: Path, open_geopackage: Any
    ) -> None:
        """`path.resolve().as_uri()`: a `?` or `#` or space in the name is the
        file's, not a URI query."""
        directory = tmp_path / "a dir#1"
        directory.mkdir()
        path = write_gpkg(directory / "x?y.gpkg", [Layer("clc", LAEA, squares(1))])
        conn = open_geopackage(path)
        try:
            assert conn.execute("SELECT COUNT(*) FROM clc").fetchone()[0] == 1
        finally:
            conn.close()


# ----------------------------------------------------------- decode_geometry


class TestDecodeGeometry:
    """R3's header, every branch."""

    @pytest.mark.parametrize("envelope", [0, 1, 2, 3, 4])
    @pytest.mark.parametrize("big_endian", [False, True])
    def test_every_envelope_code_in_either_byte_order(
        self, gpkg: ModuleType, envelope: int, big_endian: bool
    ) -> None:
        shape = Polygon([(0.5, 0.25), (7.0, 0.5), (3.0, 9.0)], [[(2, 2), (4, 2), (3, 4)]])
        srs_id, geometry = gpkg.decode_geometry(
            blob(shape, 25833, envelope=envelope, big_endian=big_endian)
        )
        assert srs_id == 25833
        assert shapely.equals_exact(geometry, shape, tolerance=0.0)

    def test_big_endian_wkb_after_a_little_endian_header(self, gpkg: ModuleType) -> None:
        shape = MultiPolygon([SQUARE, box(20, 20, 30, 30)])
        _, geometry = gpkg.decode_geometry(blob(shape, LAEA, wkb_big_endian=True))
        assert shapely.equals_exact(geometry, shape, tolerance=0.0)

    def test_coordinates_come_out_bit_for_bit(self, gpkg: ModuleType) -> None:
        shape = LineString([(4_321_987.123456789, 5_432_109.987654321), (0.1, 0.2)])
        _, geometry = gpkg.decode_geometry(blob(shape, LAEA, envelope=1))
        assert list(geometry.coords) == list(shape.coords)

    def test_the_empty_flag_gives_no_geometry(self, gpkg: ModuleType) -> None:
        srs_id, geometry = gpkg.decode_geometry(blob(None, LAEA, empty=True))
        assert srs_id == LAEA
        assert geometry is None or geometry.is_empty

    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({"extended": True}, id="extended-type"),
            pytest.param({"magic": b"GQ"}, id="bad-magic"),
            pytest.param({"version": 1}, id="version-1"),
            pytest.param({"envelope": 5}, id="envelope-5"),
            pytest.param({"envelope": 6}, id="envelope-6"),
            pytest.param({"envelope": 7}, id="envelope-7"),
        ],
    )
    def test_a_header_it_cannot_read_is_refused(
        self, gpkg: ModuleType, kwargs: dict[str, Any]
    ) -> None:
        with pytest.raises(gpkg.GeoPackageError):
            gpkg.decode_geometry(blob(SQUARE, LAEA, **kwargs))

    @pytest.mark.parametrize("length", [0, 1, 4, 7])
    def test_a_blob_shorter_than_the_header_is_refused(self, gpkg: ModuleType, length: int) -> None:
        with pytest.raises(gpkg.GeoPackageError):
            gpkg.decode_geometry(blob(SQUARE, LAEA)[:length])

    def test_an_envelope_longer_than_the_blob_is_refused(self, gpkg: ModuleType) -> None:
        with pytest.raises(gpkg.GeoPackageError):
            gpkg.decode_geometry(blob(SQUARE, LAEA, envelope=4)[:40])

    def test_corrupt_wkb_is_refused_as_a_geopackage_error(self, gpkg: ModuleType) -> None:
        with pytest.raises(gpkg.GeoPackageError):
            gpkg.decode_geometry(blob(SQUARE, LAEA)[:-5])

    def test_the_error_is_a_value_error(self, gpkg: ModuleType) -> None:
        assert issubclass(gpkg.GeoPackageError, ValueError)


# ---------------------------------------------------------------- layer_info


class TestLayerInfo:
    def test_the_only_features_table_is_the_default(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        path = one_layer(tmp_path, squares(1, 2))
        conn = open_geopackage(path)
        try:
            layer = gpkg.layer_info(conn, None)
        finally:
            conn.close()
        assert (layer.table, layer.column, layer.pk, layer.srs_id) == (
            "clc",
            "Shape",
            "OBJECTID",
            LAEA,
        )
        assert CRS.from_user_input(layer.crs) == CRS.from_epsg(LAEA)
        assert layer.rtree == ("rtree_clc_Shape" if HAS_RTREE else None)

    def test_the_layer_is_frozen(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        conn = open_geopackage(one_layer(tmp_path, squares(1)))
        try:
            layer = gpkg.layer_info(conn, "clc")
        finally:
            conn.close()
        with pytest.raises((TypeError, ValueError, AttributeError)):
            layer.table = "other"

    def test_several_tables_need_a_name_and_the_refusal_lists_them(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        path = write_gpkg(
            tmp_path / "f.gpkg",
            [Layer("europe", LAEA, squares(1)), Layer("reunion", 2975, squares(2))],
        )
        conn = open_geopackage(path)
        try:
            with pytest.raises(gpkg.GeoPackageError) as caught:
                gpkg.layer_info(conn, None)
            assert "europe" in str(caught.value) and "reunion" in str(caught.value)
            assert gpkg.layer_info(conn, "reunion").srs_id == 2975
        finally:
            conn.close()

    def test_an_unknown_name_is_refused_naming_it(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        conn = open_geopackage(one_layer(tmp_path, squares(1)))
        try:
            with pytest.raises(gpkg.GeoPackageError, match="nosuch"):
                gpkg.layer_info(conn, "nosuch")
        finally:
            conn.close()

    def test_a_table_that_is_not_features_is_refused(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        path = write_gpkg(
            tmp_path / "f.gpkg",
            [Layer("clc", LAEA, squares(1)), Layer("pixels", LAEA, [], data_type="tiles")],
        )
        conn = open_geopackage(path)
        try:
            with pytest.raises(gpkg.GeoPackageError, match="pixels"):
                gpkg.layer_info(conn, "pixels")
            assert gpkg.layer_info(conn, None).table == "clc"  # tiles are not candidates
        finally:
            conn.close()

    def test_a_non_epsg_organization_takes_the_definition_wkt(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        wkt = CRS.from_epsg(25833).to_wkt()
        path = write_gpkg(
            tmp_path / "f.gpkg",
            [Layer("clc", 990001, squares(1))],
            srs=[Srs(990001, organization="rasputin", code=1, definition=wkt)],
        )
        conn = open_geopackage(path)
        try:
            layer = gpkg.layer_info(conn, None)
        finally:
            conn.close()
        assert CRS.from_user_input(layer.crs) == CRS.from_epsg(25833)

    def test_the_epsg_organization_is_case_insensitive(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        path = write_gpkg(
            tmp_path / "f.gpkg",
            [Layer("clc", 25833, squares(1))],
            srs=[Srs(25833, organization="epsg", definition="undefined")],
        )
        conn = open_geopackage(path)
        try:
            layer = gpkg.layer_info(conn, None)
        finally:
            conn.close()
        assert CRS.from_user_input(layer.crs) == CRS.from_epsg(25833)

    def test_a_file_that_is_not_a_geopackage_is_refused(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        path = write_gpkg(tmp_path / "f.gpkg", [Layer("clc", LAEA, squares(1))], application_id=0)
        conn = open_geopackage(path)
        try:
            with pytest.raises(gpkg.GeoPackageError):
                gpkg.layer_info(conn, None)
        finally:
            conn.close()

    def test_no_rtree_is_reported_as_none(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        conn = open_geopackage(one_layer(tmp_path, squares(1), rtree=False))
        try:
            assert gpkg.layer_info(conn, None).rtree is None
        finally:
            conn.close()


# ------------------------------------------------------------ query_features


class TestQueryFeatures:
    @needs_rtree
    def test_rows_come_in_primary_key_order_whatever_the_rtree_holds(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        """R3: `ORDER BY <pk>`. The table's rows are inserted out of key order
        and the R-tree's in reverse, so neither rowid nor R-tree order is it."""
        rows = squares(5, 3, 9, 1, 7)
        path = one_layer(tmp_path, rows, insert_order=[9, 1, 5, 7, 3], rtree_order=[7, 3, 5, 1, 9])
        got = read_all(gpkg, open_geopackage, path, bounds=(-1.0, -1.0, 1000.0, 20.0))
        assert [f.fid for f in got] == [1, 3, 5, 7, 9]
        assert [f.value for f in got] == ["101", "103", "105", "107", "109"]

    def test_without_an_rtree_rows_also_come_in_primary_key_order(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        path = one_layer(tmp_path, squares(5, 3, 9), rtree=False, insert_order=[9, 3, 5])
        got = read_all(gpkg, open_geopackage, path, bounds=(-1.0, -1.0, 1000.0, 20.0))
        assert [f.fid for f in got] == [3, 5, 9]

    def test_the_geometry_is_the_rows(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        rows = squares(1, 2)
        got = read_all(gpkg, open_geopackage, one_layer(tmp_path, rows))
        for feature, row in zip(got, rows, strict=True):
            assert shapely.equals_exact(feature.geometry, row.geometry, tolerance=0.0)

    @pytest.mark.parametrize("rtree", [pytest.param(True, marks=needs_rtree), False])
    def test_the_box_keeps_every_row_meeting_it_closed(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any, rtree: bool
    ) -> None:
        """Squares at x 0-10, 100-110, 200-210, 300-310. A box whose min x is
        exactly 110 meets the second square's edge: closed, so it is kept.
        Through the R-tree the far squares do not come back (its bounds are
        rounded outward, by far less than 90 m); without one the table is
        scanned (R3), and a candidate superset is all a query promises."""
        path = one_layer(tmp_path, squares(1, 2, 3, 4), rtree=rtree)
        got = read_all(gpkg, open_geopackage, path, bounds=(110.0, 2.0, 205.0, 3.0))
        fids = [f.fid for f in got]
        assert 2 in fids and 3 in fids
        if rtree:
            assert 1 not in fids and 4 not in fids

    def test_no_box_reads_every_row(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        got = read_all(gpkg, open_geopackage, one_layer(tmp_path, squares(1, 2, 3)))
        assert [f.fid for f in got] == [1, 2, 3]

    def test_an_empty_flag_row_is_yielded_empty(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        """R1: skipped and counted, which is the caller's; the reader must not
        drop it silently or refuse it."""
        rows = [*squares(1), Row(2, None, {"Code_18": "311"}, raw=blob(None, LAEA, empty=True))]
        got = read_all(gpkg, open_geopackage, one_layer(tmp_path, rows, rtree=False))
        assert [f.fid for f in got] == [1, 2]
        assert got[1].geometry is None or got[1].geometry.is_empty

    def test_a_row_in_another_srs_is_refused_naming_the_table_and_key(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        rows = [*squares(1), Row(4242, SQUARE, {"Code_18": "311"}, raw=blob(SQUARE, 4326))]
        path = one_layer(tmp_path, rows)
        with pytest.raises(gpkg.GeoPackageError) as caught:
            read_all(gpkg, open_geopackage, path)
        assert "4242" in str(caught.value) and "clc" in str(caught.value)

    def test_a_row_with_a_bad_header_is_refused_naming_its_key(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        rows = [
            *squares(1),
            Row(777, SQUARE, {"Code_18": "311"}, raw=blob(SQUARE, LAEA, envelope=5)),
        ]
        with pytest.raises(gpkg.GeoPackageError, match="777"):
            read_all(gpkg, open_geopackage, one_layer(tmp_path, rows, rtree=False))

    def test_a_missing_attribute_column_is_refused_naming_it(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        conn = open_geopackage(one_layer(tmp_path, squares(1)))
        try:
            layer = gpkg.layer_info(conn, None)
            with pytest.raises(gpkg.GeoPackageError, match="Nope_18"):
                list(gpkg.query_features(conn, layer, None, "Nope_18"))
        finally:
            conn.close()

    @needs_rtree
    def test_sqlite_without_the_rtree_module_is_refused_naming_the_cause(
        self, tmp_path: Path, gpkg: ModuleType, open_geopackage: Any
    ) -> None:
        """R3: `no such module: rtree` becomes a refusal naming the cause, never
        a silent full scan. Simulated: a connection whose every statement that
        touches the R-tree fails as SQLite without the module fails."""
        conn = open_geopackage(one_layer(tmp_path, squares(1, 2)))
        try:
            layer = gpkg.layer_info(conn, None)
            with pytest.raises(gpkg.GeoPackageError, match=r"(?i)r-?tree"):
                list(gpkg.query_features(_NoRtree(conn), layer, (0.0, 0.0, 5.0, 5.0), "Code_18"))
        finally:
            conn.close()


class _NoRtree:
    """A connection whose R-tree does not exist as a module."""

    def __init__(self, conn: sqlite3.Connection) -> None:
        self._conn = conn

    def _check(self, sql: str) -> None:
        if "rtree_" in sql:
            raise sqlite3.OperationalError("no such module: rtree")

    def execute(self, sql: str, *args: Any) -> sqlite3.Cursor:
        self._check(sql)
        return self._conn.execute(sql, *args)

    def cursor(self) -> _NoRtree:
        return self

    def __getattr__(self, name: str) -> Any:
        return getattr(self._conn, name)


# ------------------------------------------------------------ the boundaries


def _calls(path: Path) -> set[str]:
    tree = ast.parse(path.read_text())
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name):
                names.add(func.id)
            elif isinstance(func, ast.Attribute):
                names.add(func.attr)
    return names


def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text())
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            found.add(("." * node.level) + (node.module or ""))
    return found


class TestBoundaries:
    """R2: `io/repository.py` stays the one module in `io/` that opens files;
    the decoder's stream is the open connection."""

    def test_the_decoder_opens_nothing(self, gpkg: ModuleType) -> None:
        path = SRC / "io" / "geopackage.py"
        calls = _calls(path)
        assert "connect" not in calls and "open" not in calls
        assert not any(m.split(".")[-1] == "pathlib" for m in _imports(path))

    def test_the_decoder_imports_no_core(self, gpkg: ModuleType) -> None:
        assert not any("_core" in m for m in _imports(SRC / "io" / "geopackage.py"))

    def test_only_the_repository_connects(self, gpkg: ModuleType) -> None:
        connecting = sorted(
            str(p.relative_to(SRC))
            for p in SRC.rglob("*.py")
            if "connect" in _calls(p) and "sqlite3" in p.read_text()
        )
        assert connecting == ["io/repository.py"]
