"""`tin_engine.io.geopackage_write` and `crs.srs_entry` (increment 36, red).

`docs/increments/36-geopackage-output.md`, section 8: the unit half of T1 to
T7 on F1, a hand-made mesh in a `:memory:` GeoPackage, and T2 on F2's arrays.
The reader `io/geopackage.py` reads every row back, and T1 unpacks each blob's
header itself (bytes 0 to 40), because the reader skips the envelope.

F1: four vertices (0,0,0), (1,0,1), (1,1,2), (0,1,3); triangles (0,1,2) and
(0,2,3); one constraint edge, the diagonal 0-2, with the `river` bit; codes
312 and 0; `srs_id` 25833. Plan areas 0.5 and 0.5, z_mean 1 and 5/3, slopes
`atan(sqrt 2)` = 54.7356 deg and `atan(sqrt 10)` = 72.4516 deg (the triangle
plane's normal (-1,-1,1) and (1,-3,1)).

Pinned beyond the design's text (listed in the red handback):

- `write_triangles` leaves the layer registered (`gpkg_contents`,
  `gpkg_geometry_columns`) and indexed; `write_lines` leaves it registered.
  `register_layer` and `build_rtree` are theirs to call.
- `codes=None` with `palette=None` means no coded land cover: `code`,
  `land_cover` and `colour` are NULL on every row. `palette` is the code
  table itself (`palettes.CORINE_NATURAL`'s type), not its name.
- `geometry_type_name` is the OGC name in capitals: `POLYGON`, `LINESTRING`.
- `write_run(conn, rows)` takes `(name, wording, value, number)` tuples.
- `srs_entry` never raises: text pyproj refuses gives the `-1` entry.

The modules under test are imported inside fixtures, so collection succeeds
and each test fails on its own with `ModuleNotFoundError` (the writer) or
`AttributeError` (`crs.srs_entry`).
"""

from __future__ import annotations

import ast
import importlib
import math
import sqlite3
import struct
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from pyproj import CRS
from pyproj.enums import WktVersion

from cli_driver import invoke
from geotiff_fixtures import KARTVERKET, needs_codecs
from gpkg_fixtures import needs_rtree
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.io.geopackage import GpkgLayer, layer_info, query_features
from tin_engine.palettes import CORINE_NATURAL
from vtkread import polygons_as_array, read_vtk

UTM33 = 25833
FIXED_TIME = "1970-01-01T00:00:00.000Z"

VERTICES = np.array([[0, 0, 0], [1, 0, 1], [1, 1, 2], [0, 1, 3]], dtype=np.float64)
TRIANGLES = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
EDGES = np.array([[0, 2]], dtype=np.int64)
RIVER = DEFAULT_VOCABULARY.mask("river")
ROAD = DEFAULT_VOCABULARY.mask("road")
EDGE_MASKS = np.array([RIVER], dtype=np.uint32)
CODES = np.array([312, 0], dtype=np.int32)
SLOPES_DEG = (math.degrees(math.atan(math.sqrt(2))), math.degrees(math.atan(math.sqrt(10))))


@pytest.fixture(scope="module")
def gw() -> ModuleType:
    return importlib.import_module("tin_engine.io.geopackage_write")


@pytest.fixture(scope="module")
def srs_entry() -> Any:
    return importlib.import_module("tin_engine.crs").srs_entry


@pytest.fixture(scope="module")
def utm33(gw: ModuleType) -> Any:
    definition = CRS.from_epsg(UTM33).to_wkt(version=WktVersion.WKT1_GDAL)
    return gw.SrsEntry(
        srs_id=UTM33,
        organization="EPSG",
        code=UTM33,
        name="ETRS89 / UTM zone 33N",
        definition=definition,
    )


@pytest.fixture
def conn(gw: ModuleType, utm33: Any) -> Iterator[sqlite3.Connection]:
    """An in-memory GeoPackage with the `gpkg_*` tables and F1's SRS."""
    con = sqlite3.connect(":memory:")
    gw.create_geopackage(con, utm33)
    yield con
    con.close()


@pytest.fixture
def f1(gw: ModuleType, conn: sqlite3.Connection) -> sqlite3.Connection:
    """F1's `triangles` layer, the CORINE palette naming its codes."""
    gw.write_triangles(conn, VERTICES, TRIANGLES, CODES, CORINE_NATURAL, UTM33)
    return conn


def bounds(xyz: np.ndarray) -> tuple[float, float, float, float]:
    """The header envelope's order: `minx, maxx, miny, maxy` (12-128 section 2.1.3)."""
    return (xyz[:, 0].min(), xyz[:, 0].max(), xyz[:, 1].min(), xyz[:, 1].max())


def blobs(conn: sqlite3.Connection, layer: GpkgLayer) -> list[bytes]:
    column, table = layer.column, layer.table
    return [row[0] for row in conn.execute(f"SELECT {column} FROM {table} ORDER BY {layer.pk}")]


def columns(conn: sqlite3.Connection, table: str) -> dict[str, tuple[str, int]]:
    """Column name -> (declared type, primary-key position) from `PRAGMA table_info`."""
    return {row[1]: (row[2], row[5]) for row in conn.execute(f"PRAGMA table_info({table})")}


def contents(conn: sqlite3.Connection, table: str) -> tuple[Any, ...]:
    row = conn.execute(
        "SELECT data_type, last_change, min_x, min_y, max_x, max_y, srs_id "
        "FROM gpkg_contents WHERE table_name = ?",
        (table,),
    ).fetchone()
    assert row is not None, f"{table} is not in gpkg_contents"
    return tuple(row)


def geometry_column(conn: sqlite3.Connection, table: str) -> tuple[str, int, int, int]:
    row = conn.execute(
        "SELECT geometry_type_name, srs_id, z, m FROM gpkg_geometry_columns WHERE table_name = ?",
        (table,),
    ).fetchone()
    assert row is not None, f"{table} is not in gpkg_geometry_columns"
    return tuple(row)


# ----------------------------------------------------------------- T1: round trip


class TestRoundTrip:
    """T1: the reader reads each triangle back; the header bytes say what the
    design says, in the envelope's own order."""

    def test_each_triangle_reads_back_as_its_polygonz_in_triangle_order(
        self, f1: sqlite3.Connection
    ) -> None:
        layer = layer_info(f1, "triangles")
        assert layer.srs_id == UTM33
        assert layer.crs == f"EPSG:{UTM33}"
        found = list(query_features(f1, layer, None, "code"))
        assert [f.fid for f in found] == [1, 2]
        for feature, corners in zip(found, VERTICES[TRIANGLES], strict=True):
            geometry = feature.geometry
            assert geometry is not None
            assert geometry.geom_type == "Polygon"
            # A type-3 (Polygon, no z) mutant: the reader raises on the WKB
            # or reads a flat polygon; either ends here.
            assert geometry.has_z
            ring = np.vstack([corners, corners[:1]])
            assert_array_equal(np.asarray(geometry.exterior.coords), ring)
        assert [f.value for f in found] == CODES.tolist()

    def test_the_header_is_gp_version_0_little_endian_with_an_xy_envelope(
        self, f1: sqlite3.Connection
    ) -> None:
        layer = layer_info(f1, "triangles")
        for blob, corners in zip(blobs(f1, layer), VERTICES[TRIANGLES], strict=True):
            assert blob[:2] == b"GP"
            assert blob[2] == 0
            assert blob[3] == 0b0011, "flags: little-endian header, 32-byte XY envelope"
            (srs_id,) = struct.unpack_from("<i", blob, 4)
            assert srs_id == UTM33
            # The reader never sees these 32 bytes: a swapped order dies here.
            assert struct.unpack_from("<4d", blob, 8) == bounds(corners)

    def test_the_wkb_is_iso_polygonz_with_one_ring_of_four_points(
        self, f1: sqlite3.Connection
    ) -> None:
        layer = layer_info(f1, "triangles")
        for blob in blobs(f1, layer):
            assert blob[40] == 1, "WKB byte order: little-endian"
            geometry_type, rings, points = struct.unpack_from("<III", blob, 41)
            assert (geometry_type, rings, points) == (1003, 1, 4)
            assert len(blob) == 40 + 1 + 4 + 4 + 4 + 4 * 3 * 8

    def test_the_layer_is_registered_as_polygon_with_z(self, f1: sqlite3.Connection) -> None:
        minx, maxx, miny, maxy = bounds(VERTICES)
        assert contents(f1, "triangles") == (
            "features",
            FIXED_TIME,
            minx,
            miny,
            maxx,
            maxy,
            UTM33,
        )
        assert geometry_column(f1, "triangles") == ("POLYGON", UTM33, 1, 0)

    @needs_rtree
    def test_the_rtree_holds_each_triangle_box_under_its_fid(self, f1: sqlite3.Connection) -> None:
        layer = layer_info(f1, "triangles")
        assert layer.rtree is not None, "the triangles layer is indexed"
        rows = f1.execute(f"SELECT id, minx, maxx, miny, maxy FROM {layer.rtree} ORDER BY id")
        expected = [(fid, *bounds(corners)) for fid, corners in enumerate(VERTICES[TRIANGLES], 1)]
        assert rows.fetchall() == expected  # exact: 0 and 1 are float32 values


# ------------------------------------------------------- T2: index matches rows


@pytest.fixture(scope="module")
def f2_arrays(tmp_path_factory: pytest.TempPathFactory) -> tuple[np.ndarray, np.ndarray]:
    """F2: the Kartverket tile at `--tolerance 10`, as (vertices, triangles)."""
    pytest.importorskip("imagecodecs", reason="the tile is LZW; needs the `codecs` extra")
    out = tmp_path_factory.mktemp("f2") / "tile.vtk"
    code, output = invoke(
        "mesh", "--dem", str(KARTVERKET), "--tolerance", "10", "--out", str(out), "--binary"
    )
    assert code == 0, output
    vtk = read_vtk(out.read_bytes())
    return vtk.points, polygons_as_array(vtk)


def within_one_float32_ulp(stored: np.ndarray, exact: np.ndarray) -> np.ndarray:
    """SQLite's R-tree keeps float32: each stored bound is within one float32
    step of the exact float64 bound at the bound's own magnitude (UTM metres,
    about 5e5: a step of 0.03 to 0.06 m). Checked at F2, 69,571 triangles."""
    step = np.spacing(np.abs(exact).astype(np.float32)).astype(np.float64)
    return np.abs(stored - exact) <= step


@needs_codecs
@needs_rtree
class TestIndexMatchesRows:
    """T2 with `batch=30_000`: 69,571 triangles cross two batch boundaries
    and leave a partial batch of 9,571."""

    def test_every_fid_has_its_own_box_and_the_counts_agree(
        self, gw: ModuleType, conn: sqlite3.Connection, f2_arrays: tuple[np.ndarray, np.ndarray]
    ) -> None:
        vertices, triangles = f2_arrays
        assert len(triangles) == 69_571, "F2's figure in the design; a changed mesh changes T2"
        gw.write_triangles(conn, vertices, triangles, None, None, UTM33, batch=30_000)
        layer = layer_info(conn, "triangles")
        assert layer.rtree is not None
        (rows,) = conn.execute("SELECT count(*) FROM triangles").fetchone()
        (indexed,) = conn.execute(f"SELECT count(*) FROM {layer.rtree}").fetchone()
        assert rows == indexed == len(triangles)
        stored = np.array(
            conn.execute(
                f"SELECT id, minx, maxx, miny, maxy FROM {layer.rtree} ORDER BY id"
            ).fetchall(),
            dtype=np.float64,
        )
        assert_array_equal(stored[:, 0], np.arange(1, len(triangles) + 1))
        xy = vertices[triangles][:, :, :2]
        exact = np.column_stack(
            [xy[:, :, 0].min(1), xy[:, :, 0].max(1), xy[:, :, 1].min(1), xy[:, :, 1].max(1)]
        )
        off = ~within_one_float32_ulp(stored[:, 1:], exact)
        wrong = np.flatnonzero(off.any(1)) + 1
        assert wrong.size == 0, f"{wrong.size} boxes differ from their triangle; fids {wrong[:5]}"


# --------------------------------------------------------------- T3: attributes


class TestAttributes:
    """T3 on F1: the eight columns against the figures computed by hand."""

    def test_code_name_colour_heights_area_and_slope(self, f1: sqlite3.Connection) -> None:
        rows = f1.execute(
            "SELECT code, land_cover, colour, z_min, z_max, z_mean, area_m2, slope_deg "
            "FROM triangles ORDER BY fid"
        ).fetchall()
        assert [r[:3] for r in rows] == [
            (312, "Coniferous forest", "#1f5a2e"),
            (0, "No polygon; constraint lines", "#3c3c3c"),
        ]
        assert [r[3:5] for r in rows] == [(0.0, 2.0), (0.0, 3.0)]
        assert_allclose([r[5] for r in rows], [1.0, 5 / 3], rtol=1e-15)
        assert [r[6] for r in rows] == [0.5, 0.5], "plan area, with the half"
        # 1e-9 deg: both sides are float64 atan of a unit-square normal.
        assert_allclose([r[7] for r in rows], SLOPES_DEG, atol=1e-9)

    def test_a_code_the_palette_lacks_is_unnamed_and_magenta(
        self, gw: ModuleType, conn: sqlite3.Connection
    ) -> None:
        codes = np.array([312, 999], dtype=np.int32)
        gw.write_triangles(conn, VERTICES, TRIANGLES, codes, CORINE_NATURAL, UTM33)
        rows = conn.execute("SELECT code, land_cover, colour FROM triangles ORDER BY fid")
        assert rows.fetchall() == [(312, "Coniferous forest", "#1f5a2e"), (999, None, "#ff00ff")]

    def test_no_coded_land_cover_gives_null_code_on_every_row(
        self, gw: ModuleType, conn: sqlite3.Connection
    ) -> None:
        gw.write_triangles(conn, VERTICES, TRIANGLES, None, None, UTM33)
        rows = conn.execute("SELECT code, land_cover, colour FROM triangles").fetchall()
        assert rows == [(None, None, None)] * 2
        tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master")}
        assert not tables & {"land_cover", "land_cover_classes", "layer_styles"}


# ---------------------------------------------------------------- T5: run table


RUN_ROWS = [
    ("rasputin_version", "rasputin version", "0.3.1", None),
    ("command", "command", "rasputin mesh --dem tile.tif --tolerance 10 --out tile.gpkg", None),
    ("crs", "CRS", "EPSG:25833", None),
    ("tolerance_m", "Tolerance (m)", "10", 10.0),
    ("max_error_m", "Largest error (m)", "9.73", 9.73),
    ("dem_source", "DEM", "7908_3_10m_z33.tif", None),
    ("vertices", "Vertices", "38921", 38921),
]


class TestRunTable:
    """T5's unit half: the record as rows, keyed by `fid`, `name` unique."""

    def test_every_entry_is_a_row_in_record_order_with_its_number(
        self, gw: ModuleType, conn: sqlite3.Connection
    ) -> None:
        gw.write_run(conn, RUN_ROWS)
        rows = conn.execute("SELECT fid, name, wording, value, number FROM run ORDER BY fid")
        assert rows.fetchall() == [(i, *row) for i, row in enumerate(RUN_ROWS, 1)]
        assert contents(conn, "run")[:2] == ("attributes", FIXED_TIME)

    def test_fid_is_the_integer_primary_key_and_name_is_unique(
        self, gw: ModuleType, conn: sqlite3.Connection
    ) -> None:
        gw.write_run(conn, RUN_ROWS)
        assert columns(conn, "run")["fid"] == ("INTEGER", 1), "Requirement 119: a rowid alias"
        with pytest.raises(sqlite3.IntegrityError, match="UNIQUE"):
            conn.execute(
                "INSERT INTO run (name, wording, value, number) VALUES (?, ?, ?, ?)",
                ("tolerance_m", "Tolerance (m)", "5", 5.0),
            )


# -------------------------------------------------------------------- T6: lines


class TestLines:
    """T6 on F1: one row per constraint edge; one 0/1 column per property set
    on at least one edge, named as `write_vtk` names its cell arrays."""

    def test_the_diagonal_reads_back_with_its_length_and_its_river_bit(
        self, gw: ModuleType, f1: sqlite3.Connection
    ) -> None:
        gw.write_lines(f1, VERTICES, EDGES, EDGE_MASKS, DEFAULT_VOCABULARY, UTM33)
        layer = layer_info(f1, "lines")
        assert layer.srs_id == UTM33
        assert layer.rtree is None, "no R-tree on lines"
        (found,) = query_features(f1, layer, None, "length_m")
        assert found.fid == 1
        assert found.geometry is not None and found.geometry.geom_type == "LineString"
        assert found.geometry.has_z
        assert_array_equal(np.asarray(found.geometry.coords), VERTICES[[0, 2]])
        assert found.value == pytest.approx(math.sqrt(2), rel=1e-15)
        attributes = set(columns(f1, "lines")) - {layer.pk, layer.column}
        assert attributes == {"length_m", "river"}
        assert f1.execute("SELECT river FROM lines").fetchall() == [(1,)]
        assert geometry_column(f1, "lines") == ("LINESTRING", UTM33, 1, 0)
        assert contents(f1, "lines") == ("features", FIXED_TIME, 0.0, 0.0, 1.0, 1.0, UTM33)

    def test_a_property_set_on_no_edge_gets_no_column(
        self, gw: ModuleType, f1: sqlite3.Connection
    ) -> None:
        edges = np.array([[0, 2], [1, 3]], dtype=np.int64)
        masks = np.array([RIVER, RIVER | ROAD], dtype=np.uint32)
        gw.write_lines(f1, VERTICES, edges, masks, DEFAULT_VOCABULARY, UTM33)
        layer = layer_info(f1, "lines")
        attributes = set(columns(f1, "lines")) - {layer.pk, layer.column}
        assert attributes == {"length_m", "river", "road"}
        rows = f1.execute("SELECT river, road FROM lines ORDER BY fid").fetchall()
        assert rows == [(1, 0), (1, 1)]

    def test_the_line_header_envelope_is_minx_maxx_miny_maxy(
        self, gw: ModuleType, f1: sqlite3.Connection
    ) -> None:
        gw.write_lines(f1, VERTICES, EDGES, EDGE_MASKS, DEFAULT_VOCABULARY, UTM33)
        (blob,) = blobs(f1, layer_info(f1, "lines"))
        assert blob[3] == 0b0011
        assert struct.unpack_from("<4d", blob, 8) == bounds(VERTICES[[0, 2]])
        assert struct.unpack_from("<II", blob, 41) == (1002, 2), "LineStringZ of two points"


# ---------------------------------------------------------------------- T7: SRS


TMERC_NO_CODE = (
    "+proj=tmerc +lat_0=0 +lon_0=15.1 +k=0.9996 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
)


class TestSrsEntry:
    """T7's unit half: `crs.srs_entry` and the SRS rows `create_geopackage` writes."""

    def test_an_epsg_code_is_its_own_srs_id_with_wkt1(self, srs_entry: Any) -> None:
        entry = srs_entry("EPSG:25833")
        assert (entry.srs_id, entry.organization, entry.code) == (UTM33, "EPSG", UTM33)
        assert entry.definition.startswith("PROJCS["), "WKT1_GDAL, what QGIS and GDAL read"
        assert "UTM zone 33N" in entry.name

    def test_a_crs_without_a_code_gets_100000_and_wkt2(self, srs_entry: Any) -> None:
        entry = srs_entry(TMERC_NO_CODE)
        assert (entry.srs_id, entry.organization) == (100000, "rasputin")
        assert entry.definition.startswith("PROJCRS[")
        assert CRS(entry.definition).equals(CRS(TMERC_NO_CODE))

    @pytest.mark.parametrize("text", [None, "free text", "not-a-crs"])
    def test_no_crs_or_text_pyproj_refuses_is_the_undefined_row(
        self, srs_entry: Any, text: str | None
    ) -> None:
        assert srs_entry(text).srs_id == -1

    def test_the_file_holds_the_reserved_rows_wgs84_and_the_run_crs(
        self, f1: sqlite3.Connection
    ) -> None:
        rows = f1.execute(
            "SELECT srs_id, organization, organization_coordsys_id FROM gpkg_spatial_ref_sys"
        ).fetchall()
        assert {r[0] for r in rows} == {-1, 0, 4326, UTM33}
        assert (UTM33, "EPSG", UTM33) in rows
        assert (4326, "EPSG", 4326) in rows

    def test_the_undefined_entry_adds_no_second_row(self, gw: ModuleType, srs_entry: Any) -> None:
        con = sqlite3.connect(":memory:")
        gw.create_geopackage(con, srs_entry(None))
        ids = [r[0] for r in con.execute("SELECT srs_id FROM gpkg_spatial_ref_sys ORDER BY 1")]
        assert ids == [-1, 0, 4326]
        gw.write_triangles(con, VERTICES, TRIANGLES, None, None, -1)
        assert layer_info(con, "triangles").srs_id == -1
        con.close()


def test_the_writer_opens_no_file_and_imports_no_core_pyproj_or_path() -> None:
    """Section 5: pure in the reader's sense. Judged on the source's imports
    and calls (`ast`), so a docstring that names the rule does not trip it."""
    src = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "io"
    tree = ast.parse((src / "geopackage_write.py").read_text(encoding="utf-8"))
    imported: set[str] = set()
    calls: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {alias.name for alias in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
            imported |= {f"{node.module}.{alias.name}" for alias in node.names}
        elif isinstance(node, ast.Call):
            calls.add(ast.unparse(node.func))
    forbidden = {"pyproj", "pathlib", "tin_engine._core", "tin_engine.crs"}
    assert not {m for m in imported if m.split(".")[0] in forbidden or m in forbidden}
    assert not calls & {"open", "sqlite3.connect", "Path"}, sorted(calls)
