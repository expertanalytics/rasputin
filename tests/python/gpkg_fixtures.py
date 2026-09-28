"""Test support for increment 16b: minimal GeoPackages and geometry blobs.

`docs/increments/16b-terrain-polygons.md`, "Tests for @tester": "a test helper
writes a minimal GeoPackage (the four `gpkg_*` tables, one feature table, its
R-tree) with `sqlite3` to `tmp_path`, so the reader opens a real file". Nothing
here imports `tin_engine`; the blobs are built from OGC 12-128r18 section
2.1.3 as the design states it (R3), never from the reader under test.

Also the committed CORINE extract and Ola's local data, by path.
"""

from __future__ import annotations

import sqlite3
import struct
from collections.abc import Iterable, Sequence
from contextlib import closing
from dataclasses import dataclass, field
from pathlib import Path

import pytest
import shapely
from pyproj import CRS

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
CORINE_DIR = FIXTURES / "corine"
#: The committed extract: CORINE 2018 over the benchmark tile, EPSG:3035,
#: written by `tests/fixtures/corine/extract.py` from Ola's GeoPackage.
EXTRACT = CORINE_DIR / "clc2018_7908_3.gpkg"
EXTRACT_TABLE = "U2018_CLC2018_V2020_20u1"
#: The legacy GML, kept by Ola's Q6 (b): 200 CORINE polygons, EPSG:4326, GML2.
LEGACY_GML = CORINE_DIR / "0000_4326_corine2018_4e6064_GML.gml"

#: Ola's local data (never in CI). Tests using them skip when absent.
RASPUTIN_DATA = Path(__file__).resolve().parents[3] / "rasputin_data"
DTM10 = RASPUTIN_DATA / "DTM10_UTM33_20260925"
OLA_EUROPE = (
    RASPUTIN_DATA
    / "corine_sql"
    / "u2018_clc2018_v2020_20u1_geoPackage"
    / "DATA"
    / "U2018_CLC2018_V2020_20u1.gpkg"
)
OLA_NORWAY = RASPUTIN_DATA / "corine2018_dtm10_utm33.gpkg"

GPKG_APPLICATION_ID = 0x47504B47  # "GPKG"

# Envelope code -> number of doubles in the header (R3).
ENVELOPE_DOUBLES = {0: 0, 1: 4, 2: 6, 3: 6, 4: 8}


def _rtree_available() -> bool:
    with closing(sqlite3.connect(":memory:")) as con:
        try:
            con.execute("CREATE VIRTUAL TABLE t USING rtree(id, minx, maxx, miny, maxy)")
        except sqlite3.OperationalError:
            return False
    return True


#: R3: SQLite's R-tree module is compiled into the venv's Python; CI's Linux
#: Python is unverified. The suite detects it and says so rather than failing
#: on a missing module that is not the reader's fault.
HAS_RTREE = _rtree_available()
needs_rtree = pytest.mark.skipif(not HAS_RTREE, reason="this SQLite has no R-tree module")


def blob(
    geometry: shapely.Geometry | None,
    srs_id: int,
    *,
    envelope: int = 0,
    big_endian: bool = False,
    empty: bool = False,
    extended: bool = False,
    magic: bytes = b"GP",
    version: int = 0,
    wkb_big_endian: bool = False,
) -> bytes:
    """A GeoPackage geometry blob: `GP`, version, flags, srs_id, envelope, WKB.

    Flags: bit 0 the header's byte order (1 little-endian), bits 1-3 the
    envelope code, bit 4 empty, bit 5 extended type. An envelope code outside
    0-4 writes no envelope bytes: the reader must refuse the code itself.
    """
    flags = (0 if big_endian else 1) | (envelope << 1) | (int(empty) << 4) | (int(extended) << 5)
    order = ">" if big_endian else "<"
    head = magic + bytes([version, flags]) + struct.pack(f"{order}i", srs_id)
    doubles = ENVELOPE_DOUBLES.get(envelope, 0)
    if doubles and geometry is not None and not geometry.is_empty:
        minx, miny, maxx, maxy = geometry.bounds
        values = [minx, maxx, miny, maxy, 0.0, 0.0, 0.0, 0.0][:doubles]
    else:
        values = [0.0] * doubles
    head += struct.pack(f"{order}{doubles}d", *values)
    if geometry is None:
        geometry = shapely.from_wkt("POLYGON EMPTY")
    return head + shapely.to_wkb(geometry, byte_order=0 if wkb_big_endian else 1)


@dataclass
class Row:
    """One feature row: its primary key, geometry (or a raw blob) and columns."""

    pk: int
    geometry: shapely.Geometry | None
    values: dict[str, object] = field(default_factory=dict)
    raw: bytes | None = None


@dataclass
class Layer:
    """One feature table and how it is described."""

    table: str
    srs_id: int
    rows: Sequence[Row]
    columns: Sequence[str] = ("Code_18",)
    column: str = "Shape"
    pk: str = "OBJECTID"
    rtree: bool = True
    data_type: str = "features"
    #: The R-tree's rows are inserted in this order of primary keys (default:
    #: reverse), so a query driven from the R-tree does not meet the table's
    #: rowid order by accident.
    rtree_order: Sequence[int] | None = None
    #: Table rows are inserted in this order of primary keys (default: as given).
    insert_order: Sequence[int] | None = None


@dataclass
class Srs:
    srs_id: int
    organization: str = "EPSG"
    code: int | None = None
    definition: str | None = None


def write_gpkg(
    path: Path,
    layers: Iterable[Layer],
    *,
    srs: Iterable[Srs] = (),
    application_id: int = GPKG_APPLICATION_ID,
    user_version: int = 10200,
) -> Path:
    """Write a minimal GeoPackage to ``path`` and return it.

    Every layer's ``srs_id`` gets an ``EPSG`` row unless ``srs`` describes it.
    """
    layers = list(layers)
    described = {s.srs_id: s for s in srs}
    for layer in layers:
        described.setdefault(layer.srs_id, Srs(layer.srs_id))
    path.unlink(missing_ok=True)
    with closing(sqlite3.connect(path)) as con:
        con.executescript(
            f"""
            PRAGMA application_id = {application_id};
            PRAGMA user_version = {user_version};
            CREATE TABLE gpkg_spatial_ref_sys (srs_name TEXT NOT NULL,
              srs_id INTEGER PRIMARY KEY, organization TEXT NOT NULL,
              organization_coordsys_id INTEGER NOT NULL, definition TEXT NOT NULL,
              description TEXT);
            CREATE TABLE gpkg_contents (table_name TEXT NOT NULL PRIMARY KEY,
              data_type TEXT NOT NULL, identifier TEXT UNIQUE, description TEXT DEFAULT '',
              last_change DATETIME NOT NULL DEFAULT '2026-09-28T00:00:00.000Z',
              min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE, srs_id INTEGER);
            CREATE TABLE gpkg_geometry_columns (table_name TEXT NOT NULL,
              column_name TEXT NOT NULL, geometry_type_name TEXT NOT NULL,
              srs_id INTEGER NOT NULL, z TINYINT NOT NULL, m TINYINT NOT NULL,
              CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name));
            CREATE TABLE gpkg_extensions (table_name TEXT, column_name TEXT,
              extension_name TEXT NOT NULL, definition TEXT NOT NULL, scope TEXT NOT NULL,
              CONSTRAINT ge_tce UNIQUE (table_name, column_name, extension_name));
            """
        )
        for s in described.values():
            code = s.srs_id if s.code is None else s.code
            definition = s.definition
            if definition is None:
                definition = (
                    CRS.from_epsg(code).to_wkt() if s.organization.upper() == "EPSG" else ""
                )
            con.execute(
                "INSERT INTO gpkg_spatial_ref_sys VALUES (?,?,?,?,?,?)",
                (f"srs {s.srs_id}", s.srs_id, s.organization, code, definition, None),
            )
        for layer in layers:
            _write_layer(con, layer)
        con.commit()
    return path


def _write_layer(con: sqlite3.Connection, layer: Layer) -> None:
    columns = "".join(f", {c}" for c in layer.columns)
    con.execute(
        f"CREATE TABLE {layer.table} ({layer.pk} INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, "
        f"{layer.column} GEOMETRY{columns})"
    )
    con.execute(
        "INSERT INTO gpkg_contents (table_name, data_type, identifier, srs_id) VALUES (?,?,?,?)",
        (layer.table, layer.data_type, layer.table, layer.srs_id),
    )
    if layer.data_type == "features":
        con.execute(
            "INSERT INTO gpkg_geometry_columns VALUES (?, ?, 'GEOMETRY', ?, 0, 0)",
            (layer.table, layer.column, layer.srs_id),
        )
    by_pk = {row.pk: row for row in layer.rows}
    order = list(layer.insert_order) if layer.insert_order is not None else list(by_pk)
    marks = ",".join("?" * (len(layer.columns) + 2))
    for pk in order:
        row = by_pk[pk]
        data = row.raw if row.raw is not None else blob(row.geometry, layer.srs_id)
        values = [row.values.get(c) for c in layer.columns]
        con.execute(f"INSERT INTO {layer.table} VALUES ({marks})", (pk, data, *values))
    if not layer.rtree or not HAS_RTREE:
        return
    name = f"rtree_{layer.table}_{layer.column}"
    con.execute(f"CREATE VIRTUAL TABLE {name} USING rtree(id, minx, maxx, miny, maxy)")
    con.execute(
        "INSERT INTO gpkg_extensions VALUES (?, ?, 'gpkg_rtree_index', "
        "'http://www.geopackage.org/spec120/#extension_rtree', 'write-only')",
        (layer.table, layer.column),
    )
    rtree_order = (
        list(layer.rtree_order) if layer.rtree_order is not None else sorted(by_pk, reverse=True)
    )
    for pk in rtree_order:
        geometry = by_pk[pk].geometry
        if geometry is None or geometry.is_empty:
            continue
        minx, miny, maxx, maxy = geometry.bounds
        con.execute(f"INSERT INTO {name} VALUES (?,?,?,?,?)", (pk, minx, maxx, miny, maxy))


def copy_reversed(source: Path, target: Path, table: str = EXTRACT_TABLE) -> Path:
    """``source`` (a GeoPackage in the committed extract's layout: flags
    ``0x01``, no envelope) rewritten with its feature rows and its R-tree rows
    inserted in reverse primary-key order: I3's shuffled GeoPackage."""
    columns = ("Code_18", "Remark", "Area_Ha", "ID")
    with closing(sqlite3.connect(source.resolve().as_uri() + "?mode=ro", uri=True)) as src:
        (srs_id,) = src.execute(
            "SELECT srs_id FROM gpkg_geometry_columns WHERE table_name = ?", (table,)
        ).fetchone()
        (definition,) = src.execute(
            "SELECT definition FROM gpkg_spatial_ref_sys WHERE srs_id = ?", (srs_id,)
        ).fetchone()
        fetched = src.execute(
            f"SELECT OBJECTID, Shape, {', '.join(columns)} FROM {table} ORDER BY OBJECTID"
        ).fetchall()
    rows = [
        Row(pk, shapely.from_wkb(data[8:]), dict(zip(columns, rest, strict=True)), raw=data)
        for pk, data, *rest in fetched
    ]
    pks = [row.pk for row in rows]
    layer = Layer(
        table,
        srs_id,
        rows,
        columns=columns,
        insert_order=pks[::-1],
        rtree_order=pks[::-1],
    )
    return write_gpkg(target, [layer], srs=[Srs(srs_id, definition=definition)])
