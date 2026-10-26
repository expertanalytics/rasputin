"""A mesh as a GeoPackage, written without GDAL.

The mirror of `io.geopackage`. The stream is an open ``sqlite3.Connection``;
this module opens nothing, knows no path, no CRS object and no ``_core``. The
blob header is OGC 12-128r18 section 2.1.3 (``GP``, version 0, flags ``0b0011``:
little-endian with a 32-byte XY envelope stored as ``minx, maxx, miny,
maxy``), the R-tree Annex F.3, attribute tables section 2.3.

Nothing here reads the clock: ``last_change`` and ``update_time`` are the
fixed text `FIXED_TIME`, so the same run gives the same bytes.
"""

from __future__ import annotations

import sqlite3
import struct
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from typing import Any
from xml.sax.saxutils import escape

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry.base import BaseGeometry

from tin_engine.features import EdgeVocabulary
from tin_engine.io.geopackage import APPLICATION_ID
from tin_engine.io.models import SHARES, UNDEFINED_SRS, CellTable, ClassFractions, SrsEntry

#: What every ``last_change`` and ``update_time`` says, so the same run gives
#: the same bytes.
FIXED_TIME = "1970-01-01T00:00:00.000Z"
#: The colour of a code the palette lacks: the ParaView preset's NaN colour.
UNNAMED_COLOUR = "#ff00ff"

Palette = Mapping[int, tuple[str, str]]
#: The column type of each esri field type a features table takes.
_COLUMN_TYPES = {
    "esriFieldTypeString": "TEXT",
    "esriFieldTypeInteger": "INTEGER",
    "esriFieldTypeSmallInteger": "INTEGER",
    "esriFieldTypeBigInteger": "INTEGER",
    "esriFieldTypeDouble": "REAL",
    "esriFieldTypeSingle": "REAL",
}
#: A layer's ``gpkg_contents`` bounds grown to hold ``?1`` to ``?4``; NULL bounds are none yet.
_GROW_BOUNDS = (
    "UPDATE gpkg_contents SET min_x = min(coalesce(min_x, ?1), ?1), "
    "min_y = min(coalesce(min_y, ?2), ?2), max_x = max(coalesce(max_x, ?3), ?3), "
    "max_y = max(coalesce(max_y, ?4), ?4) WHERE table_name = ?5"
)
Float64 = npt.NDArray[np.float64]


_RESERVED = (
    UNDEFINED_SRS,
    SrsEntry(
        srs_id=0,
        organization="NONE",
        code=0,
        name="Undefined geographic SRS",
        definition="undefined",
    ),
    SrsEntry(
        srs_id=4326,
        organization="EPSG",
        code=4326,
        name="WGS 84",
        definition='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,'
        'AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,'
        'AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,'
        'AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]',
    ),
)

_SCHEMA = r"""
CREATE TABLE gpkg_spatial_ref_sys (srs_name TEXT NOT NULL, srs_id INTEGER PRIMARY KEY,
  organization TEXT NOT NULL, organization_coordsys_id INTEGER NOT NULL,
  definition TEXT NOT NULL, description TEXT);
CREATE TABLE gpkg_contents (table_name TEXT PRIMARY KEY, data_type TEXT NOT NULL,
  identifier TEXT UNIQUE, description TEXT DEFAULT '', last_change DATETIME NOT NULL,
  min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE, srs_id INTEGER,
  CONSTRAINT fk_gc_r_srs_id FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys(srs_id));
CREATE TABLE gpkg_geometry_columns (table_name TEXT NOT NULL, column_name TEXT NOT NULL,
  geometry_type_name TEXT NOT NULL, srs_id INTEGER NOT NULL, z TINYINT NOT NULL,
  m TINYINT NOT NULL, CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name));
CREATE TABLE gpkg_extensions (table_name TEXT, column_name TEXT, extension_name TEXT NOT NULL,
  definition TEXT NOT NULL, scope TEXT NOT NULL,
  CONSTRAINT ge_tce UNIQUE (table_name, column_name, extension_name));
"""

_STYLES = r"""
CREATE TABLE layer_styles (id INTEGER PRIMARY KEY, f_table_catalog TEXT(256),
  f_table_schema TEXT(256), f_table_name TEXT(256), f_geometry_column TEXT(256),
  styleName TEXT(30), styleQML TEXT, styleSLD TEXT, useAsDefault BOOLEAN, description TEXT,
  owner TEXT(30), ui TEXT(30), update_time DATETIME)
"""

# QGIS's categorised style on `code` (QGIS 3.28+), with the `<prop>` pairs
# older QGIS reads beside the `<Option>` ones.
_QML = r"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.28.0" styleCategories="Symbology">
<renderer-v2 type="categorizedSymbol" attr="code" symbollevels="0" enableorderby="0"
  forceraster="0">
<categories>{categories}</categories>
<symbols>{symbols}</symbols>
</renderer-v2>
</qgis>
"""
_CATEGORY = '<category value="{code}" label="{code} {name}" symbol="{i}" render="true"/>\n'
_SYMBOL = r"""<symbol type="fill" name="{i}" alpha="1" clip_to_extent="1" force_rhr="0">
<layer class="SimpleFill" pass="0" enabled="1" locked="0">
<Option type="Map">
<Option name="color" value="{rgb},255" type="QString"/>
<Option name="style" value="solid" type="QString"/>
<Option name="outline_style" value="no" type="QString"/>
</Option>
<prop k="color" v="{rgb},255"/>
<prop k="style" v="solid"/>
<prop k="outline_style" v="no"/>
</layer>
</symbol>
"""

#: The 8-byte header and the envelope, then the WKB byte order and type.
_HEADER: list[Any] = [
    ("gp", "S2"),
    ("version", "u1"),
    ("flags", "u1"),
    ("srs_id", "<i4"),
    ("envelope", "<f8", 4),
    ("byte_order", "u1"),
    ("type", "<u4"),
]


@contextmanager
def _transaction(conn: sqlite3.Connection) -> Iterator[None]:
    """One transaction around a bulk insert, unless the caller already holds one."""
    own = not conn.in_transaction
    if own:
        conn.execute("BEGIN")
    yield
    if own:
        conn.commit()


def _envelopes(xyz: Float64) -> Float64:
    """``minx, maxx, miny, maxy`` per geometry of ``xyz`` (n, points, 3)."""
    x, y = xyz[:, :, 0], xyz[:, :, 1]
    return np.column_stack([x.min(1), x.max(1), y.min(1), y.max(1)])


def _records(xyz: Float64, srs_id: int, kind: int, tail: list[Any]) -> np.ndarray:
    """One structured record per geometry: header, envelope, ISO WKB head, ``tail``."""
    rec = np.zeros(len(xyz), np.dtype([*_HEADER, *tail]))
    rec["gp"], rec["flags"], rec["srs_id"] = b"GP", 0b0011, srs_id
    rec["envelope"], rec["byte_order"], rec["type"] = _envelopes(xyz), 1, kind
    return rec


def _split(rec: np.ndarray) -> tuple[list[memoryview], Float64]:
    """Each record as one blob (a view into one buffer), and the envelopes."""
    view, width = memoryview(rec.tobytes()), rec.itemsize
    return [view[i * width : (i + 1) * width] for i in range(len(rec))], rec["envelope"].copy()


def polygon_blobs(xyz: Float64, srs_id: int) -> tuple[list[memoryview], Float64]:
    """PolygonZ blobs (type 1003, one ring of four points) for triangles ``xyz``
    of shape (n, 3, 3), and their envelopes."""
    rec = _records(xyz, srs_id, 1003, [("rings", "<u4"), ("points", "<u4"), ("xyz", "<f8", (4, 3))])
    rec["rings"], rec["points"] = 1, 4
    rec["xyz"] = np.concatenate([xyz, xyz[:, :1]], axis=1)
    return _split(rec)


def line_blobs(xyz: Float64, srs_id: int) -> tuple[list[memoryview], Float64]:
    """LineStringZ blobs (type 1002, two points) for edges ``xyz`` of shape (n, 2, 3)."""
    rec = _records(xyz, srs_id, 1002, [("points", "<u4"), ("xyz", "<f8", (2, 3))])
    rec["points"], rec["xyz"] = 2, xyz
    return _split(rec)


def _bounds(boxes: Float64) -> tuple[float, float, float, float]:
    """``min_x, min_y, max_x, max_y`` of envelopes ``minx, maxx, miny, maxy``."""
    return boxes[:, 0].min(), boxes[:, 2].min(), boxes[:, 1].max(), boxes[:, 3].max()


def create_geopackage(conn: sqlite3.Connection, srs: SrsEntry) -> None:
    """An empty GeoPackage 1.3 on ``conn``: the four ``gpkg_*`` tables, the two
    reserved SRS rows, EPSG:4326 and ``srs`` (once, if it is none of those)."""
    for pragma in ("journal_mode=OFF", "synchronous=OFF", "user_version=10300"):
        conn.execute(f"PRAGMA {pragma}")
    conn.execute(f"PRAGMA application_id={APPLICATION_ID}")
    conn.executescript(_SCHEMA)
    rows = [*_RESERVED, *([srs] if srs.srs_id not in {r.srs_id for r in _RESERVED} else [])]
    conn.executemany(
        "INSERT INTO gpkg_spatial_ref_sys VALUES (?, ?, ?, ?, ?, NULL)",
        [(r.name, r.srs_id, r.organization, r.code, r.definition) for r in rows],
    )


def register_layer(
    conn: sqlite3.Connection,
    name: str,
    kind: str,
    srs_id: int,
    bounds: Sequence[float | None],
    z: bool,
) -> None:
    """The ``gpkg_contents`` and ``gpkg_geometry_columns`` rows of a features
    table ``name`` whose geometry column ``geom`` holds ``kind`` (``POLYGON``,
    ``LINESTRING``, ``MULTIPOLYGON``); ``bounds`` is ``min_x, min_y, max_x, max_y``
    (None for an empty table)."""
    conn.execute(
        "INSERT INTO gpkg_contents VALUES (?, 'features', ?, '', ?, ?, ?, ?, ?, ?)",
        (name, name, FIXED_TIME, *bounds, srs_id),
    )
    conn.execute(
        "INSERT INTO gpkg_geometry_columns VALUES (?, 'geom', ?, ?, ?, 0)",
        (name, kind, srs_id, int(z)),
    )


def _register_attributes(conn: sqlite3.Connection, name: str) -> None:
    conn.execute(
        "INSERT INTO gpkg_contents (table_name, data_type, identifier, last_change) "
        "VALUES (?, 'attributes', ?, ?)",
        (name, name, FIXED_TIME),
    )


def create_rtree(conn: sqlite3.Connection, table: str) -> None:
    """The empty R-tree of ``table``'s ``geom``, registered as the write-only extension."""
    conn.execute(f"CREATE VIRTUAL TABLE rtree_{table}_geom USING rtree(id, minx, maxx, miny, maxy)")
    conn.execute(
        "INSERT INTO gpkg_extensions VALUES (?, 'geom', 'gpkg_rtree_index', "
        "'http://www.geopackage.org/spec120/#extension_rtree', 'write-only')",
        (table,),
    )


def fill_rtree(conn: sqlite3.Connection, table: str, ids: Iterable[int], boxes: Float64) -> None:
    """R-tree rows of ``table`` for rows ``ids``, ``boxes`` ``minx, maxx, miny, maxy`` each."""
    rows = zip(ids, *boxes.T.tolist(), strict=True)
    conn.executemany(f"INSERT INTO rtree_{table}_geom VALUES (?, ?, ?, ?, ?)", rows)


def build_rtree(conn: sqlite3.Connection, table: str, boxes: Float64, batch: int = 100_000) -> None:
    """The R-tree of ``table``'s ``geom`` from ``boxes`` (``minx, maxx, miny,
    maxy`` per row, in ``fid`` order), fed straight into the virtual table in
    batches under a 256 MB page cache; registered as the write-only extension."""
    create_rtree(conn, table)
    (cache,) = conn.execute("PRAGMA cache_size").fetchone()
    conn.execute("PRAGMA cache_size=-262144")
    with _transaction(conn):
        for start in range(0, len(boxes), batch):
            chunk = boxes[start : start + batch]
            fill_rtree(conn, table, range(start + 1, start + len(chunk) + 1), chunk)
    conn.execute(f"PRAGMA cache_size={cache}")


def _looked_up(codes: Sequence[int], palette: Palette) -> tuple[list[str | None], list[str]]:
    """Each code's name (None when the palette lacks it) and colour (magenta then)."""
    found = [palette.get(code) for code in codes]
    return [f[0] if f else None for f in found], [f[1] if f else UNNAMED_COLOUR for f in found]


def write_triangles(
    conn: sqlite3.Connection,
    vertices: Any,
    triangles: Any,
    codes: Any | None,
    palette: Palette | None,
    srs_id: int,
    batch: int = 100_000,
) -> None:
    """The ``triangles`` layer: one PolygonZ row per triangle, in triangle order
    (``fid - 1`` is the index), with ``code``, ``land_cover``, ``colour``
    (NULL on every row when ``codes`` is None), the three corner heights, the
    plan area and the slope; registered and indexed. Rows and blobs are built
    per ``batch`` from array slices, in one transaction."""
    points = np.asarray(vertices, dtype=np.float64)
    faces = np.asarray(triangles)
    conn.execute(
        "CREATE TABLE triangles (fid INTEGER PRIMARY KEY, geom POLYGON, code INTEGER, "
        "land_cover TEXT, colour TEXT, z_min REAL, z_max REAL, z_mean REAL, area_m2 REAL, "
        "slope_deg REAL)"
    )
    boxes = np.empty((len(faces), 4))
    with _transaction(conn):
        for start in range(0, len(faces), batch):
            stop = min(start + batch, len(faces))
            xyz = points[faces[start:stop]]
            blobs, boxes[start:stop] = polygon_blobs(xyz, srs_id)
            e1, e2 = xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]
            normal = np.cross(e1, e2)
            area = 0.5 * np.abs(normal[:, 2])
            slope = np.arctan2(np.hypot(normal[:, 0], normal[:, 1]), np.abs(normal[:, 2]))
            z = xyz[:, :, 2]
            if codes is None:
                code: list[Any] = [None] * (stop - start)
                name, colour = code, code
            else:
                code = np.asarray(codes)[start:stop].tolist()
                name, colour = _looked_up(code, palette or {})
            rows = zip(
                range(start + 1, stop + 1),
                blobs,
                code,
                name,
                colour,
                z.min(1).tolist(),
                z.max(1).tolist(),
                z.mean(1).tolist(),
                area.tolist(),
                np.degrees(slope).tolist(),
                strict=True,
            )
            conn.executemany("INSERT INTO triangles VALUES (?,?,?,?,?,?,?,?,?,?)", rows)
    register_layer(conn, "triangles", "POLYGON", srs_id, _bounds(boxes), z=True)
    build_rtree(conn, "triangles", boxes, batch)


def write_cells(
    conn: sqlite3.Connection, cells: CellTable, srs_id: int, catchment_ids: Any = None
) -> None:
    """The ``cells`` layer: one PointZ row per triangle at its centroid, in
    triangle order (``fid - 1`` is the index), its ``catchment_ids`` entry (1
    without them), the plan area, the
    height, slope, aspect (NULL where NaN) and the shares; registered and indexed."""
    xyz = cells.midpoints[:, None, :]
    rec = _records(xyz, srs_id, 1001, [("xyz", "<f8", 3)])
    rec["xyz"] = cells.midpoints
    blobs, boxes = _split(rec)
    names = "".join(f", {name} REAL" for name in SHARES)
    conn.execute(
        "CREATE TABLE cells (fid INTEGER PRIMARY KEY, geom POINT, catchment_id INTEGER, "
        f"area_m2 REAL, z_m REAL, slope_deg REAL, aspect_deg REAL{names})"
    )
    aspect = [None if np.isnan(a) else a for a in cells.aspect_deg.tolist()]
    columns = [cells.area_m2, cells.midpoints[:, 2], cells.slope_deg]
    rows = zip(
        range(1, len(blobs) + 1), blobs,
        [1] * len(blobs) if catchment_ids is None else np.asarray(catchment_ids).tolist(),
        *(c.tolist() for c in columns),
        aspect, *cells.shares.T.tolist(), strict=True,
    )  # fmt: skip
    with _transaction(conn):
        conn.executemany(f"INSERT INTO cells VALUES ({', '.join('?' * 13)})", rows)
    register_layer(conn, "cells", "POINT", srs_id, _bounds(boxes), z=True)
    build_rtree(conn, "cells", boxes)


def write_fractions(conn: sqlite3.Connection, fractions: ClassFractions) -> None:
    """The ``land_cover_fractions`` attributes table: one row per (triangle,
    class) entry, ``triangle`` the triangle's ``fid``, indexed on it."""
    rows = np.repeat(np.arange(1, fractions.triangles + 1), np.diff(fractions.offsets))
    conn.execute(
        "CREATE TABLE land_cover_fractions (fid INTEGER PRIMARY KEY, triangle INTEGER, "
        "code INTEGER, fraction REAL)"
    )
    with _transaction(conn):
        conn.executemany(
            "INSERT INTO land_cover_fractions VALUES (NULL, ?, ?, ?)",
            zip(
                rows.tolist(), fractions.classes.tolist(), fractions.fractions.tolist(), strict=True
            ),
        )
    conn.execute("CREATE INDEX land_cover_fractions_triangle ON land_cover_fractions (triangle)")
    _register_attributes(conn, "land_cover_fractions")


def write_lines(
    conn: sqlite3.Connection,
    vertices: Any,
    edges: Any,
    edge_masks: Any,
    vocabulary: EdgeVocabulary,
    srs_id: int,
) -> None:
    """The ``lines`` layer: one LineStringZ row per constraint edge, in edge
    order, with its plan length and one 0/1 column per vocabulary property set
    on at least one edge, named as ``write_vtk`` names its cell arrays.
    Registered; no R-tree."""
    xyz = np.asarray(vertices, dtype=np.float64)[np.asarray(edges)]
    masks = np.asarray(edge_masks).astype(np.int64)
    used = {name for mask in np.unique(masks) for name in vocabulary.names(int(mask))}
    names = [name for _, name in vocabulary.table() if name in used]
    blobs, boxes = line_blobs(xyz, srs_id)
    length = np.hypot(*(xyz[:, 1, :2] - xyz[:, 0, :2]).T)
    flags = [((masks & vocabulary.mask(name)) != 0).astype(int).tolist() for name in names]
    columns = "".join(f', "{name}" INTEGER' for name in names)
    conn.execute(
        f"CREATE TABLE lines (fid INTEGER PRIMARY KEY, geom LINESTRING, length_m REAL{columns})"
    )
    marks = ", ".join("?" * (3 + len(names)))
    with _transaction(conn):
        rows = zip(range(1, len(xyz) + 1), blobs, length.tolist(), *flags, strict=True)
        conn.executemany(f"INSERT INTO lines VALUES ({marks})", rows)
    register_layer(conn, "lines", "LINESTRING", srs_id, _bounds(boxes), z=True)


def multipolygon_blobs(
    geometries: Sequence[BaseGeometry], srs_id: int
) -> tuple[list[bytes], Float64]:
    """2-D MultiPolygon blobs of polygons or multipolygons, a polygon wrapped
    as a one-part multipolygon, and their envelopes ``minx, maxx, miny, maxy``,
    as flag ``0b0011`` declares (12-128 section 2.1.3)."""
    multi = [g if g.geom_type == "MultiPolygon" else shapely.MultiPolygon([g]) for g in geometries]
    boxes = shapely.bounds(multi).reshape(-1, 4)[:, [0, 2, 1, 3]]
    wkb = shapely.to_wkb(multi, byte_order=1, flavor="iso", output_dimension=2)
    head = struct.Struct("<2sBBi4d")
    blobs = [
        head.pack(b"GP", 0, 0b0011, srs_id, *box) + bytes(w)
        for box, w in zip(boxes, wkb, strict=True)
    ]
    return blobs, boxes


def create_features_table(
    conn: sqlite3.Connection,
    name: str,
    srs_id: int,
    kind: str,
    pk: str,
    columns: Sequence[tuple[str, str]],
) -> None:
    """An empty features table ``name``: integer key ``pk``, ``geom`` of
    ``kind``, and one column per ``(name, esri field type)``; registered and
    with its empty R-tree. A field type with no column type is a `ValueError`."""
    for column, esri in columns:
        if esri not in _COLUMN_TYPES:
            raise ValueError(f"{column}: the field type {esri} has no GeoPackage column type")
    rest = "".join(f', "{c}" {_COLUMN_TYPES[t]}' for c, t in columns)
    conn.execute(f'CREATE TABLE "{name}" ("{pk}" INTEGER PRIMARY KEY, geom {kind}{rest})')
    register_layer(conn, name, kind, srs_id, (None,) * 4, z=False)
    create_rtree(conn, name)


def insert_features(
    conn: sqlite3.Connection,
    name: str,
    srs_id: int,
    rows: Sequence[tuple[int, BaseGeometry, Sequence[Any]]],
) -> None:
    """``(key, polygon, values)`` rows into the table ``name`` as
    multipolygons, their R-tree entries keyed by ``key``, and the layer's
    bounds grown to hold them: one transaction."""
    if not rows:
        return
    blobs, boxes = multipolygon_blobs([g for _, g, _ in rows], srs_id)
    marks = ", ".join("?" * (2 + len(rows[0][2])))
    with _transaction(conn):
        values = ((k, b, *v) for (k, _, v), b in zip(rows, blobs, strict=True))
        conn.executemany(f'INSERT INTO "{name}" VALUES ({marks})', values)
        fill_rtree(conn, name, [k for k, _, _ in rows], boxes)
        conn.execute(_GROW_BOUNDS, (*_bounds(boxes), name))


def write_land_cover(
    conn: sqlite3.Connection,
    polygons: Sequence[tuple[int, BaseGeometry]],
    palette: Palette | None,
    srs_id: int,
) -> None:
    """The ``land_cover`` layer: one 2-D MultiPolygon row per ``(code,
    polygon)``, with its name, colour and plan area. The header envelope is
    ``minx, maxx, miny, maxy``, as flag ``0b0011`` declares (12-128 section
    2.1.3). Registered; no R-tree."""
    blobs, boxes = multipolygon_blobs([g for _, g in polygons], srs_id)
    codes = [int(code) for code, _ in polygons]
    name, colour = _looked_up(codes, palette or {})
    conn.execute(
        "CREATE TABLE land_cover (fid INTEGER PRIMARY KEY, geom MULTIPOLYGON, code INTEGER, "
        "land_cover TEXT, colour TEXT, area_m2 REAL)"
    )
    areas = shapely.area([g for _, g in polygons]).tolist()
    rows = zip(blobs, codes, name, colour, areas, strict=True)
    with _transaction(conn):
        conn.executemany("INSERT INTO land_cover VALUES (NULL, ?, ?, ?, ?, ?)", rows)
    register_layer(conn, "land_cover", "MULTIPOLYGON", srs_id, _bounds(boxes), z=False)


def write_classes(conn: sqlite3.Connection) -> None:
    """The ``land_cover_classes`` attributes table: one row per code in
    ``triangles`` with its name, colour, triangle count and plan area."""
    conn.execute(
        "CREATE TABLE land_cover_classes (fid INTEGER PRIMARY KEY, code INTEGER, "
        "land_cover TEXT, colour TEXT, triangles INTEGER, area_m2 REAL)"
    )
    conn.execute(
        "INSERT INTO land_cover_classes (code, land_cover, colour, triangles, area_m2) "
        "SELECT code, land_cover, colour, count(*), sum(area_m2) FROM triangles "
        "GROUP BY code ORDER BY code"
    )
    _register_attributes(conn, "land_cover_classes")


def write_run(conn: sqlite3.Connection, rows: Sequence[tuple[str, str, str, float | None]]) -> None:
    """The ``run`` attributes table: one row per ``(name, wording, value,
    number)`` of the run record, in the given order; ``name`` is unique."""
    conn.execute(
        "CREATE TABLE run (fid INTEGER PRIMARY KEY, name TEXT UNIQUE, wording TEXT, "
        "value TEXT, number REAL)"
    )
    with _transaction(conn):
        conn.executemany("INSERT INTO run VALUES (NULL, ?, ?, ?, ?)", rows)
    _register_attributes(conn, "run")


def _rgb(colour: str) -> str:
    return ",".join(str(int(colour[i : i + 2], 16)) for i in (1, 3, 5))


def write_styles(conn: sqlite3.Connection, palette: Palette, layers: Sequence[str]) -> None:
    """QGIS's ``layer_styles``: for each of ``layers`` the default style, a
    categorised fill on ``code`` with the palette's colours and labels."""
    conn.execute(_STYLES)
    codes = sorted(palette)
    categories = "".join(
        _CATEGORY.format(code=code, name=escape(palette[code][0], {'"': "&quot;"}), i=i)
        for i, code in enumerate(codes)
    )
    symbols = "".join(
        _SYMBOL.format(i=i, rgb=_rgb(palette[code][1])) for i, code in enumerate(codes)
    )
    qml = _QML.format(categories=categories, symbols=symbols)
    conn.executemany(
        "INSERT INTO layer_styles (f_table_catalog, f_table_schema, f_table_name, "
        "f_geometry_column, styleName, styleQML, styleSLD, useAsDefault, description, owner, "
        "ui, update_time) VALUES ('', '', ?, 'geom', 'land cover', ?, '', 1, "
        "'land-cover colours, categorised on code', '', '', ?)",
        [(layer, qml, FIXED_TIME) for layer in layers],
    )
