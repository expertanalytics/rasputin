"""Cut increment 16b's CORINE extract out of Ola's CORINE 2018 GeoPackage.

    python tests/fixtures/corine/extract.py \
        ../rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg

Test support, not production code (`docs/increments/16b-terrain-polygons.md`,
"Test data"), as `tests/fixtures/dtm10/extract.py` is. Writes
`clc2018_7908_3.gpkg` next to this file: the features of the Europe layer
`U2018_CLC2018_V2020_20u1` whose R-tree box meets the committed benchmark
tile's (`tests/fixtures/dem_archive/7908_3_10m_z33.tif`) image in EPSG:3035,
each clipped in EPSG:3035 to that image buffered by 1 km, so the partition is
exact over the tile and every artificial edge the clip makes lies at least 1 km
outside it.

The output keeps the source's layout: the table name, its columns (`OBJECTID`
the primary key, `Shape`, `Code_18`, `Remark`, `Area_Ha`, `ID`), the SRS row
for 3035, an R-tree `rtree_U2018_CLC2018_V2020_20u1_Shape` listed in
`gpkg_extensions`, `application_id` GPKG and `user_version` 10200, and blobs in
the source's header layout: flags `0x01` (little-endian, no envelope, not
empty), then little-endian WKB. A clipped geometry is written as a
`MULTIPOLYGON`, like every source row. Rows are written in `OBJECTID` order.

Written GDAL-free with `sqlite3`, `shapely` and `pyproj`. See `NOTICE` for the
attribution the data needs.
"""

from __future__ import annotations

import sqlite3
import struct
import sys
from contextlib import closing
from pathlib import Path

import numpy as np
import shapely
from pyproj import Transformer
from shapely.geometry import MultiPolygon, Polygon, box

HERE = Path(__file__).resolve().parent
OUT = HERE / "clc2018_7908_3.gpkg"
TABLE = "U2018_CLC2018_V2020_20u1"
RTREE = f"rtree_{TABLE}_Shape"
SRS = 3035
# The benchmark tile's node rectangle (7908_3_10m_z33.tif: x_min 799 750,
# y_max 7 950 250, 5051 x 5051 nodes at 10 m), EPSG:25833.
TILE = (799_750.0, 7_899_750.0, 850_250.0, 7_950_250.0)
BUFFER = 1_000.0


def tile_image() -> Polygon:
    """The tile's rectangle in EPSG:3035, its edges densified to 100 m first,
    buffered by 1 km."""
    ring = np.asarray(shapely.segmentize(box(*TILE), 100.0).exterior.coords)
    to_3035 = Transformer.from_crs(25833, SRS, always_xy=True)
    x, y = to_3035.transform(ring[:, 0], ring[:, 1])
    return Polygon(np.column_stack([x, y])).buffer(BUFFER)


def decode(blob: bytes) -> shapely.Geometry:
    """The source's own header layout (M1): 8 bytes, no envelope."""
    assert blob[:2] == b"GP" and blob[3] == 0x01, blob[:8]
    return shapely.from_wkb(blob[8:])


def encode(geometry: MultiPolygon) -> bytes:
    return (
        b"GP" + bytes([0, 0x01]) + struct.pack("<i", SRS) + shapely.to_wkb(geometry, byte_order=1)
    )


def as_multipolygon(geometry: shapely.Geometry) -> MultiPolygon | None:
    parts = [g for g in shapely.get_parts(geometry) if isinstance(g, Polygon) and not g.is_empty]
    return MultiPolygon(parts) if parts else None


def create(con: sqlite3.Connection, srs_row: tuple[object, ...]) -> None:
    con.executescript(
        f"""
        PRAGMA application_id = 1196444487;
        PRAGMA user_version = 10200;
        CREATE TABLE gpkg_spatial_ref_sys (srs_name TEXT NOT NULL, srs_id INTEGER PRIMARY KEY,
          organization TEXT NOT NULL, organization_coordsys_id INTEGER NOT NULL,
          definition TEXT NOT NULL, description TEXT);
        CREATE TABLE gpkg_contents (table_name TEXT NOT NULL PRIMARY KEY, data_type TEXT NOT NULL,
          identifier TEXT UNIQUE, description TEXT DEFAULT '',
          last_change DATETIME NOT NULL DEFAULT '2026-09-28T00:00:00.000Z',
          min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE, srs_id INTEGER);
        CREATE TABLE gpkg_geometry_columns (table_name TEXT NOT NULL, column_name TEXT NOT NULL,
          geometry_type_name TEXT NOT NULL, srs_id INTEGER NOT NULL, z TINYINT NOT NULL,
          m TINYINT NOT NULL, CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name));
        CREATE TABLE gpkg_extensions (table_name TEXT, column_name TEXT,
          extension_name TEXT NOT NULL, definition TEXT NOT NULL, scope TEXT NOT NULL,
          CONSTRAINT ge_tce UNIQUE (table_name, column_name, extension_name));
        CREATE TABLE {TABLE} (OBJECTID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
          Shape MULTIPOLYGON, Code_18 TEXT(3), Remark TEXT(20), Area_Ha DOUBLE, ID TEXT(18));
        CREATE VIRTUAL TABLE {RTREE} USING rtree(id, minx, maxx, miny, maxy);
        """
    )
    con.executemany(
        "INSERT INTO gpkg_spatial_ref_sys VALUES (?,?,?,?,?,?)",
        [
            ("Undefined Cartesian", -1, "NONE", -1, "undefined", None),
            ("Undefined Geographic", 0, "NONE", 0, "undefined", None),
            srs_row,
        ],
    )
    con.execute(
        "INSERT INTO gpkg_contents (table_name, data_type, identifier, srs_id) "
        "VALUES (?, 'features', ?, ?)",
        (TABLE, TABLE, SRS),
    )
    con.execute(
        "INSERT INTO gpkg_geometry_columns VALUES (?, 'Shape', 'MULTIPOLYGON', ?, 0, 0)",
        (TABLE, SRS),
    )
    con.execute(
        "INSERT INTO gpkg_extensions VALUES (?, 'Shape', 'gpkg_rtree_index', "
        "'F.3 RTree Spatial Index', 'write-only')",
        (TABLE,),
    )


def main(source: Path) -> None:
    image = tile_image()
    minx, miny, maxx, maxy = image.bounds
    with closing(sqlite3.connect(source.resolve().as_uri() + "?mode=ro", uri=True)) as src:
        srs_row = src.execute(
            "SELECT * FROM gpkg_spatial_ref_sys WHERE srs_id = ?", (SRS,)
        ).fetchone()
        rows = src.execute(
            f"SELECT t.OBJECTID, t.Shape, t.Code_18, t.Remark, t.Area_Ha, t.ID FROM {TABLE} t "
            f"JOIN {RTREE} r ON r.id = t.OBJECTID "
            "WHERE r.maxx >= ? AND r.minx <= ? AND r.maxy >= ? AND r.miny <= ? "
            "ORDER BY t.OBJECTID",
            (minx, maxx, miny, maxy),
        ).fetchall()
    kept = []
    for objectid, blob, code, remark, area, ident in rows:
        clipped = as_multipolygon(shapely.intersection(decode(blob), image))
        if clipped is not None:
            kept.append((objectid, clipped, code, remark, area, ident))
    OUT.unlink(missing_ok=True)
    with closing(sqlite3.connect(OUT)) as con:
        create(con, srs_row)
        for objectid, geometry, code, remark, area, ident in kept:
            con.execute(
                f"INSERT INTO {TABLE} VALUES (?,?,?,?,?,?)",
                (objectid, encode(geometry), code, remark, area, ident),
            )
            gx0, gy0, gx1, gy1 = geometry.bounds
            con.execute(f"INSERT INTO {RTREE} VALUES (?,?,?,?,?)", (objectid, gx0, gx1, gy0, gy1))
        con.commit()
        con.execute("VACUUM")
    vertices = sum(len(shapely.get_coordinates(g)) for _, g, *_ in kept)
    classes = sorted({code for _, _, code, *_ in kept})
    print(
        f"{OUT.name}: {len(kept)} of {len(rows)} candidates, {vertices} vertices, "
        f"classes {' '.join(classes)}, {OUT.stat().st_size} bytes"
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]))
