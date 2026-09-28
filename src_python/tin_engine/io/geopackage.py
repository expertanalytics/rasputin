"""GeoPackage features, read without GDAL (increment 16b R2, R3).

`docs/increments/16b-terrain-polygons.md` R3; the blob header is OGC 12-128r18
section 2.1.3, the R-tree Annex F.3. The decoder's stream is an open
``sqlite3.Connection`` (``io.repository.open_geopackage``): this module opens
nothing, knows no path and no CRS object, and imports no ``_core``. The CRS it
reports is text for ``crs.parse_crs``.

Rows come in primary-key order, whatever order SQLite or the R-tree holds
them in, because the mesh depends on input order (R7).
"""

from __future__ import annotations

import sqlite3
import struct
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Any

import shapely
from shapely.geometry.base import BaseGeometry

#: `PRAGMA application_id` of a GeoPackage: "GPKG".
APPLICATION_ID = 0x47504B47
#: Envelope code -> bytes after the 8-byte header (R3); 5-7 are refused.
_ENVELOPE_BYTES = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}

Box = tuple[float, float, float, float]


class GeoPackageError(ValueError):
    """A GeoPackage, layer or row this reader refuses, in words."""


@dataclass(frozen=True, slots=True)
class GpkgLayer:
    """One features table: its geometry column, primary key, SRS and R-tree."""

    table: str
    column: str
    pk: str
    srs_id: int
    crs: str
    rtree: str | None


@dataclass(frozen=True, slots=True)
class RawFeature:
    """One row: its primary key, geometry (None for the empty flag) and value."""

    fid: int
    geometry: BaseGeometry | None
    value: Any


def _q(name: str) -> str:
    """An SQL identifier, quoted."""
    return '"' + name.replace('"', '""') + '"'


def decode_geometry(blob: bytes) -> tuple[int, BaseGeometry | None]:
    """A geometry blob's ``srs_id`` and geometry; None for the empty flag."""
    if len(blob) < 8 or blob[:2] != b"GP" or blob[2] != 0:
        raise GeoPackageError("not a GeoPackage geometry blob (magic GP, version 0)")
    flags = blob[3]
    envelope = (flags >> 1) & 0b111
    if envelope not in _ENVELOPE_BYTES:
        raise GeoPackageError(f"envelope code {envelope} is not defined")
    if flags & 0b10_0000:
        raise GeoPackageError("an extended geometry type is not standard WKB")
    (srs_id,) = struct.unpack("<i" if flags & 1 else ">i", blob[4:8])
    start = 8 + _ENVELOPE_BYTES[envelope]
    if len(blob) < start:
        raise GeoPackageError(f"the blob ends inside its {start - 8}-byte envelope")
    if flags & 0b1_0000:
        return srs_id, None
    try:
        return srs_id, shapely.from_wkb(bytes(blob[start:]))
    except (shapely.errors.GEOSException, ValueError) as exc:
        raise GeoPackageError(f"unreadable WKB: {exc}") from exc


def layer_info(conn: sqlite3.Connection, table: str | None) -> GpkgLayer:
    """The layer ``table``, or the file's only features table for None."""
    try:
        if conn.execute("PRAGMA application_id").fetchone()[0] != APPLICATION_ID:
            raise GeoPackageError("not a GeoPackage (application_id is not GPKG)")
        contents = "SELECT table_name, data_type FROM gpkg_contents"
        kinds = dict(conn.execute(contents).fetchall())
        if table is None:
            found = sorted(t for t, kind in kinds.items() if kind == "features")
            if len(found) != 1:
                raise GeoPackageError(f"name one layer of the features tables {found}")
            table = found[0]
        if table not in kinds:
            raise GeoPackageError(f"no layer {table!r}; the file has {sorted(kinds)}")
        if kinds[table] != "features":
            raise GeoPackageError(f"layer {table!r} holds {kinds[table]}, not features")
        column, srs_id = conn.execute(
            "SELECT column_name, srs_id FROM gpkg_geometry_columns WHERE table_name = ?", (table,)
        ).fetchone()
        org, code, definition = conn.execute(
            "SELECT organization, organization_coordsys_id, definition "
            "FROM gpkg_spatial_ref_sys WHERE srs_id = ?",
            (srs_id,),
        ).fetchone()
        info = conn.execute(f"PRAGMA table_info({_q(table)})").fetchall()
        pk = next((row[1] for row in info if row[5] == 1), "rowid")
        rtree = f"rtree_{table}_{column}"
        indexed = (
            conn.execute(
                "SELECT 1 FROM gpkg_extensions WHERE table_name = ? AND column_name = ? "
                "AND extension_name = 'gpkg_rtree_index'",
                (table, column),
            ).fetchone()
            and conn.execute("SELECT 1 FROM sqlite_master WHERE name = ?", (rtree,)).fetchone()
        )
    except (sqlite3.Error, TypeError) as exc:
        raise GeoPackageError(f"not a readable GeoPackage: {exc}") from exc
    crs = f"EPSG:{code}" if str(org).upper() == "EPSG" else str(definition)
    return GpkgLayer(
        table=table, column=column, pk=pk, srs_id=srs_id, crs=crs, rtree=rtree if indexed else None
    )


def query_features(
    conn: sqlite3.Connection,
    layer: GpkgLayer,
    box: Box | None,
    attribute: str,
    scale: float = 100_000.0,
) -> Iterator[RawFeature]:
    """Every row whose R-tree box, widened by ``s * max(1, s / scale)`` with
    ``s`` its width plus height (R5, "Long edges"), meets the closed ``box``
    ``(minx, miny, maxx, maxy)``; every row with no box or no R-tree. A
    superset of the rows meeting ``box``, in primary-key order."""
    columns = {
        row[1].lower(): row[1] for row in conn.execute(f"PRAGMA table_info({_q(layer.table)})")
    }
    if attribute.lower() not in columns:
        raise GeoPackageError(f"layer {layer.table!r} has no column {attribute!r}")
    pk, table = _q(layer.pk), _q(layer.table)
    select = (
        f"SELECT t.{pk}, t.{_q(layer.column)}, t.{_q(columns[attribute.lower()])} FROM {table} t"
    )
    params: dict[str, float] = {}
    if box is not None and layer.rtree is not None:
        # Fetch by primary key from the R-tree's hits: a join is planned as a
        # full scan of the table with one R-tree probe per row.
        s = "(r.maxx - r.minx + r.maxy - r.miny)"
        w = f"({s} * max(1.0, {s} / :scale))"
        select += (
            f" WHERE t.{pk} IN (SELECT r.id FROM {_q(layer.rtree)} r WHERE r.maxx + {w} >= :minx"
            f" AND r.minx - {w} <= :maxx AND r.maxy + {w} >= :miny AND r.miny - {w} <= :maxy)"
        )
        params = dict(zip(("minx", "miny", "maxx", "maxy"), box, strict=True), scale=scale)
    try:
        rows = conn.execute(f"{select} ORDER BY t.{pk}", params)
    except sqlite3.OperationalError as exc:
        if "rtree" in str(exc).lower():
            raise GeoPackageError(
                f"{layer.table}: this SQLite has no R-tree module ({exc})"
            ) from exc
        raise
    for fid, data, value in rows:
        try:
            srs_id, geometry = (layer.srs_id, None) if data is None else decode_geometry(data)
        except GeoPackageError as exc:
            raise GeoPackageError(f"layer {layer.table!r} row {fid}: {exc}") from exc
        if srs_id != layer.srs_id:
            raise GeoPackageError(
                f"layer {layer.table!r} row {fid}: srs_id {srs_id}, the layer's is {layer.srs_id}"
            )
        yield RawFeature(fid=fid, geometry=geometry, value=value)
