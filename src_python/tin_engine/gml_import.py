"""A GML 3.2 document of surfaces into one GeoPackage table, once.

`gml_to_geopackage` streams the features `io.gml.iter_gml32` decodes into the
connection it is given, through the writer `io.geopackage_write` (the R-tree
and the layer bounds grow with each chunk). It opens no file and knows no
path; it is blocking, so an async caller runs it in `asyncio.to_thread`.
Nothing is reprojected: the table keeps the document's CRS.
"""

from __future__ import annotations

import re
import sqlite3
from dataclasses import dataclass
from typing import IO, Any

from shapely.geometry.base import BaseGeometry

from tin_engine.crs import parse_crs, same_crs, srs_entry
from tin_engine.io.geopackage_write import (
    create_features_table,
    create_geopackage,
    insert_features,
)
from tin_engine.io.gml import Gml32Feature, GmlError, iter_gml32

#: Features per transaction, a chosen value: 1000 kept the peak resident
#: memory at 96 MiB on 16,140 features and 126 MiB on 287,730.
CHUNK = 1000

_INTEGER = re.compile(r"[+-]?[0-9]+")
_RESERVED = ("fid", "geom", "gml_id")


@dataclass(frozen=True, slots=True)
class GmlImport:
    """What a conversion wrote: the table (the feature type), the CRS as
    stored (the ``srsName`` text, or the given one), the rows written, the
    features skipped for having no surface, and the columns after ``fid`` and
    ``geom``."""

    layer: str
    crs: str
    features: int
    without_surface: int
    columns: tuple[str, ...]


def gml_to_geopackage(
    stream: IO[bytes], conn: sqlite3.Connection, *, crs: str | None = None
) -> GmlImport:
    """Every feature with a surface in ``stream`` as a row of one table on
    ``conn``, a new GeoPackage. The first feature with a surface names the
    table, its columns (``gml_id``, then its values: INTEGER when the text is
    a decimal integer, else TEXT) and the CRS (its ``srsName``, else ``crs``).
    Refused, a `GmlError`: a later feature of another type, or with a value
    the columns lack; a value named ``fid``, ``geom`` or ``gml_id``; a second
    ``srsName``; a CRS that does not parse, lists latitude first, or differs
    from ``crs``; a document with no surface at all."""
    first: Gml32Feature | None = None
    fixed = ""
    columns: tuple[str, ...] = ()
    srs_id = written = skipped = 0
    rows: list[tuple[int, BaseGeometry, list[Any]]] = []
    for feature in iter_gml32(stream):
        if feature.geometry is None:
            skipped += 1
            continue
        if first is None:
            first, fixed = feature, _fixed_crs(feature, crs)
            columns = ("gml_id", *_columns(feature))
            srs_id = _create(conn, feature, fixed)
        _check(feature, first, fixed, columns)
        written += 1
        row = [feature.gml_id, *map(feature.values.get, columns[1:])]
        rows.append((written, feature.geometry, row))
        if len(rows) == CHUNK:
            insert_features(conn, first.type_name, srs_id, rows)
            rows = []
    if first is None:
        raise GmlError("no feature has a surface: nothing to write")
    insert_features(conn, first.type_name, srs_id, rows)
    conn.commit()
    return GmlImport(first.type_name, fixed, written, skipped, columns)


def _columns(feature: Gml32Feature) -> tuple[str, ...]:
    for name in feature.values:
        if name in _RESERVED:
            raise GmlError(f"feature {feature.gml_id}: a value named {name}, a reserved column")
    return tuple(feature.values)


def _create(conn: sqlite3.Connection, feature: Gml32Feature, crs: str) -> int:
    """The GeoPackage and its one table, typed by ``feature``'s values; its srs_id."""
    srs = srs_entry(crs)
    create_geopackage(conn, srs)
    typed = [
        (k, "esriFieldTypeInteger" if _INTEGER.fullmatch(v or "") else "esriFieldTypeString")
        for k, v in feature.values.items()
    ]
    columns = [("gml_id", "esriFieldTypeString"), *typed]
    create_features_table(conn, feature.type_name, srs.srs_id, "MULTIPOLYGON", "fid", columns)
    conn.commit()
    return srs.srs_id


def _fixed_crs(feature: Gml32Feature, given: str | None) -> str:
    """The CRS the table takes from its first feature: its ``srsName``, else ``given``."""
    name = feature.srs_name or given
    if name is None:
        raise GmlError(f"feature {feature.gml_id} names no srsName; give --crs")
    try:
        parsed = parse_crs(name)
    except ValueError as exc:
        raise GmlError(f"feature {feature.gml_id}: {exc}") from exc
    if parsed.axis_info and parsed.axis_info[0].direction != "east":
        raise GmlError(
            f"feature {feature.gml_id}: {name} lists {parsed.axis_info[0].name} first;"
            " only a CRS whose first axis points east is read"
        )
    if feature.srs_name and given is not None and not same_crs(feature.srs_name, given):
        raise GmlError(f"feature {feature.gml_id}: its srsName {feature.srs_name} is not {given}")
    return name


def _check(feature: Gml32Feature, first: Gml32Feature, crs: str, columns: tuple[str, ...]) -> None:
    where = f"feature {feature.gml_id}"
    if feature.type_name != first.type_name:
        raise GmlError(f"{where} is a {feature.type_name}, not a {first.type_name}")
    if feature.srs_name is not None and feature is not first and feature.srs_name != crs:
        raise GmlError(f"{where}: srsName {feature.srs_name}, after {crs}")
    extra = [name for name in feature.values if name not in columns[1:]]
    if extra:
        raise GmlError(f"{where}: {', '.join(extra)}, not a column of {first.type_name}")
