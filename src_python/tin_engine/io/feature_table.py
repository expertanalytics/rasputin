"""A whole feature layer, every attribute column with its field type.

`docs/increments/116-feature-table.md`, *2. Design*. A GeoPackage's columns
take their field type from the declared type (OGC 12-128r18, Table 1); a
GeoJSON file's from its values (RFC 7946). Rows are read by
``feature_input.read_source``, so its refusals and box query hold here too.

Blocking (sqlite3, GEOS): an async caller runs the reader in
``asyncio.to_thread``.
"""

from __future__ import annotations

import re
import sqlite3
from collections.abc import Callable, Sequence
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from shapely.geometry.base import BaseGeometry

from tin_engine.feature_input import FeatureError, read_source
from tin_engine.io.geopackage import layer_info
from tin_engine.io.repository import open_geopackage, read_json

INTEGER, DOUBLE, STRING = "esriFieldTypeInteger", "esriFieldTypeDouble", "esriFieldTypeString"
#: A GeoPackage declared column type -> its field type; `TEXT(n)` is `TEXT`.
DECLARED = dict.fromkeys(("INTEGER", "INT", "MEDIUMINT", "SMALLINT", "TINYINT", "BOOLEAN"), INTEGER)
DECLARED |= dict.fromkeys(("REAL", "DOUBLE", "FLOAT"), DOUBLE)
DECLARED |= dict.fromkeys(("TEXT", "DATE", "DATETIME"), STRING)


@dataclass(frozen=True, slots=True)
class FeatureTable:
    """A feature layer as read: its CRS text, its GeoPackage layer (None for
    GeoJSON), its columns as (name, field type) in file order, and per row
    (fid, geometry or None, values in column order), in file order."""

    crs: str
    layer: str | None
    columns: tuple[tuple[str, str], ...]
    rows: tuple[tuple[int, BaseGeometry | None, tuple[Any, ...]], ...]


def read_feature_table(
    path: Path,
    layer: str | None,
    box_for: Callable[[str], tuple[float, float, float, float]] | None = None,
) -> FeatureTable:
    """The layer at ``path`` by suffix (`.gpkg`, `.geojson`, `.json`); with
    ``box_for(own CRS)``, a GeoPackage's rows its R-tree finds near that box.
    Every refusal is a :class:`FeatureError` naming the file."""
    suffix = path.suffix.lower()
    if suffix == ".gpkg":
        try:
            with closing(open_geopackage(path)) as conn:
                info = layer_info(conn, layer)
                declared = "SELECT name, type FROM pragma_table_info(?)"
                found = [
                    (n, t)
                    for n, t in conn.execute(declared, (info.table,))
                    if n not in (info.pk, info.column)
                ]
        except (ValueError, sqlite3.Error) as exc:
            raise FeatureError(f"{path.name}: {exc}") from exc
        columns = tuple((n, _declared(path, n, t)) for n, t in found)
        names = tuple(n for n, _ in columns)
        source = read_source(path, info.table, None, box_for or (lambda _: None), extra=names)
        values = source.extra or tuple(() for _ in source.rows)  # () when no column is asked
        rows = tuple((int(f), g, v) for (f, g, _), v in zip(source.rows, values, strict=True))
        return FeatureTable(source.crs, info.table, columns, rows)
    if suffix not in (".geojson", ".json"):
        raise FeatureError(f"{path.name}: unknown suffix {suffix or '(none)'}; use .gpkg, .geojson")
    try:
        doc = read_json(path)
        features = doc["features"] if doc.get("type") == "FeatureCollection" else [doc]
        keys = tuple(dict.fromkeys(k for f in features for k in f.get("properties") or {}))
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise FeatureError(f"{path.name}: not a GeoJSON FeatureCollection ({exc})") from exc
    source = read_source(path, layer, None, lambda _: None, extra=keys)
    values = source.extra or tuple(() for _ in source.rows)  # () when no key is asked
    columns = tuple((k, _inferred(path, k, [v[i] for v in values])) for i, k in enumerate(keys))
    rows = tuple(
        (fid, g, tuple(int(x) if isinstance(x, bool) else x for x in v))
        for fid, ((_, g, _), v) in enumerate(zip(source.rows, values, strict=True), 1)
    )
    return FeatureTable(source.crs, None, columns, rows)


def _declared(path: Path, name: str, declared: str) -> str:
    """A GeoPackage column's field type from its declared type."""
    upper = declared.strip().upper()
    field = DECLARED.get("TEXT" if re.fullmatch(r"TEXT\(\d+\)", upper) else upper)
    if field is None:
        raise FeatureError(f"{path.name}: column {name!r} has type {declared!r}, not read")
    return field


def _inferred(path: Path, key: str, values: Sequence[Any]) -> str:
    """A GeoJSON key's field type from its non-null values."""
    present = [v for v in values if v is not None]
    if all(isinstance(v, int) for v in present):
        return INTEGER
    if all(isinstance(v, int | float) and not isinstance(v, bool) for v in present):
        return DOUBLE
    if all(isinstance(v, str) for v in present):
        return STRING
    kinds = sorted({type(v).__name__ for v in present})
    raise FeatureError(f"{path.name}: property {key!r} mixes {kinds}, not read")
