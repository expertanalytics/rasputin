"""A sub-catchment network's file: its links, and its units' polygons on demand.

A GeoPackage layer is read with one ``SELECT`` of the named columns, its
geometries only for the units asked for (by ID, or through its R-tree at a
point); a GeoJSON file is parsed whole on each read. Every refusal is a
:class:`FeatureError` naming the file.
"""

from __future__ import annotations

import sqlite3
from collections.abc import Collection, Sequence
from contextlib import closing
from pathlib import Path
from typing import Any

import shapely
from shapely.geometry import Point, shape
from shapely.geometry.base import BaseGeometry

from tin_engine.crs import reprojector
from tin_engine.feature_input import FeatureError, read_source
from tin_engine.io.geojson import GEOJSON_SUFFIXES, RFC7946_CRS, read_collection
from tin_engine.io.geopackage import _q, decode_geometry, layer_info
from tin_engine.io.repository import open_geopackage, read_json
from tin_engine.network import Link, unit_id

#: IDs per ``IN`` list: under SQLite's oldest limit on bound parameters.
CHUNK = 900


def _table(
    path: Path,
    layer: str | None,
    attributes: Sequence[str],
    ids: Collection[str] | None = None,
) -> tuple[str, list[tuple[Any, ...]]]:
    """The file's CRS and, per feature, the values of `attributes`, then its
    geometry when `ids` is given, only for the features whose first
    attribute is one of `ids`."""
    suffix = path.suffix.lower()
    if suffix not in (*GEOJSON_SUFFIXES, ".gpkg"):
        raise FeatureError(f"{path.name}: use .gpkg, .geojson or .json")
    if layer is not None and suffix != ".gpkg":
        raise FeatureError(f"{path.name}: a layer applies to a GeoPackage only")
    try:
        if suffix != ".gpkg":
            features, crs = read_collection(read_json(path), default_crs=RFC7946_CRS)
            props = [f.get("properties") or {} for f in features]
            for name in attributes:
                if not any(name in p for p in props):
                    raise FeatureError(f"{path.name}: no feature has the attribute {name!r}")
            if ids is None:
                return crs, [tuple(p.get(a) for a in attributes) for p in props]
            wanted = [k for k, p in enumerate(props) if unit_id(p.get(attributes[0])) in ids]
            return crs, [
                (*(props[k].get(a) for a in attributes), shape(features[k]["geometry"]))
                for k in wanted
            ]
        with closing(open_geopackage(path)) as conn:
            info = layer_info(conn, layer)
            names = {
                r[1].lower(): r[1] for r in conn.execute(f"PRAGMA table_info({_q(info.table)})")
            }
            for name in attributes:
                if name.lower() not in names:
                    raise FeatureError(f"{path.name}: layer {info.table!r} has no {name!r}")
            columns = [_q(names[a.lower()]) for a in attributes]
            select = f"SELECT {', '.join(columns)} FROM {_q(info.table)}"
            if ids is None:
                return info.crs, conn.execute(select).fetchall()
            select = f"SELECT {', '.join(columns)}, {_q(info.column)} FROM {_q(info.table)}"
            rows, keys = [], list(ids)
            for k in range(0, len(keys), CHUNK):
                chunk = keys[k : k + CHUNK]
                where = f" WHERE {columns[0]} IN ({', '.join('?' * len(chunk))})"
                rows += conn.execute(select + where, chunk).fetchall()
            return info.crs, [(*r[:-1], r[-1] and decode_geometry(r[-1])[1]) for r in rows]
    except FeatureError:
        raise
    except (OSError, ValueError, KeyError, TypeError, sqlite3.Error) as exc:
        raise FeatureError(f"{path.name}: {exc}") from exc


def read_links(
    path: Path,
    layer: str | None,
    id_attr: str,
    down_attr: str,
    area_attr: str | None,
    area_unit: str | None,
    name_attr: str | None,
) -> tuple[str, list[Link]]:
    """The file's CRS and every unit's link; the stated area in m2."""
    named = [id_attr, down_attr, *(a for a in (area_attr, name_attr) if a is not None)]
    crs, rows = _table(path, layer, named)
    scale = 1e6 if area_unit == "km2" else 1.0
    links = []
    for row in rows:
        values = dict(zip(named, row, strict=True))
        area = values.get(area_attr) if area_attr else None
        name = values.get(name_attr) if name_attr else None
        links.append(
            Link(
                id=unit_id(values[id_attr]),
                down=unit_id(values[down_attr]),
                upstream_area_m2=None if area is None else float(area) * scale,
                name=None if name is None else str(name),
            )
        )
    return crs, links


def units_at(
    path: Path, layer: str | None, id_attr: str, point: Point
) -> list[tuple[str, BaseGeometry]]:
    """The units whose box holds `point` (in the file's CRS): candidates for
    the unit under it, through a GeoPackage's R-tree; every unit of a GeoJSON file."""
    box = (point.x, point.y, point.x, point.y)
    found = read_source(path, layer, id_attr, lambda _own: box)
    return [(unit_id(v), g) for _, g, v in found.rows if g is not None]


def unit_geometries(
    path: Path, layer: str | None, id_attr: str, ids: Collection[str]
) -> dict[str, BaseGeometry]:
    """The polygons of the units `ids`, by ID."""
    _crs, rows = _table(path, layer, [id_attr], set(ids))
    return {unit_id(own): g for own, g in rows if g is not None}


def outlet_point(path: Path, layer: str | None, id_attr: str, unit: str, crs: str) -> Point | None:
    """The outlet point of `unit` moved into `crs`; None when the file has
    none, refused when it has two."""
    own, rows = _table(path, layer, [id_attr], {unit})
    if len(rows) > 1:
        raise FeatureError(f"{path.name}: {len(rows)} outlet points for the unit {unit}")
    if not rows or rows[0][1] is None:
        return None
    ((x, y),) = reprojector(own, crs)([shapely.get_coordinates(rows[0][1])[0]])
    return Point(float(x), float(y))
