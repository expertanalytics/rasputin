"""Gauging stations and their reference catchments, read back (increment 29, PR 3).

`docs/increments/29-nve-reference-catchments.md`, "The station set". Reads
what `rasputin fetch-stations` writes, or a user's own GeoJSON of points with
a `station` property. Every file must carry a `crs` member; duplicate station
numbers and the wrong geometry type are refused with a `ValueError` naming
them. The file is opened by `io/repository.py`, the one module in `io/` that
opens files.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field
from shapely.geometry import MultiPolygon, Polygon, shape

from tin_engine.crs import parse_crs

from .repository import read_json


class Station(BaseModel):
    """One gauging station, in its file's CRS. `series` holds the discharge
    series keys (`"1001.0"`), the join to NVE's discharge later."""

    model_config = ConfigDict(frozen=True)

    station: str = Field(pattern=r"^\d+\.\d+\.\d+$")
    name: str | None = None
    x: float
    y: float
    series: tuple[str, ...] = ()
    nve_area_km2: float | None = None
    watercourse: str | None = None
    river: str | None = None


def features_of(path: Path) -> tuple[list[dict[str, Any]], str]:
    """The features of the GeoJSON FeatureCollection in `path` and the text
    of its `crs` member, which is required and must name a CRS pyproj reads."""
    doc = read_json(path)
    member = doc.get("crs") if isinstance(doc, dict) else None
    if not member:
        raise ValueError(f"{path.name}: no crs member; the file must name its CRS")
    name = (member.get("properties") or {}).get("name") if isinstance(member, dict) else None
    if name is None:
        raise ValueError(f"{path.name}: the crs member has no name; it must name the CRS")
    parse_crs(str(name))
    if not isinstance(doc.get("features"), list):
        raise ValueError(f"{path.name}: no features list; the file is not a FeatureCollection")
    return list(doc["features"]), str(name)


def required(path: Path, feature: Mapping[str, Any], key: str) -> Any:
    """The property `key` of `feature`, refused naming it when missing."""
    props = feature.get("properties") or {}
    if key not in props:
        raise ValueError(f"{path.name}: a feature has no {key} property")
    return props[key]


def _unique(path: Path, keys: list[Any], what: str) -> None:
    seen: set[Any] = set()
    for key in keys:
        if key in seen:
            raise ValueError(f"{path.name}: {what} {key} appears twice")
        seen.add(key)


def _geometry_type(path: Path, feature: Mapping[str, Any], allowed: tuple[str, ...]) -> str:
    kind = str((feature.get("geometry") or {}).get("type"))
    if kind not in allowed:
        props = feature.get("properties") or {}
        label = props.get("station", props.get("objectid"))
        raise ValueError(f"{path.name}: {label} is a {kind}, not a {' or '.join(allowed)}")
    return kind


def read_stations(path: Path) -> tuple[tuple[Station, ...], str]:
    """The stations in `path`, in file order, and the file's CRS."""
    features, crs = features_of(path)
    stations = []
    for f in features:
        _geometry_type(path, f, ("Point",))
        props = f.get("properties") or {}
        try:
            x, y = f["geometry"]["coordinates"][:2]
        except (TypeError, ValueError) as exc:
            label = props.get("station")
            raise ValueError(f"{path.name}: station {label}: the point has no x and y") from exc
        fields = ("station", "name", "series", "nve_area_km2", "watercourse", "river")
        stations.append(Station(x=x, y=y, **{k: props[k] for k in fields if k in props}))
    _unique(path, [s.station for s in stations], "station")
    return tuple(stations), crs


def read_references(path: Path) -> tuple[Mapping[str, Polygon | MultiPolygon], str]:
    """The reference polygon of each station in `path`, and the file's CRS."""
    features, crs = features_of(path)
    for f in features:
        _geometry_type(path, f, ("Polygon", "MultiPolygon"))
    numbers = [str(required(path, f, "station")) for f in features]
    _unique(path, numbers, "station")
    polygons = {n: shape(f["geometry"]) for n, f in zip(numbers, features, strict=True)}
    for n, g in polygons.items():
        # The agreement divides by the polygon's area and perimeter.
        if g.is_empty or g.area == 0:
            raise ValueError(f"{path}: the reference polygon of station {n} has no area")
    return polygons, crs


__all__ = ["Station", "features_of", "read_references", "read_stations", "required"]
