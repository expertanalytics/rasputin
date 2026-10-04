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
    crs = str(member["properties"]["name"])
    parse_crs(crs)
    return list(doc["features"]), crs


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
        x, y = f["geometry"]["coordinates"][:2]
        fields = ("station", "name", "series", "nve_area_km2", "watercourse", "river")
        stations.append(Station(x=x, y=y, **{k: props[k] for k in fields if k in props}))
    _unique(path, [s.station for s in stations], "station")
    return tuple(stations), crs


def read_references(path: Path) -> tuple[Mapping[str, Polygon | MultiPolygon], str]:
    """The reference polygon of each station in `path`, and the file's CRS."""
    features, crs = features_of(path)
    for f in features:
        _geometry_type(path, f, ("Polygon", "MultiPolygon"))
    numbers = [str(f["properties"]["station"]) for f in features]
    _unique(path, numbers, "station")
    return {n: shape(f["geometry"]) for n, f in zip(numbers, features, strict=True)}, crs


__all__ = ["Station", "features_of", "read_references", "read_stations"]
