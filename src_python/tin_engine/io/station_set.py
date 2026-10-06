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

from shapely.geometry import MultiPolygon, Polygon, shape

from tin_engine.hydrography import Lake, Station

from .geojson import read_collection
from .repository import read_json


def features_of(path: Path) -> tuple[list[dict[str, Any]], str]:
    """The features of the GeoJSON in `path` and the text of its `crs` member,
    which is required and must name a CRS pyproj reads (`read_collection`
    with no default); a refusal names the file, and a feature whose geometry
    is missing, not an object, or has no or an empty `type` is refused by its
    index."""
    try:
        features, crs = read_collection(read_json(path), default_crs=None)
    except ValueError as exc:
        raise ValueError(f"{path.name}: {exc}") from exc
    for i, f in enumerate(features):
        g = f.get("geometry")
        if not isinstance(g, dict) or not g.get("type"):
            raise ValueError(f"{path.name}: feature {i} has no geometry")
    return features, crs


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


def read_nve_lakes(path: Path) -> tuple[tuple[Lake, ...], str]:
    """The lakes in `path`, in file order, a MultiPolygon split into its parts,
    and the file's CRS. A geometry that is not a Polygon or MultiPolygon, or
    has no area, is refused naming the feature's `objectid` (else its index)."""
    features, crs = features_of(path)
    lakes: list[Lake] = []
    for k, f in enumerate(features):
        props = f.get("properties") or {}
        label = f"lake {props['objectid']}" if "objectid" in props else f"feature {k}"
        kind = (f.get("geometry") or {}).get("type")
        g = shape(f["geometry"]) if kind in ("Polygon", "MultiPolygon") else None
        if g is None or g.is_empty or g.area == 0:
            raise ValueError(f"{path.name}: {label} is not a polygon with an area ({kind})")
        number, name = props.get("vatnlnr") or None, props.get("navn")
        parts = g.geoms if isinstance(g, MultiPolygon) else [g]
        lakes += [Lake(None if number is None else int(number), name, p) for p in parts]
    return tuple(lakes), crs


__all__ = [
    "features_of",
    "read_nve_lakes",
    "read_references",
    "read_stations",
    "required",
]
