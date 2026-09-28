"""Terrain features: polygons and lines read, mapped to bits, and clipped (16b-1).

`docs/increments/16b-terrain-polygons.md` R1, R4-R7. A source is a GeoJSON
``FeatureCollection``, a GeoPackage layer or OGR's GML2, by suffix. Each
feature's attribute goes through a :class:`ClassMap` to vocabulary names and so
to a mask; its rings and lines are pre-clipped to a region around the domain
(R5), moved into the DEM's CRS vertex by vertex, and clipped to the domain as
linework, never as areas (R6). What comes out is frozen data in the DEM's CRS.

**The pre-clip keeps whole edges** (R5, Ola 2026-09-28): an edge is kept when
it lies within its widening ``w(e)`` of the region, and dropped whole
otherwise, so every kept edge reaches the engine between its own two vertices.
``w(e)`` bounds how far the edge's straight line in the DEM's CRS strays from
its straight line in the source's, so a dropped edge misses the domain (I4).
The bound is claimed only for a geographic source over a Transverse Mercator
DEM; any other pair is moved first and pre-clipped in the DEM's CRS with
``w = 0``, which is exact.

Blocking (sqlite3, GEOS): an async caller runs :func:`open_features` in
``asyncio.to_thread``. No ``_core``.
"""

from __future__ import annotations

import json
import math
import sqlite3
import time
from collections.abc import Callable, Iterable, Mapping
from contextlib import closing
from pathlib import Path
from typing import Any, Literal

import numpy as np
import numpy.typing as npt
import shapely
from pydantic import BaseModel, ConfigDict
from pyproj import CRS
from shapely.geometry import LineString, Polygon, shape
from shapely.geometry.base import BaseGeometry

from tin_engine.crs import parse_crs, reprojector, transform_definition
from tin_engine.domain import GEOJSON_DEFAULT_CRS, DomainPolygon
from tin_engine.features import DEFAULT_VOCABULARY, EdgeVocabulary
from tin_engine.io.geopackage import layer_info, query_features
from tin_engine.io.gml import read_gml
from tin_engine.io.repository import open_geopackage

GEOJSON_SUFFIXES = (".geojson", ".json")
SUFFIXES = (*GEOJSON_SUFFIXES, ".gpkg", ".gml")
ACCEPTED = frozenset({"Polygon", "MultiPolygon", "LineString", "MultiLineString"})
#: R5: the region is the domain plus this many metres in the DEM's CRS, its
#: boundary densified to this spacing before it is moved.
MARGIN, DENSIFY = 100.0, 1_000.0
#: WGS 84's smallest and largest radii of curvature, metres (R5, "Long edges").
R_MIN, R_MAX = 6_335_439.0, 6_399_594.0
WATER_CODES = ("511", "512", "521", "522", "523")
CORINE_NOTICE = (
    "Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, "
    "Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and "
    "re-encoded, produced with funding by the European Union, not endorsed by the EU"
)

Widening = Callable[[tuple[float, float], tuple[float, float]], float]


class FeatureError(ValueError):
    """A source or feature refused, naming the file and, for a feature, its fid."""


class ClassMap(BaseModel):
    """An attribute's values to vocabulary names (R4). ``otherwise`` is what an
    unlisted value gets: names, ``"drop"`` (no constraint) or ``"refuse"``."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    attribute: str
    classes: Mapping[str, tuple[str, ...]]
    otherwise: Literal["refuse", "drop"] | tuple[str, ...] = "refuse"
    notice: str = ""


def _corine(name: str, attribute: str) -> ClassMap:
    water = dict.fromkeys(WATER_CODES, ("land_cover", "water"))
    return ClassMap(
        name=name,
        attribute=attribute,
        classes=water,
        otherwise=("land_cover",),
        notice=CORINE_NOTICE,
    )


CLASS_MAPS: dict[str, ClassMap] = {
    "property": ClassMap(
        name="property",
        attribute="property",
        classes={p.name: (p.name,) for p in DEFAULT_VOCABULARY.properties},
    ),
    "corine": _corine("corine", "Code_18"),
    "corine-water": ClassMap(
        name="corine-water",
        attribute="Code_18",
        classes=dict.fromkeys(WATER_CODES, ("water",)),
        otherwise="drop",
        notice=CORINE_NOTICE,
    ),
    "clc18_kode": _corine("clc18_kode", "clc18_kode"),
}


class FeatureSource(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    path: Path
    class_map: ClassMap
    layer: str | None = None
    crs: str | None = None


class FeatureRequest(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    sources: tuple[FeatureSource, ...]
    vocabulary: EdgeVocabulary = DEFAULT_VOCABULARY


class TerrainFeature(BaseModel):
    """One feature: its fid, mask, and clipped lines in the DEM's CRS; a
    closed line (first point repeated) is an unclipped ring."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    fid: int | str
    mask: int
    lines: tuple[LineString, ...]


class FeatureSet(BaseModel):
    """The features in source order; ``outside`` counts those clipped away,
    ``empty`` empty geometries skipped. ``crs`` and ``layers`` are per source
    (a layer for a GeoPackage only); ``clip_seconds`` is the time spent after
    reading, for the ``features clip`` row."""

    model_config = ConfigDict(frozen=True)

    features: tuple[TerrainFeature, ...]
    outside: int = 0
    empty: int = 0
    crs: tuple[str, ...] = ()
    layers: tuple[str | None, ...] = ()
    clip_seconds: float = 0.0


def source_region(domain: DomainPolygon, dem_crs: str | CRS, source_crs: str | CRS) -> Polygon:
    """R5: the domain, in the DEM's CRS, buffered by 100 m, densified to 1 km,
    moved into ``source_crs``; the convex hull of those points."""
    ring = shapely.segmentize(domain.polygon.buffer(MARGIN).exterior, DENSIFY)
    xy = shapely.get_coordinates(ring)
    if parse_crs(source_crs) != parse_crs(dem_crs):
        xy = reprojector(dem_crs, source_crs)(xy)
    hull = shapely.convex_hull(shapely.multipoints(xy))
    assert isinstance(hull, Polygon)
    return hull


def pre_clip(
    geometry: BaseGeometry, region: Polygon, widening: Widening | None = None
) -> tuple[LineString, ...]:
    """R5: whole edges within ``widening(a, b)`` (0 if None) of ``region``, as
    chains of the original vertices: a polygon's exterior's, then each hole's.
    A ring with no edge dropped stays closed with its start vertex."""
    shapely.prepare(region)
    out: list[LineString] = []
    for part in shapely.get_parts(geometry):
        if isinstance(part, Polygon):
            rings = [(r, True) for r in (part.exterior, *part.interiors)]
        else:
            rings = [(part, False)]
        for ring, closed in rings:
            xy = shapely.get_coordinates(ring)
            if len(xy) < 2:
                continue
            segments = shapely.linestrings(np.stack([xy[:-1], xy[1:]], axis=1))
            keep = np.asarray(shapely.intersects(segments, region))
            if widening is not None:
                for k in np.flatnonzero(~keep):
                    w = widening((xy[k, 0], xy[k, 1]), (xy[k + 1, 0], xy[k + 1, 1]))
                    keep[k] = shapely.distance(segments[k], region) <= w
            out += _runs(xy, keep, closed)
    return tuple(out)


def _runs(
    xy: npt.NDArray[np.float64], keep: npt.NDArray[np.bool_], closed: bool
) -> list[LineString]:
    """Each maximal run of kept edges as a chain; a ring's run through its
    start vertex is one chain, joined across it."""
    if keep.all():
        return [LineString(xy)]
    n = len(keep)
    first = int(np.flatnonzero(~keep)[0]) + 1 if closed else 0
    chains: list[LineString] = []
    run: list[int] = []
    for k in [(first + i) % n for i in range(n)] + [-1]:
        if k >= 0 and keep[k]:
            run.append(k)
        elif run:
            chains.append(LineString([xy[run[0]], *(xy[j + 1] for j in run)]))
            run = []
    return chains


def open_features(request: FeatureRequest, domain: DomainPolygon, dem_crs: str | CRS) -> FeatureSet:
    """Every source's features, in the DEM's CRS, clipped to ``domain`` (in the
    DEM's CRS already). Raises :class:`FeatureError` on any refusal."""
    vocabulary = request.vocabulary
    for source in request.sources:
        cmap = source.class_map
        named = [n for names in cmap.classes.values() for n in names]
        try:
            vocabulary.mask(*named, *(() if isinstance(cmap.otherwise, str) else cmap.otherwise))
        except ValueError as exc:
            raise FeatureError(f"class map {cmap.name}: {exc}") from exc
    tally = _Tally(domain, parse_crs(dem_crs), vocabulary)
    crss: list[str] = []
    layers: list[str | None] = []
    for source in request.sources:
        own, layer = tally.source(source)
        crss.append(own)
        layers.append(layer)
    return FeatureSet(
        features=tuple(tally.features),
        outside=tally.outside,
        empty=tally.empty,
        crs=tuple(crss),
        layers=tuple(layers),
        clip_seconds=tally.seconds,
    )


class _Tally:
    """What :func:`open_features` accumulates over its sources."""

    def __init__(self, domain: DomainPolygon, dem: CRS, vocabulary: EdgeVocabulary) -> None:
        self.domain, self.dem, self.vocabulary = domain, dem, vocabulary
        self.features: list[TerrainFeature] = []
        self.outside = self.empty = 0
        self.seconds = 0.0
        shapely.prepare(domain.polygon)

    def source(self, source: FeatureSource) -> tuple[str, str | None]:
        """Read one source and take in its features; its CRS text and layer."""
        path, attribute = source.path, source.class_map.attribute
        suffix = path.suffix.lower()
        if suffix not in SUFFIXES:
            raise FeatureError(f"{path.name}: unknown suffix {suffix or '(none)'}; use {SUFFIXES}")
        if source.layer is not None and suffix != ".gpkg":
            raise FeatureError(f"{path.name}: a layer applies to a GeoPackage only")
        try:
            if suffix == ".gpkg":
                with closing(open_geopackage(path)) as conn:
                    layer = layer_info(conn, source.layer)
                    own = self._own(source, layer.crs)
                    box = source_region(self.domain, self.dem, own).bounds
                    scale = 1.0 if parse_crs(own).is_geographic else 100_000.0
                    found = query_features(conn, layer, box, attribute, scale)
                    self._take(source, own, ((r.fid, r.geometry, r.value) for r in found))
                return own, layer.table
            if suffix == ".gml":
                with path.open("rb") as stream:
                    doc = read_gml(stream, attribute)
                if doc.crs is None and source.crs is None:
                    raise FeatureError(f"{path.name}: no geometry names its srsName; give its CRS")
                own = self._own(source, doc.crs or str(source.crs))
                self._take(source, own, ((f.fid, f.geometry, f.value) for f in doc.features))
                return own, None
            doc = json.loads(path.read_text())
            member = doc.get("crs")
            own = self._own(source, member["properties"]["name"] if member else GEOJSON_DEFAULT_CRS)
            rows = [
                (f.get("id", k), f["geometry"] and shape(f["geometry"]), f["properties"] or {})
                for k, f in enumerate(doc["features"])
            ]
            self._take(source, own, ((k, g, p.get(attribute)) for k, g, p in rows))
            return own, None
        except FeatureError:
            raise
        except (OSError, ValueError, KeyError, TypeError, AttributeError, sqlite3.Error) as exc:
            raise FeatureError(f"{path.name}: {exc}") from exc

    def _own(self, source: FeatureSource, text: str) -> str:
        """The file's CRS text, checked against the one given, if any."""
        if source.crs is not None and parse_crs(source.crs) != parse_crs(text):
            raise FeatureError(f"{source.path.name} is in {text} but the given CRS is {source.crs}")
        parse_crs(text)
        return text

    def _take(self, source: FeatureSource, own: str, rows: Iterable[tuple[Any, Any, Any]]) -> None:
        name, cmap = source.path.name, source.class_map
        src = parse_crs(own)
        move = None if src == self.dem else reprojector(src, self.dem)
        bound, region = None, source_region(self.domain, self.dem, self.dem)
        if move is not None and src.is_geographic:
            steps = transform_definition(src, self.dem)
            if "gridshift" not in steps and "deformation" not in steps:
                bound = _bound(self.dem, source_region(self.domain, self.dem, src))
        for fid, geometry, value in rows:
            t0 = time.perf_counter()
            if geometry is None or geometry.is_empty:
                self.empty += 1
                continue
            geometry = shapely.force_2d(geometry)
            if geometry.geom_type not in ACCEPTED:
                raise FeatureError(
                    f"{name}: feature {fid} is a {geometry.geom_type}; input is polygons and lines"
                )
            mask = _mask(cmap, self.vocabulary, value, f"{name}: feature {fid}")
            if mask is None:
                continue
            moved = None
            if move is not None and bound is not None and bound[1](geometry):
                # Pre-clipped where its edges are straight, then moved (R5).
                chains = pre_clip(geometry, bound[2], bound[0])
                kept = [shapely.transform(c, move) for c in chains]
            else:  # moved first, then pre-clipped exactly in the DEM's CRS
                moved = shapely.transform(geometry, move) if move is not None else geometry
                kept = [moved]
            if not all(_finite(g) for g in kept):
                raise FeatureError(f"{name}: feature {fid} has a vertex with no image in the DEM")
            if moved is not None:
                kept = list(pre_clip(moved, region))
            lines = tuple(
                piece
                for line in kept
                for piece in shapely.get_parts(shapely.intersection(line, self.domain.polygon))
                if isinstance(piece, LineString) and piece.length > 0
            )
            if lines:
                self.features.append(TerrainFeature(fid=fid, mask=mask, lines=lines))
            else:
                self.outside += 1
            self.seconds += time.perf_counter() - t0


def _finite(geometry: BaseGeometry) -> bool:
    return bool(np.isfinite(shapely.get_coordinates(geometry)).all())


def _mask(cmap: ClassMap, vocabulary: EdgeVocabulary, value: Any, what: str) -> int | None:
    """R4: the mask ``value`` maps to, None for ``drop``; a list is a union."""
    names: list[str] = []
    for item in value if isinstance(value, list) else [value]:
        key = None if item is None else str(item)
        if key is not None and key in cmap.classes:
            names += cmap.classes[key]
        elif cmap.otherwise == "drop":
            return None
        elif cmap.otherwise == "refuse":
            raise FeatureError(f"{what}: {cmap.attribute} {item!r} is not in the map {cmap.name}")
        else:
            names += cmap.otherwise
    return vocabulary.mask(*names)


def _bound(
    dem: CRS, region: Polygon
) -> tuple[Widening, Callable[[BaseGeometry], bool], Polygon] | None:
    """R5, "Long edges": for a Transverse Mercator DEM, ``w(e)`` in degrees, the
    test of where it applies, and the region in the source CRS; else None."""
    op = dem.coordinate_operation
    if op is None or op.method_name != "Transverse Mercator":
        return None
    params = {p.name: float(p.value) for p in op.params}
    lam0, k0 = params["Longitude of natural origin"], params["Scale factor at natural origin"]
    top = max(abs(region.bounds[1]), abs(region.bounds[3]))

    def widening(a: tuple[float, float], b: tuple[float, float]) -> float:
        (la, pa), (lb, pb) = a, b
        star, low = max(abs(pa), abs(pb)), 0.0 if pa * pb < 0 else min(abs(pa), abs(pb))
        far = math.radians(max(abs(la - lam0), abs(lb - lam0)))
        dphi, dlam = math.radians(pb - pa), math.radians(lb - la)
        ell = k0 * R_MAX * math.hypot(dphi, math.cos(math.radians(low)) * dlam) / math.cos(far)
        bend = (
            1.25 * ell**2 * (1.09 * math.tan(math.radians(star)) + math.tan(far)) / (8 * k0 * R_MIN)
        )
        return math.degrees(bend / (k0 * R_MIN * math.cos(math.radians(max(star, top) + 0.1))))

    def applies(geometry: BaseGeometry) -> bool:
        xy = shapely.get_coordinates(geometry)
        return bool((np.abs(xy[:, 1]) <= 80).all() and (np.abs(xy[:, 0] - lam0) <= 60).all())

    return widening, applies, region
