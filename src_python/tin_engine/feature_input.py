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

import math
import sqlite3
import time
from collections.abc import Callable, Iterable, Mapping
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import numpy.typing as npt
import shapely
from pydantic import BaseModel, ConfigDict
from pyproj import CRS
from shapely.geometry import LineString, MultiPolygon, Polygon, shape
from shapely.geometry.base import BaseGeometry

from tin_engine.crs import parse_crs, reprojector, same_crs, transform_definition
from tin_engine.domain import DomainPolygon
from tin_engine.features import DEFAULT_VOCABULARY, EdgeVocabulary, TerrainFeature
from tin_engine.io.geojson import GEOJSON_SUFFIXES, RFC7946_CRS, read_collection
from tin_engine.io.geopackage import layer_info, query_features
from tin_engine.io.gml import read_gml
from tin_engine.io.repository import open_geopackage, read_json

SUFFIXES = (*GEOJSON_SUFFIXES, ".gpkg", ".gml")
ACCEPTED = frozenset({"Polygon", "MultiPolygon", "LineString", "MultiLineString"})
#: R5: the region is the domain plus this many metres in the DEM's CRS, its
#: boundary densified to this spacing before it is moved.
MARGIN, DENSIFY = 100.0, 1_000.0
#: WGS 84's smallest and largest radii of curvature, metres (R5, "Long edges").
R_MIN, R_MAX = 6_335_439.0, 6_399_594.0
WATER_CODES = ("511", "512", "521", "522", "523")
#: A map whose values are integer class codes names its code system (16c, R5).
CORINE_CODES = "CORINE Land Cover level-3 code"
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
    unlisted value gets: names, ``"drop"`` (no constraint) or ``"refuse"``.
    ``codes``, when not empty, names the code system the values are integer
    codes of; such a map labels triangles (16c, R5)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    attribute: str
    classes: Mapping[str, tuple[str, ...]]
    otherwise: Literal["refuse", "drop"] | tuple[str, ...] = "refuse"
    notice: str = ""
    codes: str = ""


def _corine(name: str, attribute: str) -> ClassMap:
    water = dict.fromkeys(WATER_CODES, ("land_cover", "water"))
    return ClassMap(
        name=name,
        attribute=attribute,
        classes=water,
        otherwise=("land_cover",),
        notice=CORINE_NOTICE,
        codes=CORINE_CODES,
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
        codes=CORINE_CODES,
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
    #: 20c-3's land-cover stage, all off by default (the CLI's are 0.05, on, 0, 5):
    #: the repair's tolerance, the same-class merge, the simplification and D.
    repair_m: float = 0.0
    merge_same_class: bool = False
    tolerance_m: float = 0.0
    outline_snap_m: float = 0.0


class FeatureSet(BaseModel):
    """The features in source order; ``outside`` counts those clipped away,
    ``clipped`` those kept that crossed the domain's boundary, ``empty`` empty
    geometries skipped, ``scanned`` the GeoPackage layers read without an
    R-tree index (R3). ``crs``, ``layers`` and ``counts`` are per source
    (a layer for a GeoPackage only; ``counts`` the features kept from each,
    16e R6/D2); ``clip_seconds`` is the time spent after reading, for the
    ``features clip`` row, ``cleanup_seconds`` the land-cover stage's share of
    it (20c-3), ``cover_vertices`` the land cover's vertices after the clip to
    the read region and after the stage (None: no stage ran), and
    ``area_changed`` the area the outline rule gave another polygon, m²."""

    model_config = ConfigDict(frozen=True)

    features: tuple[TerrainFeature, ...]
    outside: int = 0
    clipped: int = 0
    empty: int = 0
    scanned: tuple[str, ...] = ()
    crs: tuple[str, ...] = ()
    layers: tuple[str | None, ...] = ()
    counts: tuple[int, ...] = ()
    clip_seconds: float = 0.0
    cleanup_seconds: float = 0.0
    cover_vertices: tuple[int, int] | None = None
    area_changed: float = 0.0


def source_region(domain: DomainPolygon, dem_crs: str | CRS, source_crs: str | CRS) -> Polygon:
    """R5: the domain's convex hull, in the DEM's CRS, buffered by 100 m (the
    hull of the domain buffered, to the arc's chord; 30d), densified to 1 km,
    moved into ``source_crs``; the convex hull of those points."""
    ring = shapely.segmentize(shapely.convex_hull(domain.polygon).buffer(MARGIN).exterior, DENSIFY)
    xy = shapely.get_coordinates(ring)
    if not same_crs(source_crs, dem_crs):
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
            inside = shapely.intersects_xy(region, xy[:, 0], xy[:, 1])
            keep = inside[:-1] | inside[1:]  # an edge with an end in the region meets it
            rest = np.flatnonzero(~keep)
            segments = shapely.linestrings(np.stack([xy[rest], xy[rest + 1]], axis=1))
            keep[rest] = shapely.intersects(region, segments)  # prepared region first
            if widening is not None:
                for i in np.flatnonzero(~keep[rest]):
                    k = rest[i]
                    w = widening((xy[k, 0], xy[k, 1]), (xy[k + 1, 0], xy[k + 1, 1]))
                    keep[k] = shapely.distance(segments[i], region) <= w
            out += _runs(xy, keep, closed)
    return tuple(out)


def _runs(
    xy: npt.NDArray[np.float64], keep: npt.NDArray[np.bool_], closed: bool
) -> list[LineString]:
    """Each maximal run of kept edges as a chain; a ring's run through its
    start vertex is one chain, joined across it."""
    if keep.all():
        return [LineString(xy)]
    if not keep.any():
        return []
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


def _linework(geometry: BaseGeometry) -> BaseGeometry:
    """Every edge ``pre_clip`` tests, as lines: a polygon's rings (each hole's
    too, wherever it lies), or the lines themselves."""
    return geometry.boundary if isinstance(geometry, Polygon | MultiPolygon) else geometry


def _untouched(line: BaseGeometry, domain: Polygon) -> bool:
    """True when ``intersection(line, domain)`` gives ``line`` back unchanged:
    two or more vertices, no repeated one, in the domain's interior, and simple
    (30b, section 3.4). ``domain`` first, so its prepared form is used."""
    xy = shapely.get_coordinates(line)
    return (
        len(xy) >= 2
        and not (xy[1:] == xy[:-1]).all(axis=1).any()
        and bool(shapely.contains_properly(domain, line))
        and bool(shapely.is_simple(line))
    )


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
    switches = (request.repair_m, request.merge_same_class, request.tolerance_m)
    tally.cleanup = request if any(switches) or request.outline_snap_m else None
    crss: list[str] = []
    layers: list[str | None] = []
    counts: list[int] = []
    for source in request.sources:
        before = len(tally.features)
        own, layer = tally.source(source)
        crss.append(own)
        layers.append(layer)
        counts.append(len(tally.features) - before)
    return FeatureSet(
        features=tuple(tally.features),
        outside=tally.outside,
        clipped=tally.clipped,
        empty=tally.empty,
        scanned=tuple(tally.scanned),
        crs=tuple(crss),
        layers=tuple(layers),
        counts=tuple(counts),
        clip_seconds=tally.seconds,
        cleanup_seconds=tally.cleanup_seconds,
        cover_vertices=tally.cover_vertices,
        area_changed=tally.area_changed,
    )


class _Tally:
    """What :func:`open_features` accumulates over its sources."""

    def __init__(self, domain: DomainPolygon, dem: CRS, vocabulary: EdgeVocabulary) -> None:
        self.domain, self.dem, self.vocabulary = domain, dem, vocabulary
        self.features: list[TerrainFeature] = []
        self.outside = self.clipped = self.empty = 0
        self.scanned: list[str] = []
        self.seconds = self.cleanup_seconds = self.area_changed = 0.0
        self.cleanup: FeatureRequest | None = None
        self.cover_vertices: tuple[int, int] | None = None
        shapely.prepare(domain.polygon)
        self.inner = domain.polygon.point_on_surface()

    def source(self, source: FeatureSource) -> tuple[str, str | None]:
        """Read one source and take in its features; its CRS text and layer."""
        read = read_source(
            source.path,
            source.layer,
            source.class_map.attribute,
            lambda own: source_region(self.domain, self.dem, own).bounds,
            source.crs,
        )
        if read.scanned:
            self.scanned.append(str(read.layer))
        try:
            self._take(source, read.crs, read.rows)
        except FeatureError:
            raise
        except (ValueError, KeyError) as exc:
            raise FeatureError(f"{source.path.name}: {exc}") from exc
        return read.crs, read.layer

    def _take(self, source: FeatureSource, own: str, rows: Iterable[tuple[Any, Any, Any]]) -> None:
        name, cmap = source.path.name, source.class_map
        src = parse_crs(own)
        move = None if same_crs(src, self.dem) else reprojector(src, self.dem)
        bound, region = None, source_region(self.domain, self.dem, self.dem)
        shapely.prepare(region)
        if move is not None and src.is_geographic:
            steps = transform_definition(src, self.dem)
            if "gridshift" not in steps and "deformation" not in steps:
                bound = _bound(self.dem, source_region(self.domain, self.dem, src))
        cover: list[tuple[Any, int, int, BaseGeometry]] = []
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
            code = _code(value, f"{name}: feature {fid}") if cmap.codes else None
            polygonal = code is not None and geometry.geom_type in ("Polygon", "MultiPolygon")
            moved = polygon = None
            if polygonal and code is not None and self.cleanup is not None:  # 20c-3: whole
                moved = shapely.transform(geometry, move) if move is not None else geometry
                if not _finite(moved):
                    raise FeatureError(
                        f"{name}: feature {fid} has a vertex with no image in the DEM"
                    )
                cover.append((fid, mask, code, moved))
                self.seconds += time.perf_counter() - t0
                continue
            if move is not None and bound is not None and bound[1](geometry):
                # Pre-clipped where its edges are straight, then moved (R5).
                chains = pre_clip(geometry, bound[2], bound[0])
                kept = [shapely.transform(c, move) for c in chains]
                polygon = shapely.transform(geometry, move) if polygonal else None
            else:  # moved first, then pre-clipped exactly in the DEM's CRS
                moved = shapely.transform(geometry, move) if move is not None else geometry
                kept = [moved]
                polygon = moved if polygonal else None
            if not all(_finite(g) for g in kept):
                raise FeatureError(f"{name}: feature {fid} has a vertex with no image in the DEM")
            if moved is not None:  # a feature whose linework misses the region keeps no edge
                kept = list(pre_clip(moved, region)) if region.intersects(_linework(moved)) else []
            self._add(fid, mask, kept, code, polygon)
            self.seconds += time.perf_counter() - t0
        if cover and self.cleanup is not None:
            self._clean(cover, region, self.cleanup)

    def _add(
        self, fid: Any, mask: int, kept: list[Any], code: int | None, polygon: BaseGeometry | None
    ) -> None:
        """One feature's linework clipped to the domain, kept, or counted outside."""
        coded: dict[str, Any] = {"code": code, "polygon": polygon}
        dom = self.domain.polygon
        lines = tuple(
            piece
            for line in kept
            for piece in (
                (line,)
                if _untouched(line, dom)
                else shapely.get_parts(shapely.intersection(line, dom))
            )
            if isinstance(piece, LineString) and piece.length > 0
        )
        if lines:
            self.features.append(TerrainFeature(fid=fid, mask=mask, lines=lines, **coded))
            # A dropped edge lies outside, so a pre-clipped chain is not covered.
            self.clipped += not all(self.domain.polygon.covers(g) for g in kept)
        elif polygon is not None and polygon.intersects(self.inner):
            # No boundary crosses the domain, and a point of it is inside: it
            # covers the domain, and labels it (R5). Kept, not clipped.
            self.features.append(TerrainFeature(fid=fid, mask=mask, lines=(), **coded))
        else:
            self.outside += 1

    def _clean(
        self, cover: list[tuple[Any, int, int, BaseGeometry]], region: Polygon, ask: FeatureRequest
    ) -> None:
        """20c-3: one source's land cover clipped to the read region, then
        repaired, merged by class, simplified and put to the outline, as asked;
        its lines and label polygons then go the usual way."""
        t0 = time.perf_counter()
        clipped = [_polygonal(shapely.intersection(g, region)) for *_, g in cover]
        self.outside += sum(g is None for g in clipped)
        if all(g is None for g in clipped):
            return
        items = [c[:3] for c, g in zip(cover, clipped, strict=True) if g is not None]
        polys = np.array([g for g in clipped if g is not None], dtype=object)
        before = int(shapely.get_num_coordinates(polys).sum())
        polys = shapely.coverage_clean(
            polys, snapping_distance=ask.repair_m, gap_width=ask.repair_m, merge_strategy="min_area"
        )
        if ask.merge_same_class:  # one polygon per class, under its first fid
            groups: dict[int, list[int]] = {}
            for i, (_, _, code) in enumerate(items):
                groups.setdefault(code, []).append(i)
            first = sorted(g[0] for g in groups.values())
            merged = [shapely.coverage_union_all(polys[groups[items[i][2]]]) for i in first]
            polys = np.array(merged, dtype=object)
            items = [items[i] for i in first]
        if ask.tolerance_m > 0:
            polys = shapely.coverage_simplify(polys, ask.tolerance_m, simplify_boundary=False)
        snapped = snap_to_outline(list(polys), self.domain.polygon, ask.outline_snap_m)
        self.area_changed += snapped.area_changed
        after = sum(int(shapely.get_num_coordinates(line)) for ls in snapped.lines for line in ls)
        old = self.cover_vertices or (0, 0)
        self.cover_vertices = (old[0] + before, old[1] + after)
        for (fid, mask, code), lines, polygon in zip(
            items, snapped.lines, snapped.polygons, strict=True
        ):
            if polygon.is_empty:  # the repair gave all of it to its neighbours
                self.empty += 1
            else:
                self._add(fid, mask, list(lines), code, polygon)
        took = time.perf_counter() - t0
        self.seconds += took
        self.cleanup_seconds += took


def _polygonal(geometry: BaseGeometry) -> BaseGeometry | None:
    """The polygonal parts of ``geometry`` as one (Multi)Polygon; None if none."""
    parts = [g for g in shapely.get_parts(geometry) if isinstance(g, Polygon) and not g.is_empty]
    return None if not parts else parts[0] if len(parts) == 1 else MultiPolygon(parts)


@dataclass(frozen=True, slots=True)
class SourceRows:
    """What :func:`read_source` read: the file's CRS text, the GeoPackage
    layer (None otherwise) and whether it was scanned without an R-tree, and
    ``(fid, geometry, attribute value)`` per feature, in file order."""

    crs: str
    layer: str | None
    scanned: bool
    rows: list[tuple[Any, Any, Any]]


def read_source(
    path: Path,
    layer: str | None,
    attribute: str | None,
    box_for: Callable[[str], tuple[float, float, float, float]],
    crs: str | None = None,
) -> SourceRows:
    """One source's features, unclipped, in its own CRS, by suffix (R1).

    ``attribute`` None reads a GeoPackage's primary key. ``box_for(own CRS)``
    is the box a GeoPackage's R-tree is queried with. ``crs``, when given,
    must be the file's; a GML file naming none takes it. Every refusal is a
    :class:`FeatureError` naming the file."""
    suffix = path.suffix.lower()
    if suffix not in SUFFIXES:
        raise FeatureError(f"{path.name}: unknown suffix {suffix or '(none)'}; use {SUFFIXES}")
    if layer is not None and suffix != ".gpkg":
        raise FeatureError(f"{path.name}: a layer applies to a GeoPackage only")
    try:
        if suffix == ".gpkg":
            with closing(open_geopackage(path)) as conn:
                info = layer_info(conn, layer)
                own = _own(path, crs, info.crs)
                scale = 1.0 if parse_crs(own).is_geographic else 100_000.0
                found = query_features(conn, info, box_for(own), attribute or info.pk, scale)
                rows = [(r.fid, r.geometry, r.value) for r in found]
            return SourceRows(own, info.table, info.rtree is None, rows)
        if suffix == ".gml":
            with path.open("rb") as stream:
                doc = read_gml(stream, attribute or "")
            if doc.crs is None and crs is None:
                raise FeatureError(f"{path.name}: no geometry names its srsName; give its CRS")
            own = _own(path, crs, doc.crs or str(crs))
            return SourceRows(
                own, None, False, [(f.fid, f.geometry, f.value) for f in doc.features]
            )
        features, text = read_collection(read_json(path), default_crs=RFC7946_CRS)
        try:  # a malformed feature's structure raises any of these
            rows = [
                (f.get("id", k), shape(g) if (g := f["geometry"]) else None, f["properties"] or {})
                for k, f in enumerate(features)
            ]
            values = [(k, g, p.get(attribute)) for k, g, p in rows]
        except (KeyError, TypeError, AttributeError, shapely.errors.ShapelyError) as exc:
            raise FeatureError(f"{path.name}: not a GeoJSON FeatureCollection ({exc})") from exc
        return SourceRows(_own(path, crs, text), None, False, values)
    except FeatureError:
        raise
    except (OSError, ValueError, KeyError, sqlite3.Error) as exc:
        raise FeatureError(f"{path.name}: {exc}") from exc


def _own(path: Path, given: str | None, text: str) -> str:
    """The file's CRS text, checked against the one given, if any."""
    if given is not None and not same_crs(given, text):
        raise FeatureError(f"{path.name} is in {text} but the given CRS is {given}")
    parse_crs(text)
    return text


def read_lake_polygons(
    path: Path, layer: str | None, point: tuple[float, float], point_crs: str
) -> tuple[tuple[BaseGeometry, ...], str]:
    """Increment 22: every polygon or multipolygon of a lake source near
    ``point`` (in ``point_crs``), in the source's own CRS, with its CRS text.
    A GeoPackage is queried with the point's box; lines and points are skipped."""

    def box_for(own: str) -> tuple[float, float, float, float]:
        ((x, y),) = reprojector(point_crs, own)([point])
        return (float(x), float(y), float(x), float(y))

    read = read_source(path, layer, None, box_for)
    lakes = tuple(
        shapely.force_2d(g)
        for _, g, _ in read.rows
        if g is not None and not g.is_empty and g.geom_type in ("Polygon", "MultiPolygon")
    )
    return lakes, read.crs


def _finite(geometry: BaseGeometry) -> bool:
    return bool(np.isfinite(shapely.get_coordinates(geometry)).all())


def _code(value: Any, what: str) -> int:
    """R5: a class code is an integer in 1 .. 2**31 - 1; anything else is refused."""
    try:
        code = int(str(value)) if not isinstance(value, list) else 0
    except ValueError:
        code = 0
    if not 1 <= code < 2**31:
        shown = value[0] if isinstance(value, list) and value else value
        raise FeatureError(f"{what}: {shown!r} is not a class code in 1 .. 2**31 - 1")
    return code


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


@dataclass(frozen=True, slots=True)
class OutlineSnap:
    """:func:`snap_to_outline`'s result, per input polygon in input order: its
    linework (not yet clipped to the domain) and its polygon after the rule;
    and ``area_changed``, the area inside the outline that changed polygon, m²."""

    lines: tuple[tuple[LineString, ...], ...]
    polygons: tuple[BaseGeometry, ...]
    area_changed: float


#: A cut point this close to the line through its neighbours is in line with
#: them, metres; a micrometre, far above a double's rounding at 1e7 m.
IN_LINE = 1e-6


def snap_to_outline(polygons: list[BaseGeometry], outline: Polygon, distance: float) -> OutlineSnap:
    """Ola's outline rule (20c-3, M6 (c)) at D = ``distance``: every part of a
    border within D of ``outline`` (every ring of it) goes onto it, and the
    stretches then on it are dropped from the linework. No buffer, no snap:
    the outline's segments are found in an ``STRtree``. D = 0 gives each
    polygon's rings unchanged; D at or above 100 m (the read region's margin)
    is refused."""
    if not (math.isfinite(distance) and 0 <= distance < MARGIN):
        raise ValueError(f"the outline snap must be >= 0 and under {MARGIN:g} m, got {distance}")
    parts = [[q for q in shapely.get_parts(p) if not q.is_empty] for p in polygons]
    if distance == 0:
        whole = [tuple(LineString(r) for q in qs for r in shapely.get_rings(q)) for qs in parts]
        return OutlineSnap(tuple(whole), tuple(polygons), 0.0)
    path = _Outline(outline, distance)
    out_lines: list[tuple[LineString, ...]] = []
    out_polys: list[BaseGeometry] = []
    loops: list[BaseGeometry] = []
    for polygon, qs in zip(polygons, parts, strict=True):
        lines: list[LineString] = []
        rebuilt: list[BaseGeometry] = []
        same: list[BaseGeometry] = []
        for q in qs:
            rings = [path.ring(shapely.get_coordinates(r)) for r in shapely.get_rings(q)]
            for ring, done in zip(shapely.get_rings(q), rings, strict=True):
                lines += [LineString(ring)] if done is None else _chains(done[0], done[1])
                loops += [] if done is None else _loops(shapely.get_coordinates(ring), *done)
            if all(r is None for r in rings):
                same.append(q)
                continue
            closed = [
                shapely.get_coordinates(r) if d is None else d[0]
                for r, d in zip(shapely.get_rings(q), rings, strict=True)
            ]
            if len(closed[0]) >= 3:
                shell = Polygon(closed[0], [h for h in closed[1:] if len(h) >= 3])
                rebuilt.append(
                    shell if shell.is_valid else shapely.make_valid(shell, method="structure")
                )
        new: BaseGeometry | None = polygon
        if len(same) < len(qs):  # ruling T2: only the rebuilt parts are unioned
            new = rebuilt[0] if len(rebuilt) == 1 else shapely.union_all(rebuilt)
            whole = [*same, *shapely.get_parts(shapely.get_parts(new))]  # make_valid nests
            new = _polygonal(shapely.geometrycollections(whole))
            if new is not None and not new.is_valid:  # argued, not proven: today's union
                new = _polygonal(shapely.union_all([*same, *rebuilt]))
        out_lines.append(tuple(lines))
        out_polys.append(Polygon() if new is None else new)
    area = shapely.intersection(shapely.union_all(loops), outline).area if loops else 0.0
    return OutlineSnap(tuple(out_lines), tuple(out_polys), float(area))


def _loops(
    old: npt.NDArray[np.float64], new: npt.NDArray[np.float64], on: list[bool], index: list[int]
) -> list[BaseGeometry]:
    """The regions one ring's rule changed (ruling T2): between two consecutive
    kept input vertices whose stretch changed, the loop of the old stretch and
    the new one; with no kept input vertex, the symmetric difference of the old
    and the new polygon (ruling P1)."""
    n, m = len(old) - 1, len(new)
    at = [(a, i) for a, i in enumerate(index) if i >= 0]
    if not at:
        was, now = (shapely.make_valid(Polygon(xy if len(xy) >= 3 else None)) for xy in (old, new))
        return [shapely.symmetric_difference(was, now)]
    out = []
    for (a, i), (b, j) in zip(at, at[1:] + at[:1], strict=True):
        span, step = ((j - i) % n, (b - a) % m) if len(at) > 1 else (n, m)
        if span > 1 or step > 1:
            xy = np.concatenate(
                [old[(i + np.arange(span + 1)) % n], new[(a + np.arange(step, -1, -1)) % m]]
            )
            out.append(shapely.make_valid(Polygon(xy)))
    return out


def _chains(xy: npt.NDArray[np.float64], on: list[bool]) -> list[LineString]:
    """A closed path cut where its segments lie on the outline (``on[k]``:
    the segment from point k to the next), those segments dropped."""
    n = len(xy)
    if not any(on):
        return [LineString([*xy, xy[0]])]
    start = on.index(True) + 1
    out, chain = [], [xy[start % n]]
    for k in range(start, start + n):
        if on[k % n]:
            out += [LineString(chain)] if len(chain) > 1 else []
            chain = [xy[(k + 1) % n]]
        else:
            chain.append(xy[(k + 1) % n])
    return out + ([LineString(chain)] if len(chain) > 1 else [])


#: A point the rule put on the outline: (outline ring, arc length along it, xy).
Placed = tuple[int, float, npt.NDArray[np.float64]]


class _Outline:
    """The outline's rings as arc-length paths, and an ``STRtree`` of their segments."""

    def __init__(self, outline: Polygon, d: float) -> None:
        self.d = d
        self.xy = [shapely.get_coordinates(r) for r in shapely.get_rings(outline)]
        self.cum: list[npt.NDArray[np.float64]] = []
        self.lap: list[npt.NDArray[np.float64]] = []  # two laps, for a way past the start
        self.owner: list[tuple[int, int]] = []
        for i, xy in enumerate(self.xy):
            lengths = np.hypot(*np.diff(xy, axis=0).T)
            cum = np.concatenate([[0.0], np.cumsum(lengths)])
            self.cum.append(cum)
            self.lap.append(np.concatenate([cum[:-1], cum[:-1] + cum[-1]]))
            self.owner += [(i, k) for k in range(len(lengths))]
        segments = np.concatenate([np.stack([xy[:-1], xy[1:]], axis=1) for xy in self.xy])
        self.tree = shapely.STRtree(shapely.linestrings(segments))

    def ring(
        self, xy: npt.NDArray[np.float64]
    ) -> tuple[npt.NDArray[np.float64], list[bool], list[int]] | None:
        """One closed ring after the rule, as points, on-outline flags and each
        point's input-vertex index (-1 for a cut, a moved point or an outline
        vertex); None if nothing in it moved."""
        edges = shapely.linestrings(np.stack([xy[:-1], xy[1:]], axis=1))
        near = np.zeros(len(edges), dtype=bool)
        near[self.tree.query(edges, predicate="dwithin", distance=self.d)[0]] = True
        if not near.any():
            return None
        raw: list[tuple[npt.NDArray[np.float64], int]] = []  # (point, input index; -1 a cut)
        for k in range(len(edges)):  # step 1: near edges cut, the same way from either side
            raw.append((xy[k], k))
            if near[k]:
                raw += [(c, -1) for c in self._cut(xy[k], xy[k + 1])]
        placed = self._place(np.array([r[0] for r in raw]))
        if all(p is None for p in placed):
            return None
        m = len(raw)
        moved = [p is not None for p in placed]
        # Step 3: a cut point that did not move is kept only next to one that did.
        seq = [
            (placed[i], raw[i][0], raw[i][1])
            for i in range(m)
            if moved[i] or raw[i][1] >= 0 or moved[i - 1] or moved[(i + 1) % m]
        ]
        seq = [e for i, e in enumerate(seq) if not _same(e[0], seq[i - 1][0])] or seq[:1]
        kept: list[Kept] = []
        for i, e in enumerate(seq):  # ... and dropped when in line with its neighbours
            a, b = (kept[-1] if kept else seq[i - 1]), seq[(i + 1) % len(seq)]
            if e[0] is None and e[2] < 0 and _off_line(e[1], a, b) <= IN_LINE:
                continue
            kept.append(e)
        own = [e[2] if e[0] is None else -1 for e in kept]
        if len(kept) < 2:
            return np.array([_xy(e) for e in kept]).reshape(-1, 2), [True] * len(kept), own
        points: list[npt.NDArray[np.float64]] = []
        on: list[bool] = []
        index: list[int] = []
        for i, e in enumerate(kept):
            nxt = kept[(i + 1) % len(kept)]
            points.append(_xy(e))
            index.append(own[i])
            if e[0] is None or nxt[0] is None:
                on.append(False)
                continue
            *steps, last = self._join(e[0], nxt[0])
            for point, flag in steps:
                on.append(flag)
                points.append(point)
                index.append(-1)
            on.append(last[1])
        return np.array(points), on, index

    def _cut(self, a: npt.NDArray[np.float64], b: npt.NDArray[np.float64]) -> Any:
        """The points cutting edge ab into pieces of at most D: the multiples
        of D along its line, counted from the foot of the CRS's origin and in
        the direction of its lexicographically last end. So both polygons
        sharing it get the same, and an edge the read region shortens keeps
        its cuts (to rounding). A multiple within D/2 of either end is dropped
        (ruling T1), so no piece is shorter than D/2."""
        flip = tuple(b) < tuple(a)
        lo, hi = (b, a) if flip else (a, b)
        u = (hi - lo) / math.hypot(*(hi - lo))
        t0, t1, h = float(np.dot(lo, u)), float(np.dot(hi, u)), self.d / 2
        k = np.arange(math.floor((t0 + h) / self.d) + 1, math.ceil((t1 - h) / self.d))
        cuts = lo + (k * self.d - t0)[:, None] * u
        return cuts[::-1] if flip else cuts

    def _place(self, xy: npt.NDArray[np.float64]) -> list[Placed | None]:
        """Step 2: each point within D to its nearest point on the outline, then
        to an outline vertex within D/2 of that, else to the nearest multiple of
        D along the outline (and to a vertex within D/2 of that)."""
        found = self.tree.query_nearest(shapely.points(xy), max_distance=self.d, all_matches=True)
        best: dict[int, int] = {}
        for i, k in zip(*found, strict=True):  # ties: the first segment, fixed
            best[int(i)] = min(int(k), best.get(int(i), int(k)))
        out: list[Placed | None] = [None] * len(xy)
        for i, k in best.items():
            r, j = self.owner[k]
            a, b = self.xy[r][j], self.xy[r][j + 1]
            ab = b - a
            t = min(max(float(np.dot(xy[i] - a, ab) / max(np.dot(ab, ab), 1e-300)), 0.0), 1.0)
            s = self.cum[r][j] + t * (self.cum[r][j + 1] - self.cum[r][j])
            m = round(s / self.d) * self.d  # within D/2 of the end only if a vertex is
            out[i] = self._vertex(r, s) or self._vertex(r, m) or (r, m, self._at(r, m))
        return out

    def _vertex(self, r: int, s: float) -> Placed | None:
        cum = self.cum[r]
        j = int(np.searchsorted(cum, s))
        j = j if j < len(cum) and (j == 0 or cum[j] - s <= s - cum[j - 1]) else j - 1
        if abs(cum[j] - s) > self.d / 2:
            return None
        j %= len(cum) - 1  # the last is the first
        return (r, float(cum[j]), self.xy[r][j])

    def _at(self, r: int, s: float) -> npt.NDArray[np.float64]:
        cum, xy = self.cum[r], self.xy[r]
        j = min(int(np.searchsorted(cum, s, "right")) - 1, len(cum) - 2)
        t = (s - cum[j]) / max(cum[j + 1] - cum[j], 1e-300)
        return xy[j] + t * (xy[j + 1] - xy[j])  # type: ignore[no-any-return]

    def _join(self, p: Placed, q: Placed) -> list[tuple[npt.NDArray[np.float64], bool]]:
        """The way from p to q, each point with whether the segment into it lies
        on the outline: along the outline's shorter way through its vertices,
        or straight where that way is an inlet (longer than twice |pq| plus 2D),
        less what of the straight line lies on the outline at either end."""
        gap = float(np.hypot(*(q[2] - p[2])))
        if p[0] != q[0]:
            return [(q[2], False)]
        cum, xy, n = self.cum[p[0]], self.xy[p[0]], len(self.cum[p[0]]) - 1
        length = cum[-1]
        ahead = (q[1] - p[1]) % length
        lap = self.lap[p[0]]
        if ahead <= length - ahead:
            lo, hi = np.searchsorted(lap, [p[1], p[1] + ahead], "right")
            js = [j % n for j in range(lo, hi) if lap[j] < p[1] + ahead]
        else:
            lo, hi = np.searchsorted(lap, [q[1], q[1] + length - ahead], "right")
            js = [j % n for j in range(lo, hi) if lap[j] < q[1] + length - ahead][::-1]
        if min(ahead, length - ahead) <= 2 * gap + 2 * self.d:
            return [*((xy[j], True) for j in js), (q[2], True)]
        on = [_off_line(xy[j], (None, p[2], -1), (None, q[2], -1)) <= IN_LINE for j in js]
        head = next((i for i, f in enumerate(on) if not f), len(js))
        tail = next((i for i, f in enumerate(on[::-1]) if not f), len(js))
        rest = js[len(js) - tail :] if tail else []
        return [
            *((xy[j], True) for j in js[:head]),
            *((xy[j], i > 0) for i, j in enumerate(rest)),
            (q[2], bool(rest)),
        ]


Kept = tuple[Placed | None, npt.NDArray[np.float64], int]  # the int: input index, -1 a cut


def _xy(e: Kept) -> npt.NDArray[np.float64]:
    return e[1] if e[0] is None else e[0][2]


def _same(p: Placed | None, q: Placed | None) -> bool:
    return p is not None and q is not None and p[0] == q[0] and p[1] == q[1]


def _off_line(x: npt.NDArray[np.float64], a: Kept, b: Kept) -> float:
    """Distance from x to segment ab, its ends in a fixed order (either direction
    gives the same answer)."""
    u, v = sorted((_xy(a), _xy(b)), key=tuple)
    uv = v - u
    t = min(max(float(np.dot(x - u, uv) / max(np.dot(uv, uv), 1e-300)), 0.0), 1.0)
    return float(np.hypot(*(x - u - t * uv)))
