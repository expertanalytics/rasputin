"""Terrain features: class polygons and lines read, and the edges chosen from them.

`docs/increments/111-one-land-cover-map.md`, *3. Design*.

A source is a GeoJSON ``FeatureCollection``, a GeoPackage layer or OGR's GML2,
by suffix. :func:`read_class_polygons` reads a land-cover file's polygons
with their class codes (a :class:`LandCoverClasses` says what each value
means), :func:`read_linestrings` a file of lines; each moves every vertex
into the DEM's CRS in one call. :func:`edges_of` turns the chosen classes'
borders, cleaned as one coverage, and the chosen lines into the edges the
mesh follows, clipped to the domain. What comes out is frozen data in the
DEM's CRS.

Blocking (sqlite3, GEOS): an async caller runs the readers in
``asyncio.to_thread``. ``_core`` only through :func:`simplify_borders`.
"""

from __future__ import annotations

import math
import sqlite3
import time
from collections.abc import Callable, Mapping, Sequence
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Protocol

import numpy as np
import numpy.typing as npt
import shapely
from pydantic import BaseModel, ConfigDict
from pyproj import CRS
from shapely.geometry import LineString, MultiPolygon, Polygon, shape
from shapely.geometry.base import BaseGeometry

from tin_engine.border_simplify import simplify_borders
from tin_engine.crs import parse_crs, reprojector, same_crs
from tin_engine.domain import DomainPolygon
from tin_engine.features import EdgeVocabulary, TerrainFeature
from tin_engine.io.geojson import GEOJSON_SUFFIXES, RFC7946_CRS, read_collection
from tin_engine.io.geopackage import layer_info, query_features
from tin_engine.io.gml import read_gml
from tin_engine.io.repository import open_geopackage, read_json
from tin_engine.legends import Legend

if TYPE_CHECKING:
    from tin_engine.preferred import Replaced

SUFFIXES = (*GEOJSON_SUFFIXES, ".gpkg", ".gml")
POLYGONAL, LINEAR = ("Polygon", "MultiPolygon"), ("LineString", "MultiLineString", "LinearRing")
#: The read region is the domain plus this many metres in the DEM's CRS, or two
#: DEM cells when that is more, its boundary densified to this
#: spacing before it is moved.
MARGIN, DENSIFY = 100.0, 1_000.0
#: A land-cover class's share in a hydrological cell; a class without one is unspecified.
Category = Literal["forest", "lake", "glacier", "savanna", "unspecified"]
#: Parts of the domain's uncovered rest narrower than this, metres, stay out of
#: the border simplification: a chosen width, assuming a metric CRS.
REST_MIN_WIDTH_M = 0.01
#: A land-cover polygon on its way to the clean-up: fid, mask, code, merge key, geometry.
Cover = tuple[Any, int, int | None, int | str, BaseGeometry]


class FeatureError(ValueError):
    """A source or feature refused, naming the file and, for a feature, its fid."""


class Split(BaseModel):
    """Class ``code`` split by its second ``attribute``: `split_key` of that
    value to a sub-class code; an unlisted or missing value keeps ``code``."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    code: int
    attribute: str
    codes: Mapping[str, int]


def split_key(value: Any) -> str | None:
    """The text a second value is matched by: a whole number's integer text
    (52, 52.0, " 52 " and "52.0" are "52"), other text as it is, and None for
    anything else (None, a bool, 52.5, NaN, a list)."""
    if isinstance(value, bool) or not isinstance(value, int | float | str):
        return None
    if isinstance(value, int):
        return str(value)
    try:
        number = float(value)
    except ValueError:
        return str(value)
    return str(int(number)) if math.isfinite(number) and number.is_integer() else None


class LandCoverClasses(BaseModel):
    """What a class polygon's attribute value means, from a land-cover class
    file (`tin_engine.land_cover_classes`): ``codes`` names the code system
    the values are integer codes of, and ``legend`` colours and names those
    codes. ``classes`` holds the listed values (each class code as text, and
    each text of ``string_codes``); ``otherwise`` is what an unlisted value
    gets: ``"keep"`` (the value is the code), ``"drop"`` (no polygon) or
    ``"refuse"``. ``categories`` holds each class code that names its cell
    category, ``water_codes`` the water classes, and ``prefer`` the one class
    a preferred source's polygons take, per category. ``accept_no_class``
    says class code 0 (no polygon) is accepted; ``accuracy_m`` None: the file
    states no positional accuracy."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    attribute: str
    codes: str
    legend: Legend
    classes: frozenset[str]
    otherwise: Literal["refuse", "drop", "keep"] = "refuse"
    notice: str = ""
    licence: str = ""
    accuracy_m: float | None = None
    string_codes: Mapping[str, int] = {}
    splits: tuple[Split, ...] = ()
    categories: Mapping[int, Category] = {}
    water_codes: frozenset[int] = frozenset()
    prefer: Mapping[Category, int] = {}
    accept_no_class: bool = False


class EdgeChoice(BaseModel):
    """Which borders and lines become edges: every class of the map
    (``every_class``), the classes ``codes``, and the lines of the
    ``--lines`` keys ``lines``."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    every_class: bool = False
    codes: frozenset[int] = frozenset()
    lines: frozenset[str] = frozenset()


class EdgeCleanup(BaseModel):
    """The land-cover clean-up: the repair's tolerance, the same-class
    merge, the simplification and the outline snap distance D (0 is off).
    No field has a default: `tin_engine.land_cover_classes.edge_cleanup`
    derives them from the sources, as the CLI does."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    repair_m: float
    merge_same_class: bool
    max_shift_m: float
    outline_snap_m: float


class FeatureSet(BaseModel):
    """The edges in order; ``outside`` counts those clipped away,
    ``clipped`` those kept that crossed the domain's boundary, ``empty`` the
    areas the clean-up emptied. ``counts`` holds the features kept from each
    source of edges; ``cleanup_seconds`` is the land-cover clean-up's time,
    ``cover_vertices`` the land cover's vertices after the clip to the read
    region and after the clean-up (None: no clean-up ran), and
    ``area_changed`` the area the outline rule gave another polygon, m²."""

    model_config = ConfigDict(frozen=True)

    features: tuple[TerrainFeature, ...]
    outside: int = 0
    clipped: int = 0
    empty: int = 0
    counts: tuple[int, ...] = ()
    cleanup_seconds: float = 0.0
    cover_vertices: tuple[int, int] | None = None
    area_changed: float = 0.0


class _Map(Protocol):
    """What :func:`edges_of` reads of a land-cover map
    (`tin_engine.land_cover_map.LandCoverMap`)."""

    @property
    def areas(self) -> Sequence[tuple[int, BaseGeometry]]: ...
    @property
    def outside(self) -> Sequence[tuple[int, BaseGeometry]]: ...
    @property
    def classes(self) -> LandCoverClasses: ...
    @property
    def replaced(self) -> Replaced | None: ...


def source_region(
    domain: DomainPolygon, dem_crs: str | CRS, source_crs: str | CRS, *, cell: float = 0.0
) -> Polygon:
    """The read region: the domain's convex hull, in the DEM's CRS, buffered by
    100 m or two DEM cells of side `cell`, whichever is more (the hull of the
    domain buffered, to the arc's chord), densified to 1 km, moved into
    ``source_crs``; the convex hull of those points."""
    margin = max(MARGIN, 2 * cell)
    ring = shapely.segmentize(shapely.convex_hull(domain.polygon).buffer(margin).exterior, DENSIFY)
    xy = shapely.get_coordinates(ring)
    if not same_crs(source_crs, dem_crs):
        xy = reprojector(dem_crs, source_crs)(xy)
    hull = shapely.convex_hull(shapely.multipoints(xy))
    assert isinstance(hull, Polygon)
    return hull


def pre_clip(geometry: BaseGeometry, region: Polygon) -> tuple[LineString, ...]:
    """Whole edges of ``geometry``'s rings that meet ``region``, as chains of
    the original vertices: a polygon's exterior's, then each hole's. A ring
    with no edge dropped stays closed with its start vertex."""
    out: list[LineString] = []
    for ring in shapely.get_rings(shapely.get_parts(geometry)):
        xy = shapely.get_coordinates(ring)
        if len(xy) < 2:
            continue
        inside = shapely.intersects_xy(region, xy[:, 0], xy[:, 1])
        keep = inside[:-1] | inside[1:]  # an edge with an end in the region meets it
        rest = np.flatnonzero(~keep)
        segments = shapely.linestrings(np.stack([xy[rest], xy[rest + 1]], axis=1))
        keep[rest] = shapely.intersects(region, segments)  # prepared region first
        out += _runs(xy, keep)
    return tuple(out)


def _runs(xy: npt.NDArray[np.float64], keep: npt.NDArray[np.bool_]) -> list[LineString]:
    """Each maximal run of a ring's kept edges as a chain, the run through
    its start vertex joined across it."""
    if keep.all():
        return [LineString(xy)]
    if not keep.any():
        return []
    n = len(keep)
    first = int(np.flatnonzero(~keep)[0]) + 1
    chains: list[LineString] = []
    run: list[int] = []
    for k in [(first + i) % n for i in range(n)] + [-1]:
        if k >= 0 and keep[k]:
            run.append(k)
        elif run:
            chains.append(LineString([xy[run[0]], *(xy[j + 1] for j in run)]))
            run = []
    return chains


def _untouched(line: BaseGeometry, domain: Polygon) -> bool:
    """True when ``intersection(line, domain)`` gives ``line`` back unchanged:
    two or more vertices, no repeated one, simple, its interior in the
    domain's interior and its ends in the closure (the DE-9IM pattern is
    ``covers`` plus "the outline meets no interior point of the line", so the
    overlay adds no node).
    ``domain`` first, so its prepared form is used."""
    xy = shapely.get_coordinates(line)
    return (
        len(xy) >= 2
        and not (xy[1:] == xy[:-1]).all(axis=1).any()
        and bool(shapely.relate_pattern(domain, line, "T**F**FF*"))
        and bool(shapely.is_simple(line))
    )


def edges_of(
    land_cover: _Map | None,
    lines: Sequence[tuple[str, LineString]],
    choice: EdgeChoice,
    cleanup: EdgeCleanup,
    domain: DomainPolygon,
    vocabulary: EdgeVocabulary,
    *,
    region: Polygon,
) -> FeatureSet:
    """The edges the mesh follows, in the DEM's CRS, clipped to ``domain``:
    the borders of the map's chosen classes (``land_cover``, and ``water`` for
    a water class), cleaned as one coverage within ``region`` (the read
    region) when a clean-up switch is on or the map was cut by ``--prefer``,
    then the chosen classes' polygons the read box misses (the map's
    ``outside``), by their rings only: without clean-up each keeps the
    edges of its rings that meet ``region`` and labels nothing; with it
    each is counted outside. Then each chosen ``--lines`` key's lines, as
    given. Pure. Raises
    :class:`FeatureError` when the border simplification refuses."""
    tally = _Tally(domain)
    counts: list[int] = []
    if land_cover is not None and (choice.every_class or choice.codes):
        land = vocabulary.mask("land_cover")
        water = land | vocabulary.mask("water")
        wet = land_cover.classes.water_codes
        n = len(land_cover.areas)

        def chosen(pairs: Sequence[tuple[int, BaseGeometry]], start: int) -> list[Cover]:
            return [
                (i, water if code in wet else land, code, code, g)
                for i, (code, g) in enumerate(pairs, start=start)
                if choice.every_class or code in choice.codes
            ]

        cover, far = chosen(land_cover.areas, 0), chosen(land_cover.outside, n)
        cut = land_cover.replaced is not None
        switches = (cleanup.repair_m, cleanup.merge_same_class, cleanup.max_shift_m)
        if any(switches) or cleanup.outline_snap_m or cut:
            snap = cleanup.repair_m or (-1.0 if cut else 0.0)
            tally.clean(cover, region, cleanup, snap)
            tally.outside += len(far)
        else:
            shapely.prepare(region)
            for fid, mask, code, _, g in cover + far:
                kept = list(pre_clip(g, region)) if region.intersects(g.boundary) else []
                tally.add(fid, mask, kept, code, g if fid < n else None)
        counts.append(len(tally.features))
    for key in dict.fromkeys(k for k, _ in lines if k in choice.lines):
        before = len(tally.features)
        mask = vocabulary.mask(key)
        for i, (k, line) in enumerate(lines):
            if k == key:
                tally.add(f"{key} {i}", mask, [line], None, None)
        counts.append(len(tally.features) - before)
    return FeatureSet(
        features=tuple(tally.features),
        outside=tally.outside,
        clipped=tally.clipped,
        empty=tally.empty,
        counts=tuple(counts),
        cleanup_seconds=tally.cleanup_seconds,
        cover_vertices=tally.cover_vertices,
        area_changed=tally.area_changed,
    )


class _Tally:
    """What :func:`edges_of` accumulates."""

    def __init__(self, domain: DomainPolygon) -> None:
        self.domain = domain
        self.features: list[TerrainFeature] = []
        self.outside = self.clipped = self.empty = 0
        self.cleanup_seconds = self.area_changed = 0.0
        self.cover_vertices: tuple[int, int] | None = None
        shapely.prepare(domain.polygon)
        self.inner = domain.polygon.point_on_surface()

    def add(
        self, fid: Any, mask: int, kept: list[Any], code: int | None, polygon: BaseGeometry | None
    ) -> None:
        """One feature's linework clipped to the domain, kept, or counted outside."""
        coded: dict[str, Any] = {"code": code, "polygon": polygon}
        dom = self.domain.polygon
        whole = [_untouched(line, dom) for line in kept]
        lines = tuple(
            piece
            for line, untouched in zip(kept, whole, strict=True)
            for piece in (
                (line,) if untouched else shapely.get_parts(shapely.intersection(line, dom))
            )
            if isinstance(piece, LineString) and piece.length > 0
        )
        if lines:
            self.features.append(TerrainFeature(fid=fid, mask=mask, lines=lines, **coded))
            # A dropped edge lies outside, so a pre-clipped chain is not covered;
            # an untouched chain lies in the closure, so it is.
            self.clipped += not all(u or dom.covers(g) for u, g in zip(whole, kept, strict=True))
        elif polygon is not None and polygon.intersects(self.inner):
            # No boundary crosses the domain, and a point of it is inside: it
            # covers the domain, and labels it. Kept, not clipped.
            self.features.append(TerrainFeature(fid=fid, mask=mask, lines=(), **coded))
        else:
            self.outside += 1

    def clean(self, cover: list[Cover], region: Polygon, ask: EdgeCleanup, snap: float) -> None:
        """The chosen land cover clipped to the read region, then repaired,
        merged by class and put to the outline, as asked; with the simplification
        on then clipped to the domain, each border moved at most `max_shift_m`
        at the repair distance's clearance, and its lines cut
        where they lie on the outline. The domain's uncovered rest, less its
        parts narrower than `REST_MIN_WIDTH_M`, is simplified with them, so no
        border is fixed for lack of a neighbour. Its lines and label polygons
        (coded polygons only) then go the usual way."""
        t0 = time.perf_counter()
        clipped = [_polygonal(shapely.intersection(g, region)) for *_, g in cover]
        self.outside += sum(g is None for g in clipped)
        if all(g is None for g in clipped):
            return
        items = [c[:4] for c, g in zip(cover, clipped, strict=True) if g is not None]
        polys = np.array([g for g in clipped if g is not None], dtype=object)
        before = int(shapely.get_num_coordinates(polys).sum())
        polys = shapely.coverage_clean(
            polys, snapping_distance=snap, gap_width=ask.repair_m, merge_strategy="min_area"
        )
        if ask.merge_same_class:  # one polygon per class, under its first fid
            groups: dict[int | str, list[int]] = {}
            for i, (*_, key) in enumerate(items):
                groups.setdefault(key, []).append(i)
            first = sorted(g[0] for g in groups.values())
            merged = [shapely.coverage_union_all(polys[groups[items[i][3]]]) for i in first]
            polys = np.array(merged, dtype=object)
            items = [items[i] for i in first]
        snapped = snap_to_outline(list(polys), self.domain.polygon, ask.outline_snap_m)
        self.area_changed += snapped.area_changed
        lines_of, polygons = snapped.lines, snapped.polygons
        if ask.max_shift_m > 0:  # simplification on: area kept inside the domain
            dom = self.domain.polygon
            inside = clip_to_domain(polygons, dom)
            rest = shapely.difference(dom, shapely.coverage_union_all(inside))
            half = REST_MIN_WIDTH_M / 2
            wide = [p for p in shapely.get_parts(rest) if not p.buffer(-half).is_empty]
            try:
                done = simplify_borders([*inside, *wide], ask.max_shift_m, ask.repair_m)
            except ValueError as exc:
                raise FeatureError(
                    f"border simplification refused: {exc}; --edges-max-shift 0 "
                    "meshes this source with its borders as read"
                ) from exc
            polygons = done.polygons[: len(inside)]
            tree = _outline_tree(dom)  # once per clean-up, not once per polygon
            lines_of = tuple(_outline_lines(g, tree) for g in polygons)
        after = sum(int(shapely.get_num_coordinates(line)) for ls in lines_of for line in ls)
        old = self.cover_vertices or (0, 0)
        self.cover_vertices = (old[0] + before, old[1] + after)
        for (fid, mask, code, _), lines, polygon, ruled in zip(
            items, lines_of, polygons, snapped.polygons, strict=True
        ):
            if ruled.is_empty:  # the repair gave all of it to its neighbours
                self.empty += 1
            else:
                self.add(fid, mask, list(lines), code, None if code is None else polygon)
        self.cleanup_seconds += time.perf_counter() - t0


def clip_to_domain(polygons: Sequence[BaseGeometry], domain: BaseGeometry) -> list[BaseGeometry]:
    """Each of ``polygons`` cut to ``domain``, an empty Polygon where none is
    left, in input order.

    Cut one at a time, neighbours that tile need not tile after: where one
    keeps a vertex lying a hair off the outline and the other gets its own
    crossing beside it (1e-11 to 3.5e-10 m apart), their shared border no
    longer meets vertex for vertex. The cut set is cleaned as one coverage at
    GEOS's own snapping distance, which moves next to no area."""
    inside = [
        Polygon() if (c := _polygonal(shapely.intersection(g, domain))) is None else c
        for g in polygons
    ]
    return list(shapely.coverage_clean(np.array(inside, dtype=object), snapping_distance=-1))


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
    extra: tuple[tuple[Any, ...], ...] = ()


def read_source(
    path: Path,
    layer: str | None,
    attribute: str | None,
    box_for: Callable[[str], tuple[float, float, float, float]],
    crs: str | None = None,
    extra: tuple[str, ...] = (),
) -> SourceRows:
    """One source's features, unclipped, in its own CRS, by suffix.

    ``attribute`` None reads a GeoPackage's primary key. ``box_for(own CRS)``
    is the box a GeoPackage's R-tree is queried with. ``crs``, when given,
    must be the file's; a GML file naming none takes it. ``extra`` names more
    attributes, read per row into ``SourceRows.extra`` (None where a feature
    or layer has none; ``()`` when none is asked). Every refusal is a
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
                found = list(
                    query_features(conn, info, box_for(own), attribute or info.pk, scale, extra)
                )
                rows = [(r.fid, r.geometry, r.value) for r in found]
            more = tuple(r.extra for r in found) if extra else ()
            return SourceRows(own, info.table, info.rtree is None, rows, more)
        if suffix == ".gml":
            with path.open("rb") as stream:
                doc = read_gml(stream, attribute or "", extra)
            if doc.crs is None and crs is None:
                raise FeatureError(f"{path.name}: no geometry names its srsName; give its CRS")
            own = _own(path, crs, doc.crs or str(crs))
            read = [(f.fid, f.geometry, f.value) for f in doc.features]
            more = tuple(f.extra for f in doc.features) if extra else ()
            return SourceRows(own, None, False, read, more)
        features, text = read_collection(read_json(path), default_crs=RFC7946_CRS)
        try:  # a malformed feature's structure raises any of these
            rows = [
                (f.get("id", k), shape(g) if (g := f["geometry"]) else None, f["properties"] or {})
                for k, f in enumerate(features)
            ]
            values = [(k, g, p.get(attribute)) for k, g, p in rows]
            more = tuple(tuple(p.get(n) for n in extra) for *_, p in rows) if extra else ()
        except (KeyError, TypeError, AttributeError, shapely.errors.ShapelyError) as exc:
            raise FeatureError(f"{path.name}: not a GeoJSON FeatureCollection ({exc})") from exc
        return SourceRows(_own(path, crs, text), None, False, values, more)
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
    """Every polygon or multipolygon of a lake source near
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


def _repaired(
    geometry: BaseGeometry, xy: npt.NDArray[np.float64], *, what: str = "a polygon"
) -> BaseGeometry | None:
    """``geometry`` if valid, else ``shapely.make_valid(geometry)`` within the
    box of the points ``xy`` widened by 100 m, as polygons, except where two
    rings share a stretch of edge: there both sides keep the same parity.
    None if empty.
    A ring that crosses itself makes GEOS's overlays throw. Only the part near
    the mesh is ever overlaid, and repairing a large polygon whole takes seconds.
    make_valid keeps a point when an odd number of rings, counted over all
    parts, wind around it. So the rings cut to the box, with the box's own
    edge, split it into faces, and a face is kept when an odd number of parts
    contain a point of it. Cutting the polygon instead of its rings can swap
    inside and outside when a ring crosses itself outside the box.
    GEOS raises on some spikes too thin to rebuild: a :class:`FeatureError`
    names ``what`` and its first vertex."""
    if geometry.is_valid:
        return geometry
    box = shapely.box(*(xy.min(axis=0) - MARGIN), *(xy.max(axis=0) + MARGIN))
    try:
        lines = shapely.clip_by_rect(shapely.boundary(geometry), *box.bounds)
        noded = shapely.union_all([lines, box.exterior])
        faces = shapely.get_parts(shapely.polygonize(shapely.get_parts(noded)))
        x, y = shapely.get_coordinates(shapely.point_on_surface(faces)).T
        odd = np.zeros(len(faces), dtype=bool)
        for part in shapely.get_parts(geometry):
            odd ^= shapely.contains_xy(part, x, y)
        kept = shapely.get_parts(shapely.coverage_union_all(faces[odd]))
    except shapely.errors.GEOSException as exc:
        x0, y0 = shapely.get_coordinates(geometry)[0]
        raise FeatureError(
            f"{what} near ({x0:.1f}, {y0:.1f}) in the DEM's CRS cannot be repaired: it has"
            " a spike too thin to rebuild; clean the file"
        ) from exc
    return MultiPolygon(list(kept)) if len(kept) else None


@dataclass(frozen=True, slots=True)
class ClassPolygons:
    """What :func:`read_class_polygons` read: ``(code, polygon)`` per row in
    file order, in the DEM's CRS and unchecked; the file's CRS text, its
    GeoPackage layer, whether it was scanned without an R-tree, and how many
    rows were skipped as empty."""

    areas: tuple[tuple[int, BaseGeometry], ...]
    crs: str
    layer: str | None
    scanned: bool
    empty: int


@dataclass(frozen=True, slots=True)
class SourceLines:
    """What :func:`read_linestrings` read: the lines in file order, in the
    DEM's CRS, a MultiLineString as its parts; the rest as :class:`ClassPolygons`."""

    lines: tuple[LineString, ...]
    crs: str
    layer: str | None
    scanned: bool
    empty: int


def _moved(
    name: str,
    found: list[tuple[Any, BaseGeometry]],
    own: str,
    mesh_crs: str,
    threads: int,
) -> list[BaseGeometry]:
    """Every geometry of ``found`` (fid, geometry) moved from ``own`` into
    ``mesh_crs`` in one transform call; one with a vertex that has no image
    is refused naming its fid."""
    geometries = [g for _, g in found]
    if found and not same_crs(own, mesh_crs):
        move = reprojector(own, mesh_crs, threads=threads)
        geometries = list(shapely.transform(np.array(geometries, dtype=object), move))
    for (fid, _), g in zip(found, geometries, strict=True):
        if not _finite(g):
            raise FeatureError(f"{name}: feature {fid} has a vertex with no image in the DEM")
    return geometries


def read_class_polygons(
    path: Path,
    layer: str | None,
    classes: LandCoverClasses,
    mesh_crs: str,
    box_for: Callable[[str], tuple[float, float, float, float]],
    *,
    crs: str | None = None,
    threads: int = 1,
) -> ClassPolygons:
    """Every polygon of a land-cover file with its class code, moved to
    ``mesh_crs`` vertex by vertex, unchecked; a GeoPackage is queried with
    ``box_for(own CRS)``, and ``crs``, when given, must be the file's. Empty
    rows are skipped and counted. A :class:`FeatureError` names file and
    feature: a value no class takes under ``otherwise``, a point, a
    collection or a line, and a vertex with no image in the DEM."""
    spelt: dict[str, str] = {}  # lower case to the first spelling, read as written
    for split in classes.splits:
        spelt.setdefault(split.attribute.lower(), split.attribute)
    read = read_source(path, layer, classes.attribute, box_for, crs, tuple(spelt.values()))
    name, empty = path.name, 0
    found: list[tuple[Any, BaseGeometry]] = []
    codes: list[int] = []
    for (fid, g, value), more in zip(read.rows, read.extra or [()] * len(read.rows), strict=True):
        what = f"{name}: feature {fid}"
        if g is None or g.is_empty:
            empty += 1
            continue
        if g.geom_type not in POLYGONAL:
            raise FeatureError(f"{what} is a {g.geom_type}; class polygons only, lines to --lines")
        if not _listed(classes, value, what):
            continue
        codes.append(_code(classes, value, what, dict(zip(spelt, more, strict=True))))
        found.append((fid, shapely.force_2d(g)))
    moved = _moved(name, found, read.crs, mesh_crs, threads)
    areas = tuple(zip(codes, moved, strict=True))
    return ClassPolygons(areas, read.crs, read.layer, read.scanned, empty)


def read_linestrings(
    path: Path,
    layer: str | None,
    mesh_crs: str,
    box_for: Callable[[str], tuple[float, float, float, float]],
    *,
    threads: int = 1,
) -> SourceLines:
    """Every line of a file, moved to ``mesh_crs`` vertex by vertex, a
    MultiLineString as its parts; read and refused as by
    :func:`read_class_polygons`, a polygon too."""
    read = read_source(path, layer, None, box_for)
    name, empty = path.name, 0
    found: list[tuple[Any, BaseGeometry]] = []
    for fid, g, _ in read.rows:
        if g is None or g.is_empty:
            empty += 1
            continue
        if g.geom_type not in LINEAR:
            where = "polygons go to --prefer or --land-cover" if g.geom_type in POLYGONAL else (
                "lines only"
            )  # fmt: skip
            raise FeatureError(f"{name}: feature {fid} is a {g.geom_type}; {where}")
        flat = shapely.force_2d(g)
        found += [(fid, LineString(p)) for p in shapely.get_parts(flat) if not p.is_empty]
    lines = tuple(
        g for g in _moved(name, found, read.crs, mesh_crs, threads) if isinstance(g, LineString)
    )
    return SourceLines(lines, read.crs, read.layer, read.scanned, empty)


def _near(
    xy: npt.NDArray[np.float64], mesh_crs: str
) -> Callable[[str], tuple[float, float, float, float]]:
    """The box of the points ``xy`` (in ``mesh_crs``) widened by 100 m, in a
    source's own CRS: what a GeoPackage is queried with."""

    def box_for(own: str) -> tuple[float, float, float, float]:
        ring = shapely.box(*xy.min(axis=0), *xy.max(axis=0)).buffer(MARGIN).exterior
        near = shapely.get_coordinates(shapely.segmentize(ring, DENSIFY))
        moved = near if same_crs(own, mesh_crs) else reprojector(mesh_crs, own)(near)
        x0, y0, x1, y1 = shapely.bounds(shapely.multipoints(moved))
        return (float(x0), float(y0), float(x1), float(y1))

    return box_for


def read_polygons(
    path: Path,
    layer: str | None,
    mesh_crs: str,
    xy: npt.NDArray[np.float64],
    tell: Callable[[int, int], None] | None = None,
    named: tuple[str, Callable[[list[Any]], list[int]]] | None = None,
    *,
    lines: str = "",
) -> tuple[list[BaseGeometry], int]:
    """Every polygon or multipolygon of a source, moved to ``mesh_crs`` vertex
    by vertex and repaired near the points ``xy`` (:func:`_repaired`); a
    GeoPackage is queried with their box widened by 100 m. The count returned
    is how many were repaired, and ``tell`` is called with it and the number
    read. With ``named`` (an attribute and a filter of its values to the
    indices kept), only those. A line is skipped, or with ``lines`` refused
    with those words, naming file and feature."""
    read = read_source(path, layer, None if named is None else named[0], _near(xy, mesh_crs))
    rows = read.rows
    if named is not None:
        rows = [rows[k] for k in named[1]([v for *_, v in rows])]
    for fid, g, _ in rows if lines else ():
        if g is not None and g.geom_type in LINEAR:
            raise FeatureError(f"{path.name}: feature {fid} is a {g.geom_type}; {lines}")
    move = None if same_crs(read.crs, mesh_crs) else reprojector(read.crs, mesh_crs)
    flat = [
        (fid, shapely.force_2d(g))
        for fid, g, _ in rows
        if g is not None and not g.is_empty and g.geom_type in POLYGONAL
    ]
    moved = [g if move is None else shapely.transform(g, move) for _, g in flat]
    fixed = [
        _repaired(g, xy, what=f"{path.name}: feature {fid}")
        for (fid, _), g in zip(flat, moved, strict=True)
    ]
    repaired = sum(a is not b for a, b in zip(moved, fixed, strict=True))
    if tell is not None:
        tell(repaired, len(moved))
    return [g for g in fixed if g is not None], repaired


def _finite(geometry: BaseGeometry) -> bool:
    return bool(np.isfinite(shapely.get_coordinates(geometry)).all())


def _code(rules: LandCoverClasses, value: Any, what: str, second: Mapping[str, Any]) -> int:
    """A class code is an integer in 1 .. 2**31 - 1, or a text the class file's
    ``string_codes`` turn into one; anything else is refused. A split class
    takes the sub-class code its second value (``second``, by lower-case
    attribute name) keys to; a first value equal to a sub-class code is refused."""
    try:
        text = "0" if isinstance(value, list) else str(value)
        code = rules.string_codes.get(text) or int(text)
    except ValueError:
        code = 0
    if not 1 <= code < 2**31:
        shown = value[0] if isinstance(value, list) and value else value
        raise FeatureError(f"{what}: {shown!r} is not a class code in 1 .. 2**31 - 1")
    for split in rules.splits:
        if code == split.code:
            key = split_key(second.get(split.attribute.lower()))
            return split.codes.get(key, code) if key is not None else code
        if code in split.codes.values():
            raise FeatureError(f"{what}: {rules.attribute} {code} is a sub-class code")
    return code


def _listed(classes: LandCoverClasses, value: Any, what: str) -> bool:
    """Whether ``value`` (a list: each item) is a class of ``classes``,
    False for one ``drop`` drops; ``refuse`` raises."""
    for item in value if isinstance(value, list) else [value]:
        if item is not None and str(item) in classes.classes:
            continue
        if classes.otherwise == "drop":
            return False
        if classes.otherwise == "refuse":
            raise FeatureError(
                f"{what}: {classes.attribute} {item!r} is not in the land-cover class file "
                f"{classes.name}"
            )
    return True


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
    """The outline rule at D = ``distance``: every part of a
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
        if len(same) < len(qs):  # only the rebuilt parts are unioned
            new = rebuilt[0] if len(rebuilt) == 1 else shapely.union_all(rebuilt)
            whole = [*same, *shapely.get_parts(shapely.get_parts(new))]  # make_valid nests
            new = _polygonal(shapely.geometrycollections(whole))
            if new is not None and not new.is_valid:  # argued, not proven: fall back
                new = _polygonal(shapely.union_all([*same, *rebuilt]))
        out_lines.append(tuple(lines))
        out_polys.append(Polygon() if new is None else new)
    area = shapely.intersection(shapely.union_all(loops), outline).area if loops else 0.0
    return OutlineSnap(tuple(out_lines), tuple(out_polys), float(area))


def _loops(
    old: npt.NDArray[np.float64], new: npt.NDArray[np.float64], on: list[bool], index: list[int]
) -> list[BaseGeometry]:
    """The regions one ring's rule changed: between two consecutive
    kept input vertices whose stretch changed, the loop of the old stretch and
    the new one; with no kept input vertex, the symmetric difference of the old
    and the new polygon."""
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


def _outline_tree(outline: Polygon) -> tuple[npt.NDArray[Any], shapely.STRtree]:
    """``outline``'s segments as lines, and an ``STRtree`` over them."""
    rims = [shapely.get_coordinates(r) for r in shapely.get_rings(outline)]
    segments = shapely.linestrings(
        np.concatenate([np.stack([r[:-1], r[1:]], axis=1) for r in rims])
    )
    return segments, shapely.STRtree(segments)


def _outline_lines(
    polygon: BaseGeometry, outline: tuple[npt.NDArray[Any], shapely.STRtree]
) -> tuple[LineString, ...]:
    """Every ring of ``polygon`` cut where an edge lies on
    the outline (one outline segment within ``IN_LINE`` of both its ends),
    those edges dropped, as :func:`snap_to_outline`'s lines. ``outline`` is
    :func:`_outline_tree`'s pair, built once per clean-up."""
    segments, tree = outline
    out: list[LineString] = []
    for ring in shapely.get_rings(shapely.get_parts(polygon)):
        xy = shapely.get_coordinates(ring)
        edges = shapely.linestrings(np.stack([xy[:-1], xy[1:]], axis=1))
        ei, si = tree.query(edges, predicate="dwithin", distance=IN_LINE)
        near = [
            shapely.distance(shapely.points(xy[ei + k]), segments[si]) <= IN_LINE for k in (0, 1)
        ]
        on = np.zeros(len(xy) - 1, dtype=bool)
        on[ei[near[0] & near[1]]] = True
        out += _chains(xy[:-1], on.tolist())
    return tuple(out)


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
        so no piece is shorter than D/2."""
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
