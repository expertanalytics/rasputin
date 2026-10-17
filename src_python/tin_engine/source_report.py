"""What each catalogue source holds, read from its headers: `rasputin sources`.

A local source is read through the readers a run uses (`repository_for`), a
remote one through its tile cache (`CacheRepository`) when cached, a vector
file through `io/geopackage.py` and `io/geojson.py`; nothing is decoded but
headers, and a GeoJSON file, which has none, is parsed whole. `fits` says
which options each source can be given to. The figures follow STAC's conventions: coverage as
a longitude/latitude box west, south, east, north; cell size in metres; the
CRS as ``EPSG:n`` where there is one. Kept apart from `sources.py`, which
imports no pyproj or shapely.
"""

from __future__ import annotations

import sqlite3
from collections.abc import Callable, Mapping
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import shapely
from pydantic import BaseModel, ConfigDict

from tin_engine.class_maps import load, shipped_names
from tin_engine.crs import crs_label, parse_crs, reprojector
from tin_engine.dem_input import repository_for
from tin_engine.feature_input import ClassMap
from tin_engine.geodesy import LocalFrame
from tin_engine.io.geojson import collection_summary
from tin_engine.io.geopackage import distinct_values, layer_info, layer_summary
from tin_engine.io.models import TileFootprint
from tin_engine.io.repository import CacheRepository, open_geopackage, read_json
from tin_engine.sources import (
    OUT_ROOT_VARIABLE,
    RemoteSource,
    ZUnit,
    catalogue,
    plugin_source_files,
    user_source_candidates,
    user_source_file,
)

#: Segments per tile edge before reprojection: a chosen value at which an
#: edge's curvature is negligible against a cell.
SEGMENTS = 8
#: The most distinct values of a coded attribute listed: a chosen value,
#: above the 44 classes of CORINE level 3.
CLASS_VALUES = 100
Box = tuple[float, float, float, float]
_FROM_ENTRY = ("id", "kind", "values", "layer", "description", "credit", "licence_note", "cite")
#: GeoJSON's type names by a GeoPackage's declared ones, which are upper case.
_TYPES = {
    n.upper(): n
    for n in ("Point", "LineString", "Polygon", "MultiPoint", "MultiLineString", "MultiPolygon")
} | {"GEOMETRY": "Geometry", "GEOMETRYCOLLECTION": "GeometryCollection"}
_POLYGONS = frozenset({"Polygon", "MultiPolygon"})


class ClassField(BaseModel):
    """A coded class map's attribute as the file spells it, its distinct
    values sorted as text (at most `CLASS_VALUES`), and how many there are."""

    model_config = ConfigDict(frozen=True)

    attribute: str
    values: tuple[str, ...]
    distinct: int


class Fit(BaseModel):
    """An option a source can be given to: the command line's `command` and
    `option`, the MCP server's `tool` and `argument`, and with `--features`
    the class map to pass as `--features-map`. A necessary condition read
    from the header, not a promise: the reader may still refuse rows the
    header does not show (a missing or repeated `objectid`, a repeated
    `station`, a polygon with no area). What the features mean is the
    description's to say: land-cover polygons fit `catchment --lakes` too."""

    model_config = ConfigDict(frozen=True)

    command: str
    option: str
    tool: str
    argument: str
    class_map: str | None = None


class SourceReport(BaseModel):
    """One source as its headers describe it; the figures are None when nothing was read."""

    model_config = ConfigDict(frozen=True)

    id: str
    kind: str
    values: Literal["heights", "classes", "features"]
    path: Path | None
    present: bool
    description: str
    credit: str
    licence_note: str
    cite: tuple[str, ...]
    tiles: int = 0
    crs: str | None = None
    crs_name: str | None = None
    cell_m: tuple[float, float] | None = None
    z_unit: ZUnit | None = None
    z_unit_from: Literal["file", "catalogue", "assumed"] | None = None
    nodata: float | None = None
    coverage: Box | None = None
    box_covered: float | None = None
    problem: str | None = None
    #: "package", "user", or "plugin <distribution>".
    origin: str
    #: A vector source's header: its GeoJSON type names ("Geometry" for a
    #: GeoPackage column declared generic), feature count, attribute names in
    #: file order, and the values of each coded attribute that fits a class map.
    layer: str | None = None
    geometry: tuple[str, ...] = ()
    features: int | None = None
    attributes: tuple[str, ...] = ()
    classes: tuple[ClassField, ...] = ()
    fits: tuple[Fit, ...] = ()


class CatalogueReport(BaseModel):
    """The listing: the places read, and every source."""

    model_config = ConfigDict(frozen=True)

    user_file: str | None
    looked_at: tuple[str, ...]
    register_at: str | None
    user_file_problem: str | None
    cache: str | None
    out_root: str | None
    note: str | None
    sources: tuple[SourceReport, ...]
    plugin_problems: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class Facts:
    """What `fits` judges a source by: what it holds; for heights whether the
    CRS is projected (None when unknown); for features the file's suffix,
    geometry type names, attributes, whether it states its CRS (a GeoPackage
    always does), and the entry's class map ("" when none)."""

    values: Literal["heights", "classes", "features"]
    projected: bool | None
    suffix: str
    geometry: frozenset[str]
    attributes: tuple[str, ...]
    crs_stated: bool
    class_map: str


def _fit(command: str, option: str, class_map: str | None = None) -> Fit:
    """A fit, its MCP spelling derived as `mcp/schema.py` derives it."""
    tool, argument = command.replace("-", "_"), option.lstrip("-").replace("-", "_")
    return Fit(command=command, option=option, tool=tool, argument=argument, class_map=class_map)


def fits(facts: Facts, maps: Mapping[str, ClassMap]) -> tuple[Fit, ...]:
    """The options a source fits, each the acceptance condition of the reader
    the option calls: a DEM's for `--dem` (not `station-catchments`' for
    one in degrees), a class raster's for `--land-cover`; for features,
    `--features` under each of `maps` (only the entry's when it names one)
    whose attribute the file has (a GeoPackage's column case-insensitively,
    a GeoJSON property exactly), `catchment --lakes` for any polygons, and
    the GeoJSON readers, which need a `crs` member, by geometry and attribute.
    Reads no file."""
    if facts.values == "heights":
        commands = ("mesh", "catchment") + ("station-catchments",) * (facts.projected is not False)
        return tuple(_fit(c, "--dem") for c in commands)
    if facts.values == "classes":
        return (_fit("mesh", "--land-cover"),)
    fold: Callable[[str], str] = str.lower if facts.suffix == ".gpkg" else str
    names, shapes = {fold(a) for a in facts.attributes}, facts.geometry
    found = [
        _fit("mesh", "--features", name)
        for name, cmap in maps.items()
        if facts.class_map in ("", name) and fold(cmap.attribute) in names
    ]
    if shapes & (_POLYGONS | {"Geometry"}):
        found.append(_fit("catchment", "--lakes"))
    if facts.suffix != ".gpkg" and facts.crs_stated:
        if shapes == {"LineString"} and "objectid" in names:
            found += [_fit("catchment", "--rivers"), _fit("station-catchments", "--rivers")]
        if shapes and shapes <= _POLYGONS:
            found.append(_fit("station-catchments", "--lakes"))
            if "station" in names:
                found.append(_fit("station-catchments", "--reference"))
        if shapes == {"Point"} and "station" in names:
            found.append(_fit("station-catchments", "--stations"))
    return tuple(found)


def _lonlat(to_lonlat: Callable[[Any], Any], box: Box) -> shapely.Polygon:
    """`box`, densified, through `to_lonlat`: its outline in longitude and latitude."""
    x0, y0, x1, y1 = box
    step = max(x1 - x0, y1 - y0) / SEGMENTS or 1.0
    ring = shapely.segmentize(shapely.box(x0, y0, x1, y1), step)
    return shapely.Polygon(to_lonlat(np.asarray(ring.exterior.coords)))


def _cover(union: shapely.Geometry, box: Box | None) -> dict[str, object]:
    """The coverage box of `union` and the share of `box` it covers."""
    covered = None
    if box is not None:
        wanted = shapely.envelope(shapely.multipoints([box[:2], box[2:]]))
        # A point or line box has no area: covered or not.
        share = union.covers(wanted) if wanted.area == 0 else union.intersection(wanted).area
        covered = round(float(share) / (wanted.area or 1.0), 4)
    return {"coverage": tuple(round(v, 6) for v in union.bounds), "box_covered": covered}


def _figures(
    entry: RemoteSource, tiles: tuple[TileFootprint, ...], box: Box | None
) -> dict[str, object]:
    """The header figures of `tiles`: CRS, cell sizes, height unit, coverage
    and the share of `box` covered."""
    meta = tiles[0].meta
    to_lonlat = reprojector(meta.crs, "EPSG:4326")
    cells, shapes = [], []
    for tile in tiles:
        m = tile.meta
        x0, y0, x1, y1 = m.node_box()
        cells.append(LocalFrame.at(m.crs, ((x0 + x1) / 2, (y0 + y1) / 2)).cell_m(m))
        shapes.append(_lonlat(to_lonlat, (x0, y0, x1, y1)))
    heights = entry.values == "heights"
    stated = "catalogue" if entry.z_unit else "assumed"
    unit_from = stated if meta.vertical_unit_assumed else "file"
    return {
        "crs": crs_label(meta.crs),
        "crs_name": parse_crs(meta.crs).name,
        "cell_m": (round(min(cells), 1), round(max(cells), 1)),
        "z_unit": meta.z_unit if heights else None,
        "z_unit_from": unit_from if heights else None,
        "nodata": meta.nodata,
        **_cover(shapely.union_all(shapes), box),
    }


def _vector(
    path: Path, entry: RemoteSource, box: Box | None, maps: Mapping[str, ClassMap]
) -> dict[str, Any]:
    """The header figures of a vector source's file `path`, its fits, and the
    values of each coded attribute a fitting class map reads; its coverage is
    the extent box, not the features."""
    suffix, coded = path.suffix.lower(), {}
    if suffix == ".gpkg":
        with closing(open_geopackage(path)) as conn:
            layer = layer_info(conn, entry.layer)
            got = layer_summary(conn, layer)
            shapes = frozenset({_TYPES.get(got.geometry.upper(), got.geometry)})
            facts = Facts("features", None, suffix, shapes, got.attributes, True, entry.class_map)
            found = fits(facts, maps)
            for a in _coded(found, maps, got.attributes):
                coded[a] = distinct_values(conn, layer, a, CLASS_VALUES)
        crs, extent, count = layer.crs, got.extent, got.features
    else:
        doc = collection_summary(read_json(path))
        facts = Facts(
            "features", None, suffix, doc.geometry, doc.attributes, doc.crs_stated, entry.class_map
        )
        found = fits(facts, maps)
        for a in _coded(found, maps, doc.attributes):
            given = [(f.get("properties") or {}).get(a) for f in doc.rows]
            seen = sorted({str(v) for v in given if v is not None})
            coded[a] = (tuple(seen[:CLASS_VALUES]), len(seen))
        crs, extent, count = doc.crs, doc.extent, doc.features
    cover = {} if extent is None else _cover(_lonlat(reprojector(crs, "EPSG:4326"), extent), box)
    return {
        "crs": crs_label(crs),
        "crs_name": parse_crs(crs).name,
        "geometry": tuple(sorted(facts.geometry)),
        "features": count,
        "attributes": facts.attributes,
        "classes": [{"attribute": a, "values": v, "distinct": n} for a, (v, n) in coded.items()],
        "fits": found,
        **cover,
    }


def _coded(
    found: tuple[Fit, ...], maps: Mapping[str, ClassMap], attributes: tuple[str, ...]
) -> list[str]:
    """The attributes, as the file spells them, that a fitting coded class map reads."""
    read = {
        maps[f.class_map].attribute.lower()
        for f in found
        if f.class_map and maps[f.class_map].codes
    }
    return [a for a in attributes if a.lower() in read]


def _one(
    entry: RemoteSource,
    origin: str,
    cache: Path | None,
    box: Box | None,
    maps: Mapping[str, ClassMap],
) -> SourceReport:
    """`entry` with what its files or its cache say, and the options it fits;
    a read that fails is its `problem`, and a vector file not read fits none."""
    found: dict[str, object] = {}
    projected = None
    if entry.path is not None:
        where: Path | None = entry.path
        present = entry.path.exists()
    else:
        where = None if cache is None else cache / entry.id
        present = where is not None and (where / "manifest.json").is_file()
        if where is not None and entry.kind == "arcgis-features":
            where = where / "features.gpkg"
    if entry.path is not None or present:
        try:
            if entry.values == "features":
                assert where is not None
                found = _vector(where, entry, box, maps)
            else:
                if entry.path is not None:
                    repository = repository_for((entry.path,), z_unit=entry.read_unit)[0]
                else:
                    assert cache is not None
                    repository = CacheRepository(cache, entry.id, z_unit=entry.read_unit)
                tiles = repository.footprints()
                found = {"tiles": len(tiles), **_figures(entry, tiles, box)}
                projected = not tiles[0].meta.geographic
        except (OSError, ValueError, sqlite3.Error) as exc:
            found["problem"] = str(exc)
    if entry.values != "features":
        facts = Facts(entry.values, projected, "", frozenset(), (), False, entry.class_map)
        found["fits"] = fits(facts, maps)
    named = entry.model_dump(include=set(_FROM_ENTRY))
    return SourceReport.model_validate(
        {**found, **named, "present": present, "path": where, "origin": origin}
    )


def listing(environ: Mapping[str, str], cache: Path | None, box: Box | None) -> CatalogueReport:
    """The listing for this environment: the user's sources file it finds
    (a refused one is named, and the other sources are still listed), the
    plugins' files skipped, the cache, and every source in catalogue order;
    `box` (W S E N degrees) adds the share of it each source covers."""
    looked, user, plugins = (
        user_source_candidates(environ),
        user_source_file(environ),
        plugin_source_files(),
    )
    problem = None
    try:
        found = catalogue(user, plugins)
    except ValueError as exc:
        found, problem = catalogue(None, plugins), str(exc)
    # The XDG file when there is one: what a server started with only HOME reads.
    register = looked[-1] if looked else None
    maps = {name: load(name) for name in (*shipped_names(), "property")}
    sources = tuple(
        _one(entry, found.origins[key], cache, box, maps) for key, entry in found.sources.items()
    )
    notes = []
    if user is None:
        notes.append(
            f"No user sources file: looked at {', '.join(map(str, looked)) or 'no place'}. "
            f'To use a local DEM, add a [source.<id>] table with kind = "local", path, credit, '
            f"licence_note and description to {register or 'a sources.toml'}."
        )
    elif user != register:
        notes.append(f"A server started without RASPUTIN_DATA reads {register} instead of {user}.")
    if any(s.present and s.kind != "local" for s in sources):
        notes.append(
            "A remote source's figures describe its cached objects, which may hold only "
            "the part fetched."
        )
    if any(s.values == "features" and s.present for s in sources):
        notes.append("A vector source's coverage is its layer's extent box, not its features.")
    return CatalogueReport(
        user_file=None if user is None else str(user),
        looked_at=tuple(map(str, looked)),
        register_at=None if register is None else str(register),
        user_file_problem=problem,
        cache=None if cache is None else str(cache),
        out_root=environ.get(OUT_ROOT_VARIABLE),
        note=" ".join(notes) or None,
        sources=sources,
        plugin_problems=found.refused,
    )


__all__ = ["CatalogueReport", "ClassField", "Facts", "Fit", "SourceReport", "fits", "listing"]
