"""Test support for increment 16b: feature documents, domains, and the bit oracle.

`docs/increments/16b-terrain-polygons.md`. The oracle here is built from the
input features and the class map, never from the producer's records (05b,
"Guarantee 15's oracle"): a noded constraint edge's expected mask is the union
of the masks of every feature whose boundary, as given in the source and
moved by pyproj directly, lies along that edge, within the noder's snap.

`tin_engine.feature_input` and `tin_engine.chains` are imported lazily, so a
suite importing this module still collects before they exist.
"""

from __future__ import annotations

import importlib
import json
import sqlite3
from collections.abc import Callable, Iterable, Mapping, Sequence
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import shapely
from pyproj import Transformer
from shapely.geometry import LineString, Polygon, mapping
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import orient

from tin_engine.domain import DomainPolygon
from tin_engine.features import DEFAULT_VOCABULARY

UTM33 = "EPSG:25833"
SNAP = 1e-3  # cli.DEFAULT_SNAP_SPACING: the noder's grid, in the DEM's CRS
#: I2's bound: a crossing point lies within this of the domain's boundary.
CROSSING = 1e-6
#: Where the synthetic partitions sit: UTM-shaped, off every round number.
X0, Y0 = 500_000.3, 6_600_000.7

Ring = Sequence[tuple[float, float]]


def feature_input() -> ModuleType:
    return importlib.import_module("tin_engine.feature_input")


def chains_module() -> ModuleType:
    return importlib.import_module("tin_engine.chains")


def at(x: float, y: float) -> tuple[float, float]:
    return (X0 + x, Y0 + y)


def square(x0: float, y0: float, x1: float, y1: float) -> Polygon:
    return Polygon([at(x0, y0), at(x1, y0), at(x1, y1), at(x0, y1)])


def domain_of(polygon: Polygon, crs: str = UTM33) -> DomainPolygon:
    return DomainPolygon(polygon=orient(polygon, sign=1.0), crs=crs)


@dataclass(frozen=True)
class Feat:
    """One GeoJSON feature: its id, geometry and properties."""

    fid: Any
    geometry: BaseGeometry | dict[str, Any] | None
    properties: Mapping[str, Any]


def write_geojson(path: Path, features: Iterable[Feat], crs: str | None = UTM33) -> Path:
    """A `FeatureCollection`; `crs` is its `crs` member (None: none, RFC 7946)."""
    out: list[dict[str, Any]] = []
    for f in features:
        geometry = f.geometry
        if isinstance(geometry, BaseGeometry):
            geometry = mapping(geometry)
        doc: dict[str, Any] = {
            "type": "Feature",
            "geometry": geometry,
            "properties": dict(f.properties),
        }
        if f.fid is not None:
            doc["id"] = f.fid
        out.append(doc)
    collection: dict[str, Any] = {"type": "FeatureCollection", "features": out}
    if crs is not None:
        collection["crs"] = {"type": "name", "properties": {"name": crs}}
    path.write_text(json.dumps(collection))
    return path


def source(path: Path, map_name: str = "property", **kwargs: Any) -> Any:
    fi = feature_input()
    return fi.FeatureSource(path=path, class_map=fi.CLASS_MAPS[map_name], **kwargs)


def open_one(
    path: Path, domain: DomainPolygon, map_name: str = "property", dem_crs: str = UTM33, **kw: Any
) -> Any:
    """`open_features` on one source."""
    fi = feature_input()
    request = fi.FeatureRequest(sources=(source(path, map_name, **kw),))
    return fi.open_features(request, domain, dem_crs)


# ----------------------------------------------------------------- the engine


@dataclass(frozen=True)
class Noded:
    """The noded constraint graph: each undirected edge as its two node
    coordinates, sorted, and its mask; and which edges came from `Outer`/`Hole`."""

    masks: dict[tuple[tuple[float, float], tuple[float, float]], int]
    boundary: set[tuple[tuple[float, float], tuple[float, float]]]
    status: str


def start(domain: DomainPolygon, features: Sequence[Any]) -> Any:
    return chains_module().start_chains(domain, features, DEFAULT_VOCABULARY)


def run_engine(started: Any) -> Noded:
    """`StartChains` through `cli._engine` (`build_pslg` -> `node` ->
    `triangulate`) at the CLI's snap; the graph as coordinates."""
    import tin_engine.cli as cli

    chains = [
        ([int(i) for i in indices], cli.ROLES[role], int(mask))
        for indices, role, mask in started.chains
    ]
    run = cli._engine(np.asarray(started.vertices, dtype=np.float64), chains, True, SNAP)
    assert run.noded is not None and run.mesh is not None, f"{run.status}: {run.message}"
    xy = np.asarray(run.noded.vertices)

    def key(a: int, b: int) -> tuple[tuple[float, float], tuple[float, float]]:
        p, q = (float(xy[a, 0]), float(xy[a, 1])), (float(xy[b, 0]), float(xy[b, 1]))
        return (p, q) if p <= q else (q, p)

    masks = {key(a, b): m for (a, b), m in cli._chain_masks(run.noded).items()}
    boundary: set[tuple[tuple[float, float], tuple[float, float]]] = set()
    for c, chain in enumerate(run.noded.chains):
        if chain.role in cli.CORE_CLOSED_ROLES:
            walk = [int(i) for i in run.noded.indices_of(c)]
            boundary |= {key(walk[k], walk[(k + 1) % len(walk)]) for k in range(len(walk))}
    return Noded(masks=masks, boundary=boundary, status=run.status)


# ----------------------------------------------------------------- the oracle


def boundary_lines(geometry: BaseGeometry) -> list[LineString]:
    """Every ring of every polygon, and every line, as a line: the linework a
    feature contributes before any clip."""
    lines: list[LineString] = []
    for part in shapely.get_parts(geometry):
        if isinstance(part, Polygon):
            lines += [LineString(r.coords) for r in (part.exterior, *part.interiors)]
        elif isinstance(part, LineString):
            lines.append(part)
    return lines


def moved(geometry: BaseGeometry, src: str, dst: str = UTM33) -> BaseGeometry:
    """`geometry` from `src` to `dst` vertex by vertex, by pyproj directly."""
    if src == dst:
        return geometry
    t = Transformer.from_crs(src, dst, always_xy=True)

    def apply(xy: np.ndarray) -> np.ndarray:
        x, y = t.transform(xy[:, 0], xy[:, 1])
        return np.column_stack([x, y])

    return shapely.transform(geometry, apply)


@dataclass(frozen=True)
class Truth:
    """A feature as the oracle sees it: its linework in the DEM's CRS, and the
    mask its class map gives it."""

    lines: BaseGeometry
    mask: int


def truths(features: Iterable[tuple[BaseGeometry, int]], src: str = UTM33) -> list[Truth]:
    return [
        Truth(shapely.MultiLineString(boundary_lines(moved(g, src))), mask)
        for g, mask in features
        if mask is not None
    ]


def expected_masks(noded: Noded, truth: Sequence[Truth], near: float = SNAP) -> dict[Any, int]:
    """I1's oracle: for each noded edge, the union of the masks of every
    feature whose linework lies along it (both ends and the midpoint within
    `near`, the snap). Built from the features, never from the noder."""
    tree = shapely.STRtree([t.lines for t in truth])
    out: dict[Any, int] = {}
    for edge in noded.masks:
        (ax, ay), (bx, by) = edge
        probe = [
            shapely.Point(ax, ay),
            shapely.Point(bx, by),
            shapely.Point((ax + bx) / 2, (ay + by) / 2),
        ]
        mask = 0
        for i in tree.query(LineString(edge).buffer(near)):
            lines = truth[int(i)].lines
            if all(lines.distance(p) <= near for p in probe):
                mask |= truth[int(i)].mask
        out[edge] = mask
    return out


def mismatches(noded: Noded, expected: Mapping[Any, int]) -> list[tuple[Any, int, int]]:
    return [(e, m, expected[e]) for e, m in noded.masks.items() if m != expected[e]]


# ---------------------------------------------------------- the committed extract


def extract_features(path: Path, table: str) -> list[tuple[int, BaseGeometry, str]]:
    """The committed extract's rows, `(OBJECTID, geometry, Code_18)`, read
    with sqlite3 and shapely only: its blobs have the 8-byte header of flags
    `0x01` (no envelope), which `extract.py` writes."""
    with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)) as con:
        rows = con.execute(
            f"SELECT OBJECTID, Shape, Code_18 FROM {table} ORDER BY OBJECTID"
        ).fetchall()
    for _, data, _ in rows:
        assert data[:2] == b"GP" and data[3] == 0x01
    return [(pk, shapely.from_wkb(data[8:]), code) for pk, data, code in rows]


def corine_mask(code: str) -> int:
    """R4's `corine` map, written here from the design's text."""
    water = code in {"511", "512", "521", "522", "523"}
    return DEFAULT_VOCABULARY.mask("land_cover", *(("water",) if water else ()))


def all_vertices(lines: Iterable[LineString]) -> np.ndarray:
    return (
        np.concatenate([np.asarray(line.coords) for line in lines]) if lines else np.zeros((0, 2))
    )


def within(polygon: Polygon, xy: np.ndarray) -> np.ndarray:
    """Distance of each point to `polygon` (0 inside or on it)."""
    return np.asarray(shapely.distance(polygon, shapely.points(xy)))


def on_boundary(polygon: Polygon, xy: np.ndarray) -> np.ndarray:
    return np.asarray(shapely.distance(polygon.boundary, shapely.points(xy)))


Check = Callable[[Any], None]
