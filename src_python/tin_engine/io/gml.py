"""GML features, read with the standard library: OGR-written GML2 whole
(`read_gml`), and GML 3.2 surfaces streamed one feature at a time
(`iter_gml32`).

The reader is strict: a document that is not well-formed XML is refused,
naming the parser's line, and nothing is recovered from it, so a damaged file
cannot lose features without a word. It takes a binary stream and opens
nothing.

A feature is an element with a ``fid`` attribute directly inside a
``gml:featureMember``. Its geometry is the first GML2 geometry inside it,
coordinates from ``gml:coordinates`` as ``x,y`` pairs (longitude first, as OGR
writes EPSG:4326), a third coordinate dropped. A document declaring the GML
3.2 namespace is refused before its tree is built.

`iter_gml32` keeps one ``gml:featureMember`` in memory at a time. A feature's
geometry is a ``gml:Polygon``, ``gml:Surface`` or ``gml:MultiSurface`` of
``gml:posList`` rings, stored as served: not repaired, reoriented or
validated. Curves, points, other surface kinds, ``gml:pos`` rings, arcs and an
``srsDimension`` that is not a whole number of 2 or more are refused in words.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import IO, Any

import numpy as np
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon
from shapely.geometry.base import BaseGeometry

GML = "{http://www.opengis.net/gml}"
_GEOMETRIES = frozenset(
    GML + name for name in ("Polygon", "MultiPolygon", "LineString", "MultiLineString", "Point")
)
GML32 = "{http://www.opengis.net/gml/3.2}"
_NIL = "{http://www.w3.org/2001/XMLSchema-instance}nil"
_SURFACES = frozenset(GML32 + n for n in ("Polygon", "Surface", "MultiSurface"))
_NOT_READ = frozenset(
    GML32 + n
    for n in ("Point", "MultiPoint", "LineString", "Curve", "MultiCurve", "CompositeCurve",
              "OrientableCurve", "LinearRing", "Ring", "MultiGeometry", "Solid", "MultiSolid",
              "CompositeSurface", "OrientableSurface", "PolyhedralSurface",
              "TriangulatedSurface", "Tin")
)  # fmt: skip


class GmlError(ValueError):
    """A GML document or feature this reader refuses."""


@dataclass(frozen=True, slots=True)
class GmlFeature:
    fid: str
    geometry: BaseGeometry | None
    value: str | None


@dataclass(frozen=True, slots=True)
class GmlDocument:
    """The geometries' one ``srsName`` (None when none gives one), and the
    features in document order."""

    crs: str | None
    features: tuple[GmlFeature, ...]


def read_gml(stream: IO[bytes], attribute: str) -> GmlDocument:
    """Every feature of the document in ``stream``; its ``value`` is the text of
    its child element named ``attribute``, None when absent. A document that
    declares the GML 3.2 namespace is refused at the declaration."""
    try:
        events = ET.iterparse(stream, events=("start-ns",))
        for _, (_, uri) in events:
            if uri == GML32[1:-1]:
                raise GmlError(
                    "a GML 3.2 file: convert it with `rasputin gml-to-gpkg` and give the .gpkg"
                )
        root = events.root  # type: ignore[attr-defined]  # documented, absent from typeshed
    except ET.ParseError as exc:
        line, _ = exc.position
        raise GmlError(f"not well-formed XML at line {line}: {exc}") from exc
    features: list[GmlFeature] = []
    names: set[str] = set()
    for member in root.iter(GML + "featureMember"):
        for element in member:
            fid = element.get("fid")
            if fid is None:
                continue
            found = next((e for e in element.iter() if e.tag in _GEOMETRIES), None)
            value = next((c.text for c in element if c.tag.rpartition("}")[2] == attribute), None)
            try:
                geometry = None if found is None else _geometry(found)
            except (ValueError, IndexError) as exc:
                raise GmlError(f"feature {fid}: unreadable coordinates ({exc})") from exc
            if found is not None and found.get("srsName") is not None:
                names.add(str(found.get("srsName")))
            features.append(GmlFeature(fid=fid, geometry=geometry, value=value))
    if len(names) > 1:
        raise GmlError(f"the geometries name {len(names)} CRSs: {sorted(names)}")
    return GmlDocument(crs=names.pop() if names else None, features=tuple(features))


def _coordinates(element: ET.Element) -> list[tuple[float, float]]:
    """The first ``gml:coordinates`` inside ``element``, as ``(x, y)`` pairs."""
    node = element.find(f".//{GML}coordinates")
    if node is None or node.text is None:
        raise ValueError("no gml:coordinates")
    points = []
    for pair in node.text.split():
        parts = pair.split(",")
        if len(parts) < 2:
            raise ValueError(f"{pair!r} is not x,y")
        points.append((float(parts[0]), float(parts[1])))
    return points


def _polygon(element: ET.Element) -> Polygon:
    outer = element.find(f"{GML}outerBoundaryIs")
    if outer is None:
        raise ValueError("a polygon without gml:outerBoundaryIs")
    holes = [_coordinates(r) for r in element.findall(f"{GML}innerBoundaryIs")]
    return Polygon(_coordinates(outer), holes)


def _geometry(element: ET.Element) -> Any:
    kind = element.tag.removeprefix(GML)
    if kind == "Polygon":
        return _polygon(element)
    if kind == "MultiPolygon":
        return MultiPolygon([_polygon(p) for p in element.iter(GML + "Polygon")])
    if kind == "MultiLineString":
        return MultiLineString([_coordinates(g) for g in element.iter(GML + "LineString")])
    points = _coordinates(element)
    return LineString(points) if kind == "LineString" else Point(points[0])


@dataclass(frozen=True, slots=True)
class Gml32Feature:
    """A GML 3.2 feature: its ``gml:id``, its element's local name, its one
    surface (None when it has none) and that surface's ``srsName``, and its
    leaf properties by local name in document order (None when empty or nil)."""

    gml_id: str
    type_name: str
    geometry: Polygon | MultiPolygon | None
    srs_name: str | None
    values: Mapping[str, str | None]


def iter_gml32(stream: IO[bytes]) -> Iterator[Gml32Feature]:
    """The feature inside each ``gml:featureMember`` of the document in
    ``stream``, in order; the tree is cleared after each, so memory holds one.
    A document with no GML 3.2 ``gml:featureMember`` is refused."""
    members = 0
    try:
        events = ET.iterparse(stream, events=("start", "end"))
        _, root = next(events)
        for event, element in events:
            if event == "end" and element.tag == GML32 + "featureMember":
                members += 1
                yield from (_feature(child) for child in element)
                root.clear()
    except ET.ParseError as exc:
        line, _ = exc.position
        raise GmlError(f"not well-formed XML at line {line}: {exc}") from exc
    if not members:
        raise GmlError("no GML 3.2 gml:featureMember: not a GML 3.2 feature collection")


def _name(element: ET.Element) -> str:
    return element.tag.rpartition("}")[2]


def _feature(element: ET.Element) -> Gml32Feature:
    gml_id = element.get(GML32 + "id", "")
    found: list[ET.Element] = []
    values: dict[str, str | None] = {}
    for prop in element:
        if len(prop) == 0:
            text = (prop.text or "").strip()
            values[_name(prop)] = None if not text or prop.get(_NIL) in ("true", "1") else text
            continue
        for child in prop:
            if child.tag in _NOT_READ:
                raise GmlError(f"feature {gml_id}: surfaces only; gml:{_name(child)} is not read")
            if child.tag in _SURFACES:
                found.append(child)
    if len(found) > 1:
        raise GmlError(f"feature {gml_id}: {len(found)} surface properties; one is read")
    geometry = srs_name = None
    if found:
        srs_name = found[0].get("srsName")
        parts = _parts(found[0], gml_id, 2)
        if found[0].tag == GML32 + "MultiSurface" or len(parts) > 1:
            geometry = MultiPolygon(parts)
        else:
            geometry = parts[0]
    return Gml32Feature(gml_id, _name(element), geometry, srs_name, values)


def _parts(element: ET.Element, gml_id: str, dimension: int) -> list[Polygon]:
    """The polygons of a ``gml:Polygon``, ``gml:Surface`` (its patches) or
    ``gml:MultiSurface`` (its members, a member surface's patches flattened)."""
    dimension = _dimension(element, gml_id, dimension)
    if element.tag == GML32 + "Polygon":
        return [_polygon32(element, gml_id, dimension)]
    if element.tag == GML32 + "Surface":
        patches = element.findall(f"{GML32}patches/{GML32}PolygonPatch")
        return [_polygon32(p, gml_id, dimension) for p in patches]
    if element.tag != GML32 + "MultiSurface":
        raise GmlError(f"feature {gml_id}: surfaces only; gml:{_name(element)} is not read")
    members = [m for m in element if _name(m) in ("surfaceMember", "surfaceMembers")]
    return [p for m in members for g in m for p in _parts(g, gml_id, dimension)]


def _polygon32(element: ET.Element, gml_id: str, dimension: int) -> Polygon:
    exterior = element.find(GML32 + "exterior")
    if exterior is None:
        raise GmlError(f"feature {gml_id}: a polygon without gml:exterior")
    holes = [_ring32(r, gml_id, dimension) for r in element.findall(GML32 + "interior")]
    return Polygon(_ring32(exterior, gml_id, dimension), holes)


def _only(node: ET.Element, want: str, gml_id: str) -> ET.Element:
    """The one child of ``node``, which must be a ``gml:<want>``."""
    children = list(node)
    if len(children) != 1 or children[0].tag != GML32 + want:
        names = ", ".join(f"gml:{_name(c)}" for c in children) or "nothing"
        raise GmlError(f"feature {gml_id}: a ring of {names}; one gml:{want} is read")
    return children[0]


def _dimension(element: ET.Element, gml_id: str, default: int) -> int:
    """The ``srsDimension`` of ``element``, ``default`` when absent."""
    text = element.get("srsDimension")
    if text is None:
        return default
    if not (text.isascii() and text.isdigit() and int(text) >= 2):
        raise GmlError(f"feature {gml_id}: srsDimension {text!r} is not a whole number >= 2")
    return int(text)


def _ring32(boundary: ET.Element, gml_id: str, dimension: int) -> np.ndarray:
    """The ``(x, y)`` positions of the one ``gml:LinearRing/gml:posList`` in ``boundary``."""
    ring = _only(boundary, "LinearRing", gml_id)
    pos_list = _only(ring, "posList", gml_id)
    dimension = _dimension(pos_list, gml_id, _dimension(ring, gml_id, dimension))
    try:
        numbers = np.array((pos_list.text or "").split(), dtype=np.float64)
    except ValueError as exc:
        raise GmlError(f"feature {gml_id}: unreadable gml:posList ({exc})") from exc
    if len(numbers) % dimension or len(numbers) < 4 * dimension:
        raise GmlError(
            f"feature {gml_id}: a gml:posList of {len(numbers)} numbers is not four or more"
            f" positions of {dimension}"
        )
    if not np.isfinite(numbers).all():
        raise GmlError(f"feature {gml_id}: a coordinate that is not finite")
    return numbers.reshape(-1, dimension)[:, :2]
