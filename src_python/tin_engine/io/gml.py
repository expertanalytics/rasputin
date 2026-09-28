"""OGR-written GML2 features, read with the standard library (increment 16b, Q6).

`docs/increments/16b-terrain-polygons.md` R1 and R2. The reader is strict (Ola,
2026-09-28): a document that is not well-formed XML is refused, naming the
parser's line, and nothing is recovered from it. It takes a binary stream and
opens nothing.

A feature is an element with a ``fid`` attribute directly inside a
``gml:featureMember``. Its geometry is the first GML2 geometry inside it,
coordinates from ``gml:coordinates`` as ``x,y`` pairs (longitude first, as OGR
writes EPSG:4326), a third coordinate dropped.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import IO, Any

from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon
from shapely.geometry.base import BaseGeometry

GML = "{http://www.opengis.net/gml}"
_GEOMETRIES = frozenset(
    GML + name for name in ("Polygon", "MultiPolygon", "LineString", "MultiLineString", "Point")
)


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
    its child element named ``attribute``, None when absent."""
    try:
        root = ET.parse(stream).getroot()
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
