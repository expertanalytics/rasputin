"""A catchment joined from a published network of sub-catchment units.

Each unit names the unit it drains to; the catchment of a unit is that unit
and every unit upstream of it, joined. Pure: no file is read here, and every
refusal is a :class:`NetworkError` in plain words.

`docs/increments/85-network-catchment.md`, *6. Constants*.
"""

from __future__ import annotations

from collections import Counter, deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import shapely
from shapely.geometry import Point, Polygon
from shapely.geometry.base import BaseGeometry

from tin_engine.geodesy import LocalFrame, area_m2

#: How far the joined area may be from the stated upstream area, as a
#: fraction of the stated one; the stated areas are rounded per unit.
AREA_TOLERANCE = 0.01
#: A unit is mostly outside the terrain catchment below this share inside it:
#: a chosen threshold, for the report only.
MOSTLY_OUTSIDE = 0.5


class NetworkError(ValueError):
    """A network, seed or join refused, in plain words."""


@dataclass(frozen=True, slots=True)
class Link:
    """One unit: its own ID, the ID it drains to, the area upstream of its
    outlet (its own included) and its name, when the file gives them."""

    id: str
    down: str
    upstream_area_m2: float | None
    name: str | None


@dataclass(frozen=True, slots=True)
class Network:
    """The links by own ID, and each down value that names no unit (the
    sea, a border, the file's edge) with the count of units draining to it."""

    links: Mapping[str, Link]
    ends: Mapping[str, int]

    @classmethod
    def build(cls, links: Sequence[Link]) -> Network:
        """The network of `links`; a null, empty or repeated own ID is refused."""
        by_id: dict[str, Link] = {}
        for link in links:
            if not link.id:
                raise NetworkError("a unit has no ID (null or empty)")
            if link.id in by_id:
                raise NetworkError(f"the ID {link.id} is given to more than one unit")
            by_id[link.id] = link
        ends = Counter(k.down for k in by_id.values() if k.down not in by_id)
        return cls(by_id, dict(ends))

    def upstream(self, unit: str) -> tuple[str, ...]:
        """`unit` and every unit draining to it, directly or not, breadth
        first. A unit reached twice means `unit` lies on a cycle: refused,
        naming the units on it."""
        above: dict[str, list[str]] = {}
        for link in self.links.values():
            above.setdefault(link.down, []).append(link.id)
        found, seen, queue = [unit], {unit}, deque([unit])
        while queue:
            for child in above.get(queue.popleft(), []):
                if child in seen:
                    cycle, step = [unit], self.links[unit].down
                    while step != unit:
                        cycle.append(step)
                        step = self.links[step].down
                    raise NetworkError(f"the units {', '.join(cycle)} drain in a cycle")
                seen.add(child)
                found.append(child)
                queue.append(child)
        return tuple(found)


def unit_id(value: object) -> str:
    """An attribute value as an ID's text: an integral float as an integer
    (7.0 is "7"), so numeric and text IDs compare equal; null is ""."""
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def home_unit(point: Point, covering: Sequence[tuple[str, BaseGeometry]]) -> str:
    """The one unit among `covering` whose polygon covers `point`; none, or
    two (the point on their border), is refused."""
    found = [own for own, geometry in covering if geometry.covers(point)]
    if not found:
        raise NetworkError("the seed lies in no unit of the network")
    if len(found) > 1:
        raise NetworkError(f"the seed lies on the border of the units {', '.join(found)}")
    return found[0]


def outlet_distance_m(seed: BaseGeometry, outlet: Point, crs: str) -> float:
    """Metres from `seed` (a point or a lake) to `outlet`, both in `crs`;
    in degrees, measured in the local frame at the outlet."""
    frame = LocalFrame.at(crs, (outlet.x, outlet.y))
    moved = shapely.transform(seed, frame.forward)
    return float(moved.distance(Point(frame.forward([(outlet.x, outlet.y)])[0])))


@dataclass(frozen=True, slots=True)
class Joined:
    """The joined polygon with its holes, its area with the holes
    subtracted, and the holes' count and area."""

    polygon: Polygon
    area_m2: float
    holes: int
    holes_m2: float


def join(units: Sequence[BaseGeometry], crs: str) -> Joined:
    """The union of `units`, which must be one polygon; parts that do not
    touch are refused with their areas."""
    parts = shapely.get_parts(shapely.union_all(list(units)))
    if len(parts) != 1:
        areas = ", ".join(f"{area_m2(p, crs) / 1e6:.6f} km2" for p in parts)
        raise NetworkError(f"the units join into {len(parts)} parts that do not touch: {areas}")
    polygon = parts[0]
    holes = [Polygon(ring) for ring in polygon.interiors]
    return Joined(polygon, area_m2(polygon, crs), len(holes), sum(area_m2(h, crs) for h in holes))


def check_area(joined_m2: float, stated_m2: float) -> None:
    """Refuse a joined area more than `AREA_TOLERANCE` off the stated one."""
    if abs(joined_m2 - stated_m2) > AREA_TOLERANCE * stated_m2:
        raise NetworkError(
            f"the joined units cover {joined_m2 / 1e6:.6f} km2 but the stated upstream area "
            f"is {stated_m2 / 1e6:.6f} km2 ({100 * (joined_m2 / stated_m2 - 1):+.2f} %): "
            "the network file lacks units, or its links are wrong"
        )


@dataclass(frozen=True, slots=True)
class TerrainCheck:
    """The network catchment against the terrain route's, in km2, and the
    units mostly outside the terrain catchment as (ID, own area in km2)."""

    terrain_area_km2: float
    shared_km2: float
    network_only_km2: float
    terrain_only_km2: float
    units_mostly_outside: tuple[tuple[str, float], ...]


def compare(
    network: Polygon, terrain: Polygon, units: Mapping[str, BaseGeometry], crs: str
) -> TerrainCheck:
    """The areas shared and apart of `network` and `terrain`, and each of
    `units` with under `MOSTLY_OUTSIDE` of its area inside `terrain`."""

    def km2(geometry: BaseGeometry) -> float:
        return area_m2(geometry, crs) / 1e6

    outside = tuple(
        (own, km2(g))
        for own, g in units.items()
        if km2(g.intersection(terrain)) < MOSTLY_OUTSIDE * km2(g)
    )
    return TerrainCheck(
        km2(terrain),
        km2(network.intersection(terrain)),
        km2(network.difference(terrain)),
        km2(terrain.difference(network)),
        outside,
    )
