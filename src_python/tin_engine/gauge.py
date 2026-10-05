"""The gauge on NVE's mapped river line (increment 29, PR 2).

`docs/increments/29-nve-reference-catchments.md`, "Placing the gauge". A
station is placed where it physically is: on the nearest line of the best
tier (its own watercourse number, a prefix of it, its river's name, any),
at the foot `P` of the station on that line. No rule here looks at an area
or a flow count. The reach is the chain of that line's `elvid` round `P`,
`reach_up` metres upstream and `U + 100` m downstream, stopped at a fork or
where the river ends. `lake_seed` (PR 4, "Lake gauges: the lake is the
seed") says which gauges are seeded with their lake instead: by containment
and one distance, never an area. Pure shapely: no DEM, no file.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from itertools import pairwise
from typing import Literal

import shapely
from pydantic import BaseModel, ConfigDict
from shapely.geometry import LineString, Point
from shapely.ops import substring

from tin_engine.io.rivers import RiverSegment
from tin_engine.io.station_set import Lake

#: Metres: an end within this of the next segment's start continues the chain.
JOIN_M = 1.0
#: Metres: `U` is never below three DTM10 cells.
U_FLOOR_M = 30.0
#: Metres: the chain runs `U + DOWN_EXTRA_M` below `P`.
DOWN_EXTRA_M = 100.0
#: Metres from `P` within which another river's line flags `confluence_near`.
CONFLUENCE_M = 100.0
#: Metres (the river file's CRS) from a station to the lake its `P` lies in,
#: within which a gauge placed on a lake line is seeded with the lake: `U`'s floor.
LAKE_GAP_M = 30.0

Tier = Literal["number", "prefix", "name", "any"]
_TIERS: tuple[Tier, ...] = ("number", "prefix", "name", "any")


class Gauge(BaseModel):
    """A station point, its watercourse number and its river's name, if known."""

    model_config = ConfigDict(frozen=True)

    x: float
    y: float
    watercourse: str | None = None
    river: str | None = None


class Reach(BaseModel):
    """The mapped reach in the request's `seed_crs`, in downstream order, `at`
    metres along it to `P`, `uncertainty` (`U`) and the burn's half-width."""

    model_config = ConfigDict(frozen=True)

    line: tuple[tuple[float, float], ...]
    at: float
    uncertainty: float
    corridor: float = 30.0


class Placement(BaseModel):
    """Where the gauge went, and why there (step 5's flags, for the report)."""

    model_config = ConfigDict(frozen=True)

    reach: Reach
    placed_on: Tier
    distance_m: float
    position: tuple[float, float]
    elvid: str | None
    objectid: int
    lake: bool
    confluence_near: bool
    reach_up_m: float
    reach_down_m: float
    reach_fork: bool


class LakeSeed(BaseModel):
    """A gauge seeded with its lake: `rule` says why, `point` is the station
    (`inside`) or `P` (`lake_line`), `lakes` all lakes containing `point`,
    and `distance_m` the station's distance to them (0 inside)."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    rule: Literal["inside", "lake_line"]
    point: tuple[float, float]
    lakes: tuple[Lake, ...]
    distance_m: float


def _tier(segment: RiverSegment, gauge: Gauge) -> int:
    number, station = segment.vassdragsnr, gauge.watercourse
    if number is not None and station is not None:
        if number == station:
            return 0
        # A leading part with its own dot shares the main number; a bare main
        # number ("002") is never a prefix (step 2, pinned 2026-10-05).
        if "." in number and station.startswith(number):
            return 1
    if gauge.river is not None and segment.name == gauge.river:
        return 2
    return 3


def _walk(
    start: RiverSegment, same: Sequence[RiverSegment], need: float, down: bool
) -> tuple[list[RiverSegment], bool]:
    """Segments joined beyond `start`, downstream or upstream, until `need`
    metres are had; stopped where none continues, or at a fork (flagged)."""
    out: list[RiverSegment] = []
    cur, have, seen = start, 0.0, {start.objectid}
    while have < need:
        end = cur.line[-1] if down else cur.line[0]
        nxt = [
            s
            for s in same
            if s.objectid not in seen
            and math.dist(end, s.line[0] if down else s.line[-1]) <= JOIN_M
        ]
        if len(nxt) != 1:
            return out, len(nxt) > 1
        cur = nxt[0]
        seen.add(cur.objectid)
        out.append(cur)
        have += LineString(cur.line).length
    return out, False


def _length(coords: Sequence[tuple[float, float]]) -> float:
    return sum(math.dist(a, b) for a, b in pairwise(coords))


def place(
    gauge: Gauge,
    segments: Sequence[RiverSegment],
    *,
    map_radius: float = 500.0,
    reach_up: float = 1000.0,
) -> Placement | None:
    """The placement of `gauge` on `segments` (same CRS), or None when no
    segment lies within `map_radius` of it."""
    station = Point(gauge.x, gauge.y)
    near = [(LineString(s.line).distance(station), s) for s in segments]
    near = [(d, s) for d, s in near if d <= map_radius]
    if not near:
        return None
    tier, d, chosen = min(((_tier(s, gauge), d, s) for d, s in near), key=lambda t: t[:2])
    line = LineString(chosen.line)
    along = line.project(station)
    p = line.interpolate(along)
    u = min(max(d, U_FLOOR_M), map_radius)
    same = [s for s in segments if s.elvid == chosen.elvid]
    ups, fork_up = _walk(chosen, same, reach_up - along, down=False)
    downs, fork_down = _walk(chosen, same, u + DOWN_EXTRA_M - (line.length - along), down=True)
    coords: list[tuple[float, float]] = []
    for s in [*reversed(ups), chosen, *downs]:
        if s is chosen:  # the metres of the chain above the chosen line's start
            offset = _length([*coords, s.line[0]])
        coords += [v for v in s.line if not coords or v != coords[-1]]
    full = LineString(coords)
    s_p = offset + along
    start, end = max(0.0, s_p - reach_up), min(full.length, s_p + u + DOWN_EXTRA_M)
    cut = substring(full, start, end)
    other = [s for s in segments if s.elvid != chosen.elvid]
    return Placement(
        reach=Reach(line=tuple((x, y) for x, y in cut.coords), at=s_p - start, uncertainty=u),
        placed_on=_TIERS[tier],
        distance_m=d,
        position=(p.x, p.y),
        elvid=chosen.elvid,
        objectid=chosen.objectid,
        lake=chosen.kind == "lake",
        confluence_near=any(shapely.distance(LineString(s.line), p) <= CONFLUENCE_M for s in other),
        reach_up_m=s_p - start,
        reach_down_m=end - s_p,
        reach_fork=fork_up or fork_down,
    )


def lake_seed(
    gauge: Gauge, placement: Placement | None, lakes: Sequence[Lake], *, gap: float = LAKE_GAP_M
) -> LakeSeed | None:
    """The lake seed of `gauge`, or None for PR 2's river path: `inside` when
    the station lies strictly inside a lake; else `lake_line` when it was
    placed on a lake line whose `P` lies in a lake at most `gap` from it."""
    if not (math.isfinite(gauge.x) and math.isfinite(gauge.y)):
        return None
    station = Point(gauge.x, gauge.y)
    holding = tuple(lk for lk in lakes if lk.polygon.contains(station))
    if holding:
        return LakeSeed(rule="inside", point=(gauge.x, gauge.y), lakes=holding, distance_m=0.0)
    if placement is None or not placement.lake:
        return None
    p = Point(placement.position)
    holding = tuple(lk for lk in lakes if lk.polygon.contains(p))
    if not holding:
        return None
    distance = min(lk.polygon.distance(station) for lk in holding)
    if distance > gap + 1e-6:  # the corridor's slack
        return None
    return LakeSeed(rule="lake_line", point=placement.position, lakes=holding, distance_m=distance)


__all__ = ["LAKE_GAP_M", "Gauge", "LakeSeed", "Placement", "Reach", "lake_seed", "place"]
