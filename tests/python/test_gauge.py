"""`gauge.place`: the gauge on NVE's mapped river line (increment 29, PR 2).

`docs/increments/29-nve-reference-catchments.md`, "Placing the gauge" (steps
1 to 5 and the `Reach` model) and "The red suites", PR 2's `test_gauge.py`.
Pure shapely, no DEM: segments are built as PR 3's `RiverSegment` directly,
or read through PR 3's `read_segments` where the reader's dropping of exact
copies is the point.

Interface pinned here (the design names `Gauge`, `Placement`, `Reach` and
`place(gauge, segments, *, map_radius=500.0, reach_up=1000.0)`; the field
names below that it leaves open are chosen here and listed in the handback):

- `Gauge(x=, y=, watercourse=None, river=None)`, frozen.
- `Placement`: `reach` (a `Reach`), `placed_on` (`"number"`, `"prefix"`,
  `"name"` or `"any"`), `distance_m` (station to line, `d`), `position` (the
  mapped position `P`, (x, y)), `elvid`, `objectid` (of the chosen segment),
  `lake`, `confluence_near`, `reach_up_m`, `reach_down_m`, `reach_fork`.
- `Reach(line=, at=, uncertainty=, corridor=30.0)`, frozen: `line` in
  downstream order, `at` the metres along it to `P`, `uncertainty` = `U`.

The prefix tier is the design's ("Placing the gauge", step 2, pinned
2026-10-05): the segment's number, not equal to the station's, is a leading
part of it as a string and the two share the main number before the first
dot; it comes before the name tier. Its own examples are pinned: `"002.A"` is
a prefix of `"002.AB"`, `"002"` is not, and `"002.AB"` is not one of
`"002.A"`.

Every distance below is in metres on a line a few kilometres long, at UTM
magnitudes (x 5e5, y 6.6e6); the absolute tolerances (1e-6 m) are about 1e3
times float64's spacing at 6.6e6 m (9.3e-10 m), checked at that scale only.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from pydantic import ValidationError
from shapely.geometry import LineString, Point

from nve_fixtures import collection, river, write
from tin_engine.hydrography import RiverSegment
from tin_engine.io.rivers import read_segments

X, Y = 500_000.0, 6_600_000.0
TOL = 1e-6


@pytest.fixture(scope="module")
def gauge() -> ModuleType:
    import tin_engine.gauge as module

    return module


def seg(
    objectid: int,
    coords: Sequence[tuple[float, float]],
    *,
    elvid: str = "1-1-1",
    vassdragsnr: str | None = "002.DCB2",
    name: str | None = "Nea",
    kind: str = "river",
) -> RiverSegment:
    return RiverSegment(
        objectid=objectid,
        elvid=elvid,
        vassdragsnr=vassdragsnr,
        name=name,
        objekttype="InnsjøMidtlinje" if kind == "lake" else "ElvBekk",
        kind=kind,
        line=tuple((float(x), float(y)) for x, y in coords),
    )


def south(x: float, y0: float, y1: float) -> tuple[tuple[float, float], ...]:
    """A straight segment at `x` from northing y0 down to y1 (digitised downstream)."""
    return ((x, y0), (x, y1))


def river_down(
    x: float = X,
    top: float = Y + 3000.0,
    bottom: float = Y - 3000.0,
    step: float = 1000.0,
    **kw: Any,
) -> list[RiverSegment]:
    """A straight river flowing south at `x`, cut into `step`-metre segments,
    objectids 1, 2, ... from the top; joined end to start exactly."""
    out, y, oid = [], top, 1
    while y > bottom + 1e-9:
        out.append(seg(oid, south(x, y, max(y - step, bottom)), **kw))
        y, oid = y - step, oid + 1
    return out


def at_station(gauge: ModuleType, x: float, y: float, **kw: Any) -> Any:
    return gauge.Gauge(x=x, y=y, **kw)


def length(reach: Any) -> float:
    return float(LineString(reach.line).length)


def assert_consistent(placement: Any) -> None:
    """What every placement must satisfy: `at` is `P` along the line, the
    metres up and down add up to the line, and `U` follows `d`."""
    reach = placement.reach
    line = LineString(reach.line)
    assert len(reach.line) >= 2
    assert reach.at == pytest.approx(placement.reach_up_m, abs=TOL)
    assert length(reach) == pytest.approx(placement.reach_up_m + placement.reach_down_m, abs=TOL)
    p = line.interpolate(reach.at)
    assert (p.x, p.y) == pytest.approx(placement.position, abs=TOL)


# ---------------------------------------------------------------------------
# The models
# ---------------------------------------------------------------------------


def test_gauge_and_reach_are_frozen(gauge: ModuleType) -> None:
    g = gauge.Gauge(x=X, y=Y)
    assert g.watercourse is None and g.river is None
    with pytest.raises((ValidationError, AttributeError, TypeError)):
        g.x = 0.0
    reach = gauge.Reach(line=((X, Y), (X, Y - 100.0)), at=50.0, uncertainty=30.0)
    assert reach.corridor == 30.0
    with pytest.raises(ValidationError):
        reach.at = 0.0


def test_the_defaults_are_500_m_and_1000_m(gauge: ModuleType) -> None:
    import inspect

    params = inspect.signature(gauge.place).parameters
    assert params["map_radius"].default == 500.0
    assert params["reach_up"].default == 1000.0
    assert params["map_radius"].kind is inspect.Parameter.KEYWORD_ONLY
    assert params["reach_up"].kind is inspect.Parameter.KEYWORD_ONLY


# ---------------------------------------------------------------------------
# Step 1: candidates within the map radius
# ---------------------------------------------------------------------------


def test_no_line_within_the_radius_gives_none(gauge: ModuleType) -> None:
    segments = river_down(x=X + 501.0)
    assert gauge.place(at_station(gauge, X, Y), segments) is None


def test_a_line_inside_the_radius_is_found(gauge: ModuleType) -> None:
    segments = river_down(x=X + 499.0)
    placement = gauge.place(at_station(gauge, X, Y), segments)
    assert placement is not None
    assert placement.distance_m == pytest.approx(499.0, abs=TOL)


def test_the_radius_is_the_argument(gauge: ModuleType) -> None:
    segments = river_down(x=X + 120.0)
    assert gauge.place(at_station(gauge, X, Y), segments, map_radius=100.0) is None
    assert gauge.place(at_station(gauge, X, Y), segments, map_radius=150.0) is not None


# ---------------------------------------------------------------------------
# Step 2: which river, by tier; the nearest of the best tier
# ---------------------------------------------------------------------------


def two_rivers(near: dict[str, Any], far: dict[str, Any]) -> list[RiverSegment]:
    """A river 20 m west of the station (objectid 1) and one 80 m east (2)."""
    return [
        seg(1, south(X - 20.0, Y + 1500.0, Y - 1500.0), elvid="near", **near),
        seg(2, south(X + 80.0, Y + 1500.0, Y - 1500.0), elvid="far", **far),
    ]


def test_its_own_watercourse_number_beats_a_nearer_line(gauge: ModuleType) -> None:
    segments = two_rivers(
        near={"vassdragsnr": "123.A", "name": "Gaula"},
        far={"vassdragsnr": "002.DCB2", "name": "Gaula"},
    )
    g = at_station(gauge, X, Y, watercourse="002.DCB2", river="Nea")
    placement = gauge.place(g, segments)
    assert placement.objectid == 2
    assert placement.elvid == "far"
    assert placement.placed_on == "number"
    assert placement.distance_m == pytest.approx(80.0, abs=TOL)


def test_the_prefix_tier_beats_a_nearer_line_of_another_number(gauge: ModuleType) -> None:
    segments = two_rivers(
        near={"vassdragsnr": "123.A", "name": "Gaula"},
        far={"vassdragsnr": "002.DC", "name": "Gaula"},
    )
    g = at_station(gauge, X, Y, watercourse="002.DCB2", river="Nea")
    placement = gauge.place(g, segments)
    assert placement.objectid == 2
    assert placement.placed_on == "prefix"


@pytest.mark.parametrize(
    ("segment", "station", "tier"),
    [
        ("002.A", "002.AB", "prefix"),  # the design's own example
        ("002", "002.AB", "any"),  # the design's example: a bare main number is not
        ("002.AB", "002.A", "any"),  # the segment's area drains into the station's
        ("00", "002.AB", "any"),  # a leading string, but not the main number
        ("002.B", "002.AB", "any"),  # same main number, not a leading part
    ],
)
def test_what_is_a_prefix(gauge: ModuleType, segment: str, station: str, tier: str) -> None:
    """The far line (80 m) has `segment`, the near one (20 m) another main
    number; neither shares the station's river name. A prefix places on the
    far line, anything else falls through to the nearest, as `any`."""
    segments = two_rivers(
        near={"vassdragsnr": "123.A", "name": "Gaula"},
        far={"vassdragsnr": segment, "name": "Gaula"},
    )
    placement = gauge.place(at_station(gauge, X, Y, watercourse=station, river="Nea"), segments)
    assert placement.placed_on == tier
    assert placement.objectid == (2 if tier == "prefix" else 1)


def test_the_prefix_tier_comes_before_the_name_tier(gauge: ModuleType) -> None:
    """The near line has the station's river name, the far line a prefix of
    its number: the prefix wins."""
    segments = two_rivers(
        near={"vassdragsnr": "123.A", "name": "Nea"},
        far={"vassdragsnr": "002.A", "name": "Gaula"},
    )
    g = at_station(gauge, X, Y, watercourse="002.AB", river="Nea")
    placement = gauge.place(g, segments)
    assert placement.objectid == 2
    assert placement.placed_on == "prefix"


def test_the_name_tier_beats_a_nearer_line_of_another_name(gauge: ModuleType) -> None:
    segments = two_rivers(
        near={"vassdragsnr": "123.A", "name": "Gaula"},
        far={"vassdragsnr": "456.B", "name": "Nea"},
    )
    g = at_station(gauge, X, Y, watercourse="002.DCB2", river="Nea")
    placement = gauge.place(g, segments)
    assert placement.objectid == 2
    assert placement.placed_on == "name"


def test_no_tier_matches_the_nearest_line_wins_as_any(gauge: ModuleType) -> None:
    segments = two_rivers(
        near={"vassdragsnr": "123.A", "name": "Gaula"},
        far={"vassdragsnr": "456.B", "name": "Orkla"},
    )
    g = at_station(gauge, X, Y, watercourse="002.DCB2", river="Nea")
    placement = gauge.place(g, segments)
    assert placement.objectid == 1
    assert placement.placed_on == "any"


def test_a_gauge_without_number_or_name_takes_the_nearest_line(gauge: ModuleType) -> None:
    """The command `catchment --rivers` places with no watercourse number."""
    segments = two_rivers(
        near={"vassdragsnr": "123.A", "name": "Gaula"},
        far={"vassdragsnr": "002.DCB2", "name": "Nea"},
    )
    placement = gauge.place(at_station(gauge, X, Y), segments)
    assert placement.objectid == 1
    assert placement.placed_on == "any"


def test_the_nearest_of_the_best_tier_wins(gauge: ModuleType) -> None:
    segments = [
        seg(1, south(X - 10.0, Y + 1500.0, Y - 1500.0), elvid="a", vassdragsnr="123.A"),
        seg(2, south(X + 70.0, Y + 1500.0, Y - 1500.0), elvid="b", vassdragsnr="002.DCB2"),
        seg(3, south(X - 40.0, Y + 1500.0, Y - 1500.0), elvid="c", vassdragsnr="002.DCB2"),
    ]
    g = at_station(gauge, X, Y, watercourse="002.DCB2")
    placement = gauge.place(g, segments)
    assert placement.objectid == 3
    assert placement.placed_on == "number"
    assert placement.distance_m == pytest.approx(40.0, abs=TOL)


# ---------------------------------------------------------------------------
# Step 3: the mapped position P, d and U
# ---------------------------------------------------------------------------


def test_p_is_the_foot_of_the_perpendicular(gauge: ModuleType) -> None:
    placement = gauge.place(at_station(gauge, X + 20.0, Y + 250.0), river_down())
    assert placement.position == pytest.approx((X, Y + 250.0), abs=TOL)
    assert placement.distance_m == pytest.approx(20.0, abs=TOL)
    assert_consistent(placement)


def test_p_is_an_end_vertex_past_the_line(gauge: ModuleType) -> None:
    """South of the river's last vertex: P is that vertex, and no river is
    left downstream of it."""
    segments = river_down(bottom=Y)
    placement = gauge.place(at_station(gauge, X + 30.0, Y - 40.0), segments)
    assert placement.position == pytest.approx((X, Y), abs=TOL)
    assert placement.distance_m == pytest.approx(50.0, abs=TOL)
    assert placement.reach_down_m == pytest.approx(0.0, abs=TOL)
    assert_consistent(placement)


@pytest.mark.parametrize(
    ("d", "radius", "expected"),
    [
        (10.0, 500.0, 30.0),  # the 30 m floor
        (100.0, 500.0, 100.0),  # d itself
        (10.0, 20.0, 20.0),  # the cap: the map radius
    ],
)
def test_u_is_d_floored_at_30_m_and_capped_at_the_radius(
    gauge: ModuleType, d: float, radius: float, expected: float
) -> None:
    placement = gauge.place(at_station(gauge, X + d, Y), river_down(), map_radius=radius)
    assert placement.reach.uncertainty == pytest.approx(expected, abs=TOL)


# ---------------------------------------------------------------------------
# Step 4: the reach
# ---------------------------------------------------------------------------


def test_the_reach_runs_reach_up_above_p_and_u_plus_100_below(gauge: ModuleType) -> None:
    placement = gauge.place(at_station(gauge, X + 10.0, Y), river_down())
    assert placement.reach_up_m == pytest.approx(1000.0, abs=TOL)
    assert placement.reach_down_m == pytest.approx(30.0 + 100.0, abs=TOL)
    reach = placement.reach
    assert reach.line[0] == pytest.approx((X, Y + 1000.0), abs=TOL)
    assert reach.line[-1] == pytest.approx((X, Y - 130.0), abs=TOL)
    assert_consistent(placement)


@pytest.mark.parametrize("reach_up", [250.0, 1000.0, 1234.5])
def test_reach_up_cuts_at_the_metre(gauge: ModuleType, reach_up: float) -> None:
    placement = gauge.place(at_station(gauge, X + 10.0, Y), river_down(), reach_up=reach_up)
    assert placement.reach_up_m == pytest.approx(reach_up, abs=TOL)
    assert placement.reach.at == pytest.approx(reach_up, abs=TOL)
    assert placement.reach.line[0] == pytest.approx((X, Y + reach_up), abs=TOL)


def test_a_river_that_starts_sooner_gives_a_shorter_reach(gauge: ModuleType) -> None:
    placement = gauge.place(at_station(gauge, X + 10.0, Y), river_down(top=Y + 300.0))
    assert placement.reach_up_m == pytest.approx(300.0, abs=TOL)
    assert placement.reach.at == pytest.approx(300.0, abs=TOL)
    assert not placement.reach_fork
    assert_consistent(placement)


def test_a_river_that_ends_sooner_gives_a_shorter_reach(gauge: ModuleType) -> None:
    placement = gauge.place(at_station(gauge, X + 10.0, Y), river_down(bottom=Y - 70.0))
    assert placement.reach_down_m == pytest.approx(70.0, abs=TOL)
    assert not placement.reach_fork
    assert_consistent(placement)


def test_the_chain_joins_ends_within_1_m(gauge: ModuleType) -> None:
    """The next segment starts 0.8 m from the chosen one's end: one chain."""
    segments = [
        seg(1, south(X, Y + 1500.0, Y - 50.0)),
        seg(2, ((X + 0.8, Y - 50.0), (X + 0.8, Y - 1500.0))),
    ]
    placement = gauge.place(at_station(gauge, X + 10.0, Y), segments)
    assert placement.reach_down_m == pytest.approx(130.0, abs=1.0)
    assert placement.reach_down_m > 50.0 + 1.0


def test_the_chain_stops_at_a_gap_of_2_m(gauge: ModuleType) -> None:
    segments = [
        seg(1, south(X, Y + 1500.0, Y - 50.0)),
        seg(2, south(X, Y - 52.0, Y - 1500.0)),
    ]
    placement = gauge.place(at_station(gauge, X + 10.0, Y), segments)
    assert placement.reach_down_m == pytest.approx(50.0, abs=TOL)
    assert not placement.reach_fork


def test_the_chain_follows_its_own_elvid_only(gauge: ModuleType) -> None:
    """A segment of another river starting at the chosen one's end does not
    continue the chain."""
    segments = [
        seg(1, south(X, Y + 1500.0, Y - 50.0), elvid="main"),
        seg(2, south(X, Y - 50.0, Y - 1500.0), elvid="other"),
    ]
    placement = gauge.place(at_station(gauge, X + 10.0, Y), segments)
    assert placement.reach_down_m == pytest.approx(50.0, abs=TOL)


def test_a_fork_downstream_stops_the_chain_and_is_flagged(gauge: ModuleType) -> None:
    """Two different segments of the elvid continue it, as at `19.80.0`: the
    chain stops at the fork and reports the metres it has."""
    q = (X, Y - 60.0)
    segments = [
        seg(1, (( X, Y + 1500.0), q)),
        seg(2, (q, (X - 300.0, Y - 1000.0))),
        seg(3, (q, (X + 300.0, Y - 1000.0))),
    ]  # fmt: skip
    placement = gauge.place(at_station(gauge, X + 10.0, Y), segments)
    assert placement.reach_fork
    assert placement.reach_down_m == pytest.approx(60.0, abs=TOL)
    assert placement.reach_up_m == pytest.approx(1000.0, abs=TOL)
    assert_consistent(placement)


def test_a_fork_upstream_stops_the_chain_and_is_flagged(gauge: ModuleType) -> None:
    q = (X, Y + 400.0)
    segments = [
        seg(1, ((X - 300.0, Y + 1500.0), q)),
        seg(2, ((X + 300.0, Y + 1500.0), q)),
        seg(3, (q, (X, Y - 1500.0))),
    ]
    placement = gauge.place(at_station(gauge, X + 10.0, Y), segments)
    assert placement.reach_fork
    assert placement.reach_up_m == pytest.approx(400.0, abs=TOL)
    assert placement.reach_down_m == pytest.approx(130.0, abs=TOL)
    assert_consistent(placement)


def test_exact_copies_read_from_a_file_do_not_make_a_fork(
    gauge: ModuleType, tmp_path: Path
) -> None:
    """Two features with one geometry under one elvid (as at `82.4.0`): the
    reader keeps one, and the chain runs on through it."""
    features = [
        river(11, [(X, Y + 1500.0), (X, Y - 60.0)], elvid="e"),
        river(12, [(X, Y - 60.0), (X, Y - 1500.0)], elvid="e"),
        river(13, [(X, Y - 60.0), (X, Y - 1500.0)], elvid="e"),
    ]
    segments, _, dropped = read_segments(write(tmp_path / "r.geojson", collection(features)))
    assert dropped == 1
    placement = gauge.place(at_station(gauge, X + 10.0, Y), segments)
    assert not placement.reach_fork
    assert placement.reach_down_m == pytest.approx(130.0, abs=TOL)


def test_the_reach_is_in_the_segments_coordinates_and_downstream_order(
    gauge: ModuleType,
) -> None:
    """A bent river: the reach follows its vertices, upstream first."""
    segments = [
        seg(1, ((X - 2000.0, Y + 500.0), (X - 500.0, Y + 500.0), (X, Y))),
        seg(2, ((X, Y), (X, Y - 2000.0))),
    ]
    placement = gauge.place(at_station(gauge, X + 15.0, Y - 300.0), segments)
    reach = placement.reach
    assert placement.position == pytest.approx((X, Y - 300.0), abs=TOL)
    assert placement.reach_up_m == pytest.approx(1000.0, abs=TOL)
    first = Point(reach.line[0])
    # 300 m of the reach lie on the second segment, so 700 m lie on the first
    # one's diagonal (707.1 m long): the first vertex is on that diagonal.
    up_diag = 1000.0 - 300.0
    expected = (X - up_diag / math.sqrt(2.0), Y + up_diag / math.sqrt(2.0))
    assert (first.x, first.y) == pytest.approx(expected, abs=1e-6)
    ys = [y for _, y in reach.line]
    assert ys[0] > ys[-1]
    assert_consistent(placement)


# ---------------------------------------------------------------------------
# Step 5: the flags
# ---------------------------------------------------------------------------


def test_a_lake_centreline_flags_lake(gauge: ModuleType) -> None:
    placement = gauge.place(at_station(gauge, X + 10.0, Y), river_down(kind="lake"))
    assert placement.lake


def test_a_river_line_does_not_flag_lake(gauge: ModuleType) -> None:
    placement = gauge.place(at_station(gauge, X + 10.0, Y), river_down())
    assert not placement.lake


def test_kind_read_from_the_file_reaches_the_lake_flag(gauge: ModuleType, tmp_path: Path) -> None:
    """A null `objekttype` with a lake number is a lake (PR 3's mapping)."""
    features = [
        river(21, [(X, Y + 1500.0), (X, Y - 1500.0)], objekttype=None, vatnlnr=495),
    ]
    segments, _, _ = read_segments(write(tmp_path / "r.geojson", collection(features)))
    placement = gauge.place(at_station(gauge, X + 10.0, Y), segments)
    assert placement.lake


def branch(x: float) -> RiverSegment:
    """A short segment of another river, centred at (x, Y)."""
    return seg(99, ((x, Y + 5.0), (x, Y - 5.0)), elvid="tributary")


@pytest.mark.parametrize(
    ("branch_x", "expected"),
    [
        (X - 90.0, True),  # 90 m from P, 110 m from the station: counts
        (X + 110.0, False),  # 110 m from P, 90 m from the station: does not
        (X - 99.0, True),
        (X - 101.0, False),
    ],
)
def test_confluence_near_is_measured_from_p_at_100_m(
    gauge: ModuleType, branch_x: float, expected: bool
) -> None:
    placement = gauge.place(at_station(gauge, X + 20.0, Y), [*river_down(), branch(branch_x)])
    assert placement.objectid != 99
    assert placement.position == pytest.approx((X, Y), abs=TOL)
    assert placement.confluence_near is expected


def test_its_own_river_is_not_a_confluence(gauge: ModuleType) -> None:
    placement = gauge.place(at_station(gauge, X + 20.0, Y), river_down())
    assert placement.confluence_near is False


# ---------------------------------------------------------------------------
# PR 4, lake gauges: `lake_seed` ("Lake gauges: the lake is the seed")
# ---------------------------------------------------------------------------
#
# Interface pinned here, from the design's own signatures: `gauge.LAKE_GAP_M`
# (30.0), `gauge.LakeSeed` (frozen; `rule`, `point`, `lakes`, `distance_m`)
# and `gauge.lake_seed(gauge, placement, lakes, *, gap=LAKE_GAP_M)`; `Lake`
# is `io/station_set.py`'s frozen dataclass (`number`, `name`, `polygon`).
# Placements come from `place` on a straight line down x = X, so `P` is
# (X, Y) for a station at (X + d, Y); lakes are axis-aligned boxes, so every
# station-to-lake distance below is a difference of two whole or half
# metres at x 5e5 and exact in float64. Before the change, `lake_seed`,
# `LakeSeed`, `LAKE_GAP_M` and `Lake` do not exist.


def lake(
    x0: float, y0: float, x1: float, y1: float, number: int | None = 7, name: str | None = "Vatnet"
) -> Any:
    from shapely.geometry import box

    from tin_engine.hydrography import Lake

    return Lake(number=number, name=name, polygon=box(x0, y0, x1, y1))


def on_lake_line(gauge: ModuleType, station_x: float) -> Any:
    """The placement of a station at (station_x, Y) on a lake centreline down x = X."""
    placement = gauge.place(at_station(gauge, station_x, Y), river_down(kind="lake"))
    assert placement.lake and placement.position == pytest.approx((X, Y), abs=TOL)
    return placement


def on_river_line(gauge: ModuleType, station_x: float) -> Any:
    placement = gauge.place(at_station(gauge, station_x, Y), river_down())
    assert not placement.lake and placement.position == pytest.approx((X, Y), abs=TOL)
    return placement


def seed_of(
    gauge: ModuleType, station_x: float, placement: Any, lakes: Sequence[Any], **kw: Any
) -> Any:
    return gauge.lake_seed(at_station(gauge, station_x, Y), placement, lakes, **kw)


def test_the_gap_is_30_m(gauge: ModuleType) -> None:
    assert gauge.LAKE_GAP_M == 30.0


def test_a_station_inside_a_lake_on_a_river_line_is_inside(gauge: ModuleType) -> None:
    """`62.18.0` Svartavatn's case: placed on a river line, inside a lake."""
    vatn = lake(X + 5.0, Y - 100.0, X + 200.0, Y + 100.0)
    seed = seed_of(gauge, X + 20.0, on_river_line(gauge, X + 20.0), [vatn])
    assert isinstance(seed, gauge.LakeSeed)
    assert seed.rule == "inside"
    assert seed.point == pytest.approx((X + 20.0, Y), abs=TOL)
    assert seed.lakes == (vatn,)
    assert seed.distance_m == 0.0


def test_a_station_inside_a_lake_with_no_placement_is_inside(gauge: ModuleType) -> None:
    vatn = lake(X + 5.0, Y - 100.0, X + 200.0, Y + 100.0)
    seed = seed_of(gauge, X + 20.0, None, [vatn])
    assert seed is not None and seed.rule == "inside"
    assert seed.point == pytest.approx((X + 20.0, Y), abs=TOL)
    assert seed.lakes == (vatn,)


def test_a_lake_line_with_p_in_a_lake_20_m_away_is_lake_line(gauge: ModuleType) -> None:
    """The distance is the station's to the lake (20 m), not to the line (25 m)."""
    vatn = lake(X - 500.0, Y - 500.0, X + 5.0, Y + 500.0)
    seed = seed_of(gauge, X + 25.0, on_lake_line(gauge, X + 25.0), [vatn])
    assert seed is not None and seed.rule == "lake_line"
    assert seed.point == pytest.approx((X, Y), abs=TOL)  # P, not the station
    assert seed.lakes == (vatn,)
    assert seed.distance_m == pytest.approx(20.0, abs=TOL)


@pytest.mark.parametrize(
    ("station_x", "expected"),
    [(X + 35.0, "lake_line"), (X + 35.5, None)],
    ids=["exactly-30-m", "30.5-m"],
)
def test_the_gap_boundary(gauge: ModuleType, station_x: float, expected: str | None) -> None:
    vatn = lake(X - 500.0, Y - 500.0, X + 5.0, Y + 500.0)
    seed = seed_of(gauge, station_x, on_lake_line(gauge, station_x), [vatn])
    assert (None if seed is None else seed.rule) == expected


def test_the_gap_is_the_argument(gauge: ModuleType) -> None:
    vatn = lake(X - 500.0, Y - 500.0, X + 5.0, Y + 500.0)
    placement = on_lake_line(gauge, X + 35.5)
    seed = seed_of(gauge, X + 35.5, placement, [vatn], gap=31.0)
    assert seed is not None and seed.distance_m == pytest.approx(30.5, abs=TOL)
    assert seed_of(gauge, X + 25.0, on_lake_line(gauge, X + 25.0), [vatn], gap=10.0) is None


def test_a_river_line_10_m_below_a_lake_is_not_seeded(gauge: ModuleType) -> None:
    """A gauge below the outlet stays on the river path (question 10's default),
    though the lake is within 30 m of the station."""
    vatn = lake(X - 200.0, Y + 10.0, X + 200.0, Y + 500.0)
    assert seed_of(gauge, X + 5.0, on_river_line(gauge, X + 5.0), [vatn]) is None


def test_a_lake_line_whose_p_is_in_no_lake_is_not_seeded(gauge: ModuleType) -> None:
    """The station is 10 m from a lake, but `P` (30 m from it) is in none."""
    vatn = lake(X + 30.0, Y - 100.0, X + 300.0, Y + 100.0)
    assert seed_of(gauge, X + 20.0, on_lake_line(gauge, X + 20.0), [vatn]) is None


def test_no_placement_and_a_lake_14_m_away_is_not_seeded(gauge: ModuleType) -> None:
    """Femundsenden (`311.4.0`), question 9's default: no river line, so no
    rule says which way the water runs past the station."""
    femunden = lake(X - 5000.0, Y - 5000.0, X + 5.0, Y + 5000.0)
    assert seed_of(gauge, X + 19.0, None, [femunden]) is None


def test_two_overlapping_lakes_both_reach_the_seed(gauge: ModuleType) -> None:
    """A user's file with overlapping polygons: the seed holds both, and
    22's `_lake` refuses it later (not `lake_seed`'s job)."""
    a = lake(X + 5.0, Y - 100.0, X + 200.0, Y + 100.0, number=1, name="A")
    b = lake(X + 10.0, Y - 50.0, X + 100.0, Y + 50.0, number=2, name="B")
    c = lake(X + 300.0, Y - 50.0, X + 400.0, Y + 50.0, number=3, name="C")
    seed = seed_of(gauge, X + 20.0, None, [a, b, c])
    assert seed is not None and seed.rule == "inside"
    assert {lk.number for lk in seed.lakes} == {1, 2}


def test_no_rule_reads_an_area_the_smaller_lake_containing_the_station_wins(
    gauge: ModuleType,
) -> None:
    small = lake(X + 10.0, Y - 10.0, X + 30.0, Y + 10.0, number=1, name="Tjørna")
    large = lake(X + 35.0, Y - 50_000.0, X + 50_000.0, Y + 50_000.0, number=2, name="Storvatnet")
    seed = seed_of(gauge, X + 20.0, on_river_line(gauge, X + 20.0), [large, small])
    assert seed is not None and seed.lakes == (small,)


def test_the_lake_line_lake_is_ps_not_the_nearest_or_the_largest(gauge: ModuleType) -> None:
    """`P` lies in a small lake 15 m from the station; a large lake 5 m from
    the station holds neither the station nor `P`, and is not chosen."""
    small = lake(X - 10.0, Y - 10.0, X + 10.0, Y + 10.0, number=1, name="Tjørna")
    large = lake(X + 30.0, Y - 50_000.0, X + 50_000.0, Y + 50_000.0, number=2, name="Storvatnet")
    seed = seed_of(gauge, X + 25.0, on_lake_line(gauge, X + 25.0), [large, small])
    assert seed is not None and seed.rule == "lake_line"
    assert seed.lakes == (small,)
    assert seed.distance_m == pytest.approx(15.0, abs=TOL)


def test_a_lake_without_a_number_is_found_by_geometry(gauge: ModuleType) -> None:
    """`97.1.0` Fetvatn: a lake line with no lake number still finds its lake."""
    vatn = lake(X - 500.0, Y - 500.0, X + 5.0, Y + 500.0, number=None, name=None)
    seed = seed_of(gauge, X + 25.0, on_lake_line(gauge, X + 25.0), [vatn])
    assert seed is not None and seed.lakes[0].number is None


def test_inside_comes_before_lake_line(gauge: ModuleType) -> None:
    """The station inside lake B, `P` inside lake A: `inside`, with B."""
    a = lake(X - 10.0, Y - 10.0, X + 10.0, Y + 10.0, number=1, name="A")
    b = lake(X + 20.0, Y - 10.0, X + 40.0, Y + 10.0, number=2, name="B")
    seed = seed_of(gauge, X + 25.0, on_lake_line(gauge, X + 25.0), [a, b])
    assert seed is not None and seed.rule == "inside"
    assert seed.lakes == (b,)
    assert seed.point == pytest.approx((X + 25.0, Y), abs=TOL)


def test_a_station_on_the_boundary_is_not_inside(gauge: ModuleType) -> None:
    """`Polygon.contains` is strict, as 22's `_lake`: a station exactly on
    the shore is not `inside`; on a lake line with `P` in that lake it is
    `lake_line` at distance 0."""
    vatn = lake(X - 500.0, Y - 500.0, X + 20.0, Y + 500.0)
    assert seed_of(gauge, X + 20.0, None, [vatn]) is None
    seed = seed_of(gauge, X + 20.0, on_lake_line(gauge, X + 20.0), [vatn])
    assert seed is not None and seed.rule == "lake_line"
    assert seed.distance_m == pytest.approx(0.0, abs=TOL)


@pytest.mark.parametrize(("offset", "inside"), [(-0.001, True), (0.001, False)])
def test_a_millimetre_at_utm_magnitudes(gauge: ModuleType, offset: float, inside: bool) -> None:
    """1 mm inside or outside a shore at x 5e5 (float64's spacing there is
    about 6e-11 m), with no placement."""
    vatn = lake(X - 500.0, Y - 500.0, X + 20.0, Y + 500.0)
    seed = seed_of(gauge, X + 20.0 + offset, None, [vatn])
    assert (seed is not None) is inside


@pytest.mark.parametrize("bad", [math.nan, math.inf], ids=["nan", "inf"])
def test_a_non_finite_station_is_not_seeded(gauge: ModuleType, bad: float) -> None:
    vatn = lake(X - 500.0, Y - 500.0, X + 20.0, Y + 500.0)
    assert gauge.lake_seed(gauge.Gauge(x=bad, y=Y), None, [vatn]) is None


def test_no_lakes_is_no_seed(gauge: ModuleType) -> None:
    assert seed_of(gauge, X + 20.0, on_lake_line(gauge, X + 20.0), []) is None


def test_the_seed_and_the_lake_are_frozen(gauge: ModuleType) -> None:
    import dataclasses

    vatn = lake(X + 5.0, Y - 100.0, X + 200.0, Y + 100.0)
    seed = seed_of(gauge, X + 20.0, None, [vatn])
    with pytest.raises(ValidationError):
        seed.rule = "lake_line"
    with pytest.raises(dataclasses.FrozenInstanceError):
        vatn.number = 8
