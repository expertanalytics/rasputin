"""`burn.burn_reach`: the mapped reach moved to the valley floor and burnt (29, PR 2).

`docs/increments/29-nve-reference-catchments.md`, "Following the river in the
DEM" (steps 1 to 6 and the paragraph after them) and "The red suites", PR 2's
`test_burn.py`. Hand-built terrains on DTM10's 10 m lattice
(`gauge_fixtures`); every flood is PR 1's `accumulate` or 22's `upstream`.

Interface pinned here (the design names `burn_reach(window, reach) -> burnt
window, GaugePath` and the path's contents; the names it leaves open are
chosen here and listed in the handback):

- `burn_reach(window: DemTile, reach: Reach) -> (DemTile, GaugePath)`: the
  burnt copy has the window's `meta`.
- `GaugePath`: `chain` ((n, 2) (row, column), downstream, extension
  included), `placed` (one Python `int`, an index into `chain`), `arc` (metres along the
  chain from the placed node, negative upstream), `end_extended_m`,
  `end_closed`, `direction_ok`, `lowered_nodes`, `lowered_max_m`,
  `node_offset_m` (P to the placed node), `dropped_m` (metres of the mapped
  reach outside the window, dropped).

THE ORACLE for every burn below is `gauge_fixtures.hand_burn` along the chain
the burn returns: the burnt array must equal it bit for bit. That one check is
the design's "`<=` the input everywhere, equal off the chain, and on the chain
equal to the input wherever the input is already below the carried level minus
0.001 m, computed in the array's dtype", extension included.

Most reaches here use `corridor=5.0` with resample points on nodes, so each
point's corridor holds exactly the node under it and the chain is known
exactly; the valley-floor tests use the default 30 m. The valley floor is
taken across the line ("Following the river", step 2): a point's
cross-section is the nodes within the corridor whose offset along the line is
at most half the 10 m step, so on these straight lines down a column it is the
point's own row, and the chosen node is abeam the point, not downstream.

Elevations are a few hundred metres (3000 m in the float32 test); the
tolerances on metres along the chain (1e-9 m) are about 1e4 times float64's
spacing at 500 m, checked only at these sizes.
"""

from __future__ import annotations

import itertools
import math
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest

from gauge_fixtures import (
    CC,
    CELL_KM2,
    DAM_CREST,
    DAM_ROWS,
    channel,
    column_line,
    drains_along,
    floor_z,
    hand_burn,
    is_taut,
    lat,
    tile_of,
    valley,
)
from tin_engine._core import accumulate, upstream
from tin_engine.io.models import DemTile
from tin_engine.raster import to_core

NODATA = -32767.0


@pytest.fixture(scope="module")
def burn() -> ModuleType:
    import tin_engine.burn as module

    return module


@pytest.fixture(scope="module")
def gauge() -> ModuleType:
    import tin_engine.gauge as module

    return module


def reach(gauge: ModuleType, line: Any, at: float, u: float = 30.0, corridor: float = 5.0) -> Any:
    return gauge.Reach(line=tuple(line), at=at, uncertainty=u, corridor=corridor)


def nodes(path: Any) -> list[tuple[int, int]]:
    return [(int(r), int(c)) for r, c in np.asarray(path.chain).reshape(-1, 2)]


def run(
    burn: ModuleType, z: npt.NDArray[np.floating], r: Any, nodata: float | None = None
) -> tuple[npt.NDArray[np.floating], Any]:
    """Burn, and check what every burn must satisfy."""
    tile = tile_of(z, nodata)
    before = np.array(tile.array, copy=True)
    out, path = burn.burn_reach(tile, r)
    assert isinstance(out, DemTile)
    assert out.meta == tile.meta
    burnt = np.asarray(out.array)
    raw = np.asarray(tile.array)
    assert np.array_equal(raw, before)  # the window is not written
    assert not np.shares_memory(burnt, raw)
    assert burnt.dtype == raw.dtype and burnt.shape == raw.shape
    chain = nodes(path)
    assert is_taut(chain), chain
    assert np.array_equal(burnt, hand_burn(raw, chain))
    assert np.all(burnt <= raw)
    values = [burnt[n] for n in chain]
    assert all(b < a for a, b in itertools.pairwise(values)), "the chain falls strictly"
    if nodata is not None:
        assert all(raw[n] != nodata for n in chain)
    placed = int(path.placed)
    arc = np.asarray(path.arc, dtype=np.float64)
    assert arc.shape == (len(chain),)
    assert isinstance(path.placed, int) and not isinstance(path.placed, bool)
    assert 0 <= placed < len(chain)
    assert arc[placed] == 0.0
    assert np.all(np.diff(arc) > 0.0)
    return burnt, path


def flow_to(z: npt.NDArray[np.floating], nodata: float | None = None) -> npt.NDArray[np.uint8]:
    return np.asarray(accumulate(to_core(tile_of(z, nodata))).flow_to)


def catchment(z: npt.NDArray[np.floating], node: tuple[int, int]) -> Any:
    seed = np.zeros(z.shape, dtype=np.uint8)
    seed[node] = 1
    return upstream(to_core(tile_of(z)), seed)


def straight(path: Any) -> list[float]:
    """The metres of each step along the chain."""
    chain = np.asarray(path.chain, dtype=np.float64).reshape(-1, 2)
    return list(np.hypot(*np.diff(chain, axis=0).T) * 10.0)


# ---------------------------------------------------------------------------
# The embankment across the valley
# ---------------------------------------------------------------------------


def test_an_embankment_diverts_the_raw_valley_and_the_burn_restores_it(
    burn: ModuleType, gauge: ModuleType
) -> None:
    z = valley(dam=True)
    burnt, path = run(burn, z, reach(gauge, column_line(CC, 40, 200), at=1100.0, corridor=30.0))
    below = (150, CC)
    head = (50, CC)  # the valley floor far above the embankment
    raw_out, burnt_out = catchment(z, below), catchment(burnt, below)
    assert np.asarray(raw_out.mask)[head] == 0
    assert np.asarray(burnt_out.mask)[head] == 1
    assert burnt_out.nodes_in > raw_out.nodes_in
    chain = nodes(path)
    assert all(c == CC for _, c in chain), "the chain is on the valley floor"
    assert {(r, CC) for r in DAM_ROWS} <= set(chain)
    assert path.lowered_nodes == len(DAM_ROWS)
    deepest = DAM_CREST - (floor_z(DAM_ROWS[0] - 1) - 0.001 * len(DAM_ROWS))
    assert path.lowered_max_m == pytest.approx(deepest, abs=1e-6)
    assert path.direction_ok
    assert path.end_closed
    assert path.end_extended_m == 0.0


def test_the_burn_is_only_on_the_chain(burn: ModuleType, gauge: ModuleType) -> None:
    z = valley(dam=True)
    burnt, path = run(burn, z, reach(gauge, column_line(CC, 40, 200), at=1100.0, corridor=30.0))
    off = np.ones(z.shape, dtype=bool)
    for n in nodes(path):
        off[n] = False
    assert np.array_equal(burnt[off], z[off])


# ---------------------------------------------------------------------------
# The descent: strict, in the array's own dtype
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_a_flat_chain_at_3000_m_falls_strictly(
    burn: ModuleType, gauge: ModuleType, dtype: type
) -> None:
    """float32's spacing at 3000 m is 0.00024 m, so 0.001 m steps stay strict."""
    z = channel(120, 41, [(0, 20), (119, 20)], top=3000.0, slope=0.0, dtype=dtype)
    burnt, path = run(burn, z, reach(gauge, column_line(20, 5, 40), at=100.0))
    chain = nodes(path)
    assert chain[:36] == [(r, 20) for r in range(5, 41)]
    assert burnt[chain[0]] == z[chain[0]]
    assert burnt[chain[1]] == (z[chain[0]] - np.asarray(0.001, dtype=dtype)).astype(dtype)


# ---------------------------------------------------------------------------
# The end of the chain (step 5)
# ---------------------------------------------------------------------------

LINE_END = 40  # the mapped line runs down column 20 from row 5 to row 40


def falling(rows: int, cols: int = 41) -> npt.NDArray[np.float64]:
    return channel(rows, cols, [(0, 20), (rows - 1, 20)])


def flat_below(rows: int, cols: int = 41) -> npt.NDArray[np.float64]:
    """Falling to row 38, then a flat floor (row 38's level) from row 39 on:
    the line's last two nodes are lowered, and so is every node below."""
    z = falling(rows, cols)
    z[39:] = z[38]
    return z


def end_reach(gauge: ModuleType) -> Any:
    return reach(gauge, column_line(20, 5, LINE_END), at=100.0)


def test_an_embankment_in_the_last_50_m_extends_the_chain_until_the_ground_falls(
    burn: ModuleType, gauge: ModuleType
) -> None:
    z = falling(110)
    z[38:43] = np.maximum(z[38:43], floor_z(30))
    burnt, path = run(burn, z, end_reach(gauge))
    chain = nodes(path)
    assert chain == [(r, 20) for r in range(5, 44)]
    assert path.end_extended_m == pytest.approx(30.0, abs=1e-9)
    assert path.end_closed
    assert burnt[chain[-1]] == z[chain[-1]], "the last node needed no lowering"
    assert drains_along(flow_to(burnt), chain) == len(chain) - 1


def test_a_flat_floor_below_the_end_hits_the_cap_after_50_straight_steps(
    burn: ModuleType, gauge: ModuleType
) -> None:
    z = flat_below(130)
    _, path = run(burn, z, end_reach(gauge))
    chain = nodes(path)
    assert chain == [(r, 20) for r in range(5, 41 + 50)]
    assert path.end_extended_m == pytest.approx(500.0, abs=1e-9)
    assert path.end_closed is False


def test_a_diagonal_flat_floor_hits_the_cap_after_35_steps(
    burn: ModuleType, gauge: ModuleType
) -> None:
    """35 diagonal steps are 494.97 m; a 36th would make 509.12 m."""
    rows, cols = 100, 80
    z = flat_below(rows, cols)
    r, c = np.indices((rows, cols)).astype(np.float64)
    trough = z[38, 20] + 3.0 * np.abs((r - 40.0) - (c - 20.0)) / math.sqrt(2.0)
    z[41:] = trough[41:]
    _, path = run(burn, z, end_reach(gauge))
    chain = nodes(path)
    assert chain[:36] == [(r, 20) for r in range(5, 41)]
    assert chain[36:] == [(40 + k, 20 + k) for k in range(1, 36)]
    assert path.end_extended_m == pytest.approx(35 * 10.0 * math.sqrt(2.0), abs=1e-9)
    assert path.end_closed is False


def test_the_window_edge_before_the_ground_falls_leaves_the_end_open(
    burn: ModuleType, gauge: ModuleType
) -> None:
    z = flat_below(56)
    _, path = run(burn, z, end_reach(gauge))
    assert path.end_closed is False
    assert 0.0 < path.end_extended_m < 500.0
    assert all(0 <= r < 56 for r, _ in nodes(path))


def test_nodata_before_the_ground_falls_leaves_the_end_open(
    burn: ModuleType, gauge: ModuleType
) -> None:
    z = flat_below(130)
    z[50] = NODATA
    _, path = run(burn, z, end_reach(gauge), nodata=NODATA)
    assert path.end_closed is False
    assert 0.0 < path.end_extended_m < 500.0
    assert max(r for r, _ in nodes(path)) < 50


def test_no_neighbour_to_extend_into_leaves_the_end_open(
    burn: ModuleType, gauge: ModuleType
) -> None:
    """The three nodes below the last one are NoData; the two beside it are
    neighbours of the node before it: no step keeps the chain taut."""
    z = flat_below(130)
    z[LINE_END + 1] = NODATA
    _, path = run(burn, z, end_reach(gauge), nodata=NODATA)
    assert path.end_closed is False
    assert path.end_extended_m == 0.0
    assert nodes(path)[-1] == (LINE_END, 20)


def test_the_extension_does_not_step_into_a_corner(burn: ModuleType, gauge: ModuleType) -> None:
    """The lowest neighbour of the last node, (40, 21), is also a neighbour of
    (39, 20): taking it would make a corner. The next lowest, (41, 21), is
    taken, and its raw ground already falls, so the end is closed."""
    z = flat_below(130)
    level = z[38, 20]
    z[LINE_END, 21] = level - 5.0
    z[LINE_END + 1, 21] = level - 0.5
    _, path = run(burn, z, end_reach(gauge))
    chain = nodes(path)
    assert (LINE_END, 21) not in chain
    assert chain[-1] == (LINE_END + 1, 21)
    assert path.end_closed
    assert path.end_extended_m == pytest.approx(10.0 * math.sqrt(2.0), abs=1e-9)


# ---------------------------------------------------------------------------
# The taut chain (step 2)
# ---------------------------------------------------------------------------


def test_a_cornered_chain_is_made_taut_and_drains(burn: ModuleType, gauge: ModuleType) -> None:
    """The round 3 case: the line turns a square corner at (30, 20). The raw
    valley-floor chain would hold (29, 20), (30, 20), (30, 21), and (29, 20)
    drains past (30, 20); the taut chain drops (30, 20) and drains node by
    node."""
    z = channel(60, 70, [(0, 20), (30, 20), (30, 69)])
    line = (lat(20, 10), lat(20, 30), lat(50, 30))
    burnt, path = run(burn, z, reach(gauge, line, at=100.0))
    taut = [(r, 20) for r in range(10, 30)] + [(30, c) for c in range(21, 51)]
    assert nodes(path) == taut
    assert drains_along(flow_to(burnt), taut) == len(taut) - 1
    # The premise: the chain with its corner does not drain.
    cornered = [(r, 20) for r in range(10, 31)] + [(30, c) for c in range(21, 51)]
    assert drains_along(flow_to(hand_burn(z, cornered)), cornered) < len(cornered) - 1


def test_a_line_that_doubles_back_drops_the_loop(burn: ModuleType, gauge: ModuleType) -> None:
    """Down column 20 from row 5 to row 20, back up to row 15, down to row 40.
    The taut pass cuts at row 14 (the first node with a later neighbour, the
    second visit of row 15): the chain is rows 5 to 40 once. The resample
    point nearest `at` = 150 m is the turn (row 20, dropped), so the placed
    node is row 14, the node the cut starts from."""
    z = channel(60, 41, [(0, 20), (59, 20)])
    line = (lat(20, 5), lat(20, 20), lat(20, 15), lat(20, 40))
    _, path = run(burn, z, reach(gauge, line, at=150.0))
    chain = nodes(path)
    assert chain == [(r, 20) for r in range(5, 41)]
    assert chain[int(path.placed)] == (14, 20)
    assert int(path.placed) == 9
    assert np.asarray(path.arc) == pytest.approx(10.0 * (np.arange(36) - 9), abs=1e-9)
    assert path.node_offset_m == pytest.approx(60.0, abs=1e-9)


# ---------------------------------------------------------------------------
# The valley floor (step 2) and the placed node (step 6)
# ---------------------------------------------------------------------------


def test_a_line_off_the_valley_floor_gives_a_chain_on_it(
    burn: ModuleType, gauge: ModuleType
) -> None:
    """The line is 28 m east of the floor. Each resample point's cross-section
    is its own row, columns 20 to 25 within 30 m, and its least elevation is
    the floor node abeam it, 28 m away. The whole 30 m disc would have chosen
    the floor node one row downstream (29.7 m away, 0.5 m lower), rows 6 to
    41: the downstream bias step 2 removes."""
    z = channel(80, 41, [(0, 20), (79, 20)])
    assert z[6, 20] < z[5, 20]  # the premise: one row down is lower
    _, path = run(burn, z, reach(gauge, column_line(22.8, 5, 40), at=100.0, corridor=30.0))
    assert nodes(path) == [(r, 20) for r in range(5, 41)]


def test_a_corridor_of_two_cells_keeps_the_chain_within_two(
    burn: ModuleType, gauge: ModuleType
) -> None:
    z = channel(80, 41, [(0, 20), (79, 20)])
    _, path = run(burn, z, reach(gauge, column_line(22.8, 5, 40), at=100.0, corridor=20.0))
    chain = nodes(path)
    assert chain == [(r, 21) for r in range(5, 41)]
    assert all(abs(c - 22.8) * 10.0 <= 20.0 for _, c in chain)


@pytest.mark.parametrize(("offset", "column"), [(0.3, 20), (0.7, 21), (0.5, 20)])
def test_ties_go_to_the_nearer_node_then_the_smaller_row_and_column(
    burn: ModuleType, gauge: ModuleType, offset: float, column: int
) -> None:
    """A floor two nodes wide (columns 20 and 21 equal), falling south: each
    point's least elevation is in its own row (its cross-section), on either
    floor column; the nearer wins, and at equal distance (offset 0.5) the
    smaller column."""
    z = channel(80, 42, [(0, 20.5), (79, 20.5)], width=0.5)
    assert z[30, 20] == z[30, 21]  # the premise: a tie in elevation
    _, path = run(burn, z, reach(gauge, column_line(20 + offset, 5, 40), at=100.0, corridor=30.0))
    assert nodes(path) == [(r, column) for r in range(5, 41)]


def test_a_node_exactly_corridor_away_is_inside(burn: ModuleType, gauge: ModuleType) -> None:
    """The line runs halfway between columns 20 and 21, so with a 5 m
    corridor each point's cross-section is (r, 20) and (r, 21), both exactly
    5 m away (to float64's rounding of the resampled northing, well inside
    the design's 1e-6 m slack). The channel is down column 21, so the
    cross-section's least elevation is (r, 21). A strict `<` would leave the
    cross-section empty and fall back to the nearest node, which by the tie
    on distance is the smaller column, (r, 20)."""
    z = channel(60, 42, [(0, 21), (59, 21)])
    _, path = run(burn, z, reach(gauge, column_line(20.5, 5, 40), at=100.0))
    assert nodes(path) == [(r, 21) for r in range(5, 41)]
    assert path.node_offset_m == pytest.approx(5.0, abs=1e-9)


def test_a_cross_section_with_no_node_takes_the_nearest_node(
    burn: ModuleType, gauge: ModuleType
) -> None:
    """The line runs halfway between columns 20 and 21 and 0.3 rows below the
    node rows, so with a 5 m corridor no node is within 5 m of any point
    (the nearest are (r, 20) and (r, 21), sqrt(3^2 + 5^2) = 5.83 m away).
    Each point then takes the node nearest it, whatever its elevation: of
    the two at equal distance (the line is vertical, so its easting is exact
    and the tie is exact), the smaller column, (r, 20), though the channel
    down column 21 is 3 m lower."""
    z = channel(60, 42, [(0, 21), (59, 21)])
    _, path = run(burn, z, reach(gauge, column_line(20.5, 5.3, 40.3), at=100.0))
    assert nodes(path) == [(r, 20) for r in range(5, 41)]
    assert nodes(path)[int(path.placed)] == (15, 20)
    # 1e-6 m: the 3 m is a difference of northings near 6.6e6 m, where
    # float64's spacing is 9.3e-10 m; checked at that northing only.
    assert path.node_offset_m == pytest.approx(math.hypot(3.0, 5.0), abs=1e-6)


def test_the_placed_node_is_the_one_the_point_nearest_at_chose(
    burn: ModuleType, gauge: ModuleType
) -> None:
    """`at` = 100 m is the resample point on row 15; it chose the floor node
    abeam, (15, 20), 28 m across the river from P."""
    z = channel(80, 41, [(0, 20), (79, 20)])
    _, path = run(burn, z, reach(gauge, column_line(22.8, 5, 40), at=100.0, corridor=30.0))
    assert nodes(path)[int(path.placed)] == (15, 20)
    assert int(path.placed) == 10
    assert path.node_offset_m == pytest.approx(28.0, abs=1e-9)
    assert np.asarray(path.arc) == pytest.approx(
        10.0 * (np.arange(36) - int(path.placed)), abs=1e-9
    )
    assert straight(path) == pytest.approx([10.0] * 35, abs=1e-9)


# ---------------------------------------------------------------------------
# The direction check (step 3)
# ---------------------------------------------------------------------------


def rising(rows: int, profile: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """A channel down column 20 whose floor follows `profile` (one value per row)."""
    z = channel(rows, 41, [(0, 20), (rows - 1, 20)], slope=0.0)
    return z + profile[:, None]


def ramp(rows: int, rise: float) -> npt.NDArray[np.float64]:
    """0 to row 10, a ramp to `rise` at row 30, `rise` after: the means of the
    first and last tenths of a chain over rows 0 to 40 are 0 and `rise`,
    whatever the tenth's rounding."""
    r = np.arange(rows, dtype=np.float64)
    return np.clip((r - 10.0) / 20.0, 0.0, 1.0) * rise


@pytest.mark.parametrize(("rise", "ok"), [(2.0, True), (2.5, False), (0.0, True)])
def test_a_line_against_the_slope_by_more_than_2_m_is_flagged(
    burn: ModuleType, gauge: ModuleType, rise: float, ok: bool
) -> None:
    z = rising(100, ramp(100, rise))
    _, path = run(burn, z, reach(gauge, column_line(20, 0, 40), at=200.0))
    assert path.direction_ok is ok


def step_profile(rows: int) -> npt.NDArray[np.float64]:
    """0 on rows 0 to 4, 2.5 m on row 5, 5 m from row 6 on: the first and last
    tenths of a chain over rows 0 to 9 or 0 to 10 differ by 5 m."""
    p = np.zeros(rows)
    p[5] = 2.5
    p[6:] = 5.0
    return p


def test_a_reach_shorter_than_100_m_is_not_checked(burn: ModuleType, gauge: ModuleType) -> None:
    z = rising(60, step_profile(60))
    _, path = run(burn, z, reach(gauge, column_line(20, 0, 9), at=40.0))
    assert path.direction_ok


def test_a_reach_of_100_m_is_checked(burn: ModuleType, gauge: ModuleType) -> None:
    z = rising(60, step_profile(60))
    _, path = run(burn, z, reach(gauge, column_line(20, 0, 10), at=40.0))
    assert path.direction_ok is False


# ---------------------------------------------------------------------------
# The drainage is checked, not assumed
# ---------------------------------------------------------------------------


def test_a_lower_node_beside_the_chain_takes_the_flow_and_drains_is_false(
    burn: ModuleType, gauge: ModuleType
) -> None:
    """The notch leaves the floor 0.01 m below it at (70, CC + 1): row 69 of
    the chain drains into the notch, not into row 70. The burn is right (the
    oracle holds); the chain is just not the drainage, and the sensitivity
    says so."""
    import tin_engine.sensitivity as sensitivity

    z = valley(dam=True, notch_drop=0.0)
    burnt, path = run(burn, z, reach(gauge, column_line(CC, 40, 200), at=290.0))
    assert nodes(path)[int(path.placed)] == (69, CC)
    assert nodes(path)[drains_along(flow_to(burnt), nodes(path))] == (69, CC)  # the premise
    a = accumulate(to_core(tile_of(burnt)))
    s = sensitivity.assess(
        np.asarray(a.count), np.asarray(a.reach), np.asarray(a.flow_to), path, CELL_KM2, 30.0, 1e6
    )
    assert s.drains is False
    assert "chain_not_draining" in s.causes
    assert s.well_posed is False


# ---------------------------------------------------------------------------
# The window
# ---------------------------------------------------------------------------


def test_a_reach_end_outside_the_window_is_dropped_and_reported(
    burn: ModuleType, gauge: ModuleType
) -> None:
    z = channel(60, 41, [(0, 20), (59, 20)])
    _, path = run(burn, z, reach(gauge, column_line(20, -10, 30), at=200.0))
    chain = nodes(path)
    assert chain[0] == (0, 20)
    assert chain[int(path.placed)] == (10, 20)
    assert path.dropped_m == pytest.approx(100.0, abs=10.0)


def test_a_reach_inside_the_window_drops_nothing(burn: ModuleType, gauge: ModuleType) -> None:
    z = channel(60, 41, [(0, 20), (59, 20)])
    _, path = run(burn, z, reach(gauge, column_line(20, 5, 30), at=100.0))
    assert path.dropped_m == 0.0


def test_a_placed_position_outside_the_window_is_refused(
    burn: ModuleType, gauge: ModuleType
) -> None:
    z = channel(60, 41, [(0, 20), (59, 20)])
    with pytest.raises(ValueError, match=r"(?i)outside"):
        burn.burn_reach(tile_of(z), reach(gauge, column_line(20, -10, 30), at=50.0))


# ---------------------------------------------------------------------------
# No NoData on the chain (step 2; PR 2's code review, round 1)
# ---------------------------------------------------------------------------

#: The two sentinels the design names, both in float32: DTM10's -32767, and
#: a positive one near float32's largest value, which the burn would carry
#: into `lowered_max_m`.
SENTINELS = [-32767.0, 3.4e38]
GAP = (20, 21)


def jog(gap: float | None) -> npt.NDArray[np.float32]:
    """A float32 floor falling 0.5 m per row down column 20 on rows 0 to 19
    and down column 22 from row 20 on, walls 3 m per column. `gap` is put at
    (20, 21), the node between the two floor runs; None leaves its data."""
    r, c = np.indices((60, 41)).astype(np.float64)
    floor_col = np.where(r < 20, 20.0, 22.0)
    z = (300.0 - 0.5 * r + 3.0 * np.abs(c - floor_col)).astype(np.float32)
    if gap is not None:
        z[GAP] = np.float32(gap)
    return z


def least_in_row(z: npt.NDArray[np.floating], row: int, nodata: float) -> tuple[int, int]:
    """The least-elevation node with data in `row` within 30 m of column 20:
    the cross-section of a resample point on the line down column 20."""
    cols = [c for c in range(17, 24) if z[row, c] != np.float32(nodata)]
    return row, min(cols, key=lambda c: (float(z[row, c]), abs(c - 20), c))


def jog_reach(gauge: ModuleType) -> Any:
    """The mapped line down column 20, rows 0 to 40; `at` = 200 m is row 20."""
    return reach(gauge, column_line(20, 0, 40), at=200.0, corridor=30.0)


@pytest.mark.parametrize("nodata", SENTINELS, ids=["minus32767", "3.4e38"])
def test_a_chain_through_nodata_is_refused_naming_the_gap(
    burn: ModuleType, gauge: ModuleType, nodata: float
) -> None:
    """The cross-sections of rows 19 and 20 choose (19, 20) and (20, 22),
    which have data; the straight 8-connected join between them steps through
    its middle node, (19.5, 21) rounded half up, (20, 21), which has none.
    The taut cut then drops (20, 22), so (20, 21) is a node of the chain the
    burn would lower and read. The station is refused, naming the gap and
    where it is in the DEM's CRS."""
    z = jog(nodata)
    # The premises, from the fixture alone.
    assert z[GAP] == np.float32(nodata)
    a, b = least_in_row(z, 19, nodata), least_in_row(z, 20, nodata)
    assert (a, b) == ((19, 20), (20, 22))
    middle = tuple(math.floor((p + q) / 2 + 0.5) for p, q in zip(a, b, strict=True))
    assert middle == GAP

    with pytest.raises(ValueError, match=r"(?i)gap \(NoData\) in the DEM") as raised:
        burn.burn_reach(tile_of(z, nodata), jog_reach(gauge))
    x, y = lat(GAP[1], GAP[0])  # (500210, 6599800)
    message = str(raised.value)
    assert f"{x:.0f}" in message and f"{y:.0f}" in message, message


def test_the_same_chain_with_data_at_the_join_burns(burn: ModuleType, gauge: ModuleType) -> None:
    """The control: with data at (20, 21) the chain is the same and burns
    without error, so the refusal is about a node on the chain, not about
    NoData anywhere in the window (the window keeps a NoData node off the
    chain, at (50, 5))."""
    z = jog(None)
    z[50, 5] = np.float32(NODATA)
    _, path = run(burn, z, jog_reach(gauge), nodata=NODATA)
    assert GAP in nodes(path)
