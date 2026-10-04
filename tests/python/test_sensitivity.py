"""`sensitivity.assess`: is the area well defined at this gauge? (increment 29, PR 2).

`docs/increments/29-nve-reference-catchments.md`, "Sensitivity" (steps 1 to
5) and "The red suites", PR 2's `test_sensitivity.py`. Arrays in, a frozen
`Sensitivity` out: the counts, flag bits and `flow_to` are written by hand
along a chain down column 1 of a small grid, one node per 10 m.

Interface pinned here (the design names `assess(count, reach_bits, flow_to,
path, cell_area_km2, uncertainty, reach_down_m)`, called positionally here,
and `SWING_MAX = 0.05`; the field names it leaves open are chosen here and
listed in the handback):

- `path` is read for `chain`, `placed`, `arc` and `end_closed` only, so a
  stand-in with those four attributes is passed, not `burn.GaugePath`.
- `Sensitivity`: `a0`, `area_up`, `area_down` (km2), `swing`,
  `largest_step` (km2) and `largest_step_at_m` (metres from the placed node,
  negative upstream), `checked_up_m`, `checked_down_m`, `drains`,
  `monotone`, `causes` (a tuple of `"swing"`, `"downstream_unread"`,
  `"chain_not_draining"`, `"chain_end_open"`) and `well_posed`.

`U` is 30 m except in the tests of a `U` between nodes (35 m on 10 m steps),
where "read to `U`" means every sample in `(0, U]` trusted and a chain node at
arc length `>= U` ("Sensitivity", step 2), not `checked_down_m >= U`.

The cell area is 2^-10 km2, so every area is a count times a power of two and
the swing's arithmetic, `(A_down - A0) / A0`, is exact up to the final
division: a swing of exactly 0.05 compares equal to `SWING_MAX`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from types import ModuleType

import numpy as np
import numpy.typing as npt
import pytest

from gauge_fixtures import DIRECTION, OUTLET

CELL = 2.0**-10
U = 30.0
STEP = 10.0
SOUTH = DIRECTION[(1, 0)]
EAST_SOUTH = DIRECTION[(1, 1)]


@pytest.fixture(scope="module")
def sensitivity() -> ModuleType:
    import tin_engine.sensitivity as module

    return module


@dataclass(frozen=True)
class Path:
    chain: npt.NDArray[np.int64]
    placed: int
    arc: npt.NDArray[np.float64]
    end_closed: bool = True
    end_extended_m: float = 0.0


@dataclass
class Case:
    count: npt.NDArray[np.uint32]
    bits: npt.NDArray[np.uint8]
    flow_to: npt.NDArray[np.uint8]
    path: Path


def case(
    counts: Sequence[int],
    placed: int,
    *,
    bits: dict[int, int] | None = None,
    off_chain: Sequence[int] = (),
    end_closed: bool = True,
) -> Case:
    """A chain down column 1, node k at row k, `counts[k]` its count; each node
    drains into the next (`flow_to` south) unless its index is in `off_chain`
    (then south-east, off the chain); `bits[k]` sets node k's flag bits."""
    n = len(counts)
    count = np.zeros((n + 1, 3), dtype=np.uint32)
    reach = np.zeros((n + 1, 3), dtype=np.uint8)
    flow = np.full((n + 1, 3), OUTLET, dtype=np.uint8)
    for k, c in enumerate(counts):
        count[k, 1] = c
        flow[k, 1] = EAST_SOUTH if k in off_chain else SOUTH
    for k, b in (bits or {}).items():
        reach[k, 1] = b
    chain = np.array([(k, 1) for k in range(n)], dtype=np.int64)
    arc = STEP * (np.arange(n, dtype=np.float64) - placed)
    return Case(count, reach, flow, Path(chain, placed, arc, end_closed))


def assess(sensitivity: ModuleType, c: Case, u: float = U, reach_down_m: float = 1000.0) -> object:
    return sensitivity.assess(c.count, c.bits, c.flow_to, c.path, CELL, u, reach_down_m)


def test_swing_max_is_the_match_bar(sensitivity: ModuleType) -> None:
    assert sensitivity.SWING_MAX == 0.05


def test_a_smooth_gain_is_well_posed(sensitivity: ModuleType) -> None:
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4))
    assert s.a0 == pytest.approx(1000 * CELL, rel=1e-15)
    assert s.area_up == pytest.approx(990 * CELL, rel=1e-15)
    assert s.area_down == pytest.approx(1009 * CELL, rel=1e-15)
    assert s.swing == pytest.approx(0.01, rel=1e-12)
    assert s.checked_up_m == pytest.approx(30.0)
    assert s.checked_down_m == pytest.approx(30.0)
    assert s.drains and s.monotone
    assert tuple(s.causes) == ()
    assert s.well_posed


def test_a_confluence_step_of_known_size_at_a_known_place(sensitivity: ModuleType) -> None:
    """A tributary of 20 % joins 20 m below the placed node."""
    counts = [994, 997, 999, 1000, 1001, 1202, 1203, 1204, 1205]
    s = assess(sensitivity, case(counts, placed=3))
    assert s.swing == pytest.approx((1203 - 1000) / 1000, rel=1e-12)
    assert s.largest_step == pytest.approx((1202 - 1001) * CELL, rel=1e-12)
    assert 10.0 <= s.largest_step_at_m <= 20.0
    assert "swing" in s.causes
    assert s.well_posed is False


def test_a_flat_floor_upstream_gives_a_step_there(sensitivity: ModuleType) -> None:
    """Where a flat's nodes join the chain the count jumps; here 20 m above."""
    counts = [590, 595, 600, 990, 1000, 1002, 1004, 1006, 1008]
    s = assess(sensitivity, case(counts, placed=4))
    assert s.largest_step == pytest.approx((990 - 600) * CELL, rel=1e-12)
    assert -20.0 <= s.largest_step_at_m <= -10.0
    assert s.swing == pytest.approx((1000 - 595) / 1000, rel=1e-12)
    assert "swing" in s.causes


@pytest.mark.parametrize(("down", "posed"), [(1050, True), (1051, False)])
def test_a_swing_of_exactly_005_is_well_posed_and_just_above_is_not(
    sensitivity: ModuleType, down: int, posed: bool
) -> None:
    counts = [985, 990, 995, 998, 1000, 1010, 1030, down, down + 5]
    s = assess(sensitivity, case(counts, placed=4))
    assert s.swing == (down - 1000) / 1000
    assert s.well_posed is posed
    assert ("swing" in s.causes) is not posed


def test_the_swing_is_one_sided(sensitivity: ModuleType) -> None:
    """0.04 below and 0.04 above: 0.04, not 0.08, and well posed."""
    counts = [950, 960, 975, 990, 1000, 1010, 1025, 1040, 1050]
    s = assess(sensitivity, case(counts, placed=4))
    assert s.swing == pytest.approx(0.04, rel=1e-12)
    assert s.well_posed


# ---------------------------------------------------------------------------
# Which counts are trusted: by their flag bits, on both sides
# ---------------------------------------------------------------------------


def test_a_flag_upstream_stops_the_read_there(sensitivity: ModuleType) -> None:
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4, bits={2: 1}))
    assert s.checked_up_m == pytest.approx(10.0)
    assert s.area_up == pytest.approx(997 * CELL, rel=1e-15)
    assert s.checked_down_m == pytest.approx(30.0)


def test_the_read_stops_at_the_first_flag_not_at_a_position(sensitivity: ModuleType) -> None:
    """Node -10 m is flagged, -20 m and -30 m are not: the read stops at the
    placed node, and the far counts (much smaller) are not used."""
    counts = [100, 200, 300, 999, 1000, 1003, 1006, 1010, 1012]
    s = assess(sensitivity, case(counts, placed=4, bits={3: 1}))
    assert s.checked_up_m == pytest.approx(0.0)
    assert s.swing == pytest.approx(0.01, rel=1e-12)


@pytest.mark.parametrize("bit", [1, 2, 3])
def test_a_flag_downstream_before_u_is_downstream_unread(sensitivity: ModuleType, bit: int) -> None:
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4, bits={6: bit}))
    assert s.checked_down_m == pytest.approx(10.0)
    assert s.area_down == pytest.approx(1003 * CELL, rel=1e-15)
    assert "downstream_unread" in s.causes
    assert s.well_posed is False


def test_a_flag_past_u_does_not_shorten_the_read(sensitivity: ModuleType) -> None:
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4, bits={8: 1, 0: 2}))
    assert s.checked_up_m == pytest.approx(30.0)
    assert s.checked_down_m == pytest.approx(30.0)
    assert s.well_posed


def test_a_mapped_river_that_ends_before_u_is_downstream_unread(
    sensitivity: ModuleType,
) -> None:
    """Every count is trusted, but the mapped river stops 25 m below."""
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4), reach_down_m=25.0)
    assert "downstream_unread" in s.causes
    assert s.well_posed is False


def test_a_mapped_river_that_reaches_u_is_read(sensitivity: ModuleType) -> None:
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4), reach_down_m=30.0)
    assert "downstream_unread" not in s.causes


def test_u_between_nodes_is_read_when_a_chain_node_lies_at_or_past_it(
    sensitivity: ModuleType,
) -> None:
    """U = 35 m: the samples downstream are at 10, 20 and 30 m, all trusted,
    and the chain's next node, at 40 m, is past U. That is read to U, though
    the farthest sample is 30 m."""
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4), u=35.0)
    assert s.checked_down_m == pytest.approx(30.0)
    assert s.checked_up_m == pytest.approx(30.0)
    assert s.area_down == pytest.approx(1009 * CELL, rel=1e-15)
    assert "downstream_unread" not in s.causes
    assert s.well_posed


def test_the_flags_of_the_node_past_u_do_not_matter(sensitivity: ModuleType) -> None:
    """U = 35 m: the node at 40 m is D, not a sample; its flag bit does not
    shorten the read."""
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4, bits={8: 1}), u=35.0)
    assert s.checked_down_m == pytest.approx(30.0)
    assert "downstream_unread" not in s.causes
    assert s.well_posed


def test_u_between_nodes_is_unread_when_the_chain_ends_short_of_it(
    sensitivity: ModuleType,
) -> None:
    """U = 35 m and the chain's last node is at 30 m: every sample is
    trusted, but no chain node lies at or past U, so the read does not reach
    it."""
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009]
    s = assess(sensitivity, case(counts, placed=4), u=35.0)
    assert s.checked_down_m == pytest.approx(30.0)
    assert "downstream_unread" in s.causes
    assert s.well_posed is False


def test_a_river_that_starts_within_u_upstream_is_not_penalised(
    sensitivity: ModuleType,
) -> None:
    counts = [997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=1))
    assert s.checked_up_m == pytest.approx(10.0)
    assert s.area_up == pytest.approx(997 * CELL, rel=1e-15)
    assert tuple(s.causes) == ()
    assert s.well_posed


# ---------------------------------------------------------------------------
# monotone and drains
# ---------------------------------------------------------------------------


def test_two_equal_neighbouring_counts_are_not_monotone(sensitivity: ModuleType) -> None:
    counts = [988, 990, 994, 997, 1000, 1000, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4))
    assert s.monotone is False
    assert s.well_posed is False


def test_a_fall_is_not_monotone(sensitivity: ModuleType) -> None:
    counts = [988, 990, 994, 997, 1000, 1003, 1001, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4))
    assert s.monotone is False
    assert s.well_posed is False


def test_an_open_end_is_not_monotone(sensitivity: ModuleType) -> None:
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4, end_closed=False))
    assert s.monotone is False
    assert "chain_end_open" in s.causes
    assert s.well_posed is False


def test_a_read_sample_draining_off_the_chain_is_not_draining(sensitivity: ModuleType) -> None:
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012]
    s = assess(sensitivity, case(counts, placed=4, off_chain=(5,)))
    assert s.drains is False
    assert "chain_not_draining" in s.causes
    assert s.well_posed is False


def test_only_read_samples_must_drain(sensitivity: ModuleType) -> None:
    """Node +40 m (index 8) drains off the chain, but it is past U."""
    counts = [988, 990, 994, 997, 1000, 1003, 1006, 1009, 1012, 1015]
    s = assess(sensitivity, case(counts, placed=4, off_chain=(8,)))
    assert s.drains
    assert s.well_posed


def test_each_cause_is_counted_apart(sensitivity: ModuleType) -> None:
    """A big step upstream, a flag downstream, an off-chain drain inside the
    read and an open end: all four causes, each once."""
    counts = [500, 600, 700, 990, 1000, 1003, 1006, 1009, 1012]
    c = case(counts, placed=4, bits={6: 2}, off_chain=(3,), end_closed=False)
    s = assess(sensitivity, c)
    assert set(s.causes) == {"swing", "downstream_unread", "chain_not_draining", "chain_end_open"}
    assert len(tuple(s.causes)) == 4
    assert s.well_posed is False
