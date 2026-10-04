"""DC0: the partition rule of increment 23c, pure (no DEM, no mesh).

`docs/increments/23-basin-scale.md`, "The partition (23c; B2, B3 and B14 as
ruled)" and "The memory estimate and the defaults (B14 as ruled)"; the
figures are the ones `docs/increments/23-probes/partition.py` prints, and the
cell count approximates the request and is not a floor (the review of B13 and
B14, round 1).

The surface, pinned at the red step and confirmed under "Settled after 23c's
red step (65e3990)", item 7: `tin_engine.decompose` holds

- `DEFAULT_MEMORY_BUDGET`, 16 GiB in bytes (2**34);
- `bytes_per_node(tolerance) -> int`, the `b(T)` table, linear between its
  columns and rounded up, 19 above 50 m;
- `partition(cols, rows, tolerance, *, pieces=1, memory_budget=DEFAULT_MEMORY_BUDGET)`
  returning an object with integer `nx`, `ny`, `dx`, `dy` (cells across and
  down, and their sides in whole lattice nodes).

Two refusals the design does not state, pinned as `ValueError`: `pieces < 1`,
and a budget under one node's bytes (`memory_budget < b(T)`), for which step
3's loop would never end (no cell of one node fits).

Tolerances in the interpolation cases are binary fractions (0.25, 0.5, 0.75,
1.5, 7.5), so the rounding up is of an exact value and a correct
implementation cannot land one above it by floating-point noise.
"""

from __future__ import annotations

import importlib
import itertools
import math
import os
from types import ModuleType

import numpy as np
import pytest

GIB = 1 << 30

#: The design's table, restated here from the record (not imported).
TABLE = ((0.0, 637), (1.0, 267), (2.0, 155), (5.0, 65), (10.0, 37), (20.0, 25), (50.0, 19))


@pytest.fixture
def dec() -> ModuleType:
    return importlib.import_module("tin_engine.decompose")


def cells(p: object) -> tuple[int, int, int, int]:
    nx, ny, dx, dy = (getattr(p, k) for k in ("nx", "ny", "dx", "dy"))
    for v in (nx, ny, dx, dy):
        assert isinstance(v, int) and not isinstance(v, bool), (nx, ny, dx, dy)
    return nx, ny, dx, dy


# ------------------------------------------------------------------ b(T)


class TestBytesPerNode:
    @pytest.mark.parametrize(("tolerance", "expected"), TABLE)
    def test_the_tables_columns(self, dec: ModuleType, tolerance: float, expected: int) -> None:
        assert dec.bytes_per_node(tolerance) == expected

    @pytest.mark.parametrize(
        ("tolerance", "expected"),
        [(0.5, 452), (0.25, 545), (0.75, 360), (1.5, 211), (7.5, 51)],
    )
    def test_linear_between_columns_rounded_up(
        self, dec: ModuleType, tolerance: float, expected: int
    ) -> None:
        """0.5 m gives 452 (the design's own example); 0.25 m is 544.5 exact,
        so 545 is the rounding up, not to nearest-even."""
        assert dec.bytes_per_node(tolerance) == expected

    @pytest.mark.parametrize("tolerance", [50.0, 50.5, 100.0, 1e6])
    def test_nineteen_from_fifty_metres_up(self, dec: ModuleType, tolerance: float) -> None:
        assert dec.bytes_per_node(tolerance) == 19

    def test_never_increases_with_the_tolerance(self, dec: ModuleType) -> None:
        values = [dec.bytes_per_node(t) for t in np.linspace(0.0, 60.0, 601)]
        assert all(a >= b for a, b in itertools.pairwise(values))

    def test_the_default_budget_is_sixteen_gib(self, dec: ModuleType) -> None:
        assert dec.DEFAULT_MEMORY_BUDGET == 16 * GIB == 2**34


# ------------------------------------------------------------------ the figures


class TestThePartitionProbesFigures:
    """`partition.py`'s rows under "What was measured for this design"."""

    def test_velhas_at_one_metre_with_pieces_64(self, dec: ModuleType) -> None:
        assert cells(dec.partition(4208, 7347, 1.0, pieces=64)) == (6, 11, 702, 668)

    def test_velhas_at_one_metre_with_pieces_16_gives_fifteen_cells(self, dec: ModuleType) -> None:
        """The count approximates the request and is not a floor: 3 x 5 = 15."""
        nx, ny, _, _ = cells(dec.partition(4208, 7347, 1.0, pieces=16))
        assert (nx, ny) == (3, 5)
        assert nx * ny == 15 < 16

    def test_velhas_is_one_piece_at_half_a_metre_and_up(self, dec: ModuleType) -> None:
        for tolerance in (0.5, 1.0, 5.0, 50.0):
            assert cells(dec.partition(4208, 7347, tolerance)) == (1, 1, 4208, 7347)

    def test_velhas_at_zero_is_cut_in_two(self, dec: ModuleType) -> None:
        assert cells(dec.partition(4208, 7347, 0.0))[:2] == (1, 2)

    def test_the_basin_at_one_metre(self, dec: ModuleType) -> None:
        assert cells(dec.partition(41332, 50297, 1.0)) == (5, 7, 8267, 7186)

    def test_the_basin_at_fifty_metres(self, dec: ModuleType) -> None:
        assert cells(dec.partition(41332, 50297, 50.0))[:2] == (2, 2)


# ------------------------------------------------------------------ the boundary


class TestTheBudgetBoundary:
    @pytest.mark.parametrize("tolerance", [0.0, 1.0, 10.0, 50.0])
    def test_one_piece_at_n_b_equal_to_the_budget(self, dec: ModuleType, tolerance: float) -> None:
        cols, rows = 300, 200
        budget = cols * rows * dec.bytes_per_node(tolerance)
        assert cells(dec.partition(cols, rows, tolerance, memory_budget=budget)) == (
            1,
            1,
            cols,
            rows,
        )

    @pytest.mark.parametrize("tolerance", [0.0, 1.0, 10.0, 50.0])
    def test_a_cut_at_n_b_one_byte_over_the_budget(self, dec: ModuleType, tolerance: float) -> None:
        cols, rows = 300, 200
        b = dec.bytes_per_node(tolerance)
        budget = cols * rows * b - 1
        nx, ny, dx, dy = cells(dec.partition(cols, rows, tolerance, memory_budget=budget))
        assert nx * ny > 1
        assert dx * dy * b <= budget

    def test_a_large_budget_gives_one_piece_above_every_former_cap(self, dec: ModuleType) -> None:
        """B14: no N_MAX (2**22 nodes) and no per-piece cap; the basin at 0 m
        is one piece under a budget that holds it."""
        cols, rows = 41332, 50297
        budget = cols * rows * 637
        assert cols * rows > 2**22
        assert cells(dec.partition(cols, rows, 0.0, memory_budget=budget)) == (1, 1, cols, rows)

    def test_pieces_one_is_no_request(self, dec: ModuleType) -> None:
        assert cells(dec.partition(64, 64, 1.0, pieces=1)) == (1, 1, 64, 64)

    def test_a_request_cuts_a_small_window(self, dec: ModuleType) -> None:
        """`--pieces` asks for more pieces than the budget needs: 4 on a
        square window is 2 x 2."""
        assert cells(dec.partition(41, 41, 1.0, pieces=4)) == (2, 2, 21, 21)


class TestRefusals:
    @pytest.mark.parametrize("pieces", [0, -1])
    def test_pieces_under_one(self, dec: ModuleType, pieces: int) -> None:
        with pytest.raises(ValueError, match="pieces"):
            dec.partition(64, 64, 1.0, pieces=pieces)

    def test_a_budget_under_one_node(self, dec: ModuleType) -> None:
        """No cell of one node fits, so step 3 cannot end: refused, not looped."""
        b = dec.bytes_per_node(1.0)
        with pytest.raises(ValueError, match="budget"):
            dec.partition(64, 64, 1.0, memory_budget=b - 1)

    def test_a_budget_of_one_node_is_every_node_a_cell(self, dec: ModuleType) -> None:
        b = dec.bytes_per_node(1.0)
        assert cells(dec.partition(5, 4, 1.0, memory_budget=b)) == (5, 4, 1, 1)


# ------------------------------------------------------------------ the invariants


def random_cases(seed: int, n: int) -> list[tuple[int, int, float, int, int]]:
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        cols = int(rng.integers(2, 60_000))
        rows = int(rng.integers(2, 60_000))
        tolerance = float(rng.choice([0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 33.0, 80.0]))
        pieces = int(rng.choice([1, 1, 2, 3, 4, 7, 16, 64, 100, 1000]))
        # Budgets from a few nodes' worth to far above the window.
        budget = int(637 * 2 ** rng.uniform(0, 42))
        out.append((cols, rows, tolerance, pieces, budget))
    return out


class TestInvariantsOverRandomWindows:
    """For windows, tolerances, `--pieces` and budgets drawn at random: no
    cell over the budget, no empty row or column of cells, whole-node sides,
    and the count at least what the budget needs."""

    def test_seeded_draws(self, dec: ModuleType) -> None:
        checked = 0
        for cols, rows, tolerance, pieces, budget in random_cases(23, 3000):
            b = dec.bytes_per_node(tolerance)
            if budget < b:
                continue
            nx, ny, dx, dy = cells(
                dec.partition(cols, rows, tolerance, pieces=pieces, memory_budget=budget)
            )
            case = (cols, rows, tolerance, pieces, budget, nx, ny, dx, dy)
            assert dx * dy * b <= budget, case
            assert 1 <= dx <= cols and 1 <= dy <= rows, case
            # No empty column or row of cells: the last starts inside the window.
            assert (nx - 1) * dx < cols <= nx * dx, case
            assert (ny - 1) * dy < rows <= ny * dy, case
            if nx * ny == 1:
                assert (dx, dy) == (cols, rows), case
                assert max(pieces, math.ceil(cols * rows * b / budget)) <= 1, case
            checked += 1
        assert checked > 2500

    def test_one_piece_exactly_when_neither_the_request_nor_the_budget_asks(
        self, dec: ModuleType
    ) -> None:
        for cols, rows, tolerance, pieces, budget in random_cases(29, 2000):
            b = dec.bytes_per_node(tolerance)
            if budget < b:
                continue
            nx, ny, _, _ = cells(
                dec.partition(cols, rows, tolerance, pieces=pieces, memory_budget=budget)
            )
            wanted = max(pieces, -(-cols * rows * b // budget))
            assert (nx * ny == 1) == (wanted <= 1), (cols, rows, tolerance, pieces, budget)


class TestNotTheMachines:
    """K5: neither the count nor the cut is the machine's."""

    @pytest.mark.parametrize(
        ("cores", "memory"), [(1, 1 * GIB), (512, 4096 * GIB), (None, 2 * GIB)]
    )
    def test_the_same_partition_on_any_machine(
        self, dec: ModuleType, monkeypatch: pytest.MonkeyPatch, cores: int | None, memory: int
    ) -> None:
        cases = [(4208, 7347, 1.0, 64), (41332, 50297, 1.0, 1), (41332, 50297, 0.0, 1)]
        before = [cells(dec.partition(c, r, t, pieces=p)) for c, r, t, p in cases]
        import tin_engine.mosaic as mosaic

        monkeypatch.setattr(os, "cpu_count", lambda: cores)
        monkeypatch.setattr(mosaic, "physical_memory", lambda: memory, raising=False)
        if hasattr(os, "sched_getaffinity"):
            monkeypatch.setattr(os, "sched_getaffinity", lambda pid: set(range(cores or 1)))
        after = [cells(dec.partition(c, r, t, pieces=p)) for c, r, t, p in cases]
        assert after == before
