"""`tin_engine.mosaic`: planning and assembly of a multi-tile DEM (increment 15a).

The design is `docs/increments/15-dem-mosaic.md` (R1, R4, R5, R7; invariants
I1-I7); names not fixed there are pinned in its "Pinned by the red suite (15a)"
subsection. Test ids follow its "Tests for @tester" list, which carries the
parked design's M1-M12 (`git show f8cedbe:docs/increments/15-dem-repository.md`)
with 15a's changes: M2 is lattice grouping, M5 the cap from physical memory,
and M13-M15 are new.

This is 15a's **invariant-critical** suite. Its tests were run against a
scratch implementation and against deliberate mutants of it (a wrong lattice
alignment check, an off-by-one in stitching offsets, an overlap check that
accepts a disagreement, a coverage check that misses a gap); the handback of
the red commit lists which test killed which.

HOW THIS FILE GOES RED. `tin_engine.mosaic` and `tin_engine.io.repository`
are imported in module-scoped fixtures, as `test_bench.py` loads its tool, so
while they are missing each test fails on its own with `ModuleNotFoundError`
and the rest of `tests/python` still collects and runs.

Tiles are built from `RasterMeta` and `DemTile` directly (`mosaic_fixtures`),
and the repository is a dict: nothing here touches a file.
"""

from __future__ import annotations

import importlib
import itertools
import math
import re
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import shapely

from mosaic_fixtures import (
    DX,
    DY,
    SENTINEL,
    X0,
    Y0,
    Loads,
    blocks,
    footprints,
    meta,
    piece,
    quadrants,
    same_array,
    values,
    whole,
)
from tin_engine.io.models import DemTile

SRC = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine"


@pytest.fixture(scope="module")
def mz() -> ModuleType:
    return importlib.import_module("tin_engine.mosaic")


@pytest.fixture(scope="module")
def footprint() -> Any:
    return importlib.import_module("tin_engine.io.repository").TileFootprint


@pytest.fixture
def plan(mz: ModuleType, footprint: Any) -> Any:
    """`plan(tiles, box=None, needed=None)`: plan a dict of tiles; `box` is a 4-tuple."""

    def run(
        tiles: dict[str, DemTile], box: tuple[float, ...] | None = None, needed: Any = None
    ) -> Any:
        bounds = None if box is None else bounds_of(mz, box)
        return mz.plan_mosaic(footprints(footprint, tiles), bounds, needed)

    return run


@pytest.fixture
def build(mz: ModuleType, plan: Any) -> Any:
    """`build(tiles, box=None, needed=None) -> (Mosaic, Loads)`: plan, then assemble."""

    def run(
        tiles: dict[str, DemTile], box: tuple[float, ...] | None = None, needed: Any = None
    ) -> Any:
        loads = Loads(tiles)
        return mz.assemble(plan(tiles, box, needed), loads), loads

    return run


@pytest.fixture
def refused(mz: ModuleType) -> Any:
    """`refused(call, *tokens)`: `call()` raises `MosaicError` naming each token."""

    def check(call: Any, *tokens: str) -> str:
        with pytest.raises(mz.MosaicError) as info:
            call()
        message = str(info.value)
        for token in tokens:
            assert token.lower() in message.lower(), f"{token!r} not named in {message!r}"
        return message

    return check


def bounds_of(mz: ModuleType, box: tuple[float, ...]) -> Any:
    x_min, y_min, x_max, y_max = box
    return mz.Bounds(x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max)


def window(mz: ModuleType, row0: int, col0: int, rows: int, cols: int) -> Any:
    return mz.IndexWindow(row0=row0, col0=col0, rows=rows, cols=cols)


def names_number(message: str, number: int) -> None:
    """`number` appears as a whole number, not inside another (`3` is not in `317`)."""
    assert re.search(rf"(?<![\d.]){number}(?![\d])", message), f"{number} not named in {message!r}"


def global_origin(plan: Any, name: str) -> tuple[int, int]:
    """The global index of tile `name`'s first node, as the plan places it."""
    (p,) = [p for p in plan.tiles if p.name == name]
    return (
        plan.window.row0 + p.canvas.row0 - p.source.row0,
        plan.window.col0 + p.canvas.col0 - p.source.col0,
    )


# ---------------------------------------------------------------------------
# The module's surface
# ---------------------------------------------------------------------------


def test_mosaic_error_is_a_value_error(mz: ModuleType) -> None:
    """R11: refusals become `typer.BadParameter` like `GeoTiffError`, one type."""
    assert issubclass(mz.MosaicError, ValueError)


def test_align_tolerance_is_a_millionth_of_a_cell(mz: ModuleType) -> None:
    assert mz.ALIGN_TOLERANCE == 1e-6


def test_an_empty_repository_is_refused(mz: ModuleType, refused: Any) -> None:
    refused(lambda: mz.plan_mosaic([], None, None))


class TestBounds:
    """R6: `Bounds` is four finite floats with `x_min < x_max`, `y_min < y_max`."""

    def test_a_valid_box(self, mz: ModuleType) -> None:
        box = bounds_of(mz, (1.0, 2.0, 3.0, 4.0))
        assert (box.x_min, box.y_min, box.x_max, box.y_max) == (1.0, 2.0, 3.0, 4.0)

    @pytest.mark.parametrize(
        "box",
        [
            (math.nan, 0.0, 1.0, 1.0),
            (0.0, 0.0, math.inf, 1.0),
            (0.0, -math.inf, 1.0, 1.0),
            (1.0, 0.0, 0.0, 1.0),  # inverted in x
            (0.0, 1.0, 1.0, 0.0),  # inverted in y
            (0.0, 0.0, 0.0, 1.0),  # zero width
            (0.0, 1.0, 1.0, 1.0),  # zero height
        ],
        ids=["nan", "inf", "-inf", "inverted-x", "inverted-y", "zero-width", "zero-height"],
    )
    def test_refused(self, mz: ModuleType, box: tuple[float, ...]) -> None:
        with pytest.raises(ValueError):
            bounds_of(mz, box)


# ---------------------------------------------------------------------------
# Planning (no pixels)
# ---------------------------------------------------------------------------


class TestM1Layout:
    """M1: canvas shape, origin and each tile's two windows."""

    def test_point_registered_2x2_sharing_a_row_and_a_column(
        self, mz: ModuleType, plan: Any
    ) -> None:
        source = whole()  # 9 x 13
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        result = plan(tiles)
        assert result.meta == source.meta
        assert result.window == window(mz, 0, 0, 9, 13)
        assert [p.name for p in result.tiles] == ["ne.tif", "nw.tif", "se.tif", "sw.tif"]
        placed = {p.name: p for p in result.tiles}
        expected = {"nw.tif": (0, 0), "ne.tif": (0, 6), "sw.tif": (4, 0), "se.tif": (4, 6)}
        for name, (row0, col0) in expected.items():
            assert placed[name].canvas == window(mz, row0, col0, 5, 7), name
            assert placed[name].source == window(mz, 0, 0, 5, 7), name
            assert placed[name].meta == tiles[name].meta, name

    def test_area_registered_2x2_abutting(self, mz: ModuleType, plan: Any) -> None:
        """Disjoint node sets one spacing apart: the canvas is contiguous and
        nothing between them is a gap (the coverage check must not see one)."""
        source = whole(rows=8, cols=12, area=True)
        result = plan(quadrants(source, row_cut=4, col_cut=6, overlap=0))
        assert result.meta == source.meta
        placed = {p.name: p.canvas for p in result.tiles}
        assert placed == {
            "nw.tif": window(mz, 0, 0, 4, 6),
            "ne.tif": window(mz, 0, 6, 4, 6),
            "sw.tif": window(mz, 4, 0, 4, 6),
            "se.tif": window(mz, 4, 6, 4, 6),
        }

    def test_one_tile_plans_as_itself(self, mz: ModuleType, plan: Any) -> None:
        source = whole()
        result = plan({"only.tif": source})
        assert result.meta == source.meta
        assert result.reference == (X0, Y0)
        assert result.window == window(mz, 0, 0, 9, 13)


class TestM2LatticeGrouping:
    """M2, changed for 15a (R4 points 1 and 3, Q5): lattices, not one lattice.

    The main lattice is the 9 x 13 point grid cut in four. `odd.tif` lies
    south-east of it, half a cell east-west off (N2's case), overlapping `se`.
    """

    @staticmethod
    def repository(**odd_changes: Any) -> dict[str, DemTile]:
        tiles = quadrants(whole(), row_cut=4, col_cut=6, overlap=1)
        fields: dict[str, Any] = {"x_min": X0 + 10.5 * DX, "y_max": Y0 - 6 * DY}
        fields.update(odd_changes)
        tiles["odd.tif"] = whole(rows=5, cols=7, **fields)
        return tiles

    def test_a_request_inside_the_main_lattice_plans(self, plan: Any) -> None:
        result = plan(self.repository(), (500012.0, 6599987.0, 500033.0, 6599996.0))
        assert [p.name for p in result.tiles] == ["nw.tif"]
        assert result.reference == (X0, Y0)

    def test_a_request_inside_the_odd_lattice_plans_on_that_lattice(
        self, mz: ModuleType, plan: Any
    ) -> None:
        result = plan(self.repository(), (500130.0, 6599952.0, 500160.0, 6599958.0))
        assert [p.name for p in result.tiles] == ["odd.tif"]
        assert result.reference == (X0 + 10.5 * DX, Y0 - 6 * DY)
        assert result.meta.x_min == 500125.0  # odd lattice: 500105 + 2 * 10
        assert result.window == window(mz, 2, 2, 3, 5)

    def test_a_request_straddling_both_is_refused_naming_one_of_each(
        self, plan: Any, refused: Any
    ) -> None:
        tiles = self.repository()
        refused(
            lambda: plan(tiles, (500100.0, 6599955.0, 500130.0, 6599965.0)),
            "se.tif", "odd.tif", "0.5 cell", "east-west",
        )  # fmt: skip

    def test_no_bounds_selects_both_lattices_and_is_refused(self, plan: Any, refused: Any) -> None:
        tiles = self.repository()
        refused(lambda: plan(tiles), "odd.tif", "0.5 cell")

    def test_a_half_cell_north_south_is_named_so(self, plan: Any, refused: Any) -> None:
        tiles = self.repository(x_min=X0 + 11 * DX, y_max=Y0 - 6.5 * DY)
        refused(
            lambda: plan(tiles, (500100.0, 6599955.0, 500130.0, 6599965.0)),
            "odd.tif",
            "0.5 cell",
            "north-south",
        )

    def test_a_thousandth_of_a_cell_is_another_lattice(self, plan: Any, refused: Any) -> None:
        tiles = self.repository(x_min=X0 + (11 + 1e-3) * DX)
        refused(
            lambda: plan(tiles, (500100.0, 6599955.0, 500130.0, 6599965.0)),
            "se.tif",
            "odd.tif",
            "cell",
        )

    @pytest.mark.parametrize(
        ("change", "tokens"),
        [
            ({"epsg": 25832}, ("25833", "25832")),
            ({"epsg": 32633}, ("25833", "32633")),
            ({"dx": 12.5}, ("spacing", "12.5")),
            ({"dy": 7.5}, ("spacing", "7.5")),
            ({"area": True}, ("registration",)),
            ({"nodata": -9999.0, "nodata_source": "tag"}, ("nodata", "-9999", "none")),
        ],
        ids=["crs-25832", "crs-32633", "spacing-dx", "spacing-dy-alone", "registration", "nodata"],
    )
    def test_other_differences_are_named_as_such(
        self, plan: Any, refused: Any, change: dict[str, Any], tokens: tuple[str, ...]
    ) -> None:
        """The parked `refuses_mixed_crs_mosaic`, `_spacing`, `_registration`, `_nodata`.

        `ne.tif` differs from `nw.tif` in exactly one field and is otherwise
        where `quadrants` puts it, so only that field can be the reason.
        """
        source = whole()
        tiles = {
            "nw.tif": piece(source, 0, 5, 0, 7),
            "ne.tif": piece(source, 0, 5, 6, 13, **_meta_change(change)),
        }
        refused(lambda: plan(tiles), "nw.tif", "ne.tif", *tokens)

    def test_mixed_nodata_sentinels_are_refused_q4(self, plan: Any, refused: Any) -> None:
        """Q4 (a): two sentinels, both named."""
        source = whole(nodata=SENTINEL)
        tiles = {
            "nw.tif": piece(source, 0, 5, 0, 7),
            "ne.tif": piece(source, 0, 5, 6, 13, nodata=-9999.0),
        }
        refused(lambda: plan(tiles), "nw.tif", "ne.tif", "-32767", "-9999")

    def test_a_different_crs_outside_the_request_is_not_a_refusal(self, plan: Any) -> None:
        """R4 changed the parked rule: only the *selected* tiles share a lattice."""
        tiles = quadrants(whole(), row_cut=4, col_cut=6, overlap=1)
        tiles["far.tif"] = whole(rows=3, cols=3, x_min=X0 + 1000 * DX, epsg=25832)
        result = plan(tiles, (500012.0, 6599987.0, 500033.0, 6599996.0))
        assert [p.name for p in result.tiles] == ["nw.tif"]


def _meta_change(change: dict[str, Any]) -> dict[str, Any]:
    """`piece` takes `RasterMeta` field names."""
    names = {"dx": "delta_x", "dy": "delta_y", "area": "pixel_is_area"}
    return {names.get(key, key): value for key, value in change.items()}


class TestM3AlignmentTolerance:
    """M3 and R4 point 1: offsets within `ALIGN_TOLERANCE` cells of an integer."""

    def test_a_decimetre_grid_with_float_noise_plans_on_integer_windows(
        self, mz: ModuleType, plan: Any
    ) -> None:
        source = whole(x_min=600_000.0, y_max=7_000_000.0, dx=0.1, dy=0.1)
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        # A producer that accumulates its tie point: 6 steps of 0.1 is not 0.6.
        accumulated = 600_000.0
        for _ in range(6):
            accumulated += 0.1
        tiles["ne.tif"] = piece(source, 0, 5, 6, 13, x_min=accumulated)
        assert (accumulated - 600_000.0) / 0.1 != 6  # the noise is there (1.4e-9 cells)
        result = plan(tiles)
        placed = {p.name: p.canvas for p in result.tiles}
        assert placed["ne.tif"] == window(mz, 0, 6, 5, 7)
        assert placed["se.tif"] == window(mz, 4, 6, 5, 7)
        assert result.meta.x_min == 600_000.0

    @pytest.mark.parametrize(("noise", "accepted"), [(5e-7, True), (2e-6, False)])
    def test_the_tolerance_boundary(
        self, plan: Any, noise: float, accepted: bool, refused: Any
    ) -> None:
        source = whole()
        tiles = {
            "nw.tif": piece(source, 0, 5, 0, 7),
            "ne.tif": piece(source, 0, 5, 6, 13, x_min=X0 + (6 + noise) * DX),
        }
        if accepted:
            assert len(plan(tiles).tiles) == 2
        else:
            refused(lambda: plan(tiles), "nw.tif", "ne.tif")


class TestM4Bounds:
    """M4 and R4 point 4: the box snapped outward to the lattice, clamped."""

    @pytest.fixture
    def tiles(self) -> dict[str, DemTile]:
        return quadrants(whole(), row_cut=4, col_cut=6, overlap=1)

    @pytest.mark.parametrize(
        ("box", "expected"),
        [
            # cols floor(1.2)=1 .. ceil(3.3)=4; rows floor(0.8)=0 .. ceil(2.6)=3
            ((500012.0, 6599987.0, 500033.0, 6599996.0), (0, 1, 4, 4)),
            # inside one cell: 2 x 2 nodes
            ((500012.0, 6599991.0, 500018.0, 6599994.0), (1, 1, 2, 2)),
            # on nodes exactly: floor and ceil of integers add nothing
            ((500010.0, 6599985.0, 500030.0, 6599995.0), (1, 1, 3, 3)),
            # past the west and north edges: clamped
            ((499000.0, 6599991.0, 500025.0, 6600500.0), (0, 0, 3, 4)),
            # at the east edge, two columns remain
            ((500115.0, 6599981.0, 500200.0, 6599984.0), (3, 11, 2, 2)),
        ],
        ids=["outward", "inside-one-cell", "on-nodes", "clamped", "east-edge"],
    )
    def test_the_window(
        self,
        mz: ModuleType,
        build: Any,
        tiles: dict[str, DemTile],
        box: tuple[float, ...],
        expected: tuple[int, ...],
    ) -> None:
        row0, col0, rows, cols = expected
        result, _ = build(tiles, box)
        assert result.plan.window == window(mz, *expected)
        assert result.plan.meta.x_min == X0 + col0 * DX
        assert result.plan.meta.y_max == Y0 - row0 * DY
        assert (result.plan.meta.rows, result.plan.meta.cols) == (rows, cols)
        assert same_array(result.tile.array, whole().array[row0 : row0 + rows, col0 : col0 + cols])

    def test_a_tile_touching_the_box_on_one_node_line_contributes_it(
        self, mz: ModuleType, build: Any
    ) -> None:
        source = whole(rows=8, cols=12, area=True)
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=0)
        # x_max lands exactly on ne's first column (500060); y on the top rows.
        result, loads = build(tiles, (500012.0, 6599991.0, 500060.0, 6599999.0))
        placed = {p.name: p for p in result.plan.tiles}
        assert set(placed) == {"nw.tif", "ne.tif"}
        assert placed["ne.tif"].source == window(mz, 0, 0, 3, 1)
        assert placed["ne.tif"].canvas == window(mz, 0, 5, 3, 1)
        assert same_array(result.tile.array, source.array[0:3, 1:7])
        assert sorted(loads.calls) == ["ne.tif", "nw.tif"]

    def test_a_box_meeting_no_tile_is_refused(
        self, plan: Any, refused: Any, tiles: dict[str, DemTile]
    ) -> None:
        refused(lambda: plan(tiles, (400000.0, 6000000.0, 400100.0, 6000100.0)))

    def test_a_box_leaving_one_column_at_the_edge_is_refused(
        self, plan: Any, refused: Any, tiles: dict[str, DemTile]
    ) -> None:
        refused(lambda: plan(tiles, (500120.0, 6599980.0, 500200.0, 6599990.0)))


class TestM5MemoryCap:
    """M5, changed for 15a (R7): the cap is half of physical memory, from a function."""

    NODES = 9 * 13

    def test_physical_memory_is_the_machines(self, mz: ModuleType) -> None:
        import os

        pages = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        assert mz.physical_memory() == pages

    @pytest.mark.parametrize(
        ("memory", "accepted"), [(2 * NODES * 4, True), (2 * NODES * 4 - 1, False)]
    )
    def test_the_boundary_at_four_bytes_a_node(
        self,
        mz: ModuleType,
        plan: Any,
        refused: Any,
        monkeypatch: pytest.MonkeyPatch,
        memory: int,
        accepted: bool,
    ) -> None:
        monkeypatch.setattr(mz, "physical_memory", lambda: memory)
        tiles = quadrants(whole(), row_cut=4, col_cut=6, overlap=1)
        if accepted:
            assert plan(tiles).meta.rows == 9
        else:
            refused(lambda: plan(tiles), "--bbox")

    def test_the_cap_is_on_the_window_not_the_union(
        self, mz: ModuleType, plan: Any, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(mz, "physical_memory", lambda: 2 * 16 * 4)
        tiles = quadrants(whole(), row_cut=4, col_cut=6, overlap=1)
        assert plan(tiles, (500012.0, 6599987.0, 500033.0, 6599996.0)).meta.cols == 4


class TestM6PlanOrderIndependence:
    """M6 and I2: every permutation of the footprints gives the same plan."""

    def test_every_permutation(self, mz: ModuleType, footprint: Any) -> None:
        tiles = quadrants(whole(), row_cut=4, col_cut=6, overlap=1)
        prints = footprints(footprint, tiles)
        first = mz.plan_mosaic(prints, None, None)
        for order in itertools.permutations(prints):
            assert mz.plan_mosaic(list(order), None, None) == first


class TestM13ReferenceNode:
    """M13 and R4 point 2: global indices count from the lattice's reference node.

    Spacing 0.3 x 0.7 and an origin picked (by search) so that stepping from
    `se`'s own origin, `(XR + 6 dx) + dx`, differs in the last bit from
    `XR + 7 dx`: a formula other than `X_ref + c0 * dx` fails bit-for-bit.
    """

    XR, YR, SX, SY = 3448.3132, 515.2302, 0.3, 0.7  # searched: se's own origin + 1 cell differs

    @pytest.fixture
    def tiles(self) -> dict[str, DemTile]:
        source = whole(x_min=self.XR, y_max=self.YR, dx=self.SX, dy=self.SY)
        return quadrants(source, row_cut=4, col_cut=6, overlap=1)

    def box(self, c0: float, r0: float, c1: float, r1: float) -> tuple[float, ...]:
        return (
            self.XR + c0 * self.SX,
            self.YR - r1 * self.SY,
            self.XR + c1 * self.SX,
            self.YR - r0 * self.SY,
        )

    def test_the_reference_is_the_north_west_node_of_the_whole_lattice(
        self, plan: Any, tiles: dict[str, DemTile]
    ) -> None:
        only_se = plan(tiles, self.box(7.5, 5.5, 10.2, 7.3))
        assert [p.name for p in only_se.tiles] == ["se.tif"]
        assert only_se.reference == (self.XR, self.YR)

    @pytest.mark.parametrize(
        "box", [None, (7.5, 5.5, 10.2, 7.3), (0.5, 0.5, 12.0, 8.0), (5.5, 3.5, 6.5, 4.5)],
        ids=["no-bounds", "se-alone", "most", "the-shared-corner"],
    )  # fmt: skip
    def test_the_origin_is_reference_plus_index_times_spacing_bit_for_bit(
        self, plan: Any, tiles: dict[str, DemTile], box: tuple[float, ...] | None
    ) -> None:
        result = plan(tiles, None if box is None else self.box(*box))
        assert result.meta.x_min == result.reference[0] + result.window.col0 * self.SX
        assert result.meta.y_max == result.reference[1] - result.window.row0 * self.SY

    def test_a_tiles_global_index_does_not_depend_on_the_request(
        self, plan: Any, tiles: dict[str, DemTile]
    ) -> None:
        requests = [
            None,
            self.box(7.5, 5.5, 10.2, 7.3),
            self.box(0.5, 0.5, 12.0, 8.0),
            self.box(5.5, 3.5, 6.5, 4.5),
        ]
        for box in requests:
            assert global_origin(plan(tiles, box), "se.tif") == (4, 6), box

    def test_the_mosaic_meta_takes_the_combined_provenance(self, plan: Any) -> None:
        """R4: `nodata_source` from the first tile by name; `vertical_unit_assumed` if any."""
        source = whole(nodata=SENTINEL, vertical_unit_assumed=False)
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        tiles["ne.tif"] = piece(source, 0, 5, 6, 13, nodata_source="caller")
        tiles["sw.tif"] = piece(source, 4, 9, 0, 7, vertical_unit_assumed=True)
        result = plan(tiles)
        assert result.meta.nodata == SENTINEL
        assert result.meta.nodata_source == "caller"
        assert result.meta.vertical_unit_assumed is True


class TestM14Coverage:
    """M14 and R4 point 5 (Ola's ruling in 18 R6): a needed node no tile covers.

    Area-registered 3 x 4-node blocks of a 6 x 8 grid, with block (1, 1)
    missing: its nodes are x 500040..500070, y 6599975..6599985.
    """

    @pytest.fixture
    def source(self) -> DemTile:
        return whole(rows=6, cols=8, area=True)

    @pytest.fixture
    def holed(self, source: DemTile) -> dict[str, DemTile]:
        return blocks(source, 3, 4, skip={(1, 1)})

    def test_a_hole_inside_the_window_is_refused_naming_it(
        self, plan: Any, refused: Any, holed: dict[str, DemTile]
    ) -> None:
        refused(lambda: plan(holed), "500040", "500070", "6599975", "6599985")

    def test_a_hole_in_the_middle_is_refused(self, plan: Any, refused: Any) -> None:
        """No corner of the window is uncovered: a check of the corners, or of
        the union's bounding box, passes this; only a check of every node fails it."""
        tiles = blocks(whole(rows=9, cols=12, area=True), 3, 4, skip={(1, 1)})
        refused(lambda: plan(tiles), "500040", "500070")

    def test_a_request_away_from_the_hole_plans(self, plan: Any, holed: dict[str, DemTile]) -> None:
        result = plan(holed, (500001.0, 6599991.0, 500069.0, 6599999.0))
        assert {p.name for p in result.tiles} == {"b00.tif", "b01.tif"}

    def test_a_request_reaching_one_hole_node_is_refused(
        self, plan: Any, refused: Any, holed: dict[str, DemTile]
    ) -> None:
        # rows 0..3 (y 6600000..6599985), cols 0..4 (x 500000..500040): node
        # (500040, 6599985) is the hole's north-west node.
        refused(
            lambda: plan(holed, (500000.0, 6599985.0, 500040.0, 6600000.0)), "500040", "6599985"
        )

    def test_the_hole_outside_the_needed_region_is_nan_filler(
        self, build: Any, source: DemTile, holed: dict[str, DemTile]
    ) -> None:
        needed = shapely.box(500000.0, 6599990.0, 500070.0, 6600000.0).union(
            shapely.box(500000.0, 6599975.0, 500030.0, 6599990.0)
        )
        result, _ = build(holed, None, needed)
        array = result.tile.array
        assert np.isnan(array[3:, 4:]).all()
        mask = np.ones(array.shape, dtype=bool)
        mask[3:, 4:] = False
        assert np.array_equal(array[mask], source.array[mask])

    def test_a_needed_region_touching_one_hole_node_is_refused(
        self, plan: Any, refused: Any, holed: dict[str, DemTile]
    ) -> None:
        needed = shapely.box(500000.0, 6599990.0, 500070.0, 6600000.0).union(
            shapely.box(500000.0, 6599975.0, 500040.0, 6599990.0)
        )  # reaches x 500040: the hole's west column, on the boundary
        refused(lambda: plan(holed, None, needed), "500040")


# ---------------------------------------------------------------------------
# Assembly (pixels, no files)
# ---------------------------------------------------------------------------


class TestM7Values:
    """M7: values at the right nodes; the result's array contract."""

    def test_values_land_at_their_nodes(self, build: Any) -> None:
        source = whole()
        result, _ = build(quadrants(source, row_cut=4, col_cut=6, overlap=1))
        assert result.tile.meta == source.meta
        assert same_array(result.tile.array, source.array)

    def test_the_array_is_read_only_and_c_contiguous(self, build: Any) -> None:
        result, _ = build(quadrants(whole(), row_cut=4, col_cut=6, overlap=1))
        array = result.tile.array
        assert not array.flags.writeable
        assert array.flags.c_contiguous
        with pytest.raises(ValueError):
            array[0, 0] = 1.0

    def test_float64_when_any_tile_is(self, build: Any) -> None:
        source = whole()
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        tiles["se.tif"] = piece(
            source, 4, 9, 6, 13, array=source.array[4:9, 6:13].astype(np.float64)
        )
        result, _ = build(tiles)
        assert result.tile.array.dtype == np.float64
        assert same_array(result.tile.array, source.array.astype(np.float64))

    def test_float32_stays_float32(self, build: Any) -> None:
        result, _ = build(quadrants(whole(), row_cut=4, col_cut=6, overlap=1))
        assert result.tile.array.dtype == np.float32


class TestM8Overlaps:
    """M8, R5 and I3: valid beats NoData; disagreement refuses (Q1).

    `w.tif` and `e.tif` cut from a 4 x 9 grid, sharing `overlap` columns.
    """

    ROWS, COLS = 4, 9

    @staticmethod
    def pair(overlap: int, nodata: float | None = None) -> tuple[DemTile, dict[str, np.ndarray]]:
        source = whole(rows=TestM8Overlaps.ROWS, cols=TestM8Overlaps.COLS, nodata=nodata)
        cut = 4
        return source, {
            "w.tif": source.array[:, : cut + overlap].copy(),
            "e.tif": source.array[:, cut:].copy(),
        }

    @staticmethod
    def tiles(source: DemTile, arrays: dict[str, np.ndarray]) -> dict[str, DemTile]:
        w, e = arrays["w.tif"], arrays["e.tif"]
        return {
            "w.tif": piece(source, 0, w.shape[0], 0, w.shape[1], array=w),
            "e.tif": piece(source, 0, e.shape[0], 4, 4 + e.shape[1], array=e),
        }

    @pytest.mark.parametrize("overlap", [1, 3])
    @pytest.mark.parametrize("holder", ["w.tif", "e.tif"])
    def test_valid_beats_nan_in_either_tile(self, build: Any, overlap: int, holder: str) -> None:
        source, arrays = self.pair(overlap)
        column = 4 if holder == "e.tif" else 4 + overlap - 1
        local = column - (4 if holder == "e.tif" else 0)
        arrays[holder][1, local] = np.nan
        result, _ = build(self.tiles(source, arrays))
        assert same_array(result.tile.array, source.array)

    @pytest.mark.parametrize("holder", ["w.tif", "e.tif"])
    def test_valid_beats_the_sentinel_in_either_tile(self, build: Any, holder: str) -> None:
        source, arrays = self.pair(3, nodata=SENTINEL)
        arrays[holder][2, 5 if holder == "w.tif" else 1] = SENTINEL  # canvas column 5
        result, _ = build(self.tiles(source, arrays))
        assert same_array(result.tile.array, source.array)

    def test_both_nan_stays_nan_and_both_sentinel_stays_sentinel(self, build: Any) -> None:
        source, arrays = self.pair(3, nodata=SENTINEL)
        arrays["w.tif"][0, 5] = arrays["e.tif"][0, 1] = np.nan
        arrays["w.tif"][3, 6] = arrays["e.tif"][3, 2] = SENTINEL
        result, _ = build(self.tiles(source, arrays))
        array = result.tile.array
        assert np.isnan(array[0, 5])
        assert array[3, 6] == SENTINEL
        expected = source.array.copy()
        expected[0, 5], expected[3, 6] = np.nan, SENTINEL
        assert same_array(array, expected)

    def test_nan_against_the_sentinel_is_nodata_in_either_order(
        self, mz: ModuleType, plan: Any
    ) -> None:
        """Which NoData wins is not ruled; that it is NoData, and the same in
        every order, is (I2, I3)."""
        source, arrays = self.pair(3, nodata=SENTINEL)
        arrays["w.tif"][1, 4] = np.nan
        arrays["e.tif"][1, 0] = SENTINEL
        tiles = self.tiles(source, arrays)
        planned = plan(tiles)
        results = [
            mz.assemble(planned.model_copy(update={"tiles": order}), Loads(tiles)).tile.array
            for order in itertools.permutations(planned.tiles)
        ]
        value = results[0][1, 4]
        assert np.isnan(value) or value == SENTINEL
        assert all(same_array(r, results[0]) for r in results)

    def test_equal_overlaps_are_accepted(self, build: Any) -> None:
        source, arrays = self.pair(3)
        result, _ = build(self.tiles(source, arrays))
        assert same_array(result.tile.array, source.array)

    @pytest.mark.parametrize("order", ["as-named", "reversed"])
    def test_a_disagreement_is_refused_naming_both_the_count_and_the_largest(
        self, mz: ModuleType, plan: Any, refused: Any, order: str
    ) -> None:
        source, arrays = self.pair(3)
        arrays["e.tif"][0, 1] += 0.5
        arrays["e.tif"][2, 1] += 2.5  # the middle column of the overlap
        arrays["e.tif"][3, 2] -= 1.0
        arrays["w.tif"][1, 4] = np.nan  # a NoData node in the overlap is not counted
        tiles = self.tiles(source, arrays)
        planned = plan(tiles)
        if order == "reversed":
            planned = planned.model_copy(update={"tiles": tuple(reversed(planned.tiles))})
        message = refused(lambda: mz.assemble(planned, Loads(tiles)), "w.tif", "e.tif", "2.5")
        names_number(message, 3)

    def test_one_ulp_is_a_disagreement(self, build: Any, refused: Any) -> None:
        """`==`, not `isclose`: DTM10's and ANADEM's overlaps agree bit for bit (N3, B4)."""
        source, arrays = self.pair(1)
        value = arrays["e.tif"][2, 0]
        arrays["e.tif"][2, 0] = np.nextafter(value, np.float32(np.inf), dtype=np.float32)
        message = refused(lambda: build(self.tiles(source, arrays)), "w.tif", "e.tif")
        names_number(message, 1)

    def test_a_disagreement_against_a_sentinel_tile_is_still_one(
        self, build: Any, refused: Any
    ) -> None:
        """The sentinel is NoData only where a node holds it; elsewhere a value is a value."""
        source, arrays = self.pair(3, nodata=SENTINEL)
        arrays["w.tif"][0, 4] = SENTINEL
        arrays["w.tif"][0, 5] += 7.0
        message = refused(lambda: build(self.tiles(source, arrays)), "w.tif", "e.tif", "7")
        names_number(message, 1)


class TestM9SplitAndRestitch:
    """M9 and I5: a grid cut into tiles reassembles to itself, meta and array."""

    @pytest.mark.parametrize(
        ("area", "overlap", "shape"),
        [(False, 1, (9, 13)), (True, 0, (8, 12)), (True, 3, (8, 12)), (False, 3, (9, 13))],
        ids=["point-shared-line", "area-abutting", "area-overlap-3", "point-overlap-3"],
    )
    def test_quadrants(self, build: Any, area: bool, overlap: int, shape: tuple[int, int]) -> None:
        source = whole(*shape, area=area)
        result, _ = build(quadrants(source, row_cut=4, col_cut=6, overlap=overlap))
        assert result.tile.meta == source.meta
        assert same_array(result.tile.array, source.array)

    def test_with_nodata_inside_and_across_the_seams(self, build: Any) -> None:
        grid = values(9, 13)
        grid[0, 0] = grid[4, 6] = grid[5, 7] = grid[8, 12] = SENTINEL  # (4, 6) is in all four
        grid[4, 2] = grid[6, 6] = grid[2, 9] = np.nan  # on the seams
        source = whole(9, 13, array=grid, nodata=SENTINEL)
        result, _ = build(quadrants(source, row_cut=4, col_cut=6, overlap=3))
        assert result.tile.meta == source.meta
        assert same_array(result.tile.array, source.array)

    def test_a_three_by_three_block_split(self, build: Any) -> None:
        source = whole(9, 12, area=True)
        result, _ = build(blocks(source, 3, 4))
        assert result.tile.meta == source.meta
        assert same_array(result.tile.array, source.array)


class TestM10AssemblyOrderIndependence:
    """M10 and I2: the assembled array does not depend on tile order."""

    def test_every_order_of_the_plans_tiles(self, mz: ModuleType, plan: Any) -> None:
        grid = values(9, 13)
        grid[5, 7] = np.nan
        grid[4, 6] = SENTINEL
        source = whole(9, 13, array=grid, nodata=SENTINEL)
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=3)
        planned = plan(tiles)
        for order in itertools.permutations(planned.tiles):
            result = mz.assemble(planned.model_copy(update={"tiles": order}), Loads(tiles))
            assert same_array(result.tile.array, source.array), [p.name for p in order]

    def test_tiles_are_loaded_once_each_in_plan_order(self, build: Any) -> None:
        result, loads = build(quadrants(whole(), row_cut=4, col_cut=6, overlap=1))
        assert loads.calls == [p.name for p in result.plan.tiles] == sorted(loads.calls)


class TestM11SingleTile:
    """M11 and I7: one tile, no bounds: the loaded tile itself, no canvas."""

    def test_the_loaded_tile_is_returned_as_is(self, build: Any) -> None:
        source = whole()
        result, loads = build({"only.tif": source})
        assert result.tile is loads.tiles["only.tif"]

    def test_with_bounds_it_is_a_window(self, build: Any) -> None:
        source = whole()
        result, _ = build({"only.tif": source}, (500012.0, 6599987.0, 500033.0, 6599996.0))
        assert same_array(result.tile.array, source.array[0:4, 1:5])


class TestM12Loads:
    """M12: only the selected tiles are loaded; a changed tile is refused."""

    def test_only_selected_tiles_are_loaded(self, build: Any) -> None:
        _, loads = build(
            quadrants(whole(), row_cut=4, col_cut=6, overlap=1),
            (500012.0, 6599987.0, 500033.0, 6599996.0),
        )
        assert loads.calls == ["nw.tif"]

    def test_a_tile_changed_since_it_was_listed(
        self, mz: ModuleType, plan: Any, refused: Any
    ) -> None:
        source = whole()
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        planned = plan(tiles)
        moved = dict(tiles)
        moved["nw.tif"] = piece(source, 0, 5, 0, 7, y_max=Y0 + DY)
        refused(lambda: mz.assemble(planned, Loads(moved)), "nw.tif", "changed since it was listed")


class TestM15Adopt:
    """M15 and R7: one canvas, handed to `DemTile` without a copy."""

    def test_the_canvas_is_adopted_not_copied(
        self, build: Any, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        adopted: list[np.ndarray] = []
        original = DemTile._adopt

        def spy(meta: Any, array: np.ndarray) -> DemTile:
            adopted.append(array)
            return original(meta, array)

        monkeypatch.setattr(DemTile, "_adopt", staticmethod(spy))
        result, _ = build(quadrants(whole(), row_cut=4, col_cut=6, overlap=1))
        assert len(adopted) == 1
        assert np.shares_memory(adopted[0], result.tile.array)
        assert not result.tile.array.flags.writeable

    def test_adopt_keeps_the_buffer_and_makes_it_read_only(self) -> None:
        buffer = values(3, 4)
        tile = DemTile._adopt(meta(rows=3, cols=4), buffer)
        assert np.shares_memory(tile.array, buffer)
        assert not tile.array.flags.writeable
        assert not buffer.flags.writeable
        assert tile.meta == meta(rows=3, cols=4)

    @pytest.mark.parametrize(
        "array",
        [
            values(4, 3),
            values(3, 4, np.int32),
            values(3, 4).reshape(1, 3, 4),
            np.asfortranarray(values(3, 4)),
        ],
        ids=["shape", "dtype", "ndim", "fortran-order"],
    )
    def test_adopt_runs_the_same_checks(self, array: np.ndarray) -> None:
        with pytest.raises(ValueError):
            DemTile._adopt(meta(rows=3, cols=4), array)

    def test_the_public_constructor_still_copies(self) -> None:
        buffer = values(3, 4)
        tile = DemTile(meta=meta(rows=3, cols=4), array=buffer)
        assert not np.shares_memory(tile.array, buffer)
        assert buffer.flags.writeable

    def test_adopt_is_called_only_from_the_mosaic(self) -> None:
        """R7: "Its one caller is `assemble`". Grep the package for the name."""
        users = sorted(
            str(path.relative_to(SRC))
            for path in SRC.rglob("*.py")
            if "_adopt" in path.read_text(encoding="utf-8")
        )
        assert users == ["io/models.py", "mosaic.py"]


def test_the_plan_never_loads(mz: ModuleType, footprint: Any) -> None:
    """I6 as a type fact: `plan_mosaic` takes footprints, and no `load`."""
    import inspect

    assert list(inspect.signature(mz.plan_mosaic).parameters) == ["footprints", "bounds", "needed"]
