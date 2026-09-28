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

AMENDED FOR OLA'S Q1 REVISED (2026-09-28). Overlaps that disagree are no
longer refused: each node takes the value of the tile it lies deepest in, ties
by name, and each disagreeing pair is reported in `Mosaic.seams`
(`TestQ1DeepestInterior`, `TestQ1SeamReport`, and four tests of
`TestM8Overlaps` that pinned the refusal and now pin the report).

AMENDED FOR OLA'S 1 MM THRESHOLD (2026-09-28, "Ignore below 1mm"). A seam
counts only nodes where |a - b| >= 0.001; the midline rule does not look at
the threshold (`TestSeamThreshold`; `TestM8Overlaps`' one-ulp test inverted).

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
from typing import Any, ClassVar

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
    deepest_interior,
    footprints,
    meta,
    piece,
    quadrants,
    same_array,
    seams,
    seams_of,
    shifted_by,
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
        """Unchanged by Ola's Q5 reading: without a box the request is every
        tile, and neither lattice covers the other's nodes (odd.tif reaches
        x 500165, past the main lattice's 500120), so none covers it (R6)."""
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


class TestB1LatticeByCoverage:
    """Ola's reading of Q5 after review ("Ruled by Ola", 2026-09-27): per
    lattice, whether **its own tiles cover every node the request needs**.

    Exactly one covers: plan on it, and the other lattices' tiles are not in
    the plan. None covers: the Q5 refusal (the straddle tests in M2). Several
    cover (a box wholly inside an overlap strip): the lattice with the most
    tiles in the repository, ties by the name of its first tile.

    The repository is M2's: the 9 x 13 main lattice in four, and `odd.tif`
    (x 500105..500165, y 6599970..6599950) half a cell east-west off it. Their
    overlap strip is x 500105..500120, y 6599970..6599960.

    The odd tile is also named `a_odd.tif`, which sorts before the main
    lattice's `ne.tif`: an implementation that takes whichever lattice comes
    first, whatever it covers, fails one of the two names in each case.
    """

    MAIN_REACHING_INTO_THE_STRIP = (500072.0, 6599962.0, 500112.0, 6599974.0)
    ODD_REACHING_INTO_THE_STRIP = (500106.0, 6599952.0, 500150.0, 6599969.0)
    INSIDE_THE_STRIP = (500106.0, 6599961.0, 500114.0, 6599969.0)
    ODD_ORIGIN = (X0 + 10.5 * DX, Y0 - 6 * DY)

    @staticmethod
    def main() -> dict[str, DemTile]:
        return quadrants(whole(), row_cut=4, col_cut=6, overlap=1)

    @staticmethod
    def odd(name: str = "odd.tif", pieces: int = 1) -> dict[str, DemTile]:
        """The odd lattice's 5 x 7 tile, whole, or in four named `<name>-<q>`."""
        x_min, y_max = TestB1LatticeByCoverage.ODD_ORIGIN
        tile = whole(rows=5, cols=7, x_min=x_min, y_max=y_max)
        if pieces == 1:
            return {name: tile}
        return {
            f"{name}-{q}": t for q, t in quadrants(tile, row_cut=2, col_cut=3, overlap=1).items()
        }

    @pytest.mark.parametrize("odd_name", ["odd.tif", "a_odd.tif"])
    def test_a_main_box_reaching_into_the_strip_plans_on_the_main_lattice(
        self, mz: ModuleType, plan: Any, build: Any, odd_name: str
    ) -> None:
        """The reviewer's case: the box selects the odd tile, and today is
        refused as mixed-lattice, but the main lattice covers every node."""
        main = self.main()
        tiles = {**main, **self.odd(odd_name)}
        result = plan(tiles, self.MAIN_REACHING_INTO_THE_STRIP)
        assert [p.name for p in result.tiles] == ["se.tif"]
        assert result.reference == (X0, Y0)
        assert result.window == window(mz, 5, 7, 4, 6)
        assert result == plan(main, self.MAIN_REACHING_INTO_THE_STRIP)
        mosaic, loads = build(tiles, self.MAIN_REACHING_INTO_THE_STRIP)
        assert loads.calls == ["se.tif"]
        assert same_array(mosaic.tile.array, whole().array[5:9, 7:13])

    @pytest.mark.parametrize("odd_name", ["odd.tif", "a_odd.tif"])
    def test_an_odd_box_reaching_into_the_strip_plans_on_the_odd_lattice(
        self, mz: ModuleType, plan: Any, odd_name: str
    ) -> None:
        """The mirror case. The main lattice's tiles cover the part of the box
        over them, but not x 500130..500150 or y 6599955..6599950, which the
        odd tile holds: the box is not clamped to the main lattice's union."""
        odd = self.odd(odd_name)
        result = plan({**self.main(), **odd}, self.ODD_REACHING_INTO_THE_STRIP)
        assert [p.name for p in result.tiles] == [odd_name]
        assert result.reference == self.ODD_ORIGIN
        assert result.window == window(mz, 0, 0, 5, 6)
        assert result == plan(odd, self.ODD_REACHING_INTO_THE_STRIP)

    def test_several_cover_the_most_tiles_in_the_repository_wins(
        self, mz: ModuleType, plan: Any
    ) -> None:
        """Main 4 tiles, odd 1; each selects exactly one tile for this box, so
        a count of *selected* tiles would tie."""
        main = self.main()
        result = plan({**main, **self.odd()}, self.INSIDE_THE_STRIP)
        assert [p.name for p in result.tiles] == ["se.tif"]
        assert result.reference == (X0, Y0)
        assert result == plan(main, self.INSIDE_THE_STRIP)

    def test_the_count_is_of_the_repository_not_of_the_selection(
        self, footprint: Any, mz: ModuleType
    ) -> None:
        """Main 1 tile (se alone), odd 2 (one far away, never selected); both
        select one. The odd names sort after `se.tif`, so a count of the
        selection, tied and broken by name, would choose the main lattice."""
        x_min, y_max = self.ODD_ORIGIN
        tiles = {
            "se.tif": self.main()["se.tif"],
            "x_odd.tif": self.odd()["odd.tif"],
            "x_far.tif": whole(rows=3, cols=3, x_min=x_min + 100 * DX, y_max=y_max),
        }
        result = mz.plan_mosaic(
            footprints(footprint, tiles), bounds_of(mz, self.INSIDE_THE_STRIP), None
        )
        assert [p.name for p in result.tiles] == ["x_odd.tif"]
        assert result.reference == (x_min, y_max)

    @pytest.mark.parametrize(
        ("prefix", "winner"),
        [("a", "odd"), ("o", "main")],
        ids=["odd-first-name-sorts-first", "main-first-name-sorts-first"],
    )
    def test_a_tie_goes_to_the_lattice_whose_first_tile_sorts_first(
        self, plan: Any, prefix: str, winner: str
    ) -> None:
        """Four tiles each. Main's first tile is `ne.tif` (last `sw.tif`); the
        odd lattice's first is `<prefix>-ne.tif` (last `<prefix>-sw.tif`).
        With prefix `o`, main's first name sorts first but the odd lattice's
        last name does: a tie broken by the last tile, or by the larger first
        name, picks the wrong one in one of the two cases."""
        odd = self.odd(prefix, pieces=4)
        result = plan({**self.main(), **odd}, self.INSIDE_THE_STRIP)
        names = {p.name for p in result.tiles}
        if winner == "odd":
            assert names <= set(odd), names
            assert result.reference == self.ODD_ORIGIN
        else:
            assert names == {"se.tif"}
            assert result.reference == (X0, Y0)

    @pytest.mark.parametrize(
        "odd_pieces", [1, 4], ids=["most-tiles", "tie-by-name"]
    )  # fmt: skip
    def test_the_choice_does_not_depend_on_footprint_order(
        self, mz: ModuleType, footprint: Any, odd_pieces: int
    ) -> None:
        tiles = {**self.main(), **self.odd("a", pieces=odd_pieces)}
        prints = footprints(footprint, tiles)
        box = bounds_of(mz, self.INSIDE_THE_STRIP)
        first = mz.plan_mosaic(prints, box, None)
        for k in range(len(prints)):  # every rotation, forwards and backwards
            rotated = prints[k:] + prints[:k]
            for order in (rotated, rotated[::-1]):
                assert mz.plan_mosaic(order, box, None) == first


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

    def test_float_noise_in_a_decimetre_spacing_does_not_add_a_node_line(
        self, mz: ModuleType, plan: Any
    ) -> None:
        """S3, R4 point 4's snap. At spacing 0.1 from a reference at (0, 1):
        0.3 / 0.1 is 2.9999999999999996 (floor 2), (1 - 0.9) / 0.1 is
        0.9999999999999998 (floor 0) and (1 - 0.7) / 0.1 is 3.0000000000000004
        (ceil 4). Each edge is a node, so each would add a line without the snap."""
        assert (0.3 / 0.1, (1.0 - 0.9) / 0.1, (1.0 - 0.7) / 0.1) == (
            2.9999999999999996,
            0.9999999999999998,
            3.0000000000000004,
        )  # the noise is there
        tiles = {"d.tif": whole(x_min=0.0, y_max=1.0, dx=0.1, dy=0.1)}
        result = plan(tiles, (0.3, 0.7, 0.6, 0.9))
        assert result.window == window(mz, 1, 3, 3, 4)

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

    @pytest.mark.parametrize(
        ("memory", "accepted"), [(2 * NODES * 8, True), (2 * NODES * 8 - 1, False)]
    )
    def test_the_boundary_at_eight_bytes_a_node_when_a_selected_tile_decodes_to_float64(
        self,
        mz: ModuleType,
        footprint: Any,
        refused: Any,
        monkeypatch: pytest.MonkeyPatch,
        memory: int,
        accepted: bool,
    ) -> None:
        """S2: the cap counts the planned decoded dtype, `np.result_type` of the
        selected footprints' `dtype`. One float64 tile makes the canvas float64
        (R5), so the window costs 8 bytes a node, not 4."""
        monkeypatch.setattr(mz, "physical_memory", lambda: memory)
        tiles = quadrants(whole(), row_cut=4, col_cut=6, overlap=1)
        prints = [
            footprint(name=name, meta=tile.meta, dtype=np.dtype(np.float64))
            if name == "se.tif"
            else footprint(name=name, meta=tile.meta)
            for name, tile in tiles.items()
        ]
        call = lambda: mz.plan_mosaic(prints, None, None)  # noqa: E731
        if accepted:
            assert call().meta.rows == 9
        else:
            refused(call, "--bbox", "float64")

    def test_a_float64_tile_outside_the_window_does_not_count(
        self, mz: ModuleType, footprint: Any, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Only the selected tiles decode, so only their dtypes size the canvas."""
        monkeypatch.setattr(mz, "physical_memory", lambda: 2 * 16 * 4)
        tiles = quadrants(whole(), row_cut=4, col_cut=6, overlap=1)
        prints = [
            footprint(name=name, meta=tile.meta, dtype=np.dtype(np.float64))
            if name == "se.tif"
            else footprint(name=name, meta=tile.meta)
            for name, tile in tiles.items()
        ]
        box = bounds_of(mz, (500012.0, 6599987.0, 500033.0, 6599996.0))  # 4 x 4, nw alone
        assert [p.name for p in mz.plan_mosaic(prints, box, None).tiles] == ["nw.tif"]

    def test_a_footprints_dtype_defaults_to_float32(self, footprint: Any) -> None:
        """S2: `TileFootprint` gains an optional decoded `dtype`, float32 unless given."""
        tile = whole()
        assert np.dtype(footprint(name="a.tif", meta=tile.meta).dtype) == np.float32
        given = footprint(name="a.tif", meta=tile.meta, dtype=np.dtype(np.float64))
        assert np.dtype(given.dtype) == np.float64

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

    def test_refusing_a_mostly_uncovered_window_costs_less_than_its_canvas(
        self, plan: Any, refused: Any
    ) -> None:
        """S1: the check runs on windows the cap allows, up to half of physical
        memory at 4 bytes a node, so what it allocates per node bounds what the
        cap can promise. Two 2 x 2 tiles at opposite corners of a 2000 x 2000
        window: nearly every node is uncovered. Index and coordinate arrays of
        the uncovered nodes cost 32 bytes a node, 8x the float32 canvas; the
        refusal must peak below the canvas itself (16 MB here). tracemalloc
        sees numpy's buffers, so this is a count of bytes, not a timing."""
        import tracemalloc

        n = 2000
        corners = {
            "nw.tif": whole(rows=2, cols=2),
            "se.tif": whole(rows=2, cols=2, x_min=X0 + (n - 2) * DX, y_max=Y0 - (n - 2) * DY),
        }
        canvas_bytes = n * n * 4
        tracemalloc.start()
        try:
            refused(lambda: plan(corners), f"{X0 + (n - 1) * DX:.0f}", f"{Y0 - (n - 1) * DY:.0f}")
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        assert peak < canvas_bytes, f"peak {peak} bytes for a {canvas_bytes}-byte canvas"

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
    """M8, R5 and I3: valid beats NoData; a disagreement is not refused but
    decided by depth and reported (Ola's Q1 revised, 2026-09-28).

    `w.tif` and `e.tif` cut from a 4 x 9 grid, sharing `overlap` columns.
    With overlap 3 the shared columns are 4, 5 and 6, and each node's own
    depth (`mosaic_fixtures.own_depth`) is, in `w.tif` then `e.tif`:

        rows 0 and 3:  every node is on both tiles' top or bottom border,
                       0 against 0, a tie, so `e.tif` (it sorts first);
        rows 1 and 2:  column 4 is 1 against 0, `w.tif`; column 5 is 1
                       against 1, `e.tif`; column 6 is 0 against 1, `e.tif`.

    The four tests below that pinned the old refusal (Q1 as first ruled) now
    pin the report; each says what it kept.
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

    @pytest.mark.parametrize(
        ("west_wins", "expected"),
        [(False, (3, 2.5, 1.0)), (True, (4, 4.0, 1.75))],
        ids=["odd-count", "even-count"],
    )
    @pytest.mark.parametrize("order", ["as-named", "reversed"])
    def test_a_disagreement_is_reported_with_the_count_the_largest_and_the_median(
        self,
        mz: ModuleType,
        plan: Any,
        order: str,
        west_wins: bool,
        expected: tuple[int, float, float],
    ) -> None:
        """Was `..._is_refused_naming_both_the_count_and_the_largest`. Same
        planted differences; the run is no longer refused. It keeps: both names,
        the count of differing nodes only (a NoData node in the overlap is not
        one, nor are the agreeing ones), the largest |difference|, in either
        order. New: the median, of the differing nodes only, and the value each
        node takes. The median of an even count is the mean of the middle two.
        """
        source, arrays = self.pair(3)
        arrays["e.tif"][0, 1] += 0.5  # row 0, column 5: a tie, e.tif's value
        arrays["e.tif"][2, 1] += 2.5  # row 2, column 5, mid-overlap: 1 against 1, e.tif
        arrays["e.tif"][3, 2] -= 1.0  # row 3, column 6: |difference| is 1
        arrays["w.tif"][1, 4] = np.nan  # a NoData node in the overlap is not counted
        if west_wins:
            arrays["e.tif"][2, 0] += 4.0  # row 2, column 4: w.tif is deeper, its value
        tiles = self.tiles(source, arrays)
        planned = plan(tiles)
        if order == "reversed":
            planned = planned.model_copy(update={"tiles": tuple(reversed(planned.tiles))})
        result = mz.assemble(planned, Loads(tiles))
        assert seams(result) == [("e.tif", "w.tif", *expected)]
        taken = source.array.copy()
        taken[0, 5] += 0.5
        taken[2, 5] += 2.5
        taken[3, 6] -= 1.0
        assert same_array(result.tile.array, taken)

    def test_the_report_names_the_tiles_whose_values_disagree(
        self, mz: ModuleType, plan: Any
    ) -> None:
        """B2, kept in the report (was `test_the_refusal_names_...`). Three
        tiles on one 3 x 3 grid, equal everywhere but node (0, 0), where `a` is
        NoData, `b` 1.0 and `c` 5.0. The seam is b | c; `a` covers the node
        too, but its value is not in the disagreement, so no seam names it.
        Node (0, 0) is on all three borders: a tie, and of the two valid
        values `b.tif`'s, whose name sorts first. In every assembly order."""
        source = whole(rows=3, cols=3)
        tiles = {}
        for name, first in (("a.tif", np.nan), ("b.tif", 1.0), ("c.tif", 5.0)):
            array = source.array.copy()
            array[0, 0] = first
            tiles[name] = piece(source, 0, 3, 0, 3, array=array)
        planned = plan(tiles)
        for order in itertools.permutations(planned.tiles):
            result = mz.assemble(planned.model_copy(update={"tiles": order}), Loads(tiles))
            names = [p.name for p in order]
            assert seams(result) == [("b.tif", "c.tif", 1, 4.0, 4.0)], names
            assert result.tile.array[0, 0] == 1.0, names

    def test_one_ulp_is_below_the_threshold_and_still_decided_by_depth(self, build: Any) -> None:
        """Was `test_one_ulp_is_a_disagreement`, which pinned `!=` (one ulp
        reported). Ola, 2026-09-28: "Ignore below 1mm", so one ulp is not
        reported. What it kept: the node's value is still decided by the
        midline rule. Overlap 1: the shared column is on both tiles' borders,
        a tie, `e.tif`'s value."""
        source, arrays = self.pair(1)
        value = arrays["e.tif"][2, 0]
        bumped = np.nextafter(value, np.float32(np.inf), dtype=np.float32)
        arrays["e.tif"][2, 0] = bumped
        result, _ = build(self.tiles(source, arrays))
        assert result.seams == ()
        assert result.tile.array[2, 4] == bumped

    def test_a_disagreement_against_a_sentinel_tile_is_still_one(self, build: Any) -> None:
        """The sentinel is NoData only where a node holds it; elsewhere a value
        is a value, reported and decided by depth like any other."""
        source, arrays = self.pair(3, nodata=SENTINEL)
        arrays["w.tif"][0, 4] = SENTINEL  # NoData: e.tif's value, not counted
        arrays["w.tif"][0, 5] += 7.0  # a tie on row 0: e.tif's value
        arrays["w.tif"][1, 4] += 3.0  # w.tif is deeper: its value
        result, _ = build(self.tiles(source, arrays))
        assert seams(result) == [("e.tif", "w.tif", 2, 7.0, 5.0)]
        taken = source.array.copy()
        taken[1, 4] += 3.0
        assert same_array(result.tile.array, taken)


def taken_from(result: Any, source: np.ndarray, offsets: dict[str, float]) -> list[str]:
    """Row strings naming, node by node, the tile each value came from.

    Every tile is `source` plus its own offset, so the offset identifies it;
    a tile is written as its name's first letter.
    """
    letter = {offset: name[0] for name, offset in offsets.items()}
    taken = np.asarray(result.tile.array, dtype=np.float64) - source
    return ["".join(letter[float(v)] for v in row) for row in taken]


def offset_tiles(tiles: dict[str, DemTile]) -> tuple[dict[str, DemTile], dict[str, float]]:
    """Each tile plus 1000 times its rank by name: every overlap disagrees."""
    offsets = {name: 1000.0 * (i + 1) for i, name in enumerate(sorted(tiles))}
    return {name: shifted_by(t, offsets[name]) for name, t in tiles.items()}, offsets


class TestQ1DeepestInterior:
    """Ola's Q1 revised (2026-09-28): each overlap node takes the value of the
    tile whose interior it lies deepest in (the largest distance, in nodes, to
    that tile's *own* nearest border); ties go to the name that sorts first;
    NoData still loses to a valid value; the result is independent of tile
    order.

    Every tile is the source plus its own offset, so every overlap disagrees
    and the offset names the tile a value came from. The expected maps are
    written out by hand, and the whole array is also checked against
    `mosaic_fixtures.deepest_interior`, a node-by-node oracle that places
    tiles by coordinates and shares no code with `tin_engine.mosaic`.
    """

    ROWS, WIDTH = 7, 8

    def two(self, overlap: int, west: str, east: str) -> tuple[DemTile, dict[str, DemTile]]:
        """Two 7 x 8 tiles side by side sharing `overlap` columns; west + 1000,
        east + 2000."""
        cols = 2 * self.WIDTH - overlap
        source = whole(self.ROWS, cols)
        return source, {
            west: shifted_by(piece(source, 0, self.ROWS, 0, self.WIDTH), 1000.0),
            east: shifted_by(piece(source, 0, self.ROWS, self.WIDTH - overlap, cols), 2000.0),
        }

    #: The overlap columns only, row by row. W: the west tile is deeper; E: the
    #: east; T: equally deep, so the name that sorts first. Rows 0 and 6 lie on
    #: both tiles' borders. An odd overlap has a tie column down the middle; an
    #: even one splits cleanly where the rows are deep enough (rows 2-4), and
    #: rows 1 and 5, 1 deep, tie the middle two columns.
    MAPS: ClassVar[dict[int, list[str]]] = {
        3: ["TTT", "WTE", "WTE", "WTE", "WTE", "WTE", "TTT"],
        4: ["TTTT", "WTTE", "WWEE", "WWEE", "WWEE", "WTTE", "TTTT"],
    }

    @pytest.mark.parametrize("overlap", [3, 4], ids=["odd-tie-column", "even-midline"])
    @pytest.mark.parametrize(
        ("west", "east"), [("w.tif", "e.tif"), ("a.tif", "b.tif")], ids=["east-first", "west-first"]
    )
    def test_two_tiles_split_down_the_middle_ties_by_name(
        self, build: Any, overlap: int, west: str, east: str
    ) -> None:
        source, tiles = self.two(overlap, west, east)
        result, _ = build(tiles)
        letter = {"W": west[0], "E": east[0], "T": min(west, east)[0]}
        own = self.WIDTH - overlap
        expected = [
            west[0] * own + "".join(letter[k] for k in row) + east[0] * own
            for row in self.MAPS[overlap]
        ]
        assert taken_from(result, source.array, {west: 1000.0, east: 2000.0}) == expected
        oracle = deepest_interior(tiles, result.tile.meta).astype(np.float32)
        assert same_array(result.tile.array, oracle)

    def test_one_shared_line_is_all_ties(self, build: Any) -> None:
        """Point-registered neighbours share one line, 0 deep in both: every
        node is a tie, and the name that sorts first takes the whole line."""
        source, tiles = self.two(1, "w.tif", "e.tif")
        result, _ = build(tiles)
        taken = taken_from(result, source.array, {"w.tif": 1000.0, "e.tif": 2000.0})
        assert [row[self.WIDTH - 1] for row in taken] == ["e"] * self.ROWS

    def test_depth_is_to_the_tiles_own_border_not_the_requests_window(self, build: Any) -> None:
        """A box of rows 2-4 and columns 4-8 cuts both tiles of the odd pair.
        The window's first and last rows are neither tile's border, so the
        answer is rows 2-4 of the full tiles' map (W T E in columns 5-7), not
        the window's (all ties, `e.tif`, on its first and last row)."""
        source, tiles = self.two(3, "w.tif", "e.tif")
        box = (X0 + 4 * DX, Y0 - 4 * DY, X0 + 8 * DX, Y0 - 2 * DY)
        result, _ = build(tiles, box)
        assert (result.tile.meta.rows, result.tile.meta.cols) == (3, 5)
        taken = taken_from(result, source.array[2:5, 4:9], {"w.tif": 1000.0, "e.tif": 2000.0})
        assert taken == ["wweee"] * 3

    #: Rows 2-6, columns 7-11 of a 9 x 12 tile holding a 5 x 5 one flush with
    #: its east border. The small tile's own depth never exceeds the big one's,
    #: so it takes only the ties, and only when its name sorts first.
    INSIDE: ClassVar[dict[str, list[str]]] = {
        "a.tif": ["bbbba", "bbbaa", "bbaaa", "bbbaa", "bbbba"],
        "s.tif": ["bbbbb"] * 5,
    }

    @pytest.mark.parametrize("small", ["a.tif", "s.tif"], ids=["small-first", "small-last"])
    def test_unequal_sizes_measure_depth_to_each_tiles_own_border(
        self, build: Any, small: str
    ) -> None:
        """Depth is to each tile's own border: measured to the mosaic's border
        (here the big tile's), both would tie at every node."""
        source = whole(9, 12)
        tiles = {
            "b.tif": shifted_by(piece(source, 0, 9, 0, 12), 1000.0),
            small: shifted_by(piece(source, 2, 7, 7, 12), 2000.0),
        }
        result, _ = build(tiles)
        taken = taken_from(result, source.array, {"b.tif": 1000.0, small: 2000.0})
        assert [row[7:] for row in taken[2:7]] == self.INSIDE[small]
        assert all(set(row) == {"b"} for row in taken[:2] + taken[7:])
        oracle = deepest_interior(tiles, result.tile.meta).astype(np.float32)
        assert same_array(result.tile.array, oracle)

    def test_a_corner_of_three_tiles(self, build: Any) -> None:
        """`a` 10 x 7 down the west, `b` 7 x 8 and `c` 6 x 8 stacked on the
        east; rows 4-6, columns 4-6 are in all three. At (6, 6) `c` is 2 deep
        and the others 0; at (5, 5) all three are 1 deep, so `a`."""
        source = whole(10, 12)
        tiles, offsets = offset_tiles(
            {
                "a.tif": piece(source, 0, 10, 0, 7),
                "b.tif": piece(source, 0, 7, 4, 12),
                "c.tif": piece(source, 4, 10, 4, 12),
            }
        )
        result, _ = build(tiles)
        assert taken_from(result, source.array, offsets) == [
            "aaaaaaabbbbb",
            "aaaaaabbbbbb",
            "aaaaaabbbbbb",
            "aaaaaabbbbbb",
            "aaaaaabbbbbb",
            "aaaaaabbbbbb",
            "aaaaaacccccb",
            "aaaaaacccccc",
            "aaaaaacccccc",
            "aaaaaaaccccc",
        ]
        oracle = deepest_interior(tiles, result.tile.meta).astype(np.float32)
        assert same_array(result.tile.array, oracle)

    #: Quadrants of a 10 x 12 grid cut at row 5 and column 6 with overlap 3:
    #: 8 x 9, 8 x 6, 5 x 9 and 5 x 6 tiles, all four sharing rows 5-7 and
    #: columns 6-8. Named so the name order is nw < ne < sw < se, and then its
    #: reverse; the ties move with it, the deeper tile does not.
    CORNERS: ClassVar[dict[str, list[str]]] = {
        "abcd": [
            "aaaaaaaaabbb",
            "aaaaaaaabbbb",
            "aaaaaaaabbbb",
            "aaaaaaaabbbb",
            "aaaaaaaabbbb",
            "aaaaaaaabbbb",
            "aaaaaaaabbbb",
            "acccccccdddb",
            "ccccccccdddd",
            "cccccccccddd",
        ],
        "dcba": [
            "ddddddcccccc",
            "dddddddccccc",
            "dddddddccccc",
            "dddddddccccc",
            "dddddddccccc",
            "bddddddcccca",
            "bbbbbbbaaaaa",
            "bbbbbbbaaaaa",
            "bbbbbbbaaaaa",
            "bbbbbbaaaaaa",
        ],
    }

    @staticmethod
    def four(
        letters: str, nodata: float | None = None
    ) -> tuple[DemTile, dict[str, DemTile], dict[str, float]]:
        source = whole(10, 12, nodata=nodata)
        cut = quadrants(source, row_cut=5, col_cut=6, overlap=3)
        named = {
            f"{letter}.tif": cut[q]
            for letter, q in zip(letters, ("nw.tif", "ne.tif", "sw.tif", "se.tif"), strict=True)
        }
        tiles, offsets = offset_tiles(named)
        return source, tiles, offsets

    @pytest.mark.parametrize("letters", ["abcd", "dcba"])
    def test_a_corner_of_four_tiles(self, build: Any, letters: str) -> None:
        source, tiles, offsets = self.four(letters)
        result, _ = build(tiles)
        assert taken_from(result, source.array, offsets) == self.CORNERS[letters]
        oracle = deepest_interior(tiles, result.tile.meta).astype(np.float32)
        assert same_array(result.tile.array, oracle)

    @pytest.mark.parametrize("nodata", [np.nan, SENTINEL], ids=["nan", "sentinel"])
    def test_nodata_in_the_deepest_tile_loses_to_a_shallower_valid_value(
        self, build: Any, nodata: float
    ) -> None:
        """At (5, 6) `a.tif` (nw) is 2 deep and the other three 0: with NoData
        there, the three valid values tie, and `b.tif` sorts first. At (4, 7),
        in `a.tif` and `b.tif` only, `a.tif` is 1 deep and `b.tif` 1: a tie
        `a.tif` would take, but its NoData leaves `b.tif`."""
        sentinel = None if np.isnan(nodata) else SENTINEL
        source, tiles, offsets = self.four("abcd", sentinel)
        array = np.asarray(tiles["a.tif"].array).copy()
        array[5, 6] = array[4, 7] = nodata
        tiles["a.tif"] = DemTile(meta=tiles["a.tif"].meta, array=array)
        result, _ = build(tiles)
        expected = [list(row) for row in self.CORNERS["abcd"]]
        assert (expected[5][6], expected[4][7]) == ("a", "a")  # with no NoData
        expected[5][6] = expected[4][7] = "b"
        assert taken_from(result, source.array, offsets) == ["".join(r) for r in expected]

    def test_every_order_gives_the_same_mosaic_and_report(self, mz: ModuleType, plan: Any) -> None:
        """I2 under the new rule: all 24 assembly orders of the four-tile
        corner, with NoData in the deepest tile at one node, give one array bit
        for bit and one seam report."""
        _, tiles, _ = self.four("dcba")
        array = np.asarray(tiles["d.tif"].array).copy()
        array[5, 6] = np.nan
        tiles["d.tif"] = DemTile(meta=tiles["d.tif"].meta, array=array)
        planned = plan(tiles)
        first = mz.assemble(planned, Loads(tiles))
        for order in itertools.permutations(planned.tiles):
            result = mz.assemble(planned.model_copy(update={"tiles": order}), Loads(tiles))
            names = [p.name for p in order]
            assert same_array(result.tile.array, first.tile.array), names
            assert seams(result) == seams(first), names
        oracle = deepest_interior(tiles, first.tile.meta).astype(np.float32)
        assert same_array(first.tile.array, oracle)


class TestQ1SeamReport:
    """Ola's Q1 revised: each disagreeing seam is reported, the two tiles (by
    name, the one that sorts first first), the number of nodes where both hold
    a valid value and the values differ, and the largest and the median
    |difference| over those nodes. Pairs that agree are not listed. A
    `Mosaic`'s `seams` is a tuple sorted by the pair's names.
    """

    def test_every_pair_of_the_four_corner(self, build: Any) -> None:
        """Every overlap of the offset quadrants disagrees by a constant, so
        each pair's count is its overlap's node count: 8 x 3 for a | b (nw,
        ne), 3 x 9 for a | c (nw, sw), the 3 x 3 corner for a | d and b | c,
        3 x 6 for b | d and 5 x 3 for c | d."""
        _, tiles, _ = TestQ1DeepestInterior.four("abcd")
        result, _ = build(tiles)
        assert seams(result) == [
            ("a.tif", "b.tif", 24, 1000.0, 1000.0),
            ("a.tif", "c.tif", 27, 2000.0, 2000.0),
            ("a.tif", "d.tif", 9, 3000.0, 3000.0),
            ("b.tif", "c.tif", 9, 1000.0, 1000.0),
            ("b.tif", "d.tif", 18, 2000.0, 2000.0),
            ("c.tif", "d.tif", 15, 1000.0, 1000.0),
        ]
        assert seams(result) == seams_of(tiles, result.tile.meta)

    @pytest.mark.parametrize("overlap", [0, 1, 3], ids=["abutting", "shared-line", "overlap-3"])
    def test_agreeing_seams_are_not_listed(self, build: Any, overlap: int) -> None:
        area = overlap == 0
        source = whole(8, 12, area=area) if area else whole(9, 13)
        result, _ = build(quadrants(source, row_cut=4, col_cut=6, overlap=overlap))
        assert result.seams == ()

    def test_one_tile_has_no_seams(self, build: Any) -> None:
        result, _ = build({"only.tif": whole()})
        assert result.seams == ()

    def test_only_the_disagreeing_pair_is_listed(self, build: Any) -> None:
        """Four quadrants that agree but for one node on the nw | ne seam
        (outside the other two tiles): one entry, not six."""
        source = whole(9, 13)
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=3)
        array = np.asarray(tiles["ne.tif"].array).copy()
        array[1, 1] -= 6.0
        tiles["ne.tif"] = DemTile(meta=tiles["ne.tif"].meta, array=array)
        result, _ = build(tiles)
        assert seams(result) == [("ne.tif", "nw.tif", 1, 6.0, 6.0)]

    def test_a_three_way_node_is_counted_in_each_disagreeing_pair(self, build: Any) -> None:
        """At node (0, 0) `a` and `b` hold 1.0 and `c` 5.0: a | c and b | c
        each count it; a | b agree and are not listed."""
        source = whole(rows=3, cols=3)
        tiles = {}
        for name, first in (("a.tif", 1.0), ("b.tif", 1.0), ("c.tif", 5.0)):
            array = source.array.copy()
            array[0, 0] = first
            tiles[name] = piece(source, 0, 3, 0, 3, array=array)
        result, _ = build(tiles)
        assert seams(result) == [
            ("a.tif", "c.tif", 1, 4.0, 4.0),
            ("b.tif", "c.tif", 1, 4.0, 4.0),
        ]

    def test_the_report_covers_the_mosaics_nodes_only(self, build: Any) -> None:
        """A difference outside `--bbox` is in no node the run reads, so it is
        not reported; one inside is."""
        source = whole(9, 13)
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=3)
        array = np.asarray(tiles["ne.tif"].array).copy()
        array[0, 0] += 2.0  # global (0, 6): outside the box
        array[3, 1] += 5.0  # global (3, 7): inside it
        tiles["ne.tif"] = DemTile(meta=tiles["ne.tif"].meta, array=array)
        box = (X0 + 5 * DX, Y0 - 6 * DY, X0 + 8 * DX, Y0 - 2 * DY)  # rows 2-6, cols 5-8
        result, _ = build(tiles, box)
        assert seams(result) == [("ne.tif", "nw.tif", 1, 5.0, 5.0)]

    def test_differences_are_taken_in_float64(self, build: Any) -> None:
        """0.5 against 2**24 + 2, both exact in float32: the difference is
        16777217.5 in float64 and would round to 16777218 in float32."""
        source = whole(rows=3, cols=3)
        tiles = {}
        for name, first in (("a.tif", 0.5), ("b.tif", 2.0**24 + 2)):
            array = source.array.copy()
            array[1, 1] = first
            tiles[name] = piece(source, 0, 3, 0, 3, array=array)
        result, _ = build(tiles)
        assert seams(result) == [("a.tif", "b.tif", 1, 16777217.5, 16777217.5)]


class TestSeamThreshold:
    """Ola, 2026-09-28: "Ignore below 1mm". A seam counts, and the report
    lists, only nodes where both tiles hold a valid value and
    |a - b| >= 1 mm, 0.001 in the DEM's units, in float64; `nodes`, `largest`
    and `median` are over those nodes only, and a pair with none is not
    listed. Which tile's value a node takes does not depend on the threshold.
    """

    #: Exactly 1 mm (the float64 0.001 the rule compares with) and the float64
    #: just below it. Planted as `gap` against 0.0, so `|a - b|` is the float64
    #: itself: no rounding stands between the planted value and the compare.
    AT = 0.001
    BELOW = float(np.nextafter(0.001, 0.0))

    @staticmethod
    def same_grid(first: float, second: float) -> dict[str, DemTile]:
        """`a.tif` and `b.tif` on one float64 3 x 3 grid, equal but at the
        centre, where they hold `first` and `second`. Both are 1 deep there:
        a tie, `a.tif`'s value."""
        source = whole(rows=3, cols=3, dtype=np.float64)
        tiles = {}
        for name, centre in (("a.tif", first), ("b.tif", second)):
            array = source.array.copy()
            array[1, 1] = centre
            tiles[name] = piece(source, 0, 3, 0, 3, array=array)
        return tiles

    @pytest.mark.parametrize("larger", ["a.tif", "b.tif"])
    def test_exactly_one_millimetre_counts(self, build: Any, larger: str) -> None:
        assert self.AT == 1e-3  # the literal the ruling names, not a neighbour
        first, second = (self.AT, 0.0) if larger == "a.tif" else (0.0, self.AT)
        result, _ = build(self.same_grid(first, second))
        assert seams(result) == [("a.tif", "b.tif", 1, self.AT, self.AT)]
        assert result.tile.array[1, 1] == first

    @pytest.mark.parametrize("larger", ["a.tif", "b.tif"])
    def test_just_below_one_millimetre_does_not(self, build: Any, larger: str) -> None:
        assert 0.0 < self.BELOW < self.AT
        first, second = (self.BELOW, 0.0) if larger == "a.tif" else (0.0, self.BELOW)
        result, _ = build(self.same_grid(first, second))
        assert result.seams == ()
        assert result.tile.array[1, 1] == first

    def test_nodes_largest_and_median_are_over_the_qualifying_nodes_only(self, build: Any) -> None:
        """`TestM8Overlaps`' pair with overlap 3; `e.tif` differs at four
        nodes, by 0.0005, 0.0009, 2 and 4. Two qualify: nodes 2, largest 4,
        median 3. Over all four the count would be 4 and the median ~1."""
        source, arrays = TestM8Overlaps.pair(3)
        arrays["e.tif"][1, 0] += np.float32(0.0005)  # global (1, 4)
        arrays["e.tif"][1, 1] += np.float32(0.0009)  # global (1, 5)
        arrays["e.tif"][2, 1] += np.float32(2.0)  # global (2, 5)
        arrays["e.tif"][2, 2] += np.float32(4.0)  # global (2, 6)
        tiles = TestM8Overlaps.tiles(source, arrays)
        result, _ = build(tiles)
        assert seams(result) == [("e.tif", "w.tif", 2, 4.0, 3.0)]
        assert seams(result) == seams_of(tiles, result.tile.meta)

    def test_a_pair_below_the_threshold_everywhere_is_not_listed(self, build: Any) -> None:
        """Quadrants with overlap 3: `ne.tif` 0.5 mm off `nw.tif` at global
        (1, 7), in those two only; `sw.tif` 3 off `nw.tif` at global (5, 1),
        in those two only. One entry, nw | sw."""
        source = whole(9, 13)
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=3)
        for name, (r, c), delta in (("ne.tif", (1, 1), 0.0005), ("sw.tif", (1, 1), 3.0)):
            array = np.asarray(tiles[name].array).copy()
            array[r, c] += np.float32(delta)
            tiles[name] = DemTile(meta=tiles[name].meta, array=array)
        result, _ = build(tiles)
        assert seams(result) == [("nw.tif", "sw.tif", 1, 3.0, 3.0)]

    def test_a_sub_millimetre_disagreement_is_still_decided_by_depth_then_name(
        self, mz: ModuleType, plan: Any
    ) -> None:
        """The midline rule does not look at the threshold. `TestM8Overlaps`'
        pair with overlap 3, 0.5 mm apart at three nodes of row 1: at column 4
        `w.tif` is deeper (1 against 0), at 5 a tie (`e.tif`, by name), at 6
        `e.tif` is deeper. So the name that sorts first does not take every
        node, and in every assembly order; nothing is reported."""
        source, arrays = TestM8Overlaps.pair(3)
        arrays["w.tif"][1, 4] += np.float32(0.0005)  # global (1, 4): w.tif deeper
        arrays["e.tif"][1, 1] += np.float32(0.0005)  # global (1, 5): a tie
        arrays["e.tif"][1, 2] -= np.float32(0.0005)  # global (1, 6): e.tif deeper
        tiles = TestM8Overlaps.tiles(source, arrays)
        expected = source.array.copy()
        expected[1, 4] = arrays["w.tif"][1, 4]
        expected[1, 5] = arrays["e.tif"][1, 1]
        expected[1, 6] = arrays["e.tif"][1, 2]
        assert (expected[1, 4:7] != source.array[1, 4:7]).all()  # all three planted
        planned = plan(tiles)
        for order in itertools.permutations(planned.tiles):
            result = mz.assemble(planned.model_copy(update={"tiles": order}), Loads(tiles))
            names = [p.name for p in order]
            assert same_array(result.tile.array, expected), names
            assert result.seams == (), names


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
