"""`assemble(plan, load, load_window=...)`: a windowed assembly equals the whole one (23a-1, W6).

`docs/increments/23-basin-scale.md`, "Windowed source reads (23a-1)",
decided 1: no second assembler; `assemble` takes each placement's source
window from `load_window` instead of slicing a whole tile. The oracle is
`assemble` without it, on 15a's fixtures (`mosaic_fixtures`): quadrants with
shared and true overlaps, a disagreeing overlap (so there are seams), mixed
dtypes, with and without a box. 15a's own suite is not edited.

HOW THIS FILE GOES RED. `window_meta` is reached through a module-scoped
fixture and `load_window` is a keyword argument, so if either is missing
(`ModuleNotFoundError`, `TypeError`) each test fails on its own and the rest
of `tests/python` still collects.
"""

from __future__ import annotations

import importlib
from typing import Any

import numpy as np
import pytest

from mosaic_fixtures import X0, Y0, piece, quadrants, same_array, whole
from tin_engine.io.models import DemTile
from tin_engine.io.repository import TileFootprint
from tin_engine.mosaic import Bounds, MosaicError, assemble, plan_mosaic


@pytest.fixture(scope="module")
def window_meta() -> Any:
    return importlib.import_module("tin_engine.io.cog").window_meta


class Windows:
    """A `load_window` over whole tiles that records what it was asked for."""

    def __init__(self, tiles: dict[str, DemTile], window_meta: Any) -> None:
        self.tiles, self.window_meta = tiles, window_meta
        self.calls: list[tuple[str, Any]] = []

    def __call__(self, name: str, window: Any) -> DemTile:
        self.calls.append((name, window))
        tile = self.tiles[name]
        rows = slice(window.row0, window.row0 + window.rows)
        cols = slice(window.col0, window.col0 + window.cols)
        return DemTile(meta=self.window_meta(tile.meta, window), array=tile.array[rows, cols])


def disagreeing(tiles: dict[str, DemTile]) -> dict[str, DemTile]:
    """`tiles` with the north-east tile 0.5 higher, so every overlap with it is a seam."""
    ne = tiles["ne.tif"]
    out = dict(tiles)
    out["ne.tif"] = DemTile(meta=ne.meta, array=np.asarray(ne.array) + np.float32(0.5))
    return out


def mixed(tiles: dict[str, DemTile]) -> dict[str, DemTile]:
    out = dict(tiles)
    se = tiles["se.tif"]
    out["se.tif"] = DemTile(meta=se.meta, array=np.asarray(se.array, dtype=np.float64))
    return out


CASES = {
    "shared_line": lambda: quadrants(whole(), row_cut=4, col_cut=6, overlap=1),
    "true_overlap": lambda: quadrants(whole(), row_cut=4, col_cut=6, overlap=3),
    "disagreeing_overlap": lambda: disagreeing(quadrants(whole(), row_cut=4, col_cut=6, overlap=3)),
    "mixed_dtypes": lambda: mixed(quadrants(whole(), row_cut=4, col_cut=6, overlap=2)),
}
BOXES: dict[str, Bounds | None] = {
    "no_box": None,
    "inner_box": Bounds(x_min=X0 + 20.0, y_min=Y0 - 35.0, x_max=X0 + 100.0, y_max=Y0 - 5.0),
}


def footprints(tiles: dict[str, DemTile]) -> list[TileFootprint]:
    return [
        TileFootprint(name=name, meta=tile.meta, dtype=tile.array.dtype)
        for name, tile in tiles.items()
    ]


@pytest.mark.parametrize("box", BOXES)
@pytest.mark.parametrize("case", CASES)
def test_w6_windowed_assembly_equals_whole_assembly(window_meta: Any, case: str, box: str) -> None:
    tiles = CASES[case]()
    plan = plan_mosaic(footprints(tiles), BOXES[box], None)
    whole_result = assemble(plan, tiles.__getitem__)
    windows = Windows(tiles, window_meta)
    windowed = assemble(plan, tiles.__getitem__, load_window=windows)
    assert windowed.tile.meta == whole_result.tile.meta
    assert same_array(windowed.tile.array, whole_result.tile.array)
    assert windowed.seams == whole_result.seams
    assert sorted(windows.calls, key=lambda c: c[0]) == sorted(
        ((p.name, p.source) for p in plan.tiles), key=lambda c: c[0]
    )
    if case == "disagreeing_overlap" and box == "no_box":
        assert whole_result.seams, "the fixture must produce seams"


def test_w6_a_window_whose_meta_is_not_window_meta_changed_since_listed(window_meta: Any) -> None:
    tiles = quadrants(whole(), row_cut=4, col_cut=6, overlap=1)
    plan = plan_mosaic(footprints(tiles), None, None)

    def shifted(name: str, window: Any) -> DemTile:
        tile = Windows(tiles, window_meta)(name, window)
        moved = tile.meta.model_copy(update={"x_min": tile.meta.x_min + 10.0})
        return DemTile(meta=moved, array=tile.array)

    with pytest.raises(MosaicError, match="changed since it was listed"):
        assemble(plan, tiles.__getitem__, load_window=shifted)


def test_w6_a_window_of_the_whole_tile_is_the_single_tile_case(window_meta: Any) -> None:
    """One tile on the mosaic's grid: the window is the whole tile."""
    source = whole()
    tiles = {"only.tif": piece(source, 0, 9, 0, 13)}
    plan = plan_mosaic(footprints(tiles), None, None)
    windows = Windows(tiles, window_meta)
    result = assemble(plan, tiles.__getitem__, load_window=windows)
    assert same_array(result.tile.array, source.array)
    assert len(windows.calls) == 1
