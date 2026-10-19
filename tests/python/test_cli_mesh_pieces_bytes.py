"""MP5 and SL1: reading each piece by itself changes no byte (K6).

`docs/increments/23-basin-scale.md`, "23c-3: each piece's memory its own",
"Red tests", MP5 and SL1. MP5 is the invariant-critical test of the step.

- MP5: `dem_input.load_piece` is replaced by a double that reads as 23c-2
  does: the run's DEM assembled once (`dem_input.assemble_dem`, today's
  `open_dem` second half), each piece's window sliced from it with 23c-2's
  meta formula (`x_min + c dx`, `y_max - r dy`); on the reprojected path the
  whole target grid resampled once, the piece's sub-grid cut from it as
  23c-2's `cut_run` does, and every piece handed every check-point block.
  Every piece file and `index.json` are byte for byte the real run's.
- SL1: with `--tolerance-slope`, each piece's `slope_nodes_tightened` is the
  same with `basin_run.MARGIN` patched to 3 as at 1 (whole files differ at
  margin 3 in the last bit of an off-node z, D5, so only the count is
  compared). That the patch reaches the run is checked on the index's piece
  windows, which grow by two nodes on every inner side.
"""

from __future__ import annotations

import dataclasses
import inspect
from pathlib import Path
from typing import Any

import pytest

import piece_spies as ps
import pieces_fixtures as pf

SLOPE = ("--tolerance-slope", "0.3", "5", "20")


def tile_window_23c2(tile: Any, r0: int, r1: int, c0: int, c1: int) -> tuple[Any, int, int]:
    """23c-2's `basin_run.tile_window`, copied: rows `r0..r1` and columns
    `c0..c1`, inclusive and cut to the tile, with its first row and column."""
    from tin_engine.io.models import DemTile

    m = tile.meta
    r0, c0, r1, c1 = max(r0, 0), max(c0, 0), min(r1, m.rows - 1), min(c1, m.cols - 1)
    meta = m.model_copy(
        update={"x_min": m.x_min + c0 * m.delta_x, "y_max": m.y_max - r0 * m.delta_y,
                "rows": r1 - r0 + 1, "cols": c1 - c0 + 1}
    )  # fmt: skip
    return DemTile(meta=meta, array=tile.array[r0 : r1 + 1, c0 : c1 + 1]), r0, c0


def install_23c2_reader(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """`load_piece` replaced by the double; the number of pieces it served."""
    from tin_engine import dem_input

    served: list[int] = []
    whole: dict[int, tuple[Any, Any, list[Any] | None]] = {}
    original = None

    def double(*args: Any, **kwargs: Any) -> Any:
        bound = inspect.signature(original).bind(*args, **kwargs)  # type: ignore[arg-type]
        plan, piece = bound.arguments["plan"], bound.arguments["piece"]
        resources = bound.arguments.get("resources")
        if id(plan) not in whole:
            run = dem_input.assemble_dem(plan, resources)
            checks = list(run.checks) if run.checks is not None else None
            whole[id(plan)] = (plan, run, checks)
        _, run, checks = whole[id(plan)]
        real = original(*args, **kwargs)  # type: ignore[misc]
        r0, r1, c0, c1 = ps.window_bounds(piece.window)
        tile, r0, c0 = tile_window_23c2(run.tile, r0, r1, c0, c1)
        grid = run.grid
        if grid is not None:
            rows, cols = tile.meta.rows, tile.meta.cols
            grid = grid.model_copy(update={"row0": grid.row0 + r0, "col0": grid.col0 + c0,
                                           "rows": rows, "cols": cols})  # fmt: skip
        served.append(len(served))
        return dataclasses.replace(
            real,
            tile=tile,
            row0=r0,
            col0=c0,
            seams=run.seams,
            grid=grid,
            checks=iter(checks) if checks is not None else None,
        )

    original = ps.patch_everywhere(monkeypatch, dem_input, "load_piece", double)
    return served


def both(where: Path, monkeypatch: pytest.MonkeyPatch, *args: str) -> tuple[dict[str, bytes], ...]:
    real = pf.cut(where / "real", *args)
    with monkeypatch.context() as mp:
        served = install_23c2_reader(mp)
        read_as_23c2 = pf.cut(where / "as-23c-2", *args)
    assert len(served) == len(read_as_23c2.index["pieces"]) > 1, "the double served no piece"
    assert real.dir is not None and read_as_23c2.dir is not None
    return ps.piece_files(real.dir), ps.piece_files(read_as_23c2.dir)


@pytest.mark.parametrize("slope", [False, True], ids=["plain", "tolerance-slope"])
@pytest.mark.parametrize("pieces", ["4", "16"])
def test_mp5_square_scene_bytes_are_23c2s(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, pieces: str, slope: bool
) -> None:
    args = (*ps.square_scene(tmp_path), "--pieces", pieces, *(SLOPE if slope else ()))
    real, as_23c2 = both(tmp_path, monkeypatch, *args)
    assert "index.json" in real
    assert real == as_23c2


def test_mp5_four_tile_scene_bytes_are_23c2s(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The overlaps agree, so the per-piece overlap report (D6) is the run's."""
    real, as_23c2 = both(tmp_path, monkeypatch, *ps.four_tile_scene(tmp_path), "--pieces", "4")
    assert real == as_23c2


@pytest.mark.parametrize("block", [512, 8])
@pytest.mark.parametrize("quality", ["off", "default"])
def test_mp5_reprojected_scene_bytes_are_23c2s(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, block: int, quality: str
) -> None:
    """DC10's scene; at `BLOCK` 8 the outline test skips blocks (49, of which
    42 meet the domain), so a block clipped to a piece before it is judged
    would drop check points the double still hands over."""
    from tin_engine import target_grid

    monkeypatch.setattr(target_grid, "BLOCK", block)
    angle = ("--start-min-angle", "0") if quality == "off" else ()
    args = (*ps.reprojected_scene(tmp_path), "--pieces", "4", *angle)
    real, as_23c2 = both(tmp_path, monkeypatch, *args)
    assert real == as_23c2


# ------------------------------------------------------------------ SL1


#: Measured on 23c-2 (`321bac56`) and in the design: pieces 0-0-0, 0-1-0,
#: 1-0-0 and 1-1-0 at margin 1.
TIGHTENED = {
    (0, 0, 0): (317, 338),
    (0, 1, 0): (276, 300),
    (1, 0, 0): (278, 312),
    (1, 1, 0): (254, 272),
}


def tightened(done: pf.Cut) -> dict[tuple[int, int, int], tuple[int, int]]:
    out = {}
    for piece in done.pieces:
        text = piece.vtk.field_data["slope_nodes_tightened"].values[0]
        n, of, total = text.split()[:3]
        assert of == "of", text
        out[piece.id] = (int(n), int(total))
    return out


def test_sl1_the_margin_does_not_move_the_slope_count(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from tin_engine import basin_run

    args = (*ps.square_scene(tmp_path), "--pieces", "4", *SLOPE, "--start-min-angle", "0")
    assert basin_run.MARGIN == 1
    at_one = pf.cut(tmp_path / "margin-1", *args)
    monkeypatch.setattr(basin_run, "MARGIN", 3)
    at_three = pf.cut(tmp_path / "margin-3", *args)
    assert tightened(at_one) == TIGHTENED
    assert tightened(at_three) == tightened(at_one)
    one = {tuple(e["id"]): e["window"] for e in at_one.index["pieces"]}
    three = {tuple(e["id"]): e["window"] for e in at_three.index["pieces"]}
    grown = [
        k
        for k in one
        if three[k]["rows"] == one[k]["rows"] + 2 and three[k]["cols"] == one[k]["cols"] + 2
    ]
    assert sorted(grown) == sorted(one), ("MARGIN 3 did not reach the windows", one, three)
