"""A cut run, piece by piece (increment 23c).

`docs/increments/23-basin-scale.md`, "The subdomain model" and "The seam
protocol": the window's partition (`cut_partition`), then each piece's slice
of the one start triangulation, refined on its cell's window with the seams
frozen, and K4 over every piece's seam records (`cut_run`).

`cli.py` parses options and hands this module a start mesh, a tile and its
refine (`refined`, the uncut run's own, so both paths refine alike and a
wrap of `cli.refine` sees every piece). No Typer here, so the run is testable
without the command line: a refusal is a `RunError` in the engine's words,
a seam two pieces wrote differently a `ConformityError`, and `cli.py` makes
them a usage error and exit 1.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, replace

import numpy as np
import numpy.typing as npt

from tin_engine._core import IndexedMesh2, RefineOutcome, SeamOutcome, indexed_mesh, refine_seam
from tin_engine.decompose import bytes_per_node, partition
from tin_engine.elevation import Trimmed
from tin_engine.features import PIECE_VOCABULARY
from tin_engine.io.mesh_index import PartitionRecord, PieceId, PieceWindow, check_conformity
from tin_engine.io.models import DemTile, RasterMeta
from tin_engine.pieces import piece_ids, seam_records, start_slice
from tin_engine.raster import to_core
from tin_engine.run_record import Value
from tin_engine.stats import PhaseClock
from tin_engine.target_grid import Block, TargetGrid


class RunError(ValueError):
    """A run the engine refuses, in its own words; no file is written."""


class Once[T]:
    """A value handed on once (15e, fix 3): `take` returns it and forgets it,
    so the holder keeps nothing alive after; a second `take` is a bug."""

    def __init__(self, value: T) -> None:
        self._value: T | None = value

    def take(self) -> T:
        value, self._value = self._value, None
        if value is None:  # explicit: `python -O` strips an assert
            raise RuntimeError("Once.take called twice")
        return value


@dataclass(frozen=True, slots=True)
class How:
    """How a start mesh is refined (15c-2's phase 2 with ``grid`` and
    ``checks``)."""

    dem: str
    tolerance: float
    min_angle: float
    feet: bool
    grid: TargetGrid | None
    checks: Iterable[Block] | None


@dataclass(frozen=True, slots=True)
class Refined:
    out: RefineOutcome
    final: RefineOutcome
    z: npt.NDArray[np.float64]
    trimmed: Trimmed
    values: Mapping[str, Value]  # the run record's entries (increment 25)
    on_frozen: int


#: A piece's refine, given by `cli.py` (its `_refined`): the held window, the
#: start, its constraint edges and masks, how, the clock, the frozen mask and
#: the seam records' heights.
Refine = Callable[
    [Once[DemTile], IndexedMesh2, npt.NDArray[np.uint32], npt.NDArray[np.uint32], How,
     PhaseClock, int, Mapping[int, float] | None],
    Refined,
]  # fmt: skip


@dataclass(frozen=True, slots=True)
class CutAsk:
    """``--pieces`` and ``--memory-budget``, and the window on its lattice:
    lattice node (0, 0) at ``origin``, the window's first node ``(row0, col0)``."""

    pieces: int
    budget: int
    origin: tuple[float, float]
    row0: int
    col0: int


@dataclass(frozen=True, slots=True)
class Partitioned:
    """A partition of more than one cell: its record and its lines' x and y."""

    record: PartitionRecord
    xs: list[float]
    ys: list[float]
    row0: int
    col0: int


@dataclass(frozen=True, slots=True)
class Piece:
    id: PieceId
    trimmed: Trimmed
    values: Mapping[str, Value]  # its own record entries, never the run's
    window: PieceWindow
    on_frozen: int


@dataclass(frozen=True, slots=True)
class CutRun:
    record: PartitionRecord
    pieces: list[Piece]


def cut_partition(meta: RasterMeta, tolerance: float, ask: CutAsk) -> Partitioned | None:
    """The window's partition ("The partition", 23c), or None for one cell,
    which is today's run (K1). Line ``i`` is at the index's
    ``origin_x + (col0 + i dx) h``, line ``j`` at ``origin_y - (row0 + j dy) h``."""
    cells = partition(meta.cols, meta.rows, tolerance, pieces=ask.pieces, memory_budget=ask.budget)
    if cells.nx * cells.ny == 1:
        return None
    h = meta.delta_x
    if meta.delta_y != h:
        raise RunError(f"the DEM's cells are {h} x {meta.delta_y}; a cut needs square ones")
    (ox, oy), r0, c0 = ask.origin, ask.row0, ask.col0
    record = PartitionRecord(
        pieces=ask.pieces, memory_budget=ask.budget, bytes_per_node=bytes_per_node(tolerance),
        nx=cells.nx, ny=cells.ny, dx=cells.dx, dy=cells.dy, cols=meta.cols, rows=meta.rows,
        row0=r0, col0=c0, spacing=h, origin=(ox, oy),
    )  # fmt: skip
    xs = [ox + (c0 + i * cells.dx) * h for i in range(1, cells.nx)]
    ys = [oy - (r0 + j * cells.dy) * h for j in range(1, cells.ny)]
    return Partitioned(record, xs, ys, r0, c0)


def tile_window(tile: DemTile, r0: int, r1: int, c0: int, c1: int) -> tuple[DemTile, int, int]:
    """Rows ``r0..r1`` and columns ``c0..c1`` of ``tile``, inclusive and cut to
    it, as a tile of their own; with its first row and column."""
    m = tile.meta
    r0, c0, r1, c1 = max(r0, 0), max(c0, 0), min(r1, m.rows - 1), min(c1, m.cols - 1)
    meta = m.model_copy(
        update={"x_min": m.x_min + c0 * m.delta_x, "y_max": m.y_max - r0 * m.delta_y,
                "rows": r1 - r0 + 1, "cols": c1 - c0 + 1}
    )  # fmt: skip
    return DemTile(meta=meta, array=tile.array[r0 : r1 + 1, c0 : c1 + 1]), r0, c0


def cut_run(
    tile: DemTile,
    start: IndexedMesh2,
    edges: npt.NDArray[np.uint32],
    masks: npt.NDArray[np.uint32],
    cut: Partitioned,
    how: How,
    clock: PhaseClock,
    refined: Refine,
) -> CutRun:
    """23c, one piece after another: each piece's slice of the one start
    triangulation, its seam edges fanned at their seam-pass points, refined
    on its cell's window (grown by a node) with the seams frozen; then K4
    over the seam records of every piece (``ConformityError`` exits 1)."""
    seam, part = PIECE_VOCABULARY.mask("seam"), cut.record
    xy = np.asarray(start.vertices)
    tris = np.asarray(start.triangles, dtype=np.int64)
    flags = np.asarray(start.constrained_edges)
    seams = edges[(masks & seam) != 0].astype(np.int64)
    labels, ids = piece_ids(
        tris, xy, seams, lambda x, y: (sum(y < v for v in cut.ys), sum(x > v for v in cut.xs))
    )
    m = tile.meta

    def seam_pass(r: int) -> SeamOutcome:
        """The seam pass on a strip that is a function of the edge alone: its
        box snapped outward to the lattice and grown by one node (step 2)."""
        (ax, ay), (bx, by) = xy[edges[r]].tolist()
        c0, c1 = ((v - m.x_min) / m.delta_x for v in (min(ax, bx), max(ax, bx)))
        r0, r1 = ((m.y_max - v) / m.delta_y for v in (max(ay, by), min(ay, by)))
        f, c = math.floor, math.ceil
        strip, _, _ = tile_window(tile, f(r0) - 1, c(r1) + 1, f(c0) - 1, c(c1) + 1)
        return refine_seam(to_core(strip), (ax, ay), (bx, by), tolerance=how.tolerance)

    checks = list(how.checks) if how.checks is not None else None
    built, records = [], []
    for p, (j, i, k) in enumerate(ids):
        piece = start_slice(xy, tris, flags, edges, masks, labels == p, seam, seam_pass)
        window, r0, c0 = tile_window(
            tile, j * part.dy - 1, (j + 1) * part.dy + 1, i * part.dx - 1, (i + 1) * part.dx + 1
        )
        wm, grid = window.meta, how.grid
        if grid is not None:  # phase 2's store on the piece's window
            grid = grid.model_copy(update={"row0": grid.row0 + r0, "col0": grid.col0 + c0,
                                           "rows": wm.rows, "cols": wm.cols})  # fmt: skip
        mine = replace(how, grid=grid, checks=checks)
        mesh = indexed_mesh(piece.vertices, piece.triangles, piece.flags)
        done = refined(Once(window), mesh, piece.edges, piece.masks, mine, clock, seam, piece.z)
        del window
        records += seam_records((j, i, k), piece, done.final.vertices, done.z, done.final.valid)
        place = PieceWindow(row0=cut.row0 + r0, col0=cut.col0 + c0, rows=wm.rows, cols=wm.cols)
        built.append(Piece((j, i, k), done.trimmed, done.values, place, done.on_frozen))
    if not any(len(p.trimmed.triangles) for p in built):
        raise RunError(f"{how.dem} has no data under any triangle; nothing to write")
    check_conformity(records)  # K4: a ConformityError fails the run
    return CutRun(part, built)
