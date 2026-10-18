"""A cut run, piece by piece.

`docs/increments/23-basin-scale.md`, *The seam protocol*: the window's
partition (`cut_partition`), then each piece's slice of the one start
triangulation, refined on its cell's window with the seams
frozen, and K4 over every piece's seam records (`cut_run`).

`cli.py` parses options and hands this module a start mesh, a `load` that
reads one piece's window of the DEM, the refine (`refined`, the uncut run's
own, so both paths refine alike and a wrap of `cli.refine` sees every piece)
and a `write` for each finished piece. One piece is loaded, refined, written
and dropped before the next is loaded, so a run holds one piece's DEM at a
time. No Typer here, and no path, so the run is testable without the command
line: a refusal is a `RunError` in the engine's words, a seam two pieces
wrote differently a `ConformityError`, and `cli.py` makes them a usage error
and exit 1.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, replace
from typing import Protocol

import numpy as np
import numpy.typing as npt

from tin_engine._core import IndexedMesh2, SeamOutcome, indexed_mesh, refine_seam
from tin_engine.decompose import Partition, bytes_per_node
from tin_engine.elevation import Trimmed
from tin_engine.features import PIECE_VOCABULARY
from tin_engine.io.mesh_index import (
    PartitionRecord,
    PieceEntry,
    PieceId,
    PieceWindow,
    SeamRecord,
    check_conformity,
)
from tin_engine.io.models import DemTile, IndexWindow, RasterMeta
from tin_engine.pieces import piece_ids, seam_records, start_slice
from tin_engine.raster import to_core
from tin_engine.refine_job import How, Once, Refined, RunError
from tin_engine.run_record import Value
from tin_engine.stats import PhaseClock
from tin_engine.target_grid import Block, TargetGrid

#: The nodes a piece's window reaches past its cell on each side. Not free to
#: change: a piece's bytes depend on how far its raster reaches.
MARGIN = 1


class Loaded(Protocol):
    """One piece's DEM, as `load` gives it: its tile, whose node (0, 0) is
    the run window's ``(row0, col0)``, and on a resampled run its sub-grid
    and check points."""

    @property
    def tile(self) -> DemTile: ...
    @property
    def row0(self) -> int: ...
    @property
    def col0(self) -> int: ...
    @property
    def grid(self) -> TargetGrid | None: ...
    @property
    def checks(self) -> Iterable[Block] | None: ...


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
    """``--pieces`` and ``--piece-budget``, and the window on its lattice:
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
    meta: RasterMeta  # the run window's


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
    entries: list[PieceEntry]


def cut_partition(meta: RasterMeta, tolerance: float, cells: Partition, ask: CutAsk) -> Partitioned:
    """The window's partition into ``cells``, of more than one cell (the
    design's "The partition"). Line ``i`` is at the index's
    ``origin_x + (col0 + i dx) h``, line ``j`` at ``origin_y - (row0 + j dy) h``."""
    h = meta.delta_x
    if meta.delta_y != h:
        raise RunError(f"the DEM's cells are {h} x {meta.delta_y}; a cut needs square ones")
    (ox, oy), r0, c0 = ask.origin, ask.row0, ask.col0
    record = PartitionRecord(
        pieces=ask.pieces, piece_budget=ask.budget, bytes_per_node=bytes_per_node(tolerance),
        nx=cells.nx, ny=cells.ny, dx=cells.dx, dy=cells.dy, cols=meta.cols, rows=meta.rows,
        row0=r0, col0=c0, spacing=h, origin=(ox, oy),
    )  # fmt: skip
    xs = [ox + (c0 + i * cells.dx) * h for i in range(1, cells.nx)]
    ys = [oy - (r0 + j * cells.dy) * h for j in range(1, cells.ny)]
    return Partitioned(record, xs, ys, r0, c0, meta)


def piece_windows(cut: Partitioned) -> dict[tuple[int, int], IndexWindow]:
    """Each cell ``(j, i)`` grown by `MARGIN` nodes and cut to the run window,
    in the run window's indices."""
    p, windows = cut.record, {}
    for j in range(p.ny):
        for i in range(p.nx):
            r0, c0 = max(j * p.dy - MARGIN, 0), max(i * p.dx - MARGIN, 0)
            r1 = min((j + 1) * p.dy + MARGIN, p.rows - 1)
            c1 = min((i + 1) * p.dx + MARGIN, p.cols - 1)
            windows[j, i] = IndexWindow(row0=r0, col0=c0, rows=r1 - r0 + 1, cols=c1 - c0 + 1)
    return windows


def strip(dem: Loaded, meta: RasterMeta, r0: int, r1: int, c0: int, c1: int) -> DemTile:
    """Rows ``r0..r1`` and columns ``c0..c1`` of the run window, inclusive and
    cut to it, from the piece's tile. Its corner is the run's ``meta`` moved by
    whole cells from the run window's corner, so two pieces build one strip
    bit for bit; a strip the piece does not hold is refused."""
    r0, c0, r1, c1 = max(r0, 0), max(c0, 0), min(r1, meta.rows - 1), min(c1, meta.cols - 1)
    i0, i1, j0, j1 = r0 - dem.row0, r1 - dem.row0, c0 - dem.col0, c1 - dem.col0
    if min(i0, j0) < 0 or i1 >= dem.tile.meta.rows or j1 >= dem.tile.meta.cols:
        raise RunError(f"a seam strip, rows {r0}-{r1}, columns {c0}-{c1}, is outside its piece")
    at = meta.model_copy(
        update={"x_min": meta.x_min + c0 * meta.delta_x, "y_max": meta.y_max - r0 * meta.delta_y,
                "rows": r1 - r0 + 1, "cols": c1 - c0 + 1}
    )  # fmt: skip
    return DemTile(meta=at, array=dem.tile.array[i0 : i1 + 1, j0 : j1 + 1])


def cut_run(
    load: Callable[[IndexWindow], Loaded],
    start: IndexedMesh2,
    edges: npt.NDArray[np.uint32],
    masks: npt.NDArray[np.uint32],
    cut: Partitioned,
    how: How,
    clock: PhaseClock,
    refined: Refine,
    write: Callable[[Piece], PieceEntry],
) -> CutRun:
    """One piece after another, in index order: its window loaded, its slice
    of the one start triangulation with its seam edges fanned at their
    seam-pass points, refined with the seams frozen, written, and dropped
    before the next is loaded; then K4 over the seam records of every piece
    (``ConformityError`` exits 1)."""
    seam, part, m = PIECE_VOCABULARY.mask("seam"), cut.record, cut.meta
    xy = np.asarray(start.vertices)
    tris = np.asarray(start.triangles, dtype=np.int64)
    flags = np.asarray(start.constrained_edges)
    seams = edges[(masks & seam) != 0].astype(np.int64)
    labels, ids = piece_ids(
        tris, xy, seams, lambda x, y: (sum(y < v for v in cut.ys), sum(x > v for v in cut.xs))
    )
    windows = piece_windows(cut)

    def one(p: int, j: int, i: int, k: int) -> tuple[list[SeamRecord], PieceEntry]:
        """Piece ``p``, in a frame of its own, so nothing of it outlives its write."""
        w = windows[j, i]
        dem = load(w)

        def seam_pass(r: int) -> SeamOutcome:
            """The seam pass on a strip that is a function of the edge alone: its
            box snapped outward to the lattice and grown by one node (step 2)."""
            (ax, ay), (bx, by) = xy[edges[r]].tolist()
            c0, c1 = ((v - m.x_min) / m.delta_x for v in (min(ax, bx), max(ax, bx)))
            r0, r1 = ((m.y_max - v) / m.delta_y for v in (max(ay, by), min(ay, by)))
            f, c = math.floor, math.ceil
            box = strip(dem, m, f(r0) - 1, c(r1) + 1, f(c0) - 1, c(c1) + 1)
            return refine_seam(to_core(box), (ax, ay), (bx, by), tolerance=how.tolerance)

        piece = start_slice(xy, tris, flags, edges, masks, labels == p, seam, seam_pass)
        mine, held = replace(how, grid=dem.grid, checks=dem.checks), Once(dem.tile)
        # The constraints travel in `edges` and `masks`; refine reads no flag.
        zeros = np.zeros(len(piece.triangles), dtype=np.uint8)
        mesh = indexed_mesh(piece.vertices, piece.triangles, zeros)
        done = refined(held, mesh, piece.edges, piece.masks, mine, clock, seam, piece.z)
        records = seam_records((j, i, k), piece, done.final.vertices, done.z, done.final.valid)
        place = PieceWindow(
            row0=cut.row0 + w.row0, col0=cut.col0 + w.col0, rows=w.rows, cols=w.cols
        )
        return records, write(Piece((j, i, k), done.trimmed, done.values, place, done.on_frozen))

    entries, records = [], []
    for p, (j, i, k) in enumerate(ids):
        found, entry = one(p, j, i, k)
        records += found
        entries.append(entry)
    if not any(e.counts.triangles for e in entries):
        raise RunError(f"{how.dem} has no data under any triangle; nothing to write")
    check_conformity(records)  # K4: a ConformityError fails the run
    return CutRun(part, entries)
