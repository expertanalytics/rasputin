"""Many DEM tiles as one node grid: planning from headers, then assembly (increment 15a).

`docs/increments/15-dem-mosaic.md` R1, R4, R5 and R7, with the two
corrections in its "Pinned by the red suite (15a)": tiles are selected by
index window, and coverage is checked per node.

`plan_mosaic` is pure: footprints and a request in, a `MosaicPlan` out, and
every refusal a header can decide fires there, before any pixel (I6).
`assemble` reads pixels through a `load` callable and never sees a path.

Positions are integers from here on. Each lattice has a reference node, its
north-west-most node over all its tiles, and a node's global index `(R, K)`
counts from it; the mosaic's first node is `X_ref + c0 * dx`, never a sum of
steps (I1, R12's global node identity).
"""

from __future__ import annotations

import math
import os
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Self

import numpy as np
import numpy.typing as npt
import shapely
from pydantic import BaseModel, ConfigDict, model_validator

from tin_engine.io.cog import window_meta
from tin_engine.io.models import DemTile, IndexWindow, RasterMeta

if TYPE_CHECKING:
    from tin_engine.io.repository import TileFootprint

#: DEM units (metres for DTM10). A seam counts a node only where the two
#: tiles differ by at least this (Ola, 2026-09-28: "Ignore below 1mm").
SEAM_THRESHOLD = 0.001

#: Cells. Two tiles share a lattice when their node offsets are integers to
#: within this. Measured noise is at most 2.2e-11 cells (B3), and 0 for DTM10.
ALIGN_TOLERANCE = 1e-6


class MosaicError(ValueError):
    """A request the tiles cannot answer. A `ValueError`, like `GeoTiffError`."""


def physical_memory() -> int:
    """Bytes of physical memory, looked up at call time (macOS and Linux)."""
    return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")


class Bounds(BaseModel):
    """A box in the DEM's CRS: finite, with `x_min < x_max` and `y_min < y_max`."""

    model_config = ConfigDict(frozen=True)

    x_min: float
    y_min: float
    x_max: float
    y_max: float

    @model_validator(mode="after")
    def _a_box(self) -> Self:
        corners = (self.x_min, self.y_min, self.x_max, self.y_max)
        if not all(map(math.isfinite, corners)) or not (
            self.x_min < self.x_max and self.y_min < self.y_max
        ):
            raise ValueError(f"need finite x_min < x_max and y_min < y_max, got {corners}")
        return self


class TilePlacement(BaseModel):
    """One selected tile: `canvas` in canvas indices, `source` in the tile's own,
    and the dtype its footprint says it decodes to."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    name: str
    meta: RasterMeta
    canvas: IndexWindow
    source: IndexWindow
    dtype: np.dtype[Any] = np.dtype(np.float32)


class MosaicPlan(BaseModel):
    """What `assemble` will build. `window` is the mosaic's first node as a
    global index, and its shape; `tiles` is sorted by name."""

    model_config = ConfigDict(frozen=True)

    meta: RasterMeta
    reference: tuple[float, float]
    window: IndexWindow
    tiles: tuple[TilePlacement, ...]


@dataclass(frozen=True, slots=True)
class Seam:
    """Two tiles whose valid values differ by at least `SEAM_THRESHOLD` at
    `nodes` of the mosaic's nodes, by `largest` and `median` there (float64);
    `first < second` by name."""

    first: str
    second: str
    nodes: int
    largest: float
    median: float

    def entry(self) -> str:
        """The `dem_seams` entry, unescaped."""
        _, _, n, largest, median = self.cells()
        return f"{self.first} | {self.second}: nodes {n}, max {largest}, median {median}"

    def cells(self) -> tuple[str, str, str, str, str]:
        """The `--stats` table row, numbers formatted as in `entry`."""
        return self.first, self.second, f"{self.nodes}", f"{self.largest:g}", f"{self.median:g}"


@dataclass(frozen=True, slots=True)
class Mosaic:
    """The assembled tile, the plan it came from, as provenance, and each pair
    of tiles whose overlap disagrees, sorted by name (Ola's Q1 revised)."""

    tile: DemTile
    plan: MosaicPlan
    seams: tuple[Seam, ...] = ()


@dataclass(frozen=True, slots=True)
class _Lattice:
    """One lattice's tiles (sorted by name), its reference node, each tile's
    global index, and the box's window on it with the tiles that meet it."""

    group: list[TileFootprint]
    reference: tuple[float, float]
    placed: list[_Placed]
    window: tuple[int, int, int, int]
    selected: list[_Placed]


@dataclass(frozen=True, slots=True)
class _Placed:
    """A footprint and its global index on its lattice."""

    footprint: TileFootprint
    row: int
    col: int

    def meets(self, r0: int, c0: int, r1: int, c1: int) -> bool:
        """Its node index range meets `[r0, r1] x [c0, c1]`, closed on both sides."""
        m = self.footprint.meta
        return (
            self.row <= r1 and self.row + m.rows > r0 and self.col <= c1 and self.col + m.cols > c0
        )


def plan_mosaic(
    footprints: Sequence[TileFootprint], bounds: Bounds | None = None, needed: Any = None
) -> MosaicPlan:
    """Select the tiles the request needs, on one lattice, or refuse (R4).

    `bounds` is the box in the DEM's CRS; without it, the lattice's union.
    `needed` is a shapely geometry in the DEM's CRS, closed: a node on its
    boundary is needed. Without it, every node of the window is.

    When the request selects tiles on several lattices, the plan is on the one
    whose own tiles cover every needed node (Ola's Q5 reading, `_covering`).
    """
    if not footprints:
        raise MosaicError("no tiles to plan a mosaic from")
    groups: list[list[TileFootprint]] = []
    for footprint in sorted(footprints, key=lambda f: f.name):
        group = next((g for g in groups if _aligned(g[0].meta, footprint.meta)), None)
        if group is None:
            groups.append([footprint])
        else:
            group.append(footprint)

    chosen: list[_Lattice] = []
    for group in groups:
        reference, placed = _placed(group)
        window = _window(bounds, reference, group[0].meta, placed)
        selected = [p for p in placed if window is not None and p.meets(*window)]
        if window is not None and selected:
            chosen.append(_Lattice(group, reference, placed, window, selected))
    if not chosen:
        raise MosaicError(f"the box {bounds} meets no tile")
    if len(chosen) > 1:
        return plan_mosaic(_covering(chosen, bounds, needed), bounds, needed)

    lattice = chosen[0]
    (x_ref, y_ref), (r0, c0, r1, c1), selected = lattice.reference, lattice.window, lattice.selected
    rows, cols = r1 - r0 + 1, c1 - c0 + 1
    if rows < 2 or cols < 2:
        raise MosaicError(f"the request is {rows} x {cols} nodes; a mesh needs at least 2 x 2")
    dtype = np.result_type(*(p.footprint.dtype for p in selected))
    cap, size = physical_memory() // 2, rows * cols * dtype.itemsize
    if size > cap:
        raise MosaicError(
            f"the request is {rows} x {cols} nodes, {size} bytes at {dtype}, over the cap "
            f"of half the physical memory ({cap} bytes); narrow it with --bbox"
        )
    vertical = any(p.footprint.meta.vertical_unit_assumed for p in selected)
    meta = _grid(selected[0].footprint.meta, lattice.reference, lattice.window, vertical)
    tiles = tuple(_placement(p, r0, c0, r1, c1) for p in selected)
    uncovered = _uncovered(meta, tiles, needed)
    if uncovered:
        raise MosaicError(uncovered)
    return MosaicPlan(
        meta=meta,
        reference=(x_ref, y_ref),
        window=IndexWindow(row0=r0, col0=c0, rows=rows, cols=cols),
        tiles=tiles,
    )


def assemble(
    plan: MosaicPlan,
    load: Callable[[str], DemTile],
    needed: Any = None,
    *,
    load_window: Callable[[str, IndexWindow], DemTile] | None = None,
) -> Mosaic:
    """Load the plan's tiles one at a time, in plan order, into one canvas (R5).

    One tile whose grid is the mosaic's is returned as loaded: no canvas, no
    copy (I7). Otherwise the canvas is allocated once, NaN-filled, in the
    planned dtype (the placements' result dtype, from the headers), and handed
    to `DemTile` without a copy (R7). A hand-built footprint may under-state
    its dtype, so a loaded tile that promotes further re-casts the canvas.

    Each tile is written whole, then dropped before the next load; what it
    holds where it meets another tile's placement is kept as a strip. Once all
    are in, every overlap is decided from the strips (Ola's Q1 revised,
    `_decide`) and each pair that disagrees is reported (`_seam`).

    `needed`, a shapely geometry in the DEM's CRS as for `plan_mosaic`, limits
    the report to the overlap nodes it covers (closed; Ola, 2026-09-28: with
    `--domain`, the needed region). Only the report: every node is decided as
    without it, and only overlap nodes are tested, never the canvas.

    Memory: the peak is the canvas, the strips, and one load's own peak, which
    is not one tile: decoding a DTM10 tile peaks at 2.0 to 3.0 tiles,
    depending on the tile (the decoder's buffers and `DemTile`'s read-only
    copy). Measured with
    tracemalloc on the 15a acceptance box (nine DTM10 tiles, 404 MB canvas,
    102 MB tiles): 720 MB, the canvas plus 3.1 tiles.

    With `load_window` (23a-1), each placement's source window is loaded
    instead of its whole tile, `load_window(name, placement.source)`, and is
    copied whole; its meta must be `window_meta` of the listed one.
    """
    if len(plan.tiles) == 1 and plan.tiles[0].meta == plan.meta:
        return Mosaic(tile=_loaded(plan.tiles[0], load, load_window), plan=plan)
    if not plan.tiles:
        raise MosaicError("the plan has no tiles")
    planned = np.result_type(*(t.dtype for t in plan.tiles))
    canvas: npt.NDArray[Any] = np.full((plan.meta.rows, plan.meta.cols), np.nan, dtype=planned)
    ordered = sorted(plan.tiles, key=lambda t: t.name)
    strips: dict[tuple[str, str], npt.NDArray[Any]] = {}
    for placement in plan.tiles:
        array = _loaded(placement, load, load_window).array
        if np.result_type(canvas.dtype, array.dtype) != canvas.dtype:
            canvas = canvas.astype(np.result_type(canvas.dtype, array.dtype))
        c, s = placement.canvas, placement.source
        rows, cols = slice(s.row0, s.row0 + s.rows), slice(s.col0, s.col0 + s.cols)
        incoming = array if load_window else array[rows, cols]  # a window comes cut
        canvas[c.row0 : c.row0 + c.rows, c.col0 : c.col0 + c.cols] = incoming
        for other in ordered:
            box = _meet(c, other.canvas)
            if other.name != placement.name and box is not None:
                strips[placement.name, other.name] = incoming[_within(box, c)].copy()
        del array, incoming  # before the next load, or two tiles outlive this one
    seams = []
    if needed is not None:
        shapely.prepare(needed)
    for i, a in enumerate(ordered):
        for b in ordered[i + 1 :]:
            box = _meet(a.canvas, b.canvas)
            if box is None:
                continue
            _decide(canvas, box, ordered, strips, (a.name, b.name), plan.meta.nodata)
            kept = _covered(box, plan.meta, needed)
            first, second = strips[a.name, b.name][kept], strips[b.name, a.name][kept]
            seam = _seam(a.name, b.name, first, second, plan)
            if seam is not None:
                seams.append(seam)
    return Mosaic(tile=DemTile._adopt(plan.meta, canvas), plan=plan, seams=tuple(seams))


def _aligned(a: RasterMeta, b: RasterMeta) -> bool:
    """Same CRS, spacing, registration and NoData, and integer node offsets (R4.1)."""
    return _key(a) == _key(b) and all(abs(v - round(v)) <= ALIGN_TOLERANCE for v in _cells(a, b))


def _key(m: RasterMeta) -> tuple[int, float, float, bool, float | None]:
    return m.epsg, m.delta_x, m.delta_y, m.pixel_is_area, m.nodata


def _cells(a: RasterMeta, b: RasterMeta) -> tuple[float, float]:
    """`b`'s first node from `a`'s, in `a`'s cells east and south."""
    return (b.x_min - a.x_min) / a.delta_x, (a.y_max - b.y_max) / a.delta_y


def _placed(group: list[TileFootprint]) -> tuple[tuple[float, float], list[_Placed]]:
    """The lattice's reference node (R4.2), and each tile's global index from it."""
    x_ref = min(f.meta.x_min for f in group)
    y_ref = max(f.meta.y_max for f in group)
    dx, dy = group[0].meta.delta_x, group[0].meta.delta_y
    return (x_ref, y_ref), [
        _Placed(f, round((y_ref - f.meta.y_max) / dy), round((f.meta.x_min - x_ref) / dx))
        for f in group
    ]


def _window(
    bounds: Bounds | None, reference: tuple[float, float], m: RasterMeta, placed: list[_Placed]
) -> tuple[int, int, int, int] | None:
    """`(r0, c0, r1, c1)`, inclusive: the box snapped outward, clamped to the
    lattice's union (R4.4); None when the box misses the union."""
    last_row = max(p.row + p.footprint.meta.rows for p in placed) - 1
    last_col = max(p.col + p.footprint.meta.cols for p in placed) - 1
    return _clamp(_snapped(bounds, reference, m), (0, 0, last_row, last_col))


def _snapped(
    bounds: Bounds | None, reference: tuple[float, float], m: RasterMeta
) -> tuple[int, int, int, int] | None:
    """The box snapped outward to `m`'s lattice, unclamped; None for no box."""
    if bounds is None:
        return None
    (x_ref, y_ref), dx, dy = reference, m.delta_x, m.delta_y
    return (
        math.floor(_snap((y_ref - bounds.y_max) / dy)),
        math.floor(_snap((bounds.x_min - x_ref) / dx)),
        math.ceil(_snap((y_ref - bounds.y_min) / dy)),
        math.ceil(_snap((bounds.x_max - x_ref) / dx)),
    )


def _clamp(
    window: tuple[int, int, int, int] | None, limits: tuple[int, int, int, int]
) -> tuple[int, int, int, int] | None:
    """`window` (None: unbounded) inside `limits`, or None when they miss."""
    r0, c0, r1, c1 = limits if window is None else window
    r0, c0 = max(r0, limits[0]), max(c0, limits[1])
    r1, c1 = min(r1, limits[2]), min(c1, limits[3])
    return (r0, c0, r1, c1) if r0 <= r1 and c0 <= c1 else None


def _covering(chosen: list[_Lattice], bounds: Bounds | None, needed: Any) -> list[TileFootprint]:
    """Ola's Q5 reading: of the lattices the request selects tiles on, the one
    whose own tiles cover every node it needs; several, the most tiles, ties by
    the first tile's name; none, the Q5 refusal.

    A lattice needs its nodes in the box snapped outward, clamped to the
    bounding box of every selected tile on any lattice, not to its own union,
    or a box running past it into another lattice's tiles would be cut short.
    """
    corners = [(p.footprint.meta, lattice) for lattice in chosen for p in lattice.selected]
    x_lo = min(m.x_min for m, _ in corners)
    x_hi = max(m.x_min + (m.cols - 1) * m.delta_x for m, _ in corners)
    y_hi = max(m.y_max for m, _ in corners)
    y_lo = min(m.y_max - (m.rows - 1) * m.delta_y for m, _ in corners)
    covering = []
    for lattice in chosen:
        (x_ref, y_ref), m = lattice.reference, lattice.group[0].meta
        inside = (
            math.ceil(_snap((y_ref - y_hi) / m.delta_y)),
            math.ceil(_snap((x_lo - x_ref) / m.delta_x)),
            math.floor(_snap((y_ref - y_lo) / m.delta_y)),
            math.floor(_snap((x_hi - x_ref) / m.delta_x)),
        )
        window = _clamp(_snapped(bounds, lattice.reference, m), inside)
        if window is None:
            continue
        tiles = tuple(_placement(p, *window) for p in lattice.placed if p.meets(*window))
        if not _uncovered(_grid(m, lattice.reference, window, False), tiles, needed):
            covering.append(lattice)
    if not covering:
        raise MosaicError(_mixed(chosen[0].selected[0].footprint, chosen[1].selected[0].footprint))
    return min(covering, key=lambda c: (-len(c.group), c.group[0].name)).group


def _grid(
    first: RasterMeta,
    reference: tuple[float, float],
    window: tuple[int, int, int, int],
    vertical_unit_assumed: bool,
) -> RasterMeta:
    """The node grid of `window` on `first`'s lattice."""
    r0, c0, r1, c1 = window
    return RasterMeta(
        x_min=reference[0] + c0 * first.delta_x,
        y_max=reference[1] - r0 * first.delta_y,
        delta_x=first.delta_x,
        delta_y=first.delta_y,
        cols=c1 - c0 + 1,
        rows=r1 - r0 + 1,
        epsg=first.epsg,
        nodata=first.nodata,
        nodata_source=first.nodata_source,
        pixel_is_area=first.pixel_is_area,
        vertical_unit_assumed=vertical_unit_assumed,
    )


def _snap(cells: float) -> float:
    """A box edge within `ALIGN_TOLERANCE` of a node is on it, so float noise
    in a spacing like 0.1 does not grow the window by a line."""
    nearest = round(cells)
    return float(nearest) if abs(cells - nearest) <= ALIGN_TOLERANCE else cells


def _placement(p: _Placed, r0: int, c0: int, r1: int, c1: int) -> TilePlacement:
    """Where tile `p` sits in the canvas, and which part of it is used."""
    m = p.footprint.meta
    rs, re = max(r0, p.row), min(r1, p.row + m.rows - 1)
    cs, ce = max(c0, p.col), min(c1, p.col + m.cols - 1)
    rows, cols = re - rs + 1, ce - cs + 1
    return TilePlacement(
        name=p.footprint.name,
        meta=m,
        canvas=IndexWindow(row0=rs - r0, col0=cs - c0, rows=rows, cols=cols),
        source=IndexWindow(row0=rs - p.row, col0=cs - p.col, rows=rows, cols=cols),
        dtype=p.footprint.dtype,
    )


def _uncovered(meta: RasterMeta, tiles: tuple[TilePlacement, ...], needed: Any) -> str | None:
    """Why a needed node no tile covers refuses the request (R4.5, per node),
    or None. A mask, a count and per-axis reductions: no per-node index or
    coordinate array unless `needed` must be tested node by node (S1)."""
    uncovered = np.ones((meta.rows, meta.cols), dtype=bool)
    for t in tiles:
        c = t.canvas
        uncovered[c.row0 : c.row0 + c.rows, c.col0 : c.col0 + c.cols] = False
    if needed is not None and uncovered.any():
        rows, cols = np.nonzero(uncovered)
        outside = ~shapely.intersects_xy(
            needed, meta.x_min + cols * meta.delta_x, meta.y_max - rows * meta.delta_y
        )
        uncovered[rows[outside], cols[outside]] = False
    count = np.count_nonzero(uncovered)
    if not count:
        return None
    rows_hit = np.flatnonzero(uncovered.any(axis=1))
    cols_hit = np.flatnonzero(uncovered.any(axis=0))
    x0, x1 = (meta.x_min + cols_hit[i] * meta.delta_x for i in (0, -1))
    y1, y0 = (meta.y_max - rows_hit[i] * meta.delta_y for i in (0, -1))
    return (
        f"{count} nodes the request needs are in no tile: x {_num(x0)} to {_num(x1)}, "
        f"y {_num(y0)} to {_num(y1)} (EPSG:{meta.epsg})"
    )


def _loaded(
    placement: TilePlacement,
    load: Callable[[str], DemTile],
    load_window: Callable[[str, IndexWindow], DemTile] | None,
) -> DemTile:
    if load_window is None:
        tile, listed = load(placement.name), placement.meta
    else:
        tile = load_window(placement.name, placement.source)
        listed = window_meta(placement.meta, placement.source)
    if tile.meta != listed:
        raise MosaicError(
            f"{placement.name} changed since it was listed: {tile.meta} is not {listed}"
        )
    return tile


def _meet(a: IndexWindow, b: IndexWindow) -> IndexWindow | None:
    """The canvas nodes both windows hold, or None."""
    r0, c0 = max(a.row0, b.row0), max(a.col0, b.col0)
    r1, c1 = min(a.row0 + a.rows, b.row0 + b.rows), min(a.col0 + a.cols, b.col0 + b.cols)
    if r0 >= r1 or c0 >= c1:
        return None
    return IndexWindow(row0=r0, col0=c0, rows=r1 - r0, cols=c1 - c0)


def _within(inner: IndexWindow, outer: IndexWindow) -> tuple[slice, slice]:
    """`inner`'s slices in an array holding `outer` (both in canvas indices)."""
    r, c = inner.row0 - outer.row0, inner.col0 - outer.col0
    return slice(r, r + inner.rows), slice(c, c + inner.cols)


def _valid(values: npt.NDArray[Any], nodata: float | None) -> npt.NDArray[np.bool_]:
    valid = ~np.isnan(values)
    if nodata is not None:
        valid &= values != nodata
    return valid


def _decide(
    canvas: npt.NDArray[Any],
    box: IndexWindow,
    ordered: list[TilePlacement],
    strips: dict[tuple[str, str], npt.NDArray[Any]],
    pair: tuple[str, str],
    nodata: float | None,
) -> None:
    """Ola's Q1 revised over `box`, the overlap of `pair`: each node takes the
    valid value of the tile it lies deepest in, depth being
    `min(r, c, rows - 1 - r, cols - 1 - c)` in the whole tile's own indices;
    ties to the name that sorts first (`ordered` is by name, and only a deeper
    tile replaces). No valid value: the sentinel if any tile holds it, else
    NaN, in every order. A third tile's values come from its strip with one
    of the pair, which holds every node of `box` it covers."""
    best = np.full((box.rows, box.cols), -1, dtype=np.int64)
    value = np.full((box.rows, box.cols), np.nan, dtype=canvas.dtype)
    sentinel = np.zeros((box.rows, box.cols), dtype=bool)
    for t in ordered:
        part = _meet(t.canvas, box)
        if part is None:
            continue
        other = next(p for p in ordered if p.name == (pair[1] if t.name == pair[0] else pair[0]))
        held = _meet(t.canvas, other.canvas)
        assert held is not None  # it holds `part`
        values = strips[t.name, other.name][_within(part, held)]
        row = np.arange(part.rows) + part.row0 - t.canvas.row0 + t.source.row0
        col = np.arange(part.cols) + part.col0 - t.canvas.col0 + t.source.col0
        depth = np.minimum.outer(
            np.minimum(row, t.meta.rows - 1 - row), np.minimum(col, t.meta.cols - 1 - col)
        )
        here = _within(part, box)
        take = _valid(values, nodata) & (depth > best[here])
        value[here][take] = values[take]
        best[here][take] = depth[take]
        if nodata is not None:
            sentinel[here] |= values == nodata
    if nodata is not None:
        value[(best < 0) & sentinel] = nodata
    canvas[box.row0 : box.row0 + box.rows, box.col0 : box.col0 + box.cols] = value


def _covered(box: IndexWindow, meta: RasterMeta, needed: Any) -> Any:
    """The nodes of `box` that `needed` covers, as a boolean mask of its shape,
    or `...` (every node) without it. A node is `x_min + col * dx`, as in
    `_uncovered`."""
    if needed is None:
        return ...
    rows = (box.row0 + np.arange(box.rows))[:, np.newaxis]
    cols = (box.col0 + np.arange(box.cols))[np.newaxis, :]
    xs, ys = np.broadcast_arrays(meta.x_min + cols * meta.delta_x, meta.y_max - rows * meta.delta_y)
    return shapely.intersects_xy(needed, xs, ys)


def _seam(
    first: str, second: str, a: npt.NDArray[Any], b: npt.NDArray[Any], plan: MosaicPlan
) -> Seam | None:
    """The pair's report over its overlap in the mosaic: nodes where both hold
    a valid value and `|a - b| >= SEAM_THRESHOLD`, and the largest and median
    `|a - b|` over those, in float64. None when no node qualifies. The
    threshold is the report's only: `_decide` never consults it."""
    both = _valid(a, plan.meta.nodata) & _valid(b, plan.meta.nodata)
    gaps = np.abs(a[both].astype(np.float64) - b[both].astype(np.float64))
    gaps = gaps[gaps >= SEAM_THRESHOLD]
    if not gaps.size:
        return None
    return Seam(first, second, int(gaps.size), float(gaps.max()), float(np.median(gaps)))


def _mixed(a: TileFootprint, b: TileFootprint) -> str:
    """Why `b` is not on `a`'s lattice (R4.3, Q5), in the terms that differ."""
    ma, mb = a.meta, b.meta
    if ma.epsg != mb.epsg:
        why = f"is in EPSG:{mb.epsg}, not EPSG:{ma.epsg}"
    elif (ma.delta_x, ma.delta_y) != (mb.delta_x, mb.delta_y):
        why = f"has spacing {mb.delta_x!r} x {mb.delta_y!r}, not {ma.delta_x!r} x {ma.delta_y!r}"
    elif ma.pixel_is_area != mb.pixel_is_area:
        why = f"has registration {_registration(mb)}, not {_registration(ma)}"
    elif ma.nodata != mb.nodata:
        why = f"has nodata {mb.nodata}, not {ma.nodata}"
    else:
        offsets = [
            f"{abs(v - round(v)):.6g} cell ({abs(v - round(v)) * h:.6g} m) {axis}"
            for v, h, axis in zip(
                _cells(ma, mb), (ma.delta_x, ma.delta_y), ("east-west", "north-south"), strict=True
            )
            if abs(v - round(v)) > ALIGN_TOLERANCE
        ]
        why = f"is {' and '.join(offsets)} off {a.name}'s lattice"
    return (
        f"the request selects tiles on two lattices, {a.name} and {b.name}: {b.name} {why}; "
        "tiles are not resampled onto one another, but a --bbox inside one lattice is meshed"
    )


def _registration(m: RasterMeta) -> str:
    return "pixel-is-area" if m.pixel_is_area else "pixel-is-point"


def _num(value: float | np.floating[Any]) -> str:
    """A coordinate written out, never in scientific notation."""
    return f"{value:f}".rstrip("0").rstrip(".")
