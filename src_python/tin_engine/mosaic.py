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

from tin_engine.io.models import DemTile, RasterMeta

if TYPE_CHECKING:
    from tin_engine.io.repository import TileFootprint

#: Cells. Two tiles share a lattice when their node offsets are integers to
#: within this. Measured noise is at most 2.2e-11 cells (B3), and 0 for DTM10.
ALIGN_TOLERANCE = 1e-6

#: The cap counts float32: headers do not give the decoded dtype, so a float64
#: mosaic is under-counted by 2x (pinned by the red suite).
BYTES_PER_NODE = 4


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


class IndexWindow(BaseModel):
    """`rows x cols` nodes starting at `(row0, col0)`."""

    model_config = ConfigDict(frozen=True)

    row0: int
    col0: int
    rows: int
    cols: int


class TilePlacement(BaseModel):
    """One selected tile: `canvas` in canvas indices, `source` in the tile's own."""

    model_config = ConfigDict(frozen=True)

    name: str
    meta: RasterMeta
    canvas: IndexWindow
    source: IndexWindow


class MosaicPlan(BaseModel):
    """What `assemble` will build. `window` is the mosaic's first node as a
    global index, and its shape; `tiles` is sorted by name."""

    model_config = ConfigDict(frozen=True)

    meta: RasterMeta
    reference: tuple[float, float]
    window: IndexWindow
    tiles: tuple[TilePlacement, ...]


@dataclass(frozen=True, slots=True)
class Mosaic:
    """The assembled tile, and the plan it came from, as provenance."""

    tile: DemTile
    plan: MosaicPlan


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
    """
    if not footprints:
        raise MosaicError("no tiles to plan a mosaic from")
    lattices: list[list[TileFootprint]] = []
    for footprint in sorted(footprints, key=lambda f: f.name):
        group = next((g for g in lattices if _aligned(g[0].meta, footprint.meta)), None)
        if group is None:
            lattices.append([footprint])
        else:
            group.append(footprint)

    chosen: list[tuple[tuple[float, float], tuple[int, int, int, int], list[_Placed]]] = []
    for group in lattices:
        reference, placed = _placed(group)
        window = _window(bounds, reference, group[0].meta, placed)
        selected = [p for p in placed if window is not None and p.meets(*window)]
        if window is not None and selected:
            chosen.append((reference, window, selected))
    if not chosen:
        raise MosaicError(f"the box {bounds} meets no tile")
    if len(chosen) > 1:
        raise MosaicError(_mixed(chosen[0][2][0].footprint, chosen[1][2][0].footprint))

    (x_ref, y_ref), (r0, c0, r1, c1), selected = chosen[0]
    rows, cols = r1 - r0 + 1, c1 - c0 + 1
    first = selected[0].footprint.meta
    if rows < 2 or cols < 2:
        raise MosaicError(f"the request is {rows} x {cols} nodes; a mesh needs at least 2 x 2")
    cap = physical_memory() // 2
    if rows * cols * BYTES_PER_NODE > cap:
        raise MosaicError(
            f"the request is {rows} x {cols} nodes, {rows * cols * BYTES_PER_NODE} bytes at "
            f"float32, over the cap of half the physical memory ({cap} bytes); "
            "narrow it with --bbox"
        )
    meta = RasterMeta(
        x_min=x_ref + c0 * first.delta_x,
        y_max=y_ref - r0 * first.delta_y,
        delta_x=first.delta_x,
        delta_y=first.delta_y,
        cols=cols,
        rows=rows,
        epsg=first.epsg,
        nodata=first.nodata,
        nodata_source=first.nodata_source,
        pixel_is_area=first.pixel_is_area,
        vertical_unit_assumed=any(p.footprint.meta.vertical_unit_assumed for p in selected),
    )
    tiles = tuple(_placement(p, r0, c0, r1, c1) for p in selected)
    _check_coverage(meta, tiles, needed)
    return MosaicPlan(
        meta=meta,
        reference=(x_ref, y_ref),
        window=IndexWindow(row0=r0, col0=c0, rows=rows, cols=cols),
        tiles=tiles,
    )


def assemble(plan: MosaicPlan, load: Callable[[str], DemTile]) -> Mosaic:
    """Load the plan's tiles one at a time, in plan order, into one canvas (R5).

    One tile whose grid is the mosaic's is returned as loaded: no canvas, no
    copy (I7). Otherwise the canvas is allocated once, NaN-filled, in the
    first tile's dtype, and handed to `DemTile` without a copy (R7). Dtypes
    are known only after a load, so a later float64 tile re-casts it once to
    the tiles' result dtype; both archives are float32 throughout.
    """
    if len(plan.tiles) == 1 and plan.tiles[0].meta == plan.meta:
        return Mosaic(tile=_loaded(plan.tiles[0], load), plan=plan)
    canvas: npt.NDArray[Any] | None = None
    merged: list[TilePlacement] = []
    for placement in plan.tiles:
        array = _loaded(placement, load).array
        if canvas is None:
            canvas = np.full((plan.meta.rows, plan.meta.cols), np.nan, dtype=array.dtype)
        elif np.result_type(canvas.dtype, array.dtype) != canvas.dtype:
            canvas = canvas.astype(np.result_type(canvas.dtype, array.dtype))
        _merge(canvas, placement, array, plan.meta.nodata, merged)
        merged.append(placement)
    if canvas is None:
        raise MosaicError("the plan has no tiles")
    return Mosaic(tile=DemTile._adopt(plan.meta, canvas), plan=plan)


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
    if bounds is None:
        return 0, 0, last_row, last_col
    (x_ref, y_ref), dx, dy = reference, m.delta_x, m.delta_y
    c0 = max(math.floor(_snap((bounds.x_min - x_ref) / dx)), 0)
    c1 = min(math.ceil(_snap((bounds.x_max - x_ref) / dx)), last_col)
    r0 = max(math.floor(_snap((y_ref - bounds.y_max) / dy)), 0)
    r1 = min(math.ceil(_snap((y_ref - bounds.y_min) / dy)), last_row)
    return (r0, c0, r1, c1) if r0 <= r1 and c0 <= c1 else None


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
    )


def _check_coverage(meta: RasterMeta, tiles: tuple[TilePlacement, ...], needed: Any) -> None:
    """Refuse a needed node no selected tile covers (R4.5, per node)."""
    uncovered = np.ones((meta.rows, meta.cols), dtype=bool)
    for t in tiles:
        c = t.canvas
        uncovered[c.row0 : c.row0 + c.rows, c.col0 : c.col0 + c.cols] = False
    rows, cols = np.nonzero(uncovered)
    xs = meta.x_min + cols * meta.delta_x
    ys = meta.y_max - rows * meta.delta_y
    if needed is not None:
        inside = shapely.intersects_xy(needed, xs, ys)
        xs, ys = xs[inside], ys[inside]
    if len(xs):
        raise MosaicError(
            f"{len(xs)} nodes the request needs are in no tile: x {_num(xs.min())} to "
            f"{_num(xs.max())}, y {_num(ys.min())} to {_num(ys.max())} (EPSG:{meta.epsg})"
        )


def _loaded(placement: TilePlacement, load: Callable[[str], DemTile]) -> DemTile:
    tile = load(placement.name)
    if tile.meta != placement.meta:
        raise MosaicError(
            f"{placement.name} changed since it was listed: {tile.meta} is not {placement.meta}"
        )
    return tile


def _merge(
    canvas: npt.NDArray[Any],
    placement: TilePlacement,
    array: npt.NDArray[Any],
    nodata: float | None,
    merged: list[TilePlacement],
) -> None:
    """R5's overlap rule, per node: valid beats NoData; two valid values must
    be equal (`==`); NoData against NoData keeps the sentinel over NaN, in
    every order."""
    c, s = placement.canvas, placement.source
    region = canvas[c.row0 : c.row0 + c.rows, c.col0 : c.col0 + c.cols]
    incoming = array[s.row0 : s.row0 + s.rows, s.col0 : s.col0 + s.cols]
    old_nan, new_nan = np.isnan(region), np.isnan(incoming)
    old_valid, new_valid = ~old_nan, ~new_nan
    if nodata is not None:
        old_valid &= region != nodata
        new_valid &= incoming != nodata
    clash = old_valid & new_valid & (region != incoming)
    if clash.any():
        row, col = (int(v) for v in np.argwhere(clash)[0])
        other = next(
            m.name
            for m in merged
            if m.canvas.row0 <= c.row0 + row < m.canvas.row0 + m.canvas.rows
            and m.canvas.col0 <= c.col0 + col < m.canvas.col0 + m.canvas.cols
        )
        largest = float(np.max(np.abs(region[clash].astype(np.float64) - incoming[clash])))
        raise MosaicError(
            f"{other} and {placement.name} disagree at {int(clash.sum())} overlapping "
            f"nodes, by up to {largest!r}"
        )
    write = (old_nan & ~new_nan) | (~old_nan & ~old_valid & new_valid)
    region[write] = incoming[write]


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
