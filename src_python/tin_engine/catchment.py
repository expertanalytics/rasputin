"""The catchment of a lake, from the DEM, as a reduced outline (increment 22).

`docs/increments/22-auto-catchment.md`, "Data flow", "The seed", "The
window" and "The fine outline". :func:`delineate` moves the seed point and the
lake into the DEM's CRS, plans and assembles a window round them (15a), floods
it (`_core.upstream`), grows the window until the catchment lies clearly
inside it or is cut by the data's edge or NoData, and traces the outline of
the in-nodes, holes filled and other pieces dropped. The fine outline is then
reduced to the request's tolerance by `_core.reduce_ring`, keeping its area.

With a river reach (increment 29, "The window loop and the catchment"), the
pour point is the gauge's placed node on the burnt reach (`burn.py`): stage A
floods from it in each window, burnt afresh, and decides every refusal;
stage B grows on until the catchment of `D`, the first chain node at or past
`U`, is clear, and the sensitivity (`sensitivity.py`) is read from one
`_core.accumulate` of the final window.

No paths: the DEM comes through a `DemRepository`, and the lakes arrive as
shapely geometries the CLI read. Blocking (the flood releases the GIL); an
async caller runs :func:`delineate` in `asyncio.to_thread`.
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable
from dataclasses import dataclass, replace
from typing import Any, Literal, Self

import numpy as np
import numpy.typing as npt
import shapely
from pydantic import BaseModel, ConfigDict, model_validator
from shapely.geometry import Point, Polygon

from tin_engine._core import ReduceStatus, UpstreamOutcome, accumulate, reduce_ring, upstream
from tin_engine.burn import BurnRefusal, burn_reach
from tin_engine.crs import crs_label, parse_crs, reprojector, same_crs, single_crs
from tin_engine.gauge import Reach
from tin_engine.io.models import RasterMeta
from tin_engine.io.repository import DemRepository
from tin_engine.mosaic import (
    Bounds,
    MixedGridError,
    MosaicError,
    MosaicPlan,
    assemble,
    physical_memory,
    plan_mosaic,
)
from tin_engine.outline import trace
from tin_engine.raster import to_core
from tin_engine.sensitivity import Sensitivity, assess

#: Metres round the seed's bounds, and round the catchment's when the window
#: grows; doubled at each growth step.
WINDOW_MARGIN_M = 2000
_SIDES = ("north", "west", "south", "east")


class CatchmentError(ValueError):
    """A catchment that cannot be delineated, in words for the person asking."""


class LakeError(CatchmentError):
    """The lakes cannot seed it: the point in no lake or in two, or a lake
    that does not move into the DEM's CRS."""


class MixedGridRefusal(CatchmentError):  # noqa: N818 -- the name pinned by the suite
    """The window selects tiles on two grids, neither covering it: the
    mosaic's `MixedGridError`, with its `tiles` and words."""

    def __init__(self, message: str, tiles: tuple[str, ...]) -> None:
        super().__init__(message)
        self.tiles = tiles


class CatchmentRequest(BaseModel):
    """The seed point in `seed_crs`, and the lakes (shapely polygons or
    multipolygons in `lakes_crs`) of which the one containing it is the seed.
    Without lakes the seed is the DEM node nearest the point."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    seed: tuple[float, float]
    seed_crs: str = "EPSG:4326"
    lakes: tuple[Any, ...] | None = None
    lakes_crs: str | None = None
    #: Metres the reduced outline may stray from the fine one; None is twice
    #: the DEM's cell, 0 keeps the fine outline less its collinear vertices.
    outline_tolerance: float | None = None
    #: The gauge's mapped reach (increment 29), in `seed_crs`, which must then
    #: be the DEM's: the pour point is the burnt chain's placed node.
    reach: Reach | None = None

    @model_validator(mode="after")
    def _lakes_have_a_crs(self) -> Self:
        if self.lakes is not None and self.lakes_crs is None:
            raise ValueError("lakes need their CRS")
        if self.lakes is not None and self.reach is not None:
            raise ValueError(
                "a river reach and lakes cannot be combined: a lake is already the seed"
            )
        t = self.outline_tolerance
        if t is not None and not (math.isfinite(t) and t >= 0.0):
            raise ValueError(f"the outline tolerance must be finite and at least 0, got {t}")
        return self


@dataclass(frozen=True, slots=True)
class Window:
    """One flood: the window's node extent, its shape, the flood's seconds,
    and the sides it grew on afterwards (none: the catchment is contained)."""

    bounds: Bounds
    rows: int
    cols: int
    seconds: float
    grown: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class GaugeResult:
    """The gauge on the river: the placed node and the burnt chain in the
    DEM's CRS, what the burn did, and the sensitivity at the placed node
    ("The window loop and the catchment"). `causes` joins the sensitivity's
    causes and, last, `direction` when the line runs against the DEM's slope."""

    node: tuple[float, float]
    chain: tuple[tuple[float, float], ...]
    node_offset_m: float
    chain_nodes: int
    lowered_nodes: int
    lowered_max_m: float
    direction_ok: bool
    end_extended_m: float
    end_closed: bool
    downstream_checked: Literal["whole", "partly", "none"]
    sensitivity: Sensitivity
    causes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class Catchment:
    """The fine outline (holes filled) in the DEM's CRS `crs`, the seed point
    there, the lake's area (None without one), node counts, the outline's
    area and what was left out of it, every window, and the last window's
    mask of in-nodes with its `meta`. Then the reduced outline `reduced` (the
    CLI's `reduced_vertices` and `reduced_area_m2` are its vertex count and
    area), the `tolerance` in metres it was reduced to (`outline_tolerance_m`;
    0 keeps the fine outline less its collinear vertices), and the seconds the
    trace and the reduction took."""

    fine: Polygon
    crs: str
    seed: tuple[float, float]
    lake_area: float | None
    nodes: int
    seed_nodes: int
    fine_area: float
    rings_dropped: int
    dropped_nodes: int
    holes_filled: int
    holes_area: float
    windows: tuple[Window, ...]
    mask: npt.NDArray[np.uint8]
    meta: RasterMeta
    trace_seconds: float
    reduced: Polygon
    tolerance: float
    reduce_seconds: float
    gauge: GaugeResult | None = None


#: A flood of one assembled window: its outcome, and what the caller keeps.
Flood = Callable[[Any], tuple[UpstreamOutcome, Any]]


def check_reach_crs(crs: str, repository: DemRepository) -> None:
    """Refuse a river file whose CRS, `crs`, is not the DEM's (its first
    tile's): the one rule `catchment --rivers` and the batch both apply."""
    dem_crs = repository.footprints()[0].meta.crs
    if not same_crs(crs, dem_crs):
        hint = _code_hint(crs, dem_crs)
        raise ValueError(f"the river file's CRS, {crs}, is not the DEM's, {dem_crs}{hint}")


def _code_hint(given: str, dem_crs: str) -> str:
    """D15 b: when the DEM's CRS is exactly an EPSG code and `given` is none,
    the code to write instead; otherwise nothing."""
    code = parse_crs(dem_crs).to_epsg(min_confidence=100)
    if code is None or parse_crs(given).to_epsg(min_confidence=100) is not None:
        return ""
    return f"; if you mean EPSG:{code}, write EPSG:{code}"


def delineate(request: CatchmentRequest, repository: DemRepository) -> Catchment:
    """The catchment of the request's seed over the repository's DEM, or a
    :class:`CatchmentError` (truncated by the data's edge or NoData, no lake
    or two under the point, over the memory cap)."""
    footprints = repository.footprints()
    dem_crs = single_crs((f.meta.crs for f in footprints), CatchmentError)
    ((x, y),) = reprojector(request.seed_crs, dem_crs)([request.seed])
    if not (math.isfinite(x) and math.isfinite(y)):
        raise CatchmentError(f"the seed {request.seed} has no image in {dem_crs}")
    if request.reach is not None:
        if not same_crs(request.seed_crs, dem_crs):
            hint = _code_hint(request.seed_crs, dem_crs)
            raise CatchmentError(f"the river reach must be in the DEM's CRS, {dem_crs}{hint}")
        return _gauged(request, repository, footprints, dem_crs)
    lake = _lake(request, dem_crs)
    x0, y0, x1, y1 = lake.bounds if lake is not None else (x, y, x, y)
    margin = float(WINDOW_MARGIN_M)
    point = (float(x), float(y))

    def flood(tile: Any) -> tuple[UpstreamOutcome, Any]:
        seed = _seed_mask(tile.meta, lake, point)
        return upstream(to_core(tile), seed), seed

    bounds = Bounds(x_min=x0 - margin, y_min=y0 - margin, x_max=x1 + margin, y_max=y1 + margin)
    m, out, seed, windows = _grow(footprints, repository, bounds, flood, False)
    fine = _outline(dem_crs, lake, point, out, m, seed, tuple(windows))
    return _reduced(fine, _tolerance(request, m))


def _tolerance(request: CatchmentRequest, m: RasterMeta) -> float:
    t = request.outline_tolerance
    return 2.0 * max(m.delta_x, m.delta_y) if t is None else t


def _grow(
    footprints: Any, repository: DemRepository, bounds: Bounds, flood: Flood, burnt: bool
) -> tuple[RasterMeta, UpstreamOutcome, Any, list[Window]]:
    """22's window loop: flood, grow until the catchment is clear of the
    window's edge, or refuse (NoData, the data's edge, the memory cap, whose
    per-node figure grows by a burnt copy and `accumulate`'s 10 bytes when
    `burnt`)."""
    plan, margin = _plan(footprints, bounds), float(WINDOW_MARGIN_M)
    windows: list[Window] = []
    while True:
        m = plan.meta
        itemsize = np.result_type(*(t.dtype for t in plan.tiles)).itemsize
        need = m.rows * m.cols * (itemsize * (1 + burnt) + 2 + 10 * burnt)
        cap = physical_memory() // 2
        if need > cap:
            raise CatchmentError(
                f"the window is {m.rows} x {m.cols} nodes, {need} bytes to flood, over the cap "
                f"of half the physical memory ({cap} bytes)"
            )
        tile = assemble(plan, repository.load).tile
        t0 = time.perf_counter()
        out, kept = flood(tile)
        seconds = time.perf_counter() - t0
        if out.touches_nodata:
            raise CatchmentError("the catchment reaches NoData in the DEM: it is truncated")
        if out.nodes_in == 0:
            raise CatchmentError("no seed node has data: the seed lies on NoData")
        # Grow while the in-nodes' bounds plus the base margin reach past the
        # window; the step itself takes the margin doubled at each growth.
        grown = _grown(plan, _plan(footprints, _joined(m, out, float(WINDOW_MARGIN_M))))
        windows.append(Window(Bounds(**_bounds_of(m)), m.rows, m.cols, seconds, grown))
        if not grown:
            break
        margin *= 2
        plan = _plan(footprints, _joined(m, out, margin))
    if out.touches_edge:
        sides = _edge_sides(m, out)
        raise CatchmentError(f"the catchment is cut by the data's {' and '.join(sides)} edge")
    return m, out, kept, windows


def _burnt_flood(reach: Reach, pick: Callable[[Any], int]) -> Flood:
    """Burn the reach into the window, and flood from the chain node `pick`
    chooses; keeps the burnt window, the path and the placed node's seed."""

    def flood(tile: Any) -> tuple[UpstreamOutcome, Any]:
        try:
            burnt, path = burn_reach(tile, reach)
        except BurnRefusal as exc:
            raise CatchmentError(str(exc)) from exc
        seed = np.zeros(burnt.array.shape, dtype=np.uint8)
        seed[tuple(path.chain[path.placed])] = 1
        start = np.zeros_like(seed)
        start[tuple(path.chain[pick(path)])] = 1
        return upstream(to_core(burnt), start), (burnt, path, seed)

    return flood


def _gauged(
    request: CatchmentRequest, repository: DemRepository, footprints: Any, dem_crs: str
) -> Catchment:
    """Stage A floods from the placed node and decides the catchment and every
    refusal; stage B grows on from it until `D` (the first chain node at or
    past `U`) has a catchment clear of the edge, and on its refusal the
    sensitivity is read in stage A's window."""
    reach = request.reach
    assert reach is not None
    line = np.asarray(reach.line, dtype=np.float64)
    pad = reach.corridor + WINDOW_MARGIN_M
    (x0, y0), (x1, y1) = line.min(axis=0) - pad, line.max(axis=0) + pad
    bounds = Bounds(x_min=x0, y_min=y0, x_max=x1, y_max=y1)
    u = reach.uncertainty

    def at_u(path: Any) -> int:
        past = np.nonzero(path.arc >= u - 1e-6)[0]
        if not past.size:
            raise CatchmentError("the burnt chain ends before the uncertainty")
        return int(past[0])

    m, out, (burnt, path, seed), windows = _grow(
        footprints, repository, bounds, _burnt_flood(reach, lambda p: int(p.placed)), True
    )
    try:
        b = _grow(footprints, repository, Bounds(**_bounds_of(m)), _burnt_flood(reach, at_u), True)
    except CatchmentError:
        pass  # stage B never refuses a catchment: read what stage A's window holds
    else:
        # Stage B starts in stage A's last window: one entry per window.
        same = b[3][0].bounds == windows[-1].bounds
        m, (burnt, path, seed), windows = b[0], b[2], [*windows, *b[3][same:]]
        out = upstream(to_core(burnt), seed)
    acc = accumulate(to_core(burnt))
    down = float(np.hypot(*np.diff(line, axis=0).T).sum()) - reach.at
    s = assess(acc.count, acc.reach, acc.flow_to, path, m.delta_x * m.delta_y / 1e6, u, down)
    checked: Literal["whole", "partly", "none"] = "whole"
    if "downstream_unread" in s.causes:
        checked = "partly" if s.checked_down_m > 0 else "none"
    xy = [(m.x_min + c * m.delta_x, m.y_max - r * m.delta_y) for r, c in path.chain.tolist()]
    gauge = GaugeResult(
        node=xy[path.placed],
        chain=tuple(xy),
        node_offset_m=path.node_offset_m,
        chain_nodes=len(xy),
        lowered_nodes=path.lowered_nodes,
        lowered_max_m=path.lowered_max_m,
        direction_ok=path.direction_ok,
        end_extended_m=path.end_extended_m,
        end_closed=path.end_closed,
        downstream_checked=checked,
        sensitivity=s,
        causes=(*s.causes, *(() if path.direction_ok else ("direction",))),
    )
    fine = _outline(dem_crs, None, xy[path.placed], out, m, seed, tuple(windows))
    return replace(_reduced(fine, _tolerance(request, m)), gauge=gauge)


def _lake(request: CatchmentRequest, dem_crs: str) -> Polygon | None:
    """The one lake part containing the seed point, in the DEM's CRS."""
    if request.lakes is None or request.lakes_crs is None:
        return None
    ((px, py),) = reprojector(request.seed_crs, request.lakes_crs)([request.seed])
    point = Point(px, py)
    found = [p for g in request.lakes for p in shapely.get_parts(g) if p.contains(point)]
    where = f"({px!r}, {py!r}) in {request.lakes_crs}"
    if not found:
        raise LakeError(f"the seed point {where} is in no lake")
    if len(found) > 1:
        raise LakeError(f"the seed point {where} is in {len(found)} lakes; give one")
    move = reprojector(request.lakes_crs, dem_crs)
    lake = shapely.transform(found[0], move)
    if not np.isfinite(shapely.get_coordinates(lake)).all():
        raise LakeError(f"the lake at {where} has a vertex with no image in {dem_crs}")
    assert isinstance(lake, Polygon)
    return lake


def _plan(footprints: Any, bounds: Bounds) -> MosaicPlan:
    try:
        return plan_mosaic(footprints, bounds)
    except MixedGridError as exc:
        raise MixedGridRefusal(str(exc), exc.tiles) from exc
    except MosaicError as exc:
        raise CatchmentError(str(exc)) from exc


def _seed_mask(m: RasterMeta, lake: Polygon | None, point: tuple[float, float]) -> Any:
    """The lake's nodes (`shapely.contains_xy`, over its bounding box only),
    or the node nearest the point."""
    seed = np.zeros((m.rows, m.cols), dtype=np.uint8)
    if lake is None:
        r = round((m.y_max - point[1]) / m.delta_y)
        c = round((point[0] - m.x_min) / m.delta_x)
        if not (0 <= r < m.rows and 0 <= c < m.cols):
            raise CatchmentError(f"the seed point {point} is outside the DEM")
        seed[r, c] = 1
        return seed
    x0, y0, x1, y1 = lake.bounds
    r0, r1 = (
        max(0, math.floor((m.y_max - y1) / m.delta_y)),
        min(m.rows - 1, math.ceil((m.y_max - y0) / m.delta_y)),
    )
    c0, c1 = (
        max(0, math.floor((x0 - m.x_min) / m.delta_x)),
        min(m.cols - 1, math.ceil((x1 - m.x_min) / m.delta_x)),
    )
    if r0 > r1 or c0 > c1:
        return seed
    r, c = np.indices((r1 - r0 + 1, c1 - c0 + 1))
    x, y = m.x_min + (c + c0) * m.delta_x, m.y_max - (r + r0) * m.delta_y
    seed[r0 : r1 + 1, c0 : c1 + 1] = shapely.contains_xy(lake, x, y)
    return seed


def _extent(m: RasterMeta) -> tuple[float, float, float, float]:
    return (
        m.x_min,
        m.y_max - (m.rows - 1) * m.delta_y,
        m.x_min + (m.cols - 1) * m.delta_x,
        m.y_max,
    )


def _bounds_of(m: RasterMeta) -> dict[str, float]:
    return dict(zip(("x_min", "y_min", "x_max", "y_max"), _extent(m), strict=True))


def _joined(m: RasterMeta, out: UpstreamOutcome, margin: float) -> Bounds:
    """The window joined with the in-nodes' bounds grown by `margin`."""
    x0, y0, x1, y1 = _extent(m)
    return Bounds(
        x_min=min(m.x_min + out.col_min * m.delta_x - margin, x0),
        y_min=min(m.y_max - out.row_max * m.delta_y - margin, y0),
        x_max=max(m.x_min + out.col_max * m.delta_x + margin, x1),
        y_max=max(m.y_max - out.row_min * m.delta_y + margin, y1),
    )


def _grown(old: MosaicPlan, new: MosaicPlan) -> tuple[str, ...]:
    """The sides on which `new`'s window reaches past `old`'s."""
    a, b = old.window, new.window
    past = (
        b.row0 < a.row0,
        b.col0 < a.col0,
        b.row0 + b.rows > a.row0 + a.rows,
        b.col0 + b.cols > a.col0 + a.cols,
    )
    return tuple(side for side, p in zip(_SIDES, past, strict=True) if p)


def _edge_sides(m: RasterMeta, out: UpstreamOutcome) -> tuple[str, ...]:
    """The window's sides with an in-node in their first two node lines."""
    near = (
        out.row_min <= 1,
        out.col_min <= 1,
        out.row_max + 2 >= m.rows,
        out.col_max + 2 >= m.cols,
    )
    return tuple(side for side, p in zip(_SIDES, near, strict=True) if p)


def _outline(
    dem_crs: str,
    lake: Polygon | None,
    point: tuple[float, float],
    out: UpstreamOutcome,
    m: RasterMeta,
    seed: Any,
    windows: tuple[Window, ...],
) -> Catchment:
    """The outer ring round the seed point (the pour node without a lake),
    holes filled and every other ring dropped, counted."""
    mask = np.asarray(out.mask)
    t0 = time.perf_counter()
    lattice = trace(mask)
    rings = [
        np.column_stack([m.x_min + r[:, 1] * m.delta_x, m.y_max - r[:, 0] * m.delta_y])
        for r in lattice
    ]
    areas = [_signed_area(r) for r in rings]
    if lake is None:
        r, c = np.argwhere(seed)[0]
        point = (m.x_min + c * m.delta_x, m.y_max - r * m.delta_y)
    inside = Point(point)
    around = [k for k, a in enumerate(areas) if a > 0 and Polygon(rings[k]).contains(inside)]
    if not around:
        raise CatchmentError(f"no ring of the catchment contains the seed point {point}")
    chosen = max(around, key=lambda k: areas[k])
    fine = Polygon(rings[chosen])
    shapely.prepare(fine)
    within = [k != chosen and fine.contains(Point(rings[k][0])) for k in range(len(rings))]
    dropped = [k for k, a in enumerate(areas) if a > 0 and k != chosen and not within[k]]
    holes = [k for k, a in enumerate(areas) if a < 0 and within[k]]
    return Catchment(
        fine=fine,
        crs=crs_label(dem_crs),
        seed=point,
        lake_area=None if lake is None else lake.area,
        nodes=int(out.nodes_in),
        seed_nodes=int(np.count_nonzero(seed)),
        fine_area=fine.area,
        rings_dropped=len(dropped),
        dropped_nodes=_nodes_in(mask, [lattice[k] for k in dropped]),
        holes_filled=len(holes),
        holes_area=-sum(areas[k] for k in holes),
        windows=windows,
        mask=mask,
        meta=m,
        trace_seconds=time.perf_counter() - t0,
        reduced=fine,
        tolerance=0.0,
        reduce_seconds=0.0,
    )


def _reduced(fine: Catchment, tolerance: float) -> Catchment:
    """The fine outline reduced to `tolerance` by `_core.reduce_ring`, in
    metres from its lower-left corner, the seed point kept inside."""
    origin = np.array(fine.fine.bounds[:2])
    ring = np.asarray(fine.fine.exterior.coords)[:-1] - origin
    t0 = time.perf_counter()
    out = reduce_ring(ring, tolerance, np.array([fine.seed]) - origin)
    if out.status != ReduceStatus.Ok:
        raise CatchmentError(f"the outline could not be reduced: {out.status.name}")
    reduced = Polygon(np.asarray(out.ring) + origin)
    return replace(
        fine, reduced=reduced, tolerance=tolerance, reduce_seconds=time.perf_counter() - t0
    )


def _signed_area(ring: npt.NDArray[np.float64]) -> float:
    """Shoelace: positive counter-clockwise."""
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def _nodes_in(mask: Any, rings: list[npt.NDArray[np.float64]]) -> int:
    """In-nodes inside the dropped outer rings (in lattice units), for the report."""
    count = 0
    for ring in rings:
        # Only the nodes in the ring's box: its vertices are on half-cells.
        (r0, c0), (r1, c1) = np.ceil(ring.min(axis=0)), np.floor(ring.max(axis=0))
        r, c = np.nonzero(mask[int(r0) : int(r1) + 1, int(c0) : int(c1) + 1])
        count += int(shapely.contains_xy(Polygon(ring), r + r0, c + c0).sum())
    return count
