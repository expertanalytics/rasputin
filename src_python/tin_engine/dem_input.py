"""The DEM a run meshes, from a declarative request.

A `DemRequest` names the sources, the box and the caller's NoData; `open_dem`
lists the tiles, plans the mosaic and assembles it, `assemble_dem` after
`plan_dem`. A cut run plans the whole and reads one piece at a time
(`piece_plan`, `load_piece`). `cli.py` parses flags into a request; a GUI
backend or an API worker builds the same request without Typer. Besides
`io/repository.py`, this is the one module below `cli.py` with paths, and it
only hands them to the repository.

With a domain, the domain is moved into the DEM's CRS first; its bounds take
the box's place, the polygon grown by one cell is the region whose nodes must
be covered, and its extent is checked against the plan, all before any tile is
loaded, so a bad domain costs no decoding.

A cached source is a catalogue key and a cache root, `CachedSource`, read by
`CacheRepository`; its missing blocks are refused before any block is
decoded, and every tile is decoded by window.
"""

from __future__ import annotations

import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

import shapely
from pydantic import BaseModel, ConfigDict, model_validator

from tin_engine.crs import (
    CrsSuggestion,
    crs_label,
    parse_crs,
    reprojector,
    same_crs,
    single_crs,
    suggest_crs,
)
from tin_engine.domain import CrsBox, DomainError, DomainPolygon, check_extent
from tin_engine.grow import grow_mitred
from tin_engine.io.models import Bounds, DemTile, IndexWindow, RasterMeta, TileFootprint, ZUnit
from tin_engine.io.repository import CacheRepository, DemRepository, TiffDemRepository
from tin_engine.mosaic import (
    Mosaic,
    MosaicError,
    MosaicPlan,
    Seam,
    assemble,
    check_canvas,
    plan_mosaic,
    restrict,
)
from tin_engine.resources import Resources
from tin_engine.target_grid import (
    Block,
    HeldWindows,
    SourceWindows,
    TargetGrid,
    TileWindows,
    check_point_blocks,
    default_spacing,
    grid_meta,
    resample,
    source_box,
    source_region,
    target_grid_for,
)


class CachedSource(BaseModel):
    """A catalogue source's key and the cache root holding `<cache>/<source>/`."""

    model_config = ConfigDict(frozen=True)

    source: str
    cache: Path


class DemRequest(BaseModel):
    """Exactly one directory of tiles, or one or more tile files, or a cached
    source, and at most one of `bounds` (a box in the DEM's CRS, or in
    `target_crs` when given), `box` (a box in its own CRS, moved into that
    one) and a domain in its own, since a domain's bounds are its box."""

    model_config = ConfigDict(frozen=True)

    sources: tuple[Path, ...] = ()
    cached: CachedSource | None = None
    bounds: Bounds | None = None
    nodata: float | None = None
    #: The caller's height unit for files that do not state one.
    z_unit: ZUnit | None = None
    domain: DomainPolygon | None = None
    target_crs: str | None = None
    box: CrsBox | None = None

    @model_validator(mode="after")
    def _one_form(self) -> Self:
        if self.cached is not None and self.sources:
            raise ValueError("give DEM sources or a cached source, not both")
        if not self.sources and self.cached is None:
            raise ValueError("no DEM source given")
        directories = sum(p.is_dir() for p in self.sources)
        if directories and len(self.sources) > 1:
            raise ValueError(
                f"give exactly one directory, or one or more files; got {len(self.sources)} "
                f"sources of which {directories} are directories"
            )
        if sum(v is not None for v in (self.bounds, self.box, self.domain)) > 1:
            raise ValueError("give one box or a domain, not two: the domain's bounds are its box")
        return self


@dataclass(frozen=True, slots=True)
class DemInput:
    """The tile to mesh, the plan it was assembled by, the domain in the DEM's
    CRS (None without one), and the mosaic's disagreeing seams, reported
    rather than refused. On the reprojected path `tile` is the target
    grid's, `domain` is in the target CRS, and `checks` yields the source's
    nodes there for the final check, with the grid and the DEM's CRS."""

    tile: DemTile
    plan: MosaicPlan
    domain: DomainPolygon | None = None
    seams: tuple[Seam, ...] = ()
    checks: Iterator[Block] | None = None
    grid: TargetGrid | None = None
    source_crs: str = ""


def repository_for(
    sources: tuple[Path, ...],
    nodata: float | None = None,
    cached: CachedSource | None = None,
    threads: int = 1,
    z_unit: ZUnit | None = None,
) -> TiffDemRepository | CacheRepository:
    """The repository over one directory of tiles, over the files given, or
    over a cached source; `threads` decode its blocks.
    `z_unit` is the caller's height unit."""
    keys: dict[str, Any] = {"nodata": nodata, "threads": threads, "z_unit": z_unit}
    if cached is not None:
        return CacheRepository(cached.cache, cached.source, **keys)
    first = sources[0]
    if first.is_dir():
        return TiffDemRepository.from_directory(first, **keys)
    return TiffDemRepository(sources, **keys)


@dataclass(frozen=True, slots=True)
class DemPlan:
    """A run's DEM from headers alone: the repository, the run's mosaic plan
    (the source's, on the reprojected path), `meta` the meshing lattice's
    grid, the domain in the DEM's CRS (the target CRS when reprojected;
    None without one), and `needed`, the region the overlap report counts.
    Reprojected: the target grid, the domain or box the check points are
    kept for, in the target CRS, and the source's CRS."""

    repository: DemRepository
    mosaic: MosaicPlan
    meta: RasterMeta
    domain: DomainPolygon | None = None
    needed: Any = None
    grid: TargetGrid | None = None
    checked: DomainPolygon | None = None
    source_crs: str = ""


@dataclass(frozen=True, slots=True)
class PiecePlan:
    """One piece's read: `window` in the run window's indices, the plan it
    decodes, and that plan's window of the run's own mosaic (the piece's
    source window when reprojected, with its sub-grid)."""

    window: IndexWindow
    mosaic: MosaicPlan
    source: IndexWindow
    grid: TargetGrid | None = None


@dataclass(frozen=True, slots=True)
class PieceDem:
    """One piece's DEM: its tile, whose node (0, 0) is the run window's
    `(row0, col0)`, the overlaps inside its window that disagree, and on the
    reprojected path its sub-grid and its check points."""

    tile: DemTile
    row0: int
    col0: int
    seams: tuple[Seam, ...] = ()
    grid: TargetGrid | None = None
    checks: Iterator[Block] | None = None


def open_dem(request: DemRequest, resources: Resources | None = None) -> DemInput:
    """List, plan and assemble. Every refusal is a `ValueError` (`GeoTiffError`,
    `MosaicError`) or, for a file that cannot be opened, an `OSError`.
    `resources` is the thread count for the decoders and the pools and
    the memory budget for the canvas; None is `Resources.fallback()`."""
    return assemble_dem(plan_dem(request, resources), resources)


def plan_dem(request: DemRequest, resources: Resources | None = None) -> DemPlan:
    """List the tiles, move the domain and plan the mosaic, and on the
    reprojected path make the target grid; headers only, no decode."""
    res = Resources.fallback() if resources is None else resources
    repository = repository_for(
        request.sources, request.nodata, request.cached, res.threads, request.z_unit
    )
    footprints = repository.footprints()
    reprojected = bool(footprints) and _reprojected(request, [f.meta for f in footprints])
    if request.box is not None and footprints:  # with no tiles the plan refuses as without one
        frame = (
            request.target_crs
            if reprojected
            else single_crs((f.meta.crs for f in footprints), MosaicError)
        )
        bounds = Bounds.of(request.box.to_crs(str(frame)))
        request = request.model_copy(update={"bounds": bounds, "box": None})
    if reprojected:
        return _plan_reprojected(request, repository, footprints)
    if request.domain is None or not footprints:
        plan = plan_mosaic(footprints, request.bounds, None, cap=False)
        domain = grown = None
    else:
        plan, domain, grown = _domain_plan(footprints, request.domain)
    # With a domain the seam report counts only the needed region's nodes,
    # the only ones the mesh can read; which value a node takes does not
    # depend on it.
    return DemPlan(repository, plan, plan.meta, domain, grown)


def assemble_dem(plan: DemPlan, resources: Resources | None = None) -> DemInput:
    """The run's DEM: the memory cap, the repository's check, then the
    decode, and on the reprojected path the resample and the check points."""
    res = Resources.fallback() if resources is None else resources
    check_canvas(plan.mosaic, res.memory_budget)
    plan.repository.check(plan.mosaic)
    mosaic = _assembled(plan, plan.mosaic)
    if plan.grid is None:
        return DemInput(tile=mosaic.tile, plan=plan.mosaic, domain=plan.domain, seams=mosaic.seams)
    source = TileWindows(mosaic.tile)
    return DemInput(
        tile=resample(plan.grid, source, res.threads),
        plan=plan.mosaic,
        domain=plan.domain,
        seams=mosaic.seams,
        checks=_checks(plan, source, res.threads),
        grid=plan.grid,
        source_crs=plan.source_crs,
    )


def piece_plan(plan: DemPlan, window: IndexWindow) -> PiecePlan:
    """What the piece of `window` (in the run window's indices) decodes:
    the run's plan cut to it, or on the reprojected path the target grid cut
    to it and the source window its nodes read, snapped outward to the
    source's lattice and cut to the run's source window."""
    if plan.grid is None:
        return PiecePlan(window, restrict(plan.mosaic, window), window)
    g, m = plan.grid, plan.mosaic.meta
    sub = g.model_copy(update={"row0": g.row0 + window.row0, "col0": g.col0 + window.col0,
                               "rows": window.rows, "cols": window.cols})  # fmt: skip
    box = source_box(sub, m)
    top, left = (m.y_max - box.y_max) / m.delta_y, (box.x_min - m.x_min) / m.delta_x
    bottom, right = (m.y_max - box.y_min) / m.delta_y, (box.x_max - m.x_min) / m.delta_x
    r0, c0 = max(math.floor(top), 0), max(math.floor(left), 0)
    r1, c1 = min(math.ceil(bottom), m.rows - 1), min(math.ceil(right), m.cols - 1)
    source = IndexWindow(row0=r0, col0=c0, rows=r1 - r0 + 1, cols=c1 - c0 + 1)
    return PiecePlan(window, restrict(plan.mosaic, source), source, sub)


def load_piece(plan: DemPlan, piece: PiecePlan, resources: Resources | None = None) -> PieceDem:
    """Decode one piece: its part of the run's mosaic, and on the reprojected
    path its sub-grid resampled from that part, with its check points."""
    res = Resources.fallback() if resources is None else resources
    mosaic = _assembled(plan, piece.mosaic)
    w = piece.window
    if piece.grid is None:
        return PieceDem(mosaic.tile, w.row0, w.col0, mosaic.seams)
    s = piece.source
    source = HeldWindows(plan.mosaic.meta, mosaic.tile, s.row0, s.col0)
    tile = resample(piece.grid, source, res.threads)
    checks = _checks(plan, source, res.threads, s)
    return PieceDem(tile, w.row0, w.col0, mosaic.seams, piece.grid, checks)


def _assembled(plan: DemPlan, mosaic: MosaicPlan) -> Mosaic:
    repository = plan.repository
    return assemble(mosaic, repository.load, plan.needed, load_window=repository.load_window)


def _checks(
    plan: DemPlan, source: SourceWindows, threads: int, within: IndexWindow | None = None
) -> Iterator[Block]:
    assert plan.grid is not None and plan.checked is not None  # the reprojected path's
    return check_point_blocks(plan.grid, source, plan.checked, threads, within=within)


class OutCrsRequiredError(ValueError):
    """A geographic DEM without a target CRS: `head` says why and
    what is suggested; `suggestion.proj` is the `--out-crs` to paste."""

    def __init__(self, head: str, suggestion: CrsSuggestion) -> None:
        super().__init__(f"{head}: --out-crs '{suggestion.proj}'")
        self.head, self.suggestion = head, suggestion


def _reprojected(request: DemRequest, metas: list[RasterMeta]) -> bool:
    """Whether the DEM is resampled onto a target grid: it is geographic, or
    `target_crs` differs from its CRS. A geographic DEM with no target is
    refused, suggesting a CRS for the box, since the mesh is computed in
    metres and the choice of projection is the caller's."""
    first = metas[0]
    if request.target_crs is not None:
        target = parse_crs(request.target_crs)
        units = sorted({a.unit_name for a in target.axis_info})
        if not target.is_projected or any(
            a.unit_conversion_factor != 1.0 for a in target.axis_info
        ):
            raise ValueError(
                f"--out-crs {request.target_crs} is a {target.type_name} with axes in {units}; "
                "the mesh is computed in it, so it must be a projected CRS in metres"
            )
        return not same_crs(target, first.crs)
    if not first.geographic:
        return False
    if request.domain is not None:
        x0, y0, x1, y1 = request.domain.to_crs(first.crs).polygon.bounds
    else:
        boxes = [m.node_box() for m in metas]
        x0, x1 = min(b[0] for b in boxes), max(b[2] for b in boxes)
        y0, y1 = min(b[1] for b in boxes), max(b[3] for b in boxes)
    s = suggest_crs((x0, x1, y0, y1), first.crs)
    head = (
        f"the DEM is geographic ({first.crs}), so --out-crs is required; suggested for this "
        f"box ({s.family}, {s.proj4} on the DEM's datum), worst scale error "
        f"{100 * s.max_scale_error:.3g} %, worst areal error {100 * s.max_areal_error:.3g} % "
        "over the box"
    )
    raise OutCrsRequiredError(head, s)


def _plan_reprojected(
    request: DemRequest, repository: DemRepository, footprints: Sequence[TileFootprint]
) -> DemPlan:
    """The reprojected path: the target grid around the domain (or box) in
    the target CRS, and the source region planned as on the projected path."""
    source_crs = single_crs((f.meta.crs for f in footprints), MosaicError)
    target = crs_label(str(request.target_crs))
    if request.domain is not None:
        domain = request.domain.to_crs(target)
    elif request.bounds is not None:
        b = request.bounds
        domain = DomainPolygon(polygon=shapely.box(b.x_min, b.y_min, b.x_max, b.y_max), crs=target)
    else:
        raise ValueError("--out-crs needs a domain or a box in the target CRS")
    meta = footprints[0].meta
    c = domain.polygon.centroid
    ((ax, ay),) = reprojector(target, meta.crs)([(c.x, c.y)])
    grid, grown = target_grid_for(domain, target, default_spacing(meta, (ax, ay)))
    box, needed = source_region(grid, meta, grown)
    plan = plan_mosaic(footprints, box, needed, cap=False)
    given = domain if request.domain is not None else None
    return DemPlan(
        repository, plan, grid_meta(grid, plan.meta), given, needed, grid, domain, source_crs
    )


def _domain_plan(
    footprints: Sequence[TileFootprint], given: DomainPolygon
) -> tuple[MosaicPlan, DomainPolygon, Any]:
    """The domain in the DEM's CRS, the plan of its bounds with the polygon
    grown by one cell needed, and that grown polygon, the needed
    region. Grown by the cell diagonal with mitred corners, a superset of the
    four nodes bilinear z reads at any point of the polygon. A refusal names
    both CRSs and keeps its type."""
    crs = single_crs((f.meta.crs for f in footprints), MosaicError)
    try:
        domain = given.to_crs(crs)
        box = Bounds.of(domain.polygon.bounds)
        # The cell is the chosen lattice's, not a listed tile's, since the
        # tiles can differ in spacing: plan on the polygon itself, grow by
        # that plan's cell diagonal, and re-plan while the lattice chosen has
        # a larger one. The reach only
        # grows, so this ends, and the region needed is never short of the
        # chosen lattice's own cell.
        plan = plan_mosaic(footprints, box, domain.polygon, cap=False)
        reach, grown = 0.0, domain.polygon
        while (diagonal := math.hypot(plan.meta.delta_x, plan.meta.delta_y)) > reach:
            reach, grown = diagonal, grow_mitred(domain.polygon, diagonal)
            plan = plan_mosaic(footprints, box, grown, cap=False)
        if past := _past(box, plan.meta):  # a vertex within 1e-6 cell past a node line
            plan = plan_mosaic(footprints, past, grown, cap=False)
        check_extent(domain, plan.meta)
    except (DomainError, MosaicError) as exc:
        raise type(exc)(f"the domain, in {given.crs}, in the DEM's {crs}: {exc}") from exc
    return plan, domain, grown


def _past(box: Bounds, m: RasterMeta) -> Bounds | None:
    """`box` with each edge past `m`'s node rectangle moved out by half a cell,
    or None when none is. The mosaic's window snaps an edge within
    `ALIGN_TOLERANCE` cell of a node line onto it; a domain vertex just past
    that line needs the next one, and moving the edge half a cell out takes
    exactly that line."""
    _, y_min, x_max, _ = m.node_box()
    out = (box.x_min < m.x_min, box.y_min < y_min, box.x_max > x_max, box.y_max > m.y_max)
    if not any(out):
        return None
    hx, hy = m.delta_x / 2, m.delta_y / 2
    return Bounds(
        x_min=box.x_min - hx * out[0],
        y_min=box.y_min - hy * out[1],
        x_max=box.x_max + hx * out[2],
        y_max=box.y_max + hy * out[3],
    )
