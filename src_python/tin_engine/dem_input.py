"""The DEM a run meshes, from a declarative request.

A `DemRequest` names the sources, the box and the caller's NoData; `open_dem`
lists the tiles, plans the mosaic and assembles it. `cli.py` parses flags
into a request and calls `open_dem`; a GUI backend or an API worker builds
the same request without Typer. Besides `io/repository.py`, this is the one
module below `cli.py` with paths, and it only hands them to the repository.

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
from tin_engine.domain import DomainError, DomainPolygon, check_extent
from tin_engine.grow import grow_mitred
from tin_engine.io.models import Bounds, DemTile, RasterMeta, TileFootprint, ZUnit
from tin_engine.io.repository import CacheRepository, DemRepository, TiffDemRepository
from tin_engine.mosaic import MosaicError, MosaicPlan, Seam, assemble, plan_mosaic
from tin_engine.resources import Resources
from tin_engine.target_grid import (
    Block,
    TargetGrid,
    TileWindows,
    check_point_blocks,
    default_spacing,
    resample,
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
    source, and at most one of a box in the DEM's CRS and a domain in its own,
    since a domain's bounds are its box."""

    model_config = ConfigDict(frozen=True)

    sources: tuple[Path, ...] = ()
    cached: CachedSource | None = None
    bounds: Bounds | None = None
    nodata: float | None = None
    #: The caller's height unit for files that do not state one.
    z_unit: ZUnit | None = None
    domain: DomainPolygon | None = None
    target_crs: str | None = None

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
        if self.bounds is not None and self.domain is not None:
            raise ValueError("give a box or a domain, not both: the domain's bounds are the box")
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


def open_dem(request: DemRequest, resources: Resources | None = None) -> DemInput:
    """List, plan and assemble. Every refusal is a `ValueError` (`GeoTiffError`,
    `MosaicError`) or, for a file that cannot be opened, an `OSError`.
    `resources` is the thread count for the decoders and the pools and
    the memory budget for the canvas; None is `Resources.fallback()`."""
    res = Resources.fallback() if resources is None else resources
    repository = repository_for(
        request.sources, request.nodata, request.cached, res.threads, request.z_unit
    )
    footprints = repository.footprints()
    if footprints and _reprojected(request, [f.meta for f in footprints]):
        return _open_reprojected(request, repository, footprints, res)
    if request.domain is None or not footprints:
        plan = plan_mosaic(footprints, request.bounds, None, memory_budget=res.memory_budget)
        domain = grown = None
    else:
        plan, domain, grown = _domain_plan(footprints, request.domain, res.memory_budget)
    # With a domain the seam report counts only the needed region's nodes,
    # the only ones the mesh can read; which value a node takes does not
    # depend on it.
    repository.check(plan)
    mosaic = assemble(plan, repository.load, grown, load_window=repository.load_window)
    return DemInput(tile=mosaic.tile, plan=plan, domain=domain, seams=mosaic.seams)


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


def _open_reprojected(
    request: DemRequest,
    repository: DemRepository,
    footprints: Sequence[TileFootprint],
    res: Resources,
) -> DemInput:
    """The reprojected path: the target grid around the domain (or box) in
    the target CRS, the source region planned and assembled as on the
    projected path, resampled; the check points are a lazy iterator over the
    source."""
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
    plan = plan_mosaic(footprints, box, needed, memory_budget=res.memory_budget)
    repository.check(plan)
    mosaic = assemble(plan, repository.load, needed, load_window=repository.load_window)
    source, threads = TileWindows(mosaic.tile), res.threads
    return DemInput(
        tile=resample(grid, source, threads),
        plan=plan,
        domain=domain if request.domain is not None else None,
        seams=mosaic.seams,
        checks=check_point_blocks(grid, source, domain, threads),
        grid=grid,
        source_crs=source_crs,
    )


def _domain_plan(
    footprints: Sequence[TileFootprint], given: DomainPolygon, budget: int | None = None
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
        plan = plan_mosaic(footprints, box, domain.polygon, memory_budget=budget)
        reach, grown = 0.0, domain.polygon
        while (diagonal := math.hypot(plan.meta.delta_x, plan.meta.delta_y)) > reach:
            reach, grown = diagonal, grow_mitred(domain.polygon, diagonal)
            plan = plan_mosaic(footprints, box, grown, memory_budget=budget)
        if past := _past(box, plan.meta):  # a vertex within 1e-6 cell past a node line
            plan = plan_mosaic(footprints, past, grown, memory_budget=budget)
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
