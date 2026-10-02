"""The DEM a run meshes, from a declarative request (increment 15a, R1).

`docs/increments/15-dem-mosaic.md` R1 and R11. A `DemRequest` names the
sources, the box and the caller's NoData; `open_dem` lists the tiles, plans the
mosaic and assembles it. `cli.py` parses flags into a request and calls
`open_dem`; a GUI backend or an API worker builds the same request without
Typer. Besides `io/repository.py`, this is the one module below `cli.py`
with paths, and it only hands them to the repository.

With a domain (increment 15b, R4 point 5, R6 and R9), the domain is moved into
the DEM's CRS first; its bounds take the box's place, the polygon grown by one
cell is the region whose nodes must be covered, and its extent is checked
against the plan, all before any tile is loaded (I6).

A cached source (increment 23a-1) is a catalogue key and a cache root,
`CachedSource`, read by `CacheRepository`; its missing blocks are refused
before any block is decoded, and every tile is decoded by window.
"""

from __future__ import annotations

import math
import os
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

import shapely
from pydantic import BaseModel, ConfigDict, model_validator

from tin_engine.crs import crs_label, parse_crs, reprojector, suggest_crs
from tin_engine.domain import DomainError, DomainPolygon, check_extent
from tin_engine.io.models import DemTile, RasterMeta
from tin_engine.io.repository import CacheRepository, TiffDemRepository
from tin_engine.mosaic import Bounds, MosaicError, MosaicPlan, Seam, assemble, plan_mosaic
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
    """Exactly one directory of tiles, or one or more tile files (R11), or a
    cached source (23a-1), and at most one of a box in the DEM's CRS and a
    domain in its own (R6)."""

    model_config = ConfigDict(frozen=True)

    sources: tuple[Path, ...] = ()
    cached: CachedSource | None = None
    bounds: Bounds | None = None
    nodata: float | None = None
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
    """The tile to mesh, the plan it was assembled by, a name for the run (the
    directory's name, or the stem of the first file as given), the domain in
    the DEM's CRS (None without one), and the mosaic's disagreeing seams (Ola's
    Q1 revised). On the reprojected path (15c-2, D1) `tile` is the target
    grid's, `domain` is in the target CRS, and `checks` yields the source's
    nodes there for the final check, with the grid and the DEM's CRS."""

    tile: DemTile
    plan: MosaicPlan
    label: str
    domain: DomainPolygon | None = None
    seams: tuple[Seam, ...] = ()
    checks: Iterator[Block] | None = None
    grid: TargetGrid | None = None
    source_crs: str = ""


def repository_for(
    sources: tuple[Path, ...], nodata: float | None = None, cached: CachedSource | None = None
) -> tuple[TiffDemRepository | CacheRepository, str]:
    """The repository over one directory of tiles, over the files given, or
    over a cached source, and the run's name: the directory's name, the first
    file's stem, or the source's key."""
    if cached is not None:
        return CacheRepository(cached.cache, cached.source, nodata=nodata), cached.source
    first = sources[0]
    if first.is_dir():
        return TiffDemRepository.from_directory(first, nodata=nodata), first.resolve().name
    return TiffDemRepository(sources, nodata=nodata), first.stem


def open_dem(request: DemRequest) -> DemInput:
    """List, plan and assemble. Every refusal is a `ValueError` (`GeoTiffError`,
    `MosaicError`) or, for a file that cannot be opened, an `OSError`."""
    repository, label = repository_for(request.sources, request.nodata, request.cached)
    footprints = repository.footprints()
    if footprints and _reprojected(request, [f.meta for f in footprints]):
        return _open_reprojected(request, repository, footprints, label)
    if request.domain is None or not footprints:
        plan, domain, grown = plan_mosaic(footprints, request.bounds, None), None, None
    else:
        plan, domain, grown = _domain_plan(footprints, request.domain)
    # With a domain the seam report counts only the needed region's nodes
    # (Ola, 2026-09-28); which value a node takes does not depend on it.
    repository.check(plan)
    mosaic = assemble(plan, repository.load, grown, load_window=repository.load_window)
    return DemInput(tile=mosaic.tile, plan=plan, label=label, domain=domain, seams=mosaic.seams)


def _reprojected(request: DemRequest, metas: list[RasterMeta]) -> bool:
    """Whether the DEM is resampled onto a target grid (D1, Q12 (a)): it is
    geographic, or `target_crs` differs from its CRS. A geographic DEM with
    no target is refused, suggesting a CRS for the box (Q11 (a), D8)."""
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
        return target != parse_crs(first.crs)
    if not first.geographic:
        return False
    if request.domain is not None:
        x0, y0, x1, y1 = request.domain.to_crs(first.crs).polygon.bounds
    else:
        x0, x1 = min(m.x_min for m in metas), max(m.x_min + (m.cols - 1) * m.delta_x for m in metas)
        y0, y1 = min(m.y_max - (m.rows - 1) * m.delta_y for m in metas), max(m.y_max for m in metas)
    s = suggest_crs((x0, x1, y0, y1), first.crs)
    raise ValueError(
        f"the DEM is geographic ({first.crs}), so --out-crs is required; suggested for this "
        f"box ({s.family}, {s.proj4} on the DEM's datum): --out-crs '{s.proj}', worst scale error "
        f"{100 * s.max_scale_error:.3g} %, worst areal error {100 * s.max_areal_error:.3g} % "
        "over the box"
    )


def _open_reprojected(
    request: DemRequest, repository: Any, footprints: Any, label: str
) -> DemInput:
    """D1's reprojected path: the target grid around the domain (or box) in
    the target CRS, the source region planned and assembled as 15a does,
    resampled; the check points are a lazy iterator over the source."""
    crss = sorted({f.meta.crs for f in footprints})
    if len(crss) > 1:
        raise MosaicError(f"the tiles are in {len(crss)} CRSs, {crss}; need one")
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
    grid = target_grid_for(domain, target, default_spacing(meta, (ax, ay)))
    grown = domain.polygon.buffer(math.sqrt(2) * grid.spacing, join_style="mitre")
    box, needed = source_region(grid, meta, grown)
    plan = plan_mosaic(footprints, box, needed)
    repository.check(plan)
    mosaic = assemble(plan, repository.load, needed, load_window=repository.load_window)
    source, threads = TileWindows(mosaic.tile), os.cpu_count() or 1
    return DemInput(
        tile=resample(grid, source, threads),
        plan=plan,
        label=label,
        domain=domain if request.domain is not None else None,
        seams=mosaic.seams,
        checks=check_point_blocks(grid, source, domain, threads),
        grid=grid,
        source_crs=crss[0],
    )


def _domain_plan(footprints: Any, given: DomainPolygon) -> tuple[MosaicPlan, DomainPolygon, Any]:
    """The domain in the DEM's CRS, the plan of its bounds with the polygon
    grown by one cell needed (R4 point 5), and that grown polygon, the needed
    region. Grown by the cell diagonal with mitred corners, a superset of the
    four nodes bilinear z reads at any point of the polygon. A refusal names
    both CRSs and keeps its type."""
    epsgs = sorted({f.meta.epsg for f in footprints})
    if len(epsgs) > 1:
        raise MosaicError(f"the tiles are in {len(epsgs)} CRSs, EPSG:{epsgs}; a domain needs one")
    try:
        domain = given.to_crs(f"EPSG:{epsgs[0]}")
        box = Bounds(
            **dict(zip(("x_min", "y_min", "x_max", "y_max"), domain.polygon.bounds, strict=True))
        )
        # The cell is the chosen lattice's, not a listed tile's (15b review, B1):
        # plan on the polygon itself, grow by that plan's cell diagonal, and
        # re-plan while the lattice chosen has a larger one. The reach only
        # grows, so this ends, and the region needed is never short of the
        # chosen lattice's own cell.
        plan, reach, grown = plan_mosaic(footprints, box, domain.polygon), 0.0, domain.polygon
        while (diagonal := math.hypot(plan.meta.delta_x, plan.meta.delta_y)) > reach:
            reach, grown = diagonal, domain.polygon.buffer(diagonal, join_style="mitre")
            plan = plan_mosaic(footprints, box, grown)
        if past := _past(box, plan.meta):  # a vertex within 1e-6 cell past a node line
            plan = plan_mosaic(footprints, past, grown)
        check_extent(domain, plan.meta)
    except (DomainError, MosaicError) as exc:
        raise type(exc)(f"the domain, in {given.crs}, in the DEM's EPSG:{epsgs[0]}: {exc}") from exc
    return plan, domain, grown


def _past(box: Bounds, m: RasterMeta) -> Bounds | None:
    """`box` with each edge past `m`'s node rectangle moved out by half a cell,
    or None when none is. 15a's window snaps an edge within `ALIGN_TOLERANCE`
    cell of a node line onto it; a domain vertex just past that line needs the
    next one, and moving the edge half a cell out takes exactly that line."""
    x_max, y_min = m.x_min + (m.cols - 1) * m.delta_x, m.y_max - (m.rows - 1) * m.delta_y
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
