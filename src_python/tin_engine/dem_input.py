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
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

from pydantic import BaseModel, ConfigDict, model_validator

from tin_engine.domain import DomainError, DomainPolygon, check_extent
from tin_engine.io.models import DemTile, RasterMeta
from tin_engine.io.repository import TiffDemRepository
from tin_engine.mosaic import Bounds, MosaicError, MosaicPlan, assemble, plan_mosaic


class DemRequest(BaseModel):
    """Exactly one directory of tiles, or one or more tile files (R11), and at
    most one of a box in the DEM's CRS and a domain in its own (R6)."""

    model_config = ConfigDict(frozen=True)

    sources: tuple[Path, ...]
    bounds: Bounds | None = None
    nodata: float | None = None
    domain: DomainPolygon | None = None

    @model_validator(mode="after")
    def _one_form(self) -> Self:
        if not self.sources:
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
    """The tile to mesh, the plan it was assembled by, and a name for the run:
    the directory's name, or the stem of the first file as given; and the
    domain in the DEM's CRS, None without one."""

    tile: DemTile
    plan: MosaicPlan
    label: str
    domain: DomainPolygon | None = None


def open_dem(request: DemRequest) -> DemInput:
    """List, plan and assemble. Every refusal is a `ValueError` (`GeoTiffError`,
    `MosaicError`) or, for a file that cannot be opened, an `OSError`."""
    first = request.sources[0]
    if first.is_dir():
        repository = TiffDemRepository.from_directory(first, nodata=request.nodata)
        label = first.resolve().name
    else:
        repository = TiffDemRepository(request.sources, nodata=request.nodata)
        label = first.stem
    footprints = repository.footprints()
    if request.domain is None or not footprints:
        plan, domain = plan_mosaic(footprints, request.bounds, None), None
    else:
        plan, domain = _domain_plan(footprints, request.domain)
    return DemInput(
        tile=assemble(plan, repository.load).tile, plan=plan, label=label, domain=domain
    )


def _domain_plan(footprints: Any, given: DomainPolygon) -> tuple[MosaicPlan, DomainPolygon]:
    """The domain in the DEM's CRS, and the plan of its bounds with the
    polygon grown by one cell needed (R4 point 5). Grown by the cell diagonal
    with mitred corners, a superset of the four nodes bilinear z reads at any
    point of the polygon. A refusal names both CRSs and keeps its type."""
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
    return plan, domain


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
