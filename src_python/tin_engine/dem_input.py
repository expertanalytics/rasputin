"""The DEM a run meshes, from a declarative request (increment 15a, R1).

`docs/increments/15-dem-mosaic.md` R1 and R11. A `DemRequest` names the
sources, the box and the caller's NoData; `open_dem` lists the tiles, plans the
mosaic and assembles it. `cli.py` parses flags into a request and calls
`open_dem`; a GUI backend or an API worker builds the same request without
Typer. Besides `io/repository.py`, this is the one module below `cli.py`
with paths, and it only hands them to the repository.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Self

from pydantic import BaseModel, ConfigDict, model_validator

from tin_engine.io.models import DemTile
from tin_engine.io.repository import TiffDemRepository
from tin_engine.mosaic import Bounds, MosaicPlan, Seam, assemble, plan_mosaic


class DemRequest(BaseModel):
    """Exactly one directory of tiles, or one or more tile files (R11)."""

    model_config = ConfigDict(frozen=True)

    sources: tuple[Path, ...]
    bounds: Bounds | None = None
    nodata: float | None = None

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
        return self


@dataclass(frozen=True, slots=True)
class DemInput:
    """The tile to mesh, the plan it was assembled by, a name for the run (the
    directory's name, or the stem of the first file as given), and the
    mosaic's disagreeing seams (Ola's Q1 revised)."""

    tile: DemTile
    plan: MosaicPlan
    label: str
    seams: tuple[Seam, ...] = ()


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
    plan = plan_mosaic(repository.footprints(), request.bounds, None)
    mosaic = assemble(plan, repository.load)
    return DemInput(tile=mosaic.tile, plan=plan, label=label, seams=mosaic.seams)
