"""What each catalogue source holds, read from its headers: `rasputin sources`.

A local source is read through the readers a run uses (`repository_for`), a
remote one through its tile cache (`CacheRepository`) when cached; nothing
is decoded but headers. The figures follow STAC's conventions: coverage as
a longitude/latitude box west, south, east, north; cell size in metres; the
CRS as ``EPSG:n`` where there is one. Kept apart from `sources.py`, which
imports no pyproj or shapely.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Literal

import numpy as np
import shapely
from pydantic import BaseModel, ConfigDict

from tin_engine.cli_kit import words
from tin_engine.crs import crs_label, parse_crs, reprojector
from tin_engine.dem_input import repository_for
from tin_engine.geodesy import LocalFrame
from tin_engine.io.models import TileFootprint
from tin_engine.io.repository import CacheRepository
from tin_engine.sources import (
    RemoteSource,
    ZUnit,
    catalogue,
    user_source_candidates,
    user_source_file,
)

#: Segments per tile edge before reprojection: a chosen value at which an
#: edge's curvature is negligible against a cell.
SEGMENTS = 8
Box = tuple[float, float, float, float]
_FROM_ENTRY = ("id", "kind", "values", "description", "credit", "licence_note", "cite")


class SourceReport(BaseModel):
    """One source as its headers describe it; the figures are None when nothing was read."""

    model_config = ConfigDict(frozen=True)

    id: str
    kind: str
    values: Literal["heights", "classes"]
    path: Path | None
    present: bool
    description: str
    credit: str
    licence_note: str
    cite: tuple[str, ...]
    tiles: int = 0
    crs: str | None = None
    crs_name: str | None = None
    cell_m: tuple[float, float] | None = None
    z_unit: ZUnit | None = None
    z_unit_from: Literal["file", "catalogue", "assumed"] | None = None
    nodata: float | None = None
    coverage: Box | None = None
    box_covered: float | None = None
    problem: str | None = None


class CatalogueReport(BaseModel):
    """The listing: the places read, and every source."""

    model_config = ConfigDict(frozen=True)

    user_file: str | None
    looked_at: tuple[str, ...]
    register_at: str | None
    user_file_problem: str | None
    cache: str | None
    note: str | None
    sources: tuple[SourceReport, ...]


def _figures(
    entry: RemoteSource, tiles: tuple[TileFootprint, ...], box: Box | None
) -> dict[str, object]:
    """The header figures of `tiles`: CRS, cell sizes, height unit, coverage
    and the share of `box` covered."""
    meta = tiles[0].meta
    to_lonlat = reprojector(meta.crs, "EPSG:4326")
    cells, shapes = [], []
    for tile in tiles:
        m = tile.meta
        x0, y0, x1, y1 = m.node_box()
        cells.append(LocalFrame.at(m.crs, ((x0 + x1) / 2, (y0 + y1) / 2)).cell_m(m))
        ring = shapely.segmentize(shapely.box(x0, y0, x1, y1), max(x1 - x0, y1 - y0) / SEGMENTS)
        shapes.append(shapely.Polygon(to_lonlat(np.asarray(ring.exterior.coords))))
    union = shapely.union_all(shapes)
    covered = None
    if box is not None:
        wanted = shapely.envelope(shapely.multipoints([box[:2], box[2:]]))
        # A point or line box has no area: covered or not.
        share = union.covers(wanted) if wanted.area == 0 else union.intersection(wanted).area
        covered = round(float(share) / (wanted.area or 1.0), 4)
    heights = entry.values == "heights"
    stated = "catalogue" if entry.z_unit else "assumed"
    unit_from = stated if meta.vertical_unit_assumed else "file"
    return {
        "crs": crs_label(meta.crs),
        "crs_name": parse_crs(meta.crs).name,
        "cell_m": (round(min(cells), 1), round(max(cells), 1)),
        "z_unit": meta.z_unit if heights else None,
        "z_unit_from": unit_from if heights else None,
        "nodata": meta.nodata,
        "coverage": tuple(round(v, 6) for v in union.bounds),
        "box_covered": covered,
    }


def _one(entry: RemoteSource, cache: Path | None, box: Box | None) -> SourceReport:
    """`entry` with what its files or its cache say; a read that fails is its `problem`."""
    found: dict[str, object] = {}
    if entry.path is not None:
        where: Path | None = entry.path
        present = entry.path.exists()
    else:
        where = None if cache is None else cache / entry.id
        present = where is not None and (where / "manifest.json").is_file()
    if entry.path is not None or present:
        try:
            if entry.path is not None:
                repository = repository_for((entry.path,), z_unit=entry.read_unit)[0]
            else:
                assert cache is not None
                repository = CacheRepository(cache, entry.id, z_unit=entry.read_unit)
            tiles = repository.footprints()
            found = {"tiles": len(tiles), **_figures(entry, tiles, box)}
        except (OSError, ValueError) as exc:
            found["problem"] = words(exc) if isinstance(exc, ValueError) else str(exc)
    named = entry.model_dump(include=set(_FROM_ENTRY))
    return SourceReport.model_validate({**found, **named, "present": present, "path": where})


def report(
    entries: Mapping[str, RemoteSource], cache: Path | None, box: Box | None
) -> tuple[SourceReport, ...]:
    """Every entry, in catalogue order, with its figures; `box` (W S E N
    degrees) adds the share of it each source covers."""
    return tuple(_one(entry, cache, box) for entry in entries.values())


def listing(environ: Mapping[str, str], cache: Path | None, box: Box | None) -> CatalogueReport:
    """The listing for this environment: the user's sources file it finds
    (a refused one is named, and the package's sources are still listed),
    the cache, and every source."""
    looked, user = user_source_candidates(environ), user_source_file(environ)
    problem = None
    try:
        entries = catalogue(user).sources
    except ValueError as exc:
        entries, problem = catalogue(None).sources, words(exc)
    register = user or (looked[0] if looked else None)
    sources = report(entries, cache, box)
    notes = []
    if user is None:
        notes.append(
            f"No user sources file: looked at {', '.join(map(str, looked)) or 'no place'}. "
            f'To use a local DEM, add a [source.<id>] table with kind = "local", path, credit, '
            f"licence_note and description to {register or 'a sources.toml'}."
        )
    if any(s.present and s.kind != "local" for s in sources):
        notes.append(
            "A remote source's figures describe its cached objects, which may hold only "
            "the part fetched."
        )
    return CatalogueReport(
        user_file=None if user is None else str(user),
        looked_at=tuple(map(str, looked)),
        register_at=None if register is None else str(register),
        user_file_problem=problem,
        cache=None if cache is None else str(cache),
        note=" ".join(notes) or None,
        sources=sources,
    )


__all__ = ["CatalogueReport", "SourceReport", "listing", "report"]
