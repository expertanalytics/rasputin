"""Our catchment against NVE's polygon: agreement, class, summary (increment 29, PR 4).

`docs/increments/29-nve-reference-catchments.md`, "Agreement and classes"
and "PR 4's red step". Both outlines are counted on the DEM's node lattice
(the final window's origin and spacing, run on as far as either polygon
reaches; no DEM is read), a station is classed from those counts and its
gauge's causes, and a batch's rows are summarised in a deterministic model.
Pure: shapely and numpy, no file.
"""

from __future__ import annotations

import math
from bisect import bisect_right
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import shapely
from pydantic import BaseModel, ConfigDict
from shapely.geometry.base import BaseGeometry

from tin_engine.io.models import RasterMeta

#: The overlap bars of `match` and `close`, both ways, inclusive.
MATCH_OVERLAP = 0.95
CLOSE_OVERLAP = 0.80
#: Cells of mean divide offset that still `match` (30 m on DTM10), inclusive.
MATCH_OFFSET_CELLS = 3.0
#: Lattice nodes tested per `contains_xy` call, so memory stays bounded.
_BAND_NODES = 1 << 20

Class = Literal["refused", "uncertain", "match", "close", "miss"]
CLASSES: tuple[Class, ...] = ("refused", "uncertain", "match", "close", "miss")
MATCH_BY = ("overlap", "offset")
UNCERTAIN_CAUSES = (
    "swing",
    "downstream_unread",
    "chain_not_draining",
    "chain_end_open",
    "direction",
)
REFUSAL_CAUSES = ("mixed_grid", "no_river", "other")
BANDS = ("under 10", "10-100", "100-1000", "over 1000")  # km2
TILE_GROUPS = ("1", "2", "3-4", "5+")
SEEDS = ("river", "lake")  # a row's `seeded_by` (PR 4, "Lake gauges")
_MEASURES = ("area_ratio", "nve_in_ours", "ours_in_nve")
_PERCENTILES = {"min": 0, "p10": 10, "p25": 25, "p50": 50, "p75": 75, "p90": 90, "max": 100}


@dataclass(frozen=True, slots=True)
class Agreement:
    """Lattice nodes strictly inside ours, NVE's and both; the two overlaps;
    our area on the lattice (nodes times the cell's area) over NVE's polygon
    area; the mean divide offset in metres (the area between
    the outlines over NVE's perimeter); and the cell side it is judged in."""

    ours: int
    ref: int
    both: int
    nve_in_ours: float
    ours_in_nve: float
    area_ratio: float
    divide_offset_m: float
    cell_m: float


def nodes_inside(polygons: Sequence[BaseGeometry], meta: RasterMeta) -> int:
    """Lattice nodes strictly inside every one of `polygons`, counted in row
    bands over the box their bounds share."""
    b = np.array([p.bounds for p in polygons])
    x0, y0, (x1, y1) = b[:, 0].max(), b[:, 1].max(), b[:, 2:].min(axis=0)
    dx, dy = meta.delta_x, meta.delta_y
    c0, c1 = math.ceil((x0 - meta.x_min) / dx), math.floor((x1 - meta.x_min) / dx)
    r0, r1 = math.ceil((meta.y_max - y1) / dy), math.floor((meta.y_max - y0) / dy)
    if c0 > c1 or r0 > r1:
        return 0
    for p in polygons:
        shapely.prepare(p)
    x = meta.x_min + np.arange(c0, c1 + 1) * dx
    step = max(1, _BAND_NODES // x.size)
    count = 0
    for r in range(r0, r1 + 1, step):
        xx, yy = np.meshgrid(x, meta.y_max - np.arange(r, min(r + step, r1 + 1)) * dy)
        inside = np.ones(xx.shape, dtype=bool)
        for p in polygons:
            inside &= shapely.contains_xy(p, xx, yy)
        count += int(np.count_nonzero(inside))
    return count


def agreement(fine: BaseGeometry, reference: BaseGeometry, meta: RasterMeta) -> Agreement:
    """Our fine outline against NVE's polygon (every part counted), on the
    lattice of `meta`'s origin and spacing."""
    ours, ref, both = (
        nodes_inside([fine], meta),
        nodes_inside([reference], meta),
        nodes_inside([fine, reference], meta),
    )
    cell = meta.delta_x * meta.delta_y
    return Agreement(
        ours=ours,
        ref=ref,
        both=both,
        nve_in_ours=both / ref if ref else 0.0,
        ours_in_nve=both / ours if ours else 0.0,
        area_ratio=ours * cell / reference.area,
        divide_offset_m=(ref + ours - 2 * both) * cell / reference.length,
        cell_m=math.sqrt(cell),
    )


def classify(
    agreement: Agreement | None, gauge: Any, refused: bool = False
) -> tuple[Class | None, str | None]:
    """The class and, for a match, the test it passed. `gauge` is a
    `catchment.GaugeResult` (its joined `causes` make a station uncertain) or
    None; without an agreement only `refused` and `uncertain` are classes."""
    if refused:
        return "refused", None
    if gauge is not None and gauge.causes:
        return "uncertain", None
    if agreement is None:
        return None, None
    a = agreement
    if min(a.nve_in_ours, a.ours_in_nve) >= MATCH_OVERLAP:
        return "match", "overlap"
    if a.divide_offset_m <= MATCH_OFFSET_CELLS * a.cell_m:
        return "match", "offset"
    return ("close" if min(a.nve_in_ours, a.ours_in_nve) >= CLOSE_OVERLAP else "miss"), None


class Measure(BaseModel):
    """The minimum, the 10th to 90th percentiles and the maximum."""

    model_config = ConfigDict(frozen=True)

    min: float
    p10: float
    p25: float
    p50: float
    p75: float
    p90: float
    max: float


class Measures(BaseModel):
    """The three numbers over the scored stations; None where there are none."""

    model_config = ConfigDict(frozen=True)

    area_ratio: Measure | None
    nve_in_ours: Measure | None
    ours_in_nve: Measure | None


class Group(Measures):
    """A size band or a tile group: its rows (refused ones included), their
    classes, and the share uncertain of those assessed (not refused)."""

    stations: int
    classes: dict[str, int]
    uncertain_share: float | None


class KnownRefusals(BaseModel):
    model_config = ConfigDict(frozen=True)

    count: int
    stations: list[str]
    line: str | None


class Summary(BaseModel):
    """A batch in numbers; fixed key lists are written in full, 0 where none."""

    model_config = ConfigDict(frozen=True)

    stations: int
    classes: dict[str, int]
    match_by: dict[str, int]
    scored: Measures
    by_size: dict[str, Group]
    by_tiles: dict[str, Group]
    by_seed: dict[str, Group]
    uncertain_causes: dict[str, int]
    refusal_causes: dict[str, int]
    known_refusals: KnownRefusals


def _counts(keys: Sequence[str], found: Sequence[Any]) -> dict[str, int]:
    tally = Counter(found)
    return {k: tally[k] for k in keys}


def _measures(rows: Sequence[Any]) -> dict[str, Measure | None]:
    scored = [r for r in rows if r.station_class in ("match", "close", "miss")]
    out: dict[str, Measure | None] = {}
    for key in _MEASURES:
        values = [v for r in scored if (v := getattr(r, key)) is not None]
        q = np.percentile(values, list(_PERCENTILES.values())) if values else None
        out[key] = (
            None if q is None else Measure(**dict(zip(_PERCENTILES, map(float, q), strict=True)))
        )
    return out


def _group(rows: Sequence[Any]) -> Group:
    assessed = [r for r in rows if r.station_class != "refused"]
    uncertain = sum(r.station_class == "uncertain" for r in assessed)
    return Group(
        stations=len(rows),
        classes=_counts(CLASSES, [r.station_class for r in rows]),
        uncertain_share=uncertain / len(assessed) if assessed else None,
        **_measures(rows),
    )


def _band(row: Any) -> str | None:
    """NVE's polygon area decides, or our fine area without a reference; a
    band's lower bound is in it."""
    area = row.reference_area_km2 if row.reference_area_km2 is not None else row.fine_area_km2
    return None if area is None else BANDS[bisect_right((10.0, 100.0, 1000.0), area)]


def _tile_group(row: Any) -> str | None:
    return None if row.tiles is None else TILE_GROUPS[bisect_right((2, 3, 5), row.tiles)]


def known_line(n: int) -> str | None:
    """The known-refusal line for `n` mixed-grid refusals; None for none."""
    if not n:
        return None
    one = n == 1
    return (
        f"{n} station{'' if one else 's'} refused because "
        f"{'its window selects' if one else 'their windows select'} tiles on two different "
        "grids, which rasputin does not combine, and neither grid covers the window alone "
        f"({'a known refusal, not a failure' if one else 'known refusals, not failures'})"
    )


def summarise(rows: Sequence[Any]) -> Summary:
    """The summary of a batch's rows, read by attribute (`StationResult`'s
    names), in row order."""
    known = [
        r.station for r in rows if r.station_class == "refused" and r.refusal_cause == "mixed_grid"
    ]
    return Summary(
        stations=len(rows),
        classes=_counts(CLASSES, [r.station_class for r in rows]),
        match_by=_counts(MATCH_BY, [r.match_by for r in rows]),
        scored=Measures(**_measures(rows)),
        by_size={b: _group([r for r in rows if _band(r) == b]) for b in BANDS},
        by_tiles={g: _group([r for r in rows if _tile_group(r) == g]) for g in TILE_GROUPS},
        by_seed={k: _group([r for r in rows if r.seeded_by == k]) for k in SEEDS},
        uncertain_causes=_counts(
            UNCERTAIN_CAUSES, [c for r in rows if r.station_class == "uncertain" for c in r.causes]
        ),
        refusal_causes=_counts(
            REFUSAL_CAUSES, [r.refusal_cause for r in rows if r.station_class == "refused"]
        ),
        known_refusals=KnownRefusals(count=len(known), stations=known, line=known_line(len(known))),
    )


__all__ = ["Agreement", "Class", "Summary", "agreement", "classify", "nodes_inside", "summarise"]
