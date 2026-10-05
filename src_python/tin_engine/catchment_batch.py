"""A catchment per station, one at a time, each compared with its reference (increment 29, PR 4).

`docs/increments/29-nve-reference-catchments.md`, "The batch" and "PR 4's
red step". Each station is placed on the river lines (`gauge.place`), its
catchment delineated in a worker thread (`asyncio.to_thread`: a large
flood takes gigabytes, so never two at once), compared with its reference
polygon when there is one (`reference.py`), and handed to the sink as a
`StationResult` row, in the file's order. Only a `CatchmentError` is a
refusal: it becomes a `refused` row and the batch goes on; any other
exception is a bug and stops it. No paths: the sink writes.
"""

from __future__ import annotations

import asyncio
import time
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from typing import Any, Literal, Protocol

from pydantic import BaseModel, ConfigDict
from shapely.geometry import box
from shapely.geometry.base import BaseGeometry

from tin_engine.catchment import (
    Catchment,
    CatchmentError,
    CatchmentRequest,
    GaugeResult,
    MixedGridRefusal,
    delineate,
)
from tin_engine.crs import reprojector
from tin_engine.gauge import Gauge, place
from tin_engine.io.repository import DemRepository
from tin_engine.io.rivers import RiverSegment
from tin_engine.io.station_set import Station
from tin_engine.reference import Class, Summary, agreement, classify, nodes_inside, summarise


class BatchRequest(BaseModel):
    """The placement's map radius and reach above the gauge, in metres, the
    outline tolerance (None: twice the DEM's cell), and the station numbers
    to run (empty: all)."""

    model_config = ConfigDict(frozen=True)

    map_radius: float = 500.0
    reach_up: float = 1000.0
    outline_tolerance: float | None = None
    only: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class StationResult:
    """One station's row ("The batch"; the CSV's columns, in this order). A
    refused row has None in every field it could not have: a refusal by
    `delineate` keeps its placement, one by `place` has none. Areas in km2;
    `tiles` counts the repository's tiles whose node box meets our fine
    outline; `grid_tiles` is one tile from each grid, for `mixed_grid`."""

    station: str
    name: str | None
    station_class: Class | None
    match_by: str | None = None
    refusal_message: str | None = None
    refusal_cause: Literal["mixed_grid", "no_river", "other"] | None = None
    grid_tiles: tuple[str, ...] | None = None
    # The placement.
    placed_on: str | None = None
    distance_m: float | None = None
    uncertainty_m: float | None = None
    elvid: str | None = None
    objectid: int | None = None
    lake: bool | None = None
    confluence_near: bool | None = None
    reach_up_m: float | None = None
    reach_down_m: float | None = None
    reach_fork: bool | None = None
    # The gauge and the burn.
    node_offset_m: float | None = None
    chain_nodes: int | None = None
    lowered_nodes: int | None = None
    lowered_max_m: float | None = None
    direction_ok: bool | None = None
    end_extended_m: float | None = None
    end_closed: bool | None = None
    downstream_checked: str | None = None
    # The sensitivity, and the joined causes.
    a0: float | None = None
    area_up: float | None = None
    area_down: float | None = None
    swing: float | None = None
    largest_step: float | None = None
    largest_step_at_m: float | None = None
    checked_up_m: float | None = None
    checked_down_m: float | None = None
    drains: bool | None = None
    monotone: bool | None = None
    causes: tuple[str, ...] = ()
    # The catchment and the agreement.
    nodes: int | None = None
    fine_area_km2: float | None = None
    reduced_area_km2: float | None = None
    reference_area_km2: float | None = None
    nve_area_km2: float | None = None
    nve_in_ours: float | None = None
    ours_in_nve: float | None = None
    area_ratio: float | None = None
    divide_offset_m: float | None = None
    tiles: int | None = None
    windows: int | None = None
    seconds: float | None = None


class BatchSink(Protocol):
    """Where the batch's products go: a station's catchment, always before its
    row and never for a refused station, then the row."""

    def catchment(self, station: Station, result: Catchment) -> None: ...

    def row(self, row: StationResult) -> None: ...


def _gauge_fields(gauge: GaugeResult) -> dict[str, Any]:
    skip = ("node", "chain", "sensitivity", "causes")
    own = {k: getattr(gauge, k) for k in GaugeResult.__slots__ if k not in skip}
    s = asdict(gauge.sensitivity)
    del s["causes"], s["well_posed"]
    return {**own, **s, "causes": gauge.causes}


async def run_batch(
    request: BatchRequest,
    repository: DemRepository,
    stations: Sequence[Station],
    stations_crs: str,
    segments: Sequence[RiverSegment],
    segments_crs: str,
    references: Mapping[str, BaseGeometry] | None,
    sink: BatchSink,
) -> Summary:
    """Run each station in `stations` (or those in `request.only`, in file
    order); `references` are in the river file's CRS. Returns the summary."""
    unknown = set(request.only) - {s.station for s in stations}
    if unknown:
        raise ValueError(f"not in the stations file: {', '.join(sorted(unknown))}")
    move = reprojector(stations_crs, segments_crs)
    boxes = [
        box(
            m.x_min, m.y_max - (m.rows - 1) * m.delta_y, m.x_min + (m.cols - 1) * m.delta_x, m.y_max
        )
        for m in (f.meta for f in repository.footprints())
    ]
    rows: list[StationResult] = []
    for station in stations:
        if request.only and station.station not in request.only:
            continue
        t0 = time.perf_counter()
        ref = None if references is None else references.get(station.station)
        fields: dict[str, Any] = {
            "station": station.station,
            "name": station.name,
            "nve_area_km2": station.nve_area_km2,
            "reference_area_km2": None if ref is None else ref.area / 1e6,
        }
        ((x, y),) = move([(station.x, station.y)])
        gauge = Gauge(x=x, y=y, watercourse=station.watercourse, river=station.river)
        placement = place(gauge, segments, map_radius=request.map_radius, reach_up=request.reach_up)
        result: Catchment | None = None
        if placement is None:
            fields["refusal_cause"] = "no_river"
            fields["refusal_message"] = (
                f"no mapped river line within {request.map_radius:g} m of the station"
            )
        else:
            keep = placement.model_dump(exclude={"reach", "position"})
            fields |= {**keep, "uncertainty_m": placement.reach.uncertainty}
            catchment_request = CatchmentRequest(
                seed=(x, y),
                seed_crs=segments_crs,
                reach=placement.reach,
                outline_tolerance=request.outline_tolerance,
            )
            try:
                result = await asyncio.to_thread(delineate, catchment_request, repository)
            except MixedGridRefusal as exc:
                fields |= {"refusal_cause": "mixed_grid", "refusal_message": str(exc)}
                fields["grid_tiles"] = tuple(exc.tiles)
            except CatchmentError as exc:
                fields |= {"refusal_cause": "other", "refusal_message": str(exc)}
        found = None
        if result is not None:
            assert result.gauge is not None
            m = result.meta
            if ref is not None:
                found = agreement(result.fine, ref, m)
                keys = ("nve_in_ours", "ours_in_nve", "area_ratio", "divide_offset_m")
                fields |= {k: getattr(found, k) for k in keys}
            # Our fine area on the lattice, as the agreement counts it: the
            # nodes inside the outline, each its cell's area.
            ours = found.ours if found is not None else nodes_inside([result.fine], m)
            fields |= _gauge_fields(result.gauge)
            fields |= {
                "nodes": result.nodes,
                "fine_area_km2": ours * m.delta_x * m.delta_y / 1e6,
                "reduced_area_km2": result.reduced.area / 1e6,
                "tiles": sum(1 for b in boxes if b.intersects(result.fine)),
                "windows": len(result.windows),
            }
            sink.catchment(station, result)
        cls, by = classify(found, None if result is None else result.gauge, result is None)
        row = StationResult(
            **fields, station_class=cls, match_by=by, seconds=time.perf_counter() - t0
        )
        rows.append(row)
        sink.row(row)
    return summarise(rows)


__all__ = ["BatchRequest", "BatchSink", "StationResult", "run_batch"]
