"""`catchment_batch.run_batch`: a catchment per station, one at a time (increment 29, PR 4).

`docs/increments/29-nve-reference-catchments.md`, "The batch", "Agreement
and classes" and "The red suites", PR 4's `test_catchment_batch.py`. The
terrain, rivers, stations and references are `batch_fixtures`' (its
docstring says why the terrain is the stage B terrain of `test_catchment.py`
rather than the valley of `test_cli_catchment.py`); the repository is in
memory and the sink a list. No file, no network.

Interface pinned here (the design names `BatchRequest` and its four fields,
`BatchSink` with `catchment(station, result)` and `row(row)`, `run_batch`'s
eight parameters in its order, `StationResult` and the field list; the
names it leaves open are chosen here and listed in the handback):

- `run_batch(request, repository, stations, stations_crs, segments,
  segments_crs, references, sink)` is a coroutine returning
  `reference.Summary`. `references` (or None) are in the river file's CRS.
- `StationResult` fields read here: `station`, `name`, `station_class` (the
  design's "class"), `match_by`, `refusal_message`, `refusal_cause`,
  `grid_tiles` (one tile from each grid, for `mixed_grid`; else None); the
  placement's `placed_on`, `distance_m`, `uncertainty_m`, `objectid`,
  `elvid`, `reach_down_m`; `causes` (the joined list of `GaugeResult`);
  `nodes`, `fine_area_km2`, `reference_area_km2` (NVE's polygon),
  `nve_area_km2` (the station layer's), `nve_in_ours`, `ours_in_nve`,
  `area_ratio`; `tiles` (how many of the repository's tiles meet our fine
  outline), `windows` (how many) and `seconds`. A refused row has None for
  every placement, gauge and agreement field it could not have.
- `sink.catchment` is called once per station that has a catchment, before
  its row, and never for a refused station.
- A station number in `only` that is not in the list is refused with a
  `ValueError` naming it, before any station runs; `only` keeps file order.

Expected values were measured on PR 2's code by `place` and `delineate`
(the five classes' causes, node counts and the mixed-grid refusal), and the
oracles are one flood of the raw terrain (`batch_fixtures.flood_mask`).
Numeric bounds: node counts exact; areas at 1e-9 relative (the largest here
is about 0.55 km2); the 16 m and 30 m placement figures at 1e-6 m.
"""

from __future__ import annotations

import asyncio
import threading
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from pydantic import ValidationError
from shapely.geometry import Polygon, box

import batch_fixtures as bf
import gauge_fixtures as gf
import tin_engine.catchment as catchment
import tin_engine.catchment_batch as cb
from catchment_fixtures import MemoryRepository, filled
from tin_engine.io.models import DemTile
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_stations


@dataclass
class ListSink:
    """A `BatchSink` that keeps what it is given, and the order of the calls."""

    catchments: list[tuple[Any, Any]] = field(default_factory=list)
    rows: list[Any] = field(default_factory=list)
    calls: list[tuple[str, str]] = field(default_factory=list)

    def catchment(self, station: Any, result: Any) -> None:
        self.catchments.append((station, result))
        self.calls.append(("catchment", station.station))

    def row(self, row: Any) -> None:
        self.rows.append(row)
        self.calls.append(("row", row.station))


@dataclass
class Inputs:
    stations: tuple[Any, ...]
    stations_crs: str
    segments: tuple[Any, ...]
    segments_crs: str


def inputs(tmp: Path, specs: Sequence[bf.Spec] = bf.FIVE, mixed: bool = False) -> Inputs:
    stations, scrs = read_stations(bf.write_stations(tmp / "stations.geojson", specs))
    segments, rcrs, _ = read_segments(bf.write_rivers(tmp / "rivers.geojson", mixed=mixed))
    return Inputs(stations, scrs, segments, rcrs)


def run(
    tmp: Path,
    *,
    specs: Sequence[bf.Spec] = bf.FIVE,
    tiles: dict[str, DemTile] | None = None,
    references: Any = "default",
    mixed: bool = False,
    repository: Any = None,
    **request: Any,
) -> tuple[ListSink, Any]:
    i = inputs(tmp, specs, mixed)
    repo = repository or MemoryRepository(bf.basin_tiles() if tiles is None else tiles)
    refs = bf.references() if references == "default" else references
    sink = ListSink()
    summary = asyncio.run(
        cb.run_batch(
            cb.BatchRequest(**request), repo, i.stations, i.stations_crs, i.segments,
            i.segments_crs, refs, sink,
        )
    )  # fmt: skip
    return sink, summary


@pytest.fixture(scope="module")
def five(tmp_path_factory: pytest.TempPathFactory) -> tuple[ListSink, Any]:
    """The five stations of `batch_fixtures`, with their references, once."""
    return run(tmp_path_factory.mktemp("five"))


def by_station(sink: ListSink) -> dict[str, Any]:
    return {r.station: r for r in sink.rows}


# ---------------------------------------------------------------------------
# The request
# ---------------------------------------------------------------------------


def test_the_request_defaults_and_is_frozen() -> None:
    r = cb.BatchRequest()
    assert (r.map_radius, r.reach_up, r.outline_tolerance, r.only) == (500.0, 1000.0, None, ())
    with pytest.raises(ValidationError):
        r.map_radius = 100.0  # type: ignore[misc]


def test_run_batch_is_a_coroutine_function() -> None:
    assert asyncio.iscoroutinefunction(cb.run_batch)


# ---------------------------------------------------------------------------
# Five stations: the rows, the classes, the order, the summary
# ---------------------------------------------------------------------------


def test_one_row_per_station_in_file_order(five: tuple[ListSink, Any]) -> None:
    sink, _ = five
    assert [r.station for r in sink.rows] == [s.station for s in bf.FIVE]
    assert [r.name for r in sink.rows] == [s.name for s in bf.FIVE]


def test_the_classes(five: tuple[ListSink, Any]) -> None:
    rows = by_station(five[0])
    assert {n: r.station_class for n, r in rows.items()} == {
        bf.TREFF.station: "match",
        bf.BOM.station: "miss",
        bf.SAMLOP.station: "uncertain",
        bf.SLUTT.station: "uncertain",
        bf.LANGT.station: "refused",
    }
    assert rows[bf.TREFF.station].match_by == "overlap"
    assert all(r.match_by is None for n, r in rows.items() if n != bf.TREFF.station)


def test_the_match_is_its_own_flood_node_for_node(five: tuple[ListSink, Any]) -> None:
    r = by_station(five[0])[bf.TREFF.station]
    mask = filled(bf.flood_mask(bf.TREFF))
    assert r.nodes == int(bf.flood_mask(bf.TREFF).sum())
    assert r.nve_in_ours == 1.0 and r.ours_in_nve == 1.0
    assert r.area_ratio == pytest.approx(1.0, rel=1e-9)
    assert r.fine_area_km2 == pytest.approx(mask.sum() * gf.CELL_KM2, rel=1e-9)
    assert r.reference_area_km2 == pytest.approx(bf.references()[r.station].area / 1e6, rel=1e-9)
    assert r.nve_area_km2 == 12.5  # the station layer's, from the stations file


def test_the_miss_shares_nothing_with_its_moved_reference(five: tuple[ListSink, Any]) -> None:
    r = by_station(five[0])[bf.BOM.station]
    assert r.nve_in_ours == 0.0 and r.ours_in_nve == 0.0
    assert tuple(r.causes) == ()


def test_a_confluence_just_below_is_uncertain_by_its_swing(five: tuple[ListSink, Any]) -> None:
    r = by_station(five[0])[bf.SAMLOP.station]
    assert "swing" in r.causes
    assert r.nve_in_ours == 1.0 and r.ours_in_nve == 1.0  # uncertain beats match


def test_a_river_ending_within_u_is_uncertain_although_its_overlap_is_100(
    five: tuple[ListSink, Any],
) -> None:
    r = by_station(five[0])[bf.SLUTT.station]
    assert tuple(r.causes) == ("downstream_unread",)
    assert r.reach_down_m == pytest.approx(10.0, abs=1e-6)
    assert r.nve_in_ours == 1.0 and r.ours_in_nve == 1.0
    assert r.objectid == 201 and r.elvid == "2-2-1"


def test_the_tiered_placement_is_in_every_placed_row(five: tuple[ListSink, Any]) -> None:
    """Each placed station has its own watercourse number on its line, and
    the A line's exact copy (objectid 102) was dropped by the reader."""
    rows = by_station(five[0])
    for spec in (bf.TREFF, bf.BOM, bf.SAMLOP, bf.SLUTT):
        r = rows[spec.station]
        assert r.placed_on == "number"
        assert r.distance_m == pytest.approx(bf.EAST_M, abs=1e-6)
        assert r.uncertainty_m == pytest.approx(30.0, abs=1e-6)
    assert {rows[s.station].objectid for s in (bf.TREFF, bf.BOM, bf.SAMLOP)} == {101}


def test_no_river_line_near_is_refused_and_the_batch_goes_on(five: tuple[ListSink, Any]) -> None:
    sink, _ = five
    r = by_station(sink)[bf.LANGT.station]
    assert r.refusal_cause == "no_river"
    assert "no mapped river line within 500 m of the station" in r.refusal_message
    assert r.grid_tiles is None
    assert r.placed_on is None and r.nodes is None and r.nve_in_ours is None
    assert sink.rows[-1].station == bf.LANGT.station  # the last row, after the others


def test_a_catchment_per_placed_station_before_its_row(five: tuple[ListSink, Any]) -> None:
    sink, _ = five
    placed = [s.station for s in bf.FIVE if s is not bf.LANGT]
    assert [st.station for st, _ in sink.catchments] == placed
    for station, result in sink.catchments:
        assert result.gauge is not None
        assert isinstance(result.fine, Polygon)
        assert sink.calls.index(("catchment", station.station)) + 1 == sink.calls.index(
            ("row", station.station)
        )
    assert ("catchment", bf.LANGT.station) not in sink.calls


def test_tiles_windows_and_seconds(five: tuple[ListSink, Any]) -> None:
    """`tiles` counts the repository's tiles whose node box meets our fine
    outline (the four quadrants of the terrain, cut at row 150 and column
    350, overlapping by one node line); `windows` counts the floods'
    windows the result lists."""
    sink, _ = five
    boxes = {}
    for name, tile in bf.basin_tiles().items():
        m = tile.meta
        boxes[name] = box(m.x_min, m.y_max - (m.rows - 1) * m.delta_y,
                          m.x_min + (m.cols - 1) * m.delta_x, m.y_max)  # fmt: skip
    results = {st.station: res for st, res in sink.catchments}
    for r in sink.rows:
        if r.station_class == "refused":
            continue
        fine = results[r.station].fine
        assert r.tiles == sum(1 for b in boxes.values() if b.intersects(fine))
        assert r.windows == len(results[r.station].windows) >= 1
        assert r.seconds >= 0.0


def test_the_summary(five: tuple[ListSink, Any]) -> None:
    _, summary = five
    s = summary.model_dump(mode="json")
    assert s["stations"] == 5
    assert s["classes"] == {"match": 1, "close": 0, "miss": 1, "uncertain": 2, "refused": 1}
    assert s["match_by"] == {"overlap": 1, "offset": 0}
    assert s["refusal_causes"]["no_river"] == 1
    assert s["uncertain_causes"]["swing"] >= 1
    assert s["uncertain_causes"]["downstream_unread"] >= 1
    assert s["known_refusals"]["count"] == 0


def test_the_batch_is_deterministic(five: tuple[ListSink, Any], tmp_path: Path) -> None:
    sink, summary = run(tmp_path)
    keep = ("station", "station_class", "nodes", "nve_in_ours", "causes", "refusal_message")
    assert [[getattr(r, k) for k in keep] for r in sink.rows] == [
        [getattr(r, k) for k in keep] for r in five[0].rows
    ]
    assert summary.model_dump_json() == five[1].model_dump_json()


# ---------------------------------------------------------------------------
# Refusals and bugs
# ---------------------------------------------------------------------------


class FailingRepository(MemoryRepository):
    """Raises `exc` on the first tile load: a bug, not a data refusal."""

    def __init__(self, tiles: dict[str, DemTile], exc: BaseException) -> None:
        super().__init__(tiles)
        self.exc = exc

    def load(self, name: str) -> DemTile:
        raise self.exc


@pytest.mark.parametrize(
    "exc",
    [RuntimeError("a bug in the tile store"), ValueError("an ordinary bug, not a refusal")],
    ids=["runtime-error", "plain-value-error"],
)
def test_a_bug_type_exception_stops_the_batch(tmp_path: Path, exc: BaseException) -> None:
    """The refused first station (no river; nothing loaded) gets its row; the
    second station's load raises, and that exception, itself, ends the batch:
    no row for it or for the third. A plain `ValueError` is not a refusal
    either: only `CatchmentError` is."""
    repo = FailingRepository(bf.basin_tiles(), exc)
    with pytest.raises(type(exc)) as raised:
        run(tmp_path, specs=(bf.LANGT, bf.TREFF, bf.BOM), repository=repo)
    assert raised.value is exc


def test_a_bug_leaves_only_the_rows_before_it(tmp_path: Path) -> None:
    repo = FailingRepository(bf.basin_tiles(), RuntimeError("a bug"))
    i = inputs(tmp_path, (bf.LANGT, bf.TREFF, bf.BOM))
    sink = ListSink()
    with pytest.raises(RuntimeError, match="a bug"):
        asyncio.run(
            cb.run_batch(cb.BatchRequest(), repo, i.stations, i.stations_crs, i.segments,
                         i.segments_crs, bf.references(), sink)
        )  # fmt: skip
    assert [r.station for r in sink.rows] == [bf.LANGT.station]
    assert sink.catchments == []


def test_nodata_in_the_catchment_is_refused_as_other(tmp_path: Path) -> None:
    """A NaN node inside the placed catchment, off the chain: `delineate`
    refuses (the catchment reaches NoData), and the row's cause is `other`,
    its message the refusal's words."""
    z = bf.basins()
    z[100, bf.COL + 2] = np.nan
    assert bf.flood_mask(bf.TREFF)[100, bf.COL + 2] == 1  # the premise: inside
    tiles = quadrants_of(z)
    sink, summary = run(tmp_path, specs=(bf.TREFF,), tiles=tiles)
    (r,) = sink.rows
    assert r.station_class == "refused"
    assert r.refusal_cause == "other"
    assert "NoData" in r.refusal_message
    assert sink.catchments == []
    assert summary.model_dump(mode="json")["refusal_causes"]["other"] == 1


def quadrants_of(z: np.ndarray) -> dict[str, DemTile]:
    from mosaic_fixtures import quadrants

    return quadrants(gf.tile_of(z), row_cut=gf.BIG_ROWS // 2, col_cut=gf.BIG_COLS // 2, overlap=1)


def test_a_window_on_two_grids_is_refused_mixed_grid_and_the_batch_goes_on(tmp_path: Path) -> None:
    """The sixth station: its window selects `a.tif` and the half-cell-shifted
    `b.tif`, and neither covers it. Refused with `refusal_cause =
    "mixed_grid"` and one tile from each grid; the next station, inside
    `a.tif`, is a match."""
    refs = {bf.GRENSE.station: box(*gf.lat(590, 160), *gf.lat(610, 140)),
            bf.TREFF.station: bf.reference_of(bf.TREFF)}  # fmt: skip
    sink, summary = run(
        tmp_path, specs=(bf.GRENSE, bf.TREFF), tiles=bf.mixed_tiles(), references=refs, mixed=True
    )
    grense, treff = sink.rows
    assert grense.station == bf.GRENSE.station
    assert grense.station_class == "refused"
    assert grense.refusal_cause == "mixed_grid"
    assert set(grense.grid_tiles) == {"a.tif", "b.tif"}
    assert "two lattices" in grense.refusal_message
    assert treff.station_class == "match"
    known = summary.model_dump(mode="json")["known_refusals"]
    assert known["count"] == 1 and known["stations"] == [bf.GRENSE.station]
    assert summary.model_dump(mode="json")["refusal_causes"]["mixed_grid"] == 1


# ---------------------------------------------------------------------------
# One station at a time, off the event loop; `only`; no reference
# ---------------------------------------------------------------------------


class RecordingRepository(MemoryRepository):
    """Records the thread of every tile load."""

    def __init__(self, tiles: dict[str, DemTile]) -> None:
        super().__init__(tiles)
        self.threads: list[int] = []

    def load(self, name: str) -> DemTile:
        self.threads.append(threading.get_ident())
        return super().load(name)


async def test_each_delineate_runs_off_the_event_loop(tmp_path: Path) -> None:
    """`delineate` in `asyncio.to_thread`: every tile load happens on a
    thread that is not the one running the event loop."""
    loop_thread = threading.get_ident()
    repo = RecordingRepository(bf.basin_tiles())
    i = inputs(tmp_path, (bf.TREFF, bf.BOM))
    await cb.run_batch(cb.BatchRequest(), repo, i.stations, i.stations_crs, i.segments,
                       i.segments_crs, None, ListSink())  # fmt: skip
    assert repo.threads
    assert loop_thread not in repo.threads


def test_only_runs_the_named_stations_in_file_order(tmp_path: Path) -> None:
    sink, summary = run(tmp_path, only=(bf.SLUTT.station, bf.TREFF.station))
    assert [r.station for r in sink.rows] == [bf.TREFF.station, bf.SLUTT.station]
    assert summary.model_dump(mode="json")["stations"] == 2


def test_an_unknown_station_in_only_is_refused_by_name(tmp_path: Path) -> None:
    i = inputs(tmp_path)
    sink = ListSink()
    with pytest.raises(ValueError, match=r"9\.9\.9"):
        asyncio.run(
            cb.run_batch(cb.BatchRequest(only=("9.9.9",)), MemoryRepository(bf.basin_tiles()),
                         i.stations, i.stations_crs, i.segments, i.segments_crs, None, sink)
        )  # fmt: skip
    assert sink.rows == []


def test_without_references_no_class_beyond_refused_and_uncertain(tmp_path: Path) -> None:
    only = (bf.TREFF.station, bf.SAMLOP.station, bf.LANGT.station)
    sink, summary = run(tmp_path, references=None, only=only)
    rows = by_station(sink)
    assert rows[bf.TREFF.station].station_class is None
    assert rows[bf.TREFF.station].match_by is None
    assert rows[bf.TREFF.station].nve_in_ours is None
    assert rows[bf.TREFF.station].reference_area_km2 is None
    assert rows[bf.TREFF.station].nodes == int(bf.flood_mask(bf.TREFF).sum())
    assert rows[bf.SAMLOP.station].station_class == "uncertain"
    assert rows[bf.LANGT.station].station_class == "refused"
    classes = summary.model_dump(mode="json")["classes"]
    assert classes["match"] == classes["close"] == classes["miss"] == 0
    assert len(sink.catchments) == 2


def test_the_map_radius_and_reach_up_are_passed_on(tmp_path: Path) -> None:
    """At a 10 m map radius no line is within reach of a station 16 m from
    its line; at reach_up 200 m the reach above `P` is 200 m."""
    sink, _ = run(tmp_path, only=(bf.TREFF.station,), map_radius=10.0)
    assert sink.rows[0].refusal_cause == "no_river"
    assert "within 10 m" in sink.rows[0].refusal_message
    sink, _ = run(tmp_path, only=(bf.TREFF.station,), reach_up=200.0)
    assert sink.rows[0].reach_up_m == pytest.approx(200.0, abs=1e-6)


def test_the_outline_tolerance_is_passed_on(tmp_path: Path) -> None:
    sink, _ = run(tmp_path, only=(bf.TREFF.station,), outline_tolerance=0.0)
    ((_, result),) = sink.catchments
    assert result.tolerance == 0.0
    assert result.reduced.area == pytest.approx(result.fine.area, rel=1e-12)


# ---------------------------------------------------------------------------
# PR 4's code review, round 1: change (d)
# ---------------------------------------------------------------------------


def test_a_river_crs_not_the_dems_is_refused_before_the_first_station(tmp_path: Path) -> None:
    """Change (d), for callers other than the command: `segments_crs`
    EPSG:32633 over the DEM's EPSG:25833 raises `ValueError` naming both,
    before any station runs, so the sink receives nothing (before change (d),
    each station was a `refused` row with cause `other`)."""
    i = inputs(tmp_path)
    sink = ListSink()
    with pytest.raises(ValueError) as caught:
        asyncio.run(
            cb.run_batch(cb.BatchRequest(), MemoryRepository(bf.basin_tiles()), i.stations,
                         i.stations_crs, i.segments, "EPSG:32633", None, sink)
        )  # fmt: skip
    assert "32633" in str(caught.value) and "25833" in str(caught.value), str(caught.value)
    assert sink.calls == []


def test_check_reach_crs_is_the_one_rule() -> None:
    """Change (d)'s one place: `catchment.check_reach_crs(crs, repository)`
    returns None when `crs` is the first tile's CRS and raises `ValueError`
    naming both CRSs, the river file's first, when it is not."""
    repository = MemoryRepository(bf.basin_tiles())
    assert catchment.check_reach_crs(gf.EPSG, repository) is None
    with pytest.raises(ValueError, match=r"(?s)river file's CRS.*32633.*DEM's.*25833"):
        catchment.check_reach_crs("EPSG:32633", repository)


# ---------------------------------------------------------------------------
# PR 4, lake gauges: `run_batch(..., lakes=None)` ("Lake gauges", "The batch")
# ---------------------------------------------------------------------------
#
# Interface pinned here, from the design: `run_batch` gains `lakes:
# Sequence[Lake] | None = None` after `sink` (passed here by keyword);
# `StationResult` gains `seeded_by`, `lake_rule`, `lake_number`,
# `lake_name`, `lake_distance_m`, in that order, right after `reach_fork`;
# `Summary` gains `by_seed` with the groups `river` and `lake`. Before the
# change, `run_batch` takes no `lakes` and the row has none of the five.

GAUGE_COLUMNS = (
    "node_offset_m", "chain_nodes", "lowered_nodes", "lowered_max_m", "direction_ok",
    "end_extended_m", "end_closed", "downstream_checked", "a0", "area_up", "area_down",
    "swing", "largest_step", "largest_step_at_m", "checked_up_m", "checked_down_m",
    "drains", "monotone",
)  # fmt: skip
LAKE_COLUMNS = ("seeded_by", "lake_rule", "lake_number", "lake_name", "lake_distance_m")


def the_lake(polygon: Polygon | None = None, number: int | None = bf.LAKE_NUMBER,
             name: str | None = bf.LAKE_NAME) -> Any:  # fmt: skip
    from tin_engine.io.station_set import Lake

    return Lake(number=number, name=name, polygon=polygon or bf.lake_polygon())


def direct(seed: tuple[float, float], polygon: Polygon) -> Any:
    """22's path run directly: `delineate` with the lake and a seed in it."""
    request = catchment.CatchmentRequest(
        seed=seed, seed_crs=gf.EPSG, lakes=(polygon,), lakes_crs=gf.EPSG
    )
    return catchment.delineate(request, MemoryRepository(bf.basin_tiles()))


def run_lakes(
    tmp: Path,
    lakes: Sequence[Any] | None,
    *,
    specs: Sequence[bf.Spec],
    references: Any = None,
    rivers: list[dict[str, Any]] | None = None,
    **request: Any,
) -> tuple[ListSink, Any]:
    from nve_fixtures import collection, write

    stations, scrs = read_stations(bf.write_stations(tmp / "stations.geojson", specs))
    river_doc = collection(bf.river_features() if rivers is None else rivers, crs=gf.EPSG)
    segments, rcrs, _ = read_segments(write(tmp / "rivers.geojson", river_doc))
    sink = ListSink()
    summary = asyncio.run(
        cb.run_batch(cb.BatchRequest(**request), MemoryRepository(bf.basin_tiles()), stations,
                     scrs, segments, rcrs, references, sink, lakes=lakes)
    )  # fmt: skip
    return sink, summary


@pytest.fixture(scope="module")
def lake_direct() -> Any:
    return direct(bf.xy(bf.SAMLOP), bf.lake_polygon())


@pytest.fixture(scope="module")
def lake_run(tmp_path_factory: pytest.TempPathFactory, lake_direct: Any) -> tuple[ListSink, Any]:
    """`SAMLOP` inside THE LAKE, with a reference drawn from 22's catchment of
    it, and `TREFF` on the river with its own flood as its reference."""
    refs = {bf.SAMLOP.station: lake_direct.fine, bf.TREFF.station: bf.reference_of(bf.TREFF)}
    return run_lakes(tmp_path_factory.mktemp("lake"), [the_lake()],
                     specs=(bf.SAMLOP, bf.TREFF), references=refs)  # fmt: skip


def test_the_five_columns_follow_reach_fork() -> None:
    import dataclasses

    names = [f.name for f in dataclasses.fields(cb.StationResult)]
    at = names.index("reach_fork")
    assert tuple(names[at + 1 : at + 6]) == LAKE_COLUMNS


def test_a_station_inside_a_lake_is_seeded_by_the_lake(lake_run: tuple[ListSink, Any]) -> None:
    r = by_station(lake_run[0])[bf.SAMLOP.station]
    assert r.seeded_by == "lake"
    assert r.lake_rule == "inside"
    assert (r.lake_number, r.lake_name) == (bf.LAKE_NUMBER, bf.LAKE_NAME)
    assert r.lake_distance_m == 0.0


def test_a_lake_row_has_no_gauge_burn_or_sensitivity(lake_run: tuple[ListSink, Any]) -> None:
    r = by_station(lake_run[0])[bf.SAMLOP.station]
    assert {k: getattr(r, k) for k in GAUGE_COLUMNS} == dict.fromkeys(GAUGE_COLUMNS)
    assert tuple(r.causes) == ()


def test_the_lake_catchment_is_22s_path_run_directly(
    lake_run: tuple[ListSink, Any], lake_direct: Any
) -> None:
    """5715 nodes (22's path) against 5540 for the placed node on the river."""
    sink, _ = lake_run
    r = by_station(sink)[bf.SAMLOP.station]
    result = dict((st.station, res) for st, res in sink.catchments)[bf.SAMLOP.station]
    assert result.gauge is None
    assert r.nodes == result.nodes == lake_direct.nodes == 5715
    assert result.fine.equals(lake_direct.fine)
    assert int(bf.flood_mask(bf.SAMLOP).sum()) == 5540  # the river path's, told apart


def test_a_lake_row_with_its_reference_is_a_match_never_uncertain(
    lake_run: tuple[ListSink, Any],
) -> None:
    """`SAMLOP`'s river row is `uncertain` (the confluence 20 m below); its
    lake row decides from the agreement alone."""
    r = by_station(lake_run[0])[bf.SAMLOP.station]
    assert r.station_class == "match" and r.match_by == "overlap"
    assert r.nve_in_ours == 1.0 and r.ours_in_nve == 1.0


def test_a_river_row_beside_it_keeps_the_river_path(lake_run: tuple[ListSink, Any]) -> None:
    r = by_station(lake_run[0])[bf.TREFF.station]
    assert r.seeded_by == "river"
    assert (r.lake_rule, r.lake_number, r.lake_name, r.lake_distance_m) == (None,) * 4
    assert r.station_class == "match"
    assert r.node_offset_m is not None and r.swing is not None


def test_by_seed_counts_one_of_each(lake_run: tuple[ListSink, Any]) -> None:
    s = lake_run[1].model_dump(mode="json")
    assert list(s["by_seed"]) == ["river", "lake"]
    assert s["by_seed"]["river"]["stations"] == 1
    assert s["by_seed"]["lake"]["stations"] == 1
    assert s["by_seed"]["lake"]["classes"]["match"] == 1
    assert s["by_seed"]["lake"]["uncertain_share"] == 0.0
    assert set(s["by_seed"]["lake"]) == set(s["by_size"]["10-100"])  # a band's contents


def test_without_lakes_the_same_station_keeps_its_river_row(five: tuple[ListSink, Any]) -> None:
    """`lakes` None (the default, as every earlier test calls it): each
    placed row is seeded by the river, the refused-by-`place` row by nothing,
    and no row has a lake."""
    rows = by_station(five[0])
    assert rows[bf.SAMLOP.station].station_class == "uncertain"
    for spec in bf.FIVE:
        r = rows[spec.station]
        assert r.seeded_by == (None if spec is bf.LANGT else "river"), spec.station
        assert (r.lake_rule, r.lake_number, r.lake_name, r.lake_distance_m) == (None,) * 4
    by_seed = five[1].model_dump(mode="json")["by_seed"]
    assert by_seed["lake"]["stations"] == 0 and by_seed["river"]["stations"] == 4


def test_a_station_inside_two_overlapping_lakes_is_refused_and_the_batch_goes_on(
    tmp_path: Path,
) -> None:
    """22's `_lake` refuses a seed point in two lakes (a `LakeError`, so a
    `CatchmentError`): a `refused` row with cause `other`."""
    two = [the_lake(), the_lake(bf.overlapping_polygon(), number=4111, name="Indre")]
    refs = {bf.TREFF.station: bf.reference_of(bf.TREFF)}
    sink, summary = run_lakes(tmp_path, two, specs=(bf.SAMLOP, bf.TREFF), references=refs)
    samlop, treff = sink.rows
    assert samlop.station_class == "refused"
    assert samlop.refusal_cause == "other"
    assert "2 lakes" in samlop.refusal_message
    assert samlop.seeded_by == "lake"
    assert treff.station_class == "match"
    assert [st.station for st, _ in sink.catchments] == [bf.TREFF.station]
    assert summary.model_dump(mode="json")["refusal_causes"]["other"] == 1


def test_a_station_with_no_placement_inside_a_lake_is_seeded_by_it(
    tmp_path: Path, lake_direct: Any
) -> None:
    """At a 10 m map radius `SAMLOP` (16 m from its line) has no placement:
    it is `no_river` only when `lake_seed` is None, and here it is not."""
    refs = {bf.SAMLOP.station: lake_direct.fine}
    sink, _ = run_lakes(tmp_path, [the_lake()], specs=(bf.SAMLOP,), references=refs,
                        map_radius=10.0)  # fmt: skip
    (r,) = sink.rows
    assert r.refusal_cause is None and r.refusal_message is None
    assert r.placed_on is None
    assert (r.seeded_by, r.lake_rule) == ("lake", "inside")
    assert r.nodes == lake_direct.nodes
    assert r.station_class == "match"


def test_without_a_reference_a_lake_row_has_no_class(tmp_path: Path) -> None:
    """`classify(None, None)`: neither scored nor `uncertain`."""
    sink, _ = run_lakes(tmp_path, [the_lake()], specs=(bf.SAMLOP,))
    (r,) = sink.rows
    assert r.seeded_by == "lake" and r.station_class is None


def test_a_lake_line_station_is_seeded_from_p(tmp_path: Path) -> None:
    """`TREFF` placed on `A` mapped as a lake centreline; `P` (column 100) is
    in the lake, the station 6 m east of its shore: `lake_line`, and the
    catchment is 22's path seeded at `P`."""
    polygon = bf.lake_line_polygon()
    sink, _ = run_lakes(tmp_path, [the_lake(polygon)], specs=(bf.TREFF,),
                        rivers=bf.lake_line_rivers())  # fmt: skip
    (r,) = sink.rows
    assert r.lake is True
    assert (r.seeded_by, r.lake_rule) == ("lake", "lake_line")
    assert r.lake_distance_m == pytest.approx(6.0, abs=1e-6)  # at x 5e5 m
    expected = direct(gf.lat(bf.COL, bf.TREFF.row), polygon)
    assert r.nodes == expected.nodes
    ((_, result),) = sink.catchments
    assert result.fine.equals(expected.fine)


def test_an_empty_lake_list_is_the_river_path(tmp_path: Path) -> None:
    sink, _ = run_lakes(tmp_path, [], specs=(bf.SAMLOP,))
    (r,) = sink.rows
    assert r.seeded_by == "river" and r.station_class == "uncertain"
