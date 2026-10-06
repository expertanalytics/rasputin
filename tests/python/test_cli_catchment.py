"""`rasputin catchment`, and its output meshed by `rasputin mesh --domain` (increment 22, PR 1).

`docs/increments/22-auto-catchment.md`, "The command", "The seed" and "The
red suites" (PR 1, `test_cli_catchment.py`), and the Bygdin numbers under
"What the data says about Bygdin".

Pinned here, from the design:

- `rasputin catchment --dem PATH --seed X Y [--seed-crs CRS] [--lakes PATH
  [--lakes-layer NAME]] --out FILE.geojson`. `--seed-crs` defaults to
  EPSG:4326 (`LON LAT`). `--lakes` is a `.gpkg` (with `--lakes-layer` when it
  holds more than one features table) or `.geojson`/`.json`, in the file's
  own CRS.
- The output is a `FeatureCollection` of one `Feature`, the polygon in the
  DEM's CRS, with a `crs` member naming it; `domain.read_domain` reads it.
  In PR 1 the polygon is the fine outline, so its area is the marching-squares
  area of the catchment's nodes, holes filled; coordinates survive the round
  trip (`repr` precision), so the area read back is that, to 1e-12.
- stderr has a line per window, and the fine outline's vertex count and area,
  worded since increment 25 as `docs/increments/25-plain-output.md` ("stderr,
  reworded") gives them: `window <k>: ..., searched in <t> s, catchment
  inside the window` (or `window widened to the <sides>`), `start: lake of
  <a> km2, <n> DEM nodes` or `start: the outlet node at (x, y) (an outlet
  must lie on the flow line; it is not moved there)`, `catchment: <n> DEM
  nodes, <a> km2`, and `outline along DEM cells: <v> vertices, <a> km2, <p>
  separate patch(es) left out (<n> nodes), <g> enclosed gap(s) filled
  (...)`. The parts the design elides (`...`) are not pinned.
- Every refusal is a non-zero exit and no file.

PR 2 (the reduction): the default output is the reduced polygon, whose area is
the fine outline's to 1e-9 relative (the design's area guarantee), so the
PR 1 area checks above compare at 1e-9 (amended in the PR 2 red step; at
1e-12 they held only for the fine ring). `--outline-tolerance 0` writes the
fine ring with only its exactly collinear vertices removed. The properties
gain `reduced_vertices`, `reduced_area_m2` and `outline_tolerance_m` (names
chosen here, beside PR 1's `fine_vertices` and `fine_area_m2`). stderr gains
a line with the word `reduced`, its vertices, its area in km and the
difference in m2.

The Bygdin test runs against Ola's data in ../rasputin_data and is skipped
when it is absent: the area must be within 2 % of NVE's 305.54 km^2 (delfelt
1187), as the design's acceptance asks.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from shapely.geometry import Point, Polygon

import gauge_fixtures as gf
from catchment_fixtures import (
    EPSG,
    X0,
    Y0,
    D,
    bowl,
    filled,
    full_flood,
    lake_box,
    lat,
    seeds_in,
    tile_of,
)
from cli_helpers import invoke, ran
from gpkg_fixtures import DTM10, OLA_NORWAY, Layer, Row, write_gpkg
from mosaic_fixtures import quadrants
from nve_fixtures import collection, river, write
from test_catchment import moved, to_4326
from test_cli_mesh_mosaic import write_tiles
from test_outline import square_area
from tin_engine._core import accumulate as core_accumulate
from tin_engine._core import upstream as core_upstream
from tin_engine.crs import parse_crs, reprojector
from tin_engine.domain import read_domain
from tin_engine.raster import to_core

SEED = lat(100, 200)  # the bowl's centre node, inside the lake
NVE_BYGDIN_KM2 = 305.54  # NVE delfelt 1187, fetched 2026-09-29 (see the design)
BYGDIN_SEED = ("8.5425", "61.3512")


def utm_seed(x: float = SEED[0], y: float = SEED[1]) -> tuple[str, ...]:
    return ("--seed", repr(x), repr(y), "--seed-crs", EPSG)


def lake_geojson(path: Path, lake: Any = None, crs: str = EPSG) -> Path:
    lake = lake_box() if lake is None else lake
    doc = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": crs}},
        "features": [
            {
                "type": "Feature",
                "properties": {"name": "the lake"},
                "geometry": {"type": "Polygon", "coordinates": [list(lake.exterior.coords)]},
            }
        ],
    }
    path.write_text(json.dumps(doc))
    return path


@pytest.fixture
def dem_dir(tmp_path: Path) -> Path:
    write_tiles(tmp_path / "tiles", quadrants(tile_of(bowl()), row_cut=200, col_cut=100, overlap=1))
    return tmp_path / "tiles"


@pytest.fixture
def lakes(tmp_path: Path) -> Path:
    return lake_geojson(tmp_path / "lakes.geojson")


def run(dem: Path, out: Path, *args: str) -> str:
    output = ran("catchment", "--dem", str(dem), *args, "--out", str(out))
    assert out.is_file()
    return output


def refused(tmp_path: Path, *args: str, out: str = "c.geojson") -> str:
    target = tmp_path / out
    code, output = invoke("catchment", *args, "--out", str(target))
    assert code != 0, output
    # Refused by the command, not by Typer for want of one.
    assert "No such command" not in output and "No such option" not in output, output
    assert not target.exists()
    return output


def expected_area() -> float:
    tile = tile_of(bowl())
    return square_area(filled(full_flood(tile, seeds_in(tile, lake_box())))) * D * D


# ---------------------------------------------------------------------------
# The round trip
# ---------------------------------------------------------------------------


def test_the_output_is_one_feature_in_the_dem_crs(
    tmp_path: Path, dem_dir: Path, lakes: Path
) -> None:
    out = tmp_path / "c.geojson"
    run(dem_dir, out, *utm_seed(), "--lakes", str(lakes))
    doc = json.loads(out.read_text())
    assert doc["type"] == "FeatureCollection"
    (feature,) = doc["features"]
    assert feature["type"] == "Feature"
    assert feature["geometry"]["type"] == "Polygon"
    assert parse_crs(doc["crs"]["properties"]["name"]) == parse_crs(EPSG)


def test_read_domain_reads_it_back_with_the_fine_area(
    tmp_path: Path, dem_dir: Path, lakes: Path
) -> None:
    out = tmp_path / "c.geojson"
    run(dem_dir, out, *utm_seed(), "--lakes", str(lakes))
    domain = read_domain(out)
    assert parse_crs(domain.crs) == parse_crs(EPSG)
    assert domain.polygon.is_valid
    assert len(domain.polygon.interiors) == 0
    assert domain.polygon.area == pytest.approx(expected_area(), rel=1e-9)


def test_mesh_accepts_it_as_a_domain(tmp_path: Path, dem_dir: Path, lakes: Path) -> None:
    out = tmp_path / "c.geojson"
    run(dem_dir, out, *utm_seed(), "--lakes", str(lakes))
    vtk = tmp_path / "m.vtk"
    code, output = invoke(
        "mesh", "--dem", str(dem_dir), "--domain", str(out), "--tolerance", "5", "--out", str(vtk)
    )
    assert code == 0, output
    assert vtk.is_file() and vtk.stat().st_size > 0


WINDOW_LINE = re.compile(
    r"\bwindow \d+: .*?searched in [0-9.]+ s, "
    r"(catchment inside the window|window widened to the [a-z]+(, [a-z]+)*)"
)
OUTLINE_LINE = re.compile(
    r"\boutline along DEM cells: \d+ vertices, [0-9.]+ km2, \d+ separate patch(es)? left out "
    r"\(\d+ nodes\), \d+ enclosed gaps? filled \("
)
#: Today's words the rewording replaces (increment 25, "stderr, reworded").
OLD_CATCHMENT_WORDS = (
    "flood ",
    "contained",
    "grown on",
    "seed:",
    "pour node",
    "snapped",
    "of node area",
    "fine outline",
    "rings dropped",
    "holes filled",
)


def test_stderr_reports_the_windows_and_the_fine_outline(
    tmp_path: Path, dem_dir: Path, lakes: Path
) -> None:
    output = run(dem_dir, tmp_path / "c.geojson", *utm_seed(), "--lakes", str(lakes))
    assert WINDOW_LINE.search(output), output
    assert re.search(r"\bstart: lake of [0-9.]+ km2, \d+ DEM nodes\b", output), output
    assert re.search(r"\bcatchment: \d+ DEM nodes, [0-9.]+ km2\b", output), output
    assert OUTLINE_LINE.search(output), output
    assert re.search(r"\breduced outline: \d+ vertices", output), output
    for old in OLD_CATCHMENT_WORDS:
        assert old not in output, old


def test_stderr_names_the_outlet_without_lakes(tmp_path: Path, dem_dir: Path) -> None:
    output = run(dem_dir, tmp_path / "c.geojson", *utm_seed(SEED[0] + 30.0, SEED[1] - 20.0))
    outlet = (
        r"\bstart: the outlet node at \(.+?\) \(an outlet must lie on the flow line; "
        r"it is not moved there\)"
    )
    assert re.search(outlet, output), output
    for old in OLD_CATCHMENT_WORDS:
        assert old not in output, old


# ---------------------------------------------------------------------------
# The seed
# ---------------------------------------------------------------------------


def test_the_seed_crs_defaults_to_wgs84(tmp_path: Path, dem_dir: Path, lakes: Path) -> None:
    utm, wgs = tmp_path / "utm.geojson", tmp_path / "wgs.geojson"
    run(dem_dir, utm, *utm_seed(), "--lakes", str(lakes))
    lon, lat_ = to_4326(*SEED)
    run(dem_dir, wgs, "--seed", repr(lon), repr(lat_), "--lakes", str(lakes))
    assert read_domain(wgs).polygon.equals_exact(read_domain(utm).polygon, 0.0)


def test_a_lake_file_in_another_crs(tmp_path: Path, dem_dir: Path, lakes: Path) -> None:
    here, there = tmp_path / "here.geojson", tmp_path / "there.geojson"
    run(dem_dir, here, *utm_seed(), "--lakes", str(lakes))
    other = lake_geojson(tmp_path / "l3035.geojson", moved(lake_box(), "EPSG:3035"), "EPSG:3035")
    run(dem_dir, there, *utm_seed(), "--lakes", str(other))
    assert read_domain(there).polygon.equals_exact(read_domain(here).polygon, 0.0)


def gpkg_lakes(path: Path) -> Path:
    """Two features tables in EPSG:3035: `lakes` holding the lake, `other` a
    box far from it. With two tables, `--lakes-layer` is needed."""
    lake = Row(1, moved(lake_box(), "EPSG:3035"), {"Code_18": "512"})
    far = Row(1, moved(lake_box(380, 20), "EPSG:3035"), {"Code_18": "512"})
    layers = [Layer("lakes", 3035, [lake], rtree=False), Layer("other", 3035, [far], rtree=False)]
    return write_gpkg(path, layers)


def test_a_geopackage_layer_seeds_like_the_geojson(
    tmp_path: Path, dem_dir: Path, lakes: Path
) -> None:
    here, there = tmp_path / "here.geojson", tmp_path / "there.geojson"
    run(dem_dir, here, *utm_seed(), "--lakes", str(lakes))
    gpkg = gpkg_lakes(tmp_path / "lakes.gpkg")
    run(dem_dir, there, *utm_seed(), "--lakes", str(gpkg), "--lakes-layer", "lakes")
    assert read_domain(there).polygon.equals_exact(read_domain(here).polygon, 0.0)


def test_a_geopackage_of_two_tables_without_a_layer_is_refused(
    tmp_path: Path, dem_dir: Path
) -> None:
    gpkg = gpkg_lakes(tmp_path / "lakes.gpkg")
    refused(tmp_path, "--dem", str(dem_dir), *utm_seed(), "--lakes", str(gpkg))


def test_without_lakes_the_seed_is_a_pour_point(tmp_path: Path, dem_dir: Path) -> None:
    out = tmp_path / "c.geojson"
    run(dem_dir, out, *utm_seed(SEED[0] + 30.0, SEED[1] - 20.0))
    tile = tile_of(bowl())
    seed = np.zeros((400, 200), dtype=np.uint8)
    seed[200, 100] = 1
    area = square_area(filled(full_flood(tile, seed))) * D * D
    assert read_domain(out).polygon.area == pytest.approx(area, rel=1e-9)


# ---------------------------------------------------------------------------
# Refusals: a non-zero exit and no file
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("out", ["c.vtk", "c.wkt", "c"])
def test_an_output_that_is_not_geojson_is_refused(
    tmp_path: Path, dem_dir: Path, lakes: Path, out: str
) -> None:
    refused(tmp_path, "--dem", str(dem_dir), *utm_seed(), "--lakes", str(lakes), out=out)


def test_lakes_layer_without_lakes_is_refused(tmp_path: Path, dem_dir: Path) -> None:
    output = refused(tmp_path, "--dem", str(dem_dir), *utm_seed(), "--lakes-layer", "lakes")
    assert "--lakes" in output


#: Typer's line for a BadParameter with param_hint "--lakes" (6757925,
#: LakeError), as `plain` flattens it. The whole phrase, because the old hint's
#: line, "Invalid value for --dem:", also contains "--lakes" wherever the
#: message quotes the flag.
LAKES_HINT = "Invalid value for --lakes:"
DEM_HINT = "Invalid value for --dem:"


def test_a_seed_in_no_lake_is_refused_under_lakes(
    tmp_path: Path, dem_dir: Path, lakes: Path
) -> None:
    output = refused(
        tmp_path, "--dem", str(dem_dir), *utm_seed(*lat(60, 60)), "--lakes", str(lakes)
    )
    assert LAKES_HINT in output, output
    assert DEM_HINT not in output, output


def test_a_seed_in_two_lakes_is_refused_under_lakes(tmp_path: Path, dem_dir: Path) -> None:
    doc = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": EPSG}},
        "features": [
            {
                "type": "Feature",
                "properties": {"name": name},
                "geometry": {"type": "Polygon", "coordinates": [list(lake.exterior.coords)]},
            }
            for name, lake in (
                ("one", lake_box()),
                ("two", lake_box().buffer(D, join_style="mitre")),
            )
        ],
    }
    path = tmp_path / "two.geojson"
    path.write_text(json.dumps(doc))
    output = refused(tmp_path, "--dem", str(dem_dir), *utm_seed(), "--lakes", str(path))
    assert LAKES_HINT in output, output
    assert DEM_HINT not in output, output


def test_a_catchment_cut_by_the_data_edge_writes_nothing(tmp_path: Path) -> None:
    write_tiles(
        tmp_path / "cut", quadrants(tile_of(bowl(rc=40)), row_cut=200, col_cut=100, overlap=1)
    )
    lake = lake_geojson(tmp_path / "l.geojson", lake_box(40, 100))
    output = refused(
        tmp_path, "--dem", str(tmp_path / "cut"), *utm_seed(*lat(100, 40)), "--lakes", str(lake)
    )
    assert re.search(r"(?i)\bnorth\b", output), output


# ---------------------------------------------------------------------------
# Bygdin, on Ola's data
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not (DTM10.is_dir() and OLA_NORWAY.is_file()), reason="Ola's DTM10 and CORINE are not here"
)
def test_bygdin_is_within_two_percent_of_nve(tmp_path: Path) -> None:
    out = tmp_path / "bygdin.geojson"
    output = run(
        DTM10,
        out,
        "--seed",
        *BYGDIN_SEED,
        "--lakes",
        str(OLA_NORWAY),
        "--lakes-layer",
        "corine2018",
    )
    domain = read_domain(out)
    assert parse_crs(domain.crs) == parse_crs(EPSG)
    km2 = domain.polygon.area / 1e6
    assert abs(km2 / NVE_BYGDIN_KM2 - 1.0) <= 0.02, (km2, output)
    # The seed is mid-lake, so the outline contains it.
    ((x, y),) = reprojector("EPSG:4326", EPSG)([[float(v) for v in BYGDIN_SEED]])
    assert domain.polygon.contains(Point(x, y))


# ---------------------------------------------------------------------------
# PR 2: the reduction
# ---------------------------------------------------------------------------


def feature(out: Path) -> dict[str, Any]:
    (f,) = json.loads(out.read_text())["features"]
    return dict(f)


def vertices(polygon: Any) -> int:
    return len(polygon.exterior.coords) - 1


def fine_corners() -> set[tuple[float, float]]:
    """The fine outline's vertices less the exactly collinear ones: the
    outer ring round the seed of the traced whole-raster catchment. Every
    coordinate is a multiple of 50 m, so the collinearity test is exact."""
    from tin_engine.outline import trace

    tile = tile_of(bowl())
    mask = filled(full_flood(tile, seeds_in(tile, lake_box())))
    rings = [np.asarray(r, dtype=np.float64) for r in trace(mask)]
    world = [np.column_stack([X0 + r[:, 1] * D, Y0 - r[:, 0] * D]) for r in rings]
    (ring,) = [w for w in world if Polygon(w).contains(Point(SEED))]
    if np.array_equal(ring[0], ring[-1]):
        ring = ring[:-1]
    a, b, c = np.roll(ring, 1, axis=0), ring, np.roll(ring, -1, axis=0)
    cross = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    return {(float(x), float(y)) for x, y in ring[cross != 0]}


def test_outline_tolerance_0_writes_the_fine_ring_less_its_collinear_vertices(
    tmp_path: Path, dem_dir: Path, lakes: Path
) -> None:
    out = tmp_path / "fine.geojson"
    run(dem_dir, out, *utm_seed(), "--lakes", str(lakes), "--outline-tolerance", "0")
    polygon = read_domain(out).polygon
    assert {(float(x), float(y)) for x, y in polygon.exterior.coords} == fine_corners()
    assert vertices(polygon) == len(fine_corners())
    assert polygon.area == pytest.approx(expected_area(), rel=1e-12)


def test_the_default_reduces_to_twice_the_cell_keeping_the_area(
    tmp_path: Path, dem_dir: Path, lakes: Path
) -> None:
    fine_out, out = tmp_path / "fine.geojson", tmp_path / "c.geojson"
    run(dem_dir, fine_out, *utm_seed(), "--lakes", str(lakes), "--outline-tolerance", "0")
    run(dem_dir, out, *utm_seed(), "--lakes", str(lakes))
    fine, reduced = read_domain(fine_out).polygon, read_domain(out).polygon
    assert vertices(reduced) < vertices(fine)
    assert abs(reduced.area - fine.area) <= 1e-9 * fine.area
    assert reduced.is_valid
    assert reduced.contains(Point(SEED))
    tolerance = 2 * D
    # Measured, not guaranteed (698b19f); expected here because the fine ring
    # is a traced lattice outline.
    assert shapely.hausdorff_distance(fine.exterior, reduced.exterior, densify=0.05) <= (
        tolerance + D
    )

    props = feature(out)["properties"]
    assert props["outline_tolerance_m"] == pytest.approx(tolerance)
    assert props["reduced_vertices"] == vertices(reduced)
    assert props["fine_vertices"] > props["reduced_vertices"]
    assert props["reduced_area_m2"] == pytest.approx(reduced.area, rel=1e-12)
    assert props["fine_area_m2"] == pytest.approx(fine.area, rel=1e-9)


def test_an_explicit_tolerance_is_kept(tmp_path: Path, dem_dir: Path, lakes: Path) -> None:
    fine_out, out = tmp_path / "fine.geojson", tmp_path / "c.geojson"
    run(dem_dir, fine_out, *utm_seed(), "--lakes", str(lakes), "--outline-tolerance", "0")
    run(dem_dir, out, *utm_seed(), "--lakes", str(lakes), "--outline-tolerance", "50")
    fine, reduced = read_domain(fine_out).polygon, read_domain(out).polygon
    distance = shapely.distance(shapely.points(np.asarray(fine.exterior.coords)), reduced.exterior)
    assert distance.max() <= 50.0 + 1e-6
    assert feature(out)["properties"]["outline_tolerance_m"] == pytest.approx(50.0)


def test_stderr_reports_the_reduced_outline(tmp_path: Path, dem_dir: Path, lakes: Path) -> None:
    output = run(dem_dir, tmp_path / "c.geojson", *utm_seed(), "--lakes", str(lakes))
    lines = [line for line in output.splitlines() if re.search(r"(?i)\breduced\b", line)]
    assert lines, output
    line = " ".join(lines)
    assert re.search(r"\bvertices\b", line), line
    assert re.search(r"km", line), line
    assert re.search(r"\d\s*m(2|\u00b2)(?![a-z])", line.replace("km", "")), line


@pytest.mark.parametrize("value", ["-1", "nan", "inf"])
def test_a_bad_outline_tolerance_is_refused(
    tmp_path: Path, dem_dir: Path, lakes: Path, value: str
) -> None:
    refused(
        tmp_path, "--dem", str(dem_dir), *utm_seed(), "--lakes", str(lakes),
        "--outline-tolerance", value,
    )  # fmt: skip


@pytest.mark.skipif(
    not (DTM10.is_dir() and OLA_NORWAY.is_file()), reason="Ola's DTM10 and CORINE are not here"
)
def test_bygdin_reduced_keeps_the_area_and_meshes_at_10_m(tmp_path: Path) -> None:
    lakes_args = ("--lakes", str(OLA_NORWAY), "--lakes-layer", "corine2018")
    fine_out, out = tmp_path / "fine.geojson", tmp_path / "bygdin.geojson"
    run(DTM10, fine_out, "--seed", *BYGDIN_SEED, *lakes_args, "--outline-tolerance", "0")
    run(DTM10, out, "--seed", *BYGDIN_SEED, *lakes_args)
    fine, reduced = read_domain(fine_out).polygon, read_domain(out).polygon
    assert abs(reduced.area - fine.area) <= 1e-9 * fine.area
    # PR 1 measured 17,812 fine vertices; the design's gauge is about 900.
    assert vertices(reduced) <= 17_812 // 5, vertices(reduced)
    assert feature(out)["properties"]["outline_tolerance_m"] == pytest.approx(20.0)
    vtk = tmp_path / "bygdin.vtk"
    code, output = invoke(
        "mesh", "--dem", str(DTM10), "--domain", str(out), "--tolerance", "10", "--out", str(vtk)
    )
    assert code == 0, output
    assert vtk.is_file() and vtk.stat().st_size > 0


# ---------------------------------------------------------------------------
# Increment 29, PR 2: `--rivers`, `--map-radius` and `--reach-up`
#
# `docs/increments/29-nve-reference-catchments.md`, "The batch" (the paragraph
# on `rasputin catchment`), "Ola's rulings" (the stderr line of copies
# dropped) and "The red suites", PR 2's `test_cli_catchment.py`. The DEM is
# `gauge_fixtures.valley(dam=True)` on a 10 m lattice, cut into four tiles;
# the river file is one line down the valley floor (objectid 8841, river Nea)
# and one exact copy of it, which the reader drops.
#
# The line runs 4 m east of the floor's column, and the station is 16 m
# east of the line, beside row 150, so `U` is 30 m. The valley floor is taken
# across the line ("Following the river", step 2): each resample point's
# cross-section is its own row, columns CC - 2 to CC + 3 within the default
# 30 m corridor (the nearest outside are 34 m and 36 m away, so no node sits
# on the corridor's boundary), and its least elevation is the floor column CC,
# also on the embankment's rows, where the crest ties every column and the
# nearer, CC, wins. So the placed node is (150, CC), 4 m across the river from
# `P`, and the embankment's three nodes are the only ones lowered. Measured on
# PR 1's flood when this was written: the chain drains node by node, no count
# within 30 m carries a flag bit, and the swing is 0.029.
#
# Pinned here (the design gives the stderr sentence as an example; these
# pieces of it are pinned, its numbers' formats are not): the line holds
# "placed on the river line 16 m from the station", "(river Nea, line 8841)",
# "moved 4 m onto the DEM's valley floor", "drain through it", "within 30 m up
# and down the river" and "well defined"; the copies line reads "rivers: 2
# segments read, 1 exact copy dropped" (segments read counts the file's
# features, copies included). The file's properties gain `placed_on`,
# `elvid`, `objectid`, `node_offset_m`, `lowered_nodes`, `reach_up_m`,
# `downstream_checked`, `swing` and `causes` (names from the design's
# `StationResult`). A river file whose CRS is not the DEM's is refused, and so
# are a reach with `--lakes`, no line within `--map-radius`, and a radius or
# reach that is not a finite positive number.
# ---------------------------------------------------------------------------


STATION = gf.lat(gf.CC + 2, 150)
PLACED = (150, gf.CC)
#: The mapped line's column: 4 m east of the valley floor.
LINE_COL = gf.CC + 0.4
NEA = 8841


def gauge_station(x: float = STATION[0], y: float = STATION[1]) -> tuple[str, ...]:
    return ("--seed", repr(x), repr(y), "--seed-crs", gf.EPSG)


@pytest.fixture
def valley_dir(tmp_path: Path) -> Path:
    tile = gf.tile_of(gf.valley(dam=True))
    write_tiles(tmp_path / "valley", quadrants(tile, row_cut=120, col_cut=30, overlap=1))
    return tmp_path / "valley"


def rivers_file(path: Path, line: Any = None, crs: str = gf.EPSG) -> Path:
    line = list(gf.column_line(LINE_COL, 20, 230)) if line is None else line
    features = [
        river(NEA, line, elvid="2-11-1", elvenavn="Nea", vassdragsnr="002.A"),
        river(NEA + 1, line, elvid="2-11-1", elvenavn="Nea", vassdragsnr="002.A"),
    ]
    return write(path, collection(features, crs=crs))


@pytest.fixture
def rivers(tmp_path: Path) -> Path:
    return rivers_file(tmp_path / "rivers.geojson")


def burnt_valley() -> np.ndarray:
    """The valley, burnt along the floor (any chain down column CC that
    crosses the embankment lowers the same three nodes)."""
    return gf.hand_burn(gf.valley(dam=True), [(r, gf.CC) for r in range(52, 166)])


def placed_mask() -> np.ndarray:
    z = burnt_valley()
    seed = np.zeros(z.shape, dtype=np.uint8)
    seed[PLACED] = 1
    return np.asarray(core_upstream(to_core(gf.tile_of(z)), seed).mask)


def expected_swing() -> float:
    count = np.asarray(core_accumulate(to_core(gf.tile_of(burnt_valley()))).count, dtype=float)
    a0 = count[PLACED]
    return float(max(a0 - count[147, gf.CC], count[153, gf.CC] - a0) / a0)


def test_rivers_places_the_gauge_and_says_where(
    tmp_path: Path, valley_dir: Path, rivers: Path
) -> None:
    output = run(valley_dir, tmp_path / "c.geojson", *gauge_station(), "--rivers", str(rivers))
    assert re.search(r"\bplaced on the river line 16 m from the station\b", output), output
    assert "(river Nea, line 8841)" in output, output
    assert re.search(r"\bmoved 4 m onto the DEM's valley floor\b", output), output
    assert "farther on" not in output, output
    assert re.search(r"\bkm(2|²) drain through it\b", output), output
    assert re.search(r"\bwithin 30 m up and down the river\b", output), output
    assert re.search(r"\bwell defined\b", output), output


def test_rivers_writes_the_placed_catchment(tmp_path: Path, valley_dir: Path, rivers: Path) -> None:
    out = tmp_path / "c.geojson"
    run(valley_dir, out, *gauge_station(), "--rivers", str(rivers))
    area = square_area(filled(placed_mask())) * gf.D * gf.D
    polygon = read_domain(out).polygon
    assert polygon.area == pytest.approx(area, rel=1e-9)
    assert polygon.contains(Point(gf.lat(gf.CC, 60)))  # above the embankment


def test_rivers_writes_the_placement_and_the_sensitivity(
    tmp_path: Path, valley_dir: Path, rivers: Path
) -> None:
    out = tmp_path / "c.geojson"
    run(valley_dir, out, *gauge_station(), "--rivers", str(rivers))
    props = feature(out)["properties"]
    assert props["placed_on"] == "any"
    assert props["objectid"] == NEA
    assert props["elvid"] == "2-11-1"
    assert props["node_offset_m"] == pytest.approx(4.0, abs=1e-6)
    assert props["lowered_nodes"] == len(gf.DAM_ROWS)
    assert props["reach_up_m"] == pytest.approx(1000.0, abs=1e-6)
    assert props["downstream_checked"] == "whole"
    assert props["swing"] == pytest.approx(expected_swing(), rel=1e-9)
    assert props["causes"] == []


def test_rivers_reports_the_copies_dropped(tmp_path: Path, valley_dir: Path, rivers: Path) -> None:
    output = run(valley_dir, tmp_path / "c.geojson", *gauge_station(), "--rivers", str(rivers))
    assert re.search(r"\brivers: 2 segments read, 1 exact copy dropped\b", output), output


def test_a_station_in_wgs84_is_placed_the_same(
    tmp_path: Path, valley_dir: Path, rivers: Path
) -> None:
    utm, wgs = tmp_path / "utm.geojson", tmp_path / "wgs.geojson"
    run(valley_dir, utm, *gauge_station(), "--rivers", str(rivers))
    lon, lat_ = to_4326(*STATION)
    run(valley_dir, wgs, "--seed", repr(lon), repr(lat_), "--rivers", str(rivers))
    assert read_domain(wgs).polygon.equals_exact(read_domain(utm).polygon, 0.0)


def test_reach_up_is_passed_on(tmp_path: Path, valley_dir: Path, rivers: Path) -> None:
    out = tmp_path / "c.geojson"
    run(valley_dir, out, *gauge_station(), "--rivers", str(rivers), "--reach-up", "200")
    assert feature(out)["properties"]["reach_up_m"] == pytest.approx(200.0, abs=1e-6)


def test_without_rivers_nothing_is_placed(tmp_path: Path, valley_dir: Path) -> None:
    out = tmp_path / "c.geojson"
    output = run(valley_dir, out, *gauge_station(*gf.lat(gf.CC, 150)))
    assert "placed on" not in output
    assert "placed_on" not in feature(out)["properties"]


def test_no_line_within_the_map_radius_is_refused(
    tmp_path: Path, valley_dir: Path, rivers: Path
) -> None:
    output = refused(
        tmp_path, "--dem", str(valley_dir), *gauge_station(), "--rivers", str(rivers),
        "--map-radius", "15",
    )  # fmt: skip
    assert re.search(r"\bno mapped river line within 15 m\b", output), output


def test_the_map_radius_defaults_to_500_m(tmp_path: Path, valley_dir: Path) -> None:
    far = rivers_file(tmp_path / "far.geojson", list(gf.column_line(gf.CC, 0, 5)))
    output = refused(tmp_path, "--dem", str(valley_dir), *gauge_station(), "--rivers", str(far))
    assert re.search(r"\bno mapped river line within 500 m\b", output), output


@pytest.mark.parametrize("option", ["--map-radius", "--reach-up"])
@pytest.mark.parametrize("value", ["0", "-5", "nan", "inf"])
def test_a_bad_radius_or_reach_is_refused(
    tmp_path: Path, valley_dir: Path, rivers: Path, option: str, value: str
) -> None:
    output = refused(
        tmp_path, "--dem", str(valley_dir), *gauge_station(), "--rivers", str(rivers),
        option, value,
    )  # fmt: skip
    assert option in output, output


def test_rivers_with_lakes_is_refused(tmp_path: Path, valley_dir: Path, rivers: Path) -> None:
    lake = lake_geojson(tmp_path / "l.geojson", lake=Point(STATION).buffer(50.0))
    output = refused(
        tmp_path, "--dem", str(valley_dir), *gauge_station(), "--rivers", str(rivers),
        "--lakes", str(lake),
    )  # fmt: skip
    assert re.search(r"(?i)lake", output), output


def test_a_river_file_in_another_crs_is_refused(tmp_path: Path, valley_dir: Path) -> None:
    to_3035 = reprojector(gf.EPSG, "EPSG:3035")
    line = [tuple(p) for p in to_3035(list(gf.column_line(LINE_COL, 20, 230)))]
    other = rivers_file(tmp_path / "r3035.geojson", line, crs="EPSG:3035")
    output = refused(tmp_path, "--dem", str(valley_dir), *gauge_station(), "--rivers", str(other))
    assert re.search(r"\bCRS\b", output), output
    assert "--rivers" in output, output


def test_map_radius_without_rivers_is_refused(tmp_path: Path, valley_dir: Path) -> None:
    output = refused(tmp_path, "--dem", str(valley_dir), *gauge_station(), "--map-radius", "100")
    assert "--rivers" in output, output
