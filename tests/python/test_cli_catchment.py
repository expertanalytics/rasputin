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
- stderr has a line per window, and the fine outline's vertex count and area.
  The wording is not pinned; these tests look for the words `window`, `fine`,
  `vertices` and `km`.
- Every refusal is a non-zero exit and no file.

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
from shapely.geometry import Point
from typer.testing import CliRunner

from catchment_fixtures import EPSG, D, bowl, filled, full_flood, lake_box, lat, seeds_in, tile_of
from gpkg_fixtures import DTM10, OLA_NORWAY, Layer, Row, write_gpkg
from mosaic_fixtures import quadrants
from test_catchment import moved, to_4326
from test_cli_mesh import plain
from test_cli_mesh_mosaic import write_tiles
from test_outline import square_area
from tin_engine.cli import app
from tin_engine.crs import parse_crs, reprojector
from tin_engine.domain import read_domain

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

SEED = lat(100, 200)  # the bowl's centre node, inside the lake
NVE_BYGDIN_KM2 = 305.54  # NVE delfelt 1187, fetched 2026-09-29 (see the design)
BYGDIN_SEED = ("8.5425", "61.3512")


def invoke(*args: str) -> tuple[int, str]:
    result = runner.invoke(app, list(args))
    return result.exit_code, plain(result.output)


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
    code, output = invoke("catchment", "--dem", str(dem), *args, "--out", str(out))
    assert code == 0, output
    assert out.is_file()
    return output


def refused(tmp_path: Path, *args: str, out: str = "c.geojson") -> str:
    target = tmp_path / out
    code, output = invoke("catchment", *args, "--out", str(target))
    assert code != 0, output
    # Refused by the command, not by Typer for want of one.
    assert "No such command" not in output, output
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
    assert domain.polygon.area == pytest.approx(expected_area(), rel=1e-12)


def test_mesh_accepts_it_as_a_domain(tmp_path: Path, dem_dir: Path, lakes: Path) -> None:
    out = tmp_path / "c.geojson"
    run(dem_dir, out, *utm_seed(), "--lakes", str(lakes))
    vtk = tmp_path / "m.vtk"
    code, output = invoke(
        "mesh", "--dem", str(dem_dir), "--domain", str(out), "--tolerance", "5", "--out", str(vtk)
    )
    assert code == 0, output
    assert vtk.is_file() and vtk.stat().st_size > 0


def test_stderr_reports_the_windows_and_the_fine_outline(
    tmp_path: Path, dem_dir: Path, lakes: Path
) -> None:
    output = run(dem_dir, tmp_path / "c.geojson", *utm_seed(), "--lakes", str(lakes))
    assert re.search(r"(?i)\bwindow", output), output
    assert re.search(r"(?i)\bfine\b.*\bvertices\b|\bvertices\b.*\bfine\b", output), output
    assert re.search(r"(?i)\bfine\b.*km|km.*\bfine\b", output), output


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
    assert read_domain(out).polygon.area == pytest.approx(area, rel=1e-12)


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


def test_a_seed_in_no_lake_is_refused(tmp_path: Path, dem_dir: Path, lakes: Path) -> None:
    refused(tmp_path, "--dem", str(dem_dir), *utm_seed(*lat(60, 60)), "--lakes", str(lakes))


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
