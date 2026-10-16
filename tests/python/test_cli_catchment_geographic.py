"""`rasputin catchment` on a DEM in degrees, end to end (increment 40).

`docs/increments/40-geographic-catchment.md`, sections 5.5 to 5.7 and the red
suites R1, R2, R8 (through the command), R10, R11, R12 and R13.

- The ANADEM fixture (`tests/fixtures/velhas/anadem_velhas.tif`, EPSG:4674)
  from the seed of its reference catchment (`tests/fixtures/velhas/NOTICE`,
  section 2: 597180 7852380 in EPSG:31983). Section 3 measured, flooding the
  native grid: 12,725 nodes, not cut by the edge or NoData; 10.7441 km² by
  the geodesic area of the traced outline (10.7445 km² row-weighted); shared
  area against the reference 0.9903 of it and 0.9864 of ours (fine outline).
- R11: a projected run writes the bytes it wrote before increment 40. Each
  file's SHA-256 was measured on the red commit's tree (master's catchment
  code with increment 38 merged) on darwin arm64, where the suite is run
  before review. The written polygon is the reduced ring, which `_core`'s
  `reduce_ring` computes in C++, and the areas come from GEOS: another
  compiler or CPU may round them differently, so on another platform the
  test is skipped, naming the one it was measured on.

Scales: degrees to five decimals (1e-5°, about 1 m); areas relative 1e-3 on
10.7 km² (the design's bound for R2) and 1e-3 on the valley's 0.5 km²; node
offsets to 0.3 m on a 10 m lattice.
"""

from __future__ import annotations

import hashlib
import json
import platform
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from pyproj import Transformer
from shapely.geometry import shape

import batch_fixtures as bf
import gauge_fixtures as gf
from cli_helpers import USAGE, invoke, mesh_to_vtk, plain, ran
from geodesy_fixtures import TWIN_DX, TWIN_DY, TWIN_LAT0, TWIN_LON0, metres_per_degree, to_4326
from geographic_fixtures import geographic_keys, tiff
from geotiff_fixtures import KARTVERKET, needs_codecs
from mosaic_fixtures import quadrants
from test_cli_catchment import STATION, gauge_station, rivers_file
from test_cli_mesh_mosaic import write_tiles
from tin_engine.io.domain_file import read_domain

VELHAS = Path(__file__).resolve().parents[1] / "fixtures" / "velhas"
ANADEM = VELHAS / "anadem_velhas.tif"
VELHAS_SEED = ("--seed", "597180", "7852380", "--seed-crs", "EPSG:31983")
#: The fixture's node box (`tests/fixtures/velhas/NOTICE`, section 1), padded
#: by a node: about 44.13 to 44.05 W, 19.40 to 19.49 S.
VELHAS_BOX = (-44.1355, -19.4878, -44.0475, -19.4000)
#: The suggested projected CRS, and the geographic mesh suite's tolerance.
TARGET, MESH_TOLERANCE = "EPSG:31983", "5"

#: R11: (platform, machine) -> {run: SHA-256 of the written file}, measured on
#: the red commit's tree. The DTM10 seed is node (row 4191, column 2456) of
#: `7908_3_10m_z33.tif`, a stream at 120 m whose catchment (10,448 nodes, two
#: windows) is contained; the valley run is `test_cli_catchment.py`'s.
R11_SHA256 = {
    ("darwin", "arm64"): {
        "dtm10": "fc11223e2c9035559f449f5edca860772a72e89a15ce1e88d9457dba48f9f5dd",
        "valley": "ea6e28b20a00326df9cc18bc542b0dc65f5f48cb54a8ac8d5cf338bff9abb0df",
    },
}
DTM10_SEED = ("--seed", "824310", "7908340", "--seed-crs", "EPSG:25833")


def feature(path: Path) -> dict[str, Any]:
    doc = json.loads(path.read_text())
    (f,) = doc["features"]
    return {"crs": doc["crs"]["properties"]["name"], **f}


@pytest.fixture(scope="module")
def velhas(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, str]:
    """R1's run: the catchment file and stderr."""
    out = tmp_path_factory.mktemp("velhas") / "c.geojson"
    code, output = invoke("catchment", "--dem", str(ANADEM), *VELHAS_SEED, "--out", str(out))
    assert code == 0, output
    return out, output


# ---------------------------------------------------------------- R1, R2


def test_r1_a_dem_in_degrees_is_flooded_on_its_own_grid(velhas: tuple[Path, str]) -> None:
    out, _ = velhas
    f = feature(out)
    assert f["crs"] == "EPSG:4674"
    assert f["properties"]["nodes"] == 12_725
    ring = np.asarray(f["geometry"]["coordinates"][0])
    x0, y0, x1, y1 = VELHAS_BOX
    assert np.all((x0 <= ring[:, 0]) & (ring[:, 0] <= x1)), "longitude first (always_xy)"
    assert np.all((y0 <= ring[:, 1]) & (ring[:, 1] <= y1))


def test_r2_the_velhas_catchment_matches_the_resampled_reference(
    velhas: tuple[Path, str],
) -> None:
    """The written outline in EPSG:31983 against the reference the resampled
    route made: 0.98 of each in the other (the design's floor; the fine
    outline measured 0.9864 and 0.9903). The fine area is the geodesic
    10.7441 km² to 1e-3."""
    out, _ = velhas
    f = feature(out)
    t = Transformer.from_crs("EPSG:4674", TARGET, always_xy=True)
    ours = shapely.transform(
        shape(f["geometry"]), lambda xy: np.column_stack(t.transform(xy[:, 0], xy[:, 1]))
    )
    reference = shape(json.loads((VELHAS / "catchment.geojson").read_text()))
    shared = ours.intersection(reference).area
    assert shared / reference.area >= 0.98
    assert shared / ours.area >= 0.98
    assert f["properties"]["fine_area_m2"] == pytest.approx(10.7441e6, rel=1e-3)
    assert f["properties"]["reduced_area_m2"] == pytest.approx(10.7441e6, rel=1e-3)


def test_r2_the_default_tolerance_is_twice_the_cell_in_metres(velhas: tuple[Path, str]) -> None:
    """5.6: `2 * cell_m`, the cell's longer side (north-south, 29.8324 m at
    19.4436° S), not twice the spacing in degrees (0.000539)."""
    out, _ = velhas
    assert feature(out)["properties"]["outline_tolerance_m"] == pytest.approx(59.6648, abs=2e-3)


# ---------------------------------------------------------------- R12, R13


def test_r12_the_file_feeds_mesh(velhas: tuple[Path, str], tmp_path: Path) -> None:
    """`mesh --dem <ANADEM> --out-crs <suggested> --domain <R1's file>`: the
    domain is read in its own CRS, EPSG:4674, and moved."""
    out, _ = velhas
    domain = read_domain(out).polygon
    assert domain.is_valid
    vtk = mesh_to_vtk(
        tmp_path,
        *("--dem", str(ANADEM), "--domain", str(out)),
        *("--out-crs", TARGET, "--tolerance", MESH_TOLERANCE),
    )
    assert len(vtk.points) > 3


def test_r13_the_report_speaks_degrees_and_km2(velhas: tuple[Path, str]) -> None:
    _, output = velhas
    deg = r"-?\d+\.\d{5}"
    assert re.search(
        rf"\bwindow 1: longitude {deg} to {deg}, latitude {deg} to {deg}, \d+ x \d+ nodes\b",
        output,
    ), output
    assert re.search(r"\bcatchment: 12725 DEM nodes, 10\.74\d* km2\b", output), output
    assert re.search(r"\boutline along DEM cells: \d+ vertices, 10\.74\d* km2\b", output), output
    assert re.search(r"\breduced outline: \d+ vertices, 10\.74\d* km2\b", output), output
    assert re.search(r"\btolerance 59\.66\d* m\b", output), output
    assert "np.float64" not in output, output
    assert "x -44" not in output, output


# ---------------------------------------------------------------- R10


def test_r10_station_catchments_refuses_a_dem_in_degrees(tmp_path: Path) -> None:
    """5.7 and question Q3, in its words; exit 2 and nothing written. The
    stations and rivers name the DEM's CRS, so no other refusal applies."""
    stations = bf.write_stations(tmp_path / "s.geojson", crs="EPSG:4674")
    rivers = bf.write_rivers(tmp_path / "r.geojson", crs="EPSG:4674")
    out = tmp_path / "out"
    code, output = invoke(
        "station-catchments", "--dem", str(ANADEM), "--stations", str(stations),
        "--rivers", str(rivers), "--out-dir", str(out),
    )  # fmt: skip
    words = " ".join(plain(output).split())
    assert code == USAGE, words
    assert "station-catchments needs a DEM in metres (a projected CRS)" in words, words
    assert "for a DEM in degrees, run catchment for each station" in words, words
    assert not out.exists()


# ---------------------------------------------------------------- R11


def measured(run: str) -> str:
    key = (sys.platform, platform.machine())
    if key not in R11_SHA256:
        pytest.skip(f"R11's figures were measured on darwin arm64, not on {key}")
    return R11_SHA256[key][run]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@needs_codecs
def test_r11_a_projected_dtm10_run_writes_todays_bytes(tmp_path: Path) -> None:
    out = tmp_path / "c.geojson"
    output = ran("catchment", "--dem", str(KARTVERKET), *DTM10_SEED, "--out", str(out))
    assert feature(out)["properties"]["nodes"] == 10_448
    assert len(feature(out)["properties"]["windows"]) == 2
    assert sha256(out) == measured("dtm10")
    # The start line's point is two plain numbers (5.6), on a projected DEM too.
    assert "np.float64" not in output, output
    assert re.search(r"\bthe outlet node at \(824310(\.0)?, 7908340(\.0)?\)", output), output


def valley_dir(root: Path) -> Path:
    """`test_cli_catchment.py`'s DEM: the dammed valley in four tiles."""
    tile = gf.tile_of(gf.valley(dam=True))
    write_tiles(root / "valley", quadrants(tile, row_cut=120, col_cut=30, overlap=1))
    return root / "valley"


def test_r11_a_projected_rivers_run_writes_todays_bytes(tmp_path: Path) -> None:
    rivers = rivers_file(tmp_path / "rivers.geojson")
    out = tmp_path / "c.geojson"
    ran(
        "catchment", "--dem", str(valley_dir(tmp_path)), *gauge_station(),
        "--rivers", str(rivers), "--out", str(out),
    )  # fmt: skip
    assert sha256(out) == measured("valley")


# ---------------------------------------------------------------- R8 through the command


def twin_tif(root: Path) -> Path:
    """The dammed valley's twin in EPSG:4326 (`geodesy_fixtures`), one tile."""
    stream = tiff(
        gf.valley(dam=True),
        x0=TWIN_LON0,
        y0=TWIN_LAT0,
        step_x=TWIN_DX,
        step_y=TWIN_DY,
        shorts=geographic_keys(4326, area=False),
        nodata=None,
    )
    path = root / "twin.tif"
    path.write_bytes(stream.getvalue())
    return path


def test_r8_a_river_file_in_metres_over_a_dem_in_degrees(tmp_path: Path) -> None:
    """The station and the river file in EPSG:25833, the DEM the valley's
    twin in degrees: the gauge is placed in metres with `U`'s floor three of
    the DEM's cells, and the catchment is the projected valley's, node for
    node.

    The cell is measured at the station (5.5), not at the grid's centre: the
    twin's cells are exactly 10 m only at row 120, and the station is at row
    150, 300 m south, where a cell is 10.0008 m east-west. The expected floor
    is three of that cell, from `metres_per_degree` at the station's latitude,
    to 1e-6 m on 30 m: a cell measured at row 120 gives 30.0000 m, and one
    measured a row (10 m) north or south of the station moves the floor by
    8e-5 m, so either fails."""
    rivers = rivers_file(tmp_path / "rivers.geojson")
    projected, twin = tmp_path / "p.geojson", tmp_path / "t.geojson"
    ran(
        "catchment", "--dem", str(valley_dir(tmp_path)), *gauge_station(),
        "--rivers", str(rivers), "--out", str(projected),
    )  # fmt: skip
    output = ran(
        "catchment", "--dem", str(twin_tif(tmp_path)), *gauge_station(),
        "--rivers", str(rivers), "--out", str(twin),
    )  # fmt: skip
    assert re.search(r"\bwithin 30 m up and down the river\b", output), output
    p, t = feature(projected), feature(twin)
    assert t["crs"] == "EPSG:4326"
    assert t["properties"]["nodes"] == p["properties"]["nodes"]
    assert t["properties"]["node_offset_m"] == pytest.approx(4.0, abs=0.3)
    assert t["properties"]["fine_area_m2"] == pytest.approx(
        p["properties"]["fine_area_m2"], rel=1e-3
    )
    kx, ky = metres_per_degree(to_4326(*STATION)[1])
    cell_m = max(TWIN_DX * kx, TWIN_DY * ky)
    assert cell_m == pytest.approx(10.0008, abs=1e-4)
    assert t["properties"]["uncertainty_m"] == pytest.approx(3 * cell_m, abs=1e-6)
