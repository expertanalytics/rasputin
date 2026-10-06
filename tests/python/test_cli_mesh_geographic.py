"""`rasputin mesh --out-crs`: a geographic DEM, or a projected one in another CRS,
resampled onto a target grid and checked against its source (increment 15c-2).

`docs/increments/15c-geographic-dem.md`, D1, D5-D8, J1, J2, J10 and "Tests for
@tester" G2, G6-G9, adapted to 23a-1 and 23a-2's read path
(`docs/increments/23-basin-scale.md`):

- G6 and G9 run the whole path through the CLI and check the written mesh
  against the **source** nodes by an independent final check (the shape of
  `basin-piece/run_sweep.py`'s `final_check`): every valid source node
  inside the domain, projected by pyproj's own `Transformer`, located by
  brute force in every output triangle, its error from the file's own
  vertices. 0 over tolerance, interior and strip. The control is the same
  pipeline without phase 2: `open_dem`'s target tile (the Python API),
  meshed by today's projected path, over which the same check finds some.
- G7 (J1): `--out-crs` equal to the DEM's own CRS writes the same bytes.
- G8 (J10): each refusal before any pixel is read (the repository's `load`
  and `load_window` replaced by ones that fail). **The memory-sum refusal of
  D6 is not tested**: B15, ruled (a) on 2026-10-02 after 15c was designed,
  deletes the refusals at half of physical memory, and `23-basin-scale.md`
  ("What survives of 15c") says 15c's memory cap goes with them.
- **Adapted to 23a (cached sources):** a catalogue source in the cache that is
  geographic (23a-1 W9 refused it "until 15c-2") meshes with `--out-crs`
  exactly as the same file by path; and a `NotCached` refusal names the run's
  `--out-crs` in the fetch command (23a-1 W5: "`--out-crs` (once 15c-2 adds
  it)").

PINNED HERE, where D1 and D7 are silent: `DemRequest.target_crs: str | None`
(the CLI's `--out-crs`), whose `open_dem` returns the **target** tile in
`DemInput.tile` (D1: "DemInput(tile = target tile, ...)"); a geographic DEM
without it is a `ValueError` whose message carries D8's suggestion. D7's
fields, as increment 25 renamed and moved them (`docs/increments/25-plain-output.md`
D2, D4): the mesh file's `crs` is the target's label; the `--stats` rows
`dem_crs` (was `source_crs`) the DEM's, `dem_transform` (was
`source_transform`) pyproj's `description` of the source-to-target
`always_xy` transformer, `resampled_grid` (was `computation_grid`) beginning
`<h> m square grid in <crs>` (the design's example), and `dem_nodes_checked`
the store's size (was `checked against <N> source nodes`). The file's
`max_error_m` is the error at the source nodes, not the resampled grid's
(which is `resampled_grid_max_error_m` in `--stats`).

HOW THIS FILE GOES RED: `--out-crs` is not an option of `mesh`, so every run
exits 2 with "No such option" (the `refused` helper rejects that reason), and
`DemRequest` ignores `target_crs`, so the Python-API control is refused by the
reader. G7's comparison and the catchment refusal are guards. The ANADEM
realism case runs on the fixture committed with the green work (Q17).

AMENDED after `@reviewer`'s round 1 on 15c-2 ("Round 1 findings", below):
B1, `--out-crs` must be a projected CRS in metres, refused in `open_dem`
before any pixel is read; and four cases that kill mutants the first suite
let live (several check-point blocks, the store's size, the refusal's
percentages, `final.ok()` honoured).
"""

from __future__ import annotations

import gc
import importlib
import json
import math
import os
import re
import shlex
import subprocess
import sys
import weakref
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from pyproj import CRS, Transformer
from shapely.geometry import Polygon
from typer.testing import CliRunner

from cli_helpers import ANSI, USAGE, Ring, geojson, invoke, mesh_to_vtk, refused, squashed
from cog_fixtures import write_cache
from geographic_fixtures import (
    ANADEM_STEP,
    LAT0,
    LON0,
    NODATA,
    SourceCheck,
    geographic_tile_tiff,
    nodes,
    project,
    project_ring,
    projected_tiff,
    rough,
    source_check,
)
from geotiff_fixtures import KARTVERKET, needs_codecs
from test_cli_mesh_domain import quarter_circle
from test_cli_mesh_mosaic import same_mesh
from test_cli_mesh_refine import file_field, stats_row
from test_cli_mesh_stats import section, table
from tin_engine.cli import app
from tin_engine.dem_input import DemRequest, open_dem
from tin_engine.domain import read_domain
from tin_engine.io.repository import TiffDemRepository
from vtkread import VtkFile

TARGET = "EPSG:31983"
ROWS = COLS = 60
TOLERANCE = 1.0
H = 30  # ANADEM's 29.8 m at 19 S, rounded (Q16)
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
VELHAS = FIXTURES / "velhas"  # Q17: the ANADEM extract and the DEM-derived catchment


def field(vtk: VtkFile, name: str) -> str:
    return file_field(vtk, name)


def report_of(tmp_path: Path, out: str = "x.vtk") -> str:
    return (tmp_path / out).with_suffix(".md").read_text(encoding="utf-8")


def inputs(tmp_path: Path, name: str, out: str = "x.vtk") -> str:
    """The ``--stats`` row ``name`` of the ``mesh_to_vtk`` that wrote ``out``."""
    return stats_row(report_of(tmp_path, out), name)


def triangles(vtk: VtkFile) -> np.ndarray:
    tri = np.array([np.asarray(p) for p in vtk.polygons])
    assert tri.ndim == 2 and tri.shape[1] == 3
    return tri


@pytest.fixture
def no_pixels(monkeypatch: pytest.MonkeyPatch) -> None:
    """J10: a repository whose `load` and `load_window` fail, so a refusal that
    reads a pixel fails the test instead of passing it."""

    def fail(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError("a pixel was read before the refusal")

    monkeypatch.setattr(TiffDemRepository, "load", fail)
    monkeypatch.setattr(TiffDemRepository, "load_window", fail)


# ------------------------------------------------------------------ the scene

#: The domain's ring in the source's (row, col) space: inside the tile with
#: margin, every vertex off-node.
DOMAIN_RC: list[tuple[float, float]] = [(8.3, 7.7), (9.1, 51.2), (50.6, 52.4), (49.2, 8.8)]


def tile_array() -> np.ndarray:
    array = rough(ROWS, COLS)
    array[[20, 33, 41], [17, 30, 44]] = NODATA  # a few voids inside the domain
    return array


def lonlat_ring(rc: list[tuple[float, float]]) -> Ring:
    return [(LON0 + c * ANADEM_STEP, LAT0 - r * ANADEM_STEP) for r, c in rc]


@pytest.fixture
def geographic_dem(tmp_path: Path) -> Path:
    path = tmp_path / "anadem_like.tif"
    path.write_bytes(geographic_tile_tiff(tile_array()).getvalue())
    return path


@pytest.fixture
def domain_4674(tmp_path: Path) -> Path:
    """The domain in a third CRS (SIRGAS 2000 over a WGS 84 DEM): Degeneracy
    policy, "the domain in a third CRS"."""
    ring = project_ring("EPSG:4326", "EPSG:4674", lonlat_ring(DOMAIN_RC))
    return geojson(tmp_path / "domain.geojson", ring, crs="EPSG:4674")


def geographic_check(vtk: VtkFile, domain: Path, tolerance: float = TOLERANCE) -> SourceCheck:
    """The independent final check of a mesh in TARGET against the source tile."""
    array = tile_array().ravel().astype(np.float64)
    lon, lat = nodes(LON0, LAT0, ANADEM_STEP, ANADEM_STEP, (ROWS, COLS))
    valid = array != NODATA
    xy = project("EPSG:4326", TARGET, lon[valid], lat[valid])
    given = json.loads(domain.read_text())
    ring = project_ring(given["crs"]["properties"]["name"], TARGET, given["coordinates"][0][:-1])
    return source_check(vtk.points, triangles(vtk), xy, array[valid], Polygon(ring), tolerance, H)


# ------------------------------------------------------------------ G6


class TestEndToEnd:
    """G6: a synthetic geographic tile with a rough surface, a domain in
    EPSG:4674, `--out-crs EPSG:31983`."""

    def test_the_mesh_is_in_the_target_crs_and_records_d7s_fields(
        self, tmp_path: Path, geographic_dem: Path, domain_4674: Path
    ) -> None:
        vtk = mesh_to_vtk(
            tmp_path,
            *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
            *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
        )
        assert field(vtk, "crs") == TARGET
        assert inputs(tmp_path, "dem_crs") == "EPSG:4326"
        expected = Transformer.from_crs("EPSG:4326", TARGET, always_xy=True).description
        assert inputs(tmp_path, "dem_transform") == expected
        assert inputs(tmp_path, "domain_crs") == "EPSG:4674"
        grid = inputs(tmp_path, "resampled_grid")
        assert grid.startswith(f"{H} m square grid in {TARGET}"), grid
        for moved in ("source_crs", "source_transform", "computation_grid", "elevation_source"):
            assert moved not in vtk.field_data, moved
        for renamed in ("dem_crs", "dem_transform", "resampled_grid", "dem_nodes_checked"):
            assert renamed not in vtk.field_data, f"{renamed} is a --stats row (D2)"
        checked = int(inputs(tmp_path, "dem_nodes_checked"))
        found = geographic_check(vtk, domain_4674)
        assert found.nodes <= checked <= ROWS * COLS - 3
        # Every vertex inside the domain, moved to the target CRS.
        ring = project_ring(
            "EPSG:4674", TARGET, json.loads(domain_4674.read_text())["coordinates"][0][:-1]
        )
        grown = Polygon(ring).buffer(1e-3)
        assert shapely.intersects_xy(grown, vtk.points[:, 0], vtk.points[:, 1]).all()

    def test_no_source_node_is_over_tolerance_interior_or_strip(
        self, tmp_path: Path, geographic_dem: Path, domain_4674: Path
    ) -> None:
        """J2 at the source nodes, by the independent check."""
        vtk = mesh_to_vtk(
            tmp_path,
            *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
            *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
        )
        found = geographic_check(vtk, domain_4674)
        assert found.nodes > 1000
        assert (found.over_interior, found.over_strip) == (0, 0), found

    def test_without_phase_2_the_check_finds_some(
        self, tmp_path: Path, geographic_dem: Path, domain_4674: Path
    ) -> None:
        """The control: `open_dem`'s target tile meshed by today's projected
        path (phase 1 alone) leaves source nodes over tolerance."""
        opened = open_dem(
            DemRequest(
                sources=(geographic_dem,),
                domain=read_domain(domain_4674),
                target_crs=TARGET,  # type: ignore[call-arg]
            )
        )
        m = opened.tile.meta
        assert (m.crs, m.delta_x, m.x_min % H, m.y_max % H) == (TARGET, H, 0.0, 0.0)  # type: ignore[attr-defined]
        grid = tmp_path / "target.tif"
        grid.write_bytes(
            projected_tiff(
                np.asarray(opened.tile.array), x0=m.x_min, y0=m.y_max, step=H, epsg=31983
            ).getvalue()
        )
        given = json.loads(domain_4674.read_text())
        moved = geojson(
            tmp_path / "domain_31983.geojson",
            project_ring("EPSG:4674", TARGET, given["coordinates"][0][:-1]),
            crs=TARGET,
        )
        vtk = mesh_to_vtk(
            tmp_path, "--dem", str(grid), "--domain", str(moved), "--tolerance", str(TOLERANCE)
        )
        found = geographic_check(vtk, domain_4674)
        assert found.over > 0, found

    def test_the_phase_2_mesh_differs_from_phase_1s_by_the_inserted_nodes(
        self, tmp_path: Path, geographic_dem: Path, domain_4674: Path
    ) -> None:
        """J8: every vertex phase 2 adds carries a source node's own z: each
        output vertex is a target-grid node, a domain vertex, or a source node
        at its projected position (to the store's 2 um) with its value."""
        vtk = mesh_to_vtk(
            tmp_path,
            *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
            *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
        )
        array = tile_array().ravel().astype(np.float64)
        lon, lat = nodes(LON0, LAT0, ANADEM_STEP, ANADEM_STEP, (ROWS, COLS))
        valid = array != NODATA
        xy = project("EPSG:4326", TARGET, lon[valid], lat[valid])
        on_grid = (vtk.points[:, 0] % H == 0) & (vtk.points[:, 1] % H == 0)
        off = vtk.points[~on_grid]
        distance = np.hypot(off[:, None, 0] - xy[None, :, 0], off[:, None, 1] - xy[None, :, 1])
        nearest = distance.argmin(axis=1)
        is_source = distance[np.arange(len(off)), nearest] < 1e-5
        assert is_source.sum() > 0, "phase 2 inserted nothing on a rough surface"
        assert (off[is_source, 2] == array[valid][nearest[is_source]]).all()


# ------------------------------------------------------------------ 15e, fix 3


def test_the_target_tile_is_dropped_before_the_final_check(
    tmp_path: Path,
    geographic_dem: Path,
    domain_4674: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Increment 15e fix 3 (`docs/increments/15e-memory-fixes.md`): no
    reference to the resampled target tile survives into
    `final_check.run`, so phase 2 does not hold it. The probe is a weakref
    to the tile's array (a NumPy array takes one), read after a full
    `gc.collect()` on entry to the final check. Went red at 9879805 because
    `mesh()` held `opened`, and `_dem_mesh` its `tile` parameter, through the call."""
    from tin_engine import cli

    tiles: list[weakref.ref[np.ndarray]] = []
    dead_on_entry: list[bool] = []
    real_open, real_run = cli.open_dem, cli.final_check.run

    def opening(*args: Any, **kwargs: Any) -> Any:
        result = real_open(*args, **kwargs)
        tiles.append(weakref.ref(result.tile.array))
        return result

    def checking(*args: Any, **kwargs: Any) -> Any:
        gc.collect()
        dead_on_entry.extend(ref() is None for ref in tiles)
        return real_run(*args, **kwargs)

    monkeypatch.setattr(cli, "open_dem", opening)
    monkeypatch.setattr(cli.final_check, "run", checking)
    mesh_to_vtk(
        tmp_path,
        *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
        *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
    )
    assert len(tiles) == 1, "open_dem was not called once"
    assert dead_on_entry == [True], "the target tile is alive when the final check starts"


# ------------------------------------------------------------------ G9


class TestAProjectedDemInAnotherCrs:
    """G9 (Q12 (a)): EPSG:25833 at 10 m, `--out-crs EPSG:25832`."""

    X0, Y0, STEP = 640_000.0, 6_650_000.0, 10.0

    def test_takes_the_reprojected_path_with_the_sources_spacing(self, tmp_path: Path) -> None:
        array = rough(ROWS, COLS, seed=9)
        dem = tmp_path / "utm33.tif"
        dem.write_bytes(
            projected_tiff(array, x0=self.X0, y0=self.Y0, step=self.STEP, epsg=25833).getvalue()
        )
        ring = [(self.X0 + c * self.STEP, self.Y0 - r * self.STEP) for r, c in DOMAIN_RC]
        domain = geojson(tmp_path / "d.geojson", ring, crs="EPSG:25833")
        vtk = mesh_to_vtk(
            tmp_path,
            *("--dem", str(dem), "--domain", str(domain)),
            *("--out-crs", "EPSG:25832", "--tolerance", str(TOLERANCE)),
        )
        assert field(vtk, "crs") == "EPSG:25832"
        assert inputs(tmp_path, "dem_crs") == "EPSG:25833"
        grid = inputs(tmp_path, "resampled_grid")
        assert grid.startswith("10 m square grid in EPSG:25832"), grid
        assert int(inputs(tmp_path, "dem_nodes_checked")) > 1000
        x, y = nodes(self.X0, self.Y0, self.STEP, self.STEP, (ROWS, COLS))
        xy = project("EPSG:25833", "EPSG:25832", x, y)
        moved = Polygon(project_ring("EPSG:25833", "EPSG:25832", ring))
        found = source_check(
            vtk.points, triangles(vtk), xy, array.ravel().astype(np.float64), moved, TOLERANCE, 10.0
        )
        assert found.nodes > 1000
        assert (found.over_interior, found.over_strip) == (0, 0), found


# ------------------------------------------------------------------ G7


@needs_codecs
def test_out_crs_equal_to_the_dems_own_writes_the_same_bytes(tmp_path: Path) -> None:
    """J1: the committed projected tile with a domain, `--out-crs` its own CRS."""
    domain = geojson(tmp_path / "quarter.geojson", quarter_circle(), crs="EPSG:25833")
    args = ("--dem", str(KARTVERKET), "--domain", str(domain), "--tolerance", "10")
    plain_code, plain_output = invoke("mesh", *args, "--out", str(tmp_path / "plain.vtk"))
    assert plain_code == 0, plain_output
    code, output = invoke(
        "mesh", *args, "--out-crs", "EPSG:25833", "--out", str(tmp_path / "same.vtk")
    )
    assert code == 0, output
    assert (tmp_path / "same.vtk").read_bytes() == (tmp_path / "plain.vtk").read_bytes()


# ------------------------------------------------------------------ G8


class TestRefusalsBeforePixels:
    def test_a_geographic_dem_without_out_crs_suggests_one(
        self, tmp_path: Path, geographic_dem: Path, domain_4674: Path, no_pixels: None
    ) -> None:
        """Q11 (a): the refusal names `--out-crs` and prints D8's suggestion for
        the domain's box in the DEM's CRS, with its worst distortion."""
        from tin_engine.crs import suggest_crs  # type: ignore[attr-defined]

        given = json.loads(domain_4674.read_text())["coordinates"][0][:-1]
        lonlat = project_ring("EPSG:4674", "EPSG:4326", given)
        lons, lats = [p[0] for p in lonlat], [p[1] for p in lonlat]
        suggestion = suggest_crs((min(lons), max(lons), min(lats), max(lats)), "EPSG:4326")
        refused(
            tmp_path,
            "mesh",
            *("--dem", str(geographic_dem), "--domain", str(domain_4674), "--tolerance", "1"),
            says=("--out-crs", suggestion.proj, suggestion.family, "worst scale error"),
            squash=True,
        )

    def test_the_python_api_refuses_it_too(
        self, geographic_dem: Path, domain_4674: Path, no_pixels: None
    ) -> None:
        request = DemRequest(sources=(geographic_dem,), domain=read_domain(domain_4674))
        with pytest.raises(ValueError, match=r"\+proj="):
            open_dem(request)

    def test_a_source_region_across_180_degrees(self, tmp_path: Path, no_pixels: None) -> None:
        """A tile reaching past 180 E, a domain straddling it, and a target CRS
        centred there: the source region's image in degrees crosses ±180."""
        dem = tmp_path / "fiji.tif"
        dem.write_bytes(
            geographic_tile_tiff(rough(40, 200), lon0=179.9, lat0=-17.0, step=0.001).getvalue()
        )
        ring = [(179.97, -17.031), (180.03, -17.032), (180.031, -17.009), (179.971, -17.008)]
        domain = geojson(tmp_path / "d.geojson", ring, crs="EPSG:4326")
        output = refused(
            tmp_path,
            "mesh",
            *("--dem", str(dem), "--domain", str(domain), "--tolerance", "1"),
            *("--out-crs", "+proj=tmerc +lat_0=0 +lon_0=180 +k=1 +ellps=WGS84 +units=m"),
            says=("180",),
            squash=True,
        )
        assert re.search(r"(?i)antimeridian|longitude", output), output

    def test_a_domain_vertex_without_an_image(
        self, tmp_path: Path, geographic_dem: Path, no_pixels: None
    ) -> None:
        """A vertex on the far side of the globe from an orthographic target."""
        ring = [(-44.01, -19.01), (-43.99, -19.01), (136.0, 19.0)]
        domain = geojson(tmp_path / "d.geojson", ring, crs="EPSG:4326")
        refused(
            tmp_path,
            "mesh",
            *("--dem", str(geographic_dem), "--domain", str(domain), "--tolerance", "1"),
            *("--out-crs", "+proj=ortho +lat_0=-19 +lon_0=-44 +ellps=WGS84 +units=m"),
            says=("no image",),
            squash=True,
        )


# ------------------------------------------------------------------ G2


def test_catchment_on_a_geographic_tile_is_a_usage_error(
    tmp_path: Path, geographic_dem: Path
) -> None:
    """The `to_core` gate keeps a geographic DEM out of `rasputin catchment`
    (increment 22), which 15c does not extend."""
    lon, lat = LON0 + 30 * ANADEM_STEP, LAT0 - 30 * ANADEM_STEP
    code, output = invoke(
        "catchment",
        *("--dem", str(geographic_dem), "--seed", repr(lon), repr(lat), "--seed-crs", "EPSG:4326"),
        *("--out", str(tmp_path / "c.geojson")),
    )
    assert code == USAGE, output
    assert "Traceback" not in output
    assert re.search(r"(?i)geographic", output), output
    assert not (tmp_path / "c.geojson").exists()


# ------------------------------------------------------------------ 23a's cached sources

KEY = "test-geographic"


@pytest.fixture
def catalogue(monkeypatch: pytest.MonkeyPatch) -> None:
    """`SOURCES` plus a geographic test entry, wherever a module bound the name."""
    sources = importlib.import_module("tin_engine.sources")
    entry = sources.RemoteSource(
        id=KEY,
        kind="one-cog",
        url="https://example.invalid/test-geographic.tif",
        crs="EPSG:4326",
        nodata=NODATA,
        credit="Test geographic elevation",
        licence_note="test fixture only",
    )
    original = sources.SOURCES
    patched = {**original, KEY: entry}
    for module in list(sys.modules.values()):
        name = getattr(module, "__name__", "") or ""
        if name.startswith("tin_engine") and getattr(module, "SOURCES", None) is original:
            monkeypatch.setattr(module, "SOURCES", patched)


def cog_bytes() -> bytes:
    return geographic_tile_tiff(
        tile_array(), compression="deflate", tile=(16, 16), overview=True
    ).getvalue()


class TestCachedGeographicSource:
    def test_meshes_as_the_same_file_by_path(
        self, tmp_path: Path, catalogue: None, domain_4674: Path
    ) -> None:
        data = cog_bytes()
        write_cache(tmp_path / "cache", KEY, {"dem": data}, crs="EPSG:4326")
        path = tmp_path / "dem.tif"
        path.write_bytes(data)
        common = ("--domain", str(domain_4674), "--out-crs", TARGET, "--tolerance", "1")
        by_key = mesh_to_vtk(
            tmp_path, "--dem", KEY, "--cache", str(tmp_path / "cache"), *common, out="k.vtk"
        )
        by_path = mesh_to_vtk(tmp_path, "--dem", str(path), *common, out="p.vtk")
        same_mesh(by_key, by_path)

    def test_not_cached_names_the_runs_out_crs(
        self, tmp_path: Path, catalogue: None, domain_4674: Path
    ) -> None:
        write_cache(
            tmp_path / "cache", KEY, {"dem": cog_bytes()}, skip={"dem": [5]}, crs="EPSG:4326"
        )
        refused(
            tmp_path,
            "mesh",
            *("--dem", KEY, "--cache", str(tmp_path / "cache"), "--domain", str(domain_4674)),
            *("--out-crs", TARGET, "--tolerance", "1"),
            says=(f"rasputin fetch {KEY}", f"--out-crs {TARGET}"),
            squash=True,
        )


# ------------------------------------------------------------------ realism (Q17)


def test_the_anadem_extract_over_its_catchment(tmp_path: Path) -> None:
    """G6's realism case: the committed ANADEM extract (EPSG:4674, the COG's
    key set) over the DEM-derived catchment (EPSG:31983), checked against
    every source node inside it."""
    import tifffile

    dem, catchment = VELHAS / "anadem_velhas.tif", VELHAS / "catchment.geojson"
    vtk = mesh_to_vtk(
        tmp_path,
        *("--dem", str(dem), "--domain", str(catchment), "--out-crs", TARGET, "--tolerance", "5"),
    )
    with tifffile.TiffFile(dem) as tif:
        page = tif.pages.first
        (_, _, _, x0, y0, _), (dx, dy, _) = page.tags[33922].value, page.tags[33550].value
        array = page.asarray().astype(np.float64)
    half = 0.5  # area-registered, as ANADEM is
    lon, lat = nodes(x0 + half * dx, y0 - half * dy, dx, dy, array.shape)
    valid = array.ravel() != NODATA
    xy = project("EPSG:4674", TARGET, lon[valid], lat[valid])
    ring = json.loads(catchment.read_text())["coordinates"][0][:-1]
    found = source_check(
        vtk.points, triangles(vtk), xy, array.ravel()[valid], Polygon(ring), 5.0, H
    )
    assert math.isfinite(found.worst) and found.nodes > 0
    assert (found.over_interior, found.over_strip) == (0, 0), found
    assert_max_error_is_the_sources(tmp_path, vtk, found, tolerance=5.0)


def assert_max_error_is_the_sources(
    tmp_path: Path, vtk: VtkFile, found: SourceCheck, tolerance: float
) -> None:
    """Increment 25, D2: on the reprojected path the file's ``max_error_m`` is
    the error at the original DEM's nodes, raised by the at-vertex figure,
    and not phase 1's figure against the resampled grid.

    ``found.worst`` is the independent check's largest error over every
    located source node, those on a vertex included; it agrees with the
    final check to the store's 2 um rounding (measured: 4e-8 m on this
    extract, while the resampled grid's figure differs by 6e-4 m)."""
    report = report_of(tmp_path)
    stated = float(field(vtk, "max_error_m"))
    resampled = float(stats_row(report, "resampled_grid_max_error_m"))
    at_vertices = float(stats_row(report, "dem_nodes_at_vertices_max_error_m"))
    assert int(stats_row(report, "dem_nodes_at_vertices")) >= 0  # measured on this path (D6)
    assert stated == pytest.approx(found.worst, abs=1e-6)
    assert abs(resampled - found.worst) > 1e-4, "the extract no longer tells the two apart"
    assert stated >= at_vertices
    assert field(vtk, "tolerance_m") == str(tolerance).removesuffix(".0")
    assert stated <= max(tolerance, at_vertices)
    sizes = table(section(report, "Sizes"))
    assert "resampled grid nodes" in sizes, sizes
    assert "DEM nodes" not in sizes, "the resampled grid is not the DEM (D4)"


def test_g6_max_error_is_the_source_dems(
    tmp_path: Path, geographic_dem: Path, domain_4674: Path
) -> None:
    """D2 on the synthetic tile: the same rule as the ANADEM extract's."""
    vtk = mesh_to_vtk(
        tmp_path,
        *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
        *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
    )
    found = geographic_check(vtk, domain_4674)
    assert_max_error_is_the_sources(tmp_path, vtk, found, tolerance=TOLERANCE)
    report = report_of(tmp_path)
    inserted = int(stats_row(report, "dem_check_points_inserted"))
    rounds = int(stats_row(report, "dem_check_rounds"))
    assert inserted > 0 and rounds > 0  # 15c-2's phase 2 has work on this tile


# ------------------------------------------------------------------ round 1 findings


class TestOutCrsMustBeProjectedInMetres:
    """B1 (`@reviewer`, round 1): the target CRS is the frame `_core` computes
    in, so it must be projected (no degrees in `_core`, J4) and in metres
    (D6: "every number is in a projected CRS, in metres"). Refused like G8,
    in `open_dem`, before any pixel is read, naming the problem."""

    FEET = "+proj=utm +zone=23 +south +ellps=GRS80 +units=us-ft"

    @pytest.mark.parametrize(
        ("out_crs", "problem"),
        [(FEET, r"(?i)metre"), ("EPSG:4674", r"(?i)projected"), ("EPSG:4326", r"(?i)projected")],
        ids=["us_feet", "geographic_4674", "geographic_4326"],
    )
    def test_refused_naming_the_problem(
        self,
        tmp_path: Path,
        geographic_dem: Path,
        domain_4674: Path,
        no_pixels: None,
        out_crs: str,
        problem: str,
    ) -> None:
        output = refused(
            tmp_path,
            "mesh",
            *("--dem", str(geographic_dem), "--domain", str(domain_4674), "--tolerance", "1"),
            *("--out-crs", out_crs),
            says=("--out-crs",),
            squash=True,
        )
        assert re.search(problem, output), output
        assert not re.search(r"(?i)antimeridian|pole", output), output

    @pytest.mark.parametrize("out_crs", [FEET, "EPSG:4674"], ids=["us_feet", "geographic"])
    def test_open_dem_refuses_before_any_pixel(
        self, geographic_dem: Path, domain_4674: Path, no_pixels: None, out_crs: str
    ) -> None:
        request = DemRequest(
            sources=(geographic_dem,), domain=read_domain(domain_4674), target_crs=out_crs
        )
        with pytest.raises(ValueError, match=r"(?i)projected|metre"):
            open_dem(request)

    def test_a_metric_proj_string_still_meshes(
        self, tmp_path: Path, geographic_dem: Path, domain_4674: Path
    ) -> None:
        """D8's own kind of suggestion: a PROJ string with no EPSG code."""
        proj = (
            "+proj=tmerc +lat_0=0 +lon_0=-44 +k=0.9996 +x_0=500000 +y_0=10000000"
            " +ellps=GRS80 +units=m"
        )
        vtk = mesh_to_vtk(
            tmp_path,
            *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
            *("--out-crs", proj, "--tolerance", str(TOLERANCE)),
        )
        assert int(inputs(tmp_path, "dem_nodes_checked")) > 1000
        assert float(field(vtk, "max_error_m")) <= TOLERANCE


class TestRoundOneMutantKillers:
    def args(self, dem: Path, domain: Path) -> tuple[str, ...]:
        return ("--dem", str(dem), "--domain", str(domain), "--out-crs", TARGET, "--tolerance", "1")

    def test_several_check_point_blocks_are_all_checked(
        self,
        tmp_path: Path,
        geographic_dem: Path,
        domain_4674: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Source blocks of 16 x 16 nodes, so the 60 x 60 tile is 16 blocks:
        a final check that stops after the first block leaves nodes over."""
        monkeypatch.setattr("tin_engine.target_grid.BLOCK", 16)
        vtk = mesh_to_vtk(tmp_path, *self.args(geographic_dem, domain_4674))
        found = geographic_check(vtk, domain_4674)
        assert found.nodes > 1000
        assert (found.over_interior, found.over_strip) == (0, 0), found

    def test_the_store_holds_every_point_yielded(
        self,
        tmp_path: Path,
        geographic_dem: Path,
        domain_4674: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """`dem_nodes_checked` (once `checked against N source nodes`) is the
        store's size, which equals the points `check_point_blocks` yields for
        the same request."""
        monkeypatch.setattr("tin_engine.target_grid.BLOCK", 16)
        opened = open_dem(
            DemRequest(
                sources=(geographic_dem,), domain=read_domain(domain_4674), target_crs=TARGET
            )
        )
        assert opened.checks is not None
        yielded = sum(len(z) for _, z in opened.checks)
        mesh_to_vtk(tmp_path, *self.args(geographic_dem, domain_4674))
        assert int(inputs(tmp_path, "dem_nodes_checked")) == yielded > 1000

    def test_the_refusal_prints_percentages(
        self, tmp_path: Path, geographic_dem: Path, domain_4674: Path, no_pixels: None
    ) -> None:
        """D8's message: the worst errors in per cent, the fields times 100."""
        from tin_engine.crs import suggest_crs

        given = json.loads(domain_4674.read_text())["coordinates"][0][:-1]
        lonlat = project_ring("EPSG:4674", "EPSG:4326", given)
        lons, lats = [p[0] for p in lonlat], [p[1] for p in lonlat]
        s = suggest_crs((min(lons), max(lons), min(lats), max(lats)), "EPSG:4326")
        output = squashed(
            refused(
                tmp_path,
                "mesh",
                *("--dem", str(geographic_dem), "--domain", str(domain_4674), "--tolerance", "1"),
                says=("worst scale error",),
                squash=True,
            )
        )
        number = r"([0-9.]+(?:e[-+]?[0-9]+)?)%"
        scale = re.search(r"worstscaleerror" + number, output)
        areal = re.search(r"worstarealerror" + number, output)
        assert scale is not None and areal is not None, output
        assert float(scale.group(1)) == pytest.approx(100 * s.max_scale_error, rel=1e-2)
        assert float(areal.group(1)) == pytest.approx(100 * s.max_areal_error, rel=1e-2)

    def test_a_failed_final_check_writes_nothing(
        self,
        tmp_path: Path,
        geographic_dem: Path,
        domain_4674: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """`final.ok()` false is a usage error in the engine's words, no file."""
        import tin_engine.final_check as final_check

        real = final_check.run

        class Failed:
            def __init__(self, outcome: Any) -> None:
                self._outcome = outcome

            def ok(self) -> bool:
                return False

            message = "planted final-check failure"

            def __getattr__(self, name: str) -> Any:
                return getattr(self._outcome, name)

        def failing(*args: Any, **kwargs: Any) -> Any:
            outcome, n = real(*args, **kwargs)
            return Failed(outcome), n

        monkeypatch.setattr(final_check, "run", failing)
        refused(
            tmp_path,
            "mesh",
            *self.args(geographic_dem, domain_4674),
            says=("planted final-check failure",),
            squash=True,
        )


# ------------------------------------------------------------------ @perf's 15c-2 acceptance


def test_the_final_checks_timing_rows_are_not_zero(
    tmp_path: Path, geographic_dem: Path, domain_4674: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`@perf`'s 15c-2 acceptance: "final check: scan (parallel)" and "final
    check: split + flip (serial)" always read 0.000. When phase 2 inserts
    points it has scanned and split, so the outcome's `scan_seconds` and
    `split_seconds`, and the clock rows D7 names, are positive. Read from
    the clock `final_check.run` was given, not from the printed table, whose
    three decimals would round a small run's real time to 0."""
    import tin_engine.final_check as final_check

    seen: list[tuple[Any, Any]] = []
    real = final_check.run

    def spying(*args: Any, **kwargs: Any) -> Any:
        outcome, n = real(*args, **kwargs)
        clock = kwargs["clock"] if "clock" in kwargs else args[4]
        seen.append((outcome, clock))
        return outcome, n

    monkeypatch.setattr(final_check, "run", spying)
    mesh_to_vtk(
        tmp_path,
        *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
        *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
    )
    ((outcome, clock),) = seen
    assert outcome.inserted > 0, "phase 2 must have work for this test to mean anything"
    assert outcome.scan_seconds > 0 and outcome.split_seconds > 0
    rows = dict(clock.phases())
    assert rows["final check: scan (parallel)"] > 0, rows
    assert rows["final check: split + flip (serial)"] > 0, rows


# ------------------------------------------------------------------ the suggestion, copied back


SUGGESTED = re.compile(r"--out-crs[\s│]*(['\"])(.*?)\1", re.DOTALL)


def copied_suggestion(output: str) -> str:
    """The CRS the refusal says to copy, exactly as printed: the quoted text
    after `--out-crs` (a line break or panel border may come before the
    opening quote), borders and line breaks inside it included, as a reader
    selecting it in a terminal gets it. Only colour codes are removed."""
    found = SUGGESTED.search(ANSI.sub("", output))
    assert found is not None, f"no quoted --out-crs suggestion in {output!r}"
    return found.group(2)


def mesh_with(tmp_path: Path, out_crs: str) -> None:
    """The suggestion passed back as `--out-crs` meshes, in that CRS."""
    dem, catchment = VELHAS / "anadem_velhas.tif", VELHAS / "catchment.geojson"
    vtk = mesh_to_vtk(
        tmp_path,
        *("--dem", str(dem), "--domain", str(catchment), "--tolerance", "5"),
        *("--out-crs", out_crs),
        out="copied.vtk",
    )
    assert CRS.from_user_input(field(vtk, "crs")) == CRS.from_user_input(out_crs)


class TestTheSuggestionPastesBack:
    """Q11 (a): the refusal prints a CRS to copy. `@reviewer` on 18a47cb: the
    WKT2 suggestion is printed inside Rich's error panel, wrapped mid-token at
    the terminal's width with `│` borders, so what a user copies never parses.
    Run as a user would, at 80 columns, copy what is printed, paste it back."""

    ARGS = (
        "mesh",
        *("--dem", str(VELHAS / "anadem_velhas.tif")),
        *("--domain", str(VELHAS / "catchment.geojson")),
        *("--tolerance", "5"),
    )

    def test_through_the_cli_runner_at_80_columns(self, tmp_path: Path) -> None:
        result = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb", "COLUMNS": "80"}).invoke(
            app, [*self.ARGS, "--out", str(tmp_path / "refused.vtk")]
        )
        assert result.exit_code == USAGE, result.output
        mesh_with(tmp_path, copied_suggestion(result.output))

    def test_in_a_real_terminal_process_at_80_columns(self, tmp_path: Path) -> None:
        env = {**os.environ, "COLUMNS": "80", "LINES": "40"}
        env.pop("NO_COLOR", None)
        code = "import sys; from tin_engine.cli import app; app(prog_name='rasputin')"
        done = subprocess.run(
            [sys.executable, "-c", code, *self.ARGS, "--out", str(tmp_path / "refused.vtk")],
            capture_output=True,
            text=True,
            env=env,
            timeout=300,
            check=False,
        )
        assert done.returncode == USAGE, done.stderr
        mesh_with(tmp_path, copied_suggestion(done.stdout + done.stderr))


class TestASuggestionWithAnApostrophePastesBack:
    """`@reviewer` round 4 on 8b7b0dd: the suggestion is printed as
    `--out-crs '<WKT>'`, and a WKT can hold an apostrophe: EPSG:4266's datum
    is M'poraloko (Gabon). The printed line, parsed as a POSIX shell parses
    it (`shlex.split`, and `bash` itself), must give back the CRS, which
    must mesh. A geographic micro tile in EPSG:4266, 0.5 degrees south of
    the equator."""

    LON, LAT = 10.5, -0.5

    @pytest.fixture
    def scene(self, tmp_path: Path) -> tuple[Path, Path]:
        dem = tmp_path / "gabon.tif"
        dem.write_bytes(
            geographic_tile_tiff(tile_array(), lon0=self.LON, lat0=self.LAT, epsg=4266).getvalue()
        )
        ring = [(self.LON + c * ANADEM_STEP, self.LAT - r * ANADEM_STEP) for r, c in DOMAIN_RC]
        return dem, geojson(tmp_path / "d.geojson", ring, crs="EPSG:4266")

    @staticmethod
    def printed_line(output: str) -> str:
        lines = [
            line.strip()
            for line in ANSI.sub("", output).splitlines()
            if line.strip().startswith("--out-crs")
        ]
        assert len(lines) == 1, f"one line to copy, beginning --out-crs, in {output!r}"
        return lines[0]

    def refusal(self, tmp_path: Path, scene: tuple[Path, Path]) -> str:
        dem, domain = scene
        result = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb", "COLUMNS": "80"}).invoke(
            app,
            [
                *("mesh", "--dem", str(dem), "--domain", str(domain), "--tolerance", "1"),
                *("--out", str(tmp_path / "refused.vtk")),
            ],
        )
        assert result.exit_code == USAGE, result.output
        return self.printed_line(result.output)

    def mesh(self, tmp_path: Path, scene: tuple[Path, Path], out_crs: str) -> None:
        dem, domain = scene
        vtk = mesh_to_vtk(
            tmp_path,
            *("--dem", str(dem), "--domain", str(domain), "--tolerance", "1"),
            *("--out-crs", out_crs),
            out="pasted.vtk",
        )
        assert CRS.from_user_input(field(vtk, "crs")) == CRS.from_user_input(out_crs)

    def test_the_datum_does_hold_an_apostrophe(self) -> None:
        """The case is real: pyproj's own name for the datum."""
        datum = CRS.from_user_input("EPSG:4266").datum
        assert datum is not None and "'" in datum.name, datum

    def test_shlex_gives_back_the_crs_and_it_meshes(
        self, tmp_path: Path, scene: tuple[Path, Path]
    ) -> None:
        words = shlex.split(self.refusal(tmp_path, scene))
        assert len(words) == 2 and words[0] == "--out-crs", words
        datum = CRS.from_user_input(words[1]).datum
        assert datum is not None and "'" in datum.name, datum
        self.mesh(tmp_path, scene, words[1])

    def test_bash_gives_back_the_same_words(self, tmp_path: Path, scene: tuple[Path, Path]) -> None:
        line = self.refusal(tmp_path, scene)
        done = subprocess.run(
            ["bash", "-c", f"printf '%s\\0' {line}"],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        assert done.returncode == 0, done.stderr
        words = done.stdout.split("\0")[:-1]
        assert len(words) == 2 and words[0] == "--out-crs", words
        self.mesh(tmp_path, scene, words[1])
