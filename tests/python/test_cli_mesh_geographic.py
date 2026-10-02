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
fields are read as: `crs` the target's label, `source_crs` the DEM's,
`source_transform` pyproj's `description` of the source-to-target
`always_xy` transformer, `computation_grid` beginning `square <h> m grid in
<crs>` and naming `resampled bilinear from <source crs>`, and
`elevation_source` holding `against the resampled grid` and `checked against
<N> source nodes`.

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

import importlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from pyproj import Transformer
from shapely.geometry import Polygon

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
from test_cli_catchment import invoke as invoke_any
from test_cli_mesh_domain import quarter_circle
from test_cli_mesh_mosaic import USAGE, invoke, same_mesh
from tin_engine.dem_input import DemRequest, open_dem
from tin_engine.domain import read_domain
from tin_engine.io.repository import TiffDemRepository
from vtkread import VtkFile, read_vtk

Ring = list[tuple[float, float]]
TARGET = "EPSG:31983"
ROWS = COLS = 60
TOLERANCE = 1.0
H = 30  # ANADEM's 29.8 m at 19 S, rounded (Q16)
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
VELHAS = FIXTURES / "velhas"  # Q17: the ANADEM extract and the DEM-derived catchment


def squashed(text: str) -> str:
    """No whitespace: Rich may break a long message anywhere inside its panel."""
    return "".join(text.split())


def write_geojson(path: Path, ring: Ring, crs: str) -> Path:
    doc = {
        "type": "Polygon",
        "coordinates": [[*ring, ring[0]]],
        "crs": {"type": "name", "properties": {"name": crs}},
    }
    path.write_text(json.dumps(doc))
    return path


def field(vtk: VtkFile, name: str) -> str:
    assert name in vtk.field_data, sorted(vtk.field_data)
    (value,) = vtk.field_data[name].values
    return str(value)


def run(tmp_path: Path, *args: str, out: str = "x.vtk") -> VtkFile:
    target = tmp_path / out
    code, output = invoke(*args, "--out", str(target))
    assert code == 0, output
    return read_vtk(target.read_bytes())


def refused(tmp_path: Path, *args: str, says: tuple[str, ...]) -> str:
    target = tmp_path / "refused.vtk"
    code, output = invoke(*args, "--out", str(target))
    assert code == USAGE, output
    assert "Traceback" not in output
    assert "No such option" not in output, "refused for the wrong reason"
    for word in says:
        assert squashed(word) in squashed(output), f"{word!r} not in {output!r}"
    assert not target.exists()
    return output


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
    return write_geojson(tmp_path / "domain.geojson", ring, "EPSG:4674")


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
        vtk = run(
            tmp_path,
            *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
            *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
        )
        assert field(vtk, "crs") == TARGET
        assert field(vtk, "source_crs") == "EPSG:4326"
        expected = Transformer.from_crs("EPSG:4326", TARGET, always_xy=True).description
        assert field(vtk, "source_transform") == expected
        assert field(vtk, "domain_crs") == "EPSG:4674"
        grid = field(vtk, "computation_grid")
        assert grid.startswith(f"square {H} m grid in {TARGET}"), grid
        assert "resampled bilinear from EPSG:4326" in grid, grid
        sentence = field(vtk, "elevation_source")
        assert "against the resampled grid" in sentence, sentence
        checked = re.search(r"checked against (\d+) source nodes", sentence)
        assert checked is not None, sentence
        found = geographic_check(vtk, domain_4674)
        assert found.nodes <= int(checked.group(1)) <= ROWS * COLS - 3
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
        vtk = run(
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
        moved = write_geojson(
            tmp_path / "domain_31983.geojson",
            project_ring("EPSG:4674", TARGET, given["coordinates"][0][:-1]),
            TARGET,
        )
        vtk = run(
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
        vtk = run(
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
        domain = write_geojson(tmp_path / "d.geojson", ring, "EPSG:25833")
        vtk = run(
            tmp_path,
            *("--dem", str(dem), "--domain", str(domain)),
            *("--out-crs", "EPSG:25832", "--tolerance", str(TOLERANCE)),
        )
        assert field(vtk, "crs") == "EPSG:25832"
        assert field(vtk, "source_crs") == "EPSG:25833"
        grid = field(vtk, "computation_grid")
        assert grid.startswith("square 10 m grid in EPSG:25832"), grid
        assert "checked against" in field(vtk, "elevation_source")
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
    domain = write_geojson(tmp_path / "quarter.geojson", quarter_circle(), "EPSG:25833")
    args = ("--dem", str(KARTVERKET), "--domain", str(domain), "--tolerance", "10")
    plain_code, plain_output = invoke(*args, "--out", str(tmp_path / "plain.vtk"))
    assert plain_code == 0, plain_output
    code, output = invoke(*args, "--out-crs", "EPSG:25833", "--out", str(tmp_path / "same.vtk"))
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
            *("--dem", str(geographic_dem), "--domain", str(domain_4674), "--tolerance", "1"),
            says=("--out-crs", suggestion.proj, suggestion.family, "worst scale error"),
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
        domain = write_geojson(tmp_path / "d.geojson", ring, "EPSG:4326")
        output = refused(
            tmp_path,
            *("--dem", str(dem), "--domain", str(domain), "--tolerance", "1"),
            *("--out-crs", "+proj=tmerc +lat_0=0 +lon_0=180 +k=1 +ellps=WGS84 +units=m"),
            says=("180",),
        )
        assert re.search(r"(?i)antimeridian|longitude", output), output

    def test_a_domain_vertex_without_an_image(
        self, tmp_path: Path, geographic_dem: Path, no_pixels: None
    ) -> None:
        """A vertex on the far side of the globe from an orthographic target."""
        ring = [(-44.01, -19.01), (-43.99, -19.01), (136.0, 19.0)]
        domain = write_geojson(tmp_path / "d.geojson", ring, "EPSG:4326")
        refused(
            tmp_path,
            *("--dem", str(geographic_dem), "--domain", str(domain), "--tolerance", "1"),
            *("--out-crs", "+proj=ortho +lat_0=-19 +lon_0=-44 +ellps=WGS84 +units=m"),
            says=("no image",),
        )


# ------------------------------------------------------------------ G2


def test_catchment_on_a_geographic_tile_is_a_usage_error(
    tmp_path: Path, geographic_dem: Path
) -> None:
    """The `to_core` gate keeps a geographic DEM out of `rasputin catchment`
    (increment 22), which 15c does not extend."""
    lon, lat = LON0 + 30 * ANADEM_STEP, LAT0 - 30 * ANADEM_STEP
    code, output = invoke_any(
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
        by_key = run(
            tmp_path, "--dem", KEY, "--cache", str(tmp_path / "cache"), *common, out="k.vtk"
        )
        by_path = run(tmp_path, "--dem", str(path), *common, out="p.vtk")
        same_mesh(by_key, by_path)

    def test_not_cached_names_the_runs_out_crs(
        self, tmp_path: Path, catalogue: None, domain_4674: Path
    ) -> None:
        write_cache(
            tmp_path / "cache", KEY, {"dem": cog_bytes()}, skip={"dem": [5]}, crs="EPSG:4326"
        )
        refused(
            tmp_path,
            *("--dem", KEY, "--cache", str(tmp_path / "cache"), "--domain", str(domain_4674)),
            *("--out-crs", TARGET, "--tolerance", "1"),
            says=(f"rasputin fetch {KEY}", f"--out-crs {TARGET}"),
        )


# ------------------------------------------------------------------ realism (Q17)


def test_the_anadem_extract_over_its_catchment(tmp_path: Path) -> None:
    """G6's realism case: the committed ANADEM extract (EPSG:4674, the COG's
    key set) over the DEM-derived catchment (EPSG:31983), checked against
    every source node inside it."""
    import tifffile

    dem, catchment = VELHAS / "anadem_velhas.tif", VELHAS / "catchment.geojson"
    vtk = run(
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
            *("--dem", str(geographic_dem), "--domain", str(domain_4674), "--tolerance", "1"),
            *("--out-crs", out_crs),
            says=("--out-crs",),
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
        vtk = run(
            tmp_path,
            *("--dem", str(geographic_dem), "--domain", str(domain_4674)),
            *("--out-crs", proj, "--tolerance", str(TOLERANCE)),
        )
        assert "checked against" in field(vtk, "elevation_source")


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
        vtk = run(tmp_path, *self.args(geographic_dem, domain_4674))
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
        """`checked against N source nodes` is the store's size, which equals
        the points `check_point_blocks` yields for the same request."""
        monkeypatch.setattr("tin_engine.target_grid.BLOCK", 16)
        opened = open_dem(
            DemRequest(
                sources=(geographic_dem,), domain=read_domain(domain_4674), target_crs=TARGET
            )
        )
        assert opened.checks is not None
        yielded = sum(len(z) for _, z in opened.checks)
        vtk = run(tmp_path, *self.args(geographic_dem, domain_4674))
        checked = re.search(r"checked against (\d+) source nodes", field(vtk, "elevation_source"))
        assert checked is not None
        assert int(checked.group(1)) == yielded > 1000

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
                *("--dem", str(geographic_dem), "--domain", str(domain_4674), "--tolerance", "1"),
                says=("worst scale error",),
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
            *self.args(geographic_dem, domain_4674),
            says=("planted final-check failure",),
        )
