"""`rasputin mesh --tolerance` with the edge strip: increment 15f-3, PY3-PY5.

`docs/increments/15f-edge-strip.md`, D6, D7, "The guarantee, and its wording"
(E1-E3, E7), L8 and "Tests for @tester" (PY3, PY4, PY5).

- **PY3, the projected path.** A synthetic GeoTIFF (non-square 10 x 5 m
  cells, rough relief), a domain with off-node corners and one polygon
  feature. The strip oracle (`strip_oracle.py`) is built from the **input
  polygons**: their crossings with the DEM's grid lines, z from the DEM array,
  measured against the written constraint lines. The midpoints between
  neighbouring crossings are ruled along the strip run's start edges (an
  edge's ends count as neighbours, default 2), which are `refine`'s output
  constraint edges, captured from the CLI's own `refine` call. Both sets are
  within the tolerance. The same oracle on `refine`'s output alone, in the
  same run, finds points over: the control. L8's other two oracles, E2 at the
  DEM nodes and constrained Delaunay, run on the written mesh.
- **PY4, the reprojected path.** 15c-2's geographic fixture (G6's scene) with
  `--out-crs`. The final check's independent source-node check is still 0
  over (E3), and the outline's crossings with the target grid's lines, against
  the resampled grid (rebuilt here through `open_dem`'s target tile, which is
  what `target_grid` produces), are within the tolerance (E1), with the same
  control on `refine`'s output.
- **PY5, `--stats`.** `edge strip: generate` on both paths; on the projected
  path also `edge strip: scan (parallel)` and `edge strip: split + flip
  (serial)`; on the reprojected path the joint loop's times stay in the
  `final check:` rows, so the strip's own scan and split rows are absent.

What the run reports is read as increment 25 names it
(`docs/increments/25-plain-output.md`, D2 and D6): the mesh file carries
`tolerance_m` and `max_error_m`, an "at most" bound on both paths; the
strip's counts are `--stats` Result rows (`line_points_checked`,
`line_max_error_m`, `line_points_on_nodata`, `line_points_refused` and its
`_max_error_m`, `line_points_inserted`, `line_check_dem_nodes_inserted`,
`line_points_duplicate`), and the DEM nodes within rounding of a vertex they
are not (15f's L14) are `dem_nodes_at_vertices` and its `_max_error_m`. 25's
D2 answers 15f's open question on L14's rounding exception: the file states
it through `max_error_m`, which includes those nodes, with no field of its
own.

CHOSEN HERE, where the design is silent: the NoData case puts the void on a
node outside the domain whose cells the outline crosses (a DEM edge void), so
`line_points_on_nodata` is at least 1 and the run still meets E1 at every
point it keeps. `dem_nodes_at_vertices` is a `--stats` row on the projected
path too, because `refine_strip` measures it there (25's D6: "where a path
produces them, they are always in `--stats`").

HOW THIS FILE GOES RED: the CLI has no strip yet, so the `line_*` rows and
the phase rows are absent from `--stats`, and the oracle finds the strip over
the tolerance on both paths. The controls pass already.

Not invariant-critical; no mutation round.
"""

from __future__ import annotations

import importlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from shapely.geometry import Polygon
from typer.testing import CliRunner

import tin_engine.cli as cli
from feature_fixtures import Feat, write_geojson
from geographic_fixtures import geographic_tile_tiff, project_ring
from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff
from recordread import stats_names
from strip_oracle import (
    Grid,
    StripFindings,
    delaunay_violations,
    edge_segments,
    exact_vertices,
    node_findings,
    polygon_segments,
    ruled_points,
    strip_findings,
)
from test_cli_mesh import plain
from test_cli_mesh_domain import geojson
from test_cli_mesh_geographic import (
    DOMAIN_RC,
    TARGET,
    geographic_check,
    lonlat_ring,
    tile_array,
    triangles,
)
from test_cli_mesh_geographic import write_geojson as write_domain
from test_cli_mesh_refine import file_field, stats_row
from test_cli_mesh_stats import seconds
from tin_engine.cli import app
from tin_engine.dem_input import DemRequest, open_dem
from tin_engine.domain import read_domain
from vtkread import VtkFile, lines_as_array, read_vtk

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})
USAGE = 2
TOLERANCE = 1.0
#: The noder's snap (`DEFAULT_SNAP_SPACING`): how far a written constraint
#: line may lie from the input polygon, in metres, with margin.
ON_INPUT = 2e-3
#: How far a point ruled on the start's own edges may lie from a written line.
ON_START = 1e-6


@dataclass(frozen=True)
class LineCheck:
    """The strip's `--stats` rows (25's D6)."""

    checked: int
    on_nodata: int
    inserted: int
    max_error: float
    refused: int
    refused_max_error: float
    dem_nodes_inserted: int | None  # refine_strip's; absent on the reprojected path (P2)
    duplicate: int


def line_check(report: str) -> LineCheck:
    def count(name: str) -> int:
        return int(stats_row(report, name))

    def metres(name: str) -> float:
        return float(stats_row(report, name))

    return LineCheck(
        checked=count("line_points_checked"),
        on_nodata=count("line_points_on_nodata"),
        inserted=count("line_points_inserted"),
        max_error=metres("line_max_error_m"),
        refused=count("line_points_refused"),
        refused_max_error=metres("line_points_refused_max_error_m"),
        dem_nodes_inserted=(
            count("line_check_dem_nodes_inserted")
            if "line_check_dem_nodes_inserted" in stats_names(report)
            else None
        ),
        duplicate=count("line_points_duplicate"),
    )


@dataclass(frozen=True)
class Run:
    """One `rasputin mesh` run: the file, the `--stats` report, and what
    the CLI's `refine` call was given and returned."""

    vtk: VtkFile
    report: str
    refine_start: Any
    refined: Any

    def written(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self.vtk.points, self.vtk.points[:, 2], lines_as_array(self.vtk)

    def refined_arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r = self.refined
        return np.asarray(r.vertices), np.asarray(r.z), np.asarray(r.edges)


def mesh(tmp: Path, *args: str) -> Run:
    """Run the CLI with `--stats`, capturing `refine`'s call through the module
    attribute the CLI calls it by."""
    seen: list[tuple[Any, Any]] = []
    real = cli.refine

    def spy(*a: Any, **k: Any) -> Any:
        out = real(*a, **k)
        seen.append((a[1], out))
        return out

    out, md = tmp / "x.vtk", tmp / "x.md"
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(cli, "refine", spy)
        result = runner.invoke(app, ["mesh", *args, "--out", str(out), "--stats", str(md)])
    assert result.exit_code == 0, plain(result.output)
    ((start, refined),) = seen
    return Run(
        read_vtk(out.read_bytes()),
        md.read_text(encoding="utf-8"),
        start,
        refined,
    )


def measure(
    run_: Run, grid: Grid, rings: list[list[tuple[float, float]]], *, written: bool
) -> tuple[StripFindings, StripFindings]:
    """E1 by the oracle, on the written mesh or on refine's output alone:
    the input polygons' crossings, then the crossings and midpoints ruled on
    the strip run's start edges."""
    vertices, z, edges = run_.written() if written else run_.refined_arrays()
    xy, zp = ruled_points(polygon_segments(rings), grid, midpoints=False)
    by_input = strip_findings(
        xy,
        zp,
        vertices,
        z,
        edges,
        tolerance=TOLERANCE,
        on_edge=ON_INPUT,
        slack=ON_INPUT * grid.slope_bound(),
    )
    start_v, _, start_e = run_.refined_arrays()
    xy, zp = ruled_points(edge_segments(start_v, start_e), grid)
    by_start = strip_findings(
        xy, zp, vertices, z, edges, tolerance=TOLERANCE, on_edge=ON_START, slack=0.0
    )
    return by_input, by_start


# ------------------------------------------------------------------ PY3: the projected path

#: micro_tiff's frame: node (row, col) at (TIE_X + 10 col, TIE_Y - 5 row).
ROWS, COLS, DX, DY = 25, 31, 10.0, 5.0


def rel(x: float, y: float) -> tuple[float, float]:
    return (TIE_X + x, TIE_Y - y)


#: Every vertex off-node; counter-clockwise.
DOMAIN = [rel(12.3, 103.7), rel(287.9, 101.1), rel(281.3, 7.3), rel(17.7, 9.9)]
FEATURE = [rel(83.3, 71.9), rel(201.7, 74.3), rel(196.1, 33.1), rel(91.9, 37.7)]


def relief(seed: int = 1) -> np.ndarray:
    """Relief at about two cells' wavelength plus noise, float32."""
    r, c = np.indices((ROWS, COLS), dtype=np.float64)
    noise = np.random.default_rng(seed).normal(0.0, 1.5, (ROWS, COLS))
    z = 100 + 12 * np.sin(c / 2.3) * np.cos(r / 1.9) + 6 * np.sin((r + 2 * c) / 1.7) + noise
    return z.astype(np.float32)


def projected_grid(array: np.ndarray) -> Grid:
    return Grid(TIE_X, TIE_Y, DX, DY, array.astype(np.float64))


def projected_scene(tmp: Path, array: np.ndarray, *, nodata: str | None = None) -> list[str]:
    dem = tmp / "dem.tif"
    dem.write_bytes(micro_tiff(array, nodata=nodata).getvalue())
    domain = geojson(tmp / "domain.geojson", DOMAIN)
    features = write_geojson(
        tmp / "features.geojson", [Feat("lake", Polygon(FEATURE), {"property": "water"})]
    )
    return ["--dem", str(dem), "--domain", str(domain), "--features", str(features)]


@pytest.fixture(scope="module")
def projected(tmp_path_factory: pytest.TempPathFactory) -> Run:
    tmp = tmp_path_factory.mktemp("projected")
    return mesh(tmp, *projected_scene(tmp, relief()), "--tolerance", str(TOLERANCE))


class TestProjectedPath:
    """PY3."""

    def test_the_stats_carry_the_line_check(self, projected: Run) -> None:
        found = line_check(projected.report)
        assert found.checked > 0 and found.inserted > 0, found
        assert found.on_nodata == 0
        assert 0.0 <= found.max_error <= TOLERANCE
        assert (found.refused, found.refused_max_error) == (0, 0.0)
        assert found.dem_nodes_inserted is not None and found.dem_nodes_inserted >= 0
        assert found.duplicate >= 0
        assert float(file_field(projected.vtk, "max_error_m")) <= TOLERANCE

    def test_the_nodes_at_vertices_are_measured_on_the_projected_path(self, projected: Run) -> None:
        """25's D2 and D6: `refine_strip` counts the DEM nodes within rounding
        of a vertex (L14), so the rows are present, and `max_error_m` is at
        least their difference."""
        report = projected.report
        assert int(stats_row(report, "dem_nodes_at_vertices")) >= 0
        at_vertices = float(stats_row(report, "dem_nodes_at_vertices_max_error_m"))
        assert float(file_field(projected.vtk, "max_error_m")) >= at_vertices

    def test_e1_every_crossing_and_midpoint_is_within_tolerance(self, projected: Run) -> None:
        rings = [DOMAIN, FEATURE]
        by_input, by_start = measure(projected, projected_grid(relief()), rings, written=True)
        assert by_input.points > 50 and by_start.points > by_input.points
        assert (by_input.unlocated, by_input.over) == (0, 0), by_input
        assert (by_start.unlocated, by_start.over) == (0, 0), by_start

    def test_refine_alone_leaves_crossings_over(self, projected: Run) -> None:
        """The control, from the same run: refine's output, before the strip."""
        rings = [DOMAIN, FEATURE]
        by_input, by_start = measure(projected, projected_grid(relief()), rings, written=False)
        assert by_input.unlocated == 0 and by_input.over > 0, by_input
        assert by_start.unlocated == 0 and by_start.over > 0, by_start

    def test_e2_and_delaunay_on_the_written_mesh(self, projected: Run) -> None:
        grid = projected_grid(relief())
        v, z, edges = projected.written()
        tri = triangles(projected.vtk)
        nodes = node_findings(grid, v, z, np.ones(len(v), bool), tri, tolerance=TOLERANCE)
        assert nodes.nodes > 300 and nodes.over == 0, nodes
        exact = exact_vertices(grid, v, np.asarray(projected.refined.vertices))
        assert delaunay_violations(grid, v, tri, edges, exact) == 0

    def test_e7_an_outline_along_grid_lines_is_left_as_refine_wrote_it(
        self, tmp_path: Path
    ) -> None:
        """Without `--domain` the start's outline runs along the DEM's border,
        a vertex every fourth node: every crossing is a node and every
        midpoint is on a cell side,
        all within the tolerance after refine, so the strip run changes
        nothing (E7 and the grid-line remark)."""
        dem = tmp_path / "dem.tif"
        dem.write_bytes(micro_tiff(relief()).getvalue())
        run_ = mesh(tmp_path, "--dem", str(dem), "--stride", "4", "--tolerance", str(TOLERANCE))
        found = line_check(run_.report)
        assert found.checked > 0 and found.inserted == 0, found
        assert found.dem_nodes_inserted == 0, found
        v, z, _ = run_.refined_arrays()
        np.testing.assert_array_equal(run_.vtk.points[:, :2], v)
        np.testing.assert_array_equal(run_.vtk.points[:, 2], z)
        np.testing.assert_array_equal(triangles(run_.vtk), np.asarray(run_.refined.triangles))

    def test_points_beside_nodata_are_counted_and_the_rest_hold(self, tmp_path: Path) -> None:
        """§3A, non-finite input: a void on node (row 21, col 15), outside the
        domain, a corner of the cells the outline's south side crosses."""
        array = relief()
        array[21, 15] = -32767
        run_ = mesh(
            tmp_path,
            *projected_scene(tmp_path, array, nodata="-32767"),
            "--tolerance",
            str(TOLERANCE),
        )
        found = line_check(run_.report)
        assert found.on_nodata >= 1, found
        void = array.astype(np.float64)
        void[21, 15] = np.nan
        by_input, by_start = measure(run_, projected_grid(void), [DOMAIN, FEATURE], written=True)
        assert by_input.over == 0, by_input
        assert by_start.over == 0, by_start

    def test_a_failed_strip_run_writes_nothing(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """D6: a refusal from the run is a usage error in the engine's words."""
        edge_strip = importlib.import_module("tin_engine.edge_strip")
        real = edge_strip.run

        class Failed:
            message = "planted strip-run failure"

            def __init__(self, outcome: Any) -> None:
                self._outcome = outcome

            def ok(self) -> bool:
                return False

            def __getattr__(self, name: str) -> Any:
                return getattr(self._outcome, name)

        monkeypatch.setattr(edge_strip, "run", lambda *a, **k: Failed(real(*a, **k)))
        out = tmp_path / "refused.vtk"
        args = projected_scene(tmp_path, relief())
        result = runner.invoke(
            app, ["mesh", *args, "--tolerance", str(TOLERANCE), "--out", str(out)]
        )
        assert result.exit_code == USAGE, result.output
        needle = "planted strip-run failure"
        assert "".join(needle.split()) in "".join(plain(result.output).split())
        assert "Traceback" not in result.output
        assert not out.exists()


# ------------------------------------------------------------------ PY4: the reprojected path


def target_domain_ring(domain: Path) -> list[tuple[float, float]]:
    given = json.loads(domain.read_text())
    return project_ring(given["crs"]["properties"]["name"], TARGET, given["coordinates"][0][:-1])


@dataclass(frozen=True)
class Reprojected:
    run: Run
    domain: Path
    grid: Grid


@pytest.fixture(scope="module")
def reprojected(tmp_path_factory: pytest.TempPathFactory) -> Reprojected:
    tmp = tmp_path_factory.mktemp("reprojected")
    dem = tmp / "anadem_like.tif"
    dem.write_bytes(geographic_tile_tiff(tile_array()).getvalue())
    ring = project_ring("EPSG:4326", "EPSG:4674", lonlat_ring(DOMAIN_RC))
    domain = write_domain(tmp / "domain.geojson", ring, "EPSG:4674")
    run_ = mesh(
        tmp,
        *("--dem", str(dem), "--domain", str(domain)),
        *("--out-crs", TARGET, "--tolerance", str(TOLERANCE)),
    )
    opened = open_dem(DemRequest(sources=(dem,), domain=read_domain(domain), target_crs=TARGET))
    m = opened.tile.meta
    array = np.asarray(opened.tile.array, np.float64).copy()
    if m.nodata is not None:
        array[array == m.nodata] = np.nan
    return Reprojected(run_, domain, Grid(m.x_min, m.y_max, m.delta_x, m.delta_y, array))


class TestReprojectedPath:
    """PY4."""

    def test_the_stats_carry_the_line_check_beside_the_final_checks(
        self, reprojected: Reprojected
    ) -> None:
        report, vtk = reprojected.run.report, reprojected.run.vtk
        found = line_check(report)
        assert found.checked > 0 and found.inserted > 0, found
        assert 0.0 <= found.max_error <= TOLERANCE
        assert found.refused == 0
        assert found.dem_nodes_inserted is None, "refine_strip's count: absent on this path"
        assert stats_row(report, "resampled_grid").startswith("30 m square grid in")
        assert int(stats_row(report, "dem_nodes_checked")) > 1000
        stated = float(file_field(vtk, "max_error_m"))
        at_vertices = float(stats_row(report, "dem_nodes_at_vertices_max_error_m"))
        assert at_vertices <= stated <= max(TOLERANCE, at_vertices)

    def test_e1_the_outlines_crossings_against_the_resampled_grid(
        self, reprojected: Reprojected
    ) -> None:
        ring = target_domain_ring(reprojected.domain)
        by_input, by_start = measure(reprojected.run, reprojected.grid, [ring], written=True)
        assert by_input.points > 50
        assert (by_input.unlocated, by_input.over) == (0, 0), by_input
        assert (by_start.unlocated, by_start.over) == (0, 0), by_start

    def test_refine_alone_leaves_crossings_over(self, reprojected: Reprojected) -> None:
        ring = target_domain_ring(reprojected.domain)
        by_input, by_start = measure(reprojected.run, reprojected.grid, [ring], written=False)
        assert by_input.over > 0, by_input
        assert by_start.over > 0, by_start

    def test_e3_the_source_nodes_and_delaunay_still_hold(self, reprojected: Reprojected) -> None:
        vtk = reprojected.run.vtk
        found = geographic_check(vtk, reprojected.domain)
        assert found.nodes > 1000
        assert (found.over_interior, found.over_strip) == (0, 0), found
        exact = exact_vertices(
            reprojected.grid, vtk.points, np.asarray(reprojected.run.refined.vertices)
        )
        tri, edges = triangles(vtk), lines_as_array(vtk)
        assert delaunay_violations(reprojected.grid, vtk.points, tri, edges, exact) == 0


# ------------------------------------------------------------------ PY5: --stats


STRIP_ROWS = ("edge strip: scan (parallel)", "edge strip: split + flip (serial)")


class TestStatsRows:
    """PY5."""

    def test_on_the_projected_path(self, projected: Run) -> None:
        phases = seconds(projected.report)
        for name in ("edge strip: generate", *STRIP_ROWS):
            assert name in phases, (name, sorted(phases))
            assert phases[name] >= 0.0

    def test_on_the_reprojected_path(self, reprojected: Reprojected) -> None:
        phases = seconds(reprojected.run.report)
        assert "edge strip: generate" in phases, sorted(phases)
        for name in ("final check: scan (parallel)", "final check: split + flip (serial)"):
            assert name in phases, (name, sorted(phases))
        for name in STRIP_ROWS:
            assert name not in phases, f"{name}: the joint loop's time is the final check's"
