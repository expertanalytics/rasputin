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

What the sentence pins is D7's, by `re.search` on each part, never the whole
sentence: the strip clause as D7 writes it, `max error at DEM nodes at most
<e> m` in place of `achieved max error` on the projected path only, the
clause after the final check's on the reprojected path, and the stderr line.
**Ola has not ruled whether the sentence states L14's rounding exception**
(E2's DEM nodes within `r(g)` of a vertex; `ASK OLA` under L14). Nothing here
pins its absence. If Ola says yes, the one test to extend is
`TestProjectedPath.test_the_sentence_carries_the_strip_clause`, with the new
clause's wording; the other searches stay as they are, unless the clause is
put inside the strip clause or between `at most` and its number.

CHOSEN HERE, where the design is silent: the NoData case puts the void on a
node outside the domain whose cells the outline crosses (a DEM edge void), so
`no_data` is at least 1 and the run still meets E1 at every point it keeps.

HOW THIS FILE GOES RED: the CLI has no strip yet, so the clause, the stderr
line, the `--stats` rows and the new wording are absent, and the oracle finds
the strip over the tolerance on both paths. The controls pass already.

Not invariant-critical; no mutation round.
"""

from __future__ import annotations

import importlib
import json
import re
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
from test_cli_mesh_refine import NUMBER, field, sentence
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

CLAUSE = re.compile(
    r"; edge strip: (\d+) check points where constraints cross grid lines and between them "
    r"\((\d+) without data\), (\d+) inserted, max error " + NUMBER + r" m at them, "
    r"(\d+) refused \(max " + NUMBER + r" m\)"
)
REPORT = re.compile(r"(\d+) strip points inserted, (\d+) nodes inserted by the strip run")


@dataclass(frozen=True)
class Clause:
    points: int
    no_data: int
    inserted: int
    max_error: float
    refused: int
    refused_max_error: float


def clause(text: str) -> Clause:
    match = CLAUSE.search(text)
    assert match is not None, f"no edge-strip clause in {text!r}"
    p, n, i, e, r, re_ = match.groups()
    return Clause(int(p), int(n), int(i), float(e), int(r), float(re_))


@dataclass(frozen=True)
class Run:
    """One `rasputin mesh` run: the file, the streams, the report, and what
    the CLI's `refine` call was given and returned."""

    vtk: VtkFile
    stderr: str
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
        result.stderr,
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

    def test_the_sentence_carries_the_strip_clause(self, projected: Run) -> None:
        text = sentence(projected.vtk)
        found = clause(text)
        assert found.points > 0 and found.inserted > 0, found
        assert found.no_data == 0
        assert 0.0 <= found.max_error <= TOLERANCE
        assert (found.refused, found.refused_max_error) == (0, 0.0)
        assert field(text, rf"max error at DEM nodes at most {NUMBER} m") <= TOLERANCE
        assert "achieved max error" not in text, text

    def test_stderr_reports_the_strip_runs_insertions(self, projected: Run) -> None:
        match = REPORT.search(" ".join(projected.stderr.split()))
        assert match is not None, projected.stderr
        assert int(match.group(1)) == clause(sentence(projected.vtk)).inserted

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
        found = clause(sentence(run_.vtk))
        assert found.points > 0 and found.inserted == 0, found
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
        found = clause(sentence(run_.vtk))
        assert found.no_data >= 1, found
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
        assert "planted strip-run failure" in "".join(plain(result.output).split())
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

    def test_the_sentence_carries_the_clause_after_the_final_checks(
        self, reprojected: Reprojected
    ) -> None:
        text = sentence(reprojected.run.vtk)
        found = clause(text)
        assert found.points > 0 and found.inserted > 0, found
        assert 0.0 <= found.max_error <= TOLERANCE
        assert found.refused == 0
        assert "against the resampled grid" in text
        assert re.search(r"checked against \d+ source nodes", text), text
        assert text.index("checked against") < text.index("; edge strip:")
        assert "achieved max error" in text, "the projected-path wording leaked across"

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
