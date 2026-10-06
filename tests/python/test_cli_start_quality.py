"""`rasputin mesh --start-min-angle DEG`: increment 20, R9 and R11.

`docs/increments/20-start-quality.md`. Built on the main session's provisional
picks, pending Ola: C1 (a), C2 (a) 25°, C3 (a), C4 (a).

Since increment 25 (``docs/increments/25-plain-output.md``, D4) the setting
and the pass's counts are ``--stats`` rows: ``start_min_angle_deg`` (``0`` is
off), ``start_quality_points_inserted`` and ``start_quality_points_skipped``,
the outcome's ``quality_inserted`` and ``quality_skipped``; stderr no longer
carries them (D7). The timing row ``refine: start quality`` stays.

The binding keyword is the design's ``min_angle_deg``; the CLI passes it on
every refined run, 25.0 by default, so the tests read it off a wrapped
``cli.refine``.

T-real is logged, not thresholded, except the tolerance and one direction:
the quarter circle's triangles under 1° at 10 m drop against
``--start-min-angle 0``. Not invariant-critical, no mutation round.
"""

from __future__ import annotations

import io
import re
import time
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import pytest
import shapely
from shapely.geometry import Polygon

import tin_engine.cli as cli
from cli_helpers import COLS, ROWS, SQUARE, USAGE, geojson, invoke, rough_dem, write_tiff
from geotiff_fixtures import KARTVERKET, micro_tiff, needs_codecs
from test_cli_mesh_dem import SENTINEL
from test_cli_mesh_domain import SNAP, quarter_circle
from test_cli_mesh_features import _delaunay_violations
from test_cli_mesh_plain_output import located_errors
from test_cli_mesh_refine import file_field, min_angles_degrees, stats_row
from test_cli_mesh_stats import mesh, seconds
from tin_engine import _core
from tin_engine.io.geotiff import decode_dem
from vtkread import VtkFile, lines_as_array, polygons_as_array, read_vtk

QUALITY = re.compile(r"(\d+) start quality nodes inserted, (\d+) start quality skips")  # pre-25
COUNTS = ("start_quality_points_inserted", "start_quality_points_skipped")


bumpy = rough_dem(20)


@pytest.fixture
def box(tmp_path: Path) -> Path:
    """SQUARE without its hole: 175 m x 66 m, so its two start triangles have
    a 20.7° angle and the pass has work to do."""
    return geojson(tmp_path / "box.geojson", SQUARE)


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    """Every keyword set ``cli.refine`` is called with, passed through."""
    seen: list[dict[str, Any]] = []
    real = cli.refine

    def spy(*args: Any, **kwargs: Any) -> Any:
        seen.append(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(cli, "refine", spy)
    return seen


def run(tmp_path: Path, *args: str, name: str = "x.vtk") -> tuple[VtkFile, str, str]:
    """(mesh, ``--stats`` report, stderr)."""
    out, md = tmp_path / name, tmp_path / f"{name}.md"
    result = mesh(*args, "--out", str(out), "--stats", str(md))
    return read_vtk(out.read_bytes()), md.read_text(encoding="utf-8"), result.stderr


class TestTheFlag:
    """R11: default 25, ``0`` is off, forwarded on both start paths."""

    def test_the_default_is_25_on_a_domain_start(
        self, tmp_path: Path, bumpy: Path, box: Path, calls: list[dict[str, Any]]
    ) -> None:
        _, report, _ = run(tmp_path, "--dem", str(bumpy), "--domain", str(box), "--tolerance", "1")
        assert [c["min_angle_deg"] for c in calls] == [25.0]
        assert stats_row(report, "start_min_angle_deg") == "25"

    def test_the_default_is_25_on_a_stride_start(
        self, tmp_path: Path, bumpy: Path, calls: list[dict[str, Any]]
    ) -> None:
        _, report, _ = run(tmp_path, "--dem", str(bumpy), "--tolerance", "1")
        assert [c["min_angle_deg"] for c in calls] == [25.0]
        assert stats_row(report, "start_min_angle_deg") == "25"

    def test_zero_is_off(
        self, tmp_path: Path, bumpy: Path, box: Path, calls: list[dict[str, Any]]
    ) -> None:
        vtk, report, stderr = run(
            tmp_path,
            "--dem",
            str(bumpy),
            "--domain",
            str(box),
            "--tolerance",
            "1",
            "--start-min-angle",
            "0",
        )
        assert [c["min_angle_deg"] for c in calls] == [0.0]
        assert stats_row(report, "start_min_angle_deg") == "0"
        assert tuple(stats_row(report, n) for n in COUNTS) == ("0", "0")
        assert QUALITY.search(stderr) is None, stderr  # D7: the counts left stderr
        assert "start_min_angle_deg" not in vtk.field_data  # D2: not a mesh-file field

    @pytest.mark.parametrize("value", ["30", "35"])
    def test_a_value_up_to_35_is_accepted_and_named(
        self, tmp_path: Path, bumpy: Path, box: Path, value: str, calls: list[dict[str, Any]]
    ) -> None:
        _, report, _ = run(
            tmp_path,
            "--dem",
            str(bumpy),
            "--domain",
            str(box),
            "--tolerance",
            "1",
            "--start-min-angle",
            value,
        )
        assert [c["min_angle_deg"] for c in calls] == [float(value)]
        assert stats_row(report, "start_min_angle_deg") == value


class TestRefusals:
    """R11: negative, non-finite, above 35, and without ``--tolerance``."""

    @pytest.mark.parametrize("value", ["-1", "-0.5", "nan", "inf", "36", "35.5"])
    def test_a_bad_value(self, tmp_path: Path, bumpy: Path, box: Path, value: str) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "mesh",
            "--dem",
            str(bumpy),
            "--domain",
            str(box),
            "--tolerance",
            "1",
            "--start-min-angle",
            value,
            "--out",
            str(target),
        )
        assert code == USAGE, output
        assert "No such option" not in output  # refused for its value, not unknown
        assert "--start-min-angle" in output
        assert not target.exists()

    def test_without_tolerance(self, tmp_path: Path, bumpy: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "mesh", "--dem", str(bumpy), "--start-min-angle", "25", "--out", str(target)
        )
        assert code == USAGE, output
        assert "No such option" not in output  # refused for its value, not unknown
        assert "--start-min-angle" in output
        assert "--tolerance" in output
        assert not target.exists()

    def test_zero_without_tolerance_is_refused_too(self, tmp_path: Path, bumpy: Path) -> None:
        """Given at all, it needs --tolerance; 0 is a value, not an absence."""
        code, output = invoke(
            "mesh", "--dem", str(bumpy), "--start-min-angle", "0", "--out", str(tmp_path / "x.vtk")
        )
        assert code == USAGE, output
        assert "No such option" not in output  # refused for its value, not unknown
        assert "--start-min-angle" in output

    def test_without_dem(self, tmp_path: Path) -> None:
        code, output = invoke(
            "mesh",
            "catchment",
            "--flat",
            "--start-min-angle",
            "25",
            "--out",
            str(tmp_path / "x.vtk"),
        )
        assert code == USAGE, output
        assert "No such option" not in output  # refused for its value, not unknown
        assert "--start-min-angle" in output


class TestReport:
    """R11, as increment 25 words it: ``--stats`` carries the two counts and a
    timing row; stderr carries neither (D7)."""

    def test_stats_counts_the_pass(self, tmp_path: Path, bumpy: Path, box: Path) -> None:
        _, report, stderr = run(
            tmp_path, "--dem", str(bumpy), "--domain", str(box), "--tolerance", "1"
        )
        inserted = int(stats_row(report, "start_quality_points_inserted"))
        assert inserted > 0  # the box's 20.7° start triangles are bad at 25°
        assert int(stats_row(report, "refinement_rounds")) >= 0
        assert QUALITY.search(stderr) is None, stderr
        assert not re.search(r"(\d+) rounds, (\d+) points inserted, (\d+) flips", stderr)

    def test_stats_has_the_outcome_counts_and_the_timing_row(
        self, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        seen: list[_core.RefineOutcome] = []
        real = cli.refine

        def spy(*args: Any, **kwargs: Any) -> _core.RefineOutcome:
            out = real(*args, **kwargs)
            seen.append(out)
            return out

        monkeypatch.setattr(cli, "refine", spy)
        _, report, _ = run(tmp_path, "--dem", str(bumpy), "--domain", str(box), "--tolerance", "1")
        (out,) = seen
        expected = (str(out.quality_inserted), str(out.quality_skipped))
        assert tuple(stats_row(report, n) for n in COUNTS) == expected
        phases = seconds(report)
        assert "refine: start quality" in phases
        assert 0.0 <= phases["refine: start quality"] <= phases["refine"]

    def test_stats_with_the_pass_off_still_has_the_row(
        self, tmp_path: Path, bumpy: Path, box: Path
    ) -> None:
        md = tmp_path / "x.md"
        mesh(
            "--dem",
            str(bumpy),
            "--domain",
            str(box),
            "--tolerance",
            "1",
            "--start-min-angle",
            "0",
            "--out",
            str(tmp_path / "x.vtk"),
            "--stats",
            str(md),
        )
        report = md.read_text(encoding="utf-8")
        assert tuple(stats_row(report, n) for n in COUNTS) == ("0", "0")
        assert "refine: start quality" in seconds(report)


class TestNoData:
    """Q10, as the fix replaces R9 (``20-start-quality.md``, "Fix: the
    start-quality pass skips NoData nodes"): the pass's node lies in the
    block, so it is skipped and counted, and the trim has nothing to drop."""

    def test_a_domain_over_a_nodata_block(self, tmp_path: Path, box: Path) -> None:
        array = np.random.default_rng(10).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
        array[5:12, 7:14] = float(SENTINEL)  # under the box's centre, where the pass looks first
        tif = write_tiff(tmp_path / "nodata.tif", micro_tiff(array, nodata=SENTINEL))
        vtk, report, _ = run(tmp_path, "--dem", str(tif), "--domain", str(box), "--tolerance", "1")
        assert (vtk.points[:, 2] != float(SENTINEL)).all()
        assert np.isfinite(vtk.points[:, 2]).all()
        assert stats_row(report, "nodata_vertices_removed") == "0"
        assert "nodata_vertices_removed" not in vtk.field_data  # D3 rule 3: a zero is not written
        assert stats_row(report, "start_quality_points_inserted") == "0"
        assert int(stats_row(report, "start_quality_points_skipped")) >= 1


class TestTheVoidIsSkipped:
    """Q-V2 of increment 20's fix (``20-start-quality.md``, "Fix: the
    start-quality pass skips NoData nodes"): the pass never inserts a NoData
    node, so it makes no vertex for the trim to remove and no hole.

    The scene is the box over ``bumpy``. Its two start triangles are bad at
    25°, and the pass adds one node; ``pass_node`` finds it from two runs
    without the void, with the pass on and off, at a tolerance refine
    meets without inserting anything. That node is then NoData, as the
    sentinel or as NaN.

    On the projected path, refine's carving (14 R6) surrounds the NoData
    vertex with its valid neighbours before the trim, so today no valid node
    is left uncovered there, measured: the hole of 23c-2's DC10 is a cut
    piece's. Here the red is the removed vertex and the carving it caused,
    and the coverage and the two oracles (section 3D) are guards that must
    hold before and after.
    """

    HOLES: ClassVar[dict[str, float]] = {"sentinel": float(SENTINEL), "nan": float("nan")}

    @pytest.fixture
    def pass_node(self, tmp_path: Path, bumpy: Path, box: Path) -> tuple[int, int]:
        common = ("--dem", str(bumpy), "--domain", str(box), "--tolerance", "1000")
        on, _, _ = run(tmp_path, *common, name="on.vtk")
        off, _, _ = run(tmp_path, *common, "--start-min-angle", "0", name="off.vtk")
        added = _nodes_of(on, bumpy) - _nodes_of(off, bumpy)
        assert len(added) == 1, added  # the scene: one bad pair, one node
        (node,) = added
        return node

    @pytest.fixture(params=sorted(HOLES))
    def voided(
        self,
        request: pytest.FixtureRequest,
        tmp_path: Path,
        bumpy: Path,
        pass_node: tuple[int, int],
    ) -> tuple[Path, np.ndarray]:
        array = np.asarray(decode_dem(io.BytesIO(bumpy.read_bytes())).array, dtype=np.float32)
        array = array.copy()
        array[pass_node] = self.HOLES[request.param]
        tif = write_tiff(tmp_path / "void.tif", micro_tiff(array, nodata=SENTINEL))
        return tif, array

    def test_no_vertex_is_removed_and_the_pass_counts_the_skip(
        self, tmp_path: Path, box: Path, voided: tuple[Path, np.ndarray]
    ) -> None:
        tif, _ = voided
        vtk, report, _ = run(tmp_path, "--dem", str(tif), "--domain", str(box), "--tolerance", "1")
        assert stats_row(report, "nodata_vertices_removed") == "0"
        assert "nodata_vertices_removed" not in vtk.field_data  # D3 rule 3: a zero is not written
        assert stats_row(report, "points_inserted_on_nodata") == "0"  # no void triangle to carve
        assert stats_row(report, "start_quality_points_inserted") == "0"
        assert int(stats_row(report, "start_quality_points_skipped")) >= 1

    def test_no_output_vertex_is_on_a_nodata_node(
        self, tmp_path: Path, box: Path, voided: tuple[Path, np.ndarray], pass_node: tuple[int, int]
    ) -> None:
        tif, array = voided
        vtk, _, _ = run(tmp_path, "--dem", str(tif), "--domain", str(box), "--tolerance", "1")
        assert pass_node not in _nodes_of(vtk, tif)
        assert np.isfinite(vtk.points[:, 2]).all()
        assert (vtk.points[:, 2] != float(SENTINEL)).all()
        assert int(_is_nodata(array).sum()) == 1  # the scene: exactly the one node

    def test_every_valid_node_is_covered_within_tolerance(
        self, tmp_path: Path, box: Path, voided: tuple[Path, np.ndarray]
    ) -> None:
        """The hole check and section 3D's tolerance oracle, from the file."""
        tif, array = voided
        vtk, _, _ = run(tmp_path, "--dem", str(tif), "--domain", str(box), "--tolerance", "1")
        xy, z = _valid_nodes_inside(tif, array, Polygon(SQUARE).buffer(-SNAP))
        assert len(z) > 100
        error = located_errors(vtk, xy, z)
        assert np.isfinite(error).all(), (
            f"{int((~np.isfinite(error)).sum())} valid nodes in no triangle"
        )
        assert error.max() <= 1.0 + 1e-9

    def test_the_mesh_is_constrained_delaunay(
        self, tmp_path: Path, box: Path, voided: tuple[Path, np.ndarray]
    ) -> None:
        """Section 3D's other oracle, exact, in the file's frame."""
        tif, _ = voided
        vtk, _, _ = run(tmp_path, "--dem", str(tif), "--domain", str(box), "--tolerance", "1")
        constrained = {tuple(sorted(map(int, e))) for e in lines_as_array(vtk)}
        assert _delaunay_violations(vtk.points[:, :2], polygons_as_array(vtk), constrained) == []


def _is_nodata(array: np.ndarray) -> np.ndarray:
    return np.isnan(array) | (array == float(SENTINEL))


def _nodes_of(vtk: VtkFile, tif: Path) -> set[tuple[int, int]]:
    """The (row, col) of every output vertex that is exactly a DEM node."""
    m = decode_dem(io.BytesIO(tif.read_bytes())).meta
    cols = (vtk.points[:, 0] - m.x_min) / m.delta_x
    rows = (m.y_max - vtk.points[:, 1]) / m.delta_y
    on = (cols == np.round(cols)) & (rows == np.round(rows))
    return {(int(r), int(c)) for r, c in zip(rows[on], cols[on], strict=True)}


def _valid_nodes_inside(
    tif: Path, array: np.ndarray, domain: Polygon
) -> tuple[np.ndarray, np.ndarray]:
    """World (x, y) and z of every node with data strictly inside `domain`."""
    m = decode_dem(io.BytesIO(tif.read_bytes())).meta
    r, c = np.indices(array.shape)
    xy = np.column_stack([(m.x_min + c * m.delta_x).ravel(), (m.y_max - r * m.delta_y).ravel()])
    z = array.astype(np.float64).ravel()
    keep = ~_is_nodata(array).ravel() & shapely.contains_xy(domain, xy[:, 0], xy[:, 1])
    return xy[keep], z[keep]


# ---------------------------------------------------------------- T-real


def _under_one_degree(vtk: VtkFile) -> int:
    return int((min_angles_degrees(vtk) < 1.0).sum())


class TestRealTile:
    """T-real: the quarter circle at 10 m and 1 m, M1's columns, logged."""

    @needs_codecs
    def test_the_quarter_circle_with_and_without_the_pass(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        domain = geojson(tmp_path / "quarter.geojson", quarter_circle())
        under: dict[tuple[str, str], int] = {}
        runs: list[tuple[str, str]] = [("10", "0"), ("10", "25"), ("1", "25")]
        for tolerance, angle in runs:
            out = tmp_path / f"q_{tolerance}_{angle}.vtk"
            md = tmp_path / f"q_{tolerance}_{angle}.md"
            began = time.perf_counter()
            mesh(
                "--dem",
                str(KARTVERKET),
                "--domain",
                str(domain),
                "--tolerance",
                tolerance,
                "--start-min-angle",
                angle,
                "--out",
                str(out),
                "--stats",
                str(md),
            )
            elapsed = time.perf_counter() - began
            vtk = read_vtk(out.read_bytes())
            achieved = float(file_field(vtk, "max_error_m"))
            assert achieved <= float(tolerance)
            tris = np.asarray(vtk.polygons, dtype=np.int64)
            degree = np.bincount(tris.ravel(), minlength=len(vtk.points))
            angles = min_angles_degrees(vtk)
            under[(tolerance, angle)] = _under_one_degree(vtk)
            quality = tuple(stats_row(md.read_text(encoding="utf-8"), n) for n in COUNTS)
            with capsys.disabled():
                print(
                    f"\nT-real {tolerance} m, start min angle {angle}: {len(tris)} triangles, "
                    f"{under[(tolerance, angle)]} under 1 deg "
                    f"({100 * np.mean(angles < 1.0):.2f} %), "
                    f"{100 * np.mean(angles < 10.0):.2f} % under 10 deg, "
                    f"{100 * np.mean(angles < 20.0):.2f} % under 20 deg, "
                    f"worst {angles.min():.4f} deg, max degree {degree.max()}, "
                    f">= 12: {(degree >= 12).sum()}, achieved {achieved} m, {elapsed:.1f} s, "
                    f"quality {quality}"
                )
        assert under[("10", "25")] < under[("10", "0")]
