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

import re
import time
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import tin_engine.cli as cli
from geotiff_fixtures import KARTVERKET, micro_tiff, needs_codecs
from test_cli_mesh_dem import SENTINEL, USAGE, invoke, write_tiff
from test_cli_mesh_domain import SQUARE, geojson, quarter_circle
from test_cli_mesh_refine import file_field, min_angles_degrees, stats_row
from test_cli_mesh_stats import mesh, seconds
from tin_engine import _core
from vtkread import VtkFile, read_vtk

ROWS, COLS = 17, 21
QUALITY = re.compile(r"(\d+) start quality nodes inserted, (\d+) start quality skips")  # pre-25
COUNTS = ("start_quality_points_inserted", "start_quality_points_skipped")


@pytest.fixture
def bumpy(tmp_path: Path) -> Path:
    array = np.random.default_rng(20).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
    return write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))


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
        code, output = invoke("--dem", str(bumpy), "--start-min-angle", "25", "--out", str(target))
        assert code == USAGE, output
        assert "No such option" not in output  # refused for its value, not unknown
        assert "--start-min-angle" in output
        assert "--tolerance" in output
        assert not target.exists()

    def test_zero_without_tolerance_is_refused_too(self, tmp_path: Path, bumpy: Path) -> None:
        """Given at all, it needs --tolerance; 0 is a value, not an absence."""
        code, output = invoke(
            "--dem", str(bumpy), "--start-min-angle", "0", "--out", str(tmp_path / "x.vtk")
        )
        assert code == USAGE, output
        assert "No such option" not in output  # refused for its value, not unknown
        assert "--start-min-angle" in output

    def test_without_dem(self, tmp_path: Path) -> None:
        code, output = invoke(
            "catchment", "--flat", "--start-min-angle", "25", "--out", str(tmp_path / "x.vtk")
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
    """Q10 (R9): the pass may insert a NoData node; refine ends and trim drops it."""

    def test_a_domain_over_a_nodata_block(self, tmp_path: Path, box: Path) -> None:
        array = np.random.default_rng(10).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
        array[5:12, 7:14] = float(SENTINEL)  # under the box's centre, where the pass looks first
        tif = write_tiff(tmp_path / "nodata.tif", micro_tiff(array, nodata=SENTINEL))
        vtk, report, _ = run(tmp_path, "--dem", str(tif), "--domain", str(box), "--tolerance", "1")
        assert (vtk.points[:, 2] != float(SENTINEL)).all()
        assert np.isfinite(vtk.points[:, 2]).all()
        assert int(file_field(vtk, "nodata_vertices_removed")) >= 1
        assert int(stats_row(report, "start_quality_points_inserted")) > 0


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
