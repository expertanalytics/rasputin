"""`rasputin mesh --dem PATH --tolerance METRES`: increment 14, T7, T9 and T10.

`docs/increments/14-adaptive-refinement.md` R1, R6 and R9, with the user's U1
(a): refinement is opt-in, and without ``--tolerance`` increment 12's mesh is
unchanged. The golden digest below covers the file from ``POINTS`` on (the
mesh and its arrays), so increment 25's header fields do not move it.

Increment 25 (``docs/increments/25-plain-output.md``) replaced R9's sentence
with named fields: the file carries ``tolerance_m``, ``max_error_m`` and
``nodata_vertices_removed`` (only when > 0, D3 rule 3); the start mesh, the
self-check ``dem_nodes_outside_mesh``, ``edge_flips`` and the vertical unit
are ``--stats`` rows (D4). ``file_field`` and ``stats_row`` below are the two
readers every suite imports (D6). ``start_mesh`` for a stride start reads
``every DEM node`` for a stride of 1, else ``every <n><ordinal> DEM node``
with the English ordinal (``every 2nd``, ``every 40th``; "Settled after the
red step", 7).

Increment 14b (``docs/increments/14b-delaunay-insertion.md``): the default
start stride follows the user's C1 (a), at most 129 nodes a side. T10 also
reports min-angle statistics, with no threshold (R10).

T10 runs the real tile under ``needs_codecs`` only and has no ``slow`` marker:
it runs only in the codecs CI job, and takes about 2 s locally at 1 m.
"""

from __future__ import annotations

import hashlib
import io
import math
import re
import time
from pathlib import Path

import numpy as np
import pytest

from cli_driver import USAGE, invoke, ran, rough_dem, runner, write_tiff
from geotiff_fixtures import KARTVERKET, elevations, micro_tiff, needs_codecs
from recordread import (
    SEAMS_AGREE,
    file_field,
    ply_fields,
    start_stride,
    stats_names,
    stats_row,
    stats_rows_named,
)
from test_cli_mesh_dem import SENTINEL, assert_z_is_the_node_value, on_nodes
from tin_engine.cli import app
from tin_engine.io.geotiff import decode_dem
from vtkread import VtkFile, read_vtk

__all__ = [
    "SEAMS_AGREE",
    "file_field",
    "ply_fields",
    "start_stride",
    "stats_names",
    "stats_row",
    "stats_rows_named",
]

#: sha256 of the file from its ``POINTS`` line on: the mesh, the cell arrays
#: and the per-feature arrays, without the dataset FieldData that increment 25
#: rewrites. Recorded at f1765ca (the CLI before increment 25; its whole-file
#: digest was 16b-1/2's ``492ee90e...1df7``) before any increment 25 change.
INCREMENT_12_LARGER_STRIDE_2 = "34c2f7ae2bec9e0a6d76826dfba415c308d555dc213a1051c76d3616f85b0bc3"

NUMBER = r"([0-9.eE+-]+|inf|nan)"


def run(tmp_path: Path, tif: Path, *extra: str) -> tuple[VtkFile, str]:
    out = tmp_path / "x.vtk"
    output = ran("mesh", "--dem", str(tif), "--out", str(out), *extra)
    return read_vtk(out.read_bytes()), output


def run_stats(tmp_path: Path, tif: Path, *extra: str) -> tuple[VtkFile, str, str]:
    """Like ``run``, with ``--stats`` to a file: (mesh, report, stderr alone)."""
    out, md = tmp_path / "x.vtk", tmp_path / "x.md"
    result = runner.invoke(
        app, ["mesh", "--dem", str(tif), "--out", str(out), "--stats", str(md), *extra]
    )
    assert result.exit_code == 0, result.output
    return read_vtk(out.read_bytes()), md.read_text(encoding="utf-8"), result.stderr


#: 17 x 21, seeded rough terrain: refinement has work to do at 1 m.
bumpy = rough_dem(14)


def min_angles_degrees(vtk: VtkFile) -> np.ndarray:
    """Each triangle's smallest angle, in degrees, from its world x and y."""
    tris = np.asarray(vtk.polygons, dtype=np.int64)
    xy = vtk.points[:, :2]
    a, b, c = xy[tris[:, 0]], xy[tris[:, 1]], xy[tris[:, 2]]

    def angle(p: np.ndarray, q: np.ndarray, r: np.ndarray) -> np.ndarray:
        u, v = q - p, r - p
        cos = np.einsum("ij,ij->i", u, v) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1))
        return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))

    return np.minimum(np.minimum(angle(a, b, c), angle(b, c, a)), angle(c, a, b))


class TestRefinedOutput:
    """T9."""

    def test_the_fields_name_the_tolerance_and_a_max_error_within_it(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        vtk, report, stderr = run_stats(tmp_path, bumpy, "--tolerance", "1")
        assert "elevation_source" not in vtk.field_data
        assert file_field(vtk, "tolerance_m") == "1"
        assert float(file_field(vtk, "max_error_m")) <= 1.0
        assert "nodata_vertices_removed" not in vtk.field_data  # a zero count (D3 rule 3)
        assert start_stride(report) == max(1, math.ceil(20 / 128))
        assert stats_row(report, "dem_nodes_outside_mesh") == "0"
        assert stats_row(report, "nodata_vertices_removed") == "0"
        assert "assumed" in stats_row(report, "dem_vertical_unit")
        assert int(stats_row(report, "edge_flips")) >= 0
        # D7: the refine report is gone; one summary line names the triangles.
        assert f"{len(vtk.polygons)} triangles" in stderr, stderr
        assert "not covered" not in stderr
        assert not re.search(r"\b\d+ flips\b", stderr), stderr

    def test_vertices_are_nodes_carrying_their_values(self, tmp_path: Path, bumpy: Path) -> None:
        vtk, report, _ = run_stats(tmp_path, bumpy, "--tolerance", "1", "--stride", "8")
        tile = decode_dem(io.BytesIO(bumpy.read_bytes()))
        assert_z_is_the_node_value(tile, vtk.points)
        rows, cols = on_nodes(tile, vtk.points)
        start = {(r, c) for r in (0, 8, 16) for c in (0, 8, 16, 20)}
        assert start <= set(zip(rows.tolist(), cols.tolist(), strict=True))
        assert len(vtk.points) > len(start), "rough terrain at 1 m must insert points"
        assert start_stride(report) == 8

    def test_a_tighter_tolerance_gives_more_triangles(self, tmp_path: Path, bumpy: Path) -> None:
        coarse, _ = run(tmp_path, bumpy, "--tolerance", "10", "--stride", "8")
        fine, report, _ = run_stats(tmp_path, bumpy, "--tolerance", "0", "--stride", "8")
        assert len(coarse.polygons) < len(fine.polygons)
        # D3 rule 3: a measured value of 0 is written, never omitted.
        assert file_field(fine, "tolerance_m") == "0"
        # 25's D2 and D6: after 15f-3 a DEM node within rounding of a strip
        # vertex is not compared and lifts max_error_m by its difference, which
        # --stats reports; so the bound is max(tolerance_m, that difference).
        at_vertices = float(stats_row(report, "dem_nodes_at_vertices_max_error_m"))
        assert 0.0 <= float(file_field(fine, "max_error_m")) <= max(0.0, at_vertices)

    @pytest.mark.parametrize(("cols", "stride"), [(70, 1), (129, 1), (130, 2), (300, 3)])
    def test_the_default_start_stride_is_at_most_129_nodes_a_side(
        self, tmp_path: Path, cols: int, stride: int
    ) -> None:
        """14b C1 (a): max(1, ceil((max(rows, cols) - 1) / 128)); 300 columns gives 3."""
        assert stride == max(1, math.ceil((cols - 1) / 128))
        array = np.random.default_rng(3).uniform(0.0, 5.0, (4, cols)).astype(np.float32)
        tif = write_tiff(tmp_path / "wide.tif", micro_tiff(array))
        _, report, _ = run_stats(tmp_path, tif, "--tolerance", "1")
        assert start_stride(report) == stride


class TestNoData:
    """T7 through the CLI."""

    def test_a_nodata_row_is_carved_and_no_sentinel_is_written(self, tmp_path: Path) -> None:
        array = np.random.default_rng(5).uniform(0.0, 50.0, (9, 11)).astype(np.float32)
        array[0, :] = float(SENTINEL)
        tif = write_tiff(tmp_path / "nodata.tif", micro_tiff(array, nodata=SENTINEL))
        vtk, report, _ = run_stats(tmp_path, tif, "--tolerance", "1")
        tile = decode_dem(io.BytesIO(tif.read_bytes()))
        assert (vtk.points[:, 2] != float(SENTINEL)).all()
        assert_z_is_the_node_value(tile, vtk.points)
        assert stats_row(report, "dem_nodes_outside_mesh") == "0"
        removed = file_field(vtk, "nodata_vertices_removed")
        assert int(removed) > 0
        assert stats_row(report, "nodata_vertices_removed") == removed

    def test_an_all_nodata_tile_is_exit_2_and_writes_nothing(self, tmp_path: Path) -> None:
        array = np.full((3, 4), float(SENTINEL), dtype=np.float32)
        tif = write_tiff(tmp_path / "void.tif", micro_tiff(array, nodata=SENTINEL))
        out = tmp_path / "x.vtk"
        code, output = invoke("mesh", "--dem", str(tif), "--tolerance", "1", "--out", str(out))
        assert code == USAGE, output
        assert "No such option" not in output, "refused for the wrong reason"
        assert not out.exists()


class TestRefusals:
    """R9: --tolerance must be finite and >= 0, and needs --dem."""

    @pytest.mark.parametrize("value", ["-1", "nan", "inf", "-inf"])
    def test_a_bad_tolerance(self, tmp_path: Path, bumpy: Path, value: str) -> None:
        out = tmp_path / "x.vtk"
        code, output = invoke("mesh", "--dem", str(bumpy), "--tolerance", value, "--out", str(out))
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--tolerance" in output
        assert not out.exists()

    def test_tolerance_with_a_fixture(self, tmp_path: Path) -> None:
        out = tmp_path / "x.vtk"
        code, output = invoke("mesh", "catchment", "--tolerance", "1", "--out", str(out))
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--tolerance" in output
        assert not out.exists()


class TestWithoutTolerance:
    """U1 (a): increment 12's output, byte for byte."""

    def test_the_uniform_mesh_is_unchanged(self, tmp_path: Path) -> None:
        tif = write_tiff(tmp_path / "larger.tif", micro_tiff(elevations(rows=7, cols=9)))
        out = tmp_path / "x.vtk"
        code, output = invoke("mesh", "--dem", str(tif), "--stride", "2", "--out", str(out))
        assert code == 0, output
        blob = out.read_bytes()
        mesh = blob[blob.index(b"\nPOINTS ") + 1 :]
        assert hashlib.sha256(mesh).hexdigest() == INCREMENT_12_LARGER_STRIDE_2
        vtk = read_vtk(blob)
        assert "tolerance_m" not in vtk.field_data and "max_error_m" not in vtk.field_data


class TestRealTile:
    """T10: 7908_3_10m_z33.tif at 5 m and 1 m. No timing or angle assertion.

    14b: min-angle statistics are reported (world coordinates, each
    triangle's smallest angle), with no threshold: M2 makes the sea's
    grading a property of the data and the stride, not of the algorithm.
    """

    @needs_codecs
    def test_refines_the_real_dem(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
        polygons: dict[str, int] = {}
        for tolerance in ("5", "1"):
            began = time.perf_counter()
            vtk, output = run(tmp_path, KARTVERKET, "--tolerance", tolerance)
            seconds = time.perf_counter() - began
            text = {k: file_field(vtk, k) for k in ("tolerance_m", "max_error_m")}
            achieved = float(text["max_error_m"])
            assert achieved <= float(tolerance)
            polygons[tolerance] = len(vtk.polygons)
            angles = min_angles_degrees(vtk)
            with capsys.disabled():
                print(
                    f"\nT10 tolerance {tolerance} m: {len(vtk.polygons)} triangles, "
                    f"{len(vtk.points)} vertices, achieved {achieved} m, {seconds:.1f} s; {text}"
                    f"\nT10 min angle: median {np.median(angles):.2f} deg, "
                    f"{100 * np.mean(angles < 1.0):.1f} % under 1 deg, "
                    f"{100 * np.mean(angles < 10.0):.1f} % under 10 deg, "
                    f"worst {angles.min():.4f} deg; report: {output.strip()}"
                )
        assert polygons["5"] < polygons["1"]
