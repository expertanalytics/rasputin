"""`rasputin mesh --dem PATH`: increment 12, tests 13-19 (design R5, R6).

Micro-TIFFs from `geotiff_fixtures.py`, written to `tmp_path`. The oracle for
node positions and values is `decode_dem`, which increment 11 shipped and
tested: each written vertex must sit exactly on a node of the decoded grid and
carry that node's value. The `.vtk` is parsed with `vtkread`, the `.ply` with
`plyread`; nothing asks the writer what it wrote.

Every refusal asserts exit 2, what it says, and that no file was left behind.

The real fixture runs only under `needs_codecs` (test 18) or, for the refusal,
only without the extra (test 19).
"""

from __future__ import annotations

import io
import math
import re
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from geotiff_fixtures import (
    KARTVERKET,
    LZW,
    METRE,
    REFUSALS,
    VERTICAL_UNITS,
    elevations,
    micro_tiff,
    needs_codecs,
    with_compression_tag,
    with_keys,
    without_codecs,
)
from plyread import read_ply, vertex_array
from recordread import file_field, start_stride, stats_row
from test_cli_mesh import plain
from tin_engine.cli import app
from tin_engine.io.geotiff import decode_dem
from tin_engine.io.models import DemTile
from vtkread import VtkFile, read_vtk

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

SENTINEL = "-32767"
USAGE = 2


def invoke(*args: str) -> tuple[int, str]:
    result = runner.invoke(app, ["mesh", *args])
    return result.exit_code, plain(result.output)


def write_tiff(path: Path, stream: io.BytesIO) -> Path:
    path.write_bytes(stream.getvalue())
    return path


def on_nodes(tile: DemTile, xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The (row, col) of each vertex, required to be exactly a grid node."""
    m = tile.meta
    cols = (xyz[:, 0] - m.x_min) / m.delta_x
    rows = (m.y_max - xyz[:, 1]) / m.delta_y
    assert_array_integral(cols)
    assert_array_integral(rows)
    return rows.astype(int), cols.astype(int)


def assert_array_integral(values: np.ndarray) -> None:
    assert (values == np.round(values)).all(), f"off-node vertices: {values}"


def assert_z_is_the_node_value(tile: DemTile, xyz: np.ndarray) -> None:
    rows, cols = on_nodes(tile, xyz)
    expected = tile.array[rows, cols].astype(np.float64)
    assert xyz[:, 2] == pytest.approx(expected, rel=1e-9, abs=1e-9)


@pytest.fixture
def baseline(tmp_path: Path) -> Path:
    return write_tiff(tmp_path / "tile.tif", micro_tiff())


@pytest.fixture
def larger(tmp_path: Path) -> Path:
    """7 x 9 so a stride can skip nodes, with distinct values."""
    return write_tiff(tmp_path / "larger.tif", micro_tiff(elevations(rows=7, cols=9)))


def run_vtk(tmp_path: Path, tif: Path, *extra: str) -> VtkFile:
    out = tmp_path / "x.vtk"
    code, output = invoke("--dem", str(tif), "--out", str(out), *extra)
    assert code == 0, output
    return read_vtk(out.read_bytes())


def stats_of(tmp_path: Path, tif: Path, *extra: str) -> tuple[VtkFile, str]:
    """``run_vtk`` with ``--stats`` beside the mesh: (mesh, report)."""
    md = tmp_path / "x.md"
    vtk = run_vtk(tmp_path, tif, "--stats", str(md), *extra)
    return vtk, md.read_text(encoding="utf-8")


class TestTheAcceptanceInvocation:
    """Test 13."""

    def test_writes_the_tile_values_at_its_nodes(self, tmp_path: Path, baseline: Path) -> None:
        vtk = run_vtk(tmp_path, baseline)
        tile = decode_dem(io.BytesIO(baseline.read_bytes()))
        assert len(vtk.points) == tile.meta.rows * tile.meta.cols  # default stride 1
        assert_z_is_the_node_value(tile, vtk.points)
        assert len(vtk.polygons) == 2 * (tile.meta.rows - 1) * (tile.meta.cols - 1)

    def test_elevation_is_a_point_array_equal_to_z(self, tmp_path: Path, baseline: Path) -> None:
        """Increment 12 amendment: ParaView's Color By -> elevation shows the heights."""
        vtk = run_vtk(tmp_path, baseline)
        assert np.array_equal(np.asarray(vtk.point_scalars["elevation"].values), vtk.points[:, 2])
        assert "elevation" not in vtk.field_data

    def test_records_the_crs_and_the_stride(self, tmp_path: Path, baseline: Path) -> None:
        """Increment 25: the file names the CRS and the DEM; the stride is
        ``start_mesh`` in ``--stats`` (assumed on the no-tolerance path too,
        the mesh being the stride grid itself)."""
        vtk, report = stats_of(tmp_path, baseline)
        assert vtk.field_data["crs"].values == ("EPSG:25833",)
        assert file_field(vtk, "dem_source") == "tile.tif"
        assert "elevation_source" not in vtk.field_data
        for absent in ("tolerance_m", "max_error_m"):  # no tolerance, no tolerance fields
            assert absent not in vtk.field_data, absent
        assert start_stride(report) == 1

    def test_an_assumed_vertical_unit_is_said(self, tmp_path: Path, baseline: Path) -> None:
        _, report = stats_of(tmp_path, baseline)
        unit = stats_row(report, "dem_vertical_unit")
        assert "metres" in unit and "assumed" in unit, unit

    def test_a_declared_vertical_unit_is_not_called_assumed(self, tmp_path: Path) -> None:
        tif = write_tiff(tmp_path / "m.tif", micro_tiff(geokeys=with_keys({VERTICAL_UNITS: METRE})))
        _, report = stats_of(tmp_path, tif)
        unit = stats_row(report, "dem_vertical_unit")
        assert "metres" in unit and "assumed" not in unit, unit

    def test_stride_picks_every_nth_node_and_the_last(self, tmp_path: Path, larger: Path) -> None:
        vtk = run_vtk(tmp_path, larger, "--stride", "3")
        tile = decode_dem(io.BytesIO(larger.read_bytes()))
        rows, cols = on_nodes(tile, vtk.points)
        assert set(rows.tolist()) == {0, 3, 6}
        assert set(cols.tolist()) == {0, 3, 6, 8}
        assert len(vtk.points) == 3 * 4
        assert_z_is_the_node_value(tile, vtk.points)
        _, report = stats_of(tmp_path, larger, "--stride", "3")
        assert start_stride(report) == 3

    def test_the_default_stride_is_the_formula(self, tmp_path: Path) -> None:
        """300 columns: the smallest stride giving at most 256 per side is 2."""
        tif = write_tiff(tmp_path / "wide.tif", micro_tiff(elevations(rows=3, cols=300)))
        _, report = stats_of(tmp_path, tif)
        assert start_stride(report) == max(1, math.ceil(299 / 255))


class TestNoData:
    """Tests 14 and 15, the user's U1 (a): drop and count."""

    def test_a_nodata_edge_row_is_dropped_and_counted(self, tmp_path: Path) -> None:
        array = elevations(rows=5, cols=6)
        array[0, :] = float(SENTINEL)
        tif = write_tiff(tmp_path / "nodata.tif", micro_tiff(array, nodata=SENTINEL))
        out = tmp_path / "x.vtk"
        code, output = invoke("--dem", str(tif), "--out", str(out))
        assert code == 0, output
        vtk = read_vtk(out.read_bytes())
        tile = decode_dem(io.BytesIO(tif.read_bytes()))
        assert len(vtk.points) == 4 * 6
        assert_z_is_the_node_value(tile, vtk.points)
        assert (vtk.points[:, 2] != float(SENTINEL)).all()
        assert file_field(vtk, "nodata_vertices_removed") == "6"
        # R3, worded by increment 25's summary line (D7), one wording on both
        # paths since increment 27. The count is 6 under either sampler rule:
        # the NoData row is row 0, and a node's bilinear cell reaches down and
        # right, never up, so no valid node had a NoData corner.
        assert re.search(r"\b6 vertices on NoData cells\b", output), output
        assert "next to" not in output, output
        assert "without data dropped" not in output

    def test_an_all_nodata_tile_writes_nothing(self, tmp_path: Path) -> None:
        array = np.full((3, 4), float(SENTINEL), dtype=np.float32)
        tif = write_tiff(tmp_path / "void.tif", micro_tiff(array, nodata=SENTINEL))
        out = tmp_path / "x.vtk"
        code, output = invoke("--dem", str(tif), "--out", str(out))
        assert code != 0, output
        assert "No such option" not in output, "refused for the wrong reason"
        assert not out.exists()


class TestUsageErrors:
    """Test 16: exit 2, the reason in words, and no file."""

    @staticmethod
    def refused(tmp_path: Path, *args: str, says: tuple[str, ...]) -> None:
        out = tmp_path / "x.vtk"
        code, output = invoke(*args, "--out", str(out))
        assert code == USAGE, output
        # `--dem` must be a known option, or every refusal here passes as
        # Typer's "No such option" before the command runs at all.
        assert "No such option" not in output, output
        for word in says:
            assert word in output, f"{word!r} not in {output!r}"
        assert not out.exists()

    def test_both_a_fixture_and_dem(self, tmp_path: Path, baseline: Path) -> None:
        self.refused(tmp_path, "catchment", "--dem", str(baseline), says=("--dem",))

    def test_neither_a_fixture_nor_dem(self, tmp_path: Path) -> None:
        self.refused(tmp_path, says=("--dem",))

    def test_dem_with_flat(self, tmp_path: Path, baseline: Path) -> None:
        self.refused(tmp_path, "--dem", str(baseline), "--flat", says=("--flat",))

    def test_dem_with_crs(self, tmp_path: Path, baseline: Path) -> None:
        self.refused(tmp_path, "--dem", str(baseline), "--crs", "EPSG:4326", says=("--crs",))

    def test_a_missing_file(self, tmp_path: Path) -> None:
        missing = tmp_path / "nowhere.tif"
        self.refused(tmp_path, "--dem", str(missing), says=("nowhere.tif",))

    def test_a_directory(self, tmp_path: Path) -> None:
        self.refused(tmp_path, "--dem", str(tmp_path), says=(tmp_path.name,))

    @pytest.mark.parametrize("stride", ["0", "-3"])
    def test_a_non_positive_stride(self, tmp_path: Path, baseline: Path, stride: str) -> None:
        self.refused(tmp_path, "--dem", str(baseline), "--stride", stride, says=("--stride",))

    @pytest.mark.parametrize(
        "refusal",
        # The codec case depends on whether the extra is installed; it has its
        # own test below. The override case needs a caller-supplied sentinel.
        [r for r in REFUSALS if not r.decode_kwargs and r.name != "refuses_missing_codec"],
        ids=lambda r: r.name,
    )
    def test_a_file_decode_dem_refuses_shows_its_message(
        self, tmp_path: Path, refusal: object
    ) -> None:
        from geotiff_fixtures import Refusal

        assert isinstance(refusal, Refusal)
        tif = write_tiff(tmp_path / "bad.tif", refusal.build())
        self.refused(tmp_path, "--dem", str(tif), says=refusal.must_name)

    def test_not_a_tiff_at_all(self, tmp_path: Path) -> None:
        junk = tmp_path / "junk.tif"
        junk.write_bytes(b"this is not a TIFF file\n" * 4)
        self.refused(tmp_path, "--dem", str(junk), says=("TIFF structure",))

    @without_codecs
    def test_a_missing_codec_points_at_the_extra(self, tmp_path: Path) -> None:
        tif = write_tiff(tmp_path / "lzw.tif", with_compression_tag(micro_tiff(), LZW))
        self.refused(tmp_path, "--dem", str(tif), says=("codecs",))


class TestPly:
    """Test 17: the same z in a `.ply` as in the `.vtk`."""

    def test_ply_carries_the_same_vertices(self, tmp_path: Path, larger: Path) -> None:
        vtk = run_vtk(tmp_path, larger, "--stride", "2")
        out = tmp_path / "x.ply"
        code, output = invoke("--dem", str(larger), "--stride", "2", "--out", str(out))
        assert code == 0, output
        _, data = read_ply(out.read_bytes())
        ply = vertex_array(data)
        order = np.lexsort(ply.T[:2])
        assert np.array_equal(ply[order], vtk.points[np.lexsort(vtk.points.T[:2])])


class TestRealFixture:
    """Tests 18 and 19, on `tests/fixtures/dem_archive/7908_3_10m_z33.tif`."""

    @needs_codecs
    def test_meshes_the_real_dem(self, tmp_path: Path) -> None:
        out, md = tmp_path / "tile.vtk", tmp_path / "tile.md"
        code, output = invoke("--dem", str(KARTVERKET), "--out", str(out), "--stats", str(md))
        assert code == 0, output
        vtk = read_vtk(out.read_bytes())
        z = vtk.points[:, 2]
        assert np.isfinite(z).all()
        assert z.min() >= -1.3 and z.max() <= 391.8
        assert vtk.field_data["crs"].values == ("EPSG:25833",)
        assert start_stride(md.read_text(encoding="utf-8")) == 20
        assert len(vtk.points) < 254 * 254, "the NoData outline drops some vertices"
        assert int(file_field(vtk, "nodata_vertices_removed")) > 0

        vtk_module = pytest.importorskip("vtk")
        reader = vtk_module.vtkPolyDataReader()
        reader.SetFileName(str(out))
        reader.Update()
        bounds = reader.GetOutput().GetBounds()
        assert bounds[4] >= -1.3 and bounds[5] <= 391.8

    @without_codecs
    def test_without_the_extra_it_is_a_usage_error(self, tmp_path: Path) -> None:
        out = tmp_path / "tile.vtk"
        code, output = invoke("--dem", str(KARTVERKET), "--out", str(out))
        assert code == USAGE, output
        assert "codecs" in output
        assert "Traceback" not in output
        assert not out.exists()
