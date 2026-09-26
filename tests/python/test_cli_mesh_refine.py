"""`rasputin mesh --dem PATH --tolerance METRES`: increment 14, T7, T9 and T10.

`docs/increments/14-adaptive-refinement.md` R1, R6 and R9, with the user's U1
(a): refinement is opt-in, and without ``--tolerance`` increment 12's output is
unchanged byte for byte. The golden digest below was taken from increment 12's
shipped CLI (merge 31629d6) on the ``larger`` fixture at ``--stride 2``.

Wording pinned from R9's example sentence: ``refined from DEM nodes``,
``tolerance <t> m``, ``achieved max error <e> m``, ``start stride <s>``,
``<n> valid DEM nodes not covered`` and ``<k> vertices without data dropped``.
How ``<t>`` is formatted is not ruled, so it is parsed as a float.

Increment 14b (``docs/increments/14b-delaunay-insertion.md``): the default
start stride follows the user's C1 (a), at most 129 nodes a side; the sentence
says ``constrained Delaunay``; the stderr report carries ``<n> flips`` (the
design says "flips in the stderr report" and leaves the wording open; this
suite picks ``<n> flips``). T10 also reports min-angle statistics, with no
threshold (R10).

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

from geotiff_fixtures import KARTVERKET, elevations, micro_tiff, needs_codecs
from test_cli_mesh_dem import (
    SENTINEL,
    USAGE,
    assert_z_is_the_node_value,
    invoke,
    on_nodes,
    write_tiff,
)
from tin_engine.io.geotiff import decode_dem
from vtkread import VtkFile, read_vtk

INCREMENT_12_LARGER_STRIDE_2 = "c677478e7bc132531d4fb9d799f5a9288658ac1a8a50317c3e78f79d376d09fa"

NUMBER = r"([0-9.eE+-]+|inf|nan)"


def sentence(vtk: VtkFile) -> str:
    (text,) = vtk.field_data["elevation_source"].values
    assert isinstance(text, str)
    return text


def field(text: str, pattern: str) -> float:
    match = re.search(pattern, text)
    assert match is not None, f"{pattern!r} not in {text!r}"
    return float(match.group(1))


def run(tmp_path: Path, tif: Path, *extra: str) -> tuple[VtkFile, str]:
    out = tmp_path / "x.vtk"
    code, output = invoke("--dem", str(tif), "--out", str(out), *extra)
    assert code == 0, output
    return read_vtk(out.read_bytes()), output


@pytest.fixture
def bumpy(tmp_path: Path) -> Path:
    """17 x 21, seeded rough terrain: refinement has work to do at 1 m."""
    array = np.random.default_rng(14).uniform(0.0, 50.0, (17, 21)).astype(np.float32)
    return write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))


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

    def test_the_sentence_names_the_tolerance_and_an_achieved_error_within_it(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        vtk, output = run(tmp_path, bumpy, "--tolerance", "1")
        text = sentence(vtk)
        assert "refined from DEM nodes" in text
        assert "bilinear" not in text
        assert field(text, rf"tolerance {NUMBER} m") == 1.0
        assert field(text, rf"achieved max error {NUMBER} m") <= 1.0
        assert field(text, rf"start stride {NUMBER}") == max(1, math.ceil(20 / 128))
        assert "constrained Delaunay" in text
        assert field(text, rf"{NUMBER} valid DEM nodes not covered") == 0
        assert field(text, rf"{NUMBER} vertices without data dropped") == 0
        assert "vertical unit assumed metres" in text
        # R9: the same counts go to stderr. The wording there is not ruled.
        assert "not covered" in output
        assert re.search(r"\b\d+ flips\b", output), output

    def test_vertices_are_nodes_carrying_their_values(self, tmp_path: Path, bumpy: Path) -> None:
        vtk, _ = run(tmp_path, bumpy, "--tolerance", "1", "--stride", "8")
        tile = decode_dem(io.BytesIO(bumpy.read_bytes()))
        assert_z_is_the_node_value(tile, vtk.points)
        rows, cols = on_nodes(tile, vtk.points)
        start = {(r, c) for r in (0, 8, 16) for c in (0, 8, 16, 20)}
        assert start <= set(zip(rows.tolist(), cols.tolist(), strict=True))
        assert len(vtk.points) > len(start), "rough terrain at 1 m must insert points"
        assert field(sentence(vtk), rf"start stride {NUMBER}") == 8

    def test_a_tighter_tolerance_gives_more_triangles(self, tmp_path: Path, bumpy: Path) -> None:
        coarse, _ = run(tmp_path, bumpy, "--tolerance", "10", "--stride", "8")
        fine, _ = run(tmp_path, bumpy, "--tolerance", "0", "--stride", "8")
        assert len(coarse.polygons) < len(fine.polygons)
        assert field(sentence(fine), rf"achieved max error {NUMBER} m") == 0.0

    @pytest.mark.parametrize(("cols", "stride"), [(70, 1), (129, 1), (130, 2), (300, 3)])
    def test_the_default_start_stride_is_at_most_129_nodes_a_side(
        self, tmp_path: Path, cols: int, stride: int
    ) -> None:
        """14b C1 (a): max(1, ceil((max(rows, cols) - 1) / 128)); 300 columns gives 3."""
        assert stride == max(1, math.ceil((cols - 1) / 128))
        array = np.random.default_rng(3).uniform(0.0, 5.0, (4, cols)).astype(np.float32)
        tif = write_tiff(tmp_path / "wide.tif", micro_tiff(array))
        vtk, _ = run(tmp_path, tif, "--tolerance", "1")
        assert field(sentence(vtk), rf"start stride {NUMBER}") == stride


class TestNoData:
    """T7 through the CLI."""

    def test_a_void_row_is_carved_and_no_sentinel_is_written(self, tmp_path: Path) -> None:
        array = np.random.default_rng(5).uniform(0.0, 50.0, (9, 11)).astype(np.float32)
        array[0, :] = float(SENTINEL)
        tif = write_tiff(tmp_path / "void.tif", micro_tiff(array, nodata=SENTINEL))
        vtk, _ = run(tmp_path, tif, "--tolerance", "1")
        tile = decode_dem(io.BytesIO(tif.read_bytes()))
        assert (vtk.points[:, 2] != float(SENTINEL)).all()
        assert_z_is_the_node_value(tile, vtk.points)
        text = sentence(vtk)
        assert field(text, rf"{NUMBER} valid DEM nodes not covered") >= 0
        assert field(text, rf"{NUMBER} vertices without data dropped") > 0

    def test_an_all_nodata_tile_is_exit_2_and_writes_nothing(self, tmp_path: Path) -> None:
        array = np.full((3, 4), float(SENTINEL), dtype=np.float32)
        tif = write_tiff(tmp_path / "void.tif", micro_tiff(array, nodata=SENTINEL))
        out = tmp_path / "x.vtk"
        code, output = invoke("--dem", str(tif), "--tolerance", "1", "--out", str(out))
        assert code == USAGE, output
        assert "No such option" not in output, "refused for the wrong reason"
        assert not out.exists()


class TestRefusals:
    """R9: --tolerance must be finite and >= 0, and needs --dem."""

    @pytest.mark.parametrize("value", ["-1", "nan", "inf", "-inf"])
    def test_a_bad_tolerance(self, tmp_path: Path, bumpy: Path, value: str) -> None:
        out = tmp_path / "x.vtk"
        code, output = invoke("--dem", str(bumpy), "--tolerance", value, "--out", str(out))
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--tolerance" in output
        assert not out.exists()

    def test_tolerance_with_a_fixture(self, tmp_path: Path) -> None:
        out = tmp_path / "x.vtk"
        code, output = invoke("catchment", "--tolerance", "1", "--out", str(out))
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "--tolerance" in output
        assert not out.exists()


class TestWithoutTolerance:
    """U1 (a): increment 12's output, byte for byte."""

    def test_the_uniform_mesh_is_unchanged(self, tmp_path: Path) -> None:
        tif = write_tiff(tmp_path / "larger.tif", micro_tiff(elevations(rows=7, cols=9)))
        out = tmp_path / "x.vtk"
        code, output = invoke("--dem", str(tif), "--stride", "2", "--out", str(out))
        assert code == 0, output
        assert hashlib.sha256(out.read_bytes()).hexdigest() == INCREMENT_12_LARGER_STRIDE_2


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
            text = sentence(vtk)
            achieved = field(text, rf"achieved max error {NUMBER} m")
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
