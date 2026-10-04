"""`rasputin mesh`'s plain output through the CLI: increment 25.

`docs/increments/25-plain-output.md`, "Tests for @tester": one vocabulary
across the mesh file (`.vtk` and `.ply` alike), `--stats`, stderr and
`--record`; the values against independent counts; the other paths; stderr's
banned words on successful runs; and `--record PATH` (D5). The record alone
is `test_run_record.py`; the reworded stderr lines of the features, land-cover
and catchment paths are pinned in their own suites, and the reprojected
`max_error_m` in `test_cli_mesh_geographic.py`.

Readers come from `recordread.py` (`file_field`, `stats_row`, `stats_names`,
`ply_fields`), whose docstring states the `--stats` layout read here.

`command` is built from `sys.argv` (`cli.py`, `_write_report`), which
`CliRunner` does not set, so a test that reads it sets `sys.argv` with
`monkeypatch` to the argument list it invokes.

PINNED HERE, where D5 is silent: `--record`'s refusals are usage errors
(exit 2) naming `--record`, as `--stats`'s are; its path is echoed on stdout
on its own line after the mesh and `--stats` paths; `PATH` outside
`--out-parent` is refused like any other output.

External input: none new. `--record` writes a file and reads nothing, so
`tester.md` section 3C does not apply.

Not invariant-critical, so no mutation round.
"""

from __future__ import annotations

import io
import json
import re
import shlex
import socket
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from typer.testing import CliRunner, Result

import tin_engine.cli as cli
from feature_fixtures import Feat, write_geojson
from geotiff_fixtures import micro_tiff
from mosaic_fixtures import quadrants, whole
from plyread import parse_header
from recordread import file_field, ply_fields, stats_names, stats_row
from test_cli_mesh import plain
from test_cli_mesh_dem import SENTINEL, write_tiff
from test_cli_mesh_domain import SQUARE, geojson
from test_cli_mesh_features import FOREST
from test_cli_mesh_mosaic import terrain, write_tiles
from tin_engine import installed_version
from tin_engine.cli import app
from tin_engine.io.geotiff import decode_dem
from vtkread import VtkFile, read_vtk

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})
USAGE = 2
ROWS, COLS = 17, 21
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
VELHAS = FIXTURES / "velhas"

#: D3 rule 6, matched case-blind as whole words (``snap`` also in
#: ``snapped``), over the whole of stderr of a successful run.
BANNED = re.compile(
    r"(?i)\buncovered\b|\bfeet\b|\bstride\b|\bchains\b|\bcoincident\b|\bcarved\b|"
    r"\bnoded\b|R-tree|\bsnap|\bDEM holes\b|\bvoid\b"
)
#: The writers' own fields, describing the file's arrays (D2, "Fields that stay").
WRITER_FIELDS = {"feature_bits", "feature_names", "feature_vocabulary"}
#: D2's mesh-file fields plus the two kept because a licence or an array needs them.
FILE_FIELDS = {
    "crs",
    "tolerance_m",
    "max_error_m",
    "dem_source",
    "dem_credit",
    "licence_note",
    "cite",
    "nodata_vertices_removed",
    "heights",
    "features_notice",
    "land_cover_codes",
}


@pytest.fixture
def holed(tmp_path: Path) -> Path:
    """17 x 21 rough terrain with a block of NoData cells inside it, so a
    refined and a strided run both lose vertices (the "projected fixture
    with NoData" of the design's tests)."""
    array = np.random.default_rng(25).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
    array[6:11, 8:13] = float(SENTINEL)
    return write_tiff(tmp_path / "holed.tif", micro_tiff(array, nodata=SENTINEL))


def mesh(monkeypatch: pytest.MonkeyPatch, *args: str) -> Result:
    """``rasputin mesh *args`` with ``sys.argv`` set as a shell would set it."""
    monkeypatch.setattr(sys, "argv", ["rasputin", "mesh", *args])
    result = runner.invoke(app, ["mesh", *args])
    assert result.exit_code == 0, result.output
    return result


def refused(monkeypatch: pytest.MonkeyPatch, *args: str) -> str:
    monkeypatch.setattr(sys, "argv", ["rasputin", "mesh", *args])
    result = runner.invoke(app, ["mesh", *args])
    output = plain(result.output)
    assert result.exit_code == USAGE, output
    assert "No such option" not in output, "refused for the wrong reason"
    assert "Traceback" not in output
    return output


def vtk_fields(vtk: VtkFile) -> dict[str, str]:
    """Every run field of a ``.vtk``, the writers' own excluded."""
    return {name: file_field(vtk, name) for name in vtk.field_data if name not in WRITER_FIELDS}


def located_errors(vtk: VtkFile, xy: np.ndarray, z: np.ndarray) -> np.ndarray:
    """``|plane - z|`` at each node in a written triangle (closed, 1e-12
    slack; a node on an edge counts in both), the plane from the file's own
    vertices; -inf for a node in no triangle. Written here, not asked of
    the engine."""
    tris = np.asarray(vtk.polygons, dtype=np.int64)
    origin = vtk.points[:, :2].min(axis=0)
    v = vtk.points[:, :2] - origin
    a, b, c = v[tris[:, 0]], v[tris[:, 1]], v[tris[:, 2]]
    za, zb, zc = (vtk.points[tris[:, k], 2] for k in range(3))
    two_a = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    x, y = xy[:, 0, None] - origin[0], xy[:, 1, None] - origin[1]

    def weight(p: np.ndarray, q: np.ndarray) -> np.ndarray:
        return ((q[:, 0] - p[:, 0]) * (y - p[:, 1]) - (q[:, 1] - p[:, 1]) * (x - p[:, 0])) / two_a

    wa, wb, wc = weight(b, c), weight(c, a), weight(a, b)
    holds = (wa >= -1e-12) & (wb >= -1e-12) & (wc >= -1e-12)
    error = np.abs(wa * za + wb * zb + wc * zc - z[:, None])
    return np.where(holds, error, -np.inf).max(axis=1)


def dem_nodes(tif: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Every node's world (x, y), its z, and whether it has data."""
    tile = decode_dem(io.BytesIO(tif.read_bytes()))
    m = tile.meta
    r, c = np.indices((m.rows, m.cols))
    xy = np.column_stack([(m.x_min + c * m.delta_x).ravel(), (m.y_max - r * m.delta_y).ravel()])
    z = np.asarray(tile.array, dtype=np.float64).ravel()
    return xy, z, z != float(SENTINEL)


def picked_without_data(tif: Path, step: int) -> tuple[int, int]:
    """(picked nodes that are NoData, nodes picked) at ``step``: every
    ``step``-th row and column, and the last. Increment 27: a node reads only
    itself, so these are exactly the vertices the no-tolerance path removes."""
    tile = decode_dem(io.BytesIO(tif.read_bytes()))
    nodata = np.asarray(tile.array) == float(SENTINEL)
    rows = sorted({*range(0, ROWS, step), ROWS - 1})
    cols = sorted({*range(0, COLS, step), COLS - 1})
    return int(nodata[np.ix_(rows, cols)].sum()), len(rows) * len(cols)


# ---------------------------------------------------------------- one vocabulary


class TestOneVocabulary:
    """D1 and D3 rule 1: one name, one value, in every place it is printed."""

    def test_every_file_field_is_in_stats_with_the_same_value(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        out, md = tmp_path / "x.vtk", tmp_path / "x.md"
        mesh(monkeypatch, "--dem", str(holed), "--tolerance", "1", "--out", str(out),
             "--stats", str(md))  # fmt: skip
        fields = vtk_fields(read_vtk(out.read_bytes()))
        assert set(fields) == {
            "crs",
            "tolerance_m",
            "max_error_m",
            "dem_source",
            "nodata_vertices_removed",
        }
        report = md.read_text(encoding="utf-8")
        for name, value in fields.items():
            assert stats_row(report, name) == value, name

    @pytest.mark.parametrize("extra", [("--tolerance", "1"), ("--stride", "2")])
    def test_the_ply_carries_the_vtks_fields_name_for_name(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, extra: tuple[str, ...]
    ) -> None:
        """D2: the `.ply` carries the same list (it lacked nine fields before)."""
        vtk, ply, edges = tmp_path / "x.vtk", tmp_path / "x.ply", tmp_path / "e.ply"
        mesh(monkeypatch, "--dem", str(holed), *extra, "--out", str(vtk))
        mesh(monkeypatch, "--dem", str(holed), *extra, "--out", str(ply), "--out-edges", str(edges))
        expected = vtk_fields(read_vtk(vtk.read_bytes()))
        assert ply_fields(parse_header(ply.read_bytes()).comments) == expected
        assert ply_fields(parse_header(edges.read_bytes()).comments) == expected

    def test_no_sentence_and_no_elevation_comment(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        vtk, ply = tmp_path / "x.vtk", tmp_path / "x.ply"
        mesh(monkeypatch, "--dem", str(holed), "--tolerance", "1", "--out", str(vtk))
        mesh(monkeypatch, "--dem", str(holed), "--tolerance", "1", "--out", str(ply))
        assert "elevation_source" not in read_vtk(vtk.read_bytes()).field_data
        comments = parse_header(ply.read_bytes()).comments
        assert not any(c.startswith(("elevation ", "elevation_source ")) for c in comments)

    def test_the_writers_fields_stay(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        out = tmp_path / "x.vtk"
        mesh(monkeypatch, "--dem", str(holed), "--tolerance", "1", "--out", str(out))
        assert set(read_vtk(out.read_bytes()).field_data) >= WRITER_FIELDS

    def test_a_value_is_a_bare_number_and_ascii(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """D3 rule 2: ``tolerance_m 1``, no unit in the value; rule 5, no sentence."""
        out = tmp_path / "x.vtk"
        mesh(monkeypatch, "--dem", str(holed), "--tolerance", "1", "--out", str(out))
        fields = vtk_fields(read_vtk(out.read_bytes()))
        for name in ("tolerance_m", "max_error_m", "nodata_vertices_removed"):
            float(fields[name])  # a bare number parses
        for value in fields.values():
            assert value.isascii() and not BANNED.search(value), value

    def test_stderr_uses_the_same_numbers(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """D7's summary: the triangle count, the tolerance, the largest
        difference to 5 significant figures, and the NoData count."""
        out = tmp_path / "x.vtk"
        result = mesh(monkeypatch, "--dem", str(holed), "--tolerance", "1", "--out", str(out))
        vtk = read_vtk(out.read_bytes())
        fields = vtk_fields(vtk)
        stderr = result.stderr
        assert f"{len(vtk.polygons)} triangles" in stderr, stderr
        assert f"within {fields['tolerance_m']} m" in stderr, stderr
        largest = f"{float(fields['max_error_m']):.5g}"
        assert f"largest difference {largest} m" in stderr, stderr
        removed = fields["nodata_vertices_removed"]
        assert f"{removed} vertices on NoData cells" in stderr, stderr


# ---------------------------------------------------------------- the values


class TestTheValuesAreRight:
    @pytest.mark.parametrize(("step", "count"), [(1, 25), (2, 9), (3, 4)])
    def test_nodata_vertices_removed_counts_the_strided_nodes_without_data(
        self,
        tmp_path: Path,
        holed: Path,
        monkeypatch: pytest.MonkeyPatch,
        step: int,
        count: int,
    ) -> None:
        """Without a tolerance every vertex is a DEM node at a known index,
        so the count is exact from the array: the picked nodes (every
        ``step``-th row and column, and the last) that are NoData.

        Increment 27 (`27-node-sampling.md`) ended increment 12's one-cell
        trim: a node reads only itself, so a valid node next to NoData keeps
        its height and its triangles. The literal counts are the 5 x 5 NoData
        block's nodes on each stride grid; at stride 1 that is all 25."""
        out = tmp_path / "x.vtk"
        mesh(monkeypatch, "--dem", str(holed), "--stride", str(step), "--out", str(out))
        expected, picked = picked_without_data(holed, step)
        assert expected == count  # the fixture is the one these literals describe
        vtk = read_vtk(out.read_bytes())
        assert file_field(vtk, "nodata_vertices_removed") == str(expected)
        assert len(vtk.points) == picked - expected

    @pytest.mark.parametrize("step", [1, 2, 3])
    def test_the_no_tolerance_summary_says_on_nodata_cells(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, step: int
    ) -> None:
        """Increment 27: one wording on both paths, as with ``--tolerance``
        (``TestOneVocabulary.test_stderr_uses_the_same_numbers``)."""
        out = tmp_path / "x.vtk"
        result = mesh(monkeypatch, "--dem", str(holed), "--stride", str(step), "--out", str(out))
        expected, _ = picked_without_data(holed, step)
        said = f"{expected} vertices on NoData cells were removed"
        assert said in result.stderr, result.stderr
        assert "next to" not in result.stderr, result.stderr

    def test_a_grid_without_nodata_writes_no_count(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        array = np.random.default_rng(3).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
        tif = write_tiff(tmp_path / "full.tif", micro_tiff(array))
        out, md = tmp_path / "x.vtk", tmp_path / "x.md"
        mesh(monkeypatch, "--dem", str(tif), "--stride", "2", "--out", str(out), "--stats", str(md))
        assert "nodata_vertices_removed" not in read_vtk(out.read_bytes()).field_data
        assert stats_row(md.read_text(encoding="utf-8"), "nodata_vertices_removed") == "0"

    @pytest.mark.parametrize("tolerance", ["2", "0.5"])
    def test_max_error_bounds_the_dem_nodes_inside_the_mesh(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, tolerance: str
    ) -> None:
        """At least the largest error at a DEM node inside a written triangle,
        computed here, and at least the at-vertex figure; at most the tolerance
        (25's D2 and D6). Since 15f-3 the projected path measures the DEM
        nodes within rounding of a vertex (15f's L14), so both at-vertex rows
        are in `--stats` there (15f's P1, S3)."""
        out, md = tmp_path / "x.vtk", tmp_path / "x.md"
        mesh(monkeypatch, "--dem", str(holed), "--tolerance", tolerance, "--out", str(out),
             "--stats", str(md))  # fmt: skip
        vtk = read_vtk(out.read_bytes())
        xy, z, valid = dem_nodes(holed)
        errors = located_errors(vtk, xy[valid], z[valid])
        assert np.isfinite(errors).sum() > 100
        stated = float(file_field(vtk, "max_error_m"))
        assert stated >= errors.max()
        report = md.read_text(encoding="utf-8")
        names = stats_names(report)
        assert {"dem_nodes_at_vertices", "dem_nodes_at_vertices_max_error_m"} <= set(names)
        assert stated >= float(stats_row(report, "dem_nodes_at_vertices_max_error_m"))
        assert stated <= float(file_field(vtk, "tolerance_m")) == float(tolerance)


# ---------------------------------------------------------------- the other paths


class TestTheOtherPaths:
    def test_without_a_tolerance(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        out = tmp_path / "x.vtk"
        mesh(monkeypatch, "--dem", str(holed), "--out", str(out))
        assert set(vtk_fields(read_vtk(out.read_bytes()))) == {
            "crs",
            "dem_source",
            "nodata_vertices_removed",
        }

    @pytest.mark.parametrize("suffix", [".vtk", ".ply"])
    @pytest.mark.parametrize("crs", [(), ("--crs", "EPSG:25833")], ids=["no-crs", "crs"])
    def test_flat(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, suffix: str, crs: tuple[str, ...]
    ) -> None:
        out = tmp_path / f"c{suffix}"
        mesh(monkeypatch, "catchment", "--flat", *crs, "--out", str(out))
        if suffix == ".vtk":
            fields = vtk_fields(read_vtk(out.read_bytes()))
        else:
            fields = ply_fields(parse_header(out.read_bytes()).comments)
        expected = {"heights": "none: every z is 0 (--flat)"}
        if crs:
            expected["crs"] = "EPSG:25833"
        assert fields == expected

    def test_tiles(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        source = whole(9, 13, array=terrain(9, 13))
        write_tiles(tmp_path / "tiles", quadrants(source, row_cut=4, col_cut=6, overlap=1))
        out, md = tmp_path / "m.vtk", tmp_path / "m.md"
        result = mesh(monkeypatch, "--dem", str(tmp_path / "tiles"), "--tolerance", "1",
                      "--out", str(out), "--stats", str(md))  # fmt: skip
        fields = vtk_fields(read_vtk(out.read_bytes()))
        assert fields["dem_source"] == "ne.tif; nw.tif; se.tif; sw.tif"
        assert not {"dem_tiles", "dem_seams"} & set(fields)
        report = md.read_text(encoding="utf-8")
        assert stats_row(report, "dem_tiles") == "ne.tif; nw.tif; se.tif; sw.tif"
        assert stats_row(report, "dem_seams") == "none: the tiles agree where they overlap"
        assert re.match(r"13 columns x 9 rows, ", stats_row(report, "dem_grid"))
        assert "DEM: 4 files, 13 columns x 9 rows" in result.stderr, result.stderr

    def test_features(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """The CORINE notice and the code system stay in the file; the source
        descriptions are ``--stats`` rows."""
        array = np.random.default_rng(16).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
        tif = write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))
        square = geojson(tmp_path / "square.geojson", SQUARE)
        corine = write_geojson(tmp_path / "c.geojson", [Feat(1, FOREST, {"Code_18": "311"})])
        out, md = tmp_path / "x.vtk", tmp_path / "x.md"
        mesh(monkeypatch, "--dem", str(tif), "--domain", str(square), "--tolerance", "1",
             "--features", str(corine), "--features-map", "corine",
             "--out", str(out), "--stats", str(md))  # fmt: skip
        fields = vtk_fields(read_vtk(out.read_bytes()))
        assert set(fields) == {
            "crs",
            "tolerance_m",
            "max_error_m",
            "dem_source",
            "features_notice",
            "land_cover_codes",
        }
        report = md.read_text(encoding="utf-8")
        for name in ("features", "features_crs", "features_transform", "domain", "domain_crs"):
            assert stats_row(report, name), name
        assert stats_row(report, "start_mesh") == "the domain outline and the feature lines"


# ---------------------------------------------------------------- stderr


def run_kinds(tmp_path: Path, holed: Path) -> dict[str, list[str]]:
    """Successful runs over every path that prints to stderr."""
    square = geojson(tmp_path / "square.geojson", SQUARE)
    corine = write_geojson(tmp_path / "c.geojson", [Feat(1, FOREST, {"Code_18": "311"})])
    source = whole(9, 13, array=terrain(9, 13))
    write_tiles(tmp_path / "tiles", quadrants(source, row_cut=4, col_cut=6, overlap=1))
    return {
        "refined": ["--dem", str(holed), "--tolerance", "1"],
        "strided": ["--dem", str(holed), "--stride", "2"],
        "tiles": ["--dem", str(tmp_path / "tiles"), "--tolerance", "1"],
        "features": [
            *("--dem", str(holed), "--domain", str(square), "--tolerance", "1"),
            *("--features", str(corine), "--features-map", "corine"),
        ],
        "flat": ["catchment", "--flat"],
        "reprojected": [
            *("--dem", str(VELHAS / "anadem_velhas.tif")),
            *("--domain", str(VELHAS / "catchment.geojson")),
            *("--out-crs", "EPSG:31983", "--tolerance", "5"),
        ],
    }


KINDS = ("refined", "strided", "tiles", "features", "flat", "reprojected")


class TestStderr:
    @pytest.mark.parametrize("kind", KINDS)
    def test_no_banned_word_on_a_successful_run(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, kind: str
    ) -> None:
        """D3 rule 6 over the whole of stderr (a refusal quotes the user's own
        input and is out of scope, D8)."""
        args = run_kinds(tmp_path, holed)[kind]
        result = mesh(monkeypatch, *args, "--out", str(tmp_path / "x.vtk"))
        assert not BANNED.search(result.stderr), result.stderr
        assert "elevation_source" not in result.stderr
        assert "without data dropped" not in result.stderr

    @pytest.mark.parametrize("kind", KINDS)
    def test_one_summary_line_names_the_triangles(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, kind: str
    ) -> None:
        args = run_kinds(tmp_path, holed)[kind]
        out = tmp_path / "x.vtk"
        result = mesh(monkeypatch, *args, "--out", str(out))
        triangles = len(read_vtk(out.read_bytes()).polygons)
        lines = [ln for ln in result.stderr.splitlines() if ln.startswith(f"{triangles} triangles")]
        assert len(lines) == 1, result.stderr

    def test_the_reprojected_summary_names_the_original_dem(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        args = run_kinds(tmp_path, holed)["reprojected"]
        out = tmp_path / "x.vtk"
        result = mesh(monkeypatch, *args, "--out", str(out))
        largest = f"{float(file_field(read_vtk(out.read_bytes()), 'max_error_m')):.5g}"
        assert "original DEM" in result.stderr, result.stderr
        assert f"largest difference {largest} m" in result.stderr, result.stderr


# ---------------------------------------------------------------- --record (D5)


def record_run(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, *args: str, name: str = "x"
) -> tuple[Result, Path, Path, Path]:
    """One run with ``--stats`` and ``--record`` beside the mesh."""
    out, md, rec = tmp_path / f"{name}.vtk", tmp_path / f"{name}.md", tmp_path / f"{name}.json"
    result = mesh(monkeypatch, *args, "--out", str(out), "--stats", str(md), "--record", str(rec))
    return result, out, md, rec


class TestRecord:
    @pytest.mark.parametrize("kind", KINDS)
    def test_the_shape_and_the_key_order(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, kind: str
    ) -> None:
        args = run_kinds(tmp_path, holed)[kind]
        _, _, md, rec = record_run(monkeypatch, tmp_path, *args)
        text = rec.read_text(encoding="ascii")
        assert text.endswith("\n") and not text.endswith("\n\n")
        obj = json.loads(text)
        report = md.read_text(encoding="utf-8")
        assert list(obj) == ["rasputin_version", "command", *stats_names(report)]
        assert text == json.dumps(obj, indent=1, ensure_ascii=True) + "\n"

    def test_the_version_and_the_command(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        _, out, md, rec = record_run(monkeypatch, tmp_path, "--dem", str(holed), "--tolerance", "1")
        obj = json.loads(rec.read_text(encoding="ascii"))
        assert obj["rasputin_version"] == installed_version()
        args = ["--dem", str(holed), "--tolerance", "1", "--out", str(out), "--stats", str(md),
                "--record", str(rec)]  # fmt: skip
        assert obj["command"] == shlex.join(["rasputin", "mesh", *args])

    @pytest.mark.parametrize("kind", KINDS)
    def test_values_follow_stats(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, kind: str
    ) -> None:
        """Counts are JSON integers; every ``_m`` and ``_deg`` value is a
        number equal to ``float`` of its ``--stats`` text; text values equal
        their ``--stats`` text."""
        args = run_kinds(tmp_path, holed)[kind]
        _, _, md, rec = record_run(monkeypatch, tmp_path, *args)
        obj = json.loads(rec.read_text(encoding="ascii"))
        report = md.read_text(encoding="utf-8")
        for name in stats_names(report):
            value, got = stats_row(report, name), obj[name]
            if name.endswith(("_m", "_deg")):
                assert type(got) is float and got == float(value), (name, got, value)
            elif isinstance(got, int):
                assert type(got) is int and str(got) == value, (name, got, value)
            else:
                assert got == value, (name, got, value)
        assert None not in obj.values()

    def test_the_counts_are_integers(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        _, _, _, rec = record_run(monkeypatch, tmp_path, "--dem", str(holed), "--tolerance", "1")
        obj = json.loads(rec.read_text(encoding="ascii"))
        for name in ("nodata_vertices_removed", "refinement_rounds", "points_inserted",
                     "edge_flips", "dem_nodes_outside_mesh"):  # fmt: skip
            assert type(obj[name]) is int, (name, obj[name])
        for name in ("tolerance_m", "max_error_m", "start_min_angle_deg"):
            assert type(obj[name]) is float, (name, obj[name])
        assert obj["tolerance_m"] == 1.0 and obj["snap_to_lines"] == "on"

    @pytest.mark.parametrize("kind", KINDS)
    def test_every_file_field_is_in_it(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, kind: str
    ) -> None:
        args = run_kinds(tmp_path, holed)[kind]
        _, out, _, rec = record_run(monkeypatch, tmp_path, *args)
        obj = json.loads(rec.read_text(encoding="ascii"))
        fields = vtk_fields(read_vtk(out.read_bytes()))
        assert set(fields) <= FILE_FIELDS, set(fields) - FILE_FIELDS
        for name, value in fields.items():
            got = obj[name]
            assert (got if isinstance(got, str) else float(got)) == (
                value if isinstance(got, str) else float(value)
            ), (name, got, value)
        for name in WRITER_FIELDS:
            assert name not in obj, name

    def test_twice_gives_the_same_bytes(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        args = ("--dem", str(holed), "--tolerance", "1")
        _, _, _, rec = record_run(monkeypatch, tmp_path, *args)
        first = rec.read_bytes()
        record_run(monkeypatch, tmp_path, *args)
        assert rec.read_bytes() == first

    def test_one_thread_gives_the_same_bytes(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Refinement is bit-identical for any thread count (increments 14 and
        21), and the record holds no thread count: ``cli.refine`` wrapped to
        pass ``threads=1`` (the binding's keyword; the CLI has none) changes
        no byte."""
        args = ("--dem", str(holed), "--tolerance", "1")
        _, _, _, rec = record_run(monkeypatch, tmp_path, *args)
        first = rec.read_bytes()
        real = cli.refine
        seen: list[int] = []

        def one_thread(*a: Any, **k: Any) -> Any:
            seen.append(1)
            return real(*a, **k, threads=1)

        monkeypatch.setattr(cli, "refine", one_thread)
        record_run(monkeypatch, tmp_path, *args)
        assert seen == [1], "the wrapper did not run"
        assert rec.read_bytes() == first

    def test_no_time_date_host_or_seconds(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        _, _, _, rec = record_run(monkeypatch, tmp_path, "--dem", str(holed), "--tolerance", "1")
        obj = json.loads(rec.read_text(encoding="ascii"))
        host = socket.gethostname()
        for key, value in obj.items():
            text = f"{key} {value}" if key != "command" else key
            assert "seconds" not in text.lower(), text
            assert not re.search(r"\b\d{4}-\d{2}-\d{2}\b|\b\d{1,2}:\d{2}(:\d{2})?\b", text), text
            assert not host or host not in text, text
        for absent in ("threads", "total", "stats_seconds"):
            assert absent not in obj

    def test_the_path_is_echoed_after_the_others(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        result, out, md, rec = record_run(monkeypatch, tmp_path, "--dem", str(holed))
        lines = result.stdout.splitlines()
        assert [Path(p).resolve() for p in lines] == [out.resolve(), md.resolve(), rec.resolve()]

    @pytest.mark.parametrize("kind", ["strided", "flat"])
    def test_the_no_tolerance_and_flat_paths_write_a_record(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, kind: str
    ) -> None:
        args = run_kinds(tmp_path, holed)[kind]
        _, _, _, rec = record_run(monkeypatch, tmp_path, *args)
        obj = json.loads(rec.read_text(encoding="ascii"))
        assert "tolerance_m" not in obj and "max_error_m" not in obj
        if kind == "flat":
            assert obj["heights"] == "none: every z is 0 (--flat)"

    def test_without_the_flag_nothing_is_written(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.chdir(tmp_path)
        out = tmp_path / "x.vtk"
        mesh(monkeypatch, "--dem", str(holed), "--out", str(out))
        assert sorted(p.name for p in tmp_path.iterdir()) == ["holed.tif", "x.vtk"]


class TestRecordRefusals:
    """D5: refused before any file is written; a refused run writes no record."""

    def assert_nothing_written(self, tmp_path: Path, before: set[str]) -> None:
        assert {p.name for p in tmp_path.iterdir()} == before

    def test_the_mesh_file(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        before = {p.name for p in tmp_path.iterdir()}
        out = tmp_path / "x.vtk"
        output = refused(monkeypatch, "--dem", str(holed), "--out", str(out),
                         "--record", str(out))  # fmt: skip
        assert "--record" in output, output
        self.assert_nothing_written(tmp_path, before)

    def test_the_mesh_file_by_another_spelling(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Resolved, as ``--stats`` is: ``./x.vtk`` and ``x.vtk`` are one file."""
        monkeypatch.chdir(tmp_path)
        before = {p.name for p in tmp_path.iterdir()}
        output = refused(monkeypatch, "--dem", str(holed), "--out", "x.vtk",
                         "--record", "./x.vtk")  # fmt: skip
        assert "--record" in output, output
        self.assert_nothing_written(tmp_path, before)

    def test_the_edges_file(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        before = {p.name for p in tmp_path.iterdir()}
        edges = tmp_path / "e.ply"
        output = refused(monkeypatch, "--dem", str(holed), "--out", str(tmp_path / "s.ply"),
                         "--out-edges", str(edges), "--record", str(edges))  # fmt: skip
        assert "--record" in output, output
        self.assert_nothing_written(tmp_path, before)

    def test_the_stats_file(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        before = {p.name for p in tmp_path.iterdir()}
        md = tmp_path / "x.md"
        output = refused(monkeypatch, "--dem", str(holed), "--out", str(tmp_path / "x.vtk"),
                         "--stats", str(md), "--record", str(md))  # fmt: skip
        assert "--record" in output, output
        self.assert_nothing_written(tmp_path, before)

    def test_standard_output(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """``-`` is ``--stats -``'s; ``--record -`` is refused."""
        before = {p.name for p in tmp_path.iterdir()}
        output = refused(monkeypatch, "--dem", str(holed), "--out", str(tmp_path / "x.vtk"),
                         "--record", "-")  # fmt: skip
        assert "--record" in output, output
        self.assert_nothing_written(tmp_path, before)

    def test_outside_out_parent(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        inside = tmp_path / "inside"
        inside.mkdir()
        before = {p.name for p in tmp_path.iterdir()}
        refused(monkeypatch, "--dem", str(holed), "--out-parent", str(inside),
                "--out", str(inside / "x.vtk"), "--record", str(tmp_path / "x.json"))  # fmt: skip
        self.assert_nothing_written(tmp_path, before)
        assert list(inside.iterdir()) == []

    @pytest.mark.parametrize("value", ["-1", "nan"])
    def test_a_refused_run_writes_no_record(
        self, tmp_path: Path, holed: Path, monkeypatch: pytest.MonkeyPatch, value: str
    ) -> None:
        before = {p.name for p in tmp_path.iterdir()}
        output = refused(monkeypatch, "--dem", str(holed), "--tolerance", value,
                         "--out", str(tmp_path / "x.vtk"),
                         "--record", str(tmp_path / "x.json"))  # fmt: skip
        assert "--tolerance" in output, output
        self.assert_nothing_written(tmp_path, before)

    def test_a_run_that_fails_late_writes_no_record(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A DEM with no data under any triangle passes every argument check
        and is refused after decoding: still no record."""
        array = np.full((3, 4), float(SENTINEL), dtype=np.float32)
        tif = write_tiff(tmp_path / "empty.tif", micro_tiff(array, nodata=SENTINEL))
        before = {p.name for p in tmp_path.iterdir()}
        refused(monkeypatch, "--dem", str(tif), "--tolerance", "1",
                "--out", str(tmp_path / "x.vtk"), "--record", str(tmp_path / "x.json"))  # fmt: skip
        self.assert_nothing_written(tmp_path, before)
