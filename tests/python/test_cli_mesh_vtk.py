"""`rasputin mesh --out x.vtk`: increment 13's half of the command.

`13-bundled-mesh.md` ruling 8: the suffix of `--out` picks the format, an
unknown suffix is refused naming both valid ones, `--out-edges` beside a `.vtk`
is refused because the edges are already in the file, and `--binary/--ascii`
is one flag pair defaulting to text for both formats (U2 (a)).

The `.vtk` bytes are parsed with `vtkread`, not with VTK; the read-back through
VTK's own readers is `test_io_vtk_readback.py`. The oracle for the lines is the
mesh's own `constrained_edges` bits, walked by `constrained_edge_set` in
`test_cli_mesh.py`, and the oracle for their masks is `cli._chain_masks`, which
that file tests against an independent rebuild. Nothing here asks the writer
what it wrote.

Every refusal asserts what it *says* and that no file was left behind,
because `exit_code != 0` alone passes against a missing option.

`TestCrsRefusals` passes before increment 13 and is a guard, not red: the
`--crs` usage error already fires before any write, whatever the suffix. It
pins that the new dispatch keeps it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from typer.testing import CliRunner

import tin_engine.cli as cli
from plyread import parse_header
from test_cli_mesh import constrained_edge_set, plain
from tin_engine.cli import app
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.viz.fixtures import GALLERY
from vtkread import VtkFile, lines_as_array, polygons_as_array, read_vtk

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

FIXTURE = "catchment"
FEATURED = "road-crosses-river"
ELEVATION_TEXT = "none (z=0, --flat)"


def invoke(*args: str) -> tuple[int, str]:
    result = runner.invoke(app, ["mesh", *args])
    return result.exit_code, plain(result.output)


def attempt_for(name: str) -> cli.Attempt:
    return cli._triangulated(GALLERY[name], delaunay=True, spacing=cli.DEFAULT_SNAP_SPACING)


@pytest.fixture
def featured(tmp_path: Path) -> VtkFile:
    out = tmp_path / "rr.vtk"
    code, output = invoke(FEATURED, "--flat", "--crs", "EPSG:25833", "--out", str(out))
    assert code == 0, output
    return read_vtk(out.read_bytes())


class TestTheSuffixPicksTheFormat:
    """Ruling 8. There is no `--format`: it could contradict the suffix."""

    def test_vtk_writes_one_vtk_file(self, tmp_path: Path) -> None:
        out = tmp_path / "mesh.vtk"
        code, output = invoke(FIXTURE, "--flat", "--out", str(out))
        assert code == 0, output
        assert [p.name for p in tmp_path.iterdir()] == ["mesh.vtk"]
        assert out.read_bytes().startswith(b"# vtk DataFile Version 4.2\n")

    def test_it_names_what_it_wrote(self, tmp_path: Path) -> None:
        out = tmp_path / "mesh.vtk"
        code, output = invoke(FIXTURE, "--flat", "--out", str(out))
        assert code == 0
        assert str(out) in output

    def test_ply_still_writes_ply(self, tmp_path: Path) -> None:
        out = tmp_path / "mesh.ply"
        code, output = invoke(FIXTURE, "--flat", "--out", str(out))
        assert code == 0, output
        assert out.read_bytes().startswith(b"ply\n")

    @pytest.mark.parametrize("name", ["mesh.vtp", "mesh.vtu", "mesh.txt", "mesh", "mesh.ply.bak"])
    def test_an_unknown_suffix_is_refused_naming_both(self, tmp_path: Path, name: str) -> None:
        out = tmp_path / name
        code, output = invoke(FIXTURE, "--flat", "--out", str(out))
        assert code != 0, output
        assert ".vtk" in output and ".ply" in output, output
        assert not out.exists()
        assert list(tmp_path.iterdir()) == [], "refused, and yet something was written"


class TestOutEdgesWithVtk:
    """Ruling 8: the edges are already in the file, so a second one is refused."""

    def test_it_is_refused_before_anything_is_written(self, tmp_path: Path) -> None:
        out, edges = tmp_path / "mesh.vtk", tmp_path / "edges.ply"
        code, output = invoke(FIXTURE, "--flat", "--out", str(out), "--out-edges", str(edges))
        assert code != 0, output
        assert "--out-edges" in output
        assert list(tmp_path.iterdir()) == []


class TestEncoding:
    """`--binary/--ascii`, one pair for both formats, text by default (U2 (a))."""

    def test_vtk_is_ascii_by_default(self, tmp_path: Path) -> None:
        out = tmp_path / "mesh.vtk"
        assert invoke(FIXTURE, "--flat", "--out", str(out))[0] == 0
        assert read_vtk(out.read_bytes()).encoding == "ASCII"

    def test_binary_writes_the_binary_header(self, tmp_path: Path) -> None:
        out = tmp_path / "mesh.vtk"
        code, output = invoke(FIXTURE, "--flat", "--binary", "--out", str(out))
        assert code == 0, output
        assert read_vtk(out.read_bytes()).encoding == "BINARY"

    def test_ascii_is_the_explicit_spelling(self, tmp_path: Path) -> None:
        out = tmp_path / "mesh.vtk"
        code, output = invoke(FIXTURE, "--flat", "--ascii", "--out", str(out))
        assert code == 0, output
        assert read_vtk(out.read_bytes()).encoding == "ASCII"

    def test_ply_is_ascii_by_default(self, tmp_path: Path) -> None:
        surface, edges = tmp_path / "s.ply", tmp_path / "e.ply"
        code, output = invoke(FIXTURE, "--flat", "--out", str(surface), "--out-edges", str(edges))
        assert code == 0, output
        assert parse_header(surface.read_bytes()).fmt == "ascii"
        assert parse_header(edges.read_bytes()).fmt == "ascii"

    def test_binary_switches_both_ply_files(self, tmp_path: Path) -> None:
        surface, edges = tmp_path / "s.ply", tmp_path / "e.ply"
        code, output = invoke(
            FIXTURE, "--flat", "--binary", "--out", str(surface), "--out-edges", str(edges)
        )
        assert code == 0, output
        assert parse_header(surface.read_bytes()).fmt == "binary_little_endian"
        assert parse_header(edges.read_bytes()).fmt == "binary_little_endian"


class TestTheBundleIsTheEnginesMesh:
    """`road-crosses-river`: 9 points, 12 triangles, and its constraint lines."""

    def test_the_counts(self, featured: VtkFile) -> None:
        assert len(featured.points) == 9
        assert len(featured.polygons) == 12
        assert featured.cell_count == len(featured.lines) + 12

    def test_the_points_are_the_meshs_own_with_z_0(self, featured: VtkFile) -> None:
        mesh = attempt_for(FEATURED).mesh
        assert mesh is not None
        assert_array_equal(featured.points[:, :2], np.asarray(mesh.vertices))
        assert_array_equal(featured.points[:, 2], np.zeros(len(featured.points)))

    def test_the_polygons_are_the_meshs_triangles(self, featured: VtkFile) -> None:
        mesh = attempt_for(FEATURED).mesh
        assert mesh is not None
        assert_array_equal(polygons_as_array(featured), np.asarray(mesh.triangles))

    def test_the_lines_are_the_constrained_edges_once_each(self, featured: VtkFile) -> None:
        mesh = attempt_for(FEATURED).mesh
        lines = lines_as_array(featured)
        pairs = [frozenset(map(int, line)) for line in lines]
        assert len(pairs) == len(set(pairs)), "an interior constraint written twice"
        assert set(pairs) == constrained_edge_set(mesh)

    def test_each_line_carries_the_mask_its_chains_gave_it(self, featured: VtkFile) -> None:
        attempt = attempt_for(FEATURED)
        joined = cli._chain_masks(attempt.source)  # type: ignore[arg-type]
        lines = lines_as_array(featured)
        expected = [joined.get(cli._undirected(int(a), int(b)), 0) for a, b in lines]
        mask = featured.cell_array("feature_mask").values.tolist()
        assert mask[: len(lines)] == expected
        assert mask[len(lines) :] == [0] * len(featured.polygons)
        assert set(expected) > {0}, "the fixture must carry real bits"

    def test_river_and_road_are_the_per_feature_arrays(self, featured: VtkFile) -> None:
        assert set(featured.cell_fields["features"]) == {"river", "road"}

    def test_the_vocabulary_is_the_default_one(self, featured: VtkFile) -> None:
        pairs = sorted((p.bit, p.name) for p in DEFAULT_VOCABULARY.properties)
        assert featured.field_data["feature_bits"].values.tolist() == [b for b, _ in pairs]
        assert featured.field_data["feature_names"].values == tuple(n for _, n in pairs)
        assert featured.field_data["feature_vocabulary"].values == (
            DEFAULT_VOCABULARY.fingerprint(),
        )

    def test_crs_and_elevation_are_carried(self, featured: VtkFile) -> None:
        assert featured.field_data["crs"].values == ("EPSG:25833",)
        assert featured.field_data["elevation"].values == (ELEVATION_TEXT,)

    def test_no_crs_means_no_crs_field(self, tmp_path: Path) -> None:
        out = tmp_path / "mesh.vtk"
        assert invoke(FIXTURE, "--flat", "--out", str(out))[0] == 0
        assert "crs" not in read_vtk(out.read_bytes()).field_data

    def test_a_fixture_with_no_bits_writes_no_features_block(self, tmp_path: Path) -> None:
        # catchment's constraints are all unclassified (mask 0).
        out = tmp_path / "mesh.vtk"
        assert invoke(FIXTURE, "--flat", "--out", str(out))[0] == 0
        assert read_vtk(out.read_bytes()).cell_fields == {}


class TestCrsRefusals:
    """`--crs` refusals surface as usage errors for `.vtk` too (ruling 7).

    Green before increment 13: the command already checks `--crs` before any
    write. Kept as guards on the new dispatch.
    """

    @pytest.mark.parametrize("crs", ["EPSG:25833\rforged", "ETRS89 60\N{DEGREE SIGN}N"])
    def test_a_bad_crs_is_a_usage_error_and_no_file(self, tmp_path: Path, crs: str) -> None:
        out = tmp_path / "mesh.vtk"
        code, output = invoke(FIXTURE, "--flat", "--crs", crs, "--out", str(out))
        assert code == 2, output
        assert "--crs" in output
        assert not out.exists()


class TestThePlyEdgeFileCarriesTheVocabulary:
    """Ruling 9 under U1 (a): the shipped edge file said nothing about its bits."""

    @pytest.fixture
    def comments(self, tmp_path: Path) -> tuple[str, ...]:
        surface, edges = tmp_path / "s.ply", tmp_path / "e.ply"
        code, output = invoke(FEATURED, "--flat", "--out", str(surface), "--out-edges", str(edges))
        assert code == 0, output
        return parse_header(edges.read_bytes()).comments

    def test_every_bit_is_named(self, comments: tuple[str, ...]) -> None:
        expected = [
            f"feature_bit {bit} {name}"
            for bit, name in sorted((p.bit, p.name) for p in DEFAULT_VOCABULARY.properties)
        ]
        assert [c for c in comments if c.startswith("feature_bit ")] == expected

    def test_the_fingerprint_is_carried(self, comments: tuple[str, ...]) -> None:
        assert f"feature_vocabulary {DEFAULT_VOCABULARY.fingerprint()}" in comments
