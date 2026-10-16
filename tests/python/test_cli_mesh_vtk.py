"""`rasputin mesh --out x.vtk`: the suffix picks the format, and what the file carries.

The suffix of `--out` picks the format, an unknown suffix is refused naming
both valid ones, `--out-edges` beside a `.vtk` is refused because the edges are
already in the file, and `--binary/--ascii` is one flag pair for both formats,
binary by default.

The `.vtk` bytes are parsed with `vtkread`, not with VTK; the read-back through
VTK's own readers is `test_io_vtk_readback.py`. The oracles are the file's own
triangles (a line is an edge of one, written once) and the features given (a
mask carries only their bits). Nothing here asks the writer what it wrote.

Every refusal asserts what it *says* and that no file was left behind,
because `exit_code != 0` alone passes against a missing option.
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import pytest

from cli_helpers import (
    SQUARE,
    crossing_lines,
    featured_args,
    help_text,
    invoke,
    polygon_file,
    rough_dem,
)
from plyread import parse_header
from tin_engine.features import DEFAULT_VOCABULARY
from vtkread import VtkFile, lines_as_array, polygons_as_array, read_vtk

bumpy = rough_dem(16)
square = polygon_file(SQUARE)
crossing = crossing_lines()

#: The bit of each feature property, from the vocabulary rather than restated.
BITS = {p.name: 1 << p.bit for p in DEFAULT_VOCABULARY.properties}


def triangle_edges(vtk: VtkFile) -> Counter[frozenset[int]]:
    """Every edge of the file's triangles, with the number of triangles it bounds."""
    return Counter(
        frozenset((int(tri[e]), int(tri[(e + 1) % 3])))
        for tri in polygons_as_array(vtk)
        for e in range(3)
    )


@pytest.fixture
def featured(tmp_path: Path, bumpy: Path, square: Path, crossing: Path) -> VtkFile:
    out = tmp_path / "rr.vtk"
    code, output = invoke("mesh", *featured_args(bumpy, square, crossing), "--out", str(out))
    assert code == 0, output
    return read_vtk(out.read_bytes())


class TestTheSuffixPicksTheFormat:
    """There is no `--format`: it could contradict the suffix."""

    def test_vtk_writes_one_vtk_file(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "mesh.vtk"
        code, output = invoke("mesh", "--dem", str(bumpy), "--out", str(out))
        assert code == 0, output
        assert sorted(p.name for p in tmp_path.iterdir()) == [bumpy.name, "mesh.vtk"]
        assert out.read_bytes().startswith(b"# vtk DataFile Version 4.2\n")

    def test_it_names_what_it_wrote(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "mesh.vtk"
        code, output = invoke("mesh", "--dem", str(bumpy), "--out", str(out))
        assert code == 0
        assert str(out) in output

    def test_ply_still_writes_ply(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "mesh.ply"
        code, output = invoke("mesh", "--dem", str(bumpy), "--out", str(out))
        assert code == 0, output
        assert out.read_bytes().startswith(b"ply\n")

    @pytest.mark.parametrize("name", ["mesh.vtp", "mesh.vtu", "mesh.txt", "mesh", "mesh.ply.bak"])
    def test_an_unknown_suffix_is_refused_naming_both(
        self, tmp_path: Path, bumpy: Path, name: str
    ) -> None:
        out = tmp_path / name
        code, output = invoke("mesh", "--dem", str(bumpy), "--out", str(out))
        assert code != 0, output
        assert ".vtk" in output and ".ply" in output, output
        assert not out.exists()
        assert list(tmp_path.iterdir()) == [bumpy], "refused, and yet something was written"


class TestOutEdgesWithVtk:
    """The edges are already in the file, so a second one is refused."""

    def test_it_is_refused_before_anything_is_written(self, tmp_path: Path, bumpy: Path) -> None:
        out, edges = tmp_path / "mesh.vtk", tmp_path / "edges.ply"
        code, output = invoke(
            "mesh", "--dem", str(bumpy), "--out", str(out), "--out-edges", str(edges)
        )
        assert code != 0, output
        assert "--out-edges" in output
        assert list(tmp_path.iterdir()) == [bumpy]


class TestEncoding:
    """`--binary/--ascii`, one pair for both formats, binary by default (increment 31)."""

    def test_vtk_is_binary_by_default(self, tmp_path: Path, bumpy: Path) -> None:
        # Increment 31 (`31-binary-default.md`, D1): binary unless `--ascii`.
        # This was `test_vtk_is_ascii_by_default` under increment 13's U2 (a).
        out = tmp_path / "mesh.vtk"
        code, output = invoke("mesh", "--dem", str(bumpy), "--out", str(out))
        assert code == 0, output
        blob = out.read_bytes()
        assert blob.split(b"\n")[2] == b"BINARY"
        assert read_vtk(blob).encoding == "BINARY"

    def test_binary_writes_the_binary_header(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "mesh.vtk"
        code, output = invoke("mesh", "--dem", str(bumpy), "--binary", "--out", str(out))
        assert code == 0, output
        assert read_vtk(out.read_bytes()).encoding == "BINARY"

    def test_ascii_is_the_explicit_spelling(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "mesh.vtk"
        code, output = invoke("mesh", "--dem", str(bumpy), "--ascii", "--out", str(out))
        assert code == 0, output
        assert read_vtk(out.read_bytes()).encoding == "ASCII"

    def test_ply_is_binary_by_default(self, tmp_path: Path, bumpy: Path) -> None:
        # Increment 31, D2: the pair's one default covers `.ply` too.
        # This was `test_ply_is_ascii_by_default` under increment 13's U2 (a).
        surface, edges = tmp_path / "s.ply", tmp_path / "e.ply"
        code, output = invoke(
            "mesh", "--dem", str(bumpy), "--out", str(surface), "--out-edges", str(edges)
        )
        assert code == 0, output
        assert parse_header(surface.read_bytes()).fmt == "binary_little_endian"
        assert parse_header(edges.read_bytes()).fmt == "binary_little_endian"

    def test_the_help_shows_binary_as_the_default(self) -> None:
        # Increment 31, D3: Typer prints the default itself, so the help and
        # the behaviour cannot drift apart. The marker is read from the
        # `--binary --ascii` row onward, so another option's marker cannot
        # stand in for it. The width is pinned (`help_text`, on CI too) so the
        # caller's terminal cannot truncate the option names; `help_text` undoes
        # Rich's wrapping.
        text = help_text("mesh", columns=100)
        marker = re.search(r"--binary --ascii\b.*?\[default: ([^\]]+)\]", text)
        assert marker is not None, text
        assert marker.group(1) == "binary", text
        # The command's own description is part of the same help page, and
        # until increment 31 said so in words
        # (`src_python/tin_engine/cli.py@2764ef71:755`).
        assert "Text is the default" not in text, text

    def test_binary_switches_both_ply_files(self, tmp_path: Path, bumpy: Path) -> None:
        surface, edges = tmp_path / "s.ply", tmp_path / "e.ply"
        code, output = invoke(
            "mesh",
            "--dem",
            str(bumpy),
            "--binary",
            "--out",
            str(surface),
            "--out-edges",
            str(edges),
        )
        assert code == 0, output
        assert parse_header(surface.read_bytes()).fmt == "binary_little_endian"
        assert parse_header(edges.read_bytes()).fmt == "binary_little_endian"


class TestTheBundle:
    """A road crossing a river inside a square: triangles, lines and their bits."""

    def test_the_counts(self, featured: VtkFile) -> None:
        assert featured.cell_count == len(featured.lines) + len(featured.polygons)
        assert len(featured.polygons) > 0 and len(featured.lines) > 0

    def test_every_line_is_a_triangle_edge_written_once(self, featured: VtkFile) -> None:
        pairs = [frozenset(map(int, line)) for line in lines_as_array(featured)]
        assert len(pairs) == len(set(pairs)), "a constraint written twice"
        assert set(pairs) <= set(triangle_edges(featured)), "a line no triangle has"

    def test_the_outline_is_among_the_lines(self, featured: VtkFile) -> None:
        outline = {edge for edge, n in triangle_edges(featured).items() if n == 1}
        assert outline <= {frozenset(map(int, line)) for line in lines_as_array(featured)}

    def test_each_line_carries_only_the_bits_given(self, featured: VtkFile) -> None:
        road, river = BITS["road"], BITS["river"]
        lines = lines_as_array(featured)
        mask = featured.cell_array("feature_mask").values.tolist()
        assert set(mask[: len(lines)]) <= {0, road, river, road | river}
        assert {road, river} <= set(mask[: len(lines)]), "the road and the river reach the file"
        assert mask[len(lines) :] == [0] * len(featured.polygons)

    def test_river_and_road_are_the_per_feature_arrays(self, featured: VtkFile) -> None:
        assert set(featured.cell_fields["features"]) == {"river", "road"}

    def test_the_vocabulary_is_the_default_one(self, featured: VtkFile) -> None:
        pairs = sorted((p.bit, p.name) for p in DEFAULT_VOCABULARY.properties)
        assert featured.field_data["feature_bits"].values.tolist() == [b for b, _ in pairs]
        assert featured.field_data["feature_names"].values == tuple(n for _, n in pairs)
        assert featured.field_data["feature_vocabulary"].values == (
            DEFAULT_VOCABULARY.fingerprint(),
        )

    def test_the_dems_crs_is_carried(self, featured: VtkFile) -> None:
        assert featured.field_data["crs"].values == ("EPSG:25833",)

    def test_a_run_with_no_bits_writes_no_features_block(self, tmp_path: Path, bumpy: Path) -> None:
        # A plain DEM run's only constraint is the outline, unclassified (mask 0).
        out = tmp_path / "mesh.vtk"
        assert invoke("mesh", "--dem", str(bumpy), "--out", str(out))[0] == 0
        assert "features" not in read_vtk(out.read_bytes()).cell_fields


class TestThePlyEdgeFileCarriesTheVocabulary:
    """A `.ply` edge file names its bits and the vocabulary they come from."""

    @pytest.fixture
    def comments(
        self, tmp_path: Path, bumpy: Path, square: Path, crossing: Path
    ) -> tuple[str, ...]:
        surface, edges = tmp_path / "s.ply", tmp_path / "e.ply"
        code, output = invoke(
            "mesh",
            *featured_args(bumpy, square, crossing),
            *("--out", str(surface), "--out-edges", str(edges)),
        )
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
