"""`rasputin mesh --dem`: the writers' side of the command, and the chain-mask join.

Every refusal asserts what it *says*, because `exit_code != 0` alone would
pass against a missing option.

The edge file's oracle is the surface file's faces (an edge of one triangle is
on the outline), or the features given, never the writer's own record. The
chain-mask join is `cli._triangulated`'s, which `draw` still runs on the
gallery; it is tested against an independent rebuild from the chains.

`cli._destination`'s refusals (traversal, symlinks, a missing parent) are
covered in `test_cli_draw.py`; `--out` and `--out-edges` are asserted to go
through it, and nothing more.

Not asserted, and not assertable here: that QGIS or ParaView opens either file.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

import tin_engine.cli as cli
from cli_helpers import (
    SQUARE,
    crossing_lines,
    featured_args,
    plain,
    polygon_file,
    rough_dem,
    runner,
)
from plyread import element_bytes, parse_header, read_ply
from tin_engine.cli import app
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.viz.fixtures import GALLERY

#: A gallery fixture that triangulates cleanly.
FIXTURE = "catchment"

#: A fixture whose constraints carry feature bits, so the edge file's scalar
#: has something other than 0 to carry. Increment 8 renamed it from `not-noded`.
FEATURED = "road-crosses-river"


@pytest.fixture
def attempt(request: pytest.FixtureRequest) -> cli.Attempt:
    """What the engine makes of a fixture, through the command's own pipeline."""
    name = getattr(request, "param", FIXTURE)
    return cli._triangulated(GALLERY[name], delaunay=True, spacing=cli.DEFAULT_SNAP_SPACING)


bumpy = rough_dem(16)
square = polygon_file(SQUARE)
crossing = crossing_lines()


def boundary_edges(triangles: np.ndarray) -> set[frozenset[int]]:
    """The edges of exactly one triangle: the outline of a mesh with no hole."""
    count: Counter[frozenset[int]] = Counter(
        frozenset((int(tri[e]), int(tri[(e + 1) % 3]))) for tri in triangles for e in range(3)
    )
    return {edge for edge, n in count.items() if n == 1}


@pytest.fixture
def written(tmp_path: Path, bumpy: Path) -> tuple[bytes, bytes]:
    """A DEM run to two destinations, exit code 0."""
    surface = tmp_path / "surface.ply"
    edges = tmp_path / "edges.ply"
    result = runner.invoke(
        app,
        ["mesh", "--dem", str(bumpy), "--out", str(surface), "--out-edges", str(edges)],
    )
    assert result.exit_code == 0, plain(result.output)
    return surface.read_bytes(), edges.read_bytes()


class TestTwoDestinations:
    """`rasputin mesh --dem bumpy.tif --out A --out-edges B`."""

    def test_the_surface_file_is_a_2d_mesh(self, written: tuple[bytes, bytes]) -> None:
        assert parse_header(written[0]).names == ("vertex", "face")

    def test_the_constraint_file_is_a_1d_mesh(self, written: tuple[bytes, bytes]) -> None:
        assert parse_header(written[1]).names == ("vertex", "edge")

    def test_both_files_are_binary_by_default(self, written: tuple[bytes, bytes]) -> None:
        assert {parse_header(blob).fmt for blob in written} == {"binary_little_endian"}

    def test_the_vertex_blocks_are_byte_identical(self, tmp_path: Path, bumpy: Path) -> None:
        # The two layers register on each other only because of this. Asked
        # for in binary: `element_bytes` locates fixed-width records, and byte
        # identity is a claim about the packed block.
        surface, edges = tmp_path / "s.ply", tmp_path / "e.ply"
        result = runner.invoke(
            app,
            [
                "mesh",
                "--dem",
                str(bumpy),
                "--binary",
                "--out",
                str(surface),
                "--out-edges",
                str(edges),
            ],
        )
        assert result.exit_code == 0, plain(result.output)
        blobs = surface.read_bytes(), edges.read_bytes()
        assert element_bytes(blobs[0], "vertex") == element_bytes(blobs[1], "vertex")
        assert element_bytes(blobs[0], "vertex"), "two empty slices are trivially equal"

    def test_it_names_what_it_wrote(self, tmp_path: Path, bumpy: Path) -> None:
        surface = tmp_path / "surface.ply"
        result = runner.invoke(app, ["mesh", "--dem", str(bumpy), "--out", str(surface)])
        assert result.exit_code == 0
        assert str(surface) in plain(result.output)


class TestTheEdgeFileIsTheOutline:
    """With no domain and no features, the only constraint is the grid's outline.

    The oracle is the surface file's own faces: an edge of exactly one
    triangle is on the outline. Nothing asks the edge writer what it chose.
    """

    def test_the_edges_are_the_surfaces_boundary(self, written: tuple[bytes, bytes]) -> None:
        _, surface = read_ply(written[0])
        _, edges = read_ply(written[1])
        pairs = {
            frozenset((int(a), int(b)))
            for a, b in zip(edges["edge"]["vertex1"], edges["edge"]["vertex2"], strict=True)
        }
        assert pairs == boundary_edges(surface["face"]["vertex_indices"])

    def test_each_constrained_edge_is_written_once(self, written: tuple[bytes, bytes]) -> None:
        # The two triangles sharing an interior constraint must not produce
        # two edge records; on the outline there is one triangle, so the count
        # is the outline's.
        header, _ = read_ply(written[1])
        _, surface = read_ply(written[0])
        assert header.element("edge").count == len(
            boundary_edges(surface["face"]["vertex_indices"])
        )


class TestTheTwoDestinationsAreDistinct:
    """One path given twice is one file, and the second write wins.

    Without this guard both writes succeed, the second overwrites the first,
    the command echoes two paths and exits 0. The caller has lost the surface
    they asked for with nothing to tell them.
    """

    def test_the_same_path_twice_is_refused_before_anything_is_written(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        target = tmp_path / "both.ply"
        result = CliRunner().invoke(
            app,
            ["mesh", "--dem", str(bumpy), "--out", str(target), "--out-edges", str(target)],
        )
        assert result.exit_code != 0, plain(result.output)
        assert "overwrite" in plain(result.output)
        assert not target.exists(), "refused, and yet something was written"


class TestTheChainMaskJoin:
    """The whole pair-to-mask map, built from the chains independently.

    `_chain_masks` walks the flat edge enumeration with a prefix sum, and its
    own docstring says getting that wrong "shifts every mask after the first
    ring onto the wrong edge". Nothing tested it. Measured on the mutant that
    docstring describes -- `count = len(walk) - 1` unconditionally, dropping
    every ring's closing edge -- the whole suite passed, 783 of 783, while the
    join produced `(3, 4): 2` where the truth is `1`: a road edge written to
    the file as a river, and one pair lost entirely.

    A subset assertion cannot see that. This rebuilds the expected map from
    `indices_of` and `chains` without touching `edge_properties`' layout, so it
    is an independent answer rather than the producer's own arithmetic read
    back.
    """

    @pytest.mark.parametrize("attempt", [FEATURED], indirect=True)
    def test_every_pair_carries_the_mask_its_chains_gave_it(self, attempt: cli.Attempt) -> None:
        pslg = attempt.source
        assert pslg is not None
        properties = np.asarray(pslg.edge_properties)

        expected: dict[tuple[int, int], int] = {}
        at = 0
        for c, chain in enumerate(pslg.chains):
            walk = [int(i) for i in pslg.indices_of(c)]
            closed = chain.role in cli.CORE_CLOSED_ROLES
            edges = [
                (walk[k], walk[(k + 1) % len(walk)])
                for k in range(len(walk) if closed else len(walk) - 1)
            ]
            for k, (a, b) in enumerate(edges):
                pair = (a, b) if a < b else (b, a)
                expected[pair] = expected.get(pair, 0) | int(properties[at + k])
            at += len(edges)

        assert cli._chain_masks(pslg) == expected
        assert set(expected.values()) > {0}, "the fixture must carry real bits"


class TestTheFeatureScalar:
    """The edge file carries the feature bits of the lines given, or 0."""

    def test_the_masks_come_from_the_features(
        self, tmp_path: Path, bumpy: Path, square: Path, crossing: Path
    ) -> None:
        edges = tmp_path / "edges.ply"
        result = runner.invoke(
            app,
            [
                "mesh",
                *featured_args(bumpy, square, crossing),
                "--out",
                str(tmp_path / "s.ply"),
                "--out-edges",
                str(edges),
            ],
        )
        assert result.exit_code == 0, plain(result.output)
        header, data = read_ply(edges.read_bytes())
        written_masks = set(data["edge"][header.element("edge").properties[2].name].tolist())
        bits = {p.name: 1 << p.bit for p in DEFAULT_VOCABULARY.properties}
        road, river = bits["road"], bits["river"]
        assert written_masks <= {0, road, river, road | river}, "a mask no input carries"
        assert {road, river} <= written_masks, "the road and the river reach the file"


class TestTheSecondFileIsNamedOrNotWritten:
    """Ruling 3's rejected alternative: two files always is a surprise."""

    def test_without_out_edges_only_one_file_appears(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "out" / "surface.ply"
        out.parent.mkdir()
        result = runner.invoke(app, ["mesh", "--dem", str(bumpy), "--out", str(out)])
        assert result.exit_code == 0, plain(result.output)
        assert [p.name for p in out.parent.iterdir()] == ["surface.ply"]


class TestTheAsciiFlag:
    """A small mesh under `--ascii` is inspectable with `head`."""

    def test_both_files_switch_together(self, tmp_path: Path, bumpy: Path) -> None:
        surface = tmp_path / "surface.ply"
        edges = tmp_path / "edges.ply"
        result = runner.invoke(
            app,
            [
                "mesh",
                "--dem",
                str(bumpy),
                "--ascii",
                "--out",
                str(surface),
                "--out-edges",
                str(edges),
            ],
        )
        assert result.exit_code == 0, plain(result.output)
        assert parse_header(surface.read_bytes()).fmt == "ascii"
        assert parse_header(edges.read_bytes()).fmt == "ascii"
