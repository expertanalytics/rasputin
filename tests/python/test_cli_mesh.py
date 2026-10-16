"""`rasputin mesh --dem`: the writers' side of the command, the chain-mask join,
`--snap-spacing`, and the role map.

Every refusal asserts what it *says*, because `exit_code != 0` alone would
pass against a missing option.

The edge file's oracle is the surface file's faces (an edge of one triangle is
on the outline), or the features given, never the writer's own record. The
chain-mask join is checked on a road crossing a river, built here and run
through `cli._engine`, against an independent rebuild from the chains.

`cli_kit.destination`'s refusals (traversal, symlinks, a missing parent) are
covered in `test_cli_kit.py`; `--out` and `--out-edges` are asserted to go
through it, and nothing more.

Not asserted, and not assertable here: that QGIS or ParaView opens either file.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

import tin_engine._core as core
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

#: The road-crosses-river shape: a square outline with a road and a river
#: crossing at its centre, at UTM 33N magnitudes so the default snap spacing
#: (1 mm) is the one a real run uses. The noder adds the crossing at
#: (350, 350); the chains carry feature bits 2 (road) and 1 (river).
ORIGIN = np.array([430_000.0, 6_900_000.0])
CROSSING_VERTICES = ORIGIN + np.array(
    [(0, 0), (700, 0), (700, 700), (0, 700), (100, 100), (600, 600), (100, 600), (600, 100)],
    dtype=float,
)
CROSSING_CHAINS = [
    ([0, 1, 2, 3], core.ChainRole.Outer, 0),
    ([4, 5], core.ChainRole.Breakline, 2),
    ([6, 7], core.ChainRole.Breakline, 1),
]

#: Fine enough that a UTM easting overflows the snap lattice (the largest grid
#: index is 2**51, so the representable coordinate at spacing s is about
#: s * 2.25e15): only the value given on the command line can make the noder
#: answer `CoordinateOutOfRange`.
OVERFLOWING_SNAP_SPACING = 1e-12


@pytest.fixture
def crossing_graph() -> core.NodedPslg:
    """The road-crosses-river shape, noded and triangulated by the command's engine."""
    run = cli._engine(
        CROSSING_VERTICES, CROSSING_CHAINS, delaunay=True, spacing=cli.DEFAULT_SNAP_SPACING
    )
    assert run.ok, (run.status, run.message)
    assert run.noded is not None
    assert len(np.asarray(run.noded.vertices)) == 9, "the crossing was not noded"
    return run.noded


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

    def test_every_pair_carries_the_mask_its_chains_gave_it(
        self, crossing_graph: core.NodedPslg
    ) -> None:
        pslg = crossing_graph
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
        assert set(expected.values()) == {0, 1, 2}, "the shape must carry both bits"


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


class TestSnapSpacing:
    """Two refusal channels for `--snap-spacing`, neither removable for the other.

    The command refuses a non-finite or non-positive value as a usage error
    before anything is read. A value it accepts reaches the noder, whose own
    refusal comes back in the engine's words.
    """

    @pytest.mark.parametrize("value", ["0", "-1", "nan", "-inf"])
    def test_a_non_positive_or_non_finite_spacing_is_a_usage_error(
        self, value: str, tmp_path: Path, bumpy: Path
    ) -> None:
        target = tmp_path / "refused.ply"
        result = runner.invoke(
            app, ["mesh", "--dem", str(bumpy), f"--snap-spacing={value}", "--out", str(target)]
        )
        message = plain(result.output)
        assert result.exit_code == 2, message
        assert "snap-spacing" in message, f"the refusal does not name the option: {message!r}"
        assert not target.exists(), "a usage error wrote a file"

    def test_the_spacing_reaches_the_engine(self, tmp_path: Path, bumpy: Path) -> None:
        # A command that accepted the option and passed the default on would
        # mesh this DEM; only the real argument makes the noder refuse.
        target = tmp_path / "overflow.ply"
        result = runner.invoke(
            app,
            [
                "mesh",
                "--dem",
                str(bumpy),
                f"--snap-spacing={OVERFLOWING_SNAP_SPACING}",
                "--out",
                str(target),
            ],
        )
        message = plain(result.output)
        assert result.exit_code != 0, message
        assert "CoordinateOutOfRange" in message
        assert core.describe(core.NodeStatus.CoordinateOutOfRange) in message
        assert not target.exists()


class TestTheRoleMap:
    def test_the_cli_maps_every_role_the_enum_has_and_nothing_else(self) -> None:
        # Total both ways: a role the enum gains and `ROLES` lacks is an input
        # vocabulary that silently cannot express it. `__members__`, because a
        # pybind11 enum type is not iterable.
        assert set(cli.ROLES.values()) == set(core.ChainRole.__members__.values())
        assert len(cli.ROLES) == len(core.ChainRole.__members__)
