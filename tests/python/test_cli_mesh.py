"""`rasputin mesh`: increment 10's command, committed RED.

`tin_engine.cli` already exists, so the intended failure here is Click's
`No such command 'mesh'` -- exit code 2, no file written. Every test that
asserts a refusal also asserts what the refusal *says*, because an assertion of
the form `exit_code != 0` alone would pass against the missing command and
measure nothing (the same trap `test_cli_draw.py` documents).

This is the composition-root half of the increment. `io/ply.py` is pure and has
no path; this command is the only thing that has one, which keeps
`06-cdt-viewer.md`'s "no file is written below `cli.py`" true for the mesh path
as it already is for the SVG path.

The oracle for the edge file is built from the **mesh**, not from the writer:
`IndexedMesh2.constrained_edges` sets bit `e` of a triangle's mask iff the edge
`(v[e], v[(e + 1) % 3])` is constrained, per `_core.pyi`, and
`constrained_edge_set` below walks that convention directly. A test that asked
the writer which edges it had chosen would agree with any internally consistent
answer, including a wrong one (`GLOSSARY.md`, "self-confirming oracle").

`cli._destination`'s refusals -- traversal, symlinks, a missing parent -- are
increment 6b-ii's boundary and are covered in `test_cli_draw.py`. They are not
re-tested here; `--out` and `--out-edges` are asserted to *go through* it, and
nothing more.

NOT ASSERTED, AND NOT ASSERTABLE HERE: that QGIS or ParaView opens either file.
No runner in this repository has them. Ruling 2 puts that in acceptance, to be
done once by a person who then writes the version down.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from typer.testing import CliRunner

import tin_engine.cli as cli
from plyread import element_bytes, parse_header, read_ply, vertex_array
from tin_engine.cli import app
from tin_engine.viz.fixtures import GALLERY

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

ANSI = re.compile(r"\x1b\[[0-9;]*m")
BOX = re.compile(r"[─-╿]")


def plain(text: str) -> str:
    """Output as a reader sees it, per `test_cli_draw.py`'s `plain`.

    Typer renders a refusal inside a Rich panel whose width is the terminal's,
    so a message asserted verbatim is flaky by construction: it passes on a
    wide terminal and fails on a narrow one. Collapsing the box rule and the
    whitespace leaves the words and numbers every assertion here looks at.
    """
    return " ".join(BOX.sub(" ", ANSI.sub("", text)).split())

#: The design's acceptance line uses this one. It triangulates cleanly.
FIXTURE = "catchment"

#: A fixture whose constraints carry feature bits, so the edge file's scalar
#: has something other than 0 to carry. Increment 8 renamed it from `not-noded`.
FEATURED = "road-crosses-river"


def constrained_edge_set(mesh: object) -> set[frozenset[int]]:
    """Every constrained edge of a mesh, deduplicated across shared triangles.

    Built from `_core.pyi`'s stated bit convention, independently of anything
    increment 10 adds. The two triangles either side of an interior constraint
    both flag it, which is why this is a set of unordered pairs rather than a
    list.
    """
    triangles = np.asarray(mesh.triangles)  # type: ignore[attr-defined]
    masks = np.asarray(mesh.constrained_edges)  # type: ignore[attr-defined]
    return {
        frozenset((int(tri[e]), int(tri[(e + 1) % 3])))
        for tri, mask in zip(triangles, masks, strict=True)
        for e in range(3)
        if mask & (1 << e)
    }


@pytest.fixture
def attempt(request: pytest.FixtureRequest) -> cli.Attempt:
    """What the engine makes of a fixture, through the command's own pipeline."""
    name = getattr(request, "param", FIXTURE)
    return cli._triangulated(GALLERY[name], delaunay=True, spacing=cli.DEFAULT_SNAP_SPACING)


@pytest.fixture
def written(tmp_path: Path) -> tuple[bytes, bytes]:
    """The acceptance invocation: one fixture, two destinations, exit code 0."""
    surface = tmp_path / "surface.ply"
    edges = tmp_path / "edges.ply"
    result = runner.invoke(
        app,
        ["mesh", FIXTURE, "--flat", "--out", str(surface), "--out-edges", str(edges)],
    )
    assert result.exit_code == 0, plain(result.output)
    return surface.read_bytes(), edges.read_bytes()


class TestTheAcceptanceInvocation:
    """`rasputin mesh catchment --flat --out A --out-edges B`."""

    def test_the_surface_file_is_a_2d_mesh(self, written: tuple[bytes, bytes]) -> None:
        assert parse_header(written[0]).names == ("vertex", "face")

    def test_the_constraint_file_is_a_1d_mesh(self, written: tuple[bytes, bytes]) -> None:
        assert parse_header(written[1]).names == ("vertex", "edge")

    def test_both_files_are_binary_by_default(self, written: tuple[bytes, bytes]) -> None:
        assert {parse_header(blob).fmt for blob in written} == {"binary_little_endian"}

    def test_the_vertex_blocks_are_byte_identical(self, written: tuple[bytes, bytes]) -> None:
        # Ruling 3. The two layers register on each other only because of this.
        assert element_bytes(written[0], "vertex") == element_bytes(written[1], "vertex")

    def test_it_names_what_it_wrote(self, tmp_path: Path) -> None:
        surface = tmp_path / "surface.ply"
        result = runner.invoke(app, ["mesh", FIXTURE, "--flat", "--out", str(surface)])
        assert result.exit_code == 0
        assert str(surface) in plain(result.output)


class TestTheGeometryIsTheEngines:
    """The file holds the mesh the engine produced, vertex for vertex."""

    def test_the_vertices_are_the_meshs_own(
        self, written: tuple[bytes, bytes], attempt: cli.Attempt
    ) -> None:
        assert attempt.mesh is not None
        _, data = read_ply(written[0])
        assert_array_equal(vertex_array(data)[:, :2], np.asarray(attempt.mesh.vertices))

    def test_the_faces_are_the_meshs_own(
        self, written: tuple[bytes, bytes], attempt: cli.Attempt
    ) -> None:
        assert attempt.mesh is not None
        _, data = read_ply(written[0])
        assert_array_equal(data["face"]["vertex_indices"], np.asarray(attempt.mesh.triangles))

    def test_the_edges_are_the_meshs_constrained_edges(
        self, written: tuple[bytes, bytes], attempt: cli.Attempt
    ) -> None:
        assert attempt.mesh is not None
        _, data = read_ply(written[1])
        pairs = {
            frozenset((int(a), int(b)))
            for a, b in zip(data["edge"]["vertex1"], data["edge"]["vertex2"], strict=True)
        }
        assert pairs == constrained_edge_set(attempt.mesh)

    def test_each_constrained_edge_is_written_once(
        self, written: tuple[bytes, bytes], attempt: cli.Attempt
    ) -> None:
        # The dedup ruling 6 names: the two triangles sharing an interior
        # constraint must not produce two edge records.
        header, _ = read_ply(written[1])
        assert attempt.mesh is not None
        assert header.element("edge").count == len(constrained_edge_set(attempt.mesh))


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
    def test_every_pair_carries_the_mask_its_chains_gave_it(
        self, attempt: cli.Attempt
    ) -> None:
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
    """The edge file carries the noded PSLG's feature bits, or 0."""

    @pytest.mark.parametrize("attempt", [FEATURED], indirect=True)
    def test_the_masks_come_from_the_noded_graph(
        self, tmp_path: Path, attempt: cli.Attempt
    ) -> None:
        edges = tmp_path / "edges.ply"
        result = runner.invoke(
            app, ["mesh", FEATURED, "--flat", "--out", str(tmp_path / "s.ply"),
                  "--out-edges", str(edges)]
        )
        assert result.exit_code == 0, plain(result.output)
        header, data = read_ply(edges.read_bytes())
        written_masks = set(data["edge"][header.element("edge").properties[2].name].tolist())
        available = set(np.asarray(attempt.source.edge_properties).tolist())  # type: ignore[union-attr]
        assert written_masks <= available | {0}, "a mask nothing in the input carries"
        assert written_masks - {0}, f"{FEATURED} carries feature bits; none reached the file"


class TestFlatIsAlwaysAWordSomebodyTyped:
    """Ruling 4. Nothing in the tree can supply an elevation yet."""

    def test_without_flat_it_refuses(self, tmp_path: Path) -> None:
        out = tmp_path / "surface.ply"
        result = runner.invoke(app, ["mesh", FIXTURE, "--out", str(out)])
        assert result.exit_code != 0
        assert "--flat" in plain(result.output)
        assert not out.exists(), "a refusal must not leave a file behind"

    def test_flat_fills_z_with_zero(self, written: tuple[bytes, bytes]) -> None:
        _, data = read_ply(written[0])
        assert_array_equal(data["vertex"]["z"], np.zeros(len(data["vertex"]["z"])))

    def test_flat_says_so_in_both_headers(self, written: tuple[bytes, bytes]) -> None:
        # Verbatim from ruling 4, because the comment is the only thing that
        # tells a person six months later that the surface is not terrain.
        for blob in written:
            assert "elevation none (z=0, --flat)" in parse_header(blob).comments


class TestTheCrsComment:
    """Ruling 5: the CRS is a comment, and no reader acts on it."""

    def test_the_crs_reaches_the_header(self, tmp_path: Path) -> None:
        out = tmp_path / "surface.ply"
        result = runner.invoke(
            app, ["mesh", FIXTURE, "--flat", "--crs", "EPSG:25833", "--out", str(out)]
        )
        assert result.exit_code == 0, plain(result.output)
        assert "crs EPSG:25833" in parse_header(out.read_bytes()).comments

    def test_it_is_not_validated(self, tmp_path: Path) -> None:
        # The writer does not import pyproj and the command does not check the
        # string: CRS authority is `raster.py`'s and `io/`'s, and a second
        # opinion inside a byte writer would be a second authority.
        out = tmp_path / "surface.ply"
        result = runner.invoke(
            app, ["mesh", FIXTURE, "--flat", "--crs", "not-a-crs", "--out", str(out)]
        )
        assert result.exit_code == 0, plain(result.output)
        assert "crs not-a-crs" in parse_header(out.read_bytes()).comments


class TestTheSecondFileIsNamedOrNotWritten:
    """Ruling 3's rejected alternative: two files always is a surprise."""

    def test_without_out_edges_only_one_file_appears(self, tmp_path: Path) -> None:
        out = tmp_path / "surface.ply"
        result = runner.invoke(app, ["mesh", FIXTURE, "--flat", "--out", str(out)])
        assert result.exit_code == 0, plain(result.output)
        assert [p.name for p in tmp_path.iterdir()] == ["surface.ply"]


class TestTheAsciiFlag:
    """A small mesh under `--ascii` is inspectable with `head`. Ruling 1."""

    def test_both_files_switch_together(self, tmp_path: Path) -> None:
        surface = tmp_path / "surface.ply"
        edges = tmp_path / "edges.ply"
        result = runner.invoke(
            app,
            ["mesh", FIXTURE, "--flat", "--ascii", "--out", str(surface),
             "--out-edges", str(edges)],
        )
        assert result.exit_code == 0, plain(result.output)
        assert parse_header(surface.read_bytes()).fmt == "ascii"
        assert parse_header(edges.read_bytes()).fmt == "ascii"


class TestRefusals:
    """What the command says no to, and in what words."""

    def test_an_unknown_fixture_names_the_gallery(self, tmp_path: Path) -> None:
        result = runner.invoke(
            app, ["mesh", "no-such-fixture", "--flat", "--out", str(tmp_path / "s.ply")]
        )
        assert result.exit_code != 0
        assert "no-such-fixture" in plain(result.output)
        assert FIXTURE in plain(result.output), "the gallery is listed, as `draw` lists it"

    def test_a_fixture_with_no_mesh_writes_nothing(self, tmp_path: Path) -> None:
        # `degenerate` is refused by the PSLG validator before the noder runs,
        # so there is no mesh to write. `draw` answers that with a picture of
        # the failure; there is no such thing as a picture of a failed file, so
        # this command fails instead, with the engine's own words.
        out = tmp_path / "surface.ply"
        result = runner.invoke(app, ["mesh", "degenerate", "--flat", "--out", str(out)])
        assert result.exit_code != 0
        assert "DegenerateRing" in plain(result.output)
        assert not out.exists()
