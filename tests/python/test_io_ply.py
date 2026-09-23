"""`tin_engine.io.ply.write_ply`: increment 10's byte writer, committed RED.

`10-mesh-output.md` names the **binary body** as this increment's
invariant-critical suite, so this file carries the mutation round. The defect
class it exists to catch is the one a count assertion cannot see: a wrong
offset, a wrong dtype or a swapped index produces a file that is structurally
valid, opens without complaint, and is geometrically wrong.

Everything here goes through `plyread`, a parser written for the suite that
reads the header it is handed and takes every size from a declaration in the
blob. It never asks the writer where anything is. That is deliberate and it is
the difference between a test and a restatement of the implementation
(`PRINCIPLES.md` B1).

The intended failure at this commit is `ModuleNotFoundError: tin_engine.io`.
The package does not exist; `write_ply` is imported at module scope so that the
absence is one collection error rather than eleven confusing ones.

TWO THINGS THIS SUITE DELIBERATELY DOES NOT ASSERT.

1. That QGIS or ParaView opens the result. Nobody here can run either, and
   ruling 2 forbids the claim until a person has opened a file and said so in
   the record (`PRINCIPLES.md` A1, A4). What the suite asserts is what the
   bytes are.
2. `cli._destination`'s refusals. They are increment 6b-ii's and are covered in
   `test_cli_draw.py`; re-testing them here would pin one boundary in two
   places.

THE PROPERTY NAMES ARE THIS SUITE'S CHOICE, NOT THE DESIGN'S. The design fixes
the element types and the coordinate width and says nothing about spelling. The
names pinned below -- `x`/`y`/`z`, `vertex_indices`, `vertex1`/`vertex2` -- are
PLY's conventional ones and the ones MDAL's driver documentation uses, so they
are the names a reader outside this repository will look for. The per-edge
feature scalar has no convention, so it is asserted *positionally*: whatever
the edge element's third property is called, its values must be the masks that
were handed in.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from tin_engine.io.ply import write_ply

from plyread import element_bytes, parse_header, read_ply, vertex_array

#: A millimetre-resolved easting at UTM 33N magnitudes -- `viz.fixtures.ORIGIN`
#: is 430_000 -- carried by ruling 2 as the probe that float32 cannot survive.
#: `TestDoublePrecision` proves the probe can fail before trusting it (A3).
EASTING = 430_000.001

#: A z the caller supplied and nothing snapped. The grid rounds x and y; ruling
#: 4 says z is the caller's, so it is the coordinate with no lattice behind it
#: and the one a fixed-decimal ASCII format loses first.
ELEVATION = 12.345_678_901_2

VERTICES = np.array(
    [
        [EASTING, 6_900_000.5, ELEVATION],
        [430_010.0, 6_900_000.5, 13.5],
        [430_010.0, 6_900_010.0, 14.0],
        [EASTING, 6_900_010.0, 11.0],
    ],
    dtype=np.float64,
)

FACES = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)

EDGES = np.array([[0, 1], [1, 2], [2, 3], [3, 0]], dtype=np.uint32)

#: Feature masks, including a 0. `_core.pyi` defines 0 as *unclassified*
#: rather than *wrong*, so it is a value the writer must carry, not skip.
EDGE_PROPERTIES = np.array([1, 0, 4, 6], dtype=np.uint32)


@pytest.fixture
def surface() -> bytes:
    """The 2D mesh file: vertices and faces, binary, no edges. Ruling 3."""
    return write_ply(VERTICES, faces=FACES)


@pytest.fixture
def constraints() -> bytes:
    """The 1D mesh file: the same vertices, and edges instead of faces."""
    return write_ply(VERTICES, edges=EDGES, edge_properties=EDGE_PROPERTIES)


class TestHeader:
    """What the header declares. Parsed as text, not compared to a golden blob."""

    def test_binary_little_endian_is_the_default(self, surface: bytes) -> None:
        header = parse_header(surface)
        assert surface.startswith(b"ply\n")
        assert header.fmt == "binary_little_endian"
        assert header.version == "1.0"

    def test_ascii_is_available_behind_the_flag(self) -> None:
        header = parse_header(write_ply(VERTICES, faces=FACES, ascii=True))
        assert header.fmt == "ascii"

    def test_coordinates_are_declared_double(self, surface: bytes) -> None:
        # Ruling 2, at the declaration. The wire width is checked separately,
        # because a header saying `double` over a float32 body is precisely the
        # mutant that a declaration test alone would miss.
        vertex = parse_header(surface).element("vertex")
        assert [(p.name, p.type_name) for p in vertex.properties] == [
            ("x", "double"),
            ("y", "double"),
            ("z", "double"),
        ]

    def test_counts_are_the_array_lengths(self, surface: bytes, constraints: bytes) -> None:
        assert parse_header(surface).element("vertex").count == len(VERTICES)
        assert parse_header(surface).element("face").count == len(FACES)
        assert parse_header(constraints).element("edge").count == len(EDGES)

    def test_faces_are_a_vertex_index_list(self, surface: bytes) -> None:
        (prop,) = parse_header(surface).element("face").properties
        assert prop.name == "vertex_indices"
        assert prop.is_list

    def test_edges_name_their_two_endpoints_and_carry_one_scalar(
        self, constraints: bytes
    ) -> None:
        properties = parse_header(constraints).element("edge").properties
        assert [p.name for p in properties[:2]] == ["vertex1", "vertex2"]
        assert len(properties) == 3, "one scalar beyond the endpoints: the feature mask"
        assert not properties[2].is_list

    def test_edges_without_properties_declare_only_the_endpoints(self) -> None:
        properties = parse_header(write_ply(VERTICES, edges=EDGES)).element("edge").properties
        assert [p.name for p in properties] == ["vertex1", "vertex2"]

    def test_comments_are_carried_in_order(self) -> None:
        given = ("crs EPSG:25833", "elevation none (z=0, --flat)")
        assert parse_header(write_ply(VERTICES, faces=FACES, comments=given)).comments == given

    def test_a_file_holds_one_element_type_beside_the_vertices(
        self, surface: bytes, constraints: bytes
    ) -> None:
        # Ruling 3 as the reader sees it: MDAL's caveat is that a host expects
        # either a 1D mesh or a 2D one, so neither file may declare both.
        assert parse_header(surface).names == ("vertex", "face")
        assert parse_header(constraints).names == ("vertex", "edge")


class TestBinaryRoundTrip:
    """The invariant-critical suite. Write arrays, parse the bytes, compare."""

    def test_vertices_survive_exactly(self, surface: bytes) -> None:
        _, data = read_ply(surface)
        assert_array_equal(vertex_array(data), VERTICES)

    def test_faces_survive_with_their_winding(self, surface: bytes) -> None:
        # assert_array_equal and not a set comparison: a writer that reversed
        # or rotated each triangle would round trip through a set unnoticed,
        # and a mesh whose faces wind inconsistently shades inside out.
        _, data = read_ply(surface)
        assert_array_equal(data["face"]["vertex_indices"], FACES)

    def test_edges_and_their_masks_survive(self, constraints: bytes) -> None:
        header, data = read_ply(constraints)
        mask_name = header.element("edge").properties[2].name
        assert_array_equal(data["edge"]["vertex1"], EDGES[:, 0])
        assert_array_equal(data["edge"]["vertex2"], EDGES[:, 1])
        assert_array_equal(data["edge"][mask_name], EDGE_PROPERTIES)

    def test_the_body_ends_where_the_declarations_say(self, surface: bytes) -> None:
        # `read_ply` raises on a body that is longer or shorter than the
        # declared elements require. Padding, a missing record and an off-by-one
        # header terminator all land here, and none of them is visible to an
        # assertion about the arrays alone.
        read_ply(surface)
        read_ply(write_ply(VERTICES, edges=EDGES, edge_properties=EDGE_PROPERTIES))

    def test_a_single_triangle_is_a_whole_file(self) -> None:
        # The smallest mesh there is. A writer that assumed more than one
        # record per block, or that wrote a separator between them, fails here
        # and nowhere else.
        _, data = read_ply(write_ply(VERTICES[:3], faces=FACES[:1]))
        assert_array_equal(vertex_array(data), VERTICES[:3])
        assert_array_equal(data["face"]["vertex_indices"], FACES[:1])


class TestDoublePrecision:
    """Ruling 2: float32 is a corruption, and this is the probe that says so."""

    def test_the_probe_can_fail(self) -> None:
        # A3. Written before the round trip is trusted: if EASTING survived a
        # float32 trip, the test below would pass against the dtype it exists
        # to refuse. The float32 step at this magnitude is 2**19 * 2**-24.
        assert float(np.float32(EASTING)) != EASTING
        assert np.spacing(np.float32(EASTING)) > 0.03

    def test_a_millimetre_survives_the_binary_body(self, surface: bytes) -> None:
        _, data = read_ply(surface)
        assert data["vertex"]["x"][0] == EASTING

    def test_a_millimetre_survives_the_ascii_body(self) -> None:
        # The ASCII path fails this with `%.6f`-style formatting or with a
        # float32 cast, which are two different bugs with one symptom: the
        # noder's millimetre lattice quietly gone.
        _, data = read_ply(write_ply(VERTICES, faces=FACES, ascii=True))
        assert data["vertex"]["x"][0] == EASTING

    def test_the_vertex_block_is_eight_bytes_per_coordinate(self, surface: bytes) -> None:
        block = element_bytes(surface, "vertex")
        assert len(block) == 3 * 8 * len(VERTICES)


class TestAsciiMatchesBinary:
    """One geometry, two encodings. This is what catches one mode drifting."""

    def test_vertices_agree(self) -> None:
        _, binary = read_ply(write_ply(VERTICES, faces=FACES))
        _, text = read_ply(write_ply(VERTICES, faces=FACES, ascii=True))
        assert_array_equal(vertex_array(text), vertex_array(binary))

    def test_faces_agree(self) -> None:
        _, binary = read_ply(write_ply(VERTICES, faces=FACES))
        _, text = read_ply(write_ply(VERTICES, faces=FACES, ascii=True))
        assert_array_equal(text["face"]["vertex_indices"], binary["face"]["vertex_indices"])

    def test_edges_and_masks_agree(self) -> None:
        kwargs = {"edges": EDGES, "edge_properties": EDGE_PROPERTIES}
        binary_header, binary = read_ply(write_ply(VERTICES, **kwargs))
        _, text = read_ply(write_ply(VERTICES, ascii=True, **kwargs))
        mask = binary_header.element("edge").properties[2].name
        assert_array_equal(text["edge"]["vertex1"], binary["edge"]["vertex1"])
        assert_array_equal(text["edge"]["vertex2"], binary["edge"]["vertex2"])
        assert_array_equal(text["edge"][mask], binary["edge"][mask])

    def test_the_headers_differ_only_in_the_format_line(self) -> None:
        binary = parse_header(write_ply(VERTICES, faces=FACES))
        text = parse_header(write_ply(VERTICES, faces=FACES, ascii=True))
        assert binary.elements == text.elements
        assert binary.fmt != text.fmt


class TestExactlyOneElementType:
    """Ruling 3, enforced at the one place it can be: the writer's signature."""

    def test_faces_and_edges_together_are_refused(self) -> None:
        with pytest.raises(ValueError, match="faces") as raised:
            write_ply(VERTICES, faces=FACES, edges=EDGES)
        assert "edges" in str(raised.value), "the message must name both, not just the first"

    def test_neither_is_refused(self) -> None:
        with pytest.raises(ValueError, match="faces") as raised:
            write_ply(VERTICES)
        assert "edges" in str(raised.value)

    def test_edge_properties_without_edges_is_refused(self) -> None:
        with pytest.raises(ValueError, match="edge_properties") as raised:
            write_ply(VERTICES, faces=FACES, edge_properties=EDGE_PROPERTIES)
        assert "edges" in str(raised.value)


class TestTheWriterNeverInventsAZ:
    """Ruling 4: the caller supplies z, and a 2D array is a caller's mistake."""

    def test_two_column_vertices_are_refused_rather_than_filled(self) -> None:
        with pytest.raises(ValueError, match=r"\(N, 3\)|shape"):
            write_ply(VERTICES[:, :2], faces=FACES)

    def test_a_supplied_zero_z_is_written_as_zero(self) -> None:
        flat = np.column_stack([VERTICES[:, :2], np.zeros(len(VERTICES))])
        _, data = read_ply(write_ply(flat, faces=FACES))
        assert_array_equal(data["vertex"]["z"], np.zeros(len(VERTICES)))


class TestTheTwoFilesRegister:
    """Ruling 3's second half, asserted rather than assumed.

    The surface file and the constraint file repeat the same vertex block
    because PLY has no cross-file reference. Vertex *i* means the same point in
    both only if those bytes are identical, and that is what makes the two
    layers line up when a person loads them side by side.
    """

    def test_the_vertex_blocks_are_byte_identical(
        self, surface: bytes, constraints: bytes
    ) -> None:
        assert element_bytes(surface, "vertex") == element_bytes(constraints, "vertex")

    def test_the_vertex_blocks_are_not_empty(self, surface: bytes) -> None:
        # Guards the assertion above from passing on two empty slices, which is
        # what it would do if `element_bytes` ever stopped finding the block.
        assert element_bytes(surface, "vertex")


class TestPurity:
    """Ruling 6: bytes out, no path in, no first-party import."""

    def test_the_return_value_is_bytes(self, surface: bytes) -> None:
        assert isinstance(surface, bytes)

    def test_the_module_imports_nothing_first_party(self) -> None:
        # A writer that reached for `_core` would need the extension built to
        # test, and a writer that reached for `viz` would couple output to the
        # renderer -- the design refuses both by name. Read from the module's
        # own import statements, so a mention in prose is not a finding.
        import ast
        import inspect

        import tin_engine.io.ply as module

        tree = ast.parse(inspect.getsource(module))
        roots = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                roots.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                roots.add(node.module.split(".")[0])
            elif isinstance(node, ast.ImportFrom) and node.level:
                roots.add("tin_engine")  # a relative import is a first-party one
        assert "tin_engine" not in roots

    def test_it_is_deterministic(self) -> None:
        assert write_ply(VERTICES, faces=FACES) == write_ply(VERTICES, faces=FACES)
