"""`tin_engine.io.vtk_legacy.write_vtk`: increment 13's byte writer, pinned without VTK.

`13-bundled-mesh.md` names **the cell order (rulings 4 and 5)** as this
increment's invariant-critical suite. The mutant it must kill is a writer that
emits `POLYGONS` before `LINES` with its cell data in that same order: the file
is internally consistent, and VTK -- which always orders lines before polygons
-- gives every value to the wrong cell in silence (finding 1). This suite kills
it through the *file order*; `test_io_vtk_readback.py` kills it again through
the reader's semantics. The second hazard is an array covering the lines only,
which makes `vtkXMLPolyDataReader` return an empty dataset (finding 2); here it
is a count mismatch against `CELL_DATA`.

Everything goes through `vtkread`, a parser written for the suite that takes
every size from a declaration in the blob and refuses whatever it does not
recognise (`PRINCIPLES.md` B1). No `vtk` import appears here: this file runs in
the plain dev environment, and the design's pure half must not depend on the
test-only extra.

It does need the compiled extension, because `tin_engine/__init__.py` imports
`_core`, so no `tin_engine` submodule is importable without it. What it pins
instead is that `vtk_legacy` itself never imports `_core` (ruling 6).

Committed red: the intended failure is `ModuleNotFoundError:
tin_engine.io.vtk_legacy` in every test that writes. The import is inside
`write_vtk` below rather than at module scope, because a collection error
interrupts the whole pytest session and would hide the rest of the suite's
state behind this file's absence.
"""

from __future__ import annotations

import struct

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from importscan import first_party_imports
from tin_engine.features import DEFAULT_VOCABULARY, EdgeProperty, EdgeVocabulary
from vtkread import VtkFile, lines_as_array, polygons_as_array, read_vtk

#: A millimetre at UTM 33N magnitudes, in both axes. float32 cannot hold either:
#: the step at 430 000 is 3.1 cm and at 6 900 000 it is 0.5 m.
EASTING = 430_000.001
NORTHING = 6_900_000.001

#: 17 significant digits in z, which is the caller's and has no lattice behind it.
ELEVATION = 12.345_678_901_234_567

#: Arbitrary doubles at UTM scale, so the round trip is not flattered by values
#: that happen to have short decimal forms.
_RNG = np.random.default_rng(13)
_RANDOM = np.column_stack(
    [
        _RNG.uniform(430_000.0, 440_000.0, 3),
        _RNG.uniform(6_900_000.0, 6_910_000.0, 3),
        _RNG.uniform(-50.0, 2_500.0, 3),
    ]
)

VERTICES = np.vstack(
    [
        [
            [EASTING, NORTHING, ELEVATION],
            [430_010.0, NORTHING, 13.5],
            [430_010.0, 6_900_010.0, 14.0],
        ],
        _RANDOM,
    ]
).astype(np.float64)

TRIANGLES = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 2], [2, 4, 5]], dtype=np.uint32)

EDGES = np.array([[0, 1], [1, 2], [2, 3]], dtype=np.uint32)

#: river; unclassified; road|railway. Distinct and non-zero on the ends, so a
#: shifted or reversed array shows, and 0 in the middle, which is a value the
#: writer must carry rather than skip.
MASKS = np.array([0b001, 0b000, 0b110], dtype=np.uint32)

#: What ruling 5 says `feature_mask` must be: lines first, then 0 per triangle.
EXPECTED_MASK = np.concatenate([MASKS, np.zeros(len(TRIANGLES), dtype=np.uint32)])

CRS = "EPSG:25833"
ELEVATION_TEXT = "none (z=0, --flat)"
FIELDS = (("crs", CRS), ("elevation", ELEVATION_TEXT))

RESERVED = ("feature_bits", "feature_names", "feature_vocabulary")


def write_vtk(vertices: object, **kwargs: object) -> bytes:
    """The writer under test, imported per call so its absence is a per-test red."""
    from tin_engine.io.vtk_legacy import write_vtk as writer

    return writer(vertices, **kwargs)  # type: ignore[arg-type]


def write(**overrides: object) -> bytes:
    """`write_vtk` over the module's mesh, with any argument replaced."""
    kwargs: dict[str, object] = {
        "triangles": TRIANGLES,
        "edges": EDGES,
        "edge_masks": MASKS,
        "vocabulary": DEFAULT_VOCABULARY,
        "fields": FIELDS,
    }
    kwargs.update(overrides)
    vertices = kwargs.pop("vertices", VERTICES)
    return write_vtk(vertices, **kwargs)  # type: ignore[arg-type]


@pytest.fixture(params=[False, True], ids=["ascii", "binary"])
def binary(request: pytest.FixtureRequest) -> bool:
    return bool(request.param)


@pytest.fixture
def parsed(binary: bool) -> VtkFile:
    return read_vtk(write(binary=binary))


class TestHeader:
    """The five fixed lines. Parsed as text, not compared to a golden blob."""

    def test_version_is_4_2(self, parsed: VtkFile) -> None:
        # 4.2 and not 5.1: the older cell layout every VTK reader accepts.
        assert parsed.version == "4.2"

    def test_the_title_is_fixed_text(self, parsed: VtkFile) -> None:
        # Ruling 7: a title line has no escaping, so it never holds user input.
        assert parsed.title == "rasputin mesh"

    def test_ascii_is_the_default(self) -> None:
        assert read_vtk(write()).encoding == "ASCII"

    def test_binary_is_behind_the_flag(self) -> None:
        assert read_vtk(write(binary=True)).encoding == "BINARY"

    def test_the_dataset_is_polydata(self, parsed: VtkFile) -> None:
        assert parsed.dataset == "POLYDATA"

    def test_points_are_declared_double(self, parsed: VtkFile) -> None:
        assert parsed.points_type == "double"

    def test_the_first_bytes_are_the_version_line(self) -> None:
        assert write().startswith(b"# vtk DataFile Version 4.2\nrasputin mesh\nASCII\n")

    def test_it_is_ascii_through_and_through(self) -> None:
        write().decode("ascii")


class TestPoints:
    """Ruling 2: double end to end, bit for bit."""

    def test_the_probe_can_fail(self) -> None:
        # A3: if float32 held these, the round trip below would not refute it.
        assert float(np.float32(EASTING)) != EASTING
        assert float(np.float32(NORTHING)) != NORTHING

    def test_points_survive_exactly(self, parsed: VtkFile) -> None:
        assert_array_equal(parsed.points, VERTICES)

    def test_the_dtype_read_back_is_float64(self, parsed: VtkFile) -> None:
        assert parsed.points.dtype.kind == "f" and parsed.points.dtype.itemsize == 8

    def test_ascii_uses_the_shortest_round_trip_text(self) -> None:
        # `repr`, per ruling 2: `%.6f` or `%g` loses the 17th digit of z.
        assert repr(ELEVATION).encode("ascii") in write()


class TestCellOrder:
    """The invariant-critical suite: rulings 4 and 5."""

    def test_lines_come_before_polygons_in_the_file(self, parsed: VtkFile) -> None:
        assert parsed.sections.index("LINES") < parsed.sections.index("POLYGONS")

    def test_no_other_cell_kind_is_written(self, parsed: VtkFile) -> None:
        assert not {"VERTICES", "TRIANGLE_STRIPS"} & set(parsed.sections)

    def test_the_lines_are_the_edges_in_order(self, parsed: VtkFile) -> None:
        assert_array_equal(lines_as_array(parsed), EDGES)

    def test_the_polygons_are_the_triangles_with_their_winding(self, parsed: VtkFile) -> None:
        assert_array_equal(polygons_as_array(parsed), TRIANGLES)

    def test_the_cell_list_sizes_are_the_legacy_ones(self, parsed: VtkFile) -> None:
        # `LINES n size` counts one length word per cell plus its point ids.
        assert parsed.lines_header == (len(EDGES), 3 * len(EDGES))
        assert parsed.polygons_header == (len(TRIANGLES), 4 * len(TRIANGLES))

    def test_cell_data_counts_every_cell(self, parsed: VtkFile) -> None:
        assert parsed.cell_count == len(EDGES) + len(TRIANGLES)

    def test_cell_data_follows_the_cells(self, parsed: VtkFile) -> None:
        assert parsed.sections[-1] == "CELL_DATA"
        assert parsed.sections.index("POLYGONS") < parsed.sections.index("CELL_DATA")

    def test_feature_mask_is_the_masks_then_zero_per_triangle(self, parsed: VtkFile) -> None:
        assert_array_equal(parsed.cell_array("feature_mask").values, EXPECTED_MASK)

    def test_feature_mask_is_the_active_scalar(self, parsed: VtkFile) -> None:
        # Ruling 5: the SCALARS block, so ParaView offers it as the cell scalar.
        # `vtkPolyDataReader` returns only the first SCALARS block by default
        # (finding 4), so it must be the only one.
        assert list(parsed.scalars) == ["feature_mask"]
        assert parsed.scalars["feature_mask"].type_name == "unsigned_int"
        assert parsed.scalars["feature_mask"].components == 1
        assert parsed.cell_sections[0] == "SCALARS"

    def test_every_cell_array_covers_every_cell(self, parsed: VtkFile) -> None:
        # Finding 2: a cell array short of the cell count empties the dataset.
        # `vtkread` refuses a short FIELD array on its own; this pins it for the
        # scalar too, and states the invariant where a reader looks for it.
        cells = len(EDGES) + len(TRIANGLES)
        assert len(parsed.scalars["feature_mask"].values) == cells
        for arrays in parsed.cell_fields.values():
            for array in arrays.values():
                assert len(array.values) == cells, array.name


class TestPerFeatureArrays:
    """Ruling 6: one 0/1 array per property set on at least one edge."""

    def test_they_sit_in_a_field_block_named_features(self, parsed: VtkFile) -> None:
        # FIELD rather than SCALARS: all three readers return FIELD arrays,
        # while `vtkPolyDataReader` drops every SCALARS block after the first.
        assert list(parsed.cell_fields) == ["features"]

    def test_exactly_the_properties_that_occur(self, parsed: VtkFile) -> None:
        # MASKS sets river (bit 0), road (bit 1) and railway (bit 2). The four
        # other DEFAULT_VOCABULARY properties are set nowhere, so they are absent.
        assert set(parsed.cell_fields["features"]) == {"river", "road", "railway"}

    def test_each_value_is_the_bits_state_and_triangles_carry_0(self, parsed: VtkFile) -> None:
        zeros = [0] * len(TRIANGLES)
        features = parsed.cell_fields["features"]
        assert features["river"].values.tolist() == [1, 0, 0, *zeros]
        assert features["road"].values.tolist() == [0, 0, 1, *zeros]
        assert features["railway"].values.tolist() == [0, 0, 1, *zeros]

    def test_they_are_unsigned_char_scalars(self, parsed: VtkFile) -> None:
        for array in parsed.cell_fields["features"].values():
            assert (array.type_name, array.components) == ("unsigned_char", 1), array.name

    def test_no_bits_set_writes_no_features_block(self, binary: bool) -> None:
        parsed = read_vtk(write(edge_masks=np.zeros(len(EDGES), dtype=np.uint32), binary=binary))
        assert parsed.cell_fields == {}
        assert_array_equal(parsed.cell_array("feature_mask").values, np.zeros(parsed.cells))

    def test_they_are_named_by_this_vocabulary(self, binary: bool) -> None:
        # Names come from the vocabulary handed in, not from DEFAULT_VOCABULARY.
        vocabulary = EdgeVocabulary(properties=(EdgeProperty(name="dyke", bit=2),))
        masks = np.array([0, 4, 0], dtype=np.uint32)
        parsed = read_vtk(write(edge_masks=masks, vocabulary=vocabulary, binary=binary))
        assert set(parsed.cell_fields["features"]) == {"dyke"}
        assert parsed.cell_fields["features"]["dyke"].values.tolist()[:3] == [0, 1, 0]


class TestDatasetFieldData:
    """Ruling 6: the vocabulary travels in the file, as names and a fingerprint."""

    def test_it_is_a_dataset_field_before_the_points(self, parsed: VtkFile) -> None:
        assert parsed.sections[0] == "FIELD"
        assert parsed.field_data_name == "FieldData"

    def test_the_table_is_the_whole_vocabulary_sorted_by_bit(self, parsed: VtkFile) -> None:
        pairs = sorted((p.bit, p.name) for p in DEFAULT_VOCABULARY.properties)
        assert parsed.field_data["feature_bits"].values.tolist() == [b for b, _ in pairs]
        assert parsed.field_data["feature_names"].values == tuple(n for _, n in pairs)

    def test_the_table_types(self, parsed: VtkFile) -> None:
        assert parsed.field_data["feature_bits"].type_name == "unsigned_int"
        assert parsed.field_data["feature_names"].type_name == "string"

    def test_a_vocabulary_declared_out_of_order_is_written_sorted(self, binary: bool) -> None:
        vocabulary = EdgeVocabulary(
            properties=(
                EdgeProperty(name="wall", bit=5),
                EdgeProperty(name="river", bit=0),
                EdgeProperty(name="road", bit=1),
            )
        )
        parsed = read_vtk(write(vocabulary=vocabulary, edge_masks=MASKS & 0b11, binary=binary))
        assert parsed.field_data["feature_bits"].values.tolist() == [0, 1, 5]
        assert parsed.field_data["feature_names"].values == ("river", "road", "wall")

    def test_the_fingerprint_is_the_vocabularys(self, parsed: VtkFile) -> None:
        array = parsed.field_data["feature_vocabulary"]
        assert array.type_name == "string"
        assert array.values == (DEFAULT_VOCABULARY.fingerprint(),)

    def test_crs_and_elevation_appear_when_given(self, parsed: VtkFile) -> None:
        assert parsed.field_data["crs"].values == (CRS,)
        assert parsed.field_data["elevation"].values == (ELEVATION_TEXT,)
        assert parsed.field_data["crs"].type_name == "string"

    def test_no_crs_field_when_none_is_given(self) -> None:
        parsed = read_vtk(write(fields=(("elevation", ELEVATION_TEXT),)))
        assert "crs" not in parsed.field_data


class TestStrings:
    """Ruling 7: one encoder, VTK's own `%XX` convention, and no forged line."""

    AWKWARD = "+proj=utm +zone=33 +units=m 50%"

    def test_space_and_percent_are_encoded_in_ascii(self) -> None:
        parsed = read_vtk(write(fields=(("crs", self.AWKWARD),)))
        assert parsed.field_data["crs"].raw == (b"+proj=utm%20+zone=33%20+units=m%2050%25",)

    def test_they_decode_back_to_the_original(self, binary: bool) -> None:
        parsed = read_vtk(write(fields=(("crs", self.AWKWARD),), binary=binary))
        assert parsed.field_data["crs"].values == (self.AWKWARD,)

    def test_the_elevation_text_is_encoded_too(self) -> None:
        parsed = read_vtk(write())
        assert parsed.field_data["elevation"].raw == (b"none%20(z=0,%20--flat)",)

    def test_a_literal_escape_is_not_decoded_twice(self, binary: bool) -> None:
        # `%` must be escaped first, or `%20` typed by a user comes back a space.
        parsed = read_vtk(write(fields=(("crs", "a%20b"),), binary=binary))
        assert parsed.field_data["crs"].values == ("a%20b",)

    @pytest.mark.parametrize("bad", ["\n", "\r", "\t", "\x00", "\x1b", "\x7f"])
    def test_a_control_character_is_refused_by_name(self, bad: str) -> None:
        with pytest.raises(ValueError, match="control character") as raised:
            write(fields=(("crs", f"EPSG:25833{bad}POINTS 1 double"),))
        assert repr(bad) in str(raised.value)

    @pytest.mark.parametrize("bad", ["\N{DEGREE SIGN}", "\N{LATIN SMALL LETTER O WITH STROKE}"])
    def test_non_ascii_is_refused_by_name(self, bad: str) -> None:
        with pytest.raises(ValueError, match="ASCII") as raised:
            write(fields=(("crs", f"ETRS89 60{bad}N"),))
        assert repr(bad) in str(raised.value)

    def test_the_elevation_value_goes_through_the_same_encoder(self) -> None:
        with pytest.raises(ValueError, match="control character"):
            write(fields=(("elevation", "none\nPOINTS 1 double"),))


class TestFieldNames:
    """A field name is a token on a header line, so it is a pattern, not free text."""

    @pytest.mark.parametrize("reserved", RESERVED)
    def test_a_reserved_name_is_refused(self, reserved: str) -> None:
        with pytest.raises(ValueError, match=reserved):
            write(fields=((reserved, "x"),))

    @pytest.mark.parametrize("name", ["", "CRS", "1crs", "c rs", "crs\n", "cr-s", "_crs"])
    def test_a_name_outside_the_pattern_is_refused(self, name: str) -> None:
        with pytest.raises(ValueError):
            write(fields=((name, "x"),))


class TestRefusals:
    """The `ValueError`s the design's signature lists."""

    def test_an_unnamed_bit_is_refused_at_write_time(self) -> None:
        # Increment 7's mechanism 3, through `vocabulary.names()`: a mask from a
        # different vocabulary is refused, not written as an unexplained number.
        with pytest.raises(ValueError, match="bit 9"):
            write(edge_masks=np.array([1, 1 << 9, 0], dtype=np.uint32))

    def test_a_bit_named_by_another_vocabulary_is_refused(self) -> None:
        only_river = EdgeVocabulary(properties=(EdgeProperty(name="river", bit=0),))
        with pytest.raises(ValueError, match="bit 1"):
            write(vocabulary=only_river)

    def test_two_column_vertices_are_refused(self) -> None:
        with pytest.raises(ValueError, match=r"\(N, 3\)|shape"):
            write(vertices=VERTICES[:, :2])

    @pytest.mark.parametrize("count", [len(EDGES) - 1, len(EDGES) + 1])
    def test_masks_must_have_one_entry_per_edge(self, count: int) -> None:
        with pytest.raises(ValueError, match="edge"):
            write(edge_masks=np.zeros(count, dtype=np.uint32))


class TestZeroEdges:
    """A mesh with no constraint edges is still a well-formed file (ruling 4)."""

    @pytest.fixture
    def empty(self, binary: bool) -> VtkFile:
        return read_vtk(
            write(
                edges=np.zeros((0, 2), dtype=np.uint32),
                edge_masks=np.zeros(0, dtype=np.uint32),
                binary=binary,
            )
        )

    def test_lines_is_written_as_0_0(self, empty: VtkFile) -> None:
        assert "LINES" in empty.sections
        assert empty.lines_header == (0, 0)

    def test_the_ascii_line_is_literally_lines_0_0(self) -> None:
        blob = write(edges=np.zeros((0, 2), dtype=np.uint32), edge_masks=np.zeros(0, np.uint32))
        assert b"\nLINES 0 0\n" in blob

    def test_the_rest_of_the_file_is_intact(self, empty: VtkFile) -> None:
        assert_array_equal(empty.points, VERTICES)
        assert_array_equal(polygons_as_array(empty), TRIANGLES)
        assert empty.cell_count == len(TRIANGLES)
        assert_array_equal(empty.cell_array("feature_mask").values, np.zeros(len(TRIANGLES)))
        assert empty.cell_fields == {}


class TestAsciiMatchesBinary:
    """One mesh, two encodings, the same arrays. Catches one mode drifting."""

    def test_everything_decodes_equal(self) -> None:
        text, packed = read_vtk(write()), read_vtk(write(binary=True))
        assert_array_equal(text.points, packed.points)
        assert_array_equal(lines_as_array(text), lines_as_array(packed))
        assert_array_equal(polygons_as_array(text), polygons_as_array(packed))
        assert text.sections == packed.sections
        assert text.cell_count == packed.cell_count
        assert_array_equal(
            text.cell_array("feature_mask").values, packed.cell_array("feature_mask").values
        )
        assert set(text.cell_fields["features"]) == set(packed.cell_fields["features"])
        for name, array in text.cell_fields["features"].items():
            assert_array_equal(array.values, packed.cell_fields["features"][name].values)
        assert set(text.field_data) == set(packed.field_data)
        for name, array in text.field_data.items():
            other = packed.field_data[name].values
            if isinstance(array.values, tuple):
                assert array.values == other, name
            else:
                assert_array_equal(array.values, other)


class TestBinaryIsBigEndian:
    """Ruling 3: big-endian by the format's definition -- the PLY writer's opposite.

    Asserted on the raw bytes, with values whose little-endian reading differs,
    because `vtkread` decoding big-endian would agree with a writer that was
    wrong in the same direction only if both were; the raw bytes are neither.
    """

    @pytest.fixture
    def blob(self) -> bytes:
        return write(binary=True)

    def test_the_probe_can_fail(self) -> None:
        assert struct.pack(">d", EASTING) != struct.pack("<d", EASTING)
        assert struct.pack(">I", 6) != struct.pack("<I", 6)

    def test_points_are_big_endian_doubles(self, blob: bytes) -> None:
        declared = f"POINTS {len(VERTICES)} double\n".encode("ascii")
        assert declared + VERTICES.astype(">f8").tobytes() in blob
        assert VERTICES.astype("<f8").tobytes() not in blob

    def test_cells_are_big_endian_ints(self, blob: bytes) -> None:
        flat = np.column_stack([np.full(len(EDGES), 2), EDGES]).ravel()
        declared = f"LINES {len(EDGES)} {flat.size}\n".encode("ascii")
        assert declared + flat.astype(">i4").tobytes() in blob

    def test_feature_mask_is_big_endian(self, blob: bytes) -> None:
        assert b"LOOKUP_TABLE default\n" + EXPECTED_MASK.astype(">u4").tobytes() in blob
        assert EXPECTED_MASK.astype("<u4").tobytes() not in blob


class TestPurity:
    """Ruling 6: bytes out, no path in, and never `_core`."""

    def test_the_return_value_is_bytes(self) -> None:
        assert isinstance(write(), bytes)

    def test_it_is_deterministic(self) -> None:
        assert write() == write()
        assert write(binary=True) == write(binary=True)

    def test_the_only_first_party_import_is_the_vocabulary(self) -> None:
        # `features.py` imports only hashlib and pydantic, so depending on it
        # adds nothing to what `io/` can reach. `_core` and `viz` stay out.
        import tin_engine.io.vtk_legacy as module

        assert first_party_imports(module) <= {"tin_engine.features"}

    def test_it_is_re_exported_from_io(self) -> None:
        from tin_engine.io.vtk_legacy import write_vtk as writer

        import tin_engine.io as io

        assert io.write_vtk is writer
        assert "write_vtk" in io.__all__

