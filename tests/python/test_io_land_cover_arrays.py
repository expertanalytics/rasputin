"""The sparse class table in the mesh file: `.vtk` FieldData arrays and `.ply`
face lists (26b, `io.vtk_legacy`, `io.ply`).

`docs/increments/26-land-cover.md`, 1.6 "What the mesh file carries" and
section 10's "26b": the arrays are read back through `vtkread` and `plyread`,
which take every count and type from the file's own declarations and share
nothing with the writers. The pinned spellings are in `fraction_fixtures`'
docstring.
"""

from __future__ import annotations

from itertools import pairwise
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from plyread import parse_header, read_ply
from test_io_vtk_legacy import EDGES, FIELDS, MASKS, TRIANGLES, VERTICES
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.io.models import ClassFractions
from tin_engine.io.ply import write_ply
from tin_engine.io.vtk_legacy import write_vtk
from vtkread import VtkFile, read_vtk

#: Four triangles: two classes, one, two, one. Class 0 (no data) among them,
#: and a fraction that is not exact in binary (0.4) to show the float32 width.
OFFSETS = np.array([0, 2, 3, 5, 6], dtype=np.uint32)
CLASSES = np.array([3, 15, 33, 0, 24, 62], dtype=np.uint8)
FRACTIONS = np.array([0.25, 0.75, 1.0, 0.4, 0.6, 1.0], dtype=np.float32)
#: 1.4 step 7: the largest fraction's class per triangle.
CODES = np.array([15, 33, 24, 62], dtype=np.int32)
NAMES = ("land_cover_offsets", "land_cover_classes", "land_cover_fractions")


def table(**changes: Any) -> Any:
    fields: dict[str, Any] = {"offsets": OFFSETS, "classes": CLASSES, "fractions": FRACTIONS}
    return ClassFractions(**(fields | changes))


def vtk_with(binary: bool, **overrides: Any) -> VtkFile:
    kwargs: dict[str, Any] = {
        "triangles": TRIANGLES,
        "edges": EDGES,
        "edge_masks": MASKS,
        "vocabulary": DEFAULT_VOCABULARY,
        "fields": FIELDS,
        "binary": binary,
        "triangle_codes": CODES,
        "land_cover_codes": "test classes; 0 = no data, and every constraint line",
    }
    return read_vtk(write_vtk(VERTICES, **(kwargs | overrides)))


@pytest.fixture(params=[False, True], ids=["ascii", "binary"])
def binary(request: pytest.FixtureRequest) -> bool:
    return bool(request.param)


class TestVtkFieldData:
    def test_the_three_arrays_are_dataset_field_data_with_their_types(self, binary: bool) -> None:
        vtk = vtk_with(binary, class_fractions=table())
        offsets, classes, fractions = (vtk.field_data[name] for name in NAMES)
        assert offsets.type_name == "unsigned_int"
        assert classes.type_name == "unsigned_char"
        assert fractions.type_name == "float"
        assert (offsets.components, classes.components, fractions.components) == (1, 1, 1)
        assert_array_equal(np.asarray(offsets.values), OFFSETS)
        assert_array_equal(np.asarray(classes.values), CLASSES)
        # float32 on the wire: equal as float32, not merely close.
        assert_array_equal(np.asarray(fractions.values).astype(np.float32), FRACTIONS)

    def test_the_offsets_count_the_triangles_and_index_the_entries(self, binary: bool) -> None:
        vtk = vtk_with(binary, class_fractions=table())
        offsets = np.asarray(vtk.field_data["land_cover_offsets"].values)
        assert len(offsets) == len(vtk.polygons) + 1
        assert offsets[-1] == len(vtk.field_data["land_cover_classes"].values)

    def test_the_code_array_stays_beside_the_table(self, binary: bool) -> None:
        vtk = vtk_with(binary, class_fractions=table())
        codes = vtk.cell_array("land_cover_code")
        assert codes.type_name == "int"
        assert_array_equal(np.asarray(codes.values)[len(vtk.lines) :], CODES)
        assert not np.asarray(codes.values)[: len(vtk.lines)].any()
        (text,) = vtk.field_data["land_cover_codes"].values
        assert "0 =" in str(text)

    def test_without_a_table_the_arrays_are_absent(self, binary: bool) -> None:
        vtk = vtk_with(binary)
        assert not set(NAMES) & set(vtk.field_data)

    def test_the_declared_field_count_covers_the_arrays(self, binary: bool) -> None:
        """The reader takes the FieldData count from the file; a writer that
        appends arrays without raising the count makes the reader refuse."""
        vtk = vtk_with(binary, class_fractions=table())
        assert set(NAMES) <= set(vtk.field_data)
        assert vtk.sections[0] == "FIELD"

    def test_the_names_are_reserved(self) -> None:
        for name in NAMES:
            with pytest.raises(ValueError, match=name):
                vtk_with(False, class_fractions=table(), fields=(*FIELDS, (name, "x")))

    def test_a_table_for_another_triangle_count_is_refused(self) -> None:
        short = table(
            offsets=np.array([0, 2, 3, 5], dtype=np.uint32),
            classes=CLASSES[:5],
            fractions=FRACTIONS[:5],
        )
        with pytest.raises(ValueError):
            vtk_with(False, class_fractions=short)


class TestPlyFaceLists:
    def surface(self, ascii: bool, **overrides: Any) -> bytes:
        kwargs: dict[str, Any] = {
            "faces": TRIANGLES,
            "ascii": ascii,
            "face_codes": CODES,
            "face_fractions": table(),
            "comments": ["land_cover_codes test classes; 0 = no data, and every constraint line"],
        }
        return write_ply(VERTICES, **(kwargs | overrides))

    @pytest.mark.parametrize("ascii", [True, False], ids=["ascii", "binary"])
    def test_the_face_element_declares_the_lists_in_the_designs_order(self, ascii: bool) -> None:
        header = parse_header(self.surface(ascii))
        face = header.element("face")
        assert [p.name for p in face.properties] == [
            "vertex_indices", "land_cover_classes", "land_cover_fractions", "land_cover_code",
        ]  # fmt: skip
        classes, fractions, code = face.properties[1:]
        assert (classes.count_type, classes.type_name) == ("uchar", "uchar")
        assert (fractions.count_type, fractions.type_name) == ("uchar", "float")
        assert (code.count_type, code.type_name) == (None, "int")

    @pytest.mark.parametrize("ascii", [True, False], ids=["ascii", "binary"])
    def test_each_faces_lists_are_its_row_of_the_table(self, ascii: bool) -> None:
        _, data = read_ply(self.surface(ascii))
        face = data["face"]
        assert len(face["land_cover_classes"]) == len(TRIANGLES)
        for t, (a, b) in enumerate(pairwise(OFFSETS)):
            assert_array_equal(np.asarray(face["land_cover_classes"][t]), CLASSES[a:b])
            assert_array_equal(
                np.asarray(face["land_cover_fractions"][t]).astype(np.float32), FRACTIONS[a:b]
            )
        assert_array_equal(np.asarray(face["land_cover_code"]), CODES)
        assert_array_equal(np.stack(list(face["vertex_indices"])), TRIANGLES)

    def test_without_a_table_the_face_has_indices_and_code_only(self) -> None:
        header = parse_header(self.surface(True, face_fractions=None))
        assert [p.name for p in header.element("face").properties] == [
            "vertex_indices", "land_cover_code",
        ]  # fmt: skip

    def test_a_table_beside_edges_is_refused(self) -> None:
        with pytest.raises(ValueError):
            write_ply(VERTICES, edges=EDGES, face_fractions=table())

    def test_a_table_for_another_face_count_is_refused(self) -> None:
        with pytest.raises(ValueError):
            self.surface(True, faces=TRIANGLES[:3], face_codes=CODES[:3])

    def test_a_face_of_256_classes_is_refused_as_the_lists_count_is_a_uchar(self) -> None:
        """`list uchar ...` counts to 255; a row of every uint8 class (the
        only way past it) is refused rather than written short."""
        crowded = table(
            offsets=np.array([0, 256, 257, 258, 259], dtype=np.uint32),
            classes=np.concatenate([np.arange(256), [0, 0, 0]]).astype(np.uint8),
            fractions=np.concatenate([np.full(256, 1 / 256), [1.0, 1.0, 1.0]]).astype(np.float32),
        )
        with pytest.raises(ValueError, match="255"):
            self.surface(True, face_fractions=crowded)
