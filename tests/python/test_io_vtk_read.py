"""Reading rasputin's own mesh back: `tin_engine.io.vtk_read.read_mesh`.

`docs/increments/51-open-core.md`, *3.1 What moves to the private package*:
the reader stays in the core, offered to plugins. `read_mesh(data:
bytes) -> MeshArrays` reads the binary legacy VTK that `rasputin mesh` writes
(header line `rasputin mesh`): the `crs` string from the dataset's field data,
POINTS and POLYGONS, skipping LINES, and the triangles' `land_cover_code`
cell array with the dataset string `land_cover_codes` beside it
(`docs/increments/81-read-land-cover.md`). It refuses ASCII ("write the mesh without
--ascii"), another program's file and a file with no `crs`. The round trip
`read_mesh(write_vtk(...))` gives the same vertices, triangles and crs bit for
bit; the three refusals are in plain words.

Interface:
`MeshArrays` has `vertices` (float64, `(n, 3)`), `triangles` (int64,
`(T, 3)`) and `crs` (str), then `land_cover_code` (int32 `(T,)`, the line
cells' entries dropped, or `None` when the file has no such array) and
`land_cover_codes` (str, `""` when absent); a refusal is a `ValueError`.
"""

from __future__ import annotations

import importlib
from types import ModuleType

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from tin_engine.features import DEFAULT_VOCABULARY, EdgeVocabulary
from tin_engine.io.models import ClassFractions
from tin_engine.io.vtk_legacy import write_vtk

RNG = np.random.default_rng(4604)
VERTICES = np.column_stack(
    [
        494000.0 + RNG.random(9) * 12000.0,
        6794000.0 + RNG.random(9) * 12000.0,
        234.0 + RNG.random(9) * 1690.0,
    ]
)
TRIANGLES = np.array([[0, 1, 2], [2, 1, 3], [3, 4, 5], [5, 6, 7], [7, 8, 0]])
CRS = "EPSG:25833"
#: A negative code and the int32 maximum catch an unsigned or wrong-width read.
CODES = np.array([142, 322, 332, -7, 2**31 - 1], dtype=np.int32)
CODE_SYSTEM = "CORINE Land Cover level-3 code"
#: Two constraint lines, so the file's cell arrays have 7 entries, lines first.
LINES = np.array([[0, 1], [1, 2]])
#: The array's header as the writer emits it for LINES and the 5 triangles; it
#: is the file's last array: header (24 bytes), 7 big-endian int32, a newline.
CODE_HEADER = b"land_cover_code 1 7 int\n"


@pytest.fixture(scope="module")
def vtk_read() -> ModuleType:
    return importlib.import_module("tin_engine.io.vtk_read")


def mesh_bytes(
    *,
    binary: bool = True,
    fields: tuple[tuple[str, str], ...] = (("crs", CRS),),
    edges: np.ndarray | None = None,
    vocabulary: EdgeVocabulary = DEFAULT_VOCABULARY,
    edge_mask: int = 1,
    triangle_codes: np.ndarray | None = None,
    land_cover_codes: str = "",
    class_fractions: ClassFractions | None = None,
) -> bytes:
    lines = np.empty((0, 2), dtype=np.int64) if edges is None else edges
    masks = np.full(len(lines), edge_mask, dtype=np.uint32)
    return write_vtk(
        VERTICES,
        triangles=TRIANGLES,
        edges=lines,
        edge_masks=masks if len(lines) else np.empty(0),
        vocabulary=vocabulary,
        fields=fields,
        binary=binary,
        triangle_codes=triangle_codes,
        land_cover_codes=land_cover_codes,
        class_fractions=class_fractions,
    )


@pytest.fixture(scope="module")
def coded() -> bytes:
    """The increment's fixture: lines, codes and the code system's string."""
    data = mesh_bytes(edges=LINES, triangle_codes=CODES, land_cover_codes=CODE_SYSTEM)
    body = np.concatenate([np.zeros(len(LINES), np.int32), CODES]).astype(">i4").tobytes()
    assert data.endswith(CODE_HEADER + body + b"\n")
    return data


class TestRoundTrip:
    def test_vertices_triangles_and_crs_come_back_bit_for_bit(self, vtk_read: ModuleType) -> None:
        mesh = vtk_read.read_mesh(mesh_bytes())
        assert mesh.vertices.dtype == np.float64
        assert_array_equal(mesh.vertices, VERTICES)
        assert mesh.triangles.dtype == np.int64
        assert_array_equal(mesh.triangles, TRIANGLES)
        assert mesh.crs == CRS

    def test_constraint_lines_are_skipped(self, vtk_read: ModuleType) -> None:
        mesh = vtk_read.read_mesh(mesh_bytes(edges=np.array([[0, 1], [1, 3], [7, 8]])))
        assert_array_equal(mesh.triangles, TRIANGLES)
        assert_array_equal(mesh.vertices, VERTICES)

    def test_other_dataset_strings_beside_the_crs_are_ignored(self, vtk_read: ModuleType) -> None:
        fields = (("dem_source", "synth.tif"), ("crs", CRS), ("dem_vertical_unit", "metre"))
        assert vtk_read.read_mesh(mesh_bytes(fields=fields)).crs == CRS

    def test_a_proj_string_crs_with_spaces_survives(self, vtk_read: ModuleType) -> None:
        crs = "+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"
        assert vtk_read.read_mesh(mesh_bytes(fields=(("crs", crs),))).crs == crs

    def test_a_mesh_written_with_an_empty_edge_vocabulary_reads_back(
        self, vtk_read: ModuleType
    ) -> None:
        """An empty vocabulary is legal; the writer emits a zero-entry string array."""
        data = mesh_bytes(vocabulary=EdgeVocabulary(properties=()))
        assert b"feature_names 1 0 string" in data
        mesh = vtk_read.read_mesh(data)
        assert_array_equal(mesh.vertices, VERTICES)
        assert_array_equal(mesh.triangles, TRIANGLES)
        assert mesh.crs == CRS
        assert mesh.land_cover_code is None


class TestRefusals:
    def test_an_ascii_mesh_is_refused_saying_how_to_write_it(self, vtk_read: ModuleType) -> None:
        with pytest.raises(ValueError, match="write the mesh without --ascii"):
            vtk_read.read_mesh(mesh_bytes(binary=False))

    def test_another_programs_file_is_refused(self, vtk_read: ModuleType) -> None:
        foreign = mesh_bytes().replace(b"\nrasputin mesh\n", b"\nvtk output\n", 1)
        with pytest.raises(ValueError, match="not a mesh written by rasputin"):
            vtk_read.read_mesh(foreign)

    def test_a_file_that_is_not_vtk_at_all_is_refused(self, vtk_read: ModuleType) -> None:
        with pytest.raises(ValueError, match="not a mesh written by rasputin"):
            vtk_read.read_mesh(b"\x89PNG\r\n\x1a\n" + bytes(64))

    def test_a_mesh_without_a_crs_is_refused(self, vtk_read: ModuleType) -> None:
        with pytest.raises(ValueError, match="no crs"):
            vtk_read.read_mesh(mesh_bytes(fields=()))

    def test_a_truncated_mesh_is_refused_not_misread(self, vtk_read: ModuleType) -> None:
        data = mesh_bytes()
        cut = data.index(b"POLYGONS") + 20
        with pytest.raises(ValueError, match="truncated"):
            vtk_read.read_mesh(data[:cut])


class TestLandCoverCodes:
    def test_codes_and_their_string_come_back_with_lines(
        self, vtk_read: ModuleType, coded: bytes
    ) -> None:
        mesh = vtk_read.read_mesh(coded)
        assert mesh.land_cover_code.dtype == np.int32
        assert mesh.land_cover_code.shape == (5,)
        assert_array_equal(mesh.land_cover_code, CODES)
        assert mesh.land_cover_codes == CODE_SYSTEM
        assert_array_equal(mesh.vertices, VERTICES)
        assert_array_equal(mesh.triangles, TRIANGLES)
        assert mesh.crs == CRS

    def test_codes_come_back_without_lines(self, vtk_read: ModuleType) -> None:
        data = mesh_bytes(triangle_codes=CODES, land_cover_codes=CODE_SYSTEM)
        assert b"LINES" not in data
        mesh = vtk_read.read_mesh(data)
        assert mesh.land_cover_code.dtype == np.int32
        assert_array_equal(mesh.land_cover_code, CODES)

    def test_a_mesh_without_codes_reads_as_none_and_empty(self, vtk_read: ModuleType) -> None:
        mesh = vtk_read.read_mesh(mesh_bytes(edges=LINES))
        assert mesh.land_cover_code is None
        assert mesh.land_cover_codes == ""

    def test_codes_without_a_string_read_with_an_empty_string(self, vtk_read: ModuleType) -> None:
        mesh = vtk_read.read_mesh(mesh_bytes(edges=LINES, triangle_codes=CODES))
        assert_array_equal(mesh.land_cover_code, CODES)
        assert mesh.land_cover_codes == ""

    def test_codes_survive_feature_arrays_and_a_fraction_table_around_them(
        self, vtk_read: ModuleType
    ) -> None:
        table = ClassFractions(
            offsets=np.array([0, 2, 3, 5, 6, 7], dtype=np.uint32),
            classes=np.array([3, 15, 33, 0, 24, 62, 7], dtype=np.uint8),
            fractions=np.array([0.25, 0.75, 1.0, 0.4, 0.6, 1.0, 1.0], dtype=np.float32),
        )
        mask = DEFAULT_VOCABULARY.mask("river", "road")
        data = mesh_bytes(
            edges=LINES, edge_mask=mask, triangle_codes=CODES,
            land_cover_codes=CODE_SYSTEM, class_fractions=table,
        )  # fmt: skip
        assert b"FIELD features 3" in data
        assert all(name.encode() in data for name in ("land_cover_offsets", "road"))
        mesh = vtk_read.read_mesh(data)
        assert_array_equal(mesh.land_cover_code, CODES)
        assert mesh.land_cover_codes == CODE_SYSTEM
        assert mesh.crs == CRS


class TestLandCoverRefusals:
    def test_a_code_array_short_of_the_cell_count_is_refused(
        self, vtk_read: ModuleType, coded: bytes
    ) -> None:
        short = coded[:-5].replace(CODE_HEADER, b"land_cover_code 1 6 int\n") + b"\n"
        with pytest.raises(ValueError, match="land_cover_code array has 6 entries for 7 cells"):
            vtk_read.read_mesh(short)

    def test_a_code_array_of_another_type_is_refused(
        self, vtk_read: ModuleType, coded: bytes
    ) -> None:
        floats = coded.replace(CODE_HEADER, b"land_cover_code 1 7 float\n")
        with pytest.raises(ValueError, match="land_cover_code array is float, not int"):
            vtk_read.read_mesh(floats)

    def test_a_file_cut_inside_the_code_array_is_refused(
        self, vtk_read: ModuleType, coded: bytes
    ) -> None:
        with pytest.raises(ValueError, match="truncated"):
            vtk_read.read_mesh(coded[:-10])
