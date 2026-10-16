"""Reading rasputin's own mesh back: `tin_engine.io.vtk_read.read_mesh`.

`docs/increments/51-open-core.md`, *3.1 What moves to the private package*:
the reader stays in the core, offered to plugins. `read_mesh(data:
bytes) -> MeshArrays` reads the binary legacy VTK that `rasputin mesh` writes
(header line `rasputin mesh`): the `crs` string from the dataset's field data,
POINTS and POLYGONS, skipping LINES. It refuses ASCII ("write the mesh without
--ascii"), another program's file and a file with no `crs`. The round trip
`read_mesh(write_vtk(...))` gives the same vertices, triangles and crs bit for
bit; the three refusals are in plain words.

Interface:
`MeshArrays` has `vertices` (float64, `(n, 3)`), `triangles` (int64,
`(T, 3)`) and `crs` (str); a refusal is a `ValueError`.
"""

from __future__ import annotations

import importlib
from types import ModuleType

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from tin_engine.features import DEFAULT_VOCABULARY, EdgeVocabulary
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


@pytest.fixture(scope="module")
def vtk_read() -> ModuleType:
    return importlib.import_module("tin_engine.io.vtk_read")


def mesh_bytes(
    *,
    binary: bool = True,
    fields: tuple[tuple[str, str], ...] = (("crs", CRS),),
    edges: np.ndarray | None = None,
    vocabulary: EdgeVocabulary = DEFAULT_VOCABULARY,
) -> bytes:
    lines = np.empty((0, 2), dtype=np.int64) if edges is None else edges
    return write_vtk(
        VERTICES,
        triangles=TRIANGLES,
        edges=lines,
        edge_masks=np.ones(len(lines), dtype=np.uint32) if len(lines) else np.empty(0),
        vocabulary=vocabulary,
        fields=fields,
        binary=binary,
    )


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
