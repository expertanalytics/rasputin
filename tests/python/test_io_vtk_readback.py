"""Increment 13's files, read back by the readers ParaView uses. Needs the `viewer` extra.

`13-bundled-mesh.md` ruling 10: the defect this increment fixes was found only
by opening the output in the real reader. Two of its hazards live *inside* the
reader and no parser of our own can see them:

- **finding 1**: VTK orders cells lines-first whatever the file says, and
  matches cell data by that index, so a triangles-first file hands every value
  to the wrong cell with no error. Only the reader that reorders can show it.
- **finding 2**: a cell array short of the cell count makes a reader return an
  empty dataset.

So this suite writes with `write_vtk` (and once through the CLI), reads with
`vtkPDataSetReader` -- what ParaView uses for `.vtk` -- and with
`vtkPolyDataReader` and `vtkDataSetReader`, in both encodings, and asserts on
what the reader hands back.

It skips when `vtk` is absent. CI runs it in its own step after installing
`.[dev,viewer]`, so a leg whose wheel disappears fails in the install and never
reaches the skip.

Every read also asserts that the reader raised no `ErrorEvent`. A VTK reader
reports a malformed file by logging and returning an empty or partial
dataset, never by raising, so a test that only counted cells could pass on the
wreckage of a file the reader refused halfway.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from typer.testing import CliRunner

vtk = pytest.importorskip("vtk")
numpy_support = pytest.importorskip("vtk.util.numpy_support")

from tin_engine.features import DEFAULT_VOCABULARY, EdgeProperty, EdgeVocabulary  # noqa: E402

VTK_LINE = 3
VTK_TRIANGLE = 5

EASTING = 430_000.001
NORTHING = 6_900_000.001

_RNG = np.random.default_rng(1313)
VERTICES = np.vstack(
    [
        [[EASTING, NORTHING, 12.345_678_901_234_567], [430_010.0, NORTHING, 13.5]],
        np.column_stack(
            [
                _RNG.uniform(430_000.0, 440_000.0, 4),
                _RNG.uniform(6_900_000.0, 6_910_000.0, 4),
                _RNG.uniform(-50.0, 2_500.0, 4),
            ]
        ),
    ]
).astype(np.float64)

TRIANGLES = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 2], [2, 4, 5]], dtype=np.uint32)
EDGES = np.array([[0, 1], [1, 2], [2, 3]], dtype=np.uint32)

#: Distinct non-zero masks on the ends and 0 in the middle. If VTK hands the
#: triangles' 0s to the lines, or the lines' masks to the triangles, both ends
#: show it.
MASKS = np.array([0b001, 0b000, 0b110], dtype=np.uint32)

#: Deliberately awkward: a space, a `%`, and a literal `%20` that must not come
#: back as a space.
CRS = "+proj=utm +zone=33 +units=m 50% a%20b"
ELEVATION_TEXT = "none (z=0, --flat)"

READERS = {
    "PDataSet": "vtkPDataSetReader",
    "PolyData": "vtkPolyDataReader",
    "DataSet": "vtkDataSetReader",
}

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})


def read_back(path: Path, reader_name: str) -> object:
    """Read `path` with the named VTK reader, refusing any reported error."""
    reader = getattr(vtk, reader_name)()
    errors: list[str] = []
    reader.AddObserver("ErrorEvent", lambda _obj, _event: errors.append(reader_name))
    reader.SetFileName(str(path))
    reader.Update()
    assert not errors, f"{reader_name} reported an error reading {path.name}"
    output = reader.GetOutput()
    assert output is not None and output.IsA("vtkPolyData"), f"{reader_name}: not PolyData"
    return output


Writer = Callable[..., object]


@pytest.fixture(params=[False, True], ids=["ascii", "binary"])
def binary(request: pytest.FixtureRequest) -> bool:
    return bool(request.param)


@pytest.fixture(params=list(READERS.values()), ids=list(READERS))
def reader_name(request: pytest.FixtureRequest) -> str:
    return str(request.param)


@pytest.fixture
def load(tmp_path: Path, binary: bool, reader_name: str) -> Writer:
    """Write with `write_vtk` (overrides allowed) and read back with the reader under test."""

    def _load(**overrides: object) -> object:
        kwargs: dict[str, object] = {
            "triangles": TRIANGLES,
            "edges": EDGES,
            "edge_masks": MASKS,
            "vocabulary": DEFAULT_VOCABULARY,
            "fields": (("crs", CRS), ("elevation_source", ELEVATION_TEXT)),
            "binary": binary,
        }
        kwargs.update(overrides)
        from tin_engine.io.vtk_legacy import write_vtk

        path = tmp_path / "mesh.vtk"
        path.write_bytes(write_vtk(VERTICES, **kwargs))  # type: ignore[arg-type]
        return read_back(path, reader_name)

    return _load


@pytest.fixture
def poly(load: Writer) -> object:
    return load()


def cell_ids(poly: object, i: int) -> list[int]:
    ids = vtk.vtkIdList()
    poly.GetCellPoints(i, ids)  # type: ignore[attr-defined]
    return [ids.GetId(k) for k in range(ids.GetNumberOfIds())]


def cell_array(poly: object, name: str) -> np.ndarray:
    array = poly.GetCellData().GetArray(name)  # type: ignore[attr-defined]
    assert array is not None, f"no cell array {name!r}"
    return np.asarray(numpy_support.vtk_to_numpy(array))


def strings(poly: object, name: str) -> list[str]:
    array = poly.GetFieldData().GetAbstractArray(name)  # type: ignore[attr-defined]
    assert array is not None, f"no field data array {name!r}"
    return [array.GetValue(k) for k in range(array.GetNumberOfValues())]


class TestPoints:
    """Ruling 2, where it failed for PLY: inside the reader."""

    def test_they_are_float64(self, poly: object) -> None:
        points = numpy_support.vtk_to_numpy(poly.GetPoints().GetData())  # type: ignore[attr-defined]
        assert points.dtype == np.float64

    def test_they_are_bit_exact(self, poly: object) -> None:
        # In binary this is also the big-endian check: a little-endian body
        # read as big-endian is not these numbers.
        points = numpy_support.vtk_to_numpy(poly.GetPoints().GetData())  # type: ignore[attr-defined]
        assert_array_equal(points, VERTICES)


class TestCells:
    """Rulings 4 and 5, seen through the reader that reorders (finding 1)."""

    def test_the_cell_count_is_lines_plus_triangles(self, poly: object) -> None:
        assert poly.GetNumberOfLines() == len(EDGES)  # type: ignore[attr-defined]
        assert poly.GetNumberOfPolys() == len(TRIANGLES)  # type: ignore[attr-defined]
        assert poly.GetNumberOfCells() == len(EDGES) + len(TRIANGLES)  # type: ignore[attr-defined]

    def test_the_first_e_cells_are_the_edges_in_order(self, poly: object) -> None:
        for i, edge in enumerate(EDGES):
            assert poly.GetCellType(i) == VTK_LINE, i  # type: ignore[attr-defined]
            assert cell_ids(poly, i) == edge.tolist(), i

    def test_the_remaining_cells_are_the_triangles_in_order(self, poly: object) -> None:
        for t, triangle in enumerate(TRIANGLES):
            i = len(EDGES) + t
            assert poly.GetCellType(i) == VTK_TRIANGLE, i  # type: ignore[attr-defined]
            assert cell_ids(poly, i) == triangle.tolist(), i

    def test_each_line_carries_its_own_mask(self, poly: object) -> None:
        # The semantic half of the invariant-critical suite. A triangles-first
        # file with data in the same order reads back with a triangle's 0 on
        # cell 0, a line.
        mask = cell_array(poly, "feature_mask")
        assert mask[: len(EDGES)].tolist() == MASKS.tolist()
        assert mask[len(EDGES) :].tolist() == [0] * len(TRIANGLES)

    def test_every_cell_array_covers_every_cell(self, poly: object) -> None:
        # Finding 2's partial array, from the reader's side.
        data = poly.GetCellData()  # type: ignore[attr-defined]
        cells = poly.GetNumberOfCells()  # type: ignore[attr-defined]
        assert cells > 0, "an empty dataset is how a reader refuses a partial array"
        for k in range(data.GetNumberOfArrays()):
            array = data.GetAbstractArray(k)
            assert array.GetNumberOfTuples() == cells, array.GetName()


class TestPerFeatureArrays:
    """Ruling 6's arrays for people, through every reader (finding 4)."""

    def test_the_arrays_present_are_the_properties_that_occur(self, poly: object) -> None:
        data = poly.GetCellData()  # type: ignore[attr-defined]
        names = {data.GetAbstractArray(k).GetName() for k in range(data.GetNumberOfArrays())}
        assert names == {"feature_mask", "river", "road", "railway"}

    def test_each_is_the_bits_state_on_lines_and_0_on_triangles(self, poly: object) -> None:
        zeros = [0] * len(TRIANGLES)
        assert cell_array(poly, "river").tolist() == [1, 0, 0, *zeros]
        assert cell_array(poly, "road").tolist() == [0, 0, 1, *zeros]
        assert cell_array(poly, "railway").tolist() == [0, 0, 1, *zeros]

    def test_feature_mask_is_the_active_scalar(self, poly: object) -> None:
        scalars = poly.GetCellData().GetScalars()  # type: ignore[attr-defined]
        assert scalars is not None and scalars.GetName() == "feature_mask"


class TestFieldData:
    """The vocabulary table, the fingerprint, CRS and elevation, as the reader sees them."""

    def test_the_table_is_the_vocabulary_sorted_by_bit(self, poly: object) -> None:
        pairs = sorted((p.bit, p.name) for p in DEFAULT_VOCABULARY.properties)
        bits = poly.GetFieldData().GetArray("feature_bits")  # type: ignore[attr-defined]
        assert bits is not None
        assert numpy_support.vtk_to_numpy(bits).tolist() == [b for b, _ in pairs]
        assert strings(poly, "feature_names") == [n for _, n in pairs]

    def test_the_fingerprint(self, poly: object) -> None:
        assert strings(poly, "feature_vocabulary") == [DEFAULT_VOCABULARY.fingerprint()]

    def test_crs_with_space_and_percent_round_trips(self, poly: object) -> None:
        assert strings(poly, "crs") == [CRS]

    def test_elevation_round_trips(self, poly: object) -> None:
        assert strings(poly, "elevation_source") == [ELEVATION_TEXT]

    def test_an_out_of_order_vocabulary_reads_back_sorted(self, load: Writer) -> None:
        vocabulary = EdgeVocabulary(
            properties=(EdgeProperty(name="road", bit=1), EdgeProperty(name="river", bit=0))
        )
        poly = load(vocabulary=vocabulary, edge_masks=MASKS & 0b11)
        assert strings(poly, "feature_names") == ["river", "road"]
        assert strings(poly, "feature_vocabulary") == [vocabulary.fingerprint()]


class TestEmptyCases:
    """`LINES 0 0`, and a mesh whose edges carry no bits."""

    def test_zero_edges_loads_the_triangles_intact(self, load: Writer) -> None:
        poly = load(edges=np.zeros((0, 2), dtype=np.uint32), edge_masks=np.zeros(0, np.uint32))
        assert poly.GetNumberOfLines() == 0  # type: ignore[attr-defined]
        assert poly.GetNumberOfCells() == len(TRIANGLES)  # type: ignore[attr-defined]
        for t, triangle in enumerate(TRIANGLES):
            assert cell_ids(poly, t) == triangle.tolist()
        assert cell_array(poly, "feature_mask").tolist() == [0] * len(TRIANGLES)
        assert strings(poly, "feature_vocabulary") == [DEFAULT_VOCABULARY.fingerprint()]

    def test_no_bits_set_loads_with_feature_mask_only(self, load: Writer) -> None:
        poly = load(edge_masks=np.zeros(len(EDGES), dtype=np.uint32))
        data = poly.GetCellData()  # type: ignore[attr-defined]
        assert data.GetNumberOfArrays() == 1
        assert cell_array(poly, "feature_mask").tolist() == [0] * (len(EDGES) + len(TRIANGLES))


class TestThroughTheCommand:
    """`rasputin mesh road-crosses-river --flat --out x.vtk`, opened by ParaView's reader."""

    FEATURED = "road-crosses-river"

    @pytest.fixture
    def written(self, tmp_path: Path, binary: bool) -> Path:
        from tin_engine.cli import app

        out = tmp_path / "rr.vtk"
        args = ["mesh", self.FEATURED, "--flat", "--crs", "EPSG:25833", "--out", str(out)]
        result = runner.invoke(app, [*args, "--binary" if binary else "--ascii"])
        assert result.exit_code == 0, result.output
        return out

    def test_it_is_the_fixtures_mesh_with_its_constraints(self, written: Path) -> None:
        poly = read_back(written, "vtkPDataSetReader")
        assert poly.GetNumberOfPoints() == 9  # type: ignore[attr-defined]
        assert poly.GetNumberOfPolys() == 12  # type: ignore[attr-defined]
        assert poly.GetNumberOfLines() == 8  # type: ignore[attr-defined]

    def test_the_river_and_road_arrays_come_through(self, written: Path) -> None:
        poly = read_back(written, "vtkPDataSetReader")
        lines = poly.GetNumberOfLines()  # type: ignore[attr-defined]
        river, road = cell_array(poly, "river"), cell_array(poly, "road")
        assert river[:lines].any() and road[:lines].any()
        assert not river[lines:].any() and not road[lines:].any()

    def test_the_metadata_comes_through(self, written: Path) -> None:
        poly = read_back(written, "vtkPDataSetReader")
        assert strings(poly, "crs") == ["EPSG:25833"]
        assert strings(poly, "elevation_source") == [ELEVATION_TEXT]
        assert strings(poly, "feature_vocabulary") == [DEFAULT_VOCABULARY.fingerprint()]


class TestPointElevation:
    """Increment 12 amendment, inside the reader: Color By -> elevation is the heights."""

    def test_elevation_is_a_float64_point_array_equal_to_z(self, poly: object) -> None:
        array = poly.GetPointData().GetArray("elevation")  # type: ignore[attr-defined]
        assert array is not None, "no point array 'elevation'"
        values = numpy_support.vtk_to_numpy(array)
        assert values.dtype == np.float64
        assert_array_equal(values, VERTICES[:, 2])
