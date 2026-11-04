"""`read_mesh` and the DEM height arrays: `MeshArrays.dem_heights`.

`docs/increments/137-dem-stats.md`, 3.6 *Reader* and test 15 of section 5.
The files are hand-made: today's `write_vtk` output with the dataset arrays
`dem_height_shares`, `dem_height_offsets`, `dem_height_first` and
`dem_height_edges_m` added to its `FieldData`, and a cell `FIELD` block
holding `dem_nodata_share` (every cell, lines first) appended, so the reader
is tested on the layout, not on the writer that will produce it.
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from height_fixtures import edges
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.io.models import HeightTable
from tin_engine.io.vtk_legacy import write_vtk
from tin_engine.io.vtk_read import read_mesh

RNG = np.random.default_rng(1370)
VERTICES = np.column_stack(
    [494000.0 + RNG.random(9) * 1200.0, 6794000.0 + RNG.random(9) * 1200.0, RNG.random(9) * 900]
)
TRIANGLES = np.array([[0, 1, 2], [2, 1, 3], [3, 4, 5], [5, 6, 7], [7, 8, 0]])
LINES = np.array([[0, 1], [1, 2]])
#: n = 2: classes 0 ... 4, edges -2.5 ... 2.5.
EDGES = edges(2, 1.0)
OFFSETS = np.array([0, 3, 3, 4, 6, 9])
FIRST = np.array([1, 0, 4, 0, 2])
_BIG = {"float": ">f4", "double": ">f8", "unsigned_int": ">u4", "int": ">i4"}


def array(name: str, kind: str, values: Any) -> bytes:
    body = np.asarray(values).astype(_BIG[kind]).tobytes()
    return f"{name} 1 {len(np.asarray(values))} {kind}\n".encode() + body + b"\n"


def mesh(
    *,
    lines: bool = True,
    shares_kind: str = "float",
    offsets: Any = OFFSETS,
    first: Any = FIRST,
    nodata_kind: str = "float",
    edges_m: Any = EDGES,
    with_arrays: bool = True,
) -> tuple[bytes, np.ndarray, np.ndarray]:
    """A binary mesh of 5 triangles (and 2 lines), and the shares and the
    triangles' NoData shares it carries."""
    edge_rows = LINES if lines else np.empty((0, 2), np.int64)
    data = write_vtk(
        VERTICES,
        triangles=TRIANGLES,
        edges=edge_rows,
        edge_masks=np.ones(len(edge_rows), np.uint32),
        vocabulary=DEFAULT_VOCABULARY,
        fields=(("crs", "EPSG:25833"),),
        binary=True,
    )
    shares = RNG.random(9).astype(np.float32)
    nodata = RNG.random(len(edge_rows) + len(TRIANGLES)).astype(np.float32)
    nodata[: len(edge_rows)] = 0.0
    if not with_arrays:
        return data, shares, nodata[len(edge_rows) :]
    found = re.search(rb"FIELD FieldData (\d+)\n", data)
    assert found is not None, data[:200]
    head = f"FIELD FieldData {int(found.group(1)) + 4}\n".encode()
    added = (
        array("dem_height_shares", shares_kind, shares)
        + array("dem_height_offsets", "unsigned_int", offsets)
        + array("dem_height_first", "int", first)
        + array("dem_height_edges_m", "double", edges_m)
    )
    data = data[: found.start()] + head + added + data[found.end() :]
    block = b"FIELD dem 1\n" + array("dem_nodata_share", nodata_kind, nodata)
    return data + block, shares, nodata[len(edge_rows) :]


class TestRead:
    @pytest.mark.parametrize("lines", [True, False])
    def test_the_table_and_the_triangle_part_of_the_nodata_share(self, lines: bool) -> None:
        data, shares, nodata = mesh(lines=lines)
        got = read_mesh(data).dem_heights
        assert isinstance(got, HeightTable)
        assert got.offsets.dtype == np.uint32 and got.first.dtype == np.int32
        assert got.shares.dtype == np.float32 and got.nodata.dtype == np.float32
        assert got.edges_m.dtype == np.float64
        assert_array_equal(got.offsets, OFFSETS)
        assert_array_equal(got.first, FIRST)
        assert_array_equal(got.shares, shares)
        assert_array_equal(got.nodata, nodata)
        assert got.edges_m.tobytes() == EDGES.tobytes()
        assert_array_equal(got.count, [3, 0, 1, 2, 3])
        assert_array_equal(read_mesh(data).triangles, TRIANGLES)

    def test_a_file_without_them_gives_none(self) -> None:
        data, _, _ = mesh(with_arrays=False)
        assert read_mesh(data).dem_heights is None


class TestRefusals:
    """Test 15: a `ValueError` naming the array."""

    @pytest.mark.parametrize(
        ("name", "change"),
        [
            ("dem_height_shares", {"shares_kind": "double"}),
            ("dem_height_offsets", {"offsets": [0, 3, 3, 4, 6, 8]}),
            ("dem_height_offsets", {"offsets": [0, 3, 2, 4, 6, 9]}),
            ("dem_height_first", {"first": [1, 0, 4, 0]}),
            ("dem_height_(first|offsets)", {"first": [1, 0, 4, 0, 3]}),
            ("dem_nodata_share", {"nodata_kind": "double"}),
            ("dem_height_edges_m", {"edges_m": EDGES[:-1]}),
            ("dem_height_edges_m", {"edges_m": EDGES[::-1]}),
        ],
        ids=[
            "shares double",
            "offsets one short",
            "offsets decreasing",
            "first one short",
            "a row past the classes",
            "nodata double",
            "edges of odd length",
            "edges decreasing",
        ],
    )
    def test_refused_naming_the_array(self, name: str, change: dict[str, Any]) -> None:
        data, _, _ = mesh(**change)
        with pytest.raises(ValueError, match=name):
            read_mesh(data)
