"""`tin_engine.grid_domain`: increment 12, test 11 (design R1).

A regular stride subsample of the DEM's own nodes, starting at node (0, 0),
always including the last row and column; the perimeter as one
counter-clockwise ring; coordinates bit-equal to `RasterGeometry::node`.

The design names the module and what it returns, not the functions. This
suite assumes `subsample(meta, stride) -> (xy, ring)` and
`default_stride(meta) -> int`; the handback names that as a gap.

The oracle below is written from R1's words, not from the module: the expected
node set is built with `range` and the design's float64 expression.
"""

from __future__ import annotations

import ast
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from tin_engine.io.models import RasterMeta

MODULE = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "grid_domain.py"


@pytest.fixture
def grid_domain() -> Any:
    import tin_engine.grid_domain as module  # type: ignore[import-not-found]

    return module


def meta(rows: int, cols: int, *, x_min: float = 799750.0, y_max: float = 7950250.0,
         delta_x: float = 10.0, delta_y: float = 5.0) -> RasterMeta:  # fmt: skip
    return RasterMeta(
        x_min=x_min, y_max=y_max, delta_x=delta_x, delta_y=delta_y, cols=cols, rows=rows,
        epsg=25833, nodata=None, nodata_source="absent", pixel_is_area=False,
        vertical_unit_assumed=False,
    )  # fmt: skip


def axis(n: int, stride: int) -> list[int]:
    """R1's index set along one axis: every stride-th, plus the last."""
    picked = list(range(0, n, stride))
    return picked if picked[-1] == n - 1 else [*picked, n - 1]


def expected_nodes(m: RasterMeta, stride: int) -> set[tuple[float, float]]:
    return {
        (m.x_min + c * m.delta_x, m.y_max - r * m.delta_y)
        for r in axis(m.rows, stride)
        for c in axis(m.cols, stride)
    }


def signed_area(xy: np.ndarray) -> float:
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


CASES = [
    pytest.param(4, 5, 1, id="stride-1"),
    pytest.param(9, 13, 4, id="divides-n-minus-1"),
    pytest.param(10, 14, 4, id="does-not-divide"),
    pytest.param(7, 3, 5, id="stride-beyond-one-side"),
]


class TestNodes:
    @pytest.mark.parametrize(("rows", "cols", "stride"), CASES)
    def test_counts_per_side(self, grid_domain: Any, rows: int, cols: int, stride: int) -> None:
        xy, _ = grid_domain.subsample(meta(rows, cols), stride)
        assert xy.shape == (len(axis(rows, stride)) * len(axis(cols, stride)), 2)
        assert xy.dtype == np.float64

    @pytest.mark.parametrize(("rows", "cols", "stride"), CASES)
    def test_coordinates_are_bit_equal_to_the_node_expression(
        self, grid_domain: Any, rows: int, cols: int, stride: int
    ) -> None:
        m = meta(rows, cols)
        xy, _ = grid_domain.subsample(m, stride)
        got = {(float(x), float(y)) for x, y in xy}
        assert len(got) == len(xy), "a node appears twice"
        assert got == expected_nodes(m, stride)

    def test_last_row_and_column_appear_when_stride_does_not_divide(
        self, grid_domain: Any
    ) -> None:
        m = meta(10, 14)
        xy, _ = grid_domain.subsample(m, 4)
        x_max = m.x_min + (m.cols - 1) * m.delta_x
        y_min = m.y_max - (m.rows - 1) * m.delta_y
        assert (xy[:, 0] == x_max).any()
        assert (xy[:, 1] == y_min).any()
        assert xy[:, 0].max() == x_max and xy[:, 1].min() == y_min

    def test_non_integer_spacing_at_utm_scale(self, grid_domain: Any) -> None:
        """A spacing with no exact binary form: any other spelling of the
        expression (`x_min + c * dx` accumulated, or `x_max - ...`) lands ulps off."""
        m = meta(31, 29, x_min=432109.87, y_max=6789012.34, delta_x=0.1, delta_y=0.3)
        xy, _ = grid_domain.subsample(m, 3)
        assert {(float(x), float(y)) for x, y in xy} == expected_nodes(m, 3)


class TestRing:
    @pytest.mark.parametrize(("rows", "cols", "stride"), CASES)
    def test_is_the_perimeter_counter_clockwise_and_implicitly_closed(
        self, grid_domain: Any, rows: int, cols: int, stride: int
    ) -> None:
        m = meta(rows, cols)
        xy, ring = grid_domain.subsample(m, stride)
        ring = np.asarray(ring)
        nr, nc = len(axis(rows, stride)), len(axis(cols, stride))
        assert len(ring) == 2 * (nr + nc) - 4
        assert len(set(ring.tolist())) == len(ring), "closed implicitly: no repeat"
        assert signed_area(xy[ring]) > 0.0

        x_max = m.x_min + (m.cols - 1) * m.delta_x
        y_min = m.y_max - (m.rows - 1) * m.delta_y
        on_extent = (
            (xy[ring, 0] == m.x_min) | (xy[ring, 0] == x_max)
            | (xy[ring, 1] == m.y_max) | (xy[ring, 1] == y_min)
        )  # fmt: skip
        assert on_extent.all()
        # Every perimeter node is on the ring, and nothing interior is.
        perimeter = {
            p for p in expected_nodes(m, stride)
            if p[0] in (m.x_min, x_max) or p[1] in (m.y_max, y_min)
        }  # fmt: skip
        assert {(float(x), float(y)) for x, y in xy[ring]} == perimeter

    def test_consecutive_ring_vertices_are_one_step_apart(self, grid_domain: Any) -> None:
        """The ring walks the perimeter in order, never jumping across the tile."""
        m = meta(10, 14)
        xy, ring = grid_domain.subsample(m, 4)
        walk = xy[np.asarray(ring)]
        steps = np.abs(np.diff(np.vstack([walk, walk[:1]]), axis=0))
        assert ((steps[:, 0] == 0) | (steps[:, 1] == 0)).all()
        assert (steps.max(axis=1) <= 4 * max(m.delta_x, m.delta_y)).all()


class TestDefaultStride:
    @pytest.mark.parametrize(
        ("rows", "cols"),
        [(2, 2), (3, 4), (256, 10), (257, 10), (10, 5051), (5051, 5051), (511, 512)],
    )
    def test_formula_and_at_most_256_per_side(
        self, grid_domain: Any, rows: int, cols: int
    ) -> None:
        stride = grid_domain.default_stride(meta(rows, cols))
        assert stride == max(1, math.ceil((max(rows, cols) - 1) / 255))
        assert len(axis(max(rows, cols), stride)) <= 256

    def test_the_real_fixture_gets_stride_20(self, grid_domain: Any) -> None:
        assert grid_domain.default_stride(meta(5051, 5051)) == 20


def test_imports_neither_core_nor_the_reader() -> None:
    """R1: a pure module. Read from source, so it holds without importing it."""
    tree = ast.parse(MODULE.read_text(encoding="utf-8"))
    imported: set[str] = set()
    for stmt in ast.walk(tree):
        if isinstance(stmt, ast.Import):
            imported |= {alias.name for alias in stmt.names}
        elif isinstance(stmt, ast.ImportFrom):
            imported.add(f"{'.' * stmt.level}{stmt.module or ''}")
            imported |= {f"{stmt.module}.{alias.name}" for alias in stmt.names}
    assert not {name for name in imported if "_core" in name or "geotiff" in name}
