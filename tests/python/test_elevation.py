"""`tin_engine.elevation.trim`: increment 12, test 12 (design R3, the user's U1 (a)).

After sampling, drop every triangle with a vertex whose `valid` is false, every
constraint edge with such an endpoint, and every vertex no remaining triangle
uses; renumber. Pure numpy, no `_core`.

The design gives the call as `elevation.trim(mesh arrays, edges, z, valid)`
and does not fix the parameter names, the edge-mask handling or the return
type. This suite assumes keyword parameters `vertices, triangles, edges,
edge_masks, z, valid` and a result with attributes `vertices` (M, 3),
`triangles`, `edges`, `edge_masks` and `dropped`. The handback names it.

Renumbering is checked through coordinates, never through indices, so any
consistent renumbering passes and any inconsistent one fails.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

MODULE = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "elevation.py"

# A 3 x 3 node grid, row-major, split along one diagonal per quad: 8 triangles.
#   6 7 8
#   3 4 5
#   0 1 2
XY = np.array([(c * 10.0, r * 10.0) for r in range(3) for c in range(3)], dtype=np.float64)
TRIANGLES = np.array(
    [(v, v + 1, v + 4) for v in (0, 1, 3, 4)] + [(v, v + 4, v + 3) for v in (0, 1, 3, 4)],
    dtype=np.uint32,
)
# The perimeter, with a distinct feature mask per edge so a mask left behind
# on the wrong edge is visible.
EDGES = np.array([(0, 1), (1, 2), (2, 5), (5, 8), (8, 7), (7, 6), (6, 3), (3, 0)], dtype=np.uint32)
MASKS = np.array([1, 2, 4, 8, 16, 32, 64, 128], dtype=np.uint32)
Z = 100.0 + np.arange(9, dtype=np.float64)


@pytest.fixture
def trim() -> Any:
    from tin_engine.elevation import trim  # type: ignore[import-not-found]

    return trim


def run(trim: Any, valid: np.ndarray) -> Any:
    return trim(
        vertices=XY, triangles=TRIANGLES, edges=EDGES, edge_masks=MASKS, z=Z, valid=valid
    )


def valid_except(*invalid: int) -> np.ndarray:
    flags = np.ones(9, dtype=bool)
    flags[list(invalid)] = False
    return flags


def triangle_set(vertices: np.ndarray, triangles: np.ndarray) -> set[frozenset[tuple[float, ...]]]:
    return {frozenset(tuple(float(v) for v in vertices[i]) for i in t) for t in triangles}


def edge_set(vertices: np.ndarray, edges: np.ndarray, masks: np.ndarray) -> set[Any]:
    return {
        (frozenset(tuple(float(v) for v in vertices[i]) for i in e), int(m))
        for e, m in zip(edges, masks, strict=True)
    }


XYZ = np.column_stack([XY, Z])


def expected_after_dropping(*invalid: int) -> tuple[set[Any], set[Any]]:
    keep_t = [t for t in TRIANGLES if not set(t.tolist()) & set(invalid)]
    keep_e = [
        (e, m) for e, m in zip(EDGES, MASKS, strict=True) if not set(e.tolist()) & set(invalid)
    ]
    tris = triangle_set(XYZ, np.array(keep_t).reshape(-1, 3))
    edges = edge_set(
        XYZ,
        np.array([e for e, _ in keep_e]).reshape(-1, 2),
        np.array([m for _, m in keep_e]),
    )
    return tris, edges


class TestDrops:
    def test_a_corner_without_data_takes_its_triangles_and_edges(self, trim: Any) -> None:
        out = run(trim, valid_except(0))
        tris, edges = expected_after_dropping(0)
        assert out.vertices.shape == (8, 3)
        assert triangle_set(out.vertices, out.triangles) == tris
        assert edge_set(out.vertices, out.edges, out.edge_masks) == edges
        assert len(out.triangles) == 6 and len(out.edges) == 6
        assert out.dropped == 1

    def test_z_travels_with_its_vertex(self, trim: Any) -> None:
        out = run(trim, valid_except(0))
        for x, y, z in out.vertices:
            (index,) = np.flatnonzero((XY[:, 0] == x) & (XY[:, 1] == y))
            assert z == Z[index]

    def test_vertices_no_remaining_triangle_uses_are_removed(self, trim: Any) -> None:
        """Losing the centre leaves only (1, 2, 5) and (3, 7, 6); corners 0 and 8
        are valid but orphaned, and go too (R3's third bullet)."""
        out = run(trim, valid_except(4))
        tris, edges = expected_after_dropping(4)
        assert triangle_set(out.vertices, out.triangles) == tris
        assert out.vertices.shape == (6, 3)
        used = {tuple(float(v) for v in XYZ[i]) for i in (1, 2, 5, 3, 7, 6)}
        assert {tuple(float(v) for v in p) for p in out.vertices} == used
        # Perimeter edges touching an orphaned corner cannot survive either:
        # their endpoint is gone from the vertex block.
        assert all(int(i) < len(out.vertices) for e in out.edges for i in e)
        assert edge_set(out.vertices, out.edges, out.edge_masks) == {
            e for e in edges if all(p in used for p in e[0])
        }

    def test_indices_are_in_range_and_dense(self, trim: Any) -> None:
        out = run(trim, valid_except(0, 8))
        referenced = set(out.triangles.ravel().tolist())
        assert referenced == set(range(len(out.vertices)))

    def test_no_invalid_z_value_survives(self, trim: Any) -> None:
        z = Z.copy()
        z[[2, 6]] = np.nan  # what a sampler might leave in an invalid slot
        out = trim(
            vertices=XY, triangles=TRIANGLES, edges=EDGES, edge_masks=MASKS,
            z=z, valid=valid_except(2, 6),
        )  # fmt: skip
        assert np.isfinite(out.vertices).all()


class TestLimits:
    def test_all_valid_changes_nothing(self, trim: Any) -> None:
        out = run(trim, valid_except())
        assert_array_equal(out.vertices, XYZ)
        assert_array_equal(out.triangles, TRIANGLES)
        assert_array_equal(out.edges, EDGES)
        assert_array_equal(out.edge_masks, MASKS)
        assert out.dropped == 0

    def test_all_invalid_leaves_nothing(self, trim: Any) -> None:
        out = run(trim, np.zeros(9, dtype=bool))
        assert out.vertices.shape == (0, 3)
        assert len(out.triangles) == 0
        assert len(out.edges) == 0 and len(out.edge_masks) == 0
        assert out.dropped == 9


def test_imports_no_core() -> None:
    """R3: pure numpy. Read from source, so it holds without importing it."""
    tree = ast.parse(MODULE.read_text(encoding="utf-8"))
    names = [
        alias.name
        for stmt in ast.walk(tree)
        if isinstance(stmt, ast.Import | ast.ImportFrom)
        for alias in stmt.names
    ] + [stmt.module or "" for stmt in ast.walk(tree) if isinstance(stmt, ast.ImportFrom)]
    assert not [name for name in names if "_core" in name]
