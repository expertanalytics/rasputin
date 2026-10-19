"""`tin_engine.pieces` on hand-built meshes: the slice's seam chains and
heights, the seam records, and the order of the piece ids.

`docs/increments/23-basin-scale.md`, "The seam protocol", steps 1, 3 and 5,
and N17. The end-to-end suites never reach some of these branches: the
noder hands every seam edge over with its ends already in the seam pass's
`a < b` order, so the reversal of a seam pass's points never runs there;
the seam pass's z and refine's bilinear z agree at every seam vertex of the
synthetic scenes, so the heights the slice hands on are not observable in a
file; and the piece ids come out sorted in creation order on every scene.
Each test here plants the case directly, with a stand-in for the seam pass
(the slice reads only `a`, `points`, `z`, `z_a` and `z_b` from it).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from tin_engine.pieces import Slice, piece_ids, seam_records, start_slice

SEAM = 2
OUTLINE = 1

#: A 4 m square, counter-clockwise, split on its diagonal 0-2.
SQUARE_XY = np.array([(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)])
SQUARE_TRIANGLES = np.array([(0, 1, 2), (0, 2, 3)], dtype=np.int64)
#: Edge bits per triangle: bit e for edge (t[e], t[e + 1]); the diagonal
#: is edge 2 of the first triangle and edge 0 of the second.
SQUARE_FLAGS = np.array([0b111, 0b111], dtype=np.uint8)


def pass_of(table: dict[int, dict[str, Any]]) -> Callable[[int], Any]:
    """A seam pass that answers edge `r` from `table[r]`, as the core's
    `SeamOutcome` would: `a` is the end the points run from."""

    def seam_pass(r: int) -> Any:
        return SimpleNamespace(**table[r])

    return seam_pass


def square_slice(diagonal: tuple[int, int], seam_pass: Callable[[int], Any]) -> Slice:
    """The first triangle of the square, its diagonal stored as `diagonal`."""
    edges = np.array([(0, 1), (1, 2), (2, 3), (3, 0), diagonal], dtype=np.uint32)
    masks = np.array([OUTLINE] * 4 + [SEAM], dtype=np.uint32)
    keep = np.array([True, False])
    return start_slice(
        SQUARE_XY, SQUARE_TRIANGLES, SQUARE_FLAGS, edges, masks, keep, SEAM, seam_pass
    )


#: The seam pass's answer for the diagonal, oriented `a < b` by (x, y).
DIAGONAL = {
    "a": (0.0, 0.0),
    "points": [(1.0, 1.0), (3.0, 3.0)],
    "z": [11.0, 13.0],
    "z_a": 10.0,
    "z_b": 14.0,
}


class TestTheChain:
    @pytest.mark.parametrize("diagonal", [(0, 2), (2, 0)], ids=["as-the-pass", "reversed"])
    def test_runs_from_the_edges_first_end_to_its_second(self, diagonal: tuple[int, int]) -> None:
        """Stored either way round, the chain runs from the edge's first end
        through the seam pass's points in order along it, so each point keeps
        its own z; reversed, the slice turns the seam pass's answer round."""
        piece = square_slice(diagonal, pass_of({4: DIAGONAL}))
        chain = piece.chains[4]
        xy = [tuple(p) for p in piece.vertices[chain].tolist()]
        heights = [piece.z[v] for v in chain]
        forward = [(0.0, 0.0), (1.0, 1.0), (3.0, 3.0), (4.0, 4.0)]
        if diagonal == (0, 2):
            assert (xy, heights) == (forward, [10.0, 11.0, 13.0, 14.0])
        else:
            assert (xy, heights) == (forward[::-1], [14.0, 13.0, 11.0, 10.0])

    @pytest.mark.parametrize("diagonal", [(0, 2), (2, 0)], ids=["as-the-pass", "reversed"])
    def test_the_fan_is_counter_clockwise_and_covers_the_triangle(
        self, diagonal: tuple[int, int]
    ) -> None:
        piece = square_slice(diagonal, pass_of({4: DIAGONAL}))
        p = piece.vertices[piece.triangles.astype(np.int64)]
        twice = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (
            p[:, 1, 1] - p[:, 0, 1]
        ) * (p[:, 2, 0] - p[:, 0, 0])
        assert (twice > 0).all(), twice
        assert twice.sum() / 2 == 8.0  # the first triangle's area: exact at this scale

    def test_a_seam_vertex_without_data_gets_no_height(self) -> None:
        piece = square_slice((0, 2), pass_of({4: {**DIAGONAL, "z_a": None}}))
        assert 0 not in piece.z
        assert sorted(piece.z.values()) == [11.0, 13.0, 14.0]


def test_a_vertex_on_two_seams_takes_the_lower_edges_height() -> None:
    """N17: a vertex at the end of several seam edges takes the z of the
    lowest-numbered one. Vertex 2 ends seam edges 1 and 2 of one triangle;
    each seam pass gives it a different z (10 + r)."""
    triangles = np.array([(0, 1, 2)], dtype=np.int64)
    edges = np.array([(0, 1), (1, 2), (2, 0)], dtype=np.uint32)
    masks = np.array([OUTLINE, SEAM, SEAM], dtype=np.uint32)
    table = {
        1: {"a": (4.0, 0.0), "points": [(4.0, 2.0)], "z": [5.0], "z_a": 1.0, "z_b": 11.0},
        2: {"a": (0.0, 0.0), "points": [(2.0, 2.0)], "z": [6.0], "z_a": 0.0, "z_b": 12.0},
    }
    piece = start_slice(
        SQUARE_XY[:3], triangles, np.array([0b111], dtype=np.uint8), edges, masks,
        np.array([True]), SEAM, pass_of(table),
    )  # fmt: skip
    assert piece.chains[1][-1] == piece.chains[2][0] == 2
    assert piece.z[2] == 11.0
    assert sorted(piece.z.items()) == [(0, 0.0), (1, 1.0), (2, 11.0), (3, 5.0), (4, 6.0)]
    assert len(piece.triangles) == 1 + 1 + 1  # each seam edge fanned at its one point


def test_seam_records_hold_every_chain_vertex_ends_included() -> None:
    """Each record is the whole chain, both ends included, in chain order;
    a vertex without data has NaN z (`valid` decides, never the z value)."""
    piece = square_slice((0, 2), pass_of({4: DIAGONAL}))
    n = len(piece.vertices)
    z = np.arange(n, dtype=np.float64) + 100.0
    valid = np.ones(n, dtype=bool)
    valid[piece.chains[4][-1]] = False
    (record,) = seam_records((0, 1, 0), piece, piece.vertices, z, valid)
    assert (record.piece, record.edge) == ((0, 1, 0), 4)
    chain = piece.chains[4]
    assert [r[:2] for r in record.xyz] == [tuple(piece.vertices[v]) for v in chain]
    assert [r[2] for r in record.xyz[:-1]] == [z[v] for v in chain[:-1]]
    assert math.isnan(record.xyz[-1][2])


def test_piece_ids_are_sorted_whatever_the_creation_order() -> None:
    """The ids come back sorted, and each triangle's label indexes them.
    The component holding triangle 0 lies in cell (1, 0) and the other in
    cell (0, 0), so creation order (by lowest triangle) is not id order."""
    labels, ids = piece_ids(
        SQUARE_TRIANGLES, SQUARE_XY, np.array([(0, 2)], dtype=np.int64),
        lambda x, y: (1, 0) if x > y else (0, 0),
    )  # fmt: skip
    assert ids == [(0, 0, 0), (1, 0, 0)]
    assert labels.tolist() == [1, 0]


class TestFixedHeights:
    """`start_slice`'s `fixed`, keyed by start-triangulation vertex: a vertex
    in it takes that z, or none, over every seam edge's `z_a` / `z_b`; a
    vertex not in it keeps N17's rule. It is how a corner where two cut lines
    cross gets its node's height in every piece (`basin_run.corner_heights`)."""

    @staticmethod
    def sliced(keep: list[bool], fixed: dict[int, float | None]) -> Slice:
        edges = np.array([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)], dtype=np.uint32)
        masks = np.array([OUTLINE] * 4 + [SEAM], dtype=np.uint32)
        return start_slice(
            SQUARE_XY, SQUARE_TRIANGLES, SQUARE_FLAGS, edges, masks, np.array(keep), SEAM,
            pass_of({4: DIAGONAL}), fixed=fixed,
        )  # fmt: skip

    def test_a_fixed_vertex_takes_its_height_over_the_seam_passs(self) -> None:
        piece = self.sliced([True, False], {2: 99.0})
        assert piece.z[piece.chains[4][-1]] == 99.0
        assert sorted(piece.z.values()) == [10.0, 11.0, 13.0, 99.0]

    def test_a_fixed_none_leaves_the_vertex_without_a_height(self) -> None:
        piece = self.sliced([True, False], {0: None})
        assert piece.chains[4][0] == 0
        assert 0 not in piece.z
        assert sorted(piece.z.values()) == [11.0, 13.0, 14.0]

    def test_the_key_is_the_start_triangulations_vertex_not_the_slices(self) -> None:
        """Keeping triangle (0, 2, 3) renumbers start vertex 2 to slice vertex 1;
        `fixed` names it 2, and slice vertex 2 (start vertex 3, on no seam) gets
        no height."""
        piece = self.sliced([False, True], {2: 99.0})
        assert piece.vertices[1].tolist() == [4.0, 4.0]
        assert piece.z[1] == 99.0
        assert 2 not in piece.z
        assert sorted(piece.z.values()) == [10.0, 11.0, 13.0, 99.0]

    def test_no_fixed_vertex_is_the_slice_without_it(self) -> None:
        with_none = self.sliced([True, False], {})
        assert with_none.z == square_slice((0, 2), pass_of({4: DIAGONAL})).z
