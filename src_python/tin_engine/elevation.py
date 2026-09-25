"""Joining sampled z to a 2-D mesh, and removing what has none.

Increment 12, R3 and the user's U1 (a). After sampling, every triangle with a
vertex whose ``valid`` is false goes, every constraint edge with such an
endpoint goes, and every vertex no remaining triangle uses goes; the rest is
renumbered. The mesh then has a ragged edge or a hole where the DEM has no
data, instead of NaN in its points (which breaks ParaView's bounds) or a refusal
(which would refuse the project's only real DEM).

Pure numpy, no ``_core``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True, slots=True)
class Trimmed:
    """A 3-D mesh with every vertex carrying a sampled z.

    ``dropped`` counts the vertices that had no data. Valid vertices orphaned by
    the triangles around them are removed too, but are not counted: they had
    data, and the count is what the output file reports as "without data".
    """

    vertices: npt.NDArray[np.float64]
    triangles: npt.NDArray[np.uint32]
    edges: npt.NDArray[np.uint32]
    edge_masks: npt.NDArray[np.uint32]
    dropped: int


def trim(
    *,
    vertices: npt.NDArray[np.float64],
    triangles: npt.ArrayLike,
    edges: npt.ArrayLike,
    edge_masks: npt.ArrayLike,
    z: npt.NDArray[np.float64],
    valid: npt.NDArray[np.bool_],
) -> Trimmed:
    """Drop what has no z, renumber, and return ``(M, 3)`` vertices."""
    valid = np.asarray(valid, dtype=bool)
    tris = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    pairs = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    masks = np.asarray(edge_masks, dtype=np.uint32).reshape(-1)

    tris = tris[valid[tris].all(axis=1)]
    keep = np.zeros(len(valid), dtype=bool)
    keep[tris.ravel()] = True
    new_index = np.cumsum(keep) - 1

    edge_ok = keep[pairs].all(axis=1)
    xyz = np.column_stack([np.asarray(vertices, dtype=np.float64)[:, :2], z])[keep]
    return Trimmed(
        vertices=xyz.reshape(-1, 3),
        triangles=new_index[tris].astype(np.uint32).reshape(-1, 3),
        edges=new_index[pairs[edge_ok]].astype(np.uint32).reshape(-1, 2),
        edge_masks=masks[edge_ok],
        dropped=int(np.count_nonzero(~valid)),
    )
