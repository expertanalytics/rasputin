"""Test support for increment 16c: the land-cover oracle and the spread check.

`docs/increments/16c-landcover-labels.md`, R1 ("The oracle (tests only)") and
the invariants I1 and I2. Everything here is built from the triangles and the
input polygons, never from the producer's records: no component ids, no chosen
points. It borrows only the producer's *predicate* (a point against the input
polygons, `intersects`) and its overlap rule (smallest area, ties to the
smaller code; Default D2).

- `landcover_oracle`: each triangle's centroid tested against the polygons,
  for every triangle whose inradius exceeds `1.5 * margin` (the centroid is
  then further than `margin` from every input boundary, so the answer is
  exact). The legacy's per-centre test
  (`legacy-archive:legacy/rasputin/gml_repository.py`, `land_cover`),
  vectorised, with 0 for a centre in no polygon.
- `spread_violations`: I1. Every interior edge that is not a constraint edge
  has the same code on both sides. Arrays only, no geometry.
- `vtk_labels`: I3, I1 and I2 on one `.vtk` file, for the CLI suites.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry.base import BaseGeometry

from vtkread import VtkFile, lines_as_array, polygons_as_array

#: The noder's snap spacing at the CLI's default (`cli.DEFAULT_SNAP_SPACING`).
SNAP = 1e-3
#: R1: `margin = 2h`.
MARGIN = 2 * SNAP

Coded = Sequence[tuple[BaseGeometry, int]]


def xy_of(vertices: npt.ArrayLike) -> np.ndarray:
    return np.asarray(vertices, dtype=np.float64)[:, :2]


def inradius(vertices: npt.ArrayLike, triangles: npt.ArrayLike) -> np.ndarray:
    """`2 * area / perimeter` per triangle, in x and y."""
    xy = xy_of(vertices)
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    a, b, c = xy[tri[:, 0]], xy[tri[:, 1]], xy[tri[:, 2]]
    area = 0.5 * np.abs(
        (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    )
    perimeter = np.hypot(*(b - c).T) + np.hypot(*(c - a).T) + np.hypot(*(a - b).T)
    return np.divide(2 * area, perimeter, out=np.zeros_like(area), where=perimeter > 0)


def centroids(vertices: npt.ArrayLike, triangles: npt.ArrayLike) -> np.ndarray:
    xy = xy_of(vertices)
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    return xy[tri].mean(axis=1)


def pick(hits: Sequence[int], polygons: Coded) -> int:
    """D2: the smallest area wins, ties to the smaller code; none is 0."""
    if not hits:
        return 0
    return min((polygons[i][0].area, polygons[i][1]) for i in hits)[1]


def codes_at(points: np.ndarray, polygons: Coded) -> tuple[np.ndarray, np.ndarray]:
    """Each point's code by D2, and how many polygons it is in."""
    codes = np.zeros(len(points), dtype=np.int64)
    counts = np.zeros(len(points), dtype=np.int64)
    if not len(polygons) or not len(points):
        return codes, counts
    tree = shapely.STRtree([p for p, _ in polygons])
    found, hit = tree.query(shapely.points(points), predicate="intersects")
    per_point: dict[int, list[int]] = {}
    for i, j in zip(found.tolist(), hit.tolist(), strict=True):
        per_point.setdefault(i, []).append(j)
    for i, hits in per_point.items():
        codes[i] = pick(hits, polygons)
        counts[i] = len(hits)
    return codes, counts


@dataclass(frozen=True)
class Oracle:
    """The expected code of each triangle, and which triangles it is exact for."""

    codes: np.ndarray
    checked: np.ndarray

    def mismatches(self, produced: npt.ArrayLike) -> list[tuple[int, int, int]]:
        got = np.asarray(produced).reshape(-1)
        bad = np.flatnonzero(self.checked & (got != self.codes))
        return [(int(t), int(got[t]), int(self.codes[t])) for t in bad[:20]]


def landcover_oracle(
    vertices: npt.ArrayLike, triangles: npt.ArrayLike, polygons: Coded, margin: float = MARGIN
) -> Oracle:
    """R1's oracle: centroids against the polygons, where `r > 1.5 * margin`."""
    codes, _ = codes_at(centroids(vertices, triangles), polygons)
    checked = inradius(vertices, triangles) > 1.5 * margin
    return Oracle(codes=codes, checked=checked)


def edge_keys(pairs: npt.ArrayLike, n: int) -> np.ndarray:
    e = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    return np.minimum(e[:, 0], e[:, 1]) * n + np.maximum(e[:, 0], e[:, 1])


def spread_violations(
    triangles: npt.ArrayLike, constraints: npt.ArrayLike, codes: npt.ArrayLike
) -> list[tuple[int, int, int, int]]:
    """I1: `(t, u, code_t, code_u)` for each unconstrained interior edge whose
    two triangles differ in code."""
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    got = np.asarray(codes).reshape(-1)
    n = int(tri.max()) + 1 if len(tri) else 0
    sides = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    keys = edge_keys(sides, n)
    owner = np.tile(np.arange(len(tri)), 3)
    free = ~np.isin(keys, edge_keys(constraints, n))
    keys, owner = keys[free], owner[free]
    order = np.argsort(keys, kind="stable")
    keys, owner = keys[order], owner[order]
    pair = np.flatnonzero(keys[1:] == keys[:-1])
    t, u = owner[pair], owner[pair + 1]
    bad = np.flatnonzero(got[t] != got[u])
    return [(int(t[k]), int(u[k]), int(got[t[k]]), int(got[u[k]])) for k in bad[:20]]


def interior_edge_count(triangles: npt.ArrayLike) -> int:
    """How many edges have two triangles: the probe I1 needs to mean anything."""
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    n = int(tri.max()) + 1
    keys = edge_keys(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), n)
    _, counts = np.unique(keys, return_counts=True)
    return int((counts == 2).sum())


def triangle_areas(vertices: npt.ArrayLike, triangles: npt.ArrayLike) -> np.ndarray:
    xy = xy_of(vertices)
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    a, b, c = xy[tri[:, 0]], xy[tri[:, 1]], xy[tri[:, 2]]
    return 0.5 * np.abs(
        (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    )


def vtk_labels(
    vtk: VtkFile, polygons: Coded, margin: float = MARGIN
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """I3, I1 and I2 on one `.vtk`; returns its lines, triangles and the
    triangles' codes."""
    lines, triangles = lines_as_array(vtk), polygons_as_array(vtk)
    array = vtk.cell_array("land_cover_code")
    values = np.asarray(array.values)
    # I3: an int scalar, one value per cell, 0 on every line.
    assert (array.type_name, array.components) == ("int", 1)
    assert len(values) == vtk.cells == len(lines) + len(triangles)
    assert not values[: len(lines)].any()
    codes = values[len(lines) :]
    # I1: the spread, on the arrays alone.
    assert interior_edge_count(triangles) > 0
    assert spread_violations(triangles, lines, codes) == []
    # I2: the oracle, from the input polygons; it must reach most triangles.
    oracle = landcover_oracle(vtk.points, triangles, polygons, margin)
    assert oracle.checked.sum() > 0.9 * len(triangles)
    assert oracle.mismatches(codes) == []
    return lines, triangles, codes
