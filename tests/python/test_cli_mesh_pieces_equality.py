"""DC6, the main session's claim as an equality: the decomposed mesh is one
`refine` over the whole start triangulation with the seams fanned and frozen.

`docs/increments/23-basin-scale.md`, "Tests @tester can write red", DC6, and
"The union argument, checked" item 1 (pieces are slices of one start
triangulation, in its order). On a node-only case (every domain vertex on a
node, grid-line seams crossing the outline at nodes, no features, no start
quality, no feet), the union of the piece files equals, as a set of
world-coordinate triangles with their heights, one run of

    refine(dem, start with each seam edge fanned at its seam-pass points,
           tolerance, frozen_mask=seam)

over the whole domain. The oracle builds that run itself:

- the start triangulation is the run's own, taken by a spy on
  `_core.triangulate` (K9 makes it the one global call), with its
  constraint arrays by `cli._constraint_arrays`, so no tie in the CDT can
  make the oracle's start differ from the pieces';
- each seam edge's points come from `_core.refine_seam` over the whole DEM
  (on a lattice line the strip window cannot change them, N17 and N19);
- the fans are built here, by index arithmetic, from the design's text
  ("The seam protocol" step 3), not by `pieces.py`;
- the fanned start reaches `refine` through `_core.indexed_mesh`
  (`test_core_mesh_arrays.py`; confirmed under "Settled after 23c's red step (65e3990)", item 1).

The octagon and the spy are confirmed there, item 13.

The domain is an asymmetric octagon on nodes rather than the design's
rectangle: its axis-parallel sides are the ones the partition lines cross,
so every crossing is a node, and its chamfers keep the cells from being
mirror images. Control: the same start refined without the seam pass and
without freezing differs, on a case where a node on a seam is over the
tolerance before the seam pass (checked here, so the control cannot pass
for want of work).
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from shapely.geometry import Polygon

import pieces_fixtures as pf
from pieces_fixtures import TIE_X, TIE_Y
from tin_engine import _core
from tin_engine.cli import _constraint_arrays

TOL = 1.0

#: (col, row) of the octagon's vertices, every one a node; chamfers 2, 3, 4, 5.
OCTAGON_RC = [(5, 3), (3, 5), (3, 33), (6, 36), (33, 36), (37, 32), (37, 8), (32, 3)]
OCTAGON = [(TIE_X + c * pf.H, TIE_Y - r * pf.H) for c, r in OCTAGON_RC]


class Captured:
    def __init__(self) -> None:
        self.calls: list[tuple[Any, Any]] = []


@pytest.fixture
def capture(monkeypatch: pytest.MonkeyPatch) -> Captured:
    """The run's start triangulation: `triangulate`'s argument and outcome."""
    seen = Captured()
    original = _core.triangulate

    def spy(noded: Any, *args: Any, **kwargs: Any) -> Any:
        outcome = original(noded, *args, **kwargs)
        seen.calls.append((noded, outcome))
        return outcome

    monkeypatch.setattr(_core, "triangulate", spy)
    for module in [m for k, m in list(sys.modules.items()) if k.startswith("tin_engine")]:
        if getattr(module, "triangulate", None) is original:
            monkeypatch.setattr(module, "triangulate", spy)
    return seen


def view_of(array: np.ndarray) -> Any:
    frozen = array.copy()
    frozen.setflags(write=False)
    return _core.raster_view(frozen, x_min=TIE_X, y_max=TIE_Y, delta_x=pf.H, delta_y=pf.H)


def fanned(
    mesh: Any, edges: np.ndarray, masks: np.ndarray, seam: int, view: Any
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Each seam edge `(a, b)` with seam-pass points `p1 .. pk`, in triangle
    `(a, b, c)`, becomes `(a, p1, c), (p1, p2, c), .., (pk, b, c)`; the
    constraint list gets the sub-edges with the edge's mask."""
    v = [(float(x), float(y)) for x, y in np.asarray(mesh.vertices)]
    tris = [[int(i) for i in t] for t in np.asarray(mesh.triangles)]
    cons = [int(c) for c in np.asarray(mesh.constrained_edges)]
    out_edges: list[tuple[int, int]] = []
    out_masks: list[int] = []
    for (a, b), m in zip(edges.tolist(), masks.tolist(), strict=True):
        if not m & seam:
            out_edges.append((a, b))
            out_masks.append(m)
            continue
        seam_out = _core.refine_seam(view, v[a], v[b], tolerance=TOL)
        points = [(float(x), float(y)) for x, y in np.asarray(seam_out.points)]
        if (float(seam_out.a[0]), float(seam_out.a[1])) != v[a]:
            points.reverse()
        chain = [a, *range(len(v), len(v) + len(points)), b]
        v += points
        out_edges += list(itertools.pairwise(chain))
        out_masks += [m] * (len(chain) - 1)
        for t in [t for t, tri in enumerate(tris) if a in tri and b in tri]:
            tri = tris[t]
            e = next(e for e in range(3) if {tri[e], tri[(e + 1) % 3]} == {a, b})
            p, c = tri[e], tri[(e + 2) % 3]
            run = chain if p == a else chain[::-1]
            qc, cp = (cons[t] >> ((e + 1) % 3)) & 1, (cons[t] >> ((e + 2) % 3)) & 1
            subs, flags = [], []
            for i in range(len(run) - 1):
                subs.append([run[i], run[i + 1], c])
                last, first = i == len(run) - 2, i == 0
                flags.append(1 | ((qc if last else 0) << 1) | ((cp if first else 0) << 2))
            tris[t], cons[t] = subs[0], flags[0]
            tris += subs[1:]
            cons += flags[1:]
    return (
        np.array(v, dtype=np.float64),
        np.array(tris, dtype=np.uint32),
        np.array(cons, dtype=np.uint8),
        np.array(out_edges, dtype=np.uint32),
        np.array(out_masks, dtype=np.uint32),
    )


def xyz_triangles(xy: np.ndarray, z: np.ndarray, triangles: Any) -> set[frozenset[pf.XYZ]]:
    return {
        frozenset((float(xy[i, 0]), float(xy[i, 1]), float(z[i])) for i in tri)
        for tri in np.asarray(triangles).tolist()
    }


def seam_node_errors(array: np.ndarray, start: np.ndarray, edges: Any) -> float:
    """The largest |z - lerp| at the nodes on a start seam edge, before the
    seam pass: what the seam pass has to do."""
    worst = 0.0
    for a, b in edges:
        (xa, ya), (xb, yb) = start[a], start[b]
        ca, ra = round((xa - TIE_X) / pf.H), round((TIE_Y - ya) / pf.H)
        cb, rb = round((xb - TIE_X) / pf.H), round((TIE_Y - yb) / pf.H)
        n = max(abs(cb - ca), abs(rb - ra))
        za, zb = float(array[ra, ca]), float(array[rb, cb])
        for k in range(1, n):
            r, c = ra + (rb - ra) * k // n, ca + (cb - ca) * k // n
            worst = max(worst, abs(float(array[r, c]) - (za + (zb - za) * k / n)))
    return worst


def test_the_decomposed_mesh_is_one_frozen_refine_over_the_fanned_start(
    tmp_path: Path, capture: Captured
) -> None:
    array = pf.hills()
    dem = pf.write_dem(tmp_path / "dem.tif", array)
    domain = pf.write_domain(tmp_path / "octagon.geojson", OCTAGON)
    done = pf.cut(
        tmp_path / "run",
        *("--dem", str(dem), "--domain", str(domain), "--tolerance", repr(TOL)),
        *("--pieces", "4", "--start-min-angle", "0", "--no-constraint-feet"),
    )
    assert len(done.pieces) == 4
    assert pf.problems(done, array, Polygon(OCTAGON), TOL) == []
    assert len(capture.calls) == 1, "not one global start triangulation (K9)"
    noded, outcome = capture.calls[0]
    assert outcome.ok(), outcome.message
    start = outcome.mesh
    edges, masks = _constraint_arrays(start, noded)
    seam = done.pieces[0].seam
    seam_edges = np.asarray(edges)[(np.asarray(masks) & seam) != 0]
    assert len(seam_edges) >= 4
    # Node-only: every start vertex is a node.
    start_xy = np.asarray(start.vertices)
    assert (((start_xy[:, 0] - TIE_X) / pf.H) % 1 == 0).all()
    assert (((TIE_Y - start_xy[:, 1]) / pf.H) % 1 == 0).all()

    view = view_of(array)
    v, t, c, e, m = fanned(start, np.asarray(edges), np.asarray(masks), seam, view)
    whole = _core.refine(
        view,
        _core.indexed_mesh(v, t, c),  # type: ignore[attr-defined]
        e,
        m,
        tolerance=TOL,
        min_angle_deg=0.0,
        constraint_feet=False,
        frozen_mask=seam,
    )
    assert whole.ok(), whole.message
    assert np.asarray(whole.valid).all()
    u = pf.union(done.pieces)
    assert xyz_triangles(u.points[:, :2], u.points[:, 2], u.triangles) == xyz_triangles(
        np.asarray(whole.vertices), np.asarray(whole.z), whole.triangles
    )


def test_control_without_the_seam_pass_and_freezing_the_mesh_differs(
    tmp_path: Path, capture: Captured
) -> None:
    array = pf.hills()
    dem = pf.write_dem(tmp_path / "dem.tif", array)
    domain = pf.write_domain(tmp_path / "octagon.geojson", OCTAGON)
    done = pf.cut(
        tmp_path / "run",
        *("--dem", str(dem), "--domain", str(domain), "--tolerance", repr(TOL)),
        *("--pieces", "4", "--start-min-angle", "0", "--no-constraint-feet"),
    )
    noded, outcome = capture.calls[0]
    edges, masks = _constraint_arrays(outcome.mesh, noded)
    seam = done.pieces[0].seam
    seam_edges = np.asarray(edges)[(np.asarray(masks) & seam) != 0].tolist()
    assert seam_node_errors(array, np.asarray(outcome.mesh.vertices), seam_edges) > TOL
    plain = _core.refine(
        view_of(array),
        outcome.mesh,
        edges,
        masks,
        tolerance=TOL,
        min_angle_deg=0.0,
        constraint_feet=False,
    )
    assert plain.ok(), plain.message
    u = pf.union(done.pieces)
    assert pf.triangle_set(u.points, u.triangles) != pf.triangle_set(
        np.asarray(plain.vertices), plain.triangles
    )
