"""The edge strip's oracles at the Python end (increment 15f-3, PY1-PY4).

`docs/increments/15f-edge-strip.md`, "The guarantee, and its wording" (E1, E2)
and "Tests for @tester" (PY3, PY4; L8 for which oracles a path carries).
Test support only; nothing in `src_python` uses it, and nothing here imports
`tin_engine`: the points are generated from the geometry the test was given,
z comes from the DEM array, and the mesh value from the output's own
constraint lines. The producer's store and records are never read (the
computational-geometry skill: "borrow the producer's predicate, never its
records").

- `ruled_points` generates the strip's points from segments: every crossing
  of a segment with a grid line (z linear between the two nodes of that cell
  side) and, with `midpoints`, the midpoint between neighbouring crossings,
  the segment's ends counting as neighbours (z bilinear in its cell).
  Crossings closer than 1e-12 in parameter are merged; the producer's
  ulp-scale arithmetic is not reproduced (L11).
- `strip_findings` measures each point against the written constraint line
  that holds it: the nearest segment within `on_edge` metres whose parameter
  is in [0, 1], interpolated linearly between its ends' z.
- `node_findings` is `tester.md` §3D's tolerance oracle (E2): every valid DEM
  node in every closed triangle with three valid vertices, against that
  triangle's plane, recomputed from the output.
- `delaunay_violations` is §3D's constrained-Delaunay oracle, decided in
  rational arithmetic in the frame the producer's Lawson flips use,
  `LatticeFrame::at`: `(col dx, -(row dy))` rounded in float64, with
  `col = (x - x_min) / dx` and `row = (y_max - y) / dy` in float64 as
  `detail::lattice_position` computes them. That is the producer's frame
  exactly for DEM nodes and for the vertices the run was given; for a vertex
  it inserted off-node it is not recoverable from the output (see the
  function).
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from itertools import pairwise

import numpy as np
import numpy.typing as npt

F64 = npt.NDArray[np.float64]


@dataclass(frozen=True)
class Grid:
    """A north-up node lattice: node `(row, col)` at `(x_min + col dx, y_max - row dy)`."""

    x_min: float
    y_max: float
    dx: float
    dy: float
    array: F64  # (rows, cols), float64, NaN for NoData

    @property
    def rows(self) -> int:
        return int(self.array.shape[0])

    @property
    def cols(self) -> int:
        return int(self.array.shape[1])

    def lattice(self, xy: F64) -> tuple[F64, F64]:
        return (xy[:, 0] - self.x_min) / self.dx, (self.y_max - xy[:, 1]) / self.dy

    def world(self, col: F64, row: F64) -> F64:
        return np.column_stack([self.x_min + col * self.dx, self.y_max - row * self.dy])

    def bilinear(self, col: F64, row: F64) -> F64:
        """The DEM's bilinear surface at lattice `(col, row)`; NaN where a corner
        of the cell is NoData (refine's `vertex_z` rule)."""
        c0 = np.clip(np.floor(col).astype(np.int64), 0, self.cols - 2)
        r0 = np.clip(np.floor(row).astype(np.int64), 0, self.rows - 2)
        fc, fr = col - c0, row - r0
        a = self.array
        return np.asarray(
            a[r0, c0] * (1 - fc) * (1 - fr)
            + a[r0, c0 + 1] * fc * (1 - fr)
            + a[r0 + 1, c0] * (1 - fc) * fr
            + a[r0 + 1, c0 + 1] * fc * fr
        )

    def slope_bound(self) -> float:
        """An upper bound on the bilinear surface's gradient, metres per metre."""
        a = self.array
        sx = np.nanmax(np.abs(np.diff(a, axis=1))) / self.dx
        sy = np.nanmax(np.abs(np.diff(a, axis=0))) / self.dy
        return float(np.hypot(sx, sy))


def ruled_points(
    segments: Iterable[tuple[F64, F64]], grid: Grid, *, midpoints: bool = True
) -> tuple[F64, F64]:
    """The strip's points on `segments` (world end pairs): positions `(P, 2)`
    and z `(P,)`, NoData points left out."""
    cols: list[float] = []
    rows: list[float] = []
    for p0, p1 in segments:
        (c0, c1), (r0, r1) = grid.lattice(np.array([p0, p1], np.float64))
        found: list[tuple[float, float, float]] = []
        for lo, hi, axis in ((c0, c1, 0), (r0, r1, 1)):
            if lo == hi:
                continue
            first, last = sorted((lo, hi))
            for k in range(int(np.floor(first)) + 1, int(np.ceil(last))):
                t = (k - lo) / (hi - lo)
                if axis == 0:
                    found.append((t, float(k), r0 + t * (r1 - r0)))
                else:
                    found.append((t, c0 + t * (c1 - c0), float(k)))
        found.sort()
        chain = [(0.0, c0, r0)]
        for t, c, r in found:
            if 1e-12 < t < 1 - 1e-12 and t - chain[-1][0] > 1e-12:
                chain.append((t, c, r))
        chain.append((1.0, c1, r1))
        for _, c, r in chain[1:-1]:
            cols.append(c)
            rows.append(r)
        if midpoints:
            for (_, ca, ra), (_, cb, rb) in pairwise(chain):
                cols.append((ca + cb) / 2)
                rows.append((ra + rb) / 2)
    col, row = np.array(cols, np.float64), np.array(rows, np.float64)
    z = grid.bilinear(col, row)
    keep = np.isfinite(z)
    return grid.world(col[keep], row[keep]), z[keep]


def polygon_segments(rings: Iterable[Iterable[tuple[float, float]]]) -> list[tuple[F64, F64]]:
    """Every edge of every closed ring, given without the closing vertex."""
    out: list[tuple[F64, F64]] = []
    for ring in rings:
        pts = np.array(list(ring), np.float64)
        out += [(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
    return out


def edge_segments(vertices: F64, edges: npt.NDArray[np.integer]) -> list[tuple[F64, F64]]:
    return [(vertices[a, :2], vertices[b, :2]) for a, b in np.asarray(edges, np.int64)]


@dataclass(frozen=True)
class StripFindings:
    points: int
    unlocated: int
    over: int
    worst: float  # largest |z - mesh value| over located points


def strip_findings(
    xy: F64,
    z: F64,
    vertices: F64,
    vertex_z: F64,
    edges: npt.NDArray[np.integer],
    *,
    tolerance: float,
    on_edge: float,
    slack: float,
) -> StripFindings:
    """Each point against the written constraint line nearest it within
    `on_edge` (parameter in [0, 1]); over when `|z - mesh value|` exceeds
    `tolerance + slack + 1e-9 max(1, |z|)`."""
    e = np.asarray(edges, np.int64)
    origin = vertices[:, :2].min(axis=0)
    a, b = vertices[e[:, 0], :2] - origin, vertices[e[:, 1], :2] - origin
    za, zb = vertex_z[e[:, 0]], vertex_z[e[:, 1]]
    d = b - a
    length2 = (d * d).sum(axis=1)
    over = unlocated = 0
    worst = 0.0
    for lo in range(0, len(xy), 256):
        p = xy[lo : lo + 256, None, :] - origin
        sigma = ((p - a) * d).sum(axis=2) / length2
        foot = a + np.clip(sigma, 0.0, 1.0)[..., None] * d
        distance = np.hypot(*(p - foot).transpose(2, 0, 1))
        holds = (distance <= on_edge) & (sigma >= -1e-12) & (sigma <= 1 + 1e-12)
        nearest = np.where(holds, distance, np.inf).argmin(axis=1)
        located = holds.any(axis=1)
        s = np.clip(sigma[np.arange(len(nearest)), nearest], 0.0, 1.0)
        value = za[nearest] + s * (zb[nearest] - za[nearest])
        zp = z[lo : lo + 256]
        error = np.abs(zp - value)
        bound = tolerance + slack + 1e-9 * np.maximum(1.0, np.abs(zp))
        unlocated += int((~located).sum())
        over += int((located & (error > bound)).sum())
        if located.any():
            worst = max(worst, float(error[located].max()))
    return StripFindings(points=len(xy), unlocated=unlocated, over=over, worst=worst)


@dataclass(frozen=True)
class NodeFindings:
    nodes: int
    over: int
    worst: float


def node_findings(
    grid: Grid,
    vertices: F64,
    vertex_z: F64,
    valid: npt.NDArray[np.bool_],
    triangles: npt.NDArray[np.integer],
    *,
    tolerance: float,
) -> NodeFindings:
    """E2: every valid DEM node in every closed triangle (barycentric, 1e-12
    slack) with three valid vertices, within `tolerance + 1e-9 max(1, |z|)`
    of the triangle's plane."""
    tri = np.asarray(triangles, np.int64)
    tri = tri[valid[tri].all(axis=1)]
    r, c = np.indices(grid.array.shape, dtype=np.float64)
    nz = grid.array.ravel()
    ok = np.isfinite(nz)
    nxy = grid.world(c.ravel()[ok], r.ravel()[ok])
    nz = nz[ok]
    origin = vertices[:, :2].min(axis=0)
    v = vertices[:, :2] - origin
    a, b, cc = v[tri[:, 0]], v[tri[:, 1]], v[tri[:, 2]]
    two_a = (b[:, 0] - a[:, 0]) * (cc[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (cc[:, 0] - a[:, 0])
    za, zb, zc = vertex_z[tri[:, 0]], vertex_z[tri[:, 1]], vertex_z[tri[:, 2]]
    over, worst, seen = 0, 0.0, 0
    for lo in range(0, len(nxy), 256):
        x = nxy[lo : lo + 256, 0, None] - origin[0]
        y = nxy[lo : lo + 256, 1, None] - origin[1]

        def weight(p: F64, q: F64, x: F64 = x, y: F64 = y) -> F64:
            cross = (q[:, 0] - p[:, 0]) * (y - p[:, 1]) - (q[:, 1] - p[:, 1]) * (x - p[:, 0])
            return np.asarray(cross / two_a)

        wa, wb, wc = weight(b, cc), weight(cc, a), weight(a, b)
        holds = (wa >= -1e-12) & (wb >= -1e-12) & (wc >= -1e-12)
        zp = nz[lo : lo + 256, None]
        error = np.where(holds, np.abs(wa * za + wb * zb + wc * zc - zp), -np.inf)
        bound = tolerance + 1e-9 * np.maximum(1.0, np.abs(zp))
        seen += int(holds.any(axis=1).sum())
        over += int((error > bound).any(axis=1).sum())
        if holds.any():
            worst = max(worst, float(error.max()))
    return NodeFindings(nodes=seen, over=over, worst=worst)


Frame = tuple[Fraction, Fraction]


def _incircle(a: Frame, b: Frame, c: Frame, d: Frame) -> Fraction:
    """Positive exactly when `d` is strictly inside the circle through the
    counter-clockwise `a`, `b`, `c`."""
    rows = [(p[0] - d[0], p[1] - d[1]) for p in (a, b, c)]
    m = [(x, y, x * x + y * y) for x, y in rows]
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def delaunay_violations(
    grid: Grid,
    vertices: F64,
    triangles: npt.NDArray[np.integer],
    edges: npt.NDArray[np.integer],
    exact: npt.NDArray[np.bool_] | None = None,
) -> int:
    """§3D: interior edges that are not constraint edges and whose opposite
    apex lies strictly inside the other triangle's circumcircle, in the
    producer's frame. Triangles are counter-clockwise there.

    A quad whose four vertices are all `exact` (default: every vertex) is
    decided exactly. A vertex the run inserted off-node is output at a
    rounded world position, from which the producer's own frame position
    cannot be recovered; a quad with such a vertex counts only when its
    incircle determinant exceeds 1e-9 of its scale to the fourth, far above
    the ulp-scale shift that rounding can cause.
    """
    col, row = grid.lattice(vertices[:, :2])
    fx, fy = col * grid.dx, -(row * grid.dy)
    frame = [(Fraction(float(x)), Fraction(float(y))) for x, y in zip(fx, fy, strict=True)]
    sure = np.ones(len(vertices), bool) if exact is None else np.asarray(exact, bool)
    constrained = {(min(int(a), int(b)), max(int(a), int(b))) for a, b in np.asarray(edges)}
    tri = np.asarray(triangles, np.int64)
    apexes: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for t, (i, j, k) in enumerate(tri):
        for u, v, w in ((i, j, k), (j, k, i), (k, i, j)):
            apexes.setdefault((min(int(u), int(v)), max(int(u), int(v))), []).append((t, int(w)))
    bad = 0
    for key, sides in apexes.items():
        if len(sides) != 2 or key in constrained:
            continue
        (t, _), (_, w) = sides
        quad = [*(int(x) for x in tri[t]), w]
        value = _incircle(*(frame[q] for q in quad))
        if value <= 0:
            continue
        if all(sure[q] for q in quad):
            bad += 1
            continue
        d = frame[w]
        scale = max(abs(frame[q][0] - d[0]) + abs(frame[q][1] - d[1]) for q in quad[:3])
        if value > Fraction(1, 10**9) * scale**4:
            bad += 1
    return bad


def exact_vertices(grid: Grid, vertices: F64, given: F64) -> npt.NDArray[np.bool_]:
    """Vertices whose frame position the oracle computes as the producer did:
    DEM nodes, and the off-node vertices the run was given (`given`, world)."""
    col, row = grid.lattice(vertices[:, :2])
    node = (col == np.round(col)) & (row == np.round(row))
    known = {(float(x), float(y)) for x, y in np.asarray(given)[:, :2]}
    start = np.array([(float(x), float(y)) in known for x, y in vertices[:, :2]], bool)
    return np.asarray(node | start)
