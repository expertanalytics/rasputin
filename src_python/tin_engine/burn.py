"""The mapped reach moved onto the DEM's valley floor and burnt in (increment 29, PR 2).

`docs/increments/29-nve-reference-catchments.md`, "Following the river in the
DEM". The reach is resampled at the node spacing; each point takes the
lowest node of its cross-section (step 2); the chosen nodes are joined into
an 8-connected chain and made taut; the direction is checked against the
DEM; the chain is lowered only where it does not already fall (step 4), and
extended along the raw valley floor past a lowered end (step 5). Numpy only:
the window's array in, a burnt copy and the :class:`GaugePath` out.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import pairwise

import numpy as np
import numpy.typing as npt

from tin_engine.gauge import Reach
from tin_engine.io.models import DemTile, RasterMeta, valid_mask

#: Metres each chain node falls below the one before it, at least (step 4).
DROP_M = 0.001
#: Metres of arc length the end's extension may run (step 5).
END_CAP_M = 500.0
#: Metres of slack on the corridor's and the cross-section's bounds (step 2).
SLACK_M = 1e-6
#: Metres the chain's last tenth may lie above its first (step 3).
DIRECTION_M = 2.0
_NEIGHBOURS = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if (dr, dc) != (0, 0)]


class BurnRefusal(ValueError):  # noqa: N818 -- the record's name, pinned by the suite
    """The burn cannot place the gauge: a gap on the chain, the placed position
    outside the window, or no node with data near the reach."""


@dataclass(frozen=True, slots=True)
class GaugePath:
    """The burnt chain, (row, column) in the window, downstream, extension
    included; the placed node's index; each node's metres from it (negative
    upstream); and what the burn found."""

    chain: npt.NDArray[np.int64]
    placed: int
    arc: npt.NDArray[np.float64]
    end_extended_m: float
    end_closed: bool
    direction_ok: bool
    lowered_nodes: int
    lowered_max_m: float
    node_offset_m: float
    dropped_m: float


def _floor_node(
    z: npt.NDArray[np.floating],
    ok: npt.NDArray[np.bool_],
    xy: tuple[float, float],
    t: tuple[float, float],
    spacing: tuple[float, ...],
    step: float,
    corridor: float,
) -> tuple[int, int]:
    """The least-elevation node of the cross-section at `xy` (direction `t`);
    with none, the node nearest `xy`. Ties: nearer, then smaller (row, col)."""
    x0, y1, dx, dy = spacing
    reach = max(corridor, math.hypot(dx, dy))
    r0, r1 = math.floor((y1 - xy[1] - reach) / dy), math.ceil((y1 - xy[1] + reach) / dy)
    c0, c1 = math.floor((xy[0] - reach - x0) / dx), math.ceil((xy[0] + reach - x0) / dx)
    r, c = np.mgrid[
        max(r0, 0) : min(r1, z.shape[0] - 1) + 1, max(c0, 0) : min(c1, z.shape[1] - 1) + 1
    ]
    r, c = r.ravel(), c.ravel()
    keep = ok[r, c]
    r, c = r[keep], c[keep]
    ex, ey = x0 + c * dx - xy[0], y1 - r * dy - xy[1]
    dist = np.hypot(ex, ey)
    inside = (dist <= corridor + SLACK_M) & (np.abs(ex * t[0] + ey * t[1]) <= step / 2 + SLACK_M)
    if inside.any():
        order = np.lexsort((c[inside], r[inside], dist[inside], z[r[inside], c[inside]]))
        return int(r[inside][order[0]]), int(c[inside][order[0]])
    if r.size == 0:
        raise BurnRefusal(f"no DEM node with data near the reach at {xy}")
    k = np.lexsort((c, r, dist))[0]
    return int(r[k]), int(c[k])


def _join(nodes: list[tuple[int, int]]) -> tuple[list[tuple[int, int]], list[int]]:
    """Consecutive nodes joined by straight 8-connected runs, and the index in
    the joined chain of each given node."""
    chain, where = [nodes[0]], [0]
    for (ra, ca), (rb, cb) in pairwise(nodes):
        n = max(abs(rb - ra), abs(cb - ca))
        for i in range(1, n + 1):
            chain.append(
                (ra + math.floor((rb - ra) * i / n + 0.5), ca + math.floor((cb - ca) * i / n + 0.5))
            )
        where.append(len(chain) - 1)
    return chain, where


def _taut(chain: list[tuple[int, int]]) -> tuple[list[tuple[int, int]], list[int]]:
    """Step 2's taut pass: the kept nodes, and each given node's kept index
    (a dropped node maps to the kept node the cut starts from)."""
    a = np.array(chain)
    kept: list[tuple[int, int]] = []
    to = [0] * len(chain)
    k = 0
    while k < len(chain):
        to[k] = len(kept)
        kept.append(chain[k])
        near = np.nonzero(np.abs(a[k + 1 :] - a[k]).max(axis=1) <= 1)[0] + k + 1
        j = int(near[-1]) if near.size else k + 1
        if j > k + 1:
            for i in range(k + 1, j):
                to[i] = len(kept) - 1
            if chain[j] == chain[k]:
                to[j] = len(kept) - 1
                j += 1
        k = max(j, k + 1)
    return kept, to


def _steps(chain: list[tuple[int, int]], m: RasterMeta) -> list[float]:
    """The metres of each step along the chain."""
    return [
        math.hypot((a - c) * m.delta_y, (b - d) * m.delta_x) for (a, b), (c, d) in pairwise(chain)
    ]


def burn_reach(window: DemTile, reach: Reach) -> tuple[DemTile, GaugePath]:
    """The window with the reach burnt in (a copy, same `meta`), and its path.
    A placed position outside the window, a chain node without data, or no
    node with data near the reach is a `BurnRefusal`."""
    m, raw = window.meta, np.asarray(window.array)
    # A NaN cell is NoData whatever the sentinel, as in the core's `is_nodata`.
    ok = valid_mask(raw, m.nodata)
    step = min(m.delta_x, m.delta_y)
    line = np.asarray(reach.line, dtype=np.float64)
    line = line[np.r_[True, np.any(np.diff(line, axis=0) != 0.0, axis=1)]]
    seg = np.diff(line, axis=0)
    lengths = np.hypot(seg[:, 0], seg[:, 1])
    cum = np.r_[0.0, np.cumsum(lengths)]
    s = np.arange(math.floor(cum[-1] / step + 1e-9) + 1) * step
    i = np.minimum(np.searchsorted(cum, s, side="right") - 1, len(seg) - 1)
    t = seg[i] / lengths[i, None]
    pts = line[i] + t * (s - cum[i])[:, None]
    rows, cols = m.index_of(pts[:, 0], pts[:, 1])
    inside = (rows >= -0.5) & (rows <= m.rows - 0.5) & (cols >= -0.5) & (cols <= m.cols - 0.5)
    at = min(int(np.floor(reach.at / step + 0.5)), len(s) - 1)
    if not inside[at]:
        raise BurnRefusal(
            f"the gauge's mapped position, {reach.at:g} m along the reach, is outside the window"
        )
    spacing = (m.x_min, m.y_max, m.delta_x, m.delta_y)
    keep = np.nonzero(inside)[0]
    chosen = [
        _floor_node(raw, ok, tuple(pts[k]), tuple(t[k]), spacing, step, reach.corridor)
        for k in keep
    ]
    joined, where = _join(chosen)
    chain, to = _taut(joined)
    # The join does not look at the nodes it steps through: one may be NoData.
    gap = next((n for n in chain if not ok[n]), None)
    if gap is not None:
        gx, gy = m.node_xy(*gap)
        raise BurnRefusal(
            f"the river line crosses a gap (NoData) in the DEM at ({gx:.0f}, {gy:.0f}): "
            "the station is refused"
        )
    placed = to[where[int(np.searchsorted(keep, at))]]
    z = [raw[n] for n in chain]
    tenth = max(1, math.ceil(len(chain) / 10))
    direction_ok = (
        sum(_steps(chain, m)) < 100.0
        or float(np.mean(z[-tenth:]) - np.mean(z[:tenth])) <= DIRECTION_M
    )
    burnt = np.array(raw, copy=True)
    drop = np.asarray(DROP_M, dtype=burnt.dtype)
    for a, b in pairwise(chain):
        burnt[b] = min(burnt[b], burnt[a] - drop)
    extended, closed = 0.0, True
    blocked = {(n[0] + dr, n[1] + dc) for n in chain[:-1] for dr, dc in [(0, 0), *_NEIGHBOURS]}
    blocked.add(chain[-1])
    while burnt[chain[-1]] < raw[chain[-1]]:
        last = chain[-1]
        free = [(last[0] + dr, last[1] + dc) for dr, dc in _NEIGHBOURS]
        free = [
            n
            for n in free
            if 0 <= n[0] < m.rows and 0 <= n[1] < m.cols and ok[n] and n not in blocked
        ]
        if not free:
            closed = False
            break
        # Ties: a straight step before a diagonal one, then the smaller (row, col).
        nxt = min(free, key=lambda n: (raw[n], abs(n[0] - last[0]) + abs(n[1] - last[1]), n))
        (more,) = _steps([last, nxt], m)
        if extended + more > END_CAP_M:
            closed = False
            break
        extended += more
        blocked |= {(last[0] + dr, last[1] + dc) for dr, dc in _NEIGHBOURS}
        blocked.add(nxt)
        chain.append(nxt)
        burnt[nxt] = min(burnt[nxt], burnt[last] - drop)
    arc = np.r_[0.0, np.cumsum(_steps(chain, m))]
    lowered = np.array([float(raw[n]) - float(burnt[n]) for n in chain])
    px, py = (np.interp(reach.at, cum, line[:, k]) for k in (0, 1))
    nx, ny = m.node_xy(*chain[placed])
    path = GaugePath(
        chain=np.array(chain, dtype=np.int64).reshape(-1, 2),
        placed=int(placed),
        arc=arc - arc[placed],
        end_extended_m=extended,
        end_closed=closed,
        direction_ok=direction_ok,
        lowered_nodes=int(np.count_nonzero(lowered > 0.0)),
        lowered_max_m=float(lowered.max(initial=0.0)),
        node_offset_m=math.hypot(nx - px, ny - py),
        dropped_m=float(np.count_nonzero(~inside)) * step,
    )
    return DemTile(meta=m, array=burnt), path


__all__ = ["DROP_M", "END_CAP_M", "BurnRefusal", "GaugePath", "burn_reach"]
