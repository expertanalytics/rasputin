"""Is the catchment area well defined at this gauge? (increment 29, PR 2).

`docs/increments/29-nve-reference-catchments.md`, "Sensitivity". The area is
read along the burnt chain within `U` of the placed node, from one
accumulation: a count is trusted only when its flag bits are clear, each side
is read outward up to the first untrusted count, and the downstream side must
reach `U`. The chain must drain along itself node by node (`flow_to`), and
the area may swing by at most `SWING_MAX` within the read stretch. Arrays in,
a frozen :class:`Sensitivity` out: no DEM, no shapely, no reference polygon.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

#: The largest one-sided relative change of the area within `U`: the `match`
#: class's own 5 % bar (Ola's ruling, 2026-10-04).
SWING_MAX = 0.05
#: Metres of slack when comparing arc lengths with `U`.
_SLACK_M = 1e-6


@dataclass(frozen=True, slots=True)
class Sensitivity:
    """Areas in km2 at the placed node and the two ends read; the one-sided
    swing; the largest step between neighbouring samples and where (metres
    from the placed node, negative upstream); the metres read each way; the
    two chain checks; the causes of an uncertain station; and the verdict."""

    a0: float
    area_up: float
    area_down: float
    swing: float
    largest_step: float
    largest_step_at_m: float
    checked_up_m: float
    checked_down_m: float
    drains: bool
    monotone: bool
    causes: tuple[str, ...]
    well_posed: bool


def assess(
    count: npt.NDArray[np.uint32],
    reach_bits: npt.NDArray[np.uint8],
    flow_to: npt.NDArray[np.uint8],
    path: Any,
    cell_area_km2: float,
    uncertainty: float,
    reach_down_m: float,
) -> Sensitivity:
    """The sensitivity at `path.placed`, from `accumulate`'s three arrays
    and the burn's `path` (`chain`, `placed`, `arc`, `end_closed`)."""
    chain = np.asarray(path.chain, dtype=np.int64).reshape(-1, 2)
    arc = np.asarray(path.arc, dtype=np.float64)
    p, n, u = int(path.placed), len(chain), uncertainty
    rows, cols = chain[:, 0], chain[:, 1]
    counts = np.asarray(count)[rows, cols].astype(np.float64)
    trusted = np.asarray(reach_bits)[rows, cols] == 0
    lo = hi = p
    if trusted[p]:
        while lo > 0 and arc[lo - 1] >= -u - _SLACK_M and trusted[lo - 1]:
            lo -= 1
        while hi + 1 < n and arc[hi + 1] <= u + _SLACK_M and trusted[hi + 1]:
            hi += 1
    # Read to U: every sample in (0, U] trusted, a chain node at or past U
    # (that node is D: exactly at U it is a sample and must be trusted; past U
    # its flags do not matter), the mapped river too.
    read = bool(
        trusted[p]
        and (hi + 1 == n or arc[hi + 1] > u + _SLACK_M)
        and arc[-1] >= u - _SLACK_M
        and reach_down_m >= u - _SLACK_M
    )
    areas = counts[lo : hi + 1] * cell_area_km2
    a0 = counts[p] * cell_area_km2
    swing = max(a0 - areas[0], areas[-1] - a0) / a0 if a0 > 0 else float("inf")
    steps = np.diff(areas)
    k = int(np.argmax(steps)) if steps.size else -1
    d = chain[lo + 1 : hi + 1] - chain[lo:hi]
    toward = 3 * (d[:, 0] + 1) + (d[:, 1] + 1)
    # A flag within U upstream cannot drain through the placed node, whose
    # bits are clear (`accumulate` ORs bits downstream): not draining.
    cut = lo > 0 and arc[lo - 1] >= -u - _SLACK_M and not trusted[lo - 1]
    flows = np.asarray(flow_to)[rows[lo:hi], cols[lo:hi]] == toward
    drains = bool(not cut and np.all(flows))
    strict = bool(np.all(steps > 0.0))
    closed = bool(path.end_closed)
    causes = tuple(
        cause
        for cause, failed in (
            ("swing", not swing <= SWING_MAX),
            ("downstream_unread", not read),
            ("chain_not_draining", not (drains and strict)),
            ("chain_end_open", not closed),
        )
        if failed
    )
    return Sensitivity(
        a0=float(a0),
        area_up=float(areas[0]),
        area_down=float(areas[-1]),
        swing=float(swing),
        largest_step=float(steps[k]) if k >= 0 else 0.0,
        largest_step_at_m=float(arc[lo + 1 + k]) if k >= 0 else 0.0,
        checked_up_m=float(-arc[lo]),
        checked_down_m=float(arc[hi]),
        drains=drains,
        monotone=strict and closed,
        causes=causes,
        well_posed=not causes,
    )


__all__ = ["SWING_MAX", "Sensitivity", "assess"]
