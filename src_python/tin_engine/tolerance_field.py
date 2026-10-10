"""The lines a vertical tolerance follows, ready for ``_core.LineTolerance``.

Increment 33 (``docs/increments/33-feature-tolerance.md``, sections 3, 4.5
and 9.1). Pure: a file in, a float64 ``(k, 4)`` array of segments in the
mesh's CRS out, with the margin they were simplified by. The C++ side
subtracts that margin from every distance, so simplifying never raises the
tolerance anywhere.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import shapely
from pydantic import BaseModel
from shapely.geometry import LineString

from tin_engine.crs import reprojector, transform_bounds
from tin_engine.feature_input import FeatureError, SourceRows, read_source

#: (x_min, y_min, x_max, y_max) in the mesh's CRS.
Box = tuple[float, float, float, float]
EVERYWHERE: Box = (-np.inf, -np.inf, np.inf, np.inf)
LINES = ("LineString", "MultiLineString")


class ToleranceLines(BaseModel, frozen=True):
    """``--tolerance-near FILE N`` and ``--tolerance-ramp START END``, metres."""

    path: Path
    crs: str | None  # --tolerance-near-crs; None: the file's own
    near_m: float
    start_m: float
    end_m: float
    margin_m: float = 1.0


class ToleranceSlope(BaseModel, frozen=True):
    """``--tolerance-slope N START END`` (increment 34): N in metres, the angles
    in degrees. The CLI checks the bounds."""

    near_m: float
    start_deg: float
    end_deg: float


def line_segments(
    spec: ToleranceLines, window: Box, mesh_crs: str
) -> tuple[npt.NDArray[np.float64], float]:
    """Every line of ``spec.path`` in ``mesh_crs``, simplified by the margin
    (Douglas-Peucker), as whole segments whose own box meets ``window`` grown
    by ``end_m + margin_m``; and the margin. A line simplification empties is
    a zero-length segment at its first vertex. Polygons and points, and a
    non-finite coordinate, refused."""
    reach = spec.end_m + spec.margin_m
    x0, y0, x1, y1 = window[0] - reach, window[1] - reach, window[2] + reach, window[3] + reach

    def box_for(own: str) -> Box:
        return transform_bounds(mesh_crs, own, (x0, y0, x1, y1))

    name = spec.path.name

    def shapes_in(box: Callable[[str], Box]) -> tuple[SourceRows, list[Any], list[Any]]:
        found = read_source(spec.path, None, None, box, spec.crs)
        shapes = [g for _, g, _ in found.rows if g is not None and not g.is_empty]
        return found, shapes, [g for g in shapes if g.geom_type in LINES]

    found, shapes, lines = shapes_in(box_for)
    if not lines and found.layer is not None:
        # The R-tree found no line in reach: read every row, to tell zero
        # segments (section 5) from a file with no lines, or a mixed one.
        found, shapes, lines = shapes_in(lambda _: EVERYWHERE)
    if not lines:
        raise FeatureError(f"{name} has no lines; polygons and points are not used here")
    if len(lines) < len(shapes):
        raise FeatureError(
            f"{name} holds lines and other shapes; polygons and points are not used here"
        )
    to_mesh = reprojector(found.crs, mesh_crs)
    rows = [np.empty((0, 4))]
    for part in shapely.get_parts(lines):
        xy = to_mesh(shapely.get_coordinates(part))
        if not np.isfinite(xy).all():  # before simplify, which drops a NaN vertex
            raise FeatureError(
                f"{name}: a line has a NaN or infinite coordinate, or one with no image"
                f" in {mesh_crs}"
            )
        kept = shapely.get_coordinates(
            shapely.simplify(LineString(xy), spec.margin_m, preserve_topology=False)
        )
        kept = kept if len(kept) >= 2 else xy[[0, 0]]
        rows.append(np.hstack([kept[:-1], kept[1:]]))
    segs = np.vstack(rows)
    lo_x, hi_x = np.minimum(segs[:, 0], segs[:, 2]), np.maximum(segs[:, 0], segs[:, 2])
    lo_y, hi_y = np.minimum(segs[:, 1], segs[:, 3]), np.maximum(segs[:, 1], segs[:, 3])
    keep = (lo_x <= x1) & (hi_x >= x0) & (lo_y <= y1) & (hi_y >= y0)
    return np.ascontiguousarray(segs[keep], dtype=np.float64), spec.margin_m
