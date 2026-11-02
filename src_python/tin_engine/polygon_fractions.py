"""Exact land-cover class fractions per triangle from class polygons.

`docs/increments/105-cell-table.md`, *4. The rule*.

Each polygon is first cut into pieces of at most `most` coordinates
(:func:`subdivide`), so a triangle's cut costs the few vertices of the pieces it
meets rather than the thousands of a whole land-cover polygon. Every triangle
is then intersected with the pieces an STRtree returns (GEOS, no snapping), the
areas summed per class; the area no polygon covers is class 0. The table is
the raster route's, through its `fractions_of` and `dominant`. Pure: arrays and
shapely geometry in, no file, no CRS.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any, NamedTuple

import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry

from tin_engine.land_cover_fractions import Assigned, dominant, fractions_of

#: Halvings of one piece before it is kept as it is: a bound on the work at a
#: point no halving separates (two rings touching), not a size anyone meets.
DEPTH = 64
#: The fewest coordinates halving always reaches: a rectangle less a notch at
#: one polygon vertex is 7 vertices and the closing point.
FEWEST = 8
#: Below this share of its triangle, the uncovered rest is taken as the
#: rounding of the cut, not area (a chosen threshold).
REST = 1e-9
#: Where the two halves' areas miss the piece's by more than this share, both
#: are cut again with `shapely.intersection` (a chosen threshold).
AREA = 1e-9


def subdivide(polygon: BaseGeometry, most: int = 64) -> list[BaseGeometry]:
    """`polygon` (or each part of a multipolygon) halved along the longer side
    of its bounding box with `shapely.clip_by_rect` until each piece has at
    most `most` coordinates, holes and closing points counted. The pieces tile
    the polygon; one at or under `most` comes back as itself. `most` below 8
    is a `ValueError`: halving cannot always get under it.

    Where a hole touches the outer ring, GEOS's rectangle clip can lose area
    or join the two into one self-touching ring. Only the returned pieces are
    checked for validity, since an invalid half that is halved again costs
    nothing. The half a bad piece came from is cut again from the part it was
    cut from (from an earlier half where that part is invalid too) with the
    slower `shapely.intersection`, whose output is valid, and its pieces are
    halved and checked again; halves whose areas miss their part's are cut
    with `shapely.intersection` at once."""
    if most < FEWEST:
        raise ValueError(f"most must be at least {FEWEST}, got {most}")
    halves: list[Half] = []
    todo = []
    for root in shapely.get_parts(polygon):
        halves.append(Half(-1, root, root.bounds, 0))
        todo.append((root, len(halves) - 1))
    done: list[BaseGeometry] = []
    tried: set[int] = set()  # cut once with intersection: kept as it comes
    while todo:
        pieces, of = _halve(todo, most, halves)
        bad = set(np.asarray(of)[~shapely.is_valid(pieces)].tolist())
        for n in list(bad):
            while halves[n].up >= 0 and not halves[n].part.is_valid:
                n = halves[n].up
            bad.add(n)
        tops = [_first_bad(n, bad, halves) for n in of]
        done += [p for p, t in zip(pieces, tops, strict=True) if t < 0 or t in tried]
        redo = set(tops) - {-1} - tried
        tried |= redo
        cut = [
            (shapely.intersection(halves[t].part, shapely.box(*halves[t].rect)), t) for t in redo
        ]
        todo = [(p, t) for half, t in cut for p in _polygons(half)]
    return done


def _first_bad(n: int, bad: set[int], halves: list[Half]) -> int:
    """The earliest half in `bad` that half `n` lies in, -1 for none."""
    top = -1
    while n >= 0:
        top = n if n in bad else top
        n = halves[n].up
    return top


class Half(NamedTuple):
    """One half: the half `part` lies in (-1 for none), the part it is cut
    from, the rectangle it is cut with, and its halvings from the input."""

    up: int
    part: BaseGeometry
    rect: tuple[float, float, float, float]
    depth: int


def _polygons(geometry: BaseGeometry) -> list[Polygon]:
    return [p for p in shapely.get_parts(geometry) if isinstance(p, Polygon) and not p.is_empty]


def _halve(
    todo: list[tuple[BaseGeometry, int]], most: int, halves: list[Half]
) -> tuple[list[BaseGeometry], list[int]]:
    """The pieces of halving each `(part, half it lies in)` in `todo`, and the
    half each lies in; new halves are appended to `halves`. Where the clipped
    halves lose area, the part is cut with `shapely.intersection`."""
    pieces: list[BaseGeometry] = []
    of: list[int] = []
    while todo:
        part, n = todo.pop()
        depth = halves[n].depth
        if shapely.get_num_coordinates(part) <= most or depth == DEPTH:
            pieces.append(part)
            of.append(n)
            continue
        x0, y0, x1, y1 = part.bounds
        if x1 - x0 >= y1 - y0:
            rects = [(x0, y0, (x0 + x1) / 2, y1), ((x0 + x1) / 2, y0, x1, y1)]
        else:
            rects = [(x0, y0, x1, (y0 + y1) / 2), (x0, (y0 + y1) / 2, x1, y1)]
        cut = [shapely.clip_by_rect(part, *r) for r in rects]
        if abs(sum(h.area for h in cut) - part.area) > AREA * part.area:
            cut = [shapely.intersection(part, shapely.box(*r)) for r in rects]
        for half, rect in zip(cut, rects, strict=True):
            halves.append(Half(n, part, rect, depth + 1))
            todo += [(p, len(halves) - 1) for p in _polygons(half)]
    return pieces, of


def assign(
    vertices: npt.ArrayLike,
    triangles: npt.ArrayLike,
    polygons: Iterable[tuple[int, BaseGeometry]],
    *,
    known: frozenset[int] | Sequence[int],
    most: int = 64,
) -> Assigned:
    """The exact class fractions of every triangle (corners `vertices[:, :2]`,
    in the polygons' CRS) from `(code, polygon)` pairs. A code not in `known`,
    or a non-finite polygon coordinate, is a `ValueError`; so is a triangle
    with no area."""
    xy = np.asarray(vertices, dtype=np.float64)[:, :2]
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    pieces: list[BaseGeometry] = []
    codes: list[int] = []
    for code, polygon in polygons:
        if code not in known:
            raise ValueError(f"class {code} is not in the land-cover class file")
        if not np.isfinite(shapely.get_coordinates(polygon)).all():
            raise ValueError(f"a polygon of class {code} has a coordinate that is not finite")
        cut = subdivide(polygon, most)
        pieces += cut
        codes += [code] * len(cut)
    shapes = shapely.polygons(xy[tri]) if len(tri) else np.zeros(0, dtype=object)
    tree: npt.NDArray[Any] = np.array(pieces, dtype=object)
    t, p = shapely.STRtree(tree).query(shapes)
    areas = shapely.area(shapely.intersection(shapes[t], tree[p]))
    whole = shapely.area(shapes)
    rest = whole - np.bincount(t, weights=areas, minlength=len(tri))
    open_ = np.flatnonzero(rest > REST * whole)
    t = np.concatenate([t, open_])
    classes = np.concatenate([np.asarray(codes, dtype=np.int32)[p], np.zeros(len(open_), np.int32)])
    areas = np.concatenate([areas, rest[open_]])
    keep = areas > 0.0
    t, classes, areas = t[keep], classes[keep], areas[keep]
    pairs, inverse = np.unique(np.column_stack([t, classes]), axis=0, return_inverse=True)
    summed = np.bincount(inverse.reshape(-1), weights=areas, minlength=len(pairs))
    offsets = np.zeros(len(tri) + 1, dtype=np.uint32)
    offsets[1:] = np.cumsum(np.bincount(pairs[:, 0], minlength=len(tri)))
    kept = pairs[:, 1].astype(np.int32)
    return Assigned(fractions_of(offsets, kept, summed), dominant(offsets, kept, summed))
