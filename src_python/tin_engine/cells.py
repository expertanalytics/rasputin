"""The cell table: one hydrological cell per triangle, as Shyft's TIN cells take it.

`docs/increments/105-cell-table.md`, *3. What a cell carries, and where each value comes from*.

A cell is its triangle's plan centroid, plan area, slope and aspect, and the
plan-area shares of five categories summed from the exact class fractions:
glacier, lake, forest and savanna from the legend's ``category`` keys, and
unspecified as ``1 - glacier - lake - reservoir - forest`` clamped at 0, so
savanna lies inside unspecified, as Shyft has no savanna share. Pure: arrays
and a legend in, no file, no CRS.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import numpy.typing as npt

from tin_engine.feature_input import Category, Legend
from tin_engine.io.models import CellTable, ClassFractions

Float64 = npt.NDArray[np.float64]
#: The categories `categorise` sums, in its column order; unspecified is the fifth.
SUMMED: tuple[Category, ...] = ("glacier", "lake", "forest", "savanna")


def categories_of(legend: Legend) -> Mapping[int, Category]:
    """The legend's class categories, refused with a `ValueError` naming every
    water class without one, so water never falls silently into unspecified."""
    missing = sorted(legend.water_codes - set(legend.categories))
    if missing:
        named = ", ".join(f"class {code}" for code in missing)
        raise ValueError(f"{named}: water without a category; give it lake or unspecified")
    return legend.categories


def categorise(fractions: ClassFractions, categories: Mapping[int, Category]) -> Float64:
    """(triangles, 5) float64 shares: glacier, lake, forest, savanna and
    unspecified. A class absent from ``categories`` is unspecified."""
    rows = np.repeat(np.arange(fractions.triangles), np.diff(fractions.offsets))
    shares = np.zeros((fractions.triangles, 5))
    for i, category in enumerate(SUMMED):
        codes = [code for code, named in categories.items() if named == category]
        chosen = np.isin(fractions.classes, codes)
        weights = fractions.fractions[chosen].astype(np.float64)
        shares[:, i] = np.bincount(rows[chosen], weights, minlength=fractions.triangles)
    # Shyft's 1 - glacier - lake - reservoir - forest; float32 rounding can overshoot by 1e-6.
    shares[:, 4] = np.maximum(0.0, 1.0 - shares[:, 0] - shares[:, 1] - shares[:, 2])
    return shares


def geometry(
    vertices: Float64, triangles: npt.NDArray[np.integer]
) -> tuple[Float64, Float64, Float64, Float64]:
    """Per triangle: the centroid with the mean corner height (T, 3), the plan
    area, the slope and the aspect in degrees, clockwise from grid north to
    the downslope direction; aspect NaN where the upward normal is vertical."""
    xyz = np.asarray(vertices, np.float64)[np.asarray(triangles)]
    normal = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    normal *= np.sign(normal[:, 2:3])
    nx, ny, nz = normal.T
    flat = (nx == 0) & (ny == 0)
    aspect = np.where(flat, np.nan, np.degrees(np.arctan2(nx, ny)) % 360)
    slope = np.degrees(np.arctan2(np.hypot(nx, ny), nz))
    return xyz.mean(1), 0.5 * nz, slope, aspect


def cell_table(
    vertices: Float64,
    triangles: npt.NDArray[np.integer],
    fractions: ClassFractions,
    categories: Mapping[int, Category],
) -> CellTable:
    """The cells of a mesh and its class fractions; reservoir 0 on every row."""
    midpoints, area, slope, aspect = geometry(vertices, triangles)
    glacier, lake, forest, savanna, unspecified = categorise(fractions, categories).T
    shares = np.column_stack([glacier, lake, np.zeros_like(lake), forest, unspecified, savanna])
    return CellTable(
        midpoints=midpoints, area_m2=area, slope_deg=slope, aspect_deg=aspect, shares=shares
    )
