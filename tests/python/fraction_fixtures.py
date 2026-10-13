"""Hand fixtures and independent oracles for exact land-cover fractions (26b).

`docs/increments/26-land-cover.md`, sections 1.1, 1.2, 1.6, 7 (F1, F2), 8 and
10's "26b". The oracles share nothing with the producer: the clip is GEOS's
(`shapely.intersection` of a triangle with each cell's box), the
triangulation GEOS's Delaunay of a seeded point set, the cell areas
`pyproj.Geod`'s own polygon area, and the move from the mesh's CRS to the
raster pyproj's own `Transformer`.

PINNED BY THESE SUITES (the design names the parts, not their spelling):
- C++ `include/terrain/land_cover/overlap.hpp`: the header comment of
  `tests/cpp/unit/test_land_cover_overlap.cpp`.
- `tin_engine._core.class_areas(corners, classes, row_area, *, threads=0) ->
  (offsets, classes, areas)`: `corners` float64 `(T, 3, 2)` in the grid's
  index space (x the column, y the row, y growing downward; cell `(r, c)` is
  `[c, c + 1] x [r, r + 1]`), `classes` a C-contiguous 2-D `uint8` array
  viewed without a copy, `row_area` float64 `(rows,)` in m2; out `offsets`
  uint32 `(T + 1,)`, `classes` uint8 and `areas` float64 of one length. A
  triangle past the window `[0, cols] x [0, rows]` is a `ValueError` whose
  message says `outside`; a zero-area triangle is a `ValueError`; another
  dtype or a non-contiguous `classes` is a `TypeError`; a wrong shape is a
  `ValueError`. Releases the GIL; `threads` 0 is all cores.
- `tin_engine.io.models.ClassFractions(offsets, classes, fractions)`: frozen;
  `offsets` uint32 `(T + 1,)` starting at 0, non-decreasing, ending at
  `len(classes)`; `classes` uint8 and `fractions` float32 of one length; the
  arrays are kept as given (another dtype is refused, not converted);
  `triangles` is T. Anything else is a `ValueError` (pydantic's).
- `tin_engine.land_cover_fractions.row_areas(window: ClassWindow, crs: str)
  -> float64 (rows,)`: the area in m2 of one cell of each row, on the
  raster's own ellipsoid (`crs` geographic; 1.1: "computed once per row by
  pyproj.Geod"), or on its grid plane for a projected raster.
- `io.vtk_legacy.write_vtk(..., class_fractions: ClassFractions | None = None)`:
  FieldData `land_cover_offsets` (`unsigned_int`), `land_cover_classes`
  (`unsigned_char`), `land_cover_fractions` (`float`); the three names are
  reserved like `land_cover_codes`; a table whose T is not the triangle
  count is a `ValueError`.
- `io.ply.write_ply(..., face_fractions: ClassFractions | None = None)`: the
  face element's properties in order `vertex_indices`, `land_cover_classes`
  (`list uchar uchar`), `land_cover_fractions` (`list uchar float`),
  `land_cover_code` (`int`); with `edges` a `ValueError`.
- `rasputin mesh --land-cover KEY [--land-cover-year Y]` (section 6): KEY a
  catalogue source of classes read from `--cache`, the year's object as
  `fetch --year` names it, the last year by default. A year outside the
  range, `--land-cover-year` without `--land-cover`, a heights source, a run
  without `--dem`, and a coded `--features-map` beside it are usage errors
  (exit 2) naming the option; an object not in the cache is refused naming
  the key and `fetch`, with no traceback; an unknown key and a run with
  neither `--cache` nor `RASPUTIN_DATA` are usage errors too. The file's
  strings: `land_cover_codes` is `<source id> class number; 0 = no data, and
  every constraint line` and `land_cover_source` is `<source id>, <year>`
  (1.6 spells both with a display name `RemoteSource` does not have; the id
  is the name `--land-cover` and `fetch` take), `land_cover_credit` the
  source's credit rendered with the fetch date and the URL,
  `land_cover_licence_note`, `land_cover_cite`; each is also a `--stats` row
  and a `--record` entry.
- `tin_engine.land_cover_fractions.assign(vertices, triangles, *, mesh_crs,
  repository, name, chunk=..., threads=0)`: the table and the codes do not
  depend on `chunk`; the default's value is not pinned.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from itertools import pairwise
from pathlib import Path

import numpy as np
import numpy.typing as npt
import shapely
from pyproj import Geod, Transformer

from class_fixtures import CLASSES, KEY, LAT0, LON0, STEP, class_tiff, object_id
from cli_driver import geojson
from cog_fixtures import write_cache
from geographic_fixtures import projected_tiff, rough

Row = dict[int, float]

#: GEOS's clip against our pieces, in cell units (coordinates below 20 here).
CLIP = 1e-9
#: The pieces of one cell over every triangle against the cell's covered area.
CELL = 1e-12


# ------------------------------------------------------------ class grids


def mixed_grid(rows: int, cols: int) -> np.ndarray:
    """Sixteen MapBiomas classes mixed over the grid, with a block of 0 (no
    data, a class like any other: 1.1) in the lower-left corner."""
    r, c = np.indices((rows, cols))
    grid = np.asarray(CLASSES, dtype=np.uint8)[(7 * r + 3 * c + 1) % len(CLASSES)]
    grid[rows - rows // 4 :, : cols // 4] = 0
    return np.ascontiguousarray(grid)


def unique_grid(rows: int, cols: int) -> np.ndarray:
    """Every cell its own class, 1 upward: at most 255 cells, so the pieces of
    a cell can be read off the output by class."""
    assert rows * cols <= 255, "a uint8 class per cell needs at most 255 cells"
    return np.ascontiguousarray((np.arange(rows * cols, dtype=np.uint8) + 1).reshape(rows, cols))


# ------------------------------------------------------------ triangulations


def square_triangulation(seed: int, extent: float, points: int = 60) -> np.ndarray:
    """GEOS's Delaunay triangles of `points` random points in the square
    `[0.5, extent - 0.5]^2` with its corners, `(T, 3, 2)`. Every fifth point
    is snapped to a grid corner and every seventh to a grid line (section 3A's
    coincidences). The triangles tile the square: their areas sum to its."""
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0.5, extent - 0.5, (points, 2))
    lo, hi = 0.5, extent - 0.5
    pts[::5] = np.clip(np.round(pts[::5]), lo + 0.5, hi - 0.5)
    pts[::7, 1] = np.clip(np.round(pts[::7, 1]), lo + 0.5, hi - 0.5)
    corners = np.array([[lo, lo], [hi, lo], [hi, hi], [lo, hi]])
    pts = np.unique(np.vstack([corners, pts]), axis=0)
    parts = shapely.get_parts(shapely.delaunay_triangles(shapely.multipoints(pts)))
    out = np.array([shapely.get_coordinates(t)[:3] for t in parts], dtype=np.float64)
    areas = shoelace(out)
    assert (areas > 1e-9).all(), "GEOS gave a zero-area triangle"
    assert abs(areas.sum() - (hi - lo) ** 2) <= 1e-9, "the triangles do not tile the square"
    return out


def shoelace(corners: npt.ArrayLike) -> np.ndarray:
    t = np.asarray(corners, dtype=np.float64).reshape(-1, 3, 2)
    a, b, c = t[:, 0], t[:, 1], t[:, 2]
    return 0.5 * np.abs(
        (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    )


# ------------------------------------------------------------ GEOS oracles


def geos_pieces(corners: npt.ArrayLike, rows: int, cols: int) -> list[dict[tuple[int, int], float]]:
    """Per triangle, the area of its overlap with each cell `(r, c)` it meets,
    by GEOS's clip of the triangle against the cell's box. Zero pieces left out."""
    out: list[dict[tuple[int, int], float]] = []
    for tri in np.asarray(corners, dtype=np.float64).reshape(-1, 3, 2):
        polygon = shapely.polygons(tri)
        x0, y0, x1, y1 = polygon.bounds
        c_lo, c_hi = max(int(np.floor(x0)), 0), min(int(np.ceil(x1)), cols)
        r_lo, r_hi = max(int(np.floor(y0)), 0), min(int(np.ceil(y1)), rows)
        rr, cc = np.meshgrid(np.arange(r_lo, r_hi), np.arange(c_lo, c_hi), indexing="ij")
        rr, cc = rr.ravel(), cc.ravel()
        boxes = shapely.box(cc, rr, cc + 1, rr + 1)
        areas = shapely.area(shapely.intersection(polygon, boxes))
        out.append(
            {(int(r), int(c)): float(a) for r, c, a in zip(rr, cc, areas, strict=True) if a > 0}
        )
    return out


def geos_class_areas(
    corners: npt.ArrayLike, grid: np.ndarray, row_area: npt.ArrayLike
) -> list[Row]:
    """F2's oracle: per triangle, class to area in m2, the pieces of GEOS's
    clip weighted by their row's cell area."""
    weights = np.asarray(row_area, dtype=np.float64)
    rows: list[Row] = []
    for pieces in geos_pieces(corners, *grid.shape):
        row: Row = {}
        for (r, c), area in pieces.items():
            k = int(grid[r, c])
            row[k] = row.get(k, 0.0) + area * float(weights[r])
        rows.append(row)
    return rows


def cell_cover(corners: npt.ArrayLike, rows: int, cols: int) -> np.ndarray:
    """The area of each cell covered by any triangle, by GEOS."""
    cover = np.zeros((rows, cols))
    for pieces in geos_pieces(corners, rows, cols):
        for (r, c), area in pieces.items():
            cover[r, c] += area
    return cover


def square_cover(extent: float, rows: int, cols: int) -> np.ndarray:
    """The exact cover of `[0.5, extent - 0.5]^2`: 1 inside, 1/2 on an edge
    cell, 1/4 on a corner cell, 0 beyond."""
    lo, hi = 0.5, extent - 0.5
    r, c = np.indices((rows, cols), dtype=np.float64)
    w = np.clip(np.minimum(c + 1, hi) - np.maximum(c, lo), 0.0, None)
    h = np.clip(np.minimum(r + 1, hi) - np.maximum(r, lo), 0.0, None)
    return w * h


def rows_of(offsets: np.ndarray, classes: np.ndarray, values: np.ndarray) -> list[Row]:
    """The sparse table as one dict per triangle."""
    return [
        {int(k): float(v) for k, v in zip(classes[a:b], values[a:b], strict=True)}
        for a, b in pairwise(offsets)
    ]


# ------------------------------------------------------------ cell areas on the ellipsoid

GEOD = Geod(ellps="WGS84")


def geod_cell_area(lat_top: float, step: float) -> float:
    """pyproj's own area of one cell whose top edge is at `lat_top`, m2. A
    cell's area does not depend on its longitude, so the cell at lon 0 is used."""
    area, _ = GEOD.polygon_area_perimeter(
        [0.0, step, step, 0.0], [lat_top, lat_top, lat_top - step, lat_top - step]
    )
    return abs(float(area))


def geod_row_areas(lat_top: float, step: float, rows: int) -> np.ndarray:
    return np.array([geod_cell_area(lat_top - r * step, step) for r in range(rows)])


# ------------------------------------------------------------ the CLI scene

#: SIRGAS 2000 / UTM 23S: the class raster (EPSG:4326 at lon -45) lies in it.
MESH_CRS = "EPSG:31983"
#: The DEM's upper-left node, in degrees: over the class raster's lower-left
#: quarter, so the mesh reaches cells of class 0 (rows 30 on, columns 0 to 9).
LON_DEM, LAT_DEM = -44.95, -15.22
DEM_ROWS = DEM_COLS = 40
DEM_STEP = 60.0
FETCHED = "2026-10-11"
#: The years the cache fixture holds; 2025 is the entry's last year.
CACHED_YEARS = (2025, 1985)


def to_mesh_crs() -> Transformer:
    return Transformer.from_crs("EPSG:4326", MESH_CRS, always_xy=True)


def dem_origin() -> tuple[float, float]:
    """The DEM's node (0, 0) in `MESH_CRS`, on a multiple of its step."""
    x, y = to_mesh_crs().transform(LON_DEM, LAT_DEM)
    return float(np.round(x / DEM_STEP) * DEM_STEP), float(np.round(y / DEM_STEP) * DEM_STEP)


def write_dem(path: Path) -> Path:
    x0, y0 = dem_origin()
    stream = projected_tiff(rough(DEM_ROWS, DEM_COLS), x0=x0, y0=y0, step=DEM_STEP, epsg=31983)
    path.write_bytes(stream.getvalue())
    return path


def domain_ring() -> list[tuple[float, float]]:
    """A quadrilateral inside the DEM with a margin, every vertex off-node."""
    x0, y0 = dem_origin()
    far = (DEM_COLS - 1) * DEM_STEP
    return [
        (x0 + 110.3, y0 - 110.1),
        (x0 + far - 110.7, y0 - 109.3),
        (x0 + far - 111.9, y0 - far + 110.9),
        (x0 + 109.1, y0 - far + 111.7),
    ]


def write_domain(path: Path) -> Path:
    return geojson(path, domain_ring(), crs=MESH_CRS)


def write_class_cache(root: Path, years: Sequence[int] = CACHED_YEARS) -> Path:
    """The cache of the test class source for `years`, as `fetch --year` leaves
    it, each object fetched on `FETCHED`."""
    directory = write_cache(
        root, KEY, {object_id(y): class_tiff(y) for y in years}, crs="EPSG:4326"
    )
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for entry in manifest["objects"].values():
        entry["fetched"] = FETCHED
    manifest_path.write_text(json.dumps(manifest))
    return root


def oracle_fractions(
    points: npt.ArrayLike, triangles: npt.ArrayLike, grid: np.ndarray
) -> list[Row]:
    """Each triangle's class fractions from the file's own vertices: moved to
    lon/lat by pyproj's `Transformer`, into the raster's index space by the
    fixture's origin and step, clipped by GEOS, weighted by `Geod`'s cell
    area per row, normalised."""
    xy = np.asarray(points, dtype=np.float64)[:, :2]
    lon, lat = to_mesh_crs().transform(xy[:, 0], xy[:, 1], direction="INVERSE")
    index = np.column_stack([(np.asarray(lon) - LON0) / STEP, (LAT0 - np.asarray(lat)) / STEP])
    corners = index[np.asarray(triangles, dtype=np.int64).reshape(-1, 3)]
    weights = geod_row_areas(LAT0, STEP, grid.shape[0])
    out: list[Row] = []
    for row in geos_class_areas(corners, grid, weights):
        total = sum(row.values())
        out.append({k: v / total for k, v in row.items()})
    return out


def dominant(row: Row) -> int:
    """1.4 step 7: the largest fraction's class, ties to the smaller class."""
    return min(row, key=lambda k: (-row[k], k))


def same_rows(got: Iterable[Row], want: Iterable[Row], tolerance: float) -> list[str]:
    """Where two tables differ: the first twenty findings, as text."""
    findings: list[str] = []
    for t, (g, w) in enumerate(zip(got, want, strict=True)):
        if set(g) != set(w):
            findings.append(f"triangle {t}: classes {sorted(g)} against {sorted(w)}")
            continue
        for k in w:
            if abs(g[k] - w[k]) > tolerance:
                findings.append(f"triangle {t}, class {k}: {g[k]!r} against {w[k]!r}")
    return findings[:20]
