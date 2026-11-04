"""The oracle and the readers for the DEM height-share suites.

Test support only; nothing in `src_python` uses it, and it imports nothing
the height shares add, so a missing module fails the tests that reach it,
not the collection of this helper.

`docs/increments/137-dem-stats.md`, 3.1 to 3.3 and 3.6. THE ORACLE
(`oracle_table`) shares no code with the producer. In the DEM's index space
(pixel (r, c) is [c, c + 1] x [r, r + 1], node (r, c) its centre), GEOS clips
each triangle, with straight edges, to every pixel square; a piece weighs its
area times its row's pixel area. Its residual is the pixel's height minus the
triangle's plane (through the three corners with their written z) at the
pixel centre, and its class index is c with edges[c] <= r < edges[c + 1],
clamped to 0 ... 2 half, the edges (j + 1/2) w for j = -half - 1 ... half. A
pixel equal to the nodata value, or NaN, weighs on the NoData share. The
span runs from the lowest to the highest class of the non-NoData nodes the
triangle holds (each edge function over its edge's length >= -1e-9); other
pieces fold to the nearer end. A triangle holding no such node, with data,
has one class: that of its weighted mean residual over its valid pieces,
clamped. One with no valid piece has count 0, first 0, nodata 1.

PINNED BY THESE SUITES BEYOND THE DESIGN (for `@architect`, in the red
commit's handback):
- the no-node class is the mean residual weighted by piece area times row
  area (3.1 says "area-weighted");
- `.gpkg`: `dem_height_shares.triangle` is the triangle's `fid` (index + 1),
  as `land_cover_fractions.triangle` is; `dem_height_classes` has one row
  per class index 0 ... 2n with `lower_m` = edges[c], `upper_m` = edges[c+1];
- `.vtk`: `dem_height_shares`, `_offsets`, `_first` and `_edges_m` are
  one-component dataset (`FieldData`) arrays; `dem_nodata_share` is a
  one-component cell array, in any cell `FIELD` block or as `SCALARS`;
- `.ply`: the comment is `dem_height_edges_m` then the edges separated by
  single spaces, each parsing back to the same double;
- the record entries `dem_height_class_m` and
  `dem_height_triangles_without_data` are read as numbers, not as spellings.
"""

from __future__ import annotations

import hashlib
import sqlite3
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pyproj
import shapely
import tifffile

from plyread import read_ply
from vtkread import VtkFile

#: The run-record entries of 3.6.
WITHOUT_DATA = "dem_height_triangles_without_data"
CLASS_M = "dem_height_class_m"
CLASS_FROM = "dem_height_class_from"
#: `dem_height_class_from` when the DEM is a path (3.6).
DEFAULT_FROM = "default: the DEM states no vertical accuracy"
#: Every stable name of 3.6.
NAMES = (
    "dem_height_shares",
    "dem_height_offsets",
    "dem_height_first",
    "dem_nodata_share",
    "dem_height_edges_m",
)
#: Run-record rows the host fills: the command line, the core count, the version.
HOST_ROWS = ("command", "threads", "rasputin_version")
#: 3.2's held-node margin, index units.
HELD = -1e-9

F64 = npt.NDArray[np.float64]


def edges(half: int, class_m: float) -> F64:
    """3.1's edges, (j + 1/2) w for j = -half - 1 ... half, one product each."""
    return np.array([(j + 0.5) * class_m for j in range(-half - 1, half + 1)])


def half_of(tolerance_m: float, class_m: float) -> int:
    """3.1's n: the least whole n >= 0 with (n + 1/2) w > tolerance."""
    n = 0
    while (n + 0.5) * class_m <= tolerance_m:
        n += 1
    return n


@dataclass(frozen=True)
class Table:
    """3.4's sparse table, as float64 for comparison."""

    offsets: npt.NDArray[np.int64]
    first: npt.NDArray[np.int64]
    shares: F64
    nodata: F64

    @property
    def count(self) -> npt.NDArray[np.int64]:
        return np.diff(self.offsets)

    def row(self, t: int) -> F64:
        return self.shares[self.offsets[t] : self.offsets[t + 1]]

    def classes(self, t: int, half: int) -> dict[int, float]:
        """Triangle `t` as `{k: share}`, k the class in units of w (index - half),
        the zeros inside the span left out."""
        start = int(self.first[t]) - half
        return {start + i: float(v) for i, v in enumerate(self.row(t)) if v != 0.0}


def table_of(offsets: Any, first: Any, shares: Any, nodata: Any) -> Table:
    return Table(
        np.asarray(offsets, dtype=np.int64),
        np.asarray(first, dtype=np.int64),
        np.asarray(shares, dtype=np.float64),
        np.asarray(nodata, dtype=np.float64),
    )


def _class(residual: float, bounds: F64) -> int:
    return int(np.clip(np.searchsorted(bounds, residual, side="right") - 1, 0, len(bounds) - 2))


def oracle_table(
    corners: npt.ArrayLike,
    z: npt.ArrayLike,
    dem: npt.ArrayLike,
    nodata: float | None,
    row_area: npt.ArrayLike,
    *,
    class_m: float,
    half: int,
) -> Table:
    """The sparse table of `corners` `(T, 3, 2)` in index space with written
    `z` `(T, 3)` over `dem` `(rows, cols)`, `row_area` one pixel's area per row."""
    xy = np.asarray(corners, dtype=np.float64).reshape(-1, 3, 2)
    zz = np.asarray(z, dtype=np.float64).reshape(-1, 3)
    grid = np.asarray(dem, dtype=np.float64)
    weights = np.asarray(row_area, dtype=np.float64)
    bounds = edges(half, class_m)
    rows, cols = grid.shape
    offsets, first, shares, nod = [0], [], [], []
    for tri, heights in zip(xy, zz, strict=True):
        coef = np.linalg.solve(np.column_stack([np.ones(3), tri]), heights)
        r0, c0 = (max(int(np.floor(tri[:, i].min())), 0) for i in (1, 0))
        r1 = min(int(np.ceil(tri[:, 1].max())), rows)
        c1 = min(int(np.ceil(tri[:, 0].max())), cols)
        rr, cc = (a.ravel() for a in np.mgrid[r0:r1, c0:c1])
        values = grid[rr, cc]
        missing = np.isnan(values) | (False if nodata is None else values == nodata)
        residual = values - (coef[0] + coef[1] * (cc + 0.5) + coef[2] * (rr + 0.5))
        pieces = shapely.intersection(shapely.box(cc, rr, cc + 1.0, rr + 1.0), shapely.Polygon(tri))
        w = shapely.area(pieces) * weights[rr]
        held = _held(tri, cc + 0.5, rr + 0.5) & ~missing
        valid = (w > 0.0) & ~missing
        total = w.sum()
        if not valid.any():
            offsets.append(offsets[-1])
            first.append(0)
            nod.append(1.0)
            continue
        if held.any():
            k = np.array([_class(float(v), bounds) for v in residual[held]])
            lo, hi = int(k.min()), int(k.max())
        else:
            mean = float((w[valid] * residual[valid]).sum() / w[valid].sum())
            lo = hi = _class(mean, bounds)
        acc = np.zeros(hi - lo + 1)
        for v, weight in zip(residual[valid], w[valid], strict=True):
            acc[min(max(_class(float(v), bounds), lo), hi) - lo] += weight
        offsets.append(offsets[-1] + len(acc))
        first.append(lo)
        shares.extend(acc / total)
        nod.append(float(w[missing].sum() / total))
    return table_of(offsets, first, shares, nod)


def _held(tri: F64, x: F64, y: F64) -> npt.NDArray[np.bool_]:
    """3.2: each edge function over its edge's length >= -1e-9, either winding."""
    a, b, c = tri
    sign = 1.0 if (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0 else -1.0
    ok = np.ones(len(x), dtype=bool)
    for p, q in ((a, b), (b, c), (c, a)):
        e = sign * ((q[0] - p[0]) * (y - p[1]) - (q[1] - p[1]) * (x - p[0]))
        ok &= e / float(np.hypot(q[0] - p[0], q[1] - p[1])) >= HELD
    return ok


def row_sums(table: Table) -> F64:
    """Each triangle's shares plus its NoData share (invariant 2)."""
    sums = np.add.reduceat(table.shares, table.offsets[:-1]) if len(table.shares) else 0.0
    return np.where(table.count > 0, sums, 0.0) + table.nodata


def spread_and_zero(table: Table, half: int, class_m: float, area: F64) -> tuple[float, float]:
    """Test 11's figures: each triangle's standard deviation of its class
    centres k w over its valid part, and its class-0 share of the valid part,
    each weighted by plan area times the valid fraction."""
    spreads, zeros, weights = [], [], []
    for t in range(len(table.first)):
        row = table.row(t)
        if len(row) == 0:
            continue
        valid = 1.0 - table.nodata[t]
        p = row / valid
        k = (np.arange(len(row)) + table.first[t] - half) * class_m
        mean = (p * k).sum()
        spreads.append(np.sqrt((p * (k - mean) ** 2).sum()))
        zero = half - int(table.first[t])
        zeros.append(p[zero] if 0 <= zero < len(row) else 0.0)
        weights.append(area[t] * valid)
    w = np.array(weights)
    return float((w * spreads).sum() / w.sum()), float((w * zeros).sum() / w.sum())


# ---------------------------------------------------------------- Velhas

VELHAS = Path(__file__).resolve().parents[1] / "fixtures" / "velhas"


@dataclass(frozen=True)
class Source:
    """A DEM as stored, read with tifffile: `corner` the upper-left pixel's
    upper-left corner, so index = ((x - x0) / dx, (y0 - y) / dy)."""

    dem: F64
    nodata: float
    x0: float
    y0: float
    dx: float
    dy: float
    crs: str

    def index(self, vertices: npt.ArrayLike, mesh_crs: str) -> F64:
        """Vertices of `mesh_crs` in this source's index space, by pyproj."""
        xy = np.asarray(vertices, dtype=np.float64)
        move = pyproj.Transformer.from_crs(mesh_crs, self.crs, always_xy=True)
        lon, lat = move.transform(xy[:, 0], xy[:, 1])
        return np.column_stack([(lon - self.x0) / self.dx, (self.y0 - lat) / self.dy])

    def row_areas(self) -> F64:
        """One pixel's ellipsoidal area per row, on the CRS's own ellipsoid."""
        e = pyproj.CRS(self.crs).ellipsoid
        assert e is not None
        geod = pyproj.Geod(a=e.semi_major_metre, rf=e.inverse_flattening)
        out = []
        for r in range(self.dem.shape[0]):
            top, bottom = self.y0 - r * self.dy, self.y0 - (r + 1) * self.dy
            area, _ = geod.polygon_area_perimeter(
                [0.0, self.dx, self.dx, 0.0], [top, top, bottom, bottom]
            )
            out.append(abs(area))
        return np.array(out)


def velhas_source() -> Source:
    """`tests/fixtures/velhas/anadem_velhas.tif`: area-registered EPSG:4674."""
    with tifffile.TiffFile(VELHAS / "anadem_velhas.tif") as tif:
        page = tif.pages.first
        (_, _, _, x0, y0, _), (dx, dy, _) = page.tags[33922].value, page.tags[33550].value
        nodata = float(page.tags[42113].value)
        dem = page.asarray().astype(np.float64)
    return Source(dem, nodata, x0, y0, dx, dy, "EPSG:4674")


def plan_areas(points: F64, triangles: npt.ArrayLike) -> F64:
    xy = np.asarray(points)[np.asarray(triangles)][:, :, :2]
    d1, d2 = xy[:, 1] - xy[:, 0], xy[:, 2] - xy[:, 0]
    return np.asarray(0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]))


# ---------------------------------------------------------------- readers


def vtk_table(vtk: VtkFile) -> tuple[Table, F64, F64]:
    """The table of a `.vtk`, its edges, and the NoData share of its line cells."""
    want = {
        "dem_height_shares": "float",
        "dem_height_offsets": "unsigned_int",
        "dem_height_first": "int",
        "dem_height_edges_m": "double",
    }
    for name, kind in want.items():
        found = vtk.field_data[name]
        assert (found.type_name, found.components) == (kind, 1), (name, found.type_name)
    nodata = vtk.cell_array("dem_nodata_share")
    assert (nodata.type_name, nodata.components) == ("float", 1), nodata.type_name
    lines = len(vtk.lines)
    values = np.asarray(nodata.values, dtype=np.float64)
    assert len(values) == vtk.cells
    table = table_of(
        vtk.field_data["dem_height_offsets"].values,
        vtk.field_data["dem_height_first"].values,
        vtk.field_data["dem_height_shares"].values,
        values[lines:],
    )
    assert len(table.offsets) == len(vtk.polygons) + 1
    edges_m = np.asarray(vtk.field_data["dem_height_edges_m"].values, dtype=np.float64)
    return table, edges_m, values[:lines]


def ply_table(blob: bytes) -> tuple[Table, F64]:
    """The table of a `.ply`, read by the header's own declarations, and the
    comment's edges."""
    header, data = read_ply(blob)
    face = header.element("face")
    by_name = {p.name: p for p in face.properties}
    shares = by_name["dem_height_shares"]
    assert (shares.count_type, shares.type_name) == ("uchar", "float"), shares
    assert (by_name["dem_height_first"].count_type, by_name["dem_height_first"].type_name) == (
        None,
        "int",
    )
    assert (by_name["dem_nodata_share"].count_type, by_name["dem_nodata_share"].type_name) == (
        None,
        "float",
    )
    found = [c for c in header.comments if c.split(" ", 1)[0] == "dem_height_edges_m"]
    assert len(found) == 1, header.comments
    edges_m = np.array([float(w) for w in found[0].split(" ")[1:]])
    # A 2-D array when every list has one length, else one object per face.
    rows = [np.asarray(r, dtype=np.float64).ravel() for r in data["face"]["dem_height_shares"]]
    assert len(rows) == face.count
    offsets = np.concatenate([[0], np.cumsum([len(r) for r in rows])])
    flat = np.concatenate(rows) if rows else np.zeros(0)
    table = table_of(
        offsets, data["face"]["dem_height_first"], flat, data["face"]["dem_nodata_share"]
    )
    return table, edges_m


def gpkg_table(path: Path) -> tuple[Table, F64]:
    """The table rebuilt from `dem_height_shares` (rows by `fid`), the
    `triangles` column `dem_nodata_share`, and the edges of `dem_height_classes`."""
    with closing(sqlite3.connect(path)) as conn:
        triangles = conn.execute("SELECT count(*) FROM triangles").fetchone()[0]
        nodata = [r[0] for r in conn.execute("SELECT dem_nodata_share FROM triangles ORDER BY fid")]
        rows = conn.execute(
            "SELECT triangle, class, share FROM dem_height_shares ORDER BY triangle, class"
        ).fetchall()
        classes = conn.execute(
            "SELECT class, lower_m, upper_m FROM dem_height_classes ORDER BY class"
        ).fetchall()
    assert [c for c, _, _ in classes] == list(range(len(classes))), classes
    lower = np.array([lo for _, lo, _ in classes], dtype=np.float64)
    upper = np.array([up for _, _, up in classes], dtype=np.float64)
    assert (lower[1:] == upper[:-1]).all(), (lower, upper)
    counts = np.zeros(triangles, dtype=np.int64)
    first = np.zeros(triangles, dtype=np.int64)
    flat = []
    for t in range(1, triangles + 1):
        mine = [(k, s) for tri, k, s in rows if tri == t]
        if mine:
            ks = [k for k, _ in mine]
            assert ks == list(range(ks[0], ks[0] + len(ks))), (t, ks)
            first[t - 1], counts[t - 1] = ks[0], len(ks)
            flat.extend(s for _, s in mine)
    offsets = np.concatenate([[0], np.cumsum(counts)])
    return table_of(offsets, first, flat, nodata), np.concatenate([lower, upper[-1:]])


# ---------------------------------------------------------------- digests


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def gpkg_digest(path: Path) -> str:
    """SHA-256 over every table's schema and rows, in name and rowid order,
    with the host-filled run rows (`HOST_ROWS`) reduced to their names: a
    GeoPackage records the command line and the core count."""
    digest = hashlib.sha256()
    with closing(sqlite3.connect(path)) as conn:
        tables = conn.execute(
            "SELECT name, sql FROM sqlite_master WHERE type = 'table' ORDER BY name"
        ).fetchall()
        for name, sql in tables:
            digest.update(repr((name, sql)).encode())
            for row in conn.execute(f'SELECT * FROM "{name}" ORDER BY rowid'):
                kept: tuple[Any, ...] = row
                if name == "run" and row[1] in HOST_ROWS:
                    kept = (row[0], row[1])
                digest.update(repr(kept).encode())
    return digest.hexdigest()
