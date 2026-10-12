"""26b acceptance: exactness and cost of the land-cover fractions on one Velhas mesh.

Run with the worktree's venv, after `run_meshes.sh`, one mesh at a time:

    .venv/bin/python docs/benchmarks/2026-10-12/26b-fractions/analyse.py \
        --mesh .../velhas_t10_lcon.bin.vtk --tag t10 [--threads 1,2,4,8] [--repeats 3]

What it does, against the mesh file as written by `rasputin mesh --land-cover`:

1. Reads the binary legacy VTK back (points, triangles, the three fraction
   arrays, `land_cover_code`, the `crs` field).
2. Builds the independent reference from the cached MapBiomas raster: the
   area of the cells whose centres fall inside the outline (`ref`), split into
   cells the outline's boundary does not touch (`interior`) and cells it
   crosses (`boundary`), so the exact class area lies in
   `[interior, interior + boundary]`. Cell areas are WGS 84 ellipsoidal, one
   per row (`pyproj.Geod`), as the producer's. The boundary cells come from
   GEOS (`shapely.intersects` of each outline segment with the cells of its box).
3. Recomputes the fractions with `land_cover_fractions.assign` (the function
   the CLI calls), checks the file's `land_cover_fractions` and
   `land_cover_code` against it bit for bit, and sums the exact per-triangle
   class areas per class. A second, independent sum uses the file's float32
   fractions times each triangle's own ellipsoidal area (`Geod` on its
   corners in longitude and latitude).
4. Times `assign` by thread count (median of `--repeats`), and its parts.
5. Per-cell check on a 1 km square around the mesh's centroid: for every cell
   of the square, `_core.class_areas` with that cell alone as class 1 against
   GEOS's clip of the same triangles with the cell's box, and against the
   cell's own area (the square is interior, so every cell is fully covered).
6. The misplaced share on 1, 5 and 25 km squares (UTM-aligned): mesh class
   areas binned by triangle centroid against raster class areas binned by cell
   centre. Both binnings are approximate at the squares' edges, so this is an
   upper bound on what the fractions misplace, not a measurement of it.

Writes `raw/analysis_<tag>.json` beside the mesh runs and prints the tables.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import shapely
from pyproj import Geod, Transformer

from tin_engine import _core
from tin_engine.io.models import ClassWindow, IndexWindow
from tin_engine.io.repository import CacheRepository
from tin_engine.land_cover_fractions import (
    CHUNK,
    assign,
    fractions_of,
    index_space,
    row_areas,
    window_of,
)

HERE = Path(__file__).resolve().parent
DATA = Path("/Users/skavhaug/projects/rasputin_data")
OUTLINE = DATA / "sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson"
CACHE = DATA / "cache"
SOURCE, NAME = "mapbiomas-c11", "brazil_coverage-col11_2025"
GEOD = Geod(ellps="WGS84")
TYPES = {"double": ">f8", "float": ">f4", "int": ">i4", "unsigned_int": ">u4", "unsigned_char": "u1"}


# ------------------------------------------------------------- the VTK reader
def read_vtk_binary(path: Path) -> dict:
    data = path.read_bytes()
    pos = 0

    def line() -> str:
        nonlocal pos
        end = data.index(b"\n", pos)
        text = data[pos:end].decode("ascii")
        pos = end + 1
        return text

    def block(dtype: str, count: int) -> np.ndarray:
        nonlocal pos
        arr = np.frombuffer(data, dtype=dtype, count=count, offset=pos)
        pos += arr.nbytes
        assert data[pos : pos + 1] == b"\n", (pos, data[pos : pos + 8])
        pos += 1
        return arr

    def strings(count: int) -> list[str]:
        nonlocal pos
        out = []
        for _ in range(count):
            tag = data[pos] >> 6
            width = {0b11: 1, 0b10: 2, 0b01: 4}.get(tag, 8)
            raw = int.from_bytes(data[pos : pos + width], "big")
            length = raw & ((1 << (8 * width - 2)) - 1) if width < 8 else raw
            pos += width
            out.append(data[pos : pos + length].decode("ascii"))
            pos += length
        assert data[pos : pos + 1] == b"\n"
        pos += 1
        return out

    def field(count: int) -> dict:
        out = {}
        for _ in range(count):
            name, _, n, kind = line().split()
            out[name] = strings(int(n)) if kind == "string" else block(TYPES[kind], int(n))
        return out

    assert line().startswith("# vtk DataFile")
    line()
    assert line() == "BINARY"
    assert line() == "DATASET POLYDATA"
    head = line().split()
    assert head[:2] == ["FIELD", "FieldData"]
    fields = field(int(head[2]))
    head = line().split()
    assert head[0] == "POINTS"
    points = block(">f8", 3 * int(head[1])).reshape(-1, 3).astype(np.float64)
    head = line().split()
    lines = np.zeros((0, 2), dtype=np.int64)
    if head[0] == "LINES":
        lines = block(">i4", int(head[2])).reshape(-1, 3)[:, 1:].astype(np.int64)
        head = line().split()
    assert head[0] == "POLYGONS", head
    polygons = block(">i4", int(head[2])).reshape(-1, 4)
    assert (polygons[:, 0] == 3).all()
    triangles = polygons[:, 1:].astype(np.int64)
    assert line().startswith("POINT_DATA")
    assert line().startswith("SCALARS elevation")
    line()
    block(">f8", len(points))
    assert line().startswith("CELL_DATA")
    assert line().startswith("SCALARS feature_mask")
    line()
    masks = block(">u4", len(lines) + len(triangles))
    cells = {}
    if pos < len(data):
        head = line().split()
        assert head[:2] == ["FIELD", "features"]
        cells = field(int(head[2]))
    assert pos == len(data), (pos, len(data))
    return {"fields": fields, "points": points, "lines": lines, "triangles": triangles,
            "masks": masks, "cells": cells}  # fmt: skip


def mesh_sha256(path: Path) -> str:
    """As tools/bench_quick.py hashes: from POINTS on, so header fields do not count."""
    import hashlib

    data = path.read_bytes()
    return hashlib.sha256(data[data.find(b"\nPOINTS ") + 1 :]).hexdigest()


# -------------------------------------------------------------- the reference
def outline_lonlat() -> np.ndarray:
    ring = json.loads(OUTLINE.read_text())["features"][0]["geometry"]["coordinates"]
    assert len(ring) == 1, "one ring, no holes"
    # EPSG:4674 (SIRGAS 2000) to the raster's EPSG:4326: pyproj's own move.
    to_raster = Transformer.from_crs("EPSG:4674", "EPSG:4326", always_xy=True)
    lon, lat = to_raster.transform(*np.asarray(ring[0], dtype=np.float64).T)
    return np.column_stack([lon, lat])


def reference(repo: CacheRepository, meta, ring_lonlat: np.ndarray) -> dict:
    """Class areas of the cells whose centres are inside, split into interior
    and boundary cells; the class window covering the outline's box."""
    corner_x, corner_y = meta.x_min - meta.delta_x / 2, meta.y_max + meta.delta_y / 2
    ring_ix = np.column_stack(
        [(ring_lonlat[:, 0] - corner_x) / meta.delta_x, (corner_y - ring_lonlat[:, 1]) / meta.delta_y]
    )
    c0, r0 = int(np.floor(ring_ix[:, 0].min())) - 1, int(np.floor(ring_ix[:, 1].min())) - 1
    c1, r1 = int(np.ceil(ring_ix[:, 0].max())) + 2, int(np.ceil(ring_ix[:, 1].max())) + 2
    window = repo.load_class_window(NAME, IndexWindow(row0=r0, col0=c0, rows=r1 - r0, cols=c1 - c0))
    classes = window.array
    rows, cols = classes.shape
    cell_area = row_areas(window)  # (rows,), m2, WGS 84 ellipsoid
    # Cell centres inside the outline, in the window's local index space.
    local = ring_ix - np.array([c0, r0], dtype=np.float64)
    poly = shapely.Polygon(local)
    assert poly.is_valid, "the outline must be a valid polygon"
    shapely.prepare(poly)  # in place; shapely 2.2 returns a numpy bool, not the geometry
    cx = np.arange(cols, dtype=np.float64) + 0.5
    cy = np.arange(rows, dtype=np.float64) + 0.5
    inside = np.zeros((rows, cols), dtype=bool)
    for r in range(rows):
        inside[r] = shapely.contains_xy(poly, cx, np.full(cols, cy[r]))
    # Boundary cells: every cell whose closed box the outline's boundary meets.
    boundary = np.zeros((rows, cols), dtype=bool)
    for p, q in zip(local[:-1], local[1:], strict=True):
        seg = shapely.LineString([p, q])
        ca, cb = int(np.floor(min(p[0], q[0]))), int(np.ceil(max(p[0], q[0])))
        ra, rb = int(np.floor(min(p[1], q[1]))), int(np.ceil(max(p[1], q[1])))
        cc, rr = np.meshgrid(np.arange(ca, cb + 1), np.arange(ra, rb + 1))
        cc, rr = cc.ravel(), rr.ravel()
        boxes = shapely.box(cc, rr, cc + 1, rr + 1)
        hit = shapely.intersects(boxes, seg)
        boundary[rr[hit], cc[hit]] = True
    interior = inside & ~boundary
    area = np.broadcast_to(cell_area[:, None], classes.shape)

    def per_class(mask: np.ndarray) -> np.ndarray:
        return np.bincount(classes[mask].astype(np.int64), weights=area[mask], minlength=256)

    return {
        "window": window, "c0": c0, "r0": r0, "cell_area": cell_area, "inside": inside,
        "boundary": boundary, "ref": per_class(inside), "interior": per_class(interior),
        "boundary_area": per_class(boundary), "cells_inside": int(inside.sum()),
        "cells_boundary": int(boundary.sum()),
        "outline_area_m2": abs(GEOD.polygon_area_perimeter(ring_lonlat[:, 0], ring_lonlat[:, 1])[0]),
    }  # fmt: skip


# ----------------------------------------------------------------- the mesh
def triangle_geod_areas(points: np.ndarray, triangles: np.ndarray, mesh_crs: str) -> np.ndarray:
    """Each triangle's WGS 84 ellipsoidal area from its corners in longitude and
    latitude: pyproj's own area, independent of the index-space construction."""
    to_lonlat = Transformer.from_crs(mesh_crs, "EPSG:4326", always_xy=True)
    lon, lat = to_lonlat.transform(points[:, 0], points[:, 1])
    out = np.empty(len(triangles), dtype=np.float64)
    for t, (a, b, c) in enumerate(triangles):
        out[t] = abs(GEOD.polygon_area_perimeter([lon[a], lon[b], lon[c]], [lat[a], lat[b], lat[c]])[0])
    return out


def timed_assign(points, triangles, mesh_crs, repo, threads: int, repeats: int):
    times = []
    result = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = assign(points, triangles, mesh_crs=mesh_crs, repository=repo, name=NAME, threads=threads)
        times.append(time.perf_counter() - t0)
    return result, times


def parts_timing(points, triangles, mesh_crs, repo, meta, threads: list[int], repeats: int) -> dict:
    """assign's three parts, timed apart: the move to index space, the window
    reads, and _core.class_areas by thread count (median of repeats)."""
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    t0 = time.perf_counter()
    corners = index_space(points, mesh_crs, meta)[tri]
    t_index = time.perf_counter() - t0
    windows, locals_, loads, shapes = [], [], [], []
    for start in range(0, len(tri), CHUNK):
        piece = corners[start : start + CHUNK]
        window = window_of(piece, meta)
        t0 = time.perf_counter()
        classes = repo.load_class_window(NAME, window)
        loads.append(time.perf_counter() - t0)
        windows.append(classes)
        shapes.append(list(classes.array.shape))
        locals_.append(piece - np.array([window.col0, window.row0], dtype=np.float64))
    areas = [row_areas(w) for w in windows]
    by_threads = {}
    for n in threads:
        samples = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            for local, w, ra in zip(locals_, windows, areas, strict=True):
                _core.class_areas(local, w.array, ra, threads=n)
            samples.append(time.perf_counter() - t0)
        by_threads[n] = samples
    return {"index_space_s": t_index, "window_loads_s": loads, "window_shapes": shapes,
            "class_areas_s": by_threads}  # fmt: skip


# ------------------------------------------------------- the per-cell check
def per_cell_check(points, triangles, mesh_crs, meta, ref: dict, side_m: float = 1000.0) -> dict:
    """Every cell of a `side_m` square around the mesh's centroid: class_areas
    with that cell alone as class 1, against GEOS's clip and the cell's area."""
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    cx, cy = points[:, 0].mean(), points[:, 1].mean()
    square = np.array([[cx - side_m / 2, cy - side_m / 2], [cx + side_m / 2, cy - side_m / 2],
                       [cx + side_m / 2, cy + side_m / 2], [cx - side_m / 2, cy + side_m / 2]])  # fmt: skip
    sq_ix = index_space(square, mesh_crs, meta)
    ca, cb = int(np.floor(sq_ix[:, 0].min())), int(np.ceil(sq_ix[:, 0].max()))
    ra, rb = int(np.floor(sq_ix[:, 1].min())), int(np.ceil(sq_ix[:, 1].max()))
    corners = index_space(points, mesh_crs, meta)[tri]
    lo, hi = corners.min(axis=1), corners.max(axis=1)
    keep = (hi[:, 0] >= ca) & (lo[:, 0] <= cb) & (hi[:, 1] >= ra) & (lo[:, 1] <= rb)
    chosen = corners[keep]
    window = window_of(chosen, meta)
    local = chosen - np.array([window.col0, window.row0], dtype=np.float64)
    # Row areas of this window, from the raster's geometry (the same Geod call).
    corner_x, corner_y = meta.x_min - meta.delta_x / 2, meta.y_max + meta.delta_y / 2
    probe = ClassWindow(
        x_min=corner_x + window.col0 * meta.delta_x, y_max=corner_y - window.row0 * meta.delta_y,
        delta_x=meta.delta_x, delta_y=meta.delta_y,
        array=np.zeros((window.rows, window.cols), dtype=np.uint8),
    )  # fmt: skip
    ra_rows = row_areas(probe)
    polys = shapely.polygons(local)
    worst_vs_geos = worst_vs_cell = 0.0
    n = 0
    for r in range(ra, rb):
        for c in range(ca, cb):
            n += 1
            lr, lc = r - window.row0, c - window.col0
            classes = np.zeros((window.rows, window.cols), dtype=np.uint8)
            classes[lr, lc] = 1
            offsets, cls, areas = _core.class_areas(local, classes, ra_rows, threads=1)
            got = float(areas[cls == 1].sum())
            geos = float(shapely.area(shapely.intersection(polys, shapely.box(lc, lr, lc + 1, lr + 1))).sum()) * ra_rows[lr]
            cell = float(ra_rows[lr])
            worst_vs_geos = max(worst_vs_geos, abs(got - geos) / cell)
            worst_vs_cell = max(worst_vs_cell, abs(got - cell) / cell)
    return {"square_side_m": side_m, "cells": n, "triangles_meeting": int(keep.sum()),
            "worst_rel_vs_geos": worst_vs_geos, "worst_rel_vs_cell_area": worst_vs_cell,
            "square_rows": [ra, rb], "square_cols": [ca, cb]}  # fmt: skip


# ----------------------------------------------------------- misplaced share
def misplaced(points, triangles, mesh_crs, areas_table, ref: dict, meta, sides_m=(1000.0, 5000.0, 25000.0)) -> dict:
    offsets, classes, areas = areas_table
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    centroid = points[tri, :2].mean(axis=1)
    counts = np.diff(offsets.astype(np.int64))
    tri_of_entry = np.repeat(np.arange(len(tri)), counts)
    window, inside, c0, r0 = ref["window"], ref["inside"], ref["c0"], ref["r0"]
    rr, cc = np.nonzero(inside)
    corner_x, corner_y = meta.x_min - meta.delta_x / 2, meta.y_max + meta.delta_y / 2
    lon = corner_x + (c0 + cc + 0.5) * meta.delta_x
    lat = corner_y - (r0 + rr + 0.5) * meta.delta_y
    to_mesh = Transformer.from_crs("EPSG:4326", mesh_crs, always_xy=True)
    x, y = to_mesh.transform(lon, lat)
    cell_class = window.array[rr, cc].astype(np.int64)
    cell_area = ref["cell_area"][rr]
    out = {}
    for side in sides_m:
        kx, ky = np.floor(x / side).astype(np.int64), np.floor(y / side).astype(np.int64)
        tx, ty = np.floor(centroid[:, 0] / side).astype(np.int64), np.floor(centroid[:, 1] / side).astype(np.int64)
        keys = np.unique(np.concatenate([kx * 100000 + ky, tx * 100000 + ty]))
        ki = np.searchsorted(keys, kx * 100000 + ky)
        ti = np.searchsorted(keys, tx * 100000 + ty)
        raster = np.zeros((len(keys), 256))
        np.add.at(raster, (ki, cell_class), cell_area)
        mesh = np.zeros((len(keys), 256))
        np.add.at(mesh, (ti[tri_of_entry], classes.astype(np.int64)), areas)
        full = raster.sum(axis=1) >= 0.5 * side * side
        share = np.abs(mesh - raster).sum(axis=1)[full] / (mesh.sum(axis=1) + raster.sum(axis=1))[full]
        out[f"{side / 1000:g} km"] = {
            "squares": int(full.sum()), "median": float(np.median(share)),
            "p95": float(np.percentile(share, 95)), "max": float(share.max()),
        }  # fmt: skip
    return out


# ---------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", type=Path, required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--threads", default="1,2,4,8")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--skip-misplaced", action="store_true")
    args = ap.parse_args()
    threads = [int(t) for t in args.threads.split(",")]
    out: dict = {"mesh": str(args.mesh), "tag": args.tag, "mesh_sha256_from_points": mesh_sha256(args.mesh)}

    mesh = read_vtk_binary(args.mesh)
    points, triangles = mesh["points"], mesh["triangles"]
    mesh_crs = mesh["fields"]["crs"][0]
    f_off = mesh["fields"]["land_cover_offsets"].astype(np.uint32)
    f_cls = mesh["fields"]["land_cover_classes"].astype(np.uint8)
    f_frac = mesh["fields"]["land_cover_fractions"].astype(np.float32)
    f_code = mesh["cells"]["land_cover_code"].astype(np.int32)[len(mesh["lines"]) :]
    counts = np.diff(f_off.astype(np.int64))
    out.update(
        crs=mesh_crs, vertices=len(points), triangles=len(triangles), entries=len(f_cls),
        classes_in_file=sorted(int(c) for c in np.unique(f_cls)),
        classes_per_triangle={"mean": float(counts.mean()), "max": int(counts.max()),
                              "histogram": {int(k): int(v) for k, v in zip(*np.unique(counts, return_counts=True))}},
        land_cover_source=mesh["fields"].get("land_cover_source"),
    )  # fmt: skip
    print(f"{args.tag}: {len(triangles)} triangles, {len(f_cls)} entries, {len(out['classes_in_file'])} classes, "
          f"classes per triangle mean {counts.mean():.3f} max {counts.max()}")  # fmt: skip

    repo = CacheRepository(CACHE, SOURCE, nodata=0.0)
    meta = {f.name: f.meta for f in repo.footprints()}[NAME]
    out["raster"] = {"cols": meta.cols, "rows": meta.rows, "delta_deg": meta.delta_x, "crs": meta.crs}

    # 3. assign, bit for bit against the file, and the thread sweep (4).
    sweep = {}
    assigned = None
    for n in threads:
        assigned, times = timed_assign(points, triangles, mesh_crs, repo, n, args.repeats)
        sweep[n] = times
        print(f"assign threads={n}: {[round(t, 3) for t in times]} s")
    assert assigned is not None
    table = assigned.table
    same = (
        np.array_equal(table.offsets, f_off) and np.array_equal(table.classes, f_cls)
        and np.array_equal(table.fractions, f_frac) and np.array_equal(assigned.codes, f_code)
    )
    out["file_equals_assign_bit_for_bit"] = bool(same)
    print("file arrays equal assign's output bit for bit:", same)
    # Exact per-triangle class areas: redo the chunks as assign does, keeping the areas.
    tri = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    corners = index_space(points, mesh_crs, meta)[tri]
    parts = []
    for start in range(0, len(tri), CHUNK):
        piece = corners[start : start + CHUNK]
        window = window_of(piece, meta)
        classes = repo.load_class_window(NAME, window)
        local = piece - np.array([window.col0, window.row0], dtype=np.float64)
        parts.append(_core.class_areas(local, classes.array, row_areas(classes), threads=0))
    offs, at = [np.zeros(1, dtype=np.int64)], 0
    for p in parts:
        offs.append(p[0][1:].astype(np.int64) + at)
        at += int(p[0][-1])
    offsets = np.concatenate(offs)
    cls_all = np.concatenate([p[1] for p in parts])
    areas_all = np.concatenate([p[2] for p in parts])
    assert np.array_equal(fractions_of(offsets.astype(np.uint32), cls_all, areas_all).fractions, f_frac)
    mesh_class = np.bincount(cls_all.astype(np.int64), weights=areas_all, minlength=256)
    tri_total = np.add.reduceat(areas_all, offsets[:-1])
    out["parts"] = parts_timing(points, triangles, mesh_crs, repo, meta, threads, args.repeats)
    print("parts:", {k: v for k, v in out["parts"].items() if k != "class_areas_s"},
          {n: round(statistics.median(v), 3) for n, v in out["parts"]["class_areas_s"].items()})

    # Independent second sum: the file's float32 fractions times Geod triangle areas.
    t0 = time.perf_counter()
    geod_tri = triangle_geod_areas(points, triangles, mesh_crs)
    out["geod_triangle_areas_s"] = time.perf_counter() - t0
    tri_of_entry = np.repeat(np.arange(len(tri)), counts)
    file_class = np.bincount(f_cls.astype(np.int64), weights=f_frac.astype(np.float64) * geod_tri[tri_of_entry], minlength=256)
    out["triangle_area_vs_geod"] = {
        "max_rel": float(np.max(np.abs(tri_total - geod_tri) / geod_tri)),
        "total_exact_m2": float(tri_total.sum()), "total_geod_m2": float(geod_tri.sum()),
    }

    # 2. The reference and the exactness table.
    t0 = time.perf_counter()
    ref = reference(repo, meta, outline_lonlat())
    out["reference_s"] = time.perf_counter() - t0
    out["outline_area_m2"] = ref["outline_area_m2"]
    out["cells_inside"] = ref["cells_inside"]
    out["cells_boundary"] = ref["cells_boundary"]
    rows = []
    all_pass = True
    for k in range(256):
        if ref["ref"][k] == 0 and mesh_class[k] == 0 and ref["boundary_area"][k] == 0:
            continue
        lo, hi = ref["interior"][k], ref["interior"][k] + ref["boundary_area"][k]
        tol = max(1e-3 * ref["ref"][k], 1e4)
        excess = max(lo - mesh_class[k], mesh_class[k] - hi, 0.0)
        ok = excess <= tol
        all_pass &= bool(ok)
        rows.append({"class": k, "mesh_m2": float(mesh_class[k]), "file_m2": float(file_class[k]),
                     "ref_m2": float(ref["ref"][k]), "interior_m2": float(lo), "touching_m2": float(hi),
                     "tol_m2": float(tol), "outside_band_m2": float(excess), "pass": bool(ok)})  # fmt: skip
    rows.sort(key=lambda r: -r["ref_m2"])
    out["exactness"] = rows
    out["exactness_all_pass"] = all_pass
    print(f"outline (Geod) {ref['outline_area_m2'] / 1e6:.3f} km2; mesh exact {tri_total.sum() / 1e6:.3f} km2; "
          f"cells inside {ref['cells_inside']}, boundary {ref['cells_boundary']}")  # fmt: skip
    print("| class | mesh km2 | ref (centres) km2 | band [interior, touching] km2 | mesh - ref ha | outside band ha | tol ha | pass |")
    for r in rows:
        print(f"| {r['class']} | {r['mesh_m2'] / 1e6:.4f} | {r['ref_m2'] / 1e6:.4f} | "
              f"[{r['interior_m2'] / 1e6:.4f}, {r['touching_m2'] / 1e6:.4f}] | {(r['mesh_m2'] - r['ref_m2']) / 1e4:+.2f} | "
              f"{r['outside_band_m2'] / 1e4:.3f} | {r['tol_m2'] / 1e4:.1f} | {'yes' if r['pass'] else 'NO'} |")  # fmt: skip
    print("exactness all pass:", all_pass)

    # 5. The per-cell check against GEOS on a 1 km square.
    out["per_cell"] = per_cell_check(points, triangles, mesh_crs, meta, ref)
    print("per cell:", out["per_cell"])

    # 6. Misplaced share.
    if not args.skip_misplaced:
        out["misplaced"] = misplaced(points, triangles, mesh_crs, (offsets, cls_all, areas_all), ref, meta)
        print("misplaced:", out["misplaced"])

    out["assign_s"] = sweep
    out["assign_median_s"] = {n: statistics.median(v) for n, v in sweep.items()}
    out["assign_s_per_million_triangles"] = {
        n: statistics.median(v) / (len(tri) / 1e6) for n, v in sweep.items()
    }
    print("assign median s by threads:", {n: round(v, 3) for n, v in out["assign_median_s"].items()})
    print("assign s per million triangles:", {n: round(v, 3) for n, v in out["assign_s_per_million_triangles"].items()})

    def plain(o):
        if isinstance(o, dict):
            return {str(k): plain(v) for k, v in o.items() if k not in ("window", "inside", "boundary", "cell_area")}
        if isinstance(o, (list, tuple)):
            return [plain(v) for v in o]
        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        return o

    (HERE / "raw" / f"analysis_{args.tag}.json").write_text(json.dumps(plain(out), indent=1) + "\n")


if __name__ == "__main__":
    sys.exit(main())
