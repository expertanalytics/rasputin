"""The independent edge-strip check for 15f-3's acceptance (@perf, 2026-10-04).

It does not import tin_engine. It reads the inputs itself: the domain GeoJSON
with shapely, CORINE from the GeoPackage with sqlite3 and shapely, DTM10 tiles
with tifffile and their .tfw world files, and the mesh by parsing rasputin's
ASCII VTK.

The points checked are where the input outline and the CORINE borders (clipped
to the domain) cross the DEM's grid lines. The grid's nodes are at pixel
centres, as the .tfw files give them. Each crossing segment's midpoint is also
checked, between neighbouring crossings, with the segment's ends counting as
neighbours. Each point gets z bilinear from the DEM in its cell, skipped when a
corner is NoData. That is compared with the written mesh's z, linear along the
nearest constraint line (VTK LINES) within --snap metres. A point with no
constraint line that near is counted as "off the lines". The error is
|mesh - DEM|, and over means error > tolerance + 1e-9 m.

It also counts the mesh's constraint-line vertices that sit at a candidate
crossing or midpoint (within --match m). Compared with the base mesh, this
answers Q1: how many midpoints the strip inserts at all.

Controls:
- alignment: mesh vertices that lie on a DEM node (lattice coordinates within
  1e-6 of integers) must carry the node's value. This checks that the grid
  convention used here is the mesher's.
- --plant-dx: shifts the DEM east by that many metres before sampling, so the
  check must then find points over.

Usage: python strip_check.py MESH.vtk --domain D.geojson --dem-dir DIR
       --tolerance T [--gpkg G --layer corine2018] [--plant-dx 10] [--json OUT]
   or: python strip_check.py MESH.vtk --resampled WINDOW.tif --out-crs CRS
       --spacing 30 --lines-from-mesh --tolerance T [--plant-dx 30] [--json OUT]
The reprojected form rebuilds the target grid (class Resampled) and checks the
mesh's own constraint lines (the outline as meshed, in the output CRS).
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
import time
from pathlib import Path

import numpy as np
import shapely
import tifffile

EPS = 1e-9


def read_vtk(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Points (n, 3) and LINES (m, 2) of a legacy ASCII POLYDATA file."""
    text = path.read_bytes().decode("latin-1").split("\n")
    if len(text) > 2 and text[2].strip() == "BINARY":
        sys.exit(f"{path}: BINARY; this check reads --ascii meshes")
    found: dict[str, list[str]] = {}
    for i, line in enumerate(text):
        head = line.split()
        if head and head[0] in ("POINTS", "LINES") and head[0] not in found:
            found[head[0]] = text[i + 1 : i + 1 + int(head[1])]
    pts = np.array(" ".join(found["POINTS"]).split(), dtype=np.float64).reshape(-1, 3)
    lines = np.array(" ".join(found["LINES"]).split(), dtype=np.int64).reshape(-1, 3)[:, 1:]
    return pts, lines


def gpkg_geometry(blob: bytes) -> shapely.Geometry:
    flags = blob[3]
    envelope = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[(flags >> 1) & 0b111]
    return shapely.from_wkb(blob[8 + envelope :])


def input_lines(domain: shapely.Geometry, gpkg: Path | None, layer: str) -> shapely.Geometry:
    """The outline and the clipped CORINE borders, noded and without duplicates."""
    parts = [domain.boundary]
    if gpkg is not None:
        x0, y0, x1, y1 = domain.bounds
        con = sqlite3.connect(f"file:{gpkg}?mode=ro", uri=True)
        rows = con.execute(
            f"select c.geom from {layer} c join rtree_{layer}_geom r on c.fid = r.id "
            "where r.maxx >= ? and r.minx <= ? and r.maxy >= ? and r.miny <= ?",
            (x0, x1, y0, y1),
        ).fetchall()
        for (g,) in rows:
            p = gpkg_geometry(g)
            if p.intersects(domain):
                parts.append(shapely.intersection(p, domain).boundary)
    return shapely.union_all(parts)


def segments(lines: shapely.Geometry) -> np.ndarray:
    out = []
    for part in shapely.get_parts(lines):
        if part.geom_type in ("LineString", "LinearRing"):
            c = np.asarray(part.coords)[:, :2]
            out.append(np.hstack([c[:-1], c[1:]]))
    s = np.vstack(out)
    return s[np.hypot(s[:, 2] - s[:, 0], s[:, 3] - s[:, 1]) > 0]


class Dem:
    """The DTM10 tiles over a box, as one array; nodes at pixel centres."""

    def __init__(self, folder: Path, box: tuple[float, float, float, float], dx_plant: float = 0.0):
        x0, y0, x1, y1 = box
        tiles = []
        for tfw in sorted(folder.glob("*.tfw")):
            a, _, _, e, cx, cy = (float(v) for v in tfw.read_text().split())
            tif = tfw.with_suffix(".tif")
            with tifffile.TiffFile(tif) as t:
                h, w = t.pages[0].shape
            if cx + a * (w - 1) < x0 or cx > x1 or cy + e * (h - 1) > y1 or cy < y0:
                continue
            tiles.append((tif, a, e, cx, cy, h, w))
        self.d = tiles[0][1]
        assert all(abs(t[1] - self.d) < 1e-9 and abs(t[2] + self.d) < 1e-9 for t in tiles)
        d = self.d
        # One lattice for all tiles: node (0, 0) at the box's north-west node.
        ox = min(t[3] for t in tiles)
        oy = max(t[4] for t in tiles)
        c_hi = max(round((t[3] - ox) / d) + t[6] for t in tiles)
        r_hi = max(round((oy - t[4]) / d) + t[5] for t in tiles)
        z = np.full((r_hi, c_hi), np.nan)
        for tif, _, _, cx, cy, h, w in tiles:
            with tifffile.TiffFile(tif) as t:
                arr = t.pages[0].asarray().astype(np.float64)
                nd = t.pages[0].tags.get("GDAL_NODATA")
            if nd is not None:
                arr[arr == float(str(nd.value).strip("\x00"))] = np.nan
            c, r = round((cx - ox) / d), round((oy - cy) / d)
            assert abs(cx - ox - c * d) < 1e-6 and abs(oy - cy - r * d) < 1e-6, "tiles off one lattice"
            sub = z[r : r + h, c : c + w]
            keep = np.isnan(sub)
            sub[keep] = arr[keep]
        self.z, self.ox, self.oy = z, ox + dx_plant, oy
        self.tiles = len(tiles)

    def lattice(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return (x - self.ox) / self.d, (self.oy - y) / self.d

    def bilinear(self, c: np.ndarray, r: np.ndarray) -> np.ndarray:
        """z at lattice (c, r); NaN when a corner of the cell is NoData. A point
        on the last grid line of the array uses the cell before it."""
        i = np.clip(np.floor(r).astype(np.int64), 0, self.z.shape[0] - 2)
        j = np.clip(np.floor(c).astype(np.int64), 0, self.z.shape[1] - 2)
        tr, tc = r - i, c - j
        z = self.z
        return (z[i, j] * (1 - tc) * (1 - tr) + z[i, j + 1] * tc * (1 - tr)
                + z[i + 1, j] * (1 - tc) * tr + z[i + 1, j + 1] * tc * tr)  # fmt: skip


class Resampled:
    """The reprojected path's target grid, rebuilt here: node (R, K) at
    (K h, -R h) in the output CRS (target_grid.py's global lattice), its value
    bilinear from the source window at the node moved into the source CRS by
    pyproj, stored as float32. Source nodes are the window's pixels, PixelIsPoint
    (tiepoint at the first node). Node values are computed only where needed."""

    def __init__(self, tif: Path, out_crs: str, spacing: float):
        import pyproj

        with tifffile.TiffFile(tif) as t:
            p = t.pages[0]
            self.src = p.asarray().astype(np.float64)
            nd = p.tags.get("GDAL_NODATA")
            tie = p.tags["ModelTiepointTag"].value
            scale = p.tags["ModelPixelScaleTag"].value
            keys = p.tags["GeoKeyDirectoryTag"].value
        if nd is not None:
            self.src[self.src == float(str(nd.value).strip("\x00"))] = np.nan
        raster_type = dict(zip(keys[4::4], keys[7::4])).get(1025, 1)
        half = 0.0 if raster_type == 2 else 0.5  # PixelIsArea: centres half a pixel in
        epsg = dict(zip(keys[4::4], keys[7::4])).get(2048) or dict(zip(keys[4::4], keys[7::4])).get(3072)
        self.lon0 = tie[3] + half * scale[0]
        self.lat0 = tie[4] - half * scale[1]
        self.s = scale[0]
        assert abs(scale[0] - scale[1]) < 1e-15
        self.h = spacing
        self.tr = pyproj.Transformer.from_crs(out_crs, f"EPSG:{epsg}", always_xy=True)
        self.cache: dict[tuple[int, int], float] = {}
        self.tiles = 1
        self.raster_type = int(raster_type)
        self.plant = 0.0

    def lattice(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return x / self.h, -y / self.h

    def nodes(self, r: np.ndarray, c: np.ndarray) -> np.ndarray:
        """Values at integer nodes (r, c), float32-rounded, NaN off data."""
        lon, lat = self.tr.transform(c * self.h - self.plant, -r * self.h)
        sc, sr = (np.asarray(lon) - self.lon0) / self.s, (self.lat0 - np.asarray(lat)) / self.s
        i, j = np.floor(sr).astype(np.int64), np.floor(sc).astype(np.int64)
        inside = (i >= 0) & (j >= 0) & (i < self.src.shape[0] - 1) & (j < self.src.shape[1] - 1)
        i, j = np.where(inside, i, 0), np.where(inside, j, 0)
        tr, tc = sr - i, sc - j
        z = self.src
        v = (z[i, j] * (1 - tc) * (1 - tr) + z[i, j + 1] * tc * (1 - tr)
             + z[i + 1, j] * (1 - tc) * tr + z[i + 1, j + 1] * tc * tr)  # fmt: skip
        v = np.where(inside, v, np.nan)
        return v.astype(np.float32).astype(np.float64)

    def node_z(self, r: np.ndarray, c: np.ndarray) -> np.ndarray:
        return self.nodes(r.astype(np.float64), c.astype(np.float64))

    def bilinear(self, c: np.ndarray, r: np.ndarray) -> np.ndarray:
        i, j = np.floor(r).astype(np.int64), np.floor(c).astype(np.int64)
        key = np.concatenate([np.stack([i + a, j + b], 1) for a in (0, 1) for b in (0, 1)])
        uniq, inv = np.unique(key, axis=0, return_inverse=True)
        vals = self.nodes(uniq[:, 0].astype(np.float64), uniq[:, 1].astype(np.float64))[inv.ravel()]
        n = len(i)
        z00, z01, z10, z11 = vals[:n], vals[n : 2 * n], vals[2 * n : 3 * n], vals[3 * n :]
        tr, tc = r - i, c - j
        return z00 * (1 - tc) * (1 - tr) + z01 * tc * (1 - tr) + z10 * (1 - tc) * tr + z11 * tc * tr


def crossings(seg_lat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per segment (in lattice coordinates): every parameter t in (0, 1) where
    it crosses an integer column or row. Returns (segment index, t, is_mid):
    the crossings, then the midpoints between neighbours (ends included)."""
    seg_ids, ts = [], []
    for a, b in ((0, 2), (1, 3)):
        lo = np.minimum(seg_lat[:, a], seg_lat[:, b])
        hi = np.maximum(seg_lat[:, a], seg_lat[:, b])
        k0 = np.floor(lo).astype(np.int64) + 1
        k1 = np.ceil(hi).astype(np.int64) - 1
        n = np.maximum(k1 - k0 + 1, 0)
        sid = np.repeat(np.arange(len(seg_lat)), n)
        k = np.repeat(k0, n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
        t = (k - seg_lat[sid, a]) / (seg_lat[sid, b] - seg_lat[sid, a])
        seg_ids.append(sid)
        ts.append(t)
    sid = np.concatenate(seg_ids)
    t = np.concatenate(ts)
    keep = (t > 0) & (t < 1)
    sid, t = sid[keep], t[keep]
    order = np.lexsort((t, sid))
    sid, t = sid[order], t[order]
    # A node crossing appears twice (a column and a row at one t); keep one.
    dup = np.zeros(len(t), bool)
    dup[1:] = (sid[1:] == sid[:-1]) & (np.abs(t[1:] - t[:-1]) < 1e-12)
    sid, t = sid[~dup], t[~dup]
    # Midpoints: between neighbours on each segment, ends 0 and 1 included.
    all_sid = np.concatenate([sid, np.arange(len(seg_lat)), np.arange(len(seg_lat))])
    all_t = np.concatenate([t, np.zeros(len(seg_lat)), np.ones(len(seg_lat))])
    o = np.lexsort((all_t, all_sid))
    all_sid, all_t = all_sid[o], all_t[o]
    same = all_sid[1:] == all_sid[:-1]
    mid_sid = all_sid[1:][same]
    mid_t = ((all_t[1:] + all_t[:-1]) / 2)[same]
    is_mid = np.concatenate([np.zeros(len(t), bool), np.ones(len(mid_t), bool)])
    return np.concatenate([sid, mid_sid]), np.concatenate([t, mid_t]), is_mid


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh", type=Path)
    ap.add_argument("--domain", type=Path)
    ap.add_argument("--dem-dir", type=Path, help="DTM10 tiles (the projected path)")
    ap.add_argument("--resampled", type=Path, help="source window GeoTIFF (the reprojected path)")
    ap.add_argument("--out-crs", help="with --resampled: the mesh's CRS")
    ap.add_argument("--spacing", type=float, default=30.0, help="with --resampled")
    ap.add_argument("--lines-from-mesh", action="store_true",
                    help="check the mesh's own constraint lines instead of the input's")
    ap.add_argument("--tolerance", type=float, required=True)
    ap.add_argument("--gpkg", type=Path)
    ap.add_argument("--layer", default="corine2018")
    ap.add_argument("--snap", type=float, default=0.01, help="metres from a constraint line")
    ap.add_argument("--match", type=float, default=2e-3, help="metres: a vertex at a candidate point")
    ap.add_argument("--plant-dx", type=float, default=0.0)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    t0 = time.perf_counter()

    pts, edges = read_vtk(a.mesh)
    if a.lines_from_mesh:
        seg = np.hstack([pts[edges[:, 0], :2], pts[edges[:, 1], :2]])
    else:
        dom = shapely.from_geojson(a.domain.read_text())
        if dom.geom_type == "GeometryCollection":
            dom = shapely.union_all(list(shapely.get_parts(dom)))
        seg = segments(input_lines(dom, a.gpkg, a.layer))
    if a.resampled:
        dem = Resampled(a.resampled, a.out_crs, a.spacing)
        dem.plant = a.plant_dx  # the grid's values taken from plant_dx m further west
    else:
        x0, y0 = pts[:, 0].min(), pts[:, 1].min()
        x1, y1 = pts[:, 0].max(), pts[:, 1].max()
        dem = Dem(a.dem_dir, (x0 - 50, y0 - 50, x1 + 50, y1 + 50), a.plant_dx)

    # Alignment control: mesh vertices on DEM nodes carry the node's value.
    vc, vr = dem.lattice(pts[:, 0], pts[:, 1])
    on_node = (np.abs(vc - np.round(vc)) < 1e-6) & (np.abs(vr - np.round(vr)) < 1e-6)
    ri, ci = np.round(vr[on_node]).astype(np.int64), np.round(vc[on_node]).astype(np.int64)
    node_z = dem.node_z(ri, ci) if a.resampled else dem.z[ri, ci]
    ok = ~np.isnan(node_z)
    align = {
        "vertices": int(len(pts)),
        "vertices_on_nodes": int(on_node.sum()),
        "max_abs_z_minus_node": float(np.max(np.abs(pts[on_node, 2][ok] - node_z[ok]))) if ok.any() else None,
    }

    # Points: crossings and midpoints on the input lines.
    sc, sr = dem.lattice(seg[:, 0], seg[:, 1])
    ec, er = dem.lattice(seg[:, 2], seg[:, 3])
    sid, t, is_mid = crossings(np.column_stack([sc, sr, ec, er]))
    px = seg[sid, 0] + t * (seg[sid, 2] - seg[sid, 0])
    py = seg[sid, 1] + t * (seg[sid, 3] - seg[sid, 1])
    pc, pr = dem.lattice(px, py)
    z_dem = dem.bilinear(pc, pr)

    # The mesh's z on its nearest constraint line.
    a_xy, b_xy = pts[edges[:, 0], :2], pts[edges[:, 1], :2]
    tree = shapely.STRtree(shapely.linestrings(np.stack([a_xy, b_xy], axis=1)))
    q = shapely.points(px, py)
    idx, dist = tree.query_nearest(q, max_distance=a.snap, return_distance=True, all_matches=False)
    hit = np.full(len(px), -1, dtype=np.int64)
    hit[idx[0]] = idx[1]
    d = np.full(len(px), np.inf)
    d[idx[0]] = dist
    on = hit >= 0
    e = edges[hit[on]]
    pa, pb = pts[e[:, 0]], pts[e[:, 1]]
    ab = pb[:, :2] - pa[:, :2]
    s = np.clip(((px[on] - pa[:, 0]) * ab[:, 0] + (py[on] - pa[:, 1]) * ab[:, 1]) / (ab ** 2).sum(1), 0, 1)
    z_mesh = np.full(len(px), np.nan)
    z_mesh[on] = pa[:, 2] + s * (pb[:, 2] - pa[:, 2])

    err = np.abs(z_mesh - z_dem)
    nodata = np.isnan(z_dem)
    valid = on & ~nodata
    over = valid & (err > a.tolerance + EPS)
    out: dict[str, object] = {
        "mesh": str(a.mesh), "tolerance": a.tolerance, "plant_dx": a.plant_dx,
        "dem_tiles": dem.tiles, "lines": "mesh" if a.lines_from_mesh else "input", "input_segments": int(len(seg)),
        "input_length_m": float(np.hypot(seg[:, 2] - seg[:, 0], seg[:, 3] - seg[:, 1]).sum()),
        "mesh_constraint_edges": int(len(edges)), "alignment": align,
    }  # fmt: skip
    for name, m in (("crossings", ~is_mid), ("midpoints", is_mid)):
        sel_over = over & m
        worst = np.argsort(-np.where(valid & m, err, -1))[:5]
        out[name] = {
            "points": int(m.sum()),
            "off_the_lines": int((m & ~on).sum()),
            "no_data": int((m & on & nodata).sum()),
            "checked": int((valid & m).sum()),
            "over": int(sel_over.sum()),
            "max_error_m": float(err[valid & m].max()) if (valid & m).any() else None,
            "max_over_m": float(err[sel_over].max()) if sel_over.any() else None,
            "worst": [
                {"x": float(px[i]), "y": float(py[i]), "err": float(err[i]), "z_mesh": float(z_mesh[i]),
                 "z_dem": float(z_dem[i])}
                for i in worst if valid[i] and m[i]
            ],
        }  # fmt: skip
    # Q1: which candidate points are vertices of the mesh (within --match m).
    cv = np.unique(edges)
    vtree = shapely.STRtree(shapely.points(pts[cv, 0], pts[cv, 1]))
    vi, _ = vtree.query_nearest(q, max_distance=a.match, return_distance=True, all_matches=False)
    is_vertex = np.zeros(len(px), bool)
    is_vertex[vi[0]] = True
    out["vertices_at_candidates"] = {
        "constraint_vertices": int(len(cv)),
        "at_crossings": int((is_vertex & ~is_mid).sum()),
        "at_midpoints": int((is_vertex & is_mid).sum()),
        "match_m": a.match,
    }
    out["max_snap_distance_m"] = float(d[on].max()) if on.any() else None
    out["seconds"] = round(time.perf_counter() - t0, 1)
    text = json.dumps(out, indent=1)
    if a.json:
        a.json.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
