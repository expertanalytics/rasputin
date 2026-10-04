"""Throwaway probe for docs/increments/26-land-cover.md: MapBiomas water bodies
in a domain, as 8-connected components, traced on cell edges and reduced with
increment 22's `reduce_ring` (pinches opened by 1e-9 first, as in
docs/research/raster-to-vector-probes/pinch_strip_probe.py).

Not production code. Needs the helper of fractions_probe.py (`label` mode) in
$WORK and the repository's venv with `_core` built:

    python docs/increments/26-probes/water_probe.py OUTLINE.geojson MESH.vtk $WORK [--classes 33]

OUTLINE is the domain in longitude/latitude; MESH gives the computation CRS
(its FieldData `crs`). Prints one JSON object.
"""

from __future__ import annotations

import argparse
import io
import json
import math
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import shapely
import tifffile
from pyproj import CRS, Transformer

sys.path.insert(0, str(Path(__file__).parent))
from fractions_probe import URL, class_window, get_range, read_mesh

from tin_engine._core import reduce_ring

E, N, W, S = (1, 0), (0, 1), (-1, 0), (0, -1)


def rings_of(mask: np.ndarray) -> list[list[tuple[int, int]]]:
    """Closed rings of cell corners (x = col, y = -row), water on the left,
    collinear corners dropped; at a corner with two ways out, the right turn
    (so water is 8-connected and the ring passes a pinch twice)."""
    m = np.pad(mask, 1)
    out: dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int]]]] = defaultdict(list)
    rr, cc = np.nonzero(m)
    for r, c in zip(rr.tolist(), cc.tolist(), strict=True):
        x, y = c, -r
        if not m[r + 1, c]:
            out[(x, y - 1)].append(((x + 1, y - 1), E))
        if not m[r - 1, c]:
            out[(x + 1, y)].append(((x, y), W))
        if not m[r, c - 1]:
            out[(x, y)].append(((x, y - 1), S))
        if not m[r, c + 1]:
            out[(x + 1, y - 1)].append(((x + 1, y), N))
    rings = []
    while out:
        start = next((p for p, e in out.items() if len(e) == 1), next(iter(out)))
        p, d = start, None
        verts: list[tuple[int, int]] = []
        dirs: list[tuple[int, int]] = []
        while True:
            cands = out[p]
            if d is None:
                k = 0
            else:
                pref = [(d[1], -d[0]), d, (-d[1], d[0])]  # right, straight, left
                k = min(range(len(cands)), key=lambda i: pref.index(cands[i][1]))
            q, d2 = cands.pop(k)
            if not cands:
                del out[p]
            verts.append(p)
            dirs.append(d2)
            p, d = q, d2
            if p == start:
                break
        # keep a vertex only where the direction changes (incoming != outgoing)
        ring = [v for i, v in enumerate(verts) if dirs[i - 1] != dirs[i]]
        if len(ring) >= 4:
            rings.append(ring)
    return rings


def signed_area(r: np.ndarray) -> float:
    """Shoelace with one cross product per edge, summed exactly (math.fsum):
    two separate dot products cancel catastrophically on a 346 km ring."""
    x, y = r[:, 0], r[:, 1]
    return 0.5 * math.fsum((x * np.roll(y, -1) - np.roll(x, -1) * y).tolist())


def open_pinches(r: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """Move each repeated vertex by eps towards the middle of its two neighbours."""
    r = r.copy()
    _, inv, cnt = np.unique(r, axis=0, return_inverse=True, return_counts=True)
    for i in np.flatnonzero(cnt[inv.ravel()] > 1):
        mid = 0.5 * (r[i - 1] + r[(i + 1) % len(r)])
        d = mid - r[i]
        n = np.hypot(*d)
        if n > 0:
            r[i] = r[i] + eps * d / n
    return r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("outline", type=Path)
    ap.add_argument("mesh", type=Path)
    ap.add_argument("work", type=Path)
    ap.add_argument("--year", type=int, default=2025)
    ap.add_argument("--classes", default="33")
    args = ap.parse_args()
    work: Path = args.work
    run = work / f"water_{args.outline.stem}_{args.classes.replace(',', '-')}"
    run.mkdir(parents=True, exist_ok=True)
    gj = json.loads(args.outline.read_text())
    outline = shapely.union_all([shapely.geometry.shape(f["geometry"]) for f in gj["features"]])
    head = work / f"c11_{args.year}_head.bin"
    if not head.exists():
        head.write_bytes(get_range(URL.format(year=args.year), 0, 16 << 20))
    page = tifffile.TiffFile(io.BytesIO(head.read_bytes())).pages[0]
    sx, sy = page.tags["ModelPixelScaleTag"].value[:2]
    x0, y0 = page.tags["ModelTiepointTag"].value[3:5]
    lon0, lat0, lon1, lat1 = outline.bounds
    c0, r0 = int((lon0 - x0) / sx) - 2, int((y0 - lat1) / sy) - 2
    c1, r1 = int((lon1 - x0) / sx) + 3, int((y0 - lat0) / sy) + 3
    win, _ = class_window(work, args.year, c0, r0, c1, r1)
    win.tofile(run / "win.u8")
    np.array(win.shape, np.int64).tofile(run / "win.shape")
    subprocess.run([str(work / "overlap_ledger"), "label", str(run), args.classes], check=True)
    lab = np.fromfile(run / "labels.i32", np.int32).reshape(win.shape)
    sizes = np.fromfile(run / "sizes.i64", np.int64)
    ids = np.arange(1, sizes.size)
    # A body belongs to the unit when any of its cells' centres is inside the
    # outline (design review, round 1: the first cell in raster order missed
    # Tres Marias, whose first cell lies north of unit 769).
    wr, wc = np.nonzero(lab)
    wl = lab[wr, wc]
    shapely.prepare(outline)
    cell_in = shapely.contains_xy(outline, x0 + (c0 + wc + 0.5) * sx, y0 - (r0 + wr + 0.5) * sy)
    cells_inside = np.bincount(wl[cell_in], minlength=sizes.size)
    inside = cells_inside[ids] > 0
    cell_m2 = 868.0  # the window's mean cell, from fractions_probe on this unit
    _, _, crs_text = read_mesh(args.mesh)
    fwd = Transformer.from_crs(CRS.from_epsg(4326), CRS.from_user_input(crs_text), always_xy=True)
    area_m2 = sizes[ids] * cell_m2
    report: dict[str, object] = {
        "outline": args.outline.name,
        "classes": args.classes,
        "year": args.year,
        "window_cells": int(win.size),
        "components_inside": int(inside.sum()),
        "water_km2_inside": round(float(area_m2[inside].sum()) / 1e6, 3),
        "water_km2_cells_inside_outline": round(float(cells_inside.sum()) * cell_m2 / 1e6, 3),
        "cell_area_m2_used": cell_m2,
        "count_at_least_ha": {
            str(ha): int(((area_m2 >= ha * 1e4) & inside).sum())
            for ha in (0.5, 1.44, 5, 10, 50, 100, 1000)
        },
        "area_share_at_least_ha": {
            str(ha): round(
                float(area_m2[(area_m2 >= ha * 1e4) & inside].sum() / area_m2[inside].sum()), 4
            )
            for ha in (0.5, 1.44, 5, 10, 50, 100, 1000)
        },
        "reduction": {},
    }
    # Trace and reduce every component of at least 4 x (60 m)^2 = 1.44 ha.
    keep = ids[(area_m2 >= 14_400) & inside]
    t0 = time.time()
    fine, pinches, outers = 0, 0, []
    nl = sizes.size
    largest = int(keep[np.argmax(sizes[keep])])
    rmin = np.full(nl, 1 << 40)
    cmin = np.full(nl, 1 << 40)
    rmax = np.full(nl, -1)
    cmax = np.full(nl, -1)
    np.minimum.at(rmin, wl, wr)
    np.minimum.at(cmin, wl, wc)
    np.maximum.at(rmax, wl, wr)
    np.maximum.at(cmax, wl, wc)
    for lid in keep.tolist():
        ra, ca, rb, cb = int(rmin[lid]), int(cmin[lid]), int(rmax[lid]) + 1, int(cmax[lid]) + 1
        mask = lab[ra:rb, ca:cb] == lid
        rings = rings_of(mask)
        if lid == largest:
            at_edge = ra == 0 or ca == 0 or rb == lab.shape[0] or cb == lab.shape[1]
            report["largest_body"] = {
                "cells": int(sizes[lid]),
                "km2": round(float(sizes[lid]) * cell_m2 / 1e6, 1),
                "box_rows_cols": [rb - ra, cb - ca],
                "box_km_approx_ns_ew": [round((rb - ra) * 0.0309, 1), round((cb - ca) * 0.030, 1)],
                "holes": sum(1 for g in rings if signed_area(np.array(g, float)) < 0),
                "touches_window_edge": bool(at_edge),
            }
        for ring in rings:
            a = np.array(ring, float)
            if signed_area(a) <= 0:
                continue  # holes (islands) are not reduced in this probe
            col = a[:, 0] - 1 + ca + c0
            row = -a[:, 1] - 1 + ra + r0
            lo, la = x0 + col * sx, y0 - row * sy
            mx, my = fwd.transform(lo, la)
            xy = np.stack([mx, my], axis=1)
            fine += len(xy)
            pinches += len(xy) - len(np.unique(xy, axis=0))
            outers.append(xy)
    t1 = time.time()
    report["traced"] = {
        "lakes": len(keep),
        "outer_rings": len(outers),
        "fine_vertices": fine,
        "pinch_vertices": pinches,
        "seconds": round(t1 - t0, 1),
    }
    big = max(range(len(outers)), key=lambda i: len(outers[i]))
    report["largest_ring"] = {"fine": len(outers[big])}
    for tol in (30.0, 60.0, 120.0):
        tot, crossings, worst_rel, worst_ring = 0, 0, 0.0, {}
        for i, xy in enumerate(outers):
            o = xy.mean(axis=0)
            r = open_pinches(xy - o)
            t_ring = time.time()
            out = reduce_ring(r, tol, np.zeros((0, 2)))
            if i == big:
                report["largest_ring"][f"{tol:g} m"] = {
                    "reduced": len(out.ring),
                    "seconds": round(time.time() - t_ring, 2),
                }
            red = np.asarray(out.ring)
            tot += len(red)
            crossings += out.rejected_crossing
            rel = abs(signed_area(red) - signed_area(r)) / signed_area(r)
            if rel > worst_rel:
                worst_rel = rel
                worst_ring = {
                    "area_m2": signed_area(r),
                    "change_m2": signed_area(red) - signed_area(r),
                    "fine": len(r),
                    "reduced": len(red),
                    "collapses": out.collapses,
                    "half_extent_m": float(np.abs(r).max()),
                }
        report["reduction"][f"{tol:g} m"] = {
            "vertices": tot,
            "rejected_crossing": crossings,
            "worst_relative_area_change": worst_rel,
            "worst_ring": worst_ring,
        }
    json.dump(report, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
