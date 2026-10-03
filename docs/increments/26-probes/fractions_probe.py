"""Throwaway probe for docs/increments/26-land-cover.md: exact MapBiomas class
areas per triangle of a real basin mesh, and the cutoff with and without a ledger.

Not production code. Run with the repository's venv (needs vtk, the `viewer`
extra, and imagecodecs, the `codecs` extra), after compiling the helper:

    c++ -std=c++20 -O2 -o $WORK/overlap_ledger docs/increments/26-probes/overlap_ledger.cpp
    python docs/increments/26-probes/fractions_probe.py MESH.vtk $WORK [--year 2025]

MESH.vtk is a rasputin mesh (its FieldData `crs` names the CRS). $WORK holds
the downloaded MapBiomas header and tiles (reused between runs) and the
intermediate arrays. Prints one JSON object per run.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import io
import json
import logging
import os
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import tifffile
from pyproj import CRS, Geod, Transformer

logging.disable(logging.CRITICAL)
URL = (
    "https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection11/"
    "lulc/coverage/brazil_coverage/brazil_coverage-col11_{year}.tif"
)
CROPS = {
    39: "soybean",
    20: "sugar cane",
    40: "rice",
    62: "cotton",
    41: "other temporary",
    46: "coffee",
    47: "citrus",
    35: "palm oil",
    48: "other perennial",
}


def read_mesh(path: Path) -> tuple[np.ndarray, np.ndarray, str]:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    r = vtk.vtkPolyDataReader()
    r.SetFileName(str(path))
    r.ReadAllFieldsOn()
    r.Update()
    pd = r.GetOutput()
    xy = vtk_to_numpy(pd.GetPoints().GetData())[:, :2].astype(np.float64)
    conn = vtk_to_numpy(pd.GetPolys().GetConnectivityArray()).astype(np.int64)
    tri = conn.reshape(-1, 3)
    crs = pd.GetFieldData().GetAbstractArray("crs").GetValue(0)
    return xy, tri, crs


def get_range(url: str, start: int, stop: int) -> bytes:
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{stop - 1}"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read()
        except OSError:
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def class_window(
    work: Path, year: int, c0: int, r0: int, c1: int, r1: int
) -> tuple[np.ndarray, object]:
    url = URL.format(year=year)
    head = work / f"c11_{year}_head.bin"
    if not head.exists():
        head.write_bytes(get_range(url, 0, 16 << 20))
    page = tifffile.TiffFile(io.BytesIO(head.read_bytes())).pages[0]
    tw, th = page.tilewidth, page.tilelength
    ncols = -(-page.shape[1] // tw)
    out = np.zeros((r1 - r0, c1 - c0), np.uint8)
    tiles = [
        (tr, tc)
        for tr in range(r0 // th, (r1 - 1) // th + 1)
        for tc in range(c0 // tw, (c1 - 1) // tw + 1)
    ]
    tdir = work / f"tiles_{year}"
    tdir.mkdir(exist_ok=True)

    def one(rc: tuple[int, int]) -> tuple[int, int, np.ndarray | None]:
        tr, tc = rc
        idx = tr * ncols + tc
        n = page.databytecounts[idx]
        if n == 0:
            return tr, tc, None
        f = tdir / f"{idx}.bin"
        if not f.exists():
            o = page.dataoffsets[idx]
            f.write_bytes(get_range(url, o, o + n))
        a = page.decode(f.read_bytes(), idx)[0]
        return tr, tc, np.asarray(a).reshape(th, tw)

    with concurrent.futures.ThreadPoolExecutor(8) as ex:
        for tr, tc, a in ex.map(one, tiles):
            if a is None:
                continue
            rr0, cc0 = tr * th, tc * tw
            ys, xs = max(rr0, r0), max(cc0, c0)
            ye, xe = min(rr0 + th, r1), min(cc0 + tw, c1)
            out[ys - r0 : ye - r0, xs - c0 : xe - c0] = a[ys - rr0 : ye - rr0, xs - cc0 : xe - cc0]
    return out, page


def hilbert_keys(x: np.ndarray, y: np.ndarray, order: int = 16) -> np.ndarray:
    """Hilbert index of integer coordinates in [0, 2^order), vectorised xy2d."""
    x = x.astype(np.int64).copy()
    y = y.astype(np.int64).copy()
    d = np.zeros_like(x)
    s = 1 << (order - 1)
    while s > 0:
        rx = ((x & s) > 0).astype(np.int64)
        ry = ((y & s) > 0).astype(np.int64)
        d += s * s * ((3 * rx) ^ ry)
        flip = ry == 0
        swap_x = np.where(flip & (rx == 1), s - 1 - x, x)
        swap_y = np.where(flip & (rx == 1), s - 1 - y, y)
        x, y = np.where(flip, swap_y, swap_x), np.where(flip, swap_x, swap_y)
        s >>= 1
    return d


def load_csr(work: Path, prefix: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.fromfile(work / f"{prefix}.csr.u32" if prefix != "csr" else work / "csr.u32", np.uint32),
        np.fromfile(work / f"{prefix}.cls.u8" if prefix != "csr" else work / "cls.u8", np.uint8),
        np.fromfile(
            work / f"{prefix}.area.f64" if prefix != "csr" else work / "area.f64", np.float64
        ),
    )


def dense_by_region(
    off: np.ndarray, cls: np.ndarray, ar: np.ndarray, region: np.ndarray
) -> dict[tuple[int, int], float]:
    tri_of = np.repeat(np.arange(off.size - 1), np.diff(off).astype(np.int64))
    key = region[tri_of].astype(np.int64) * 256 + cls
    u, inv = np.unique(key, return_inverse=True)
    s = np.zeros(u.size)
    np.add.at(s, inv, ar)
    return dict(zip(u.tolist(), s.tolist(), strict=True))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh", type=Path)
    ap.add_argument("work", type=Path)
    ap.add_argument("--year", type=int, default=2025)
    ap.add_argument(
        "--cap", type=int, default=0, help="at most this many classes per triangle (0: none)"
    )
    ap.add_argument(
        "--chosen",
        action="store_true",
        help="only the design's setting: 5 %%, 1 ha, renorm and present",
    )
    args = ap.parse_args()
    work: Path = args.work
    run = work / args.mesh.stem
    run.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    xy, tri, crs_text = read_mesh(args.mesh)
    inv = Transformer.from_crs(CRS.from_user_input(crs_text), CRS.from_epsg(4326), always_xy=True)
    lon, lat = inv.transform(xy[:, 0], xy[:, 1])
    head = work / f"c11_{args.year}_head.bin"
    if not head.exists():
        head.write_bytes(get_range(URL.format(year=args.year), 0, 16 << 20))
    page = tifffile.TiffFile(io.BytesIO(head.read_bytes())).pages[0]
    sx, sy = page.tags["ModelPixelScaleTag"].value[:2]
    x0, y0 = page.tags["ModelTiepointTag"].value[3:5]
    col = (np.asarray(lon) - x0) / sx
    row = (y0 - np.asarray(lat)) / sy
    c0, r0 = int(np.floor(col.min())) - 1, int(np.floor(row.min())) - 1
    c1, r1 = int(np.ceil(col.max())) + 2, int(np.ceil(row.max())) + 2
    win, _ = class_window(work, args.year, c0, r0, c1, r1)
    geod = Geod(ellps="WGS84")
    rows = np.arange(r0, r1)
    lat_top = y0 - rows * sy
    roww = np.array(
        [
            abs(geod.polygon_area_perimeter([0, sx, sx, 0], [a, a, a - sy, a - sy])[0])
            for a in lat_top
        ]
    )
    tv = np.stack([col[tri] - c0, row[tri] - r0], axis=-1).reshape(-1, 6)
    tv.astype(np.float64).tofile(run / "tri.f64")
    win.tofile(run / "win.u8")
    np.array(win.shape, np.int64).tofile(run / "win.shape")
    roww.tofile(run / "roww.f64")
    helper = str(work / "overlap_ledger")
    t1 = time.time()
    subprocess.run([helper, "overlap", str(run)], check=True)
    t2 = time.time()
    off, cls, ar = load_csr(run, "csr")
    T = tri.shape[0]  # noqa: N806 (the triangle count, as in the increment's text)
    A = np.add.reduceat(ar, off[:-1].astype(np.int64)) if ar.size else np.zeros(T)  # noqa: N806
    nnz = np.diff(off.astype(np.int64))
    # Hilbert order of centroids in the mesh CRS, on a 2^16 grid over the mesh's box.
    cx, cy = xy[tri, 0].mean(axis=1), xy[tri, 1].mean(axis=1)
    side = max(cx.max() - cx.min(), cy.max() - cy.min()) * (1 + 1e-12)
    gx = np.floor((cx - cx.min()) / side * 65535).astype(np.int64)
    gy = np.floor((cy - cy.min()) / side * 65535).astype(np.int64)
    key = hilbert_keys(gx, gy)
    order = np.lexsort((cy, cx, key)).astype(np.uint32)
    order.tofile(run / "order.u32")
    exact_total = dense_by_region(off, cls, ar, np.zeros(T, np.int64))
    class_area = {k % 256: v for k, v in exact_total.items()}
    starts = off[:-1].astype(np.int64)
    seg_max = np.maximum.reduceat(ar, starts)
    tri_of = np.repeat(np.arange(T), nnz)
    first_max = np.full(T, -1, np.int64)
    hit = np.flatnonzero(ar == seg_max[tri_of])[::-1]
    first_max[tri_of[hit]] = hit
    dom = cls[first_max]
    crop_mask = np.isin(dom, list(CROPS))
    report: dict[str, object] = {
        "mesh": args.mesh.name,
        "triangles": T,
        "year": args.year,
        "seconds_overlap": round(t2 - t1, 1),
        "seconds_prepare": round(t1 - t0, 1),
        "window_cells": int(win.size),
        "mean_cell_m2": float(roww.mean()),
        "mean_triangle_m2": float(A.mean()),
        "mean_triangle_cells": float(A.mean() / roww.mean()),
        "triangle_m2_p50_p90_p99_max": [float(np.percentile(A, q)) for q in (50, 90, 99, 100)],
        "median_triangle_m2_dominant_crop": float(np.median(A[crop_mask]))
        if crop_mask.any()
        else None,
        "median_triangle_m2_dominant_other": float(np.median(A[~crop_mask])),
        "triangles_dominant_crop": int(crop_mask.sum()),
        "exact_nnz_mean": float(nnz.mean()),
        "exact_nnz_max": int(nnz.max()),
        "exact_nnz_hist": np.bincount(nnz, minlength=1)[:12].tolist(),
        "class_area_km2": {str(k): round(v / 1e6, 3) for k, v in sorted(class_area.items())},
        "variants": [],
    }
    tri_all = np.repeat(np.arange(T), nnz)
    regions = {}
    for s_km in (1, 5, 25):
        s = s_km * 1000.0
        rx = np.floor((cx - cx.min()) / s).astype(np.int64)
        ry = np.floor((cy - cy.min()) / s).astype(np.int64)
        _, rinv = np.unique(rx * 100000 + ry, return_inverse=True)
        reg_area = np.bincount(rinv, weights=A)
        full = reg_area >= 0.95 * s * s
        ex = np.zeros((reg_area.size, 256))
        np.add.at(ex, (rinv[tri_all], cls), ar)
        regions[s_km] = (rinv, reg_area, full, ex)
    env = {**os.environ, "RASPUTIN_PROBE_CAP": str(args.cap)}
    report["cap"] = args.cap
    for c in (0.05,) if args.chosen else (0.01, 0.05, 0.10):
        for m in (1e4,) if args.chosen else (float("inf"), 1e4, 1e5):
            for rule in ("renorm", "present") if args.chosen else ("renorm", "present", "any"):
                res = subprocess.run(
                    [helper, "ledger", str(run), str(c), str(m), rule],
                    check=True,
                    capture_output=True,
                    text=True,
                    env=env,
                )
                led = json.loads(res.stdout)
                ooff, ocls, oar = load_csr(run, f"out_{rule}_{c:g}_{m:g}")
                onnz = np.diff(ooff.astype(np.int64))
                out_total = dense_by_region(ooff, ocls, oar, np.zeros(T, np.int64))
                cls_err = {}
                for k, v in class_area.items():
                    o = out_total.get(k, 0.0)
                    cls_err[str(k)] = [
                        round((o - v) / 1e6, 4),
                        round(100 * (o - v) / v, 3) if v > 0 else None,
                    ]
                worst_cls = max(cls_err.items(), key=lambda kv: abs(kv[1][1] or 0))
                otri = np.repeat(np.arange(T), onnz)
                regional = {}
                for s_km, (rinv, reg_area, full, ex) in regions.items():
                    od = np.zeros_like(ex)
                    np.add.at(od, (rinv[otri], ocls), oar)
                    mis = (
                        50.0 * np.abs(od - ex).sum(axis=1) / reg_area
                    )  # % of the region in a wrong class
                    mf = mis[full]
                    regional[f"{s_km}km"] = (
                        [
                            round(float(np.percentile(mf, 50)), 3),
                            round(float(np.percentile(mf, 95)), 3),
                            round(float(mf.max()), 3),
                            int(full.sum()),
                        ]
                        if full.any()
                        else None
                    )
                ledger_worst = max((w for w, _ in led["worst"].values()), default=0.0)
                remainder = max((abs(r) for _, r in led["worst"].values()), default=0.0)
                vanished = [
                    k for k, v in class_area.items() if v > 0 and out_total.get(k, 0.0) < 0.5 * v
                ]
                report["variants"].append(
                    {
                        "rule": rule,
                        "cutoff": c,
                        "min_m2": m,
                        "nnz_mean": round(float(onnz.mean()), 3),
                        "nnz_max": int(onnz.max()),
                        "worst_class_error": [worst_cls[0], *worst_cls[1]],
                        "crop_errors_pct": {
                            CROPS[k]: cls_err[str(k)][1] for k in CROPS if str(k) in cls_err
                        },
                        "water_error_pct": cls_err.get("33", [None, None])[1],
                        "classes_below_half": vanished,
                        "regional_misplaced_pct_p50_p95_max_n": regional,
                        "ledger_max_abs_m2": round(ledger_worst, 1),
                        "ledger_max_in_mean_triangles": round(ledger_worst / float(A.mean()), 3),
                        "end_remainder_max_m2": round(remainder, 3),
                    }
                )
    json.dump(report, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
