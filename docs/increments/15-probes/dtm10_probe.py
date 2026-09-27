"""Increment 15's measurements on Ola's DTM10 archive (254 UTM33 tiles).

    python docs/increments/15-probes/dtm10_probe.py headers DIR   # N1, N2, N4: headers, lattice
    python docs/increments/15-probes/dtm10_probe.py seams DIR     # N3: do overlaps agree?

Header-only for `headers`, through io/geotiff.py's own private checks, so a
refusal here is the refusal `decode_dem` would raise. A measurement script, not
production code; nothing imports it.
"""

from __future__ import annotations

import collections
import sys
from pathlib import Path

import numpy as np
import tifffile

from tin_engine.io import geotiff as g


def metas(directory: Path) -> dict[str, tuple]:
    out, refused = {}, collections.Counter()
    for path in sorted(directory.glob("*.tif")):
        with path.open("rb") as stream, tifffile.TiffFile(stream) as tif:
            page, pages = tif.pages.first, tuple(tif.pages)
            try:
                g._single_page(pages)
                dtype = g._check_page(page)
                tie, scale = g._georeferencing(page)
                keys = tif.geotiff_metadata or {}
                x, y, dx, dy, area = g._placement(tie, scale, keys)
                epsg = g._projected_epsg(keys)
                nodata, source = g._nodata(page, dtype, None)
            except g.GeoTiffError as error:
                refused[str(error)[:90]] += 1
                continue
            out[path.name] = (x, y, dx, dy, page.imagelength, page.imagewidth, epsg, area,
                              nodata, source, str(dtype), page.compression.name, page.is_tiled)
    for text, n in refused.items():
        print(f"REFUSED x{n}: {text}")
    return out


def headers(directory: Path) -> None:
    m = metas(directory)
    print(f"{len(m)} tiles read")
    for i, label in ((6, "epsg"), (2, "dx"), (3, "dy"), (7, "area"), (8, "nodata"),
                     (9, "nodata source"), (10, "dtype"), (11, "compression"),
                     (12, "tiled"), (4, "rows"), (5, "cols")):
        print(f"{label}: {collections.Counter(v[i] for v in m.values()).most_common(6)}")
    x0 = min(v[0] for v in m.values())
    y0 = max(v[1] for v in m.values())
    dx = next(iter(m.values()))[2]
    off = [((v[0] - x0) / dx, (y0 - v[1]) / dx) for v in m.values()]
    worst = max(max(abs(a - round(a)), abs(b - round(b))) for a, b in off)
    print(f"reference node x {x0}, y {y0}; worst off-integer offset {worst:.1e} cells")
    boxes = {k: (round((y0 - v[1]) / dx), round((y0 - v[1]) / dx) + v[4],
                 round((v[0] - x0) / dx), round((v[0] - x0) / dx) + v[5]) for k, v in m.items()}
    rows = collections.Counter()
    names = sorted(boxes)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            ra, rb = boxes[a], boxes[b]
            h = min(ra[1], rb[1]) - max(ra[0], rb[0])
            w = min(ra[3], rb[3]) - max(ra[2], rb[2])
            if h > 0 and w > 0:
                rows[(min(h, w), "corner" if h < 1000 and w < 1000 else "edge")] += 1
    print("overlapping pairs by (overlap width in nodes, kind):", sorted(rows.items()))
    span = (max(b[1] for b in boxes.values()), max(b[3] for b in boxes.values()))
    print(f"union bbox {span[0]} x {span[1]} nodes = {span[0] * span[1] / 1e9:.2f} G, "
          f"{span[0] * span[1] * 4 / 2**30:.1f} GiB float32")
    covered = sum(v[4] * v[5] for v in m.values())
    print(f"sum of tile nodes {covered / 1e9:.2f} G ({covered * 4 / 2**30:.1f} GiB float32)")


def seams(directory: Path) -> None:
    m = metas(directory)
    x0 = min(v[0] for v in m.values())
    y0 = max(v[1] for v in m.values())
    dx = next(iter(m.values()))[2]
    boxes = {k: (round((y0 - v[1]) / dx), round((y0 - v[1]) / dx) + v[4],
                 round((v[0] - x0) / dx), round((v[0] - x0) / dx) + v[5]) for k, v in m.items()}
    names, done = sorted(boxes), 0
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            ra, rb = boxes[a], boxes[b]
            r0, r1 = max(ra[0], rb[0]), min(ra[1], rb[1])
            c0, c1 = max(ra[2], rb[2]), min(ra[3], rb[3])
            if r1 - r0 <= 0 or c1 - c0 <= 0 or min(r1 - r0, c1 - c0) > 1000 or done >= 4:
                continue
            za = tifffile.imread(directory / a)[r0 - ra[0] : r1 - ra[0], c0 - ra[2] : c1 - ra[2]]
            zb = tifffile.imread(directory / b)[r0 - rb[0] : r1 - rb[0], c0 - rb[2] : c1 - rb[2]]
            nd = m[a][8]
            va, vb = za != nd, zb != nd
            both = va & vb
            diff = np.abs(za[both].astype(float) - zb[both])
            print(f"{a} | {b}: {za.shape} overlap, both valid {int(both.sum())}, one valid "
                  f"{int((va ^ vb).sum())}, equal {int((za[both] == zb[both]).sum())}, "
                  f"max |diff| {diff.max() if diff.size else None}")
            if za.shape[0] > 1 and za.shape[1] > 1:  # the probe can fail: shift by one
                s = np.abs(za[1:, :].astype(float) - zb[:-1, :]) if za.shape[0] <= za.shape[1] \
                    else np.abs(za[:, 1:].astype(float) - zb[:, :-1])
                print(f"   shifted by one node: max |diff| {s.max():.2f}, "
                      f"equal {(s == 0).mean():.1%}")
            done += 1


if __name__ == "__main__":
    {"headers": headers, "seams": seams}[sys.argv[1]](Path(sys.argv[2]))
