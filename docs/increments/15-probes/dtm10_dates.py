"""Do DTM10 overlaps agree when the two tiles share an export date?

Evidence for "Ruled by Ola", Q1 revised (docs/increments/15-dem-mosaic.md).
The export date of a tile is the modification date of its side files
(`.tif.aux.xml`, else `.tfw`; the two agree on all 254 tiles). 143 of the
254 `.tif` files share one date, 2021-12-10, later than their side files, so
the `.tif` date does not tell the exports apart. The TIFF tags carry no date.

Usage: python dtm10_dates.py [ARCHIVE] [N_PAIRS] [SEED]
Prints, for N random neighbour pairs on one lattice, whether they share an
export date and how their jointly valid overlap nodes compare: exact, within
1 mm, or beyond.
"""

import datetime
import glob
import os
import random
import sys

import numpy as np
import tifffile

NODATA = -32767.0


def export_date(tif):
    for side in (tif + ".aux.xml", tif[:-4] + ".tfw"):
        if os.path.exists(side):
            return datetime.date.fromtimestamp(os.path.getmtime(side))
    return None


def origin(tif):
    with tifffile.TiffFile(tif) as f:
        page = f.pages[0]
        tie = page.tags["ModelTiepointTag"].value
        return tie[3], tie[4], page.shape


def overlap(a, b, info):
    (xa, ya, sa), (xb, yb, sb) = info[a], info[b]
    x0, x1 = max(xa, xb), min(xa + sa[1] * 10, xb + sb[1] * 10)
    y0, y1 = max(ya - sa[0] * 10, yb - sb[0] * 10), min(ya, yb)
    if x1 <= x0 or y1 <= y0:
        return None
    arr_a, arr_b = tifffile.imread(a), tifffile.imread(b)

    def window(arr, xm, ym):
        r0, r1 = round((ym - y1) / 10), round((ym - y0) / 10)
        c0, c1 = round((x0 - xm) / 10), round((x1 - xm) / 10)
        return arr[r0:r1, c0:c1].astype(np.float64)

    wa, wb = window(arr_a, xa, ya), window(arr_b, xb, yb)
    ok = (wa != NODATA) & (wb != NODATA)
    return np.abs(wa[ok] - wb[ok])


def main():
    archive = sys.argv[1] if len(sys.argv) > 1 else "../rasputin_data/DTM10_UTM33_20220924"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    tiles = sorted(glob.glob(f"{archive}/*_10m_z33.tif"))
    info = {t: origin(t) for t in tiles}
    pairs = [
        (a, b)
        for i, a in enumerate(tiles)
        for b in tiles[i + 1 :]
        if (info[a][0] - info[b][0]) % 10 == 0
        and (info[a][1] - info[b][1]) % 10 == 0
        and abs(info[a][0] - info[b][0]) <= 50510
        and abs(info[a][1] - info[b][1]) <= 50510
    ]
    random.seed(seed)
    rows = []
    for a, b in random.sample(pairs, min(n, len(pairs))):
        d = overlap(a, b, info)
        if d is None or d.size == 0:
            continue
        same = export_date(a) == export_date(b)
        rows.append((same, int((d > 0).sum()), int((d >= 1e-3).sum()), float(d.max()), a, b))
    for same in (True, False):
        group = [r for r in rows if r[0] == same]
        exact = sum(r[1] == 0 for r in group)
        within = sum(r[2] == 0 for r in group)
        worst = max((r[3] for r in group), default=0.0)
        label = "same date" if same else "different dates"
        print(
            f"{label:16s} pairs {len(group):3d}  exact {exact:3d}"
            f"  within 1 mm {within:3d}  worst {worst:.3f} m"
        )
    same_bad = [r for r in rows if r[0] and r[2] > 0]
    worst_rows = sorted(rows, key=lambda r: -r[3])[:8]
    for same, _n_any, n_mm, worst, a, b in same_bad + worst_rows:
        tag = "same" if same else "diff"
        name_a, name_b = os.path.basename(a)[:6], os.path.basename(b)[:6]
        print(f"  {tag} {name_a} | {name_b}  >1mm {n_mm:7d}  max {worst:.3f}")


if __name__ == "__main__":
    main()
