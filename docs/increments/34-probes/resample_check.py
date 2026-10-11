"""Increment 34 design probe: the slope on the resampled path against the
source's own (README item 1c; design review round 1, B2).

With --out-crs, refine meshes a square target grid whose spacing is the
source's north-south spacing rounded to whole metres
(src_python/tin_engine/target_grid.py, default_spacing). For a geographic
source at latitude phi the east-west spacing is cos(phi) times finer, so the
target grid is coarser east-west than the source. This probe takes a steep
window of Copernicus GLO-30 (the Wetterstein and Karwendel, N47 E011), and
compares the shares of steep nodes by Horn's slope (section 3's rule) on the
source grid (its own metre spacings, per row) and on a square grid in
EPSG:25832 at the rounded north-south spacing, bilinear from the source as the
resampled path samples it.
"""

import numpy as np
import pyproj
import tifffile
from slope_stats import horn

GLO30 = "/Users/skavhaug/projects/rasputin_data/germany_glo30"
SRC = f"{GLO30}/Copernicus_DSM_COG_10_N47_00_E011_00_DEM.tif"
LON0, LON1, LAT0, LAT1 = 11.0, 11.5, 47.3, 47.55


def main():
    with tifffile.TiffFile(SRC) as t:
        page = t.pages[0]
        z = page.asarray().astype(np.float64)
        tie = page.tags["ModelTiepointTag"].value
        scale = page.tags["ModelPixelScaleTag"].value
    lon_min, lat_max, dlon, dlat = tie[3], tie[4], scale[0], scale[1]
    c0, c1 = int((LON0 - lon_min) / dlon), int((LON1 - lon_min) / dlon)
    r0, r1 = int((lat_max - LAT1) / dlat), int((lat_max - LAT0) / dlat)
    win = z[r0:r1, c0:c1]
    geod = pyproj.Geod(ellps="WGS84")
    lat_mid = lat_max - (r0 + r1) / 2 * dlat
    _, _, dy = geod.inv(11.25, lat_mid - dlat / 2, 11.25, lat_mid + dlat / 2)
    _, _, dx = geod.inv(11.25 - dlon / 2, lat_mid, 11.25 + dlon / 2, lat_mid)
    s_src = horn(win, dx, dy)
    h = max(1, round(dy))
    fwd = pyproj.Transformer.from_crs(4326, 25832, always_xy=True)
    inv = pyproj.Transformer.from_crs(25832, 4326, always_xy=True)
    xs, ys = fwd.transform([LON0 + 0.02, LON1 - 0.02], [LAT0 + 0.02, LAT1 - 0.02])
    gx = np.arange(np.ceil(xs[0] / h) * h, xs[1], h)
    gy = np.arange(np.floor(ys[1] / h) * h, ys[0], -h)
    X, Y = np.meshgrid(gx, gy)  # noqa: N806
    lon, lat = inv.transform(X, Y)
    fc, fr = (lon - lon_min) / dlon, (lat_max - lat) / dlat
    i0, j0 = np.floor(fr).astype(int), np.floor(fc).astype(int)
    ty, tx = fr - i0, fc - j0
    zt = (z[i0, j0] * (1 - tx) * (1 - ty) + z[i0, j0 + 1] * tx * (1 - ty)
          + z[i0 + 1, j0] * (1 - tx) * ty + z[i0 + 1, j0 + 1] * tx * ty)  # fmt: skip
    s_tgt = horn(zt, float(h), float(h))
    print(
        f"source window {win.shape[0]} x {win.shape[1]} nodes, {dx:.2f} m (east-west) "
        f"by {dy:.2f} m; target {zt.shape[0]} x {zt.shape[1]} nodes at {h} m square"
    )
    for a in (15, 25, 30, 35, 40, 45):
        print(
            f"  at {a} deg or more: source {100 * np.mean(s_src >= a):5.1f}%  "
            f"target {100 * np.mean(s_tgt >= a):5.1f}%"
        )
    print(
        f"  median slope: source {np.median(s_src):.1f} deg, target {np.median(s_tgt):.1f} deg; "
        f"p99 source {np.percentile(s_src, 99):.1f}, target {np.percentile(s_tgt, 99):.1f}"
    )


if __name__ == "__main__":
    main()
