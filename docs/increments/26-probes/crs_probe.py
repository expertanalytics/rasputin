"""Throwaway probe for docs/increments/26-land-cover.md, not production code:
(1) MapBiomas cell area over the basin's latitudes; (2) how far a straight
edge in the mesh CRS bends when drawn straight in longitude/latitude; (3) the
variation of the area ratio across a triangle.

    python docs/increments/26-probes/crs_probe.py
"""

import numpy as np
from pyproj import CRS, Geod, Transformer

S = 0.00026949458523585647
g = Geod(ellps="WGS84")
for lat in (-7.0, -10.0, -15.0, -18.0, -21.0):
    lons = [-45.0, -45.0 + S, -45.0 + S, -45.0]
    lats = [lat, lat, lat - S, lat - S]
    area, _ = g.polygon_area_perimeter(lons, lats)
    print(f"cell area at {lat:6.1f}: {abs(area):7.1f} m2")

# The level-3 runs' CRS: TM on SIRGAS 2000, lon_0 -42, k 0.997548.
tm = CRS.from_proj4(
    "+proj=tmerc +lat_0=0 +lon_0=-42 +k=0.997548 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
)
for name, crs in (("UTM 23S", CRS.from_epsg(31983)), ("basin TM", tm)):
    inv = Transformer.from_crs(crs, CRS.from_epsg(4326), always_xy=True)
    fwd = Transformer.from_crs(CRS.from_epsg(4326), crs, always_xy=True)
    worst = {}
    for lon0, lat0 in ((-45.0, -10.0), (-47.5, -20.0), (-38.5, -9.0), (-44.0, -15.0)):
        x0, y0 = fwd.transform(lon0, lat0)
        for L in (100.0, 1_000.0, 5_000.0, 20_000.0):
            for ang in np.linspace(0, np.pi, 7):
                xa, ya = x0, y0
                xb, yb = x0 + L * np.cos(ang), y0 + L * np.sin(ang)
                t = np.linspace(0, 1, 101)
                lon, lat = inv.transform(xa + t * (xb - xa), ya + t * (yb - ya))
                # chord in lon/lat between the ends, mapped back to the mesh CRS
                lc = lon[0] + t * (lon[-1] - lon[0])
                la = lat[0] + t * (lat[-1] - lat[0])
                xc, yc = fwd.transform(lc, la)
                d = np.hypot(xc - (xa + t * (xb - xa)), yc - (ya + t * (yb - ya))).max()
                worst[L] = max(worst.get(L, 0.0), d)
    print(name, {int(k): f"{v:.2e} m" for k, v in worst.items()})

    # area ratio (projected area / ellipsoidal area) over 1 km and 20 km squares
    for lon0, lat0 in ((-38.5, -9.0), (-47.5, -20.0)):
        x0, y0 = fwd.transform(lon0, lat0)
        ratios = []
        for dx, dy in ((0, 0), (1_000, 0), (0, 1_000), (20_000, 0), (0, 20_000)):
            xs = np.array([0, 100, 100, 0]) + x0 + dx
            ys = np.array([0, 0, 100, 100]) + y0 + dy
            lo, la = inv.transform(xs, ys)
            a, _ = g.polygon_area_perimeter(lo, la)
            ratios.append(1e4 / abs(a))
        r = np.array(ratios)
        near, far = abs(r[1:3] - r[0]).max(), abs(r[3:5] - r[0]).max()
        print(
            f"  {name} at ({lon0},{lat0}): k^2 = {r[0]:.6f}; "
            f"change over 1 km {near:.2e}, over 20 km {far:.2e}"
        )
