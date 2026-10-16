"""M8: where the pieces 1 m moved (thicker than 25 cm) lie: their length, and on Lagan their
distance in EPSG:3035 to M2's seam N 3 772 265.53 and easting pair."""
import sys, numpy as np, shapely, pyproj
cover, crs = sys.argv[1], sys.argv[2]
a, b = np.load(cover, allow_pickle=True)
A, B = shapely.from_wkb(list(a)), shapely.from_wkb(list(b))
X = shapely.symmetric_difference(A, B); d = shapely.area(X)
parts = shapely.get_parts(X[d > 1e-9]); parts = parts[shapely.area(parts) > 1e-6]
w = np.array([2 * shapely.maximum_inscribed_circle(p, 1e-3).length for p in parts])
T = parts[w > 0.25]; L = shapely.length(T) / 2
c = shapely.get_coordinates(shapely.centroid(T))
x, y = pyproj.Transformer.from_crs(crs, "EPSG:3035", always_xy=True).transform(c[:, 0], c[:, 1])
print(f"{cover.split('/')[-1]}: thick pieces {len(T)}, length sum {L.sum():.0f} m, median {np.median(L):.1f} m, max {L.max():.0f} m; mean thickness area/length {shapely.area(T).sum()/L.sum():.3f} m")
seam = np.abs(y - 3772265.53) < 15; east = (np.abs(x - 4537573.71) < 18) | (np.abs(x - 4537579.53) < 18)
print(f"  within 15 m of N 3772265.53: {seam.sum()} ({shapely.area(T[seam]).sum():.0f} m2); within 18 m of the easting pair: {east.sum()}")
# straight-line seams in EPSG:3035: piece centroids sharing a rounded northing or easting
for name, v in (("northing", y), ("easting", x)):
    r = np.round(v / 5) * 5; u, n = np.unique(r, return_counts=True); top = np.argsort(-n)[:4]
    print(f"  most common {name} (5 m bins): {[(float(u[i]), int(n[i])) for i in top]}")
