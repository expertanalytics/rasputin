"""M8: the land cover moved, inside the domain only; and the longest thick pieces' widths of the polygons they belong to."""
import sys, numpy as np, shapely
cover, dom = sys.argv[1], sys.argv[2]
a, b = np.load(cover, allow_pickle=True); D = shapely.from_wkb(np.load(dom, allow_pickle=True)[0])
A, B = shapely.from_wkb(list(a)), shapely.from_wkb(list(b))
X = shapely.symmetric_difference(A, B); d = shapely.area(X); idx = np.nonzero(d > 1e-9)[0]
shapely.prepare(D)
inside = shapely.area(shapely.intersection(X[idx], D)).sum()
print(f"{cover.split('/')[-1]}: moved {d.sum():.1f} m2 in all, {inside:.1f} m2 inside the domain (each m2 counted in both polygons); domain {D.area/1e6:.1f} km2")
parts, owner = shapely.get_parts(X[idx], return_index=True); owner = idx[owner]
w = np.array([2 * shapely.maximum_inscribed_circle(p, 1e-3).length for p in parts])
L = shapely.length(parts) / 2; top = np.argsort(-L * (w > 0.25))[:5]
for i in top:
    o = owner[i]; c = parts[i].centroid
    print(f"  piece length {L[i]:.0f} m thick {w[i]:.2f} m at ({c.x:.0f}, {c.y:.0f}); inside domain {D.intersects(parts[i])}; owner polygon width {2*shapely.maximum_inscribed_circle(A[o], 1e-2).length:.1f} m, area {A[o].area:.0f} m2, distance to its own other border: min width at piece? ")
