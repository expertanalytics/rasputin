"""M8: where and how thick the land cover the repair moved is, from m8drive's saved cover,
and how many slivers (smallest angle < 1 deg) of a mesh lie within 2 m of it."""
import sys, re, numpy as np, shapely
cover, vtk = sys.argv[1], sys.argv[2:]
a, b = np.load(cover, allow_pickle=True)
A, B = shapely.from_wkb(list(a)), shapely.from_wkb(list(b))
X = shapely.symmetric_difference(A, B); d = shapely.area(X)
parts = shapely.get_parts(X[d > 1e-9]); parts = parts[shapely.area(parts) > 1e-6]
w = np.array([2 * shapely.maximum_inscribed_circle(p, 1e-3).length for p in parts])
em = [(i, round(float(shapely.area(A[i])), 2), round(2 * shapely.maximum_inscribed_circle(A[i], 1e-3).length, 3)) for i, y in enumerate(B) if y.is_empty]
print(f"{cover.split('/')[-1]}: polygons changed {(d > 1e-9).sum()} of {len(A)}; area {d.sum():.2f} m2 (each m2 counted in both polygons); pieces {len(parts)}; "
      f"piece thickness max {w.max() if len(w) else 0:.3f} m, p99 {np.percentile(w, 99) if len(w) else 0:.3f}, median {np.median(w) if len(w) else 0:.3f}; "
      f"pieces thicker than 0.25 m {(w > 0.25).sum()} holding {shapely.area(parts[w > 0.25]).sum() if len(w) else 0:.1f} m2; emptied (index, area m2, width m) {em}")
U = shapely.union_all(parts) if len(parts) else shapely.Polygon(); shapely.prepare(U)
for f in vtk:
    bb = open(f, "rb").read()
    m = re.search(rb"\nPOINTS (\d+) (\w+)\n", bb); n = int(m.group(1))
    p = np.frombuffer(bb, dtype={"double": ">f8", "float": ">f4"}[m.group(2).decode()], count=3*n, offset=m.end()).reshape(n, 3)[:, :2].astype(float)
    m2 = re.search(rb"\nPOLYGONS (\d+) (\d+)\n", bb); k = int(m2.group(1))
    t = np.frombuffer(bb, dtype=">i4", count=4*k, offset=m2.end()).reshape(k, 4)[:, 1:]
    q, r, s = p[t[:, 0]], p[t[:, 1]], p[t[:, 2]]
    def ang(u, v, z):
        e1, e2 = v-u, z-u
        return np.degrees(np.arctan2(np.abs(e1[:, 0]*e2[:, 1]-e1[:, 1]*e2[:, 0]), (e1*e2).sum(1)))
    mn = np.minimum(np.minimum(ang(q, r, s), ang(r, s, q)), ang(s, q, r)); sl = mn < 1
    tri = shapely.polygons(np.stack([q[sl], r[sl], s[sl]], 1))
    near = shapely.dwithin(tri, U, 2.0)
    print(f"  {f.split('/')[-1]}: slivers {sl.sum()}, of them within 2 m of land cover this repair moved {near.sum()}")
