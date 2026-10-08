"""R3 (c): triangles with smallest angle under 1 deg that have a corner within 1 mm of the slit's far corner."""
import re, sys, numpy as np
X, Y = 413677.295, 6331664.081
for f in sys.argv[1:]:
    b = open(f, "rb").read()
    m = re.search(rb"\nPOINTS (\d+) (\w+)\n", b); n = int(m.group(1))
    p = np.frombuffer(b, dtype={"double": ">f8", "float": ">f4"}[m.group(2).decode()], count=3*n, offset=m.end()).reshape(n, 3)[:, :2].astype(float)
    m2 = re.search(rb"\nPOLYGONS (\d+) (\d+)\n", b); k = int(m2.group(1))
    t = np.frombuffer(b, dtype=">i4", count=4*k, offset=m2.end()).reshape(k, 4)[:, 1:]
    a, bb, c = p[t[:, 0]], p[t[:, 1]], p[t[:, 2]]
    def ang(u, v, w):
        e1, e2 = v-u, w-u
        return np.degrees(np.arctan2(np.abs(e1[:, 0]*e2[:, 1]-e1[:, 1]*e2[:, 0]), (e1*e2).sum(1)))
    mn = np.minimum(np.minimum(ang(a, bb, c), ang(bb, c, a)), ang(c, a, bb))
    near = (np.hypot(*(p - (X, Y)).T) < 1e-3)[t].any(1)
    d = np.hypot(*(p - (X, Y)).T).min()
    print(f"{f.split('/')[-1]}: nearest vertex {d:.4f} m; triangles at the corner {near.sum()}; of them under 1 deg {(near & (mn < 1)).sum()}")
