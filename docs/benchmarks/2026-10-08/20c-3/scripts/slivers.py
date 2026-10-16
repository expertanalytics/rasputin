import re, sys, numpy as np
for f in sys.argv[1:]:
    b = open(f, "rb").read()
    m = re.search(rb"\nPOINTS (\d+) (\w+)\n", b); n = int(m.group(1))
    dt = {"double": ">f8", "float": ">f4"}[m.group(2).decode()]
    p = np.frombuffer(b, dtype=dt, count=3 * n, offset=m.end()).reshape(n, 3)[:, :2].astype(float)
    m2 = re.search(rb"\nPOLYGONS (\d+) (\d+)\n", b); k = int(m2.group(1))
    t = np.frombuffer(b, dtype=">i4", count=4 * k, offset=m2.end()).reshape(k, 4)
    assert (t[:, 0] == 3).all(); t = t[:, 1:]
    a, bb, c = p[t[:, 0]], p[t[:, 1]], p[t[:, 2]]
    def ang(u, v, w):
        e1, e2 = v - u, w - u
        return np.degrees(np.arctan2(np.abs(e1[:, 0]*e2[:, 1]-e1[:, 1]*e2[:, 0]), (e1*e2).sum(1)))
    mn = np.minimum(np.minimum(ang(a, bb, c), ang(bb, c, a)), ang(c, a, bb))
    print(f"{f.split('/')[-1]}: triangles {k}, under 1 deg {(mn<1).sum()} ({100*(mn<1).mean():.4f} %), worst {mn.min():.6f} deg")
