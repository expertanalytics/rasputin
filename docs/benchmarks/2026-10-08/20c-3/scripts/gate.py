import re, sys, json, numpy as np, shapely
from shapely.geometry import shape
dom = sys.argv[1]; g = json.load(open(dom)); gg = g["features"][0]["geometry"] if "features" in g else g
outline = shape(gg).boundary
for f in sys.argv[2:]:
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
    s = mn < 1
    side = np.minimum(np.minimum(np.hypot(*(a-bb).T), np.hypot(*(bb-c).T)), np.hypot(*(c-a).T))
    cen = (a+bb+c)[s]/3
    near = int((shapely.distance(shapely.points(cen), outline) <= 20).sum())
    print(f"{f.split('/')[-1]}: tri {k} | <1deg {s.sum()} ({100*s.mean():.4f} %) | worst {mn.min():.6f} | <1deg with side<10cm {(s & (side<0.1)).sum()} | <1deg centre within 20 m of outline {near}")
