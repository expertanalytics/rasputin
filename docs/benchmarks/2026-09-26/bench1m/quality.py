"""quality.py MESH.vtk -> one JSON line of plan-view quality measures.

Reads rasputin's legacy ASCII POLYDATA (POINTS, LINES = constraint edges,
POLYGONS). Same code for every run.
- min angle per triangle in x/y: median, % < 1 deg, count < 0.1 deg, worst
- vertex degree = incident triangles: max, count >= 12, count >= 20
- constrained Delaunay: for every interior edge that is not a LINES edge, the
  apex of one side tested against the other side's circumcircle. Coordinates are
  taken relative to the edge's first endpoint (exact for lattice nodes); a
  float determinant decides unless within 4x Shewchuk's incircle error bound,
  then Fractions on the exact doubles decide. 'Inside' (strict) is a violation.
  Counted per edge (the tests count per side, i.e. up to 2x this).
"""
import sys, json
from fractions import Fraction
import numpy as np

def read(path):
    with open(path) as fh:
        lines = fh.read().split("\n")
    i = 0
    pts = tri = edges = None
    while i < len(lines):
        s = lines[i].split()
        if s and s[0] == "POINTS":
            n = int(s[1]); pts = np.loadtxt(lines[i+1:i+1+n], ndmin=2)[:, :2]; i += n
        elif s and s[0] == "LINES":
            n = int(s[1]); edges = np.loadtxt(lines[i+1:i+1+n], dtype=np.int64, ndmin=2)[:, 1:3]; i += n
        elif s and s[0] == "POLYGONS":
            n = int(s[1]); tri = np.loadtxt(lines[i+1:i+1+n], dtype=np.int64, ndmin=2)[:, 1:4]; i += n
        i += 1
    if edges is None: edges = np.zeros((0, 2), np.int64)
    return pts, tri, edges

def incircle_exact(a, b, c, d):
    F = [ (Fraction(p[0]) - Fraction(d[0]), Fraction(p[1]) - Fraction(d[1])) for p in (a, b, c)]
    (ax, ay), (bx, by), (cx, cy) = F
    al, bl, cl = ax*ax+ay*ay, bx*bx+by*by, cx*cx+cy*cy
    return ax*(by*cl-bl*cy) - ay*(bx*cl-bl*cx) + al*(bx*cy-by*cx)

def main(path):
    P, T, E = read(path)
    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    def ang(p, q, r):  # angle at p
        u, v = q - p, r - p
        return np.degrees(np.arctan2(np.abs(u[:, 0]*v[:, 1]-u[:, 1]*v[:, 0]), (u*v).sum(1)))
    mn = np.minimum(np.minimum(ang(a, b, c), ang(b, c, a)), ang(c, a, b))
    deg = np.bincount(T.ravel(), minlength=len(P))
    # orientation: make every triangle CCW for the incircle sign convention
    o = (b[:, 0]-a[:, 0])*(c[:, 1]-a[:, 1]) - (b[:, 1]-a[:, 1])*(c[:, 0]-a[:, 0])
    T2 = T.copy(); T2[o < 0] = T2[o < 0][:, [0, 2, 1]]
    nt = len(T2)
    # half-edges (u, v, tri, apex)
    u = T2.ravel(); v = T2[:, [1, 2, 0]].ravel(); w = T2[:, [2, 0, 1]].ravel()
    t = np.repeat(np.arange(nt), 3)
    key = np.minimum(u, v) * (len(P) + 1) + np.maximum(u, v)
    order = np.argsort(key, kind="stable"); ks = key[order]
    same = np.nonzero(ks[1:] == ks[:-1])[0]
    h1, h2 = order[same], order[same + 1]
    if len(E):
        ck = set((np.minimum(E[:, 0], E[:, 1]) * (len(P) + 1) + np.maximum(E[:, 0], E[:, 1])).tolist())
        keep = np.array([k not in ck for k in key[h1].tolist()], bool)
        h1, h2 = h1[keep], h2[keep]
    tri = T2[t[h1]]; apex = w[h2]
    A, B, C, D = P[tri[:, 0]], P[tri[:, 1]], P[tri[:, 2]], P[apex]
    ad, bd, cd = A - D, B - D, C - D
    al = (ad**2).sum(1); bl = (bd**2).sum(1); cl = (cd**2).sum(1)
    x1 = bd[:, 0]*cd[:, 1] - bd[:, 1]*cd[:, 0]
    x2 = cd[:, 0]*ad[:, 1] - cd[:, 1]*ad[:, 0]
    x3 = ad[:, 0]*bd[:, 1] - ad[:, 1]*bd[:, 0]
    det = al*x1 + bl*x2 + cl*x3
    perm = al*np.abs(x1) + bl*np.abs(x2) + cl*np.abs(x3)
    eps = np.finfo(float).eps / 2
    bound = 4 * (10 + 96*eps) * eps * perm
    inside = det > bound
    amb = np.nonzero(np.abs(det) <= bound)[0]
    for k in amb.tolist():
        if incircle_exact(A[k], B[k], C[k], D[k]) > 0: inside[k] = True
    r = dict(triangles=int(nt), vertices=int(len(P)), constraint_edges=int(len(E)),
             minang_median=float(np.median(mn)), pct_lt1=float(100*np.mean(mn < 1)),
             n_lt1=int((mn < 1).sum()), n_lt01=int((mn < 0.1).sum()), worst=float(mn.min()),
             deg_max=int(deg.max()), deg_ge12=int((deg >= 12).sum()), deg_ge20=int((deg >= 20).sum()),
             cdt_edges_checked=int(len(h1)), cdt_ambiguous=int(len(amb)), cdt_violations=int(inside.sum()))
    print(json.dumps(r))

main(sys.argv[1])
