"""Fixtures for increment 32's tests 17-20b: two 100 m squares sharing the
border A=(100,0) .. D=(100,100) with inner vertices B, C (one collapse
possible), optionally a triangle island in the right square. Placements and
the anchored deviation re-enacted from border_collapse.hpp@633808b9; the
kernel (the worktree's _core, clearance 0) confirms the winner."""
import math, numpy as np
from tin_engine import _core

def cross(u, v): return u[0]*v[1]-u[1]*v[0]
def sub(a, b): return (a[0]-b[0], a[1]-b[1])
def add(a, b): return (a[0]+b[0], a[1]+b[1])
def mul(a, s): return (a[0]*s, a[1]*s)
def dot(u, v): return u[0]*v[0]+u[1]*v[1]
def segd(p, a, b):
    ab = sub(b, a); L = dot(ab, ab)
    t = 0.0 if L == 0 else min(max(dot(sub(p, a), ab)/L, 0.0), 1.0)
    q = add(a, mul(ab, t)); return math.hypot(p[0]-q[0], p[1]-q[1])

def placements(A, B, C, D):
    vb, vc, vd = sub(B, A), sub(C, A), sub(D, A)
    area = cross(vb, vc) + cross(vc, vd); out = []
    bd = cross(vb, vd)
    if bd != 0: out.append(("on A-B", add(A, mul(vb, area/bd))))
    cd = cross(vc, vd)
    if cd != 0: out.append(("on C-D", add(A, add(vc, mul(sub(vd, vc), 1.0 - area/cd)))))
    return out

def deviation(src, e):  # anchored(), fa=(0,0), fd=(last,1), open border
    a, d = src[0], src[-1]; span = len(src) - 1; best = None; prefix = 0.0
    suffix = [0.0]*(span+2)
    for r in range(span, 0, -1): suffix[r] = max(suffix[r+1], segd(src[r], e, d))
    for r in range(span+1):
        if r > 0: prefix = max(prefix, segd(src[r], a, e))
        if r >= span: continue
        p, q = src[r], sub(src[r+1], src[r]); L = dot(q, q)
        t = min(max(dot(sub(e, p), q)/L, 0.0), 1.0) if L > 0 else 0.0
        f = add(p, mul(q, t))
        v = max(prefix, suffix[r+1], math.hypot(e[0]-f[0], e[1]-f[1]))
        best = v if best is None else min(best, v)
    return best

def rings(B, C, island=None):
    A, D = (100.0, 0.0), (100.0, 100.0)
    left = [(0.0, 0.0), A, B, C, D, (0.0, 100.0)]
    right = [A, (200.0, 0.0), (200.0, 100.0), D, C, B]
    rs = [left, right]
    if island:
        rs.append(list(island))            # island, ccw
        rs.append(list(island)[::-1])      # the hole in the right square (as its own ring of the coverage)
    return rs

def quantities(B, C, E, island=None, frame=None):
    """The quantities of 15.2 for placement E (one collapse; checked over everything,
    a superset of what the grid looks at): |E-A|, |E-D|; every vertex other than B, C
    and A's and D's coordinates to A-E and E-D; A to E-D and D to A-E (each end is
    excluded only from its own edge); E to every edge other than A-B, B-C, C-D.
    ``frame`` = (A, D, rings) replaces the two squares."""
    A, D, rs = frame if frame else ((100.0, 0.0), (100.0, 100.0), rings(B, C, island))
    verts = {p for r in rs for p in r} - {B, C, A, D}
    edges = {tuple(sorted((r[k], r[(k+1) % len(r)]))) for r in rs for k in range(len(r))}
    edges -= {tuple(sorted(e)) for e in ((A, B), (B, C), (C, D))}
    ends = (math.dist(E, A), math.dist(E, D))
    vq = min(min(segd(v, A, E), segd(v, E, D)) for v in verts)
    a_ed, d_ae = segd(A, E, D), segd(D, A, E)
    eq = min(segd(E, *e) for e in edges)
    return dict(EA=ends[0], ED=ends[1], vertex_to_new=vq, A_to_ED=a_ed, D_to_AE=d_ae, E_to_edge=eq,
                min=min(*ends, vq, a_ed, d_ae, eq))

def run(B, C, island=None, band=50.0, frame=None):
    rs = frame[2] if frame else rings(B, C, island)
    pts = np.array([p for r in rs for p in r], dtype=np.float64)
    starts = np.cumsum([0, *map(len, rs)]).astype(np.uint64)
    return _core.simplify_borders(pts, starts, band), rs

def report(name, B, C, island=None, band=50.0, frame=None):
    A, D = frame[:2] if frame else ((100.0, 0.0), (100.0, 100.0))
    src = [A, B, C, D]
    pl = [(lab, e, deviation(src, e), quantities(B, C, e, island, frame)) for lab, e in placements(A, B, C, D)]
    win = min(pl, key=lambda x: x[2])
    out, rs = run(B, C, island, band, frame)
    c = out.counts
    newpts = {tuple(p) for p in out.points.tolist()} - {p for r in rs for p in r}
    print(f"{name}: B={B} C={C} island={island}")
    for lab, e, dev, q in pl:
        tag = "WINNER" if (lab, e) == win[:2] else "other "
        print(f"   {tag} {lab}: E=({e[0]:.4f},{e[1]:.4f}) deviation {dev:.3f}; |E-A| {q['EA']:.3f} |E-D| {q['ED']:.3f} vertex-to-new-edge {q['vertex_to_new']:.3f} A-to-ED {q['A_to_ED']:.3f} D-to-AE {q['D_to_AE']:.3f} E-to-edge {q['E_to_edge']:.3f}; min {q['min']:.3f}")
    ok = len(newpts) == 1 and math.dist(next(iter(newpts)), win[1]) < 1e-9 if newpts else False
    print(f"   kernel, clearance 0: collapses {c.collapses}, rejected crossing {c.rejected_crossing}, side {c.rejected_side}; new vertices {sorted(newpts)}; equals the re-enacted winner: {ok}")
    return pl, win
