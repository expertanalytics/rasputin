"""Throwaway prototype (increment 32 design probe): APSC per shared border of a
land-cover coverage, ends fixed, with a band against the source border, a
crossing test and a swept-region vertex test across all borders.

    python border_apsc.py BAND OUT.geojson [unanchored]

"unanchored" uses increment 22's band check (every source vertex within BAND
of the new edges, the new point within BAND of any source segment) instead
of the anchored one the design adopts.
"""

UNANCHORED = False

import heapq
import json
import math
import sys
import time
from collections import defaultdict

import numpy as np
import shapely
from pyproj import Transformer
from shapely.geometry import MultiPolygon, Polygon, mapping, shape

OUTLINE = "/Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/fused_catchment_glo30.geojson"
CORINE = "/Users/skavhaug/projects/rasputin_data/germany_corine/isar_loisach_clc2018_3035.geojson"


def orient(a, b, c):
    v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (v > 0) - (v < 0)


def seg_dist(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    l2 = dx * dx + dy * dy
    t = 0.0 if l2 == 0 else max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2))
    return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)


def on_seg(p, a, b):
    return min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def meet(a, b, c, d):
    """Closed segments ab and cd share a point."""
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    if o1 != o2 and o3 != o4:
        return True
    return (o1 == 0 and on_seg(c, a, b)) or (o2 == 0 and on_seg(d, a, b)) or (
        o3 == 0 and on_seg(a, c, d)) or (o4 == 0 and on_seg(b, c, d))


def winding(poly, k):
    w = 0
    for i in range(len(poly)):
        p, q = poly[i], poly[(i + 1) % len(poly)]
        o = orient(p, q, k)
        if p[1] <= k[1] < q[1] and o > 0:
            w += 1
        elif q[1] <= k[1] < p[1] and o < 0:
            w -= 1
    return w


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


class Chain:
    def __init__(self, pts, closed, fixed_all, fixed_first):
        self.fine = pts  # list of tuples; closed: first not repeated
        self.closed = closed
        n = len(pts)
        self.p = list(pts)
        self.alive = [True] * n
        self.nxt = [(i + 1) % n if closed else (i + 1 if i + 1 < n else -1) for i in range(n)]
        self.prv = [(i - 1) % n if closed else (i - 1) for i in range(n)]
        self.start = list(range(n))
        self.anc = [(i, 0.0) for i in range(n)]
        self.count = [1] * n
        self.ver = [0] * n
        self.fixed = [fixed_all] * n
        if not closed and n:
            self.fixed[0] = self.fixed[-1] = True
        if closed and fixed_first:
            self.fixed[0] = True
        self.size = n

    def at(self, k):
        return self.fine[k % len(self.fine)] if self.closed else self.fine[k]


def deviation_of(ch, s, n, A, E, D):
    at = ch.at
    suffix = [0.0] * (n + 1)
    for k in range(n - 1, -1, -1):
        suffix[k] = max(suffix[k + 1], seg_dist(at(s + k), E, D))
    best_split, best = 0, suffix[0]
    prefix = 0.0
    for m in range(1, n + 1):
        prefix = max(prefix, seg_dist(at(s + m - 1), A, E))
        v = max(prefix, suffix[m])
        if v < best:
            best, best_split = v, m
    to_fine = min(seg_dist(E, at(s + k), at(s + k + 1)) for k in range(max(n, 1)))
    return max(best, to_fine), best_split


def rel(ch, k, s):
    return (k - s) % len(ch.fine) if ch.closed else k - s


def deviation_anchored(ch, s, n, A, E, D, ancA, ancD):
    at = ch.at
    suffix = [0.0] * (n + 1)
    for k in range(n - 1, -1, -1):
        suffix[k] = max(suffix[k + 1], seg_dist(at(s + k), E, D))
    ra = (rel(ch, ancA[0], s), ancA[1])
    if ch.closed and ra[0] > n:  # A's anchor just before s
        ra = (ra[0] - len(ch.fine), ra[1])
    rd = (rel(ch, ancD[0], s), ancD[1])
    best = None
    prefix = 0.0
    for m in range(1, n + 1):
        prefix = max(prefix, seg_dist(at(s + m - 1), A, E))
        P, Q = at(s + m - 1), at(s + m)
        dx, dy = Q[0] - P[0], Q[1] - P[1]
        l2 = dx * dx + dy * dy
        t = 0.0 if l2 == 0 else max(0.0, min(1.0, ((E[0] - P[0]) * dx + (E[1] - P[1]) * dy) / l2))
        anc = (m - 1, t)
        if anc < ra or anc > rd:
            continue
        v = max(prefix, suffix[m], math.hypot(P[0] + t * dx - E[0], P[1] + t * dy - E[1]))
        if best is None or v < best[0]:
            best = (v, m, (s + m - 1, t))
    return best


def simplify(chains, band):
    cell = band
    grid = defaultdict(list)

    def cells(*pts):
        x0 = math.floor(min(p[0] for p in pts) / cell); x1 = math.floor(max(p[0] for p in pts) / cell)
        y0 = math.floor(min(p[1] for p in pts) / cell); y1 = math.floor(max(p[1] for p in pts) / cell)
        for ix in range(x0, x1 + 1):
            for iy in range(y0, y1 + 1):
                yield (ix, iy)

    def insert(ci, u):
        ch = chains[ci]
        v = ch.nxt[u]
        for c in cells(ch.p[u], ch.p[v]):
            grid[c].append((ci, u))

    for ci, ch in enumerate(chains):
        for u in range(len(ch.p)):
            if ch.nxt[u] != -1:
                insert(ci, u)

    def candidate(ci, b):
        ch = chains[ci]
        a = ch.prv[b]; c = ch.nxt[b]
        if a == -1 or c == -1:
            return None
        d = ch.nxt[c]
        if d == -1 or ch.fixed[b] or ch.fixed[c] or len({a, b, c, d}) < 4:
            return None
        A = ch.p[a]
        vb = (ch.p[b][0] - A[0], ch.p[b][1] - A[1]); vc = (ch.p[c][0] - A[0], ch.p[c][1] - A[1])
        vd = (ch.p[d][0] - A[0], ch.p[d][1] - A[1])
        area = cross(vb, vc) + cross(vc, vd)
        s = ch.start[a]; n = ch.count[a] + ch.count[b] + ch.count[c]
        best = None
        cands = []
        bd = cross(vb, vd)
        if bd != 0:
            f = area / bd; cands.append((vb[0] * f, vb[1] * f))
        cd = cross(vc, vd)
        if cd != 0:
            t = 1.0 - area / cd; cands.append((vc[0] + (vd[0] - vc[0]) * t, vc[1] + (vd[1] - vc[1]) * t))
        for e in cands:
            E = (A[0] + e[0], A[1] + e[1])
            if UNANCHORED:  # increment 22's check: vertices to edges, E to any segment
                dev0, split0 = deviation_of(ch, s, n, A, E, ch.p[d])
                r = (dev0, split0, (s + max(split0, 1) - 1, 0.0))
            else:
                r = deviation_anchored(ch, s, n, A, E, ch.p[d], ch.anc[a], ch.anc[d])
            if r is None:
                continue
            dev, split, anc = r
            if best is None or dev < best[0]:
                best = (dev, E, split, anc)
        return best

    heap = []
    pending = {}

    def evaluate(ci, b):
        ch = chains[ci]
        ch.ver[b] += 1
        cand = candidate(ci, b)
        pending[(ci, b)] = cand
        if cand and cand[0] <= band:
            heapq.heappush(heap, (cand[0], ci, b, ch.ver[b]))

    for ci, ch in enumerate(chains):
        for b in range(len(ch.p)):
            evaluate(ci, b)
    rej = defaultdict(int)
    done = 0
    while heap:
        dev, ci, b, ver = heapq.heappop(heap)
        ch = chains[ci]
        if not ch.alive[b] or ch.ver[b] != ver:
            continue
        minsize = 4 if ch.closed else 4
        if ch.size < minsize:
            continue
        a = ch.prv[b]; c = ch.nxt[b]; d = ch.nxt[c]
        _, E, split, anc = pending[(ci, b)]
        A, B, C, D = ch.p[a], ch.p[b], ch.p[c], ch.p[d]
        # crossing: AE and ED against all live edges near
        bad = orient(A, E, D) == 0 and not on_seg(E, A, D)
        seen = set()
        loop = [A, B, C, D, E]
        for cc in cells(A, B, C, D, E):
            if bad:
                break
            for (cj, u) in grid.get(cc, ()):
                if (cj, u) in seen:
                    continue
                seen.add((cj, u))
                o = chains[cj]
                if not o.alive[u] or o.nxt[u] == -1:
                    continue
                v = o.nxt[u]
                if cj == ci and u in (a, b, c):
                    continue
                P, Q = o.p[u], o.p[v]
                if (P[0], P[1]) != (o.p[u][0], o.p[u][1]):
                    continue
                # stale entry: edge u->v must still be the one inserted; we accept (grid holds both)
                for (S0, S1, end) in ((A, E, A), (E, D, D)):
                    if meet(S0, S1, P, Q):
                        shared = (P == end or Q == end)
                        if not shared:
                            bad = True; break
                        # touching at the shared end only: the other end must not lie on the segment
                        other = Q if P == end else P
                        far = S1 if end == S0 else S0
                        if orient(S0, S1, other) == 0 and on_seg(other, S0, S1):
                            bad = True; break
                        if orient(P, Q, far) == 0 and on_seg(far, P, Q):
                            bad = True; break
                if bad:
                    break
                for K in (P, Q):
                    if K in (A, D) or (cj == ci and (u in (a, b, c) or v in (b, c, d))):
                        continue
                    if winding(loop, K) != 0:
                        bad = True; break
        if bad:
            rej["topology"] += 1
            continue
        n = ch.count[a] + ch.count[b] + ch.count[c]
        e = len(ch.p)
        ch.p.append(E); ch.alive.append(True); ch.ver.append(0); ch.fixed.append(False)
        ch.nxt.append(d); ch.prv.append(a); ch.start.append(ch.start[a] + split); ch.count.append(n - split); ch.anc.append(anc)
        ch.alive[b] = ch.alive[c] = False
        ch.nxt[a] = e; ch.prv[d] = e
        ch.count[a] = split
        insert(ci, a); insert(ci, e)
        ch.size -= 1
        done += 1
        for u in (ch.prv[a], a, e, d):
            if u != -1:
                evaluate(ci, u)
    return done, rej


def walk(ch):
    if ch.closed:
        u = next(i for i in range(len(ch.p)) if ch.alive[i] and ch.fixed[i]) if any(
            ch.alive[i] and ch.fixed[i] for i in range(len(ch.p))) else next(i for i in range(len(ch.p)) if ch.alive[i])
        out = [ch.p[u]]
        v = ch.nxt[u]
        while v != u:
            out.append(ch.p[v]); v = ch.nxt[v]
        return out
    u = 0; out = [ch.p[0]]
    while ch.nxt[u] != -1:
        u = ch.nxt[u]; out.append(ch.p[u])
    return out


def main():
    global UNANCHORED
    band = float(sys.argv[1]); outpath = sys.argv[2]
    UNANCHORED = len(sys.argv) > 3 and sys.argv[3] == "unanchored"
    t = time.perf_counter()
    dom = shape(json.load(open(OUTLINE))["features"][0]["geometry"])
    to = Transformer.from_crs(3035, 25832, always_xy=True)
    feats = json.load(open(CORINE))["features"]
    geoms = [shapely.transform(shape(f["geometry"]), lambda xy: np.column_stack(to.transform(xy[:, 0], xy[:, 1]))) for f in feats]
    codes = [int(f["properties"]["Code_18"]) for f in feats]
    clipped = [shapely.intersection(g, dom) for g in geoms]
    keep = [i for i, g in enumerate(clipped) if not g.is_empty and g.area > 0]
    polys = np.array([clipped[i] for i in keep], dtype=object); codes = [codes[i] for i in keep]
    polys = shapely.coverage_clean(polys, snapping_distance=1.0, gap_width=1.0, merge_strategy="min_area")
    groups = {}
    for i, c in enumerate(codes):
        groups.setdefault(c, []).append(i)
    classes = list(groups)
    merged = [shapely.coverage_union_all(polys[groups[c]]) for c in classes]
    print(f"prep {time.perf_counter() - t:.1f} s; valid coverage: {shapely.coverage_is_valid(np.array(merged, dtype=object))}")
    parts = [(ci, q) for ci, m in enumerate(merged) for q in shapely.get_parts(m) if not q.is_empty]
    small = [(ci, q) for ci, q in parts if q.area < 25e4]
    ob = dom.boundary
    print(f"pieces {len(parts)}; under 25 ha {len(small)} ({sum(q.area for _, q in small)/1e6:.2f} km2), of which touch outline "
          f"{sum(q.distance(ob) < 1e-6 for _, q in small)}; vertices of small pieces {sum(shapely.get_num_coordinates(q) for _, q in small)}")
    # rings
    rings = []  # (class, part, ringidx, coords list closed-not-repeated)
    for pi, (ci, q) in enumerate(parts):
        for ri, r in enumerate(shapely.get_rings(q)):
            xy = [(float(p[0]), float(p[1])) for p in shapely.get_coordinates(r)[:-1]]
            rings.append((ci, pi, ri, xy))
    use = defaultdict(int); nbr = defaultdict(set)
    for *_, xy in rings:
        n = len(xy)
        for i in range(n):
            p, q = xy[i], xy[(i + 1) % n]
            key = (p, q) if p < q else (q, p)
            use[key] += 1
    for (p, q) in use:
        nbr[p].add(q); nbr[q].add(p)
    node = {v for v, s in nbr.items() if len(s) != 2}
    print(f"vertices {sum(len(r[3]) for r in rings)}, distinct edges {len(use)}, nodes {len(node)}")
    chains = []; chain_key = {}; ring_refs = []
    for ci, pi, ri, xy in rings:
        n = len(xy)
        idx = [i for i in range(n) if xy[i] in node]
        refs = []
        if not idx:
            k0 = min(range(n), key=lambda i: xy[i])
            key = ("ring", xy[k0])
            if key not in chain_key:
                pts = xy[k0:] + xy[:k0]
                e = (pts[0], pts[1]); e = e if e[0] < e[1] else (e[1], e[0])
                chain_key[key] = len(chains)
                chains.append(Chain(pts, True, use[e] == 1, False))
                chains[-1].anchor = pts
            j = chain_key[key]
            pts = chains[j].fine
            # direction: same as ring if the successor of pts[0] in xy equals pts[1]
            fw = xy[(k0 + 1) % n] == pts[1]
            refs.append((j, fw))
        else:
            for t_, s in enumerate(idx):
                e_ = idx[(t_ + 1) % len(idx)]
                seg = [xy[(s + k) % n] for k in range(((e_ - s) % n or n) + 1)]
                kf = (seg[0], seg[1]); kb = (seg[-1], seg[-2])
                if kf in chain_key:
                    refs.append((chain_key[kf], True))
                elif kb in chain_key:
                    refs.append((chain_key[kb], False))
                else:
                    e = (seg[0], seg[1]); e = e if e[0] < e[1] else (e[1], e[0])
                    j = len(chains)
                    chain_key[kf] = j
                    if seg[0] == seg[-1]:
                        chains.append(Chain(seg[:-1], True, use[e] == 1, True))
                    else:
                        chains.append(Chain(seg, False, use[e] == 1, False))
                    refs.append((j, True))
        ring_refs.append(refs)
    fixed = sum(c.fixed[0] and all(c.fixed) for c in chains)
    print(f"chains {len(chains)} ({fixed} fixed, outline); closed {sum(c.closed for c in chains)}")
    t = time.perf_counter()
    done, rej = simplify(chains, band) if band > 0 else (0, {})
    print(f"simplify {band} m: {time.perf_counter() - t:.1f} s, collapses {done}, rejected {dict(rej)}")
    out = [walk(c) for c in chains]
    # rebuild
    polys_new = defaultdict(lambda: defaultdict(list))
    for (ci, pi, ri, xy), refs in zip(rings, ring_refs):
        coords = []
        for j, fw in refs:
            seq = out[j] if fw else out[j][::-1]
            if chains[j].closed:
                # rotate to start at the ring's own start: fine, a closed chain stands alone
                coords = seq + [seq[0]] if fw else [seq[-1]] + seq
                if not fw:
                    coords = seq[::-1][::-1]
                    coords = list(reversed(out[j])); coords = coords + [coords[0]]
                else:
                    coords = out[j] + [out[j][0]]
                break
            coords += seq if not coords else seq[1:]
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        polys_new[ci][pi].append(coords)
    newgeoms = []
    for ci in range(len(classes)):
        ps = [Polygon(r[0], r[1:]) for pi, r in sorted(polys_new[ci].items())]
        newgeoms.append(MultiPolygon(ps))
    arr = np.array(newgeoms, dtype=object)
    print(f"valid polygons {all(shapely.is_valid(arr))}; valid coverage {shapely.coverage_is_valid(arr)}")
    va = int(shapely.get_num_coordinates(np.array(merged, dtype=object)).sum()); vb = int(shapely.get_num_coordinates(arr).sum())
    print(f"polygon vertices {va} -> {vb}")
    worst = 0.0
    for ci in range(len(classes)):
        da = newgeoms[ci].area - merged[ci].area
        worst = max(worst, abs(da) / merged[ci].area)
    print(f"largest class area change {100 * worst:.2e} %")
    def directed(a, b):
        pts = shapely.points(shapely.get_coordinates(shapely.segmentize(a, 0.5)))
        return float(shapely.distance(pts, b).max())
    fwd, bwd = [], []
    for c, o in zip(chains, out):
        if len(o) < 2 or len(c.fine) < 2 or all(c.fixed):
            continue
        L1 = shapely.LineString(c.fine + ([c.fine[0]] if c.closed else []))
        L2 = shapely.LineString(o + ([o[0]] if c.closed else []))
        fwd.append(directed(L1, L2)); bwd.append(directed(L2, L1))
    fwd, bwd = np.array(fwd), np.array(bwd)
    print(f"source->simplified: max {fwd.max():.1f} m, over band {int((fwd > band + 1e-6).sum())} borders; "
          f"simplified->source: max {bwd.max():.1f} m, over band {int((bwd > band + 1e-6).sum())}; median Hausdorff {np.median(np.maximum(fwd,bwd)):.1f} m")
    changed = sum(shapely.difference(newgeoms[ci], merged[ci]).area for ci in range(len(classes)))
    print(f"ground with a different class: {changed/1e6:.2f} km2 ({100*changed/dom.area:.2f} %)")
    fc = {"type": "FeatureCollection", "crs": {"type": "name", "properties": {"name": "EPSG:25832"}},
          "features": [{"type": "Feature", "properties": {"Code_18": str(classes[ci])}, "geometry": mapping(g)} for ci, g in enumerate(newgeoms)]}
    json.dump(fc, open(outpath, "w"))


if __name__ == "__main__":
    main()
