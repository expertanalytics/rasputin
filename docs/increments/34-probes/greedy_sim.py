"""Increment 34 design probe: greedy insertion to a slope-dependent tolerance,
simulated in Python on a DTM10 window, to compare the rules of section 3
(docs/increments/34-probes/README.md, item 2).

Not rasputin's loop: each round re-triangulates the vertex set (scipy's
Delaunay, unconstrained), assigns every DEM node to one triangle
(find_simplex; a node on a shared edge goes to one of the two), and inserts the
worst node of every triangle that has not converged. rasputin inserts in the
same per-round fashion but legalises locally; the counts are compared with
rasputin's own on the uniform cases (item 3) to show how close the simulation is.

Rules, with t(s) the tolerance at slope s (degrees): F below START, N above
END, linear between (a step when START == END):
  uniform-F, uniform-N  one tolerance;
  node      every node held to t(s(n)): converge when max err(n) / t(s(n)) <= 1,
            insert the node with the largest ratio;
  tri       the triangle held to t(max s over its nodes and corners);
            insert the node of largest error (refine's choice);
  tri+halo  as tri, with s replaced by its 3x3 maximum (the cells meeting T);
  normal    the triangle held to t(slope of its own plane), the rule Ola warned of.
Reported per rule: triangles, rounds, nodes whose error exceeds t(s(n)) after
convergence (the guarantee; must be 0 for node, tri and tri+halo), and the
share of (triangle, round) evaluations whose error lay between N and F, the
evaluations where refine's lazy test would ask the policy.
"""

import argparse
import time

import numpy as np
from scipy.spatial import Delaunay
from slope_stats import horn, window


def tol_of(s, near, far, start, end):
    if end <= start:
        return np.where(s >= start, near, far)
    w = np.clip((s - start) / (end - start), 0.0, 1.0)
    return far + (near - far) * w


def dilate(s):
    p = np.pad(s, 1, mode="edge")
    return np.maximum.reduce(
        [p[i : i + s.shape[0], j : j + s.shape[1]] for i in range(3) for j in range(3)]
    )


def group_first(keys, order_vals):
    """For each key, the index with the largest order_val (ties: first)."""
    o = np.lexsort((-order_vals, keys))
    k = keys[o]
    first = np.ones(len(o), bool)
    first[1:] = k[1:] != k[:-1]
    return k[first], o[first]


def run(z, s, rule, near, far, start, end, max_rounds=200):
    rows, cols = z.shape
    rr, cc = np.mgrid[0:rows, 0:cols]
    nodes = np.column_stack([cc.ravel() * 10.0, -rr.ravel() * 10.0])
    zf, sf = z.ravel(), s.ravel()
    tn = tol_of(sf, near, far, start, end)
    sh = dilate(s).ravel()
    corner_ids = [0, cols - 1, (rows - 1) * cols, rows * cols - 1]
    verts = list(corner_ids)
    is_vert = np.zeros(len(zf), bool)
    is_vert[corner_ids] = True
    evals = between = 0
    seen = set()
    for rnd in range(1, max_rounds + 1):  # noqa: B007 (returned after the loop)
        v = np.array(verts)
        tri = Delaunay(nodes[v])
        simp = tri.find_simplex(nodes)
        xform = tri.transform[simp, :2]
        r = nodes - tri.transform[simp, 2]
        b = np.einsum("ijk,ik->ij", xform, r)
        bary = np.column_stack([b, 1 - b.sum(1)])
        corners = v[tri.simplices]  # (ntri, 3) node ids
        plane = (bary * zf[corners[simp]]).sum(1)
        err = np.abs(zf - plane)
        err[is_vert] = 0.0
        ntri = len(tri.simplices)
        maxerr = np.zeros(ntri)
        np.maximum.at(maxerr, simp, err)
        if rule == "node":
            ratio = err / np.maximum(tn, 1e-12)
            keys, idx = group_first(simp, ratio)
            worst = np.zeros(ntri)
            worst[keys] = ratio[idx]
            split = worst > 1.0
        else:
            keys, idx = group_first(simp, err)
            if rule.startswith("uniform"):
                allowed = np.full(ntri, near if rule == "uniform-N" else far)
            elif rule in ("tri", "tri+halo"):
                src = sh if rule == "tri+halo" else sf
                smax = np.zeros(ntri)
                np.maximum.at(smax, simp, src)
                smax = np.maximum(smax, src[corners].max(1))
                allowed = tol_of(smax, near, far, start, end)
            elif rule == "normal":
                p = np.column_stack([nodes[corners[:, k]] for k in range(3)]).reshape(ntri, 3, 2)
                zc = zf[corners]
                a = np.column_stack([p[:, 1] - p[:, 0], (zc[:, 1] - zc[:, 0])[:, None]])
                bb = np.column_stack([p[:, 2] - p[:, 0], (zc[:, 2] - zc[:, 0])[:, None]])
                nrm = np.cross(a, bb)
                slope = np.degrees(np.arctan2(np.hypot(nrm[:, 0], nrm[:, 1]), np.abs(nrm[:, 2])))
                allowed = tol_of(slope, near, far, start, end)
            split = maxerr > allowed
        # The lazy test's share: evaluations of new triangles whose error is in (N, F].
        for t in range(ntri):
            key = tuple(sorted(corners[t]))
            if key not in seen:
                seen.add(key)
                evals += 1
                between += near < maxerr[t] <= far
        new = idx[split[keys]]
        if len(new) == 0:
            break
        verts.extend(new.tolist())
        is_vert[new] = True
    over = np.mean(err > tn + 1e-9) * 100
    return ntri, rnd, over, 100.0 * between / max(evals, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="romsdalen")
    ap.add_argument("--size", type=int, default=401, help="nodes a side, from the window's corner")
    ap.add_argument("--near", type=float, default=2.0)
    ap.add_argument("--far", type=float, default=10.0)
    ap.add_argument("--start", type=float, default=30.0)
    ap.add_argument("--end", type=float, default=30.0)
    ap.add_argument("--rules", default="uniform-F,uniform-N,node,tri,tri+halo,normal")
    a = ap.parse_args()
    from slope_stats import CASES

    z = window(*CASES[a.case])[: a.size, : a.size]
    s = horn(z)
    print(
        f"{a.case} {z.shape[0]}x{z.shape[1]} nodes; N {a.near} F {a.far} "
        f"ramp {a.start}..{a.end} deg; "
        f"steep (>= end) {100 * np.mean(s >= a.end):.1f}%"
    )
    for rule in a.rules.split(","):
        t0 = time.perf_counter()
        ntri, rounds, over, between = run(z, s, rule, a.near, a.far, a.start, a.end)
        print(
            f"  {rule:10s} triangles {ntri:8d} rounds {rounds:3d} nodes over t(s(n)) {over:6.3f}% "
            f"lazy-asked {between:5.1f}%  ({time.perf_counter() - t0:.0f} s)"
        )


if __name__ == "__main__":
    main()
