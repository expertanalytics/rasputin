"""Increment 34 design probe: the slope together with increment 33's lines, on
the Geilo-Ål window (README item 6; design review round 1, S10).

The simulation of greedy_sim.py with 33's per-triangle distance rule added:
a triangle is allowed the ramp at its distance to Bergensbanen, here the
smallest distance of its nodes and corners (at least 33's exact triangle
distance, so it errs low on triangles). Rules: uniform 1 m; the lines alone
(the quick check's geilo-al-ramp flags: F 20, N 1, ramp 0 to 3000 m); the
lines with the slope (N 2 m, 25 to 35 degrees), both conditions, the node
inserted by section 4.3's rule.

Usage: python lines_sim.py PYTHON RASPUTIN_PY PKG OUT_DIR. Also runs rasputin
master on the window's box (boxes.py's geilo-al domain, written to OUT_DIR)
with the lines alone and at uniform 1 m.
"""

import json
import subprocess
import sys

import numpy as np
import shapely
from greedy_sim import group_first, tol_of
from scipy.spatial import Delaunay
from slope_stats import CASES, horn, window

LINES = "/Users/skavhaug/projects/rasputin_data/banenor_banenettverk/bergensbanen.geojson"
X0, Y1 = 127000.0, 6745000.0  # the window's north-west node
F, NL, S, E, NS, A0, A1 = 20.0, 1.0, 0.0, 3000.0, 2.0, 25.0, 35.0


def ramp(d):
    return np.where(d <= S, NL, np.where(d >= E, F, NL + (F - NL) * (d - S) / max(E - S, 1e-12)))


def run(z, s, dist, rule):
    rows, cols = z.shape
    rr, cc = np.mgrid[0:rows, 0:cols]
    nodes = np.column_stack([cc.ravel() * 10.0, -rr.ravel() * 10.0])
    zf, df = z.ravel(), dist.ravel()
    tn = tol_of(s.ravel(), NS, F, A0, A1)
    corner_ids = [0, cols - 1, (rows - 1) * cols, rows * cols - 1]
    verts, is_vert = list(corner_ids), np.zeros(len(zf), bool)
    is_vert[corner_ids] = True
    for _ in range(300):
        v = np.array(verts)
        tri = Delaunay(nodes[v])
        simp = tri.find_simplex(nodes)
        xform = tri.transform[simp, :2]
        b = np.einsum("ijk,ik->ij", xform, nodes - tri.transform[simp, 2])
        bary = np.column_stack([b, 1 - b.sum(1)])
        corners = v[tri.simplices]
        err = np.abs(zf - (bary * zf[corners[simp]]).sum(1))
        err[is_vert] = 0.0
        ntri = len(tri.simplices)
        maxerr = np.zeros(ntri)
        np.maximum.at(maxerr, simp, err)
        keys, idx = group_first(simp, err)
        if rule == "uniform-1":
            allowed = np.full(ntri, NL)
        else:
            dmin = np.full(ntri, np.inf)
            np.minimum.at(dmin, simp, df)
            allowed = ramp(np.minimum(dmin, df[corners].min(1)))
        split = maxerr > allowed
        pick = np.full(ntri, -1)
        pick[keys] = idx
        if rule == "lines+slope":
            over = np.zeros(ntri, bool)
            np.logical_or.at(over, simp, err > tn)
            o = np.lexsort((-err, -(err / tn), simp))  # largest ratio, ties to the larger error
            first = np.ones(len(o), bool)
            first[1:] = simp[o][1:] != simp[o][:-1]
            rk, ri = simp[o][first], o[first]
            pick[rk[over[rk]]] = ri[over[rk]]
            split |= over
        new = pick[np.nonzero(split)[0]]
        new = new[new >= 0]
        if len(new) == 0:
            break
        verts.extend(new.tolist())
        is_vert[new] = True
    return ntri


def main():
    py, rasputin, pkg, out = sys.argv[1:5]
    z = window(*CASES["geilo-al"])[:1001, :1001]
    s = horn(z)
    with open(LINES) as f:
        geoms = [shapely.geometry.shape(g["geometry"]) for g in json.load(f)["features"]]
    rows, cols = z.shape
    rr, cc = np.mgrid[0:rows, 0:cols]
    pts = shapely.points(X0 + cc.ravel() * 10.0, Y1 - rr.ravel() * 10.0)
    dist = shapely.distance(pts, shapely.union_all(geoms)).reshape(z.shape)
    print(
        f"geilo-al window: nodes within 3000 m of the line {100 * np.mean(dist < 3000):.1f}%, "
        f"at 35 deg or more {100 * np.mean(s >= 35):.1f}%"
    )
    for rule in ("uniform-1", "lines", "lines+slope"):
        print(f"  simulation {rule:11s} {run(z, s, dist, rule)} triangles")
    dom = f"{out}/geilo-al_box.geojson"  # written by boxes.py
    dem = "/Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925"
    for name, extra in (("uniform-1", ["--tolerance", "1"]),
                        ("lines", ["--tolerance", "20", "--tolerance-near", LINES, "1",
                                   "--tolerance-ramp", "0", "3000"])):  # fmt: skip
        r = subprocess.run([py, rasputin, "mesh", "--dem", dem, "--domain", dom, *extra,
                            "--out", f"{out}/geilo_{name}.vtk"],
                           env={"PYTHONPATH": pkg}, capture_output=True, text=True)  # fmt: skip
        said = [x for x in (r.stdout + r.stderr).splitlines() if "triangles." in x or "rror" in x]
        print(f"  rasputin {name:11s} {said}")


if __name__ == "__main__":
    main()
