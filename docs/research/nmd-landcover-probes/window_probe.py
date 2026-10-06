"""THROWAWAY probe, not production code: NMD 2018 against CORINE 2018 in one window of Lagan.

Backs the figures in docs/research/nmd-landcover.md. Run with the repository's
venv from a directory holding Kronoberg's county file, unzipped:

  curl -O https://geodata.naturvardsverket.se/nedladdning/marktacke/NMD2018/\
bas_lan_ogen/G_lan_nmd2018bas_ogeneraliserad_v1_1.zip && unzip G_lan_*.zip
  python window_probe.py ../rasputin_data 430000 6295000 440000 6305000 3000

The arguments are Ola's data folder, the window in SWEREF 99 TM (EPSG:3006) and the side in
metres of the sub-window (its north-west corner) on which polygons are built
and simplified. The CORINE extract and the Lagan outline are read from the data
folder (sweden_corine/, sweden_smhi_svar/).

The sieve relabels each region under the minimum area to the class it shares
most boundary with, among neighbours already at least that large, raising the
threshold in doublings (GRASS r.reclass.area's rmarea rule, roughly).
"""

import json
import sys
import time
from collections import Counter, deque
from pathlib import Path

import numpy as np
import shapely
import tifffile
from pyproj import Transformer
from shapely.geometry import box, shape

TIF = "G_lan_nmd2018bas_ogeneraliserad_v1_1/G_lan_nmd2018bas_ogeneraliserad_v1_1.tif"
ORX, ORY = 396130.0, 6343950.0  # the file's tie point; cell corners on 10 m multiples


def label(w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """4-connected regions of equal class: a label per cell, a size per region."""
    h, wd = w.shape
    lab = -np.ones((h, wd), np.int64)
    sizes = []
    for r in range(h):
        for c in range(wd):
            if lab[r, c] >= 0:
                continue
            v, n = w[r, c], len(sizes)
            lab[r, c] = n
            q, s = deque([(r, c)]), 0
            while q:
                y, x = q.popleft()
                s += 1
                for yy, xx in ((y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1)):
                    if 0 <= yy < h and 0 <= xx < wd and lab[yy, xx] < 0 and w[yy, xx] == v:
                        lab[yy, xx] = n
                        q.append((yy, xx))
            sizes.append(s)
    return lab, np.array(sizes)


def crack_km(w: np.ndarray) -> float:
    """Length of class boundary inside the window: cell edges between two classes."""
    n = (w[:, 1:] != w[:, :-1]).sum() + (w[1:, :] != w[:-1, :]).sum()
    return float(n) * 10 / 1000


def corners(w: np.ndarray) -> int:
    """Lattice points where a boundary turns or three classes meet: the polygon vertices."""
    p = np.pad(w, 1, mode="edge")
    a, b, c, d = p[:-1, :-1], p[:-1, 1:], p[1:, :-1], p[1:, 1:]
    flat = (a == b) & (b == c) & (c == d)
    straight = ((a == c) & (b == d) & (a != b)) | ((a == b) & (c == d) & (a != c))
    return int((~(flat | straight)).sum())


def _sieve_round(w: np.ndarray, min_cells: int) -> bool:
    lab, sizes = label(w)
    small = sizes < min_cells
    if not small.any():
        return False
    h, wd = w.shape
    votes: dict[int, Counter[int]] = {}
    for dy, dx in ((0, 1), (1, 0)):
        la, lb = lab[: h - dy, : wd - dx], lab[dy:, dx:]
        va, vb = w[: h - dy, : wd - dx], w[dy:, dx:]
        m = la != lb
        for x, y, cv in ((la[m], lb[m], vb[m]), (lb[m], la[m], va[m])):
            sel = small[x] & ~small[y]
            for comp, cls in zip(x[sel].tolist(), cv[sel].tolist(), strict=True):
                votes.setdefault(comp, Counter())[cls] += 1
    new = -np.ones(len(sizes), np.int64)
    for comp, cnt in votes.items():
        new[comp] = cnt.most_common(1)[0][0]
    m = small[lab] & (new[lab] >= 0)
    w[m] = new[lab][m]
    return bool(m.any())


def sieve(w: np.ndarray, min_cells: int) -> np.ndarray:
    w = w.copy()
    t, steps = 2, []
    while t < min_cells:
        steps.append(t)
        t *= 2
    for th in [*steps, min_cells]:
        for _ in range(20):
            if not _sieve_round(w, th):
                break
    return w


def polygons(w: np.ndarray, ox: float, oy: float) -> list[shapely.Geometry]:
    """The window as polygons: the union of each class's row runs."""
    h, wd = w.shape
    per: dict[int, list[shapely.Geometry]] = {}
    for r in range(h):
        row = w[r]
        cut = np.flatnonzero(np.diff(row)) + 1
        for s, e in zip(np.r_[0, cut], np.r_[cut, wd], strict=True):
            cell = box(ox + 10 * s, oy - 10 * (r + 1), ox + 10 * e, oy - 10 * r)
            per.setdefault(int(row[s]), []).append(cell)
    out: list[shapely.Geometry] = []
    for boxes in per.values():
        u = shapely.union_all(boxes)
        out.extend(getattr(u, "geoms", [u]))
    return out


def nverts(polys: list[shapely.Geometry]) -> int:
    return int(sum(shapely.get_num_coordinates(p) for p in polys))


def main() -> None:
    data = Path(sys.argv[1])
    x0, y0, x1, y1, sub = map(float, sys.argv[2:7])
    win = box(x0, y0, x1, y1)
    outline = data / "sweden_smhi_svar" / "lagan_mouth_svar2022_3006.geojson"
    dom = shape(json.loads(outline.read_text())["features"][0]["geometry"])
    print("share of the window inside Lagan:", round(win.intersection(dom).area / win.area, 3))
    a = tifffile.imread(TIF)
    c0, r0, c1, r1 = (
        int(v) for v in ((x0 - ORX) / 10, (ORY - y1) / 10, (x1 - ORX) / 10, (ORY - y0) / 10)
    )
    w = a[r0:r1, c0:c1].copy()
    _, sizes = label(w)
    print(
        f"NMD raw: {len(sizes)} regions, median {int(np.median(sizes))} cells, "
        f"{int((sizes == 1).sum())} of one cell; "
        f"boundary {crack_km(w):.0f} km, {corners(w)} vertices"
    )
    for mc in (25, 100, 2500):
        t = time.time()
        s = sieve(w, mc)
        _, z = label(s)
        print(
            f"NMD sieved at {mc / 100} ha: {len(z)} regions, smallest {int(z.min())} cells; "
            f"boundary {crack_km(s):.0f} km, {corners(s)} vertices ({time.time() - t:.0f} s)"
        )
    shares = {int(k): round(v / w.size, 3) for k, v in sorted(Counter(w.ravel().tolist()).items())}
    print("NMD shares", shares)

    tr = Transformer.from_crs(3035, 3006, always_xy=True)
    fc = json.loads((data / "sweden_corine" / "lagan_clc2018_3035.geojson").read_text())
    sub_box = box(x0, y1 - sub, x0 + sub, y1)
    clipped, inside, in_sub, codes = [], [], [], []
    for f in fc["features"]:
        g = shapely.transform(
            shape(f["geometry"]), lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1]))
        )
        if win.contains(g):
            inside.append(g.area / 1e4)
        for target, acc in ((win, clipped), (sub_box, in_sub)):
            if g.intersects(target):
                gi = g.intersection(target)
                parts = [p for p in getattr(gi, "geoms", [gi]) if p.area > 0]
                acc.extend(parts)
                if acc is clipped:
                    codes.extend([f["properties"]["Code_18"]] * len(parts))
    areas = np.array([p.area for p in clipped]) / 1e4
    print(
        f"CORINE: {len(clipped)} polygons after the clip, {nverts(clipped)} vertices, median "
        f"{np.median(areas):.1f} ha; {len(inside)} wholly inside, smallest {min(inside):.1f} ha"
    )
    # Each shared edge is in two clipped perimeters, the window's frame in one.
    edge_km = (sum(p.length for p in clipped) - win.length) / 2 / 1000
    area: Counter[str] = Counter()
    for p, code in zip(clipped, codes, strict=True):
        area[code] += p.area
    cshares = {k: round(v / win.area, 3) for k, v in sorted(area.items())}
    print(f"CORINE: boundary {edge_km:.0f} km, {len(cshares)} classes; shares {cshares}")
    print(
        f"CORINE in the {sub / 1000:.0f} km sub-window: "
        f"{len(in_sub)} polygons, {nverts(in_sub)} vertices"
    )

    n = int(sub / 10)
    cases = (("raw", 1), ("sieve 0.25 ha", 25), ("sieve 1 ha", 100), ("sieve 25 ha", 2500))
    for name, mc in cases:
        ww = w[:n, :n] if mc == 1 else sieve(w[:n, :n], mc)
        polys = polygons(ww, x0, y1)
        line = f"NMD {sub / 1000:.0f} km sub-window, {name}: "
        line += f"{len(polys)} polygons, {nverts(polys)} vertices"
        for tol in (5, 10, 20, 40):
            simple = shapely.coverage_simplify(
                np.array(polys, dtype=object), tol, simplify_boundary=False
            )
            line += f" | {tol} m: {nverts(list(simple))}"
        print(line)


if __name__ == "__main__":
    main()
