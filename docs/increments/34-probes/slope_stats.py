"""Increment 34 design probe: how much ground is steep, by four slope measures,
and how DEM noise moves it (docs/increments/34-probes/README.md, item 1).

Reads DTM10 windows straight from the tiles (tifffile), so nothing of rasputin
runs here. Measures, all in degrees, per node:
  horn     Horn (1981) 3x3 weighted differences, by section 3's rule at the
           border and next to NoData (reflection through the node);
  central  centred differences over the 4 neighbours (Zevenbergen-Thorne's gradient);
  cellmax  the largest bilinear-cell gradient over the 4 cells sharing the node,
           the bound refine.hpp's foot_epsilon uses (|dz| along each cell edge);
  horn_s   Horn after a 3x3 mean filter.
Noise: Gaussian noise of sigma metres added to z, Horn recomputed.
"""

import sys
import time

import numpy as np
import tifffile

DTM = "/Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925"
CASES = {
    # name: (tile, x0, x1, y0, y1) in EPSG:25833, 10 m nodes
    "romsdalen": ("6901_3", 118000, 128000, 6938000, 6948000),
    "geilo-al": ("6701_3", 127000, 150000, 6720000, 6745000),
}


def window(tile, x0, x1, y0, y1):
    with open(f"{DTM}/{tile}_10m_z33.tfw") as f:
        tfw = [float(v) for v in f]
    ox, oy = tfw[4], tfw[5]
    c0, c1 = int((x0 - ox) / 10), int((x1 - ox) / 10) + 1
    r0, r1 = int((oy - y1) / 10), int((oy - y0) / 10) + 1
    z = tifffile.imread(f"{DTM}/{tile}_10m_z33.tif")[r0:r1, c0:c1].astype(np.float64)
    return z


def horn(z, dx=10.0, dy=None, valid=None):
    """Section 3's slope in degrees: Horn's 3x3 differences with every missing
    neighbour (outside the grid, or NoData) filled so that a plane is exact:
    an edge neighbour by its reflection through the node, 2 z(node) -
    z(opposite), or z(node) when the opposite is missing too; a corner
    neighbour by its reflection when the opposite corner is there, else by
    z(row neighbour) + z(column neighbour) - z(node), from the filled edge
    neighbours. A NoData node gets 0."""
    dy = dx if dy is None else dy
    z = np.asarray(z, dtype=np.float64)
    ok = np.ones(z.shape, bool) if valid is None else np.asarray(valid, bool)
    p = np.pad(np.where(ok, z, np.nan), 1, mode="constant", constant_values=np.nan)
    rows, cols = z.shape

    def at(di, dj):  # the neighbour at row offset di, column offset dj, NaN if missing
        return p[1 + di : 1 + di + rows, 1 + dj : 1 + dj + cols]

    def edge(di, dj):
        v, o = at(di, dj), at(-di, -dj)
        return np.where(np.isnan(v), np.where(np.isnan(o), z, 2.0 * z - o), v)

    e = {k: edge(*k) for k in ((-1, 0), (1, 0), (0, -1), (0, 1))}

    def corner(di, dj):
        v, o = at(di, dj), at(-di, -dj)
        return np.where(
            np.isnan(v), np.where(np.isnan(o), e[(di, 0)] + e[(0, dj)] - z, 2.0 * z - o), v
        )

    a, b, c = corner(-1, -1), e[(-1, 0)], corner(-1, 1)
    d, f = e[(0, -1)], e[(0, 1)]
    g, h, i = corner(1, -1), e[(1, 0)], corner(1, 1)
    gx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * dx)
    gy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * dy)
    return np.where(ok, np.degrees(np.arctan(np.hypot(gx, gy))), 0.0)


def horn_edge(z, d=10.0):
    """Round 1's rule, kept for the comparison of item 1b: a missing neighbour
    takes the node's own z (numpy's edge padding)."""
    p = np.pad(z, 1, mode="edge")
    a, b, c = p[:-2, :-2], p[:-2, 1:-1], p[:-2, 2:]
    dd, f = p[1:-1, :-2], p[1:-1, 2:]
    g, h, i = p[2:, :-2], p[2:, 1:-1], p[2:, 2:]
    gx = ((c + 2 * f + i) - (a + 2 * dd + g)) / (8 * d)
    gy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * d)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def central(z, d=10.0):
    gy, gx = np.gradient(z, d)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def cellmax(z, d=10.0):
    z00, z01, z10, z11 = z[:-1, :-1], z[:-1, 1:], z[1:, :-1], z[1:, 1:]
    gx = np.maximum(abs(z01 - z00), abs(z11 - z10)) / d
    gy = np.maximum(abs(z10 - z00), abs(z11 - z01)) / d
    cell = np.hypot(gx, gy)  # one value per cell
    p = np.pad(cell, 1, mode="constant", constant_values=0.0)
    node = np.maximum.reduce([p[:-1, :-1], p[:-1, 1:], p[1:, :-1], p[1:, 1:]])
    return np.degrees(np.arctan(node))


def smooth(z):
    p = np.pad(z, 1, mode="edge")
    return sum(p[i : i + z.shape[0], j : j + z.shape[1]] for i in range(3) for j in range(3)) / 9.0


def share(s, angles):
    return " ".join(f"{a}:{100 * np.mean(s >= a):5.1f}%" for a in angles)


def main():
    angles = (15, 20, 25, 30, 35, 40, 45)
    rng = np.random.default_rng(34)
    for name, spec in CASES.items():
        z = window(*spec)
        print(f"{name}: {z.shape[0]} x {z.shape[1]} nodes, z {z.min():.0f}..{z.max():.0f} m")
        t = time.perf_counter()
        sh = horn(z)
        print(
            f"  horn    {share(sh, angles)}  "
            f"({(time.perf_counter() - t) * 1e9 / z.size:.1f} ns/node)"
        )
        print(f"  central {share(central(z), angles)}")
        print(f"  cellmax {share(cellmax(z), angles)}")
        print(f"  horn_s  {share(horn(smooth(z)), angles)}")
        for sigma in (0.1, 0.3, 1.0):
            sn = horn(z + rng.normal(0.0, sigma, z.shape))
            moved = np.mean(np.abs(sn - sh))
            print(f"  noise {sigma} m: {share(sn, angles)}  mean |change| {moved:.2f} deg")
        flat = sh < 5
        print(
            f"  ground under 5 deg: {100 * flat.mean():.1f}%; p50/p90/p99 horn "
            f"{np.percentile(sh, 50):.1f}/{np.percentile(sh, 90):.1f}/"
            f"{np.percentile(sh, 99):.1f} deg"
        )


if __name__ == "__main__":
    sys.exit(main())
