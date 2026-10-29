"""Figures for increment 114's design review fixes: how far pyproj's
`transform_bounds` falls inside the true image of a box's edges, by
densify_pts, on a São Francisco-sized box (degrees to EPSG:32723) and on the
test fixtures; a box whose geographic image crosses ±180°; and the share of a
projected box a coverage holds, by area and by bounds.

    python docs/increments/114-probes/edge_probe.py
"""

import timeit

import numpy as np
import shapely
from pyproj import Transformer

from tin_engine.crs import reprojector, transform_bounds
from tin_engine.source_report import _lonlat


def true_bounds(src: str, dst: str, box: tuple[float, float, float, float], n: int = 200_001):
    """The image's bounds from `n` points per edge: the reference."""
    t = Transformer.from_crs(src, dst, always_xy=True)
    x0, y0, x1, y1 = box
    s = np.linspace(0.0, 1.0, n)
    xs = np.concatenate([x0 + (x1 - x0) * s, np.full(n, x1), x1 - (x1 - x0) * s, np.full(n, x0)])
    ys = np.concatenate([np.full(n, y0), y0 + (y1 - y0) * s, np.full(n, y1), y1 - (y1 - y0) * s])
    u, v = t.transform(xs, ys)
    return float(np.min(u)), float(np.min(v)), float(np.max(u)), float(np.max(v))


def shortfall(src, dst, box, pts):
    t = Transformer.from_crs(src, dst, always_xy=True)
    got = t.transform_bounds(*box, densify_pts=pts)
    ref = true_bounds(src, dst, box)
    # positive: the edge lies inside the true image by that much
    return (got[0] - ref[0], got[1] - ref[1], ref[2] - got[2], ref[3] - got[3])


SF = (-48.0, -21.0, -36.0, -7.0)  # São Francisco basin's box, degrees, about 1300 x 1550 km
for pts in (21, 101, 1001):
    s = shortfall("EPSG:4326", "EPSG:32723", SF, pts)
    t = Transformer.from_crs("EPSG:4326", "EPSG:32723", always_xy=True)
    ms = (
        1e3
        * timeit.timeit(lambda t=t, pts=pts: t.transform_bounds(*SF, densify_pts=pts), number=50)
        / 50
    )
    print(f"SF 4326->32723 pts {pts}: shortfall m {tuple(round(v, 4) for v in s)} {ms:.3f} ms")
# the other way: a UTM box of that size to degrees, shortfall in metres (x 111 km/deg)
SFU = transform_bounds("EPSG:4326", "EPSG:32723", SF)
for pts in (21, 101, 1001):
    s = shortfall("EPSG:32723", "EPSG:4326", SFU, pts)
    print(f"SF 32723->4326 pts {pts}: shortfall m {tuple(round(v * 111_320, 4) for v in s)}")
# the test fixtures (quad box, fetch box), 25833 <-> 4326, at 21 points
QUAD = (500_022.0, 6_599_969.0, 500_087.0, 6_599_991.0)
FETCH = (500_101.3, 6_599_812.3, 500_333.9, 6_599_938.9)
for name, b in (("quad", QUAD), ("fetch", FETCH)):
    ll = transform_bounds("EPSG:25833", "EPSG:4326", b)
    s = shortfall("EPSG:4326", "EPSG:25833", ll, 21)
    print(f"{name} 4326->25833 pts 21: shortfall m {tuple(round(v, 6) for v in s)}")

# a UTM zone 60 box across 180 degrees east
AM = (700_000.0, 5_000_000.0, 900_000.0, 5_100_000.0)
print(
    "EPSG:32660 box",
    AM,
    "-> 4326",
    tuple(round(v, 2) for v in transform_bounds("EPSG:32660", "EPSG:4326", AM)),
)

# area share against bounds share: the sources COVERAGE, a 25833 box over its east edge
COVERAGE = shapely.box(15.0, 59.179575, 15.351879, 59.538349)
EDGE = transform_bounds("EPSG:4326", "EPSG:25833", (15.30, 59.30, 15.40, 59.35))
poly = _lonlat(reprojector("EPSG:25833", "EPSG:4326"), EDGE)
env = shapely.box(*transform_bounds("EPSG:25833", "EPSG:4326", EDGE))
print("EDGE in 25833", tuple(round(v, 1) for v in EDGE))
print(
    "EDGE share by area",
    round(COVERAGE.intersection(poly).area / poly.area, 4),
    "by bounds",
    round(COVERAGE.intersection(env).area / env.area, 4),
)


# The design's move: `n` equal steps along each side of the box, through
# `reprojector`, its bounds; and the jump test for ±180°.
def outline(src, dst, box, n):
    x0, y0, x1, y1 = box
    s = np.linspace(0.0, 1.0, n + 1)
    xs = np.concatenate(
        [x0 + (x1 - x0) * s, np.full_like(s, x1), x1 - (x1 - x0) * s, np.full_like(s, x0)]
    )
    ys = np.concatenate(
        [np.full_like(s, y0), y0 + (y1 - y0) * s, np.full_like(s, y1), y1 - (y1 - y0) * s]
    )
    return reprojector(src, dst)(np.column_stack([xs, ys]))


def jumps(xy):
    return float(np.abs(np.diff(xy[:, 0])).max())


for src, dst, box, scale in (
    ("EPSG:4326", "EPSG:32723", SF, 1.0),
    ("EPSG:32723", "EPSG:4326", SFU, 111_320.0),
):
    ref = true_bounds(src, dst, box)
    for n in (1000, 2000, 4000):
        xy = outline(src, dst, box, n)
        got = (xy[:, 0].min(), xy[:, 1].min(), xy[:, 0].max(), xy[:, 1].max())
        s = (got[0] - ref[0], got[1] - ref[1], ref[2] - got[2], ref[3] - got[3])
        ms = 1e3 * timeit.timeit(lambda a=(src, dst, box, n): outline(*a), number=20) / 20
        print(
            f"outline {src}->{dst} n {n} points {len(xy)}: shortfall m "
            f"{tuple(round(float(v) * scale, 5) for v in s)} {ms:.3f} ms"
        )
print("AM outline max lon jump", round(jumps(outline("EPSG:32660", "EPSG:4326", AM, 2000)), 2))
POLE = (-500_000.0, -500_000.0, 500_000.0, 500_000.0)  # EPSG:3413 box around the north pole
print("pole outline max lon jump", round(jumps(outline("EPSG:3413", "EPSG:4326", POLE, 2000)), 2))
print("SF outline max lon jump", round(jumps(outline("EPSG:32723", "EPSG:4326", SFU, 2000)), 5))
print(
    "pyproj pole box", tuple(round(v, 2) for v in transform_bounds("EPSG:3413", "EPSG:4326", POLE))
)
