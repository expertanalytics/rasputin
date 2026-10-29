"""Figures for increment 114's tests, with the design's move: the box's ring
segmentized to EDGE_STEPS steps along its longest side, mapped point by
point (`reprojector`), and the bounds of that outline. Fixtures: the quad
(25833 nodes x 500000..500120 dx 10, y 6599960..6600000 dy 5), the sources
fixture (dataset_fixtures.py) and the fetch fixture (test_cli_fetch.py BBOX).

    python docs/increments/114-probes/box_probe.py
"""

import math
import timeit

import numpy as np
import shapely

from tin_engine.crs import reprojector, same_crs

EDGE_STEPS = 2000


def outline(src: str, dst: str, box: tuple[float, float, float, float]) -> np.ndarray:
    x0, y0, x1, y1 = box
    s = np.linspace(0.0, 1.0, EDGE_STEPS + 1)
    xs = np.concatenate(
        [x0 + (x1 - x0) * s, np.full_like(s, x1), x1 - (x1 - x0) * s, np.full_like(s, x0)]
    )
    ys = np.concatenate(
        [np.full_like(s, y0), y0 + (y1 - y0) * s, np.full_like(s, y1), y1 - (y1 - y0) * s]
    )
    return reprojector(src, dst)(np.column_stack([xs, ys]))


def to_crs(src: str, dst: str, box: tuple[float, float, float, float]):
    if same_crs(src, dst):
        return box
    xy = outline(src, dst, box)
    return tuple(float(v) for v in (*xy.min(axis=0), *xy.max(axis=0)))


X0, Y0 = 500_000.0, 6_600_000.0
utm = (X0 + 22.0, Y0 - 31.0, X0 + 87.0, Y0 - 9.0)  # inside, off-node
lonlat = to_crs("EPSG:25833", "EPSG:4326", utm)
back = to_crs("EPSG:4326", "EPSG:25833", lonlat)
print("utm box", utm)
print("lon/lat box", tuple(round(v, 9) for v in lonlat))
print("back in 25833", tuple(round(v, 4) for v in back))


def snapped(b: tuple[float, float, float, float]) -> tuple[int, int, int, int]:
    """Node columns/rows of the box snapped outward to dx 10, dy 5."""
    c0, c1 = math.floor((b[0] - X0) / 10), math.ceil((b[2] - X0) / 10)
    r0, r1 = math.floor((Y0 - b[3]) / 5), math.ceil((Y0 - b[1]) / 5)
    return c0, c1, r0, r1


for name, b in (("utm", utm), ("back", back)):
    c0, c1, r0, r1 = snapped(b)
    print(name, "cols", c0, c1, "rows", r0, r1, "nodes", (c1 - c0 + 1) * (r1 - r0 + 1))
print("growth m", tuple(round(a - b, 4) for a, b in zip(back, utm, strict=True)))
print(
    "point box 25833->4326",
    tuple(round(v, 9) for v in to_crs("EPSG:25833", "EPSG:4326", (X0, Y0, X0, Y0))),
)
print("CRS84 same as 4326", same_crs("OGC:CRS84", "EPSG:4326"))
n = 200
t = timeit.timeit(lambda: to_crs("EPSG:4326", "EPSG:25833", lonlat), number=n)
print(f"to_crs {1e3 * t / n:.3f} ms per call")

# The sources fixture: coverage in degrees; INSIDE, OUTSIDE and EDGE given in 25833.
COVERAGE = shapely.box(15.0, 59.179575, 15.351879, 59.538349)
for name, b in (
    ("INSIDE", (15.15, 59.33, 15.20, 59.38)),
    ("OUTSIDE", (16.0, 60.0, 16.1, 60.1)),
    ("EDGE", (15.30, 59.30, 15.40, 59.35)),
):
    u = tuple(round(v, 1) for v in to_crs("EPSG:4326", "EPSG:25833", b))
    area = shapely.Polygon(outline("EPSG:25833", "EPSG:4326", u))
    env = shapely.box(*to_crs("EPSG:25833", "EPSG:4326", u))
    print(
        name,
        "in 25833",
        u,
        "share by area",
        round(COVERAGE.intersection(area).area / area.area, 4),
        "by bounds",
        round(COVERAGE.intersection(env).area / env.area, 4),
    )
p = to_crs("EPSG:25833", "EPSG:4326", (508_522.7, 6_576_808.3, 508_522.7, 6_576_808.3))
print(
    "UTM point in 4326",
    tuple(round(v, 6) for v in p),
    "covered",
    COVERAGE.covers(shapely.Point(p[:2])),
)

# The fetch fixture (EPSG:25833).
FETCH = (500_101.3, 6_599_812.3, 500_333.9, 6_599_938.9)
f4326 = to_crs("EPSG:25833", "EPSG:4326", FETCH)
fback = to_crs("EPSG:4326", "EPSG:25833", f4326)
print("fetch box in 4326", tuple(round(v, 9) for v in f4326))
print("fetch box back", tuple(round(v, 4) for v in fback))
print("fetch growth m", tuple(round(a - b, 4) for a, b in zip(fback, FETCH, strict=True)))

# Refusals: the outline's largest longitude step between neighbours.
for name, src, b in (
    ("UTM 60 across 180", "EPSG:32660", (700_000.0, 5_000_000.0, 900_000.0, 5_100_000.0)),
    ("3413 round the pole", "EPSG:3413", (-500_000.0, -500_000.0, 500_000.0, 500_000.0)),
    ("UTM 33 quad", "EPSG:25833", utm),
):
    xy = outline(src, "EPSG:4326", b)
    print(name, "largest longitude step", round(float(np.abs(np.diff(xy[:, 0])).max()), 6))
