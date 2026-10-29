"""Figures for increment 114's tests: a box in another CRS, moved into the
frame by transform_bounds, on the quad fixture (25833 nodes x 500000..500120
dx 10, y 6599960..6600000 dy 5) and on a point box for sources.

    python docs/increments/114-probes/box_probe.py
"""

import math
import timeit

from tin_engine.crs import same_crs, transform_bounds

X0, Y0 = 500_000.0, 6_600_000.0
utm = (X0 + 22.0, Y0 - 31.0, X0 + 87.0, Y0 - 9.0)  # inside, off-node
lonlat = transform_bounds("EPSG:25833", "EPSG:4326", utm)
back = transform_bounds("EPSG:4326", "EPSG:25833", lonlat)
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
point = (10.5, 59.5, 10.5, 59.5)
print("point box 4326->4326", transform_bounds("EPSG:4326", "EPSG:4326", point))
print("CRS84 same as 4326", same_crs("OGC:CRS84", "EPSG:4326"))
n = 200
t = timeit.timeit(lambda: transform_bounds("EPSG:4326", "EPSG:25833", lonlat), number=n)
print(f"transform_bounds {1e3 * t / n:.3f} ms per call")

# The sources fixture (tests/python/dataset_fixtures.py): coverage and boxes in degrees.
COVERAGE = (15.0, 59.179575, 15.351879, 59.538349)
INSIDE = (15.15, 59.33, 15.20, 59.38)
OUTSIDE = (16.0, 60.0, 16.1, 60.1)
for name, b in (("INSIDE", INSIDE), ("OUTSIDE", OUTSIDE)):
    u = transform_bounds("EPSG:4326", "EPSG:25833", b)
    again = transform_bounds("EPSG:25833", "EPSG:4326", u)
    inside = (again[0] >= COVERAGE[0] and again[1] >= COVERAGE[1]
              and again[2] <= COVERAGE[2] and again[3] <= COVERAGE[3])  # fmt: skip
    print(name, "in 25833", tuple(round(v, 1) for v in u))
    print(name, "back in 4326", tuple(round(v, 6) for v in again), "inside coverage", inside)

# The fetch fixture (tests/python/test_cli_fetch.py BBOX, EPSG:25833).
FETCH = (500_101.3, 6_599_812.3, 500_333.9, 6_599_938.9)
f4326 = transform_bounds("EPSG:25833", "EPSG:4326", FETCH)
fback = transform_bounds("EPSG:4326", "EPSG:25833", f4326)
print("fetch box in 4326", tuple(round(v, 9) for v in f4326))
print("fetch box back", tuple(round(v, 4) for v in fback))
print("fetch growth m", tuple(round(a - b, 4) for a, b in zip(fback, FETCH, strict=True)))
