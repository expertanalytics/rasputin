"""The decode phase of the reprojected path (dem_input._open_reprojected), step by step,
on master's code, plus the cheap superset candidates for the grown region.
Layout as run.py. Usage: python decode_probe.py <domain.geojson> <cache root> <target CRS>
"""

import math
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.meta_path[:] = [f for f in sys.meta_path if "skbc" not in type(f).__module__]
sys.path.insert(0, str(HERE / "pkg"))

import shapely  # noqa: E402

from tin_engine import dem_input as di  # noqa: E402
from tin_engine import target_grid as tg  # noqa: E402
from tin_engine.crs import reprojector  # noqa: E402
from tin_engine.io.domain_file import read_domain  # noqa: E402
from tin_engine.mosaic import assemble, plan_mosaic  # noqa: E402

print("tin_engine from", di.__file__)
domain_path, cache, target = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]


def step(name, f):
    t0 = time.perf_counter()
    r = f()
    print(f"  {name:<52} {time.perf_counter() - t0:8.3f} s", flush=True)
    return r


given = step("read_domain", lambda: read_domain(domain_path, None))
cached = di.CachedSource(source="glo30", cache=cache)
req = di.DemRequest(cached=cached, domain=given, target_crs=target)
repo, label = di.repository_for(req.sources, req.nodata, req.cached)
footprints = step("footprints (headers)", repo.footprints)
domain = step("domain.to_crs(target)", lambda: given.to_crs(target))
meta = footprints[0].meta
c = domain.polygon.centroid
((ax, ay),) = reprojector(target, meta.crs)([(c.x, c.y)])
spacing = tg.default_spacing(meta, (ax, ay))
grid = step("target_grid_for (buffer #1)", lambda: tg.target_grid_for(domain, target, spacing))
d = math.sqrt(2) * grid.spacing
grown = step("grown = buffer #2 (same as #1)", lambda: domain.polygon.buffer(d, join_style="mitre"))
box, needed = step(
    "source_region (transform + buffer #3)", lambda: tg.source_region(grid, meta, grown)
)
plan = step("plan_mosaic", lambda: plan_mosaic(footprints, box, needed))
step("repository.check", lambda: repo.check(plan))
mosaic = step(
    "assemble (decode + seams)",
    lambda: assemble(plan, repo.load, needed, load_window=repo.load_window),
)
source = tg.TileWindows(mosaic.tile)
tile = step("resample", lambda: tg.resample(grid, source, os.cpu_count() or 1))
nodes = grid.rows * grid.cols / 1e6
print(f"  grid {grid.rows} x {grid.cols} = {nodes:.1f} M nodes; canvas {mosaic.tile.array.shape}")
print(
    f"  grown vertices {len(shapely.get_coordinates(grown))}, "
    f"needed vertices {len(shapely.get_coordinates(needed))}"
)

# Cheap supersets of the grown region.
for eps in (float(grid.spacing), 10.0):

    def wider(e: float = eps) -> shapely.Geometry:
        return shapely.simplify(domain.polygon, e).buffer(d + e, join_style="mitre")

    s = step(f"simplify({eps:g}) + buffer(d + {eps:g}, mitre)", wider)
    extra = (s.area - grown.area) / 1e6
    print(
        f"    covers the exact grown region: {s.buffer(1e-6).covers(grown)}; "
        f"vertices {len(shapely.get_coordinates(s))}; extra area {extra:.2f} km2"
    )
hull = domain.polygon.convex_hull
h = step("buffer of convex hull (mitre)", lambda: hull.buffer(d, join_style="mitre"))
print(f"    extra area {(h.area - grown.area) / 1e6:.1f} km2 of {grown.area / 1e6:.0f}")
