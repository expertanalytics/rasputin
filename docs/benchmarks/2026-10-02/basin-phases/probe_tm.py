"""Probe for README.md (basin-phases); a measurement script, not production code.
Run from this directory with the venv that has build-bench/pkg on its path."""
import ctypes, gc, json, math, os, sys, tracemalloc
from pathlib import Path
import shapely.geometry
from tin_engine.crs import crs_label, reprojector
from tin_engine.domain import DomainPolygon
from tin_engine.io.repository import CacheRepository
from tin_engine.mosaic import assemble, plan_mosaic
from tin_engine.target_grid import default_spacing, source_region, target_grid_for
_libc, _buf = ctypes.CDLL("/usr/lib/libSystem.B.dylib"), ctypes.create_string_buffer(512)
def fp():
    _libc.proc_pid_rusage(os.getpid(), 4, _buf); return int.from_bytes(_buf.raw[72:80], "little")/1e9
D=Path("/Users/skavhaug/projects/rasputin_data")
target = crs_label(Path("runs/basin_out_crs.wkt").read_text().strip())
geo = json.loads((D/"sao_francisco_piece/bho2017_50k_level3/761_epsg4674.geojson").read_text())["features"][0]["geometry"]
domain = DomainPolygon(polygon=shapely.geometry.shape(geo), crs="EPSG:4674").to_crs(target)
repo = CacheRepository(D/"cache", "anadem-v1"); fps = repo.footprints(); meta = fps[0].meta
c = domain.polygon.centroid; ((ax, ay),) = reprojector(target, meta.crs)([(c.x, c.y)])
grid = target_grid_for(domain, target, default_spacing(meta, (ax, ay)))
grown = domain.polygon.buffer(math.sqrt(2) * grid.spacing, join_style="mitre")
box, needed = source_region(grid, meta, grown); plan = plan_mosaic(fps, box, needed)
print("planned", fp(), "plan tiles", len(plan.tiles), plan.tiles[0].meta == plan.meta, plan.tiles[0].source)
tracemalloc.start()
t = repo.load_window(plan.tiles[0].name, plan.tiles[0].source)
print("loaded", fp(), "traced cur/peak GB", [x/1e9 for x in tracemalloc.get_traced_memory()])
for s in tracemalloc.take_snapshot().statistics("lineno")[:6]: print(s)
del t; gc.collect(); print("del tile", fp(), [x/1e9 for x in tracemalloc.get_traced_memory()])
import subprocess
out = subprocess.run(["vmmap", "-summary", str(os.getpid())], capture_output=True, text=True).stdout
print("\n".join(l for l in out.splitlines() if "MALLOC" in l or "Physical footprint" in l or "REGION TYPE" in l or "TOTAL" in l))
