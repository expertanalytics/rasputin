"""Each buffer call of `decode` and `features read`, in isolation, with the
arguments master f81b20b7 passes (read off the traces in raw/trace/), and the
decode and parse work around them, step by step. Every step is run three times
and the median printed with all three; nothing else runs.

Usage: PKG=<dir> python isolated.py <lagan|ljungan_flasjo|numedalslagen>
Prints JSON lines: {"catchment", "step", "runs": [s, s, s], "median", ...}.
"""

import json
import math
import os
import statistics
import sys
import time
from pathlib import Path

sys.meta_path[:] = [f for f in sys.meta_path if "skbc" not in type(f).__module__]
sys.path.insert(0, str(Path(os.environ["PKG"]).resolve()))

import shapely  # noqa: E402

import tin_engine  # noqa: E402
from tin_engine import dem_input as di  # noqa: E402
from tin_engine import feature_input as fi  # noqa: E402
from tin_engine import target_grid as tg  # noqa: E402
from tin_engine.crs import reprojector  # noqa: E402
from tin_engine.io.domain_file import read_domain  # noqa: E402
from tin_engine.io.geopackage import layer_info, query_features  # noqa: E402
from tin_engine.io.repository import open_geopackage, read_json  # noqa: E402
from tin_engine.mosaic import assemble, plan_mosaic  # noqa: E402

DATA = Path("/Users/skavhaug/projects/rasputin_data")
EU = DATA / "corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg"
c = sys.argv[1]
print(json.dumps({"tin_engine": tin_engine.__file__, "shapely": shapely.__version__,
                  "geos": shapely.geos_version_string}), flush=True)  # fmt: skip


def step(name, f, repeat=3, **extra):  # type: ignore[no-untyped-def]
    runs, r = [], None
    for _ in range(repeat):
        t0 = time.perf_counter()
        r = f()
        runs.append(time.perf_counter() - t0)
    row = {"catchment": c, "step": name, "runs": runs, "median": statistics.median(runs), **extra}
    print(json.dumps(row), flush=True)
    return r


def nverts(rows):  # type: ignore[no-untyped-def]
    return int(sum(shapely.get_num_coordinates(g) for _, g, _ in rows if g is not None))


if c == "numedalslagen":
    dem = DATA / "DTM10_UTM33_20260925"
    given = read_domain(
        Path("/Users/skavhaug/projects/rasputin_scratch/norway/numedalslagen/numedalslagen_outline_nve.geojson")
    )
    repo, _ = di.repository_for((dem,))
    footprints = step("footprints (tile headers)", repo.footprints)
    plan, domain, grown = step("_domain_plan (dem_input.py:250's buffer inside)",
                               lambda: di._domain_plan(footprints, given))  # fmt: skip
    diag = math.hypot(plan.meta.delta_x, plan.meta.delta_y)
    step("buffer dem_input.py:250", lambda: domain.polygon.buffer(diag, join_style="mitre"),
         distance=diag, vertices=int(shapely.get_num_coordinates(domain.polygon)))  # fmt: skip
    step("repository.check", lambda: repo.check(plan))
    step("assemble (decode + seams)", lambda: assemble(plan, repo.load, grown, load_window=repo.load_window))
    dem_crs = plan.meta.crs
    sources = [("gpkg33", DATA / "corine2018_dtm10_utm33.gpkg", "corine2018", "Code_18"),
               ("gpkg", EU, "U2018_CLC2018_V2020_20u1", "Code_18")]  # fmt: skip
    geojson = None
else:
    cache = DATA / "sweden_glo30_cache"
    outline = DATA / f"sweden_smhi_svar/{'lagan_mouth' if c == 'lagan' else c}_svar2022_3006.geojson"
    target = "EPSG:3006"
    given = read_domain(outline, None)
    req = di.DemRequest(cached=di.CachedSource(source="glo30", cache=cache), domain=given, target_crs=target)
    repo, _ = di.repository_for(req.sources, req.nodata, req.cached)
    footprints = step("footprints (tile headers)", repo.footprints)
    domain = step("domain.to_crs(target)", lambda: given.to_crs(target))
    meta = footprints[0].meta
    cen = domain.polygon.centroid
    ((ax, ay),) = reprojector(target, meta.crs)([(cen.x, cen.y)])
    spacing = tg.default_spacing(meta, (ax, ay))
    grid = step("target_grid_for (target_grid.py:107's buffer inside)",
                lambda: tg.target_grid_for(domain, target, spacing))  # fmt: skip
    d = math.sqrt(2) * grid.spacing
    grown = step("buffer target_grid.py:107 = dem_input.py:213 (the same call)",
                 lambda: domain.polygon.buffer(d, join_style="mitre"), distance=d,
                 vertices=int(shapely.get_num_coordinates(domain.polygon)))  # fmt: skip
    box, needed = step("target_grid.source_region (target_grid.py:145's buffer inside)",
                       lambda: tg.source_region(grid, meta, grown))  # fmt: skip
    moved = shapely.transform(grown, reprojector(grid.crs, meta.crs))
    cell = math.hypot(meta.delta_x, meta.delta_y)
    step("buffer target_grid.py:145", lambda: moved.buffer(cell, join_style="mitre"),
         distance=cell, vertices=int(shapely.get_num_coordinates(moved)))  # fmt: skip
    plan = step("plan_mosaic", lambda: plan_mosaic(footprints, box, needed))
    step("repository.check", lambda: repo.check(plan))
    mosaic = step("assemble (decode + seams)",
                  lambda: assemble(plan, repo.load, needed, load_window=repo.load_window))  # fmt: skip
    threads = os.cpu_count() or 1
    step("resample", lambda: tg.resample(grid, tg.TileWindows(mosaic.tile), threads))
    dem_crs = target
    sources = [("gpkg", EU, "U2018_CLC2018_V2020_20u1", "Code_18")]
    geojson = DATA / f"sweden_corine/{c}_clc2018_3035.geojson"

# features read: the region's buffer, then each source's own reading.
step("buffer feature_input.py:161", lambda: domain.polygon.buffer(fi.MARGIN), distance=fi.MARGIN,
     vertices=int(shapely.get_num_coordinates(domain.polygon)))  # fmt: skip
step("feature_input.source_region(domain, dem, dem)", lambda: fi.source_region(domain, dem_crs, dem_crs))
if geojson is not None:
    doc = step("read_json (GeoJSON window)", lambda: read_json(geojson), mb=geojson.stat().st_size / 1e6)
    rows = step("read_source (GeoJSON window: read_json + shape)",
                lambda: fi.read_source(geojson, None, "Code_18", lambda own: (0, 0, 0, 0)).rows)  # fmt: skip
    print(json.dumps({"catchment": c, "source": "geojson", "rows": len(rows), "vertices": nverts(rows)}))
for kind, path, layer, attr in sources:
    own = "EPSG:3035" if kind == "gpkg" else dem_crs
    box = step(f"feature_input.source_region(domain, dem, {own}) [{kind} box]",
               lambda own=own: fi.source_region(domain, dem_crs, own).bounds)  # fmt: skip

    def query(path=path, layer=layer, attr=attr, box=box):  # type: ignore[no-untyped-def]
        conn = open_geopackage(path)
        try:
            info = layer_info(conn, layer)
            return [(r.fid, r.geometry, r.value) for r in query_features(conn, info, box, attr, 100_000.0)]
        finally:
            conn.close()

    rows = step(f"query_features ({kind}, {path.stat().st_size / 1e9:.2f} GB)", query)
    print(json.dumps({"catchment": c, "source": kind, "rows": len(rows), "vertices": nverts(rows)}))
