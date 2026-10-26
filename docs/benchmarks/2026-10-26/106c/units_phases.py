"""106c acceptance: the `catchments` phase split into reading and selecting,
repairing, and simplifying, with the calls `cli._units` makes, one run each.
Run from the worktree's venv:
    .venv/bin/python docs/benchmarks/2026-10-26/106c/units_phases.py \
        UNITS.gpkg DOMAIN.geojson TRIANGLES
"""

import sys
import time
from pathlib import Path

import shapely

from tin_engine import sub_catchments
from tin_engine.feature_input import read_source
from tin_engine.io.domain_file import read_domain
from tin_engine.network import join

units_path, domain_path, triangles = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
t0 = time.perf_counter()
given = read_domain(domain_path)
read = read_source(units_path, None, None, lambda own: given.to_crs(own).polygon.bounds)
rows = [(f, shapely.force_2d(g)) for f, g, _ in read.rows if g is not None and not g.is_empty]
inside = given.to_crs(read.crs).polygon
picked = [rows[k] for k in sub_catchments.select([g for _, g in rows], inside)]
t1 = time.perf_counter()
repaired = sub_catchments.repair([g for _, g in picked])
area = join(repaired, read.crs).polygon.area
band = sub_catchments.band_m(area, triangles)
t2 = time.perf_counter()
kept = sub_catchments.simplify(repaired, band)
t3 = time.perf_counter()
print(
    f"{domain_path.stem}: {len(rows)} read, {len(picked)} selected, band {band:.3f} m; "
    f"read and select {t1 - t0:.3f} s, repair {t2 - t1:.3f} s, simplify {t3 - t2:.3f} s"
)
