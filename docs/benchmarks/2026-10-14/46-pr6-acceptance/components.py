"""Increment 46 PR 6 acceptance: the shadow parts of `rasputin snow`, timed
one by one, the way the CLI calls them.

    python components.py <mesh.vtk> <dem path> <days> <day stride> <threads> <out.json>

Times: the near raster's read (`near_raster`), the horizon search
(`search_horizons`), and the evaluator (`evaluate_instant`) on all 24 hourly
steps of every <day stride>-th day of the spring from 2027-04-01 (4 samples
per step), each step timed on its own; the searched geometry is what the
evaluator gets, as in the run.
Also ru_maxrss of this process at the end.
"""

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

from tin_engine.cli import _per_triangle
from tin_engine.crs import site_of  # noqa: F401  (imported as the CLI does)
from tin_engine.io.vtk_read import read_mesh
from tin_engine.radiation.evaluate import evaluate_instant
from tin_engine.radiation.frame import receivers_from_mesh
from tin_engine.radiation.request import RadiationRequest
from tin_engine.radiation.search import near_raster, search_horizons
from tin_engine.snow.forcing import sample_times
from tin_engine.snow.geometry import triangle_geometry
from tin_engine.snow.params import SyntheticSpring

mesh, dem, days, stride, threads, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
arrays = read_mesh(Path(mesh).read_bytes())
v, t = arrays.vertices, arrays.triangles
geometry = triangle_geometry(v, t, arrays.crs)
ok = ~geometry.zero_area
frame = receivers_from_mesh(v, t[ok], arrays.crs)
box = (*v[:, :2].min(axis=0), *v[:, :2].max(axis=0))
request = RadiationRequest()

c = time.perf_counter()
near = near_raster((Path(dem),), box, request)
t_near = time.perf_counter() - c
c = time.perf_counter()
sky = _per_triangle(search_horizons(frame, near, request, threads=threads), ok)
t_search = time.perf_counter() - c

spring = SyntheticSpring(days=days)
steps = []
for k in (24 * d + h for d in range(0, days, stride) for h in range(24)):
    times = sample_times(spring, k)
    c = time.perf_counter()
    evaluate_instant(sky, times, threads=threads)
    steps.append(time.perf_counter() - c)
s = np.array(steps)
record = {
    "mesh": mesh, "dem": dem, "triangles": int(len(t)), "zero_area": int((~ok).sum()),
    "threads": threads, "days": days, "day_stride": stride, "near_raster_shape": list(near.array.shape),
    "near_raster_read_s": t_near, "horizon_search_s": t_search,
    "evaluator_steps": int(len(s)), "evaluator_total_s": float(s.sum()),
    "evaluator_step_median_ms": float(np.median(s) * 1e3),
    "evaluator_step_p05_ms": float(np.percentile(s, 5) * 1e3),
    "evaluator_step_p95_ms": float(np.percentile(s, 95) * 1e3),
    "evaluator_ns_per_triangle_sample": float(s.sum() / (len(s) * 4 * len(t)) * 1e9),
    "ru_maxrss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    "reach_255_share": float((sky.reach[ok] == 255).mean()),
}
Path(out).write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
