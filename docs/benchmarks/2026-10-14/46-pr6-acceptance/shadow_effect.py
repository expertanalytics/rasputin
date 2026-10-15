"""Increment 46 PR 6 acceptance: what the cast shadows take away on the mesh,
on 2027-04-15, 05-15 and 06-21 (section 10.6's dates), for comparison with
section 10.6's figures on DEM nodes.

    python shadow_effect.py <mesh.vtk> <dem path> <out.json>

As `rasputin snow` computes it: one receiver per triangle, the horizon
searched on the near raster to 20 km, the evaluator at the 96 midpoints of
the day's 15-minute intervals (the run's 4 samples per hourly step), the beam
seen from the mesh's centre at its mean vertex height. "off" is the same
receivers under an open sky (self-shading only). Per triangle, the day's mean
direct beam on the slope, beam·b, and the day's mean dI change,
beam·(b_on - b_off) averaged over the 96 samples.

Weights: "plan" is the triangle's map area, the analogue of section 10.6's
evenly spaced DEM nodes and of test E2; "surface" is its inclined area (the
energy on the ground); "triangle" counts each triangle once.
"""

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

from tin_engine.cli import _per_triangle
from tin_engine.crs import site_of
from tin_engine.io.vtk_read import read_mesh
from tin_engine.radiation.evaluate import evaluate_instant
from tin_engine.radiation.frame import receivers_from_mesh
from tin_engine.radiation.geometry import RadiationGeometry
from tin_engine.radiation.request import RadiationRequest
from tin_engine.radiation.search import near_raster, search_horizons
from tin_engine.snow.forcing import sun_and_beam
from tin_engine.snow.geometry import triangle_geometry

mesh, dem, out = sys.argv[1], sys.argv[2], sys.argv[3]
a = read_mesh(Path(mesh).read_bytes())
v, t = a.vertices, a.triangles
geo = triangle_geometry(v, t, a.crs)
ok = ~geo.zero_area
frame = receivers_from_mesh(v, t[ok], a.crs)
box = (*v[:, :2].min(axis=0), *v[:, :2].max(axis=0))
req = RadiationRequest()
near = near_raster((Path(dem),), box, req)
on_g = _per_triangle(search_horizons(frame, near, req, threads=8), ok)
off_g = _per_triangle(RadiationGeometry.open_sky(frame, req), ok)
centre = (v[:, :2].min(axis=0) + v[:, :2].max(axis=0)) / 2
site = site_of(a.crs, float(centre[0]), float(centre[1]))
height = float(np.mean(v[:, 2]))
weights = {"plan": geo.area, "surface": geo.area / np.maximum(geo.normals[:, 2], 1e-12),
           "triangle": np.ones(len(t))}


def wq(x: np.ndarray, w: np.ndarray, q: float) -> float:
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    return float(x[o][np.searchsorted(c, q)])


record = {"mesh": mesh, "triangles": int(len(t)), "site": [site.lat_deg, site.lon_deg, height],
          "dates": {}}
for day in ("2027-04-15", "2027-05-15", "2027-06-21"):
    d0 = int(dt.datetime.fromisoformat(day + "T00:00:00+00:00").timestamp() * 1000)
    times = d0 + (2 * np.arange(96, dtype=np.int64) + 1) * 450_000
    _, beam = sun_and_beam(times, site.lat_deg, site.lon_deg, height)
    e_on, e_off = np.zeros(len(t)), np.zeros(len(t))
    for chunk in np.array_split(np.arange(96), 12):
        f_on = evaluate_instant(on_g, times[chunk], threads=8)
        f_off = evaluate_instant(off_g, times[chunk], threads=8)
        assert (f_on.b <= f_off.b).all()
        e_on += beam[chunk] @ f_on.b
        e_off += beam[chunk] @ f_off.b
    e_on, e_off = e_on / 96, e_off / 96  # day-mean beam on the slope, W/m^2
    change = e_on - e_off  # = the day-mean dI change: the reference cancels
    row = {}
    for name, w in weights.items():
        row[name] = {
            "beam_removed_pct": float(100 * (1 - (w @ e_on) / (w @ e_off))),
            "mean_dI_change_wm2": float(w @ change / w.sum()),
            "p05_dI_change_wm2": wq(change, w, 0.05),
            "losing_more_than_20_wm2_pct": float(100 * w[change < -20].sum() / w.sum()),
        }
    record["dates"][day] = row
    print(day, json.dumps(row["plan"]), flush=True)
Path(out).write_text(json.dumps(record, indent=2) + "\n")
