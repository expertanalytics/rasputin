"""Triangle density over the whole basin, from random boxes: the extrapolation's sample.

A measurement script, not production code (ROADMAP basin item 2.1); nothing
imports it. Run from the repository root with the venv's python, after
``tools/bench.py``'s Release build exists in ``build-bench/pkg``:

    python docs/benchmarks/2026-10-02/basin-piece-anadem/sample_boxes.py SOURCE BASIN_OUTLINE \
        DATA_DIR OUT_DIR SCRATCH [--n 100] [--side 10000] [--seed 1] \
        [--tolerances 1,2,5,10,20,50]

Box centres are drawn uniformly by area over the basin (rejection sampling in
an Albers equal-area projection on GRS80, seed fixed), and a box is kept only
if its square lies wholly inside the basin, so the sample under-represents a
band of half a box along the basin's edge. Each box is a ``--side`` metre
square in the SIRGAS 2000 / UTM zone of its centre (EPSG:31960 + zone); its
DEM is fetched, projected and resampled to 30 m by ``prep_dem.py`` exactly as
the piece's was, and it is meshed with ``--domain`` = the square at each
tolerance, once (timing is not the point here), with the final-check
measurement of ``run_sweep.py``. Per-box results go to ``OUT_DIR/boxes.json``.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import run_sweep  # noqa: E402  (puts the Release package first on sys.path)
import prep_dem  # noqa: E402
import pyproj  # noqa: E402
import shapely  # noqa: E402
import shapely.ops  # noqa: E402

AEA = "+proj=aea +lat_1=-10 +lat_2=-18.5 +lat_0=-14 +lon_0=-42 +ellps=GRS80 +units=m +no_defs"


def draw(basin: Path, n: int, side: float, seed: int) -> list[tuple[float, float]]:
    poly, crs = prep_dem._outline(basin)
    to_aea = pyproj.Transformer.from_crs(crs, AEA, always_xy=True)
    to_ll = pyproj.Transformer.from_crs(AEA, "EPSG:4326", always_xy=True)
    p = shapely.ops.transform(to_aea.transform, poly)
    shapely.prepare(p)
    rng = np.random.default_rng(seed)
    x0, y0, x1, y1 = p.bounds
    out, tried = [], 0
    while len(out) < n:
        x, y = rng.uniform(x0, x1), rng.uniform(y0, y1)
        tried += 1
        if p.contains(shapely.box(x - side / 2, y - side / 2, x + side / 2, y + side / 2)):
            out.append(to_ll.transform(x, y))
    print(f"{n} boxes from {tried} draws", file=sys.stderr)
    return out


def box_outline(k: int, lon: float, lat: float, side: float, folder: Path) -> tuple[Path, int]:
    zone = math.floor((lon + 180.0) / 6.0) + 1
    epsg = 31960 + zone  # SIRGAS 2000 / UTM zone NS
    x, y = pyproj.Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True).transform(lon, lat)
    x, y = round(x / 30.0) * 30.0, round(y / 30.0) * 30.0
    h = side / 2
    ring = [[x - h, y - h], [x + h, y - h], [x + h, y + h], [x - h, y + h], [x - h, y - h]]
    path = folder / f"box{k:03d}_outline_epsg{epsg}.geojson"
    path.write_text(json.dumps({
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": f"urn:ogc:def:crs:EPSG::{epsg}"}},
        "features": [{"type": "Feature", "properties": {"box": k, "lon": lon, "lat": lat},
                      "geometry": {"type": "Polygon", "coordinates": [ring]}}],
    }))  # fmt: skip
    return path, epsg


def main() -> None:
    source, basin, data, out, scratch = sys.argv[1:6]
    rest = sys.argv[6:]
    opt = lambda k, d: rest[rest.index(k) + 1] if k in rest else d  # noqa: E731
    n, side, seed = int(opt("--n", 100)), float(opt("--side", 10000)), int(opt("--seed", 1))
    tols = [float(t) for t in opt("--tolerances", "1,2,5,10,20,50").split(",")]
    folder = Path(data) / "boxes"
    folder.mkdir(parents=True, exist_ok=True)
    out_dir, scr = Path(out), Path(scratch)
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)
    scr.mkdir(parents=True, exist_ok=True)
    centres = draw(Path(basin), n, side, seed)
    path = out_dir / "boxes.json"
    rows = json.loads(path.read_text())["boxes"] if path.exists() else []
    done = {r["box"] for r in rows}
    for k, (lon, lat) in enumerate(centres):
        if k in done:
            continue
        outline, epsg = box_outline(k, lon, lat, side, folder)
        win, grid = prep_dem._names(source, outline, folder, epsg, 30.0)
        if not grid.exists():
            prep_dem.fetch(source, outline, folder, epsg=epsg)
            prep_dem.resample(source, outline, folder, epsg=epsg)
        z = prep_dem.read_geotiff(grid)[0]
        rec = {"box": k, "lon": lon, "lat": lat, "epsg": epsg, "z_min": float(z[z != -9999].min()),
               "z_max": float(z[z != -9999].max()), "tolerances": {}}  # fmt: skip
        for tol in tols:
            vtk = scr / f"{out_dir.name}_box{k:03d}_t{tol:g}.ascii.vtk"  # unique per sample
            r = run_sweep._run(["mesh", "--dem", str(grid), "--domain", str(outline), "--tolerance",
                                f"{tol:g}", "--ascii", "--out", str(vtk), "--stats",
                                str(out_dir / "logs" / f"box{k:03d}_t{tol:g}.stats.md")],
                               out_dir / "logs" / f"box{k:03d}_t{tol:g}.log")  # fmt: skip
            r.update(run_sweep._sizes(out_dir / "logs" / f"box{k:03d}_t{tol:g}.stats.md"))
            mesh = run_sweep.bench.read_vtk_ascii(vtk)
            fc = run_sweep.final_check(mesh, grid, win, outline, tol)
            r["final_check"] = fc
            if fc["control_resampled_grid"]["over_tolerance"]:
                raise SystemExit(f"box {k} t {tol}: the control failed: {fc}")
            rec["tolerances"][f"{tol:g}"] = r
            vtk.unlink()
        rows.append(rec)
        path.write_text(json.dumps({"source": source, "basin": basin, "n": n, "side_m": side,
                                    "seed": seed, "updated_utc": time.strftime(
                                        "%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                    "boxes": rows}, separators=(",", ":")))  # fmt: skip
        t = rec["tolerances"]
        print(f"box {k:3d} ({lon:.3f}, {lat:.3f}) z {rec['z_min']:.0f}-{rec['z_max']:.0f}: "
              + ", ".join(f"t{k2} {v['triangles']}" for k2, v in t.items()), flush=True)


if __name__ == "__main__":
    main()
