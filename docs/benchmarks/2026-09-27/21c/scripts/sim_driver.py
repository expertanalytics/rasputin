"""sim_driver.py PKG MODE DOMAIN OUTDIR [--foot] [--seed N] [--keep]

One simulated refine of the 1 m benchmark (tolerance 1, CLI defaults) with the
scratch build PKG (scripts/sim.patch applied), RASPUTIN_SIM=MODE, through
tools/bench.py's child, writing an ASCII VTK. Then bench.py's quality() on that
mesh (worst angle, max degree, tolerance, constrained Delaunay check), merged
with the simulation's own JSON, is written to OUTDIR/<domain>_<mode>.json.
DOMAIN is `tile` or a path to a GeoJSON polygon.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
DEM = REPO / "tests/fixtures/dem_archive/7908_3_10m_z33.tif"

pkg, mode, domain, outdir = sys.argv[1:5]
foot = "--foot" in sys.argv
seed = sys.argv[sys.argv.index("--seed") + 1] if "--seed" in sys.argv else None
out = Path(outdir)
out.mkdir(parents=True, exist_ok=True)
name = "tile" if domain == "tile" else Path(domain).stem
stem = f"{name}_{mode}{'_foot' if foot else ''}{'_seed' + seed if seed else ''}"
vtk = out / f"{stem}.vtk"
env = dict(os.environ, RASPUTIN_SIM=mode, RASPUTIN_SIM_OUT=str(out / f"{stem}.sim.json"))
if foot:
    env["RASPUTIN_SIM_FOOT"] = "1"
if seed:
    env["RASPUTIN_SIM_SEED"] = seed
argv = [sys.executable, str(REPO / "tools/bench.py"), "_child", "--pkg", pkg, "--threads", "0",
        "--", "mesh", "--dem", str(DEM), "--tolerance", "1",
        *([] if domain == "tile" else ["--domain", domain]), "--out", str(vtk), "--ascii"]
t0 = time.perf_counter()
done = subprocess.run(argv, env=env, capture_output=True, text=True)
wall = time.perf_counter() - t0
if done.returncode != 0:
    sys.exit(done.stderr[-3000:])
child = json.loads([l for l in done.stderr.splitlines() if l.startswith("BENCH ")][-1][6:])
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, pkg)
sys.path.insert(0, str(REPO / "tools"))
import bench  # noqa: E402

m = bench.read_vtk_ascii(vtk)
q = bench.quality(m.points, m.triangles, m.edges, 1.0, child["max_error"])
sim = json.loads((out / f"{stem}.sim.json").read_text())
rec = dict(domain=name, mode=mode, seed=seed, foot=foot, wall_s=wall, child=child, quality=q.model_dump(),
           vtk_sha256=m.sha256, vtk_triangles=int(len(m.triangles)), sim=sim)
(out / f"{stem}.json").write_text(json.dumps(rec, indent=1))
if "--keep" not in sys.argv:
    vtk.unlink()
qq = rec["quality"]
print(f"{stem}: tris {len(m.triangles)} rounds {child['rounds']} ins {child['inserted']} flips "
      f"{child['flips']} maxerr {child['max_error']:.4f} full {sim['full_rescan_max_error']:.4f}/"
      f"{sim['full_rescan_needs_split']} worst {qq['worst_angle']:.4f} deg {qq['max_degree']} "
      f"tol {qq['within_tolerance']} delaunay {qq['delaunay_violations']}/{qq['delaunay_checked']}")
