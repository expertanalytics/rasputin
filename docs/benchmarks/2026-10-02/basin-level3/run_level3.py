"""Mesh each BHO level-3 unit of the São Francisco separately. A measurement script.

    python run_level3.py RUNS_DIR MESH_DIR UNITS_JSON OUTLINE_DIR CRS_WKT CACHE \
        TOLERANCE [--binary=P,P...] PREFIX...

Per unit, in order: estimate the peak footprint from the basin-phases figures
and skip the unit if it exceeds ``LIMIT_GB``; ``rasputin fetch`` the unit's
blocks; ``MallocLargeCache=0 /usr/bin/time -l rasputin mesh`` it into
``MESH_DIR``, watching swap from outside (kill at +3 GB over its level at the
start), as binary VTK for the units named in ``--binary`` and as text for the
rest. Writes ``RUNS_DIR/<prefix>.{fetch.out,log,stats.md,json}``.

``rasputin`` is the CLI of the interpreter running this script
(``python -c "from tin_engine.cli import app; app()"``).
"""

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

LIMIT_GB = 24.0
SWAP_KILL_MB = 3000.0
THREADS = 10

# basin-phases (docs/benchmarks/2026-10-02/basin-phases/README.md), measured
# with MallocLargeCache=0 on unit 761: bytes per source node while loading,
# per source/target node held, the resample transient per row x col x thread
# over 256 rows, and the store charged per check point (0.5 per box node).
LOAD_B, HELD_B, RESAMPLE_B, STORE_B, POINTS_PER_NODE, SRC_PER_NODE = (
    12.0, 4.0, 135.0, 30.0, 0.5, 1.04)  # fmt: skip


def estimate_gb(unit: dict) -> float:
    """The largest of the load, resample and store phases, per basin-phases."""
    n, cols = unit["grid_nodes"], unit["grid_cols"]
    src = SRC_PER_NODE * n
    load = LOAD_B * src
    resample = HELD_B * (src + n) + RESAMPLE_B * 256 * cols * THREADS
    store = HELD_B * (src + n) + STORE_B * POINTS_PER_NODE * n
    return max(load, resample, store) / 1e9


def swap_used_mb() -> float:
    out = subprocess.run(["sysctl", "-n", "vm.swapusage"], capture_output=True, text=True).stdout
    return float(re.search(r"used = ([\d.]+)M", out).group(1))


def power() -> str:
    return subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout


def main() -> None:
    runs, meshes, units_json, outlines, crs_file, cache = (Path(a) for a in sys.argv[1:7])
    tol, rest = sys.argv[7], sys.argv[8:]
    binary = {u for a in rest if a.startswith("--binary=") for u in a[9:].split(",")}
    prefixes = [a for a in rest if not a.startswith("--")]
    crs = crs_file.read_text().strip()
    units = {u["prefix"]: u for u in json.loads(units_json.read_text())}
    cli = [sys.executable, "-c", "from tin_engine.cli import app; app()"]
    runs.mkdir(parents=True, exist_ok=True)
    meshes.mkdir(parents=True, exist_ok=True)
    for p in prefixes:
        unit = units[p]
        rec: dict = {"prefix": p, "estimate_gb": round(estimate_gb(unit), 2)}
        domain = outlines / unit["file"]
        if rec["estimate_gb"] > LIMIT_GB:
            rec["skipped"] = f"estimate {rec['estimate_gb']} GB > {LIMIT_GB} GB"
            (runs / f"{p}.json").write_text(json.dumps(rec, indent=1))
            print(p, rec["skipped"], flush=True)
            continue
        fetch = subprocess.run([*cli, "fetch", "anadem-v1", "--cache", str(cache),
                                "--domain", str(domain), "--out-crs", crs],
                               capture_output=True, text=True)  # fmt: skip
        (runs / f"{p}.fetch.out").write_text(fetch.stdout + fetch.stderr)
        rec["fetch_exit"] = fetch.returncode
        if fetch.returncode:
            (runs / f"{p}.json").write_text(json.dumps(rec, indent=1))
            print(p, "fetch failed", flush=True)
            continue
        vtk = meshes / f"sub_basin_{p}_anadem_tol{tol}m.vtk"
        stats = runs / f"{p}.stats.md"
        cmd = ["/usr/bin/time", "-l", *cli, "mesh", "--dem", "anadem-v1", "--cache", str(cache),
               "--domain", str(domain), "--out-crs", crs, "--tolerance", tol,
               "--out", str(vtk), "--stats", str(stats),
               "--binary" if p in binary else "--ascii"]  # fmt: skip
        rec["power_before"], rec["swap_before_mb"] = power(), swap_used_mb()
        env = {**os.environ, "MallocLargeCache": "0"}
        t0 = time.monotonic()
        with open(runs / f"{p}.log", "w") as log:
            proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
            max_swap = rec["swap_before_mb"]
            while proc.poll() is None:
                max_swap = max(max_swap, swap_used_mb())
                if max_swap - rec["swap_before_mb"] > SWAP_KILL_MB:
                    proc.kill()
                    rec["killed"] = f"swap +{max_swap - rec['swap_before_mb']:.0f} MB"
                    break
                time.sleep(1.0)
            proc.wait()
        rec["exit"], rec["outer_wall_s"] = proc.returncode, round(time.monotonic() - t0, 1)
        rec["swap_max_mb"], rec["power_after"] = max_swap, power()
        rec["mesh"] = str(vtk)
        (runs / f"{p}.json").write_text(json.dumps(rec, indent=1))
        print(p, "exit", proc.returncode, rec["outer_wall_s"], "s", flush=True)


if __name__ == "__main__":
    main()
