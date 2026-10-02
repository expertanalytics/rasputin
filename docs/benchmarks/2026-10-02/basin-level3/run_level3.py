"""Mesh each BHO level-3 unit of the São Francisco separately. A measurement script.

    python run_level3.py RUNS_DIR MESH_DIR UNITS_JSON OUTLINE_DIR CRS_WKT CACHE \
        TOLERANCE [--binary=P,P...] [--estimates=JSON] [--driver=PY]
        [--swap-kill-mb=MB] PREFIX...

Per unit, in order: estimate the peak footprint from the basin-phases figures
and skip the unit if it exceeds ``LIMIT_GB``; ``rasputin fetch`` the unit's
blocks; ``MallocLargeCache=0 /usr/bin/time -l rasputin mesh`` it into
``MESH_DIR``, as binary VTK for the units named in ``--binary`` and as text
for the rest. Writes ``RUNS_DIR/<prefix>.{fetch.out,log,stats.md,json,mem.csv}``.

Swap is watched from outside: if it grows by more than ``--swap-kill-mb``
(default 3000) over its level at the start, the run is killed and no further
unit is started. The mesh process's footprint (``phys_footprint``, what
``time -l`` reports) is sampled every 0.25 s into ``<prefix>.mem.csv``.

``--estimates`` names a JSON object ``{prefix: GB}`` that replaces the
basin-phases estimate (the 2 m runs, whose mesh may set the peak).
``--driver`` runs basin-phases' ``phase_driver.py`` (the CLI in-process, with
phase markers into ``<prefix>.markers.jsonl``) instead of the bare CLI; the
record then holds ``late_peak_gb``, the largest sampled footprint from the
start of refine phase 1 to the end of the run.

``rasputin`` is the CLI of the interpreter running this script
(``python -c "from tin_engine.cli import app; app()"``).
"""

import ctypes
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


_libc = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
_buf = ctypes.create_string_buffer(512)


def footprint(pid: int) -> int | None:
    """``phys_footprint`` of ``pid`` (``struct rusage_info_v4``, offset 72)."""
    if _libc.proc_pid_rusage(pid, 4, _buf) != 0:
        return None
    return int.from_bytes(_buf.raw[72:80], "little")


def child_of(pid: int) -> int | None:
    out = subprocess.run(["pgrep", "-P", str(pid)], capture_output=True, text=True).stdout
    return int(out.split()[0]) if out.split() else None


def late_peak(markers: Path, samples: list[tuple[float, int]]) -> float | None:
    """The largest sample from refine phase 1's start on, GB."""
    if not markers.exists():
        return None
    lines = [json.loads(x) for x in markers.read_text().splitlines()]
    t = next((m["t"] for m in lines if m["name"] == "refine phase 1" and m["ev"] == "B"), None)
    late = [f for (ts, f) in samples if t is not None and ts >= t]
    return max(late) / 1e9 if late else None


def power() -> str:
    return subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout


def main() -> None:
    runs, meshes, units_json, outlines, crs_file, cache = (Path(a) for a in sys.argv[1:7])
    tol, rest = sys.argv[7], sys.argv[8:]
    opt = dict(a[2:].split("=", 1) for a in rest if a.startswith("--"))
    binary = set(opt.get("binary", "").split(","))
    given = json.loads(Path(opt["estimates"]).read_text()) if "estimates" in opt else {}
    swap_kill = float(opt.get("swap-kill-mb", SWAP_KILL_MB))
    prefixes = [a for a in rest if not a.startswith("--")]
    crs = crs_file.read_text().strip()
    units = {u["prefix"]: u for u in json.loads(units_json.read_text())}
    cli = [sys.executable, "-c", "from tin_engine.cli import app; app()"]
    stop = False
    runs.mkdir(parents=True, exist_ok=True)
    meshes.mkdir(parents=True, exist_ok=True)
    for p in prefixes:
        if stop:
            break
        unit = units[p]
        est = given.get(p, estimate_gb(unit))
        rec: dict = {"prefix": p, "estimate_gb": round(est, 2)}
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
        markers = runs / f"{p}.markers.jsonl"
        run = [sys.executable, opt["driver"], str(markers), "--"] if "driver" in opt else cli
        cmd = ["/usr/bin/time", "-l", *run, "mesh", "--dem", "anadem-v1", "--cache", str(cache),
               "--domain", str(domain), "--out-crs", crs, "--tolerance", tol,
               "--out", str(vtk), "--stats", str(stats),
               "--binary" if p in binary else "--ascii"]  # fmt: skip
        rec["power_before"], rec["swap_before_mb"] = power(), swap_used_mb()
        env = {**os.environ, "MallocLargeCache": "0"}
        t0 = time.monotonic()
        samples: list[tuple[float, int]] = []
        with open(runs / f"{p}.log", "w") as log, open(runs / f"{p}.mem.csv", "w") as mem:
            proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
            max_swap, child, n = rec["swap_before_mb"], None, 0
            mem.write("t,footprint_bytes,swap_mb\n")
            while proc.poll() is None:
                child = child or child_of(proc.pid)
                f = footprint(child) if child else None
                swap = swap_used_mb() if n % 2 == 0 else None
                n += 1
                if f is not None:
                    samples.append((time.time(), f))
                    mem.write(f"{time.time():.3f},{f},{'' if swap is None else swap}\n")
                max_swap = max(max_swap, swap or 0.0)
                if max_swap - rec["swap_before_mb"] > swap_kill:
                    for pid in (child, proc.pid):
                        if pid:
                            subprocess.run(["kill", "-KILL", str(pid)])
                    rec["killed"] = f"swap +{max_swap - rec['swap_before_mb']:.0f} MB"
                    stop = True
                    break
                time.sleep(0.25)
            proc.wait()
        rec["sampled_peak_gb"] = max((f for _, f in samples), default=0) / 1e9
        rec["late_peak_gb"] = late_peak(markers, samples)
        rec["exit"], rec["outer_wall_s"] = proc.returncode, round(time.monotonic() - t0, 1)
        rec["swap_max_mb"], rec["power_after"] = max_swap, power()
        rec["mesh"] = str(vtk)
        (runs / f"{p}.json").write_text(json.dumps(rec, indent=1))
        print(p, "exit", proc.returncode, rec["outer_wall_s"], "s", flush=True)


if __name__ == "__main__":
    main()
