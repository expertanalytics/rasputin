"""One whole-basin ``rasputin mesh`` run, timed, with its memory watched.

A measurement script, not production code; nothing imports it.

    python run_basin.py RASPUTIN OUT_JSON LOG_DIR TAG -- MESH_ARGS...

Runs ``/usr/bin/time -l RASPUTIN MESH_ARGS`` and, every 5 s, samples the
child's RSS (``ps``) and the machine's swap in use (``sysctl vm.swapusage``)
into ``LOG_DIR/TAG.mem.csv``. If swap in use grows past its value at the start
by more than ``SWAP_GROWTH_LIMIT_GB`` (env, default 16) the run is killed and recorded as stopped (the machine has 32 GiB; a run deep in swap
measures the disk, and could wedge the session). ``pmset -g batt`` is recorded
before and after. The record (wall, exit status, ``time -l``'s maximum
resident set size and peak memory footprint, the memory trace's maxima) is
appended to OUT_JSON under ``runs``.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

SWAP_GROWTH_LIMIT_GB = float(os.environ.get("SWAP_GROWTH_LIMIT_GB", "16"))


def _swap_gb() -> float:
    out = subprocess.run(["sysctl", "-n", "vm.swapusage"], capture_output=True, text=True).stdout
    m = re.search(r"used = ([\d.]+)M", out)
    return float(m.group(1)) / 1024 if m else 0.0


def _rss_gb(pid: int) -> float:
    """RSS of `pid`'s children (``time`` runs the CLI as its child)."""
    kids = subprocess.run(["pgrep", "-P", str(pid)], capture_output=True, text=True).stdout.split()
    total = 0
    for k in kids:
        out = subprocess.run(["ps", "-o", "rss=", "-p", k], capture_output=True, text=True).stdout
        total += int(out.strip() or 0)
    return total / 2**20


def _pmset() -> str:
    return subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout.strip()


def main() -> None:
    rasputin, out_json, log_dir, tag = sys.argv[1:5]
    args = sys.argv[sys.argv.index("--") + 1 :]
    logs = Path(log_dir)
    logs.mkdir(parents=True, exist_ok=True)
    log, mem = logs / f"{tag}.log", logs / f"{tag}.mem.csv"
    rec: dict = {"tag": tag, "args": args, "pmset_before": _pmset(),
                 "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                 "swap_gb_before": _swap_gb()}  # fmt: skip
    t0 = time.time()
    with log.open("w") as fh, mem.open("w") as mf:
        mf.write("t_s,rss_gb,swap_gb\n")
        p = subprocess.Popen(["/usr/bin/time", "-l", rasputin, *args], stdout=fh, stderr=subprocess.STDOUT)
        max_rss = max_swap = 0.0
        stopped = None
        while p.poll() is None:
            time.sleep(5)
            rss, swap = _rss_gb(p.pid), _swap_gb()
            max_rss, max_swap = max(max_rss, rss), max(max_swap, swap)
            mf.write(f"{time.time() - t0:.0f},{rss:.3f},{swap:.3f}\n")
            mf.flush()
            if swap - rec["swap_gb_before"] > SWAP_GROWTH_LIMIT_GB:
                stopped = f"swap grew to {swap:.1f} GB (+{SWAP_GROWTH_LIMIT_GB} GB limit) at {time.time() - t0:.0f} s"
                subprocess.run(["pkill", "-TERM", "-P", str(p.pid)])
                p.wait()
                break
    rec.update({"wall_s": round(time.time() - t0, 2), "returncode": p.returncode,
                "stopped": stopped, "sampled_max_rss_gb": round(max_rss, 3),
                "sampled_max_swap_gb": round(max_swap, 3), "pmset_after": _pmset()})  # fmt: skip
    text = log.read_text(errors="replace")
    for key, pat in (("max_rss_bytes", r"(\d+)\s+maximum resident set size"),
                     ("peak_footprint_bytes", r"(\d+)\s+peak memory footprint"),
                     ("time_real_s", r"([\d.]+) real")):  # fmt: skip
        m = re.search(pat, text)
        rec[key] = float(m.group(1)) if m else None
    path = Path(out_json)
    doc = json.loads(path.read_text()) if path.exists() else {"runs": []}
    doc["runs"].append(rec)
    path.write_text(json.dumps(doc, indent=1))
    print(json.dumps({k: rec[k] for k in ("tag", "wall_s", "returncode", "stopped", "max_rss_bytes",
                                          "peak_footprint_bytes", "sampled_max_swap_gb")}))  # fmt: skip


if __name__ == "__main__":
    main()
