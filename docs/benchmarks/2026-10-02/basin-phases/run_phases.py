"""One ``rasputin mesh`` run through `phase_driver.py`, its memory sampled every
0.25 s and cut by phase. A measurement script, not production code.

    python run_phases.py OUT_DIR TAG -- MESH_ARGS...

Runs ``/usr/bin/time -l PYTHON phase_driver.py OUT_DIR/TAG.markers.jsonl --
MESH_ARGS`` and samples the CLI process (``time``'s child) with
``proc_pid_rusage(RUSAGE_INFO_V4)``: resident size, ``phys_footprint`` (what
``time -l`` reports as peak memory footprint, compressed and swapped pages
included) and the lifetime maximum footprint, into ``TAG.mem.csv``. Swap in
use is sampled too; if it grows by more than 8 GB the run is killed and
recorded as stopped (this run is meant to fit). ``pmset -g batt`` before and
after. Writes ``TAG.json``: the record, and per marker span the peak
footprint and resident size of the samples inside it.
"""

from __future__ import annotations

import ctypes
import json
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SWAP_GROWTH_LIMIT_GB = 8.0
_libc = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
_buf = ctypes.create_string_buffer(512)


def _usage(pid: int) -> tuple[int, int, int] | None:
    """(resident, phys_footprint, lifetime max phys_footprint) bytes, or None.
    Offsets in ``struct rusage_info_v4`` (``<sys/resource.h>``)."""
    if _libc.proc_pid_rusage(pid, 4, _buf) != 0:
        return None
    q = lambda off: int.from_bytes(_buf.raw[off : off + 8], "little")  # noqa: E731
    return q(64), q(72), q(240)


def _swap_gb() -> float:
    out = subprocess.run(["sysctl", "-n", "vm.swapusage"], capture_output=True, text=True).stdout
    m = re.search(r"used = ([\d.]+)M", out)
    return float(m.group(1)) / 1024 if m else 0.0


def _sh(*argv: str) -> str:
    return subprocess.run(argv, capture_output=True, text=True).stdout.strip()


def spans(markers: list[dict], samples: list[tuple[float, int, int]]) -> list[dict]:
    """Each B/E pair of a name (in order) with the peak footprint and resident
    size over the samples whose time falls inside it; and, from the markers'
    own readings, the exact peak when the lifetime peak footprint rose inside
    the span (a sample every 0.25 s can miss a short transient)."""
    open_: dict[str, list[dict]] = {}
    out = []
    for m in markers:
        if m["ev"] == "B":
            open_.setdefault(m["name"], []).append(m)
        else:
            b = open_[m["name"]].pop()
            inside = [s for s in samples if b["t"] <= s[0] <= m["t"]]
            out.append({"name": m["name"], "t0": round(b["t"], 3), "t1": round(m["t"], 3),
                        "seconds": round(m["t"] - b["t"], 3), "samples": len(inside),
                        "peak_footprint": max((s[2] for s in inside), default=None),
                        "peak_resident": max((s[1] for s in inside), default=None),
                        "fp_end": m.get("fp"), "life_max_end": m.get("life_max"),
                        **{f"end_{k}": v for k, v in m.items() if k not in ("t", "ev", "name", "fp", "life_max")},
                        # the exact peak over the span when the lifetime peak rose in it
                        "exact_peak": m["life_max"] if m.get("life_max", 0) > b.get("life_max", 0) else None,
                        **{k: v for k, v in b.items() if k not in ("t", "ev", "name")}})  # fmt: skip
    return sorted(out, key=lambda s: s["t0"])


def main() -> None:
    out_dir, tag = Path(sys.argv[1]), sys.argv[2]
    args = sys.argv[sys.argv.index("--") + 1 :]
    out_dir.mkdir(parents=True, exist_ok=True)
    markers_path, mem_path = out_dir / f"{tag}.markers.jsonl", out_dir / f"{tag}.mem.csv"
    rec: dict = {"tag": tag, "args": args, "pmset_before": _sh("pmset", "-g", "batt"),
                 "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                 "swap_gb_before": _swap_gb()}  # fmt: skip
    argv = ["/usr/bin/time", "-l", sys.executable, str(HERE / "phase_driver.py"), str(markers_path), "--", *args]
    log = (out_dir / f"{tag}.log").open("w")
    t0 = time.time()
    p = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT)
    child = None
    while child is None and p.poll() is None:
        kids = _sh("pgrep", "-P", str(p.pid)).split()
        child = int(kids[0]) if kids else None
        time.sleep(0.02)
    samples: list[tuple[float, int, int]] = []
    stopped, lifetime_max, n = None, 0, 0
    with mem_path.open("w") as mf:
        mf.write("t,resident_bytes,footprint_bytes,swap_gb\n")
        while p.poll() is None:
            u = _usage(child) if child else None
            if u:
                t = time.time()
                samples.append((t, u[0], u[1]))
                lifetime_max = max(lifetime_max, u[2])
                swap = _swap_gb() if n % 4 == 0 else None
                mf.write(f"{t:.3f},{u[0]},{u[1]},{'' if swap is None else f'{swap:.3f}'}\n")
                n += 1
                if swap is not None and swap - rec["swap_gb_before"] > SWAP_GROWTH_LIMIT_GB:
                    stopped = f"swap grew to {swap:.1f} GB at {t - t0:.0f} s"
                    subprocess.run(["kill", "-TERM", str(child)])
                    p.wait()
                    break
            time.sleep(0.25)
    log.close()
    text = (out_dir / f"{tag}.log").read_text(errors="replace")
    for key, pat in (("time_max_rss_bytes", r"(\d+)\s+maximum resident set size"),
                     ("time_peak_footprint_bytes", r"(\d+)\s+peak memory footprint"),
                     ("time_real_s", r"([\d.]+) real")):  # fmt: skip
        m = re.search(pat, text)
        rec[key] = float(m.group(1)) if m else None
    markers = [json.loads(line) for line in markers_path.read_text().splitlines()]
    rec.update({"wall_s": round(time.time() - t0, 2), "returncode": p.returncode, "stopped": stopped,
                "sampled_lifetime_max_footprint": lifetime_max, "n_samples": len(samples),
                "pmset_after": _sh("pmset", "-g", "batt"), "spans": spans(markers, samples)})  # fmt: skip
    (out_dir / f"{tag}.json").write_text(json.dumps(rec, indent=1))
    for s in rec["spans"]:
        pf = s["peak_footprint"]
        ex = s["exact_peak"]
        print(f"{s['name']:40s} {s['seconds']:8.2f} s  sampled {pf / 1e9 if pf else float('nan'):6.2f} GB"
              f"  exact {ex / 1e9 if ex else float('nan'):6.2f} GB")
    print({k: rec[k] for k in ("tag", "wall_s", "returncode", "stopped", "time_peak_footprint_bytes")})


if __name__ == "__main__":
    main()
