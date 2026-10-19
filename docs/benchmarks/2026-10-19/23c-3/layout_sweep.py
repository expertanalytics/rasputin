"""The bench's child at many process layouts: each sample runs `bench.py _child`
with an environment variable PAD of a different length (0 to 8 * (N - 1)
bytes), which moves the stack and the heap's start. Trees interleaved per
pad, domains tile and quarter, the given thread count. One TSV line per
sample to stdout: tree, domain, threads, pad, refine_s, app_s, proc_s (@perf).
python layout_sweep.py W N THREADS name=TREE ..."""
import json
import os
import subprocess
import sys
import time

w, n, threads = sys.argv[1], int(sys.argv[2]), sys.argv[3]
trees = [a.split("=", 1) for a in sys.argv[4:]]
dem = f"{w}/tests/fixtures/dem_archive/7908_3_10m_z33.tif"
quarter = f"{w}/docs/benchmarks/2026-09-26/quarter.geojson"
out = os.environ.get("SWEEP_OUT", "/tmp/sweep.vtk")
for pad in range(0, 8 * n, 8):
    for domain, extra in (("tile", []), ("quarter", ["--domain", quarter])):
        for name, tree in trees:
            argv = [f"{w}/.venv/bin/python", f"{w}/tools/bench.py", "_child", "--pkg",
                    f"{tree}/build-bench/pkg", "--threads", threads, "--", "mesh", "--dem", dem,
                    "--tolerance", "1", *extra, "--out", out, "--binary"]  # fmt: skip
            env = dict(os.environ, PAD="x" * pad)
            t0 = time.perf_counter()
            done = subprocess.run(argv, env=env, capture_output=True, text=True, cwd=w)
            proc = time.perf_counter() - t0
            line = next(x for x in done.stderr.splitlines() if x.startswith("BENCH "))
            d = json.loads(line[6:])
            print(f"{name}\t{domain}\t{threads}\t{pad}\t{d['refine_s']:.6f}\t{d['app_s']:.6f}\t{proc:.6f}",
                  flush=True)  # fmt: skip
