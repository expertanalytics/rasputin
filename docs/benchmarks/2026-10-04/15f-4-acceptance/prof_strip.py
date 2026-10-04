"""Where refine_strip's untimed time goes (@perf, 2026-10-04): a driver for
macOS `sample`. It runs `rasputin mesh` from a package directory, here the
RelWithDebInfo build at -O3 -g in build-prof/pkg. edge_strip.run is wrapped so
that, after the real call, `_core.refine_strip` runs K more times on the same
inputs, then K times with an empty strip (no points). Each call's wall time
and its two timed phases (scan_seconds, split_seconds) are printed as JSON on
stderr. The empty-strip calls give the per-call fixed cost.

python prof_strip.py PKG K -- <rasputin mesh argv>
"""

import json
import sys
import time

split = sys.argv.index("--")
pkg, k = sys.argv[1], int(sys.argv[2])
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, pkg)

import numpy as np  # noqa: E402

import tin_engine._core as core  # noqa: E402
import tin_engine.cli as cli  # noqa: E402
import tin_engine.edge_strip as edge_strip  # noqa: E402

assert core.__file__.startswith(pkg), core.__file__
real_run = edge_strip.run


def run(view, strip, start, tolerance, clock):
    out = real_run(view, strip, start, tolerance, clock)
    arrays = [np.asarray(a) for a in (start.vertices, start.triangles, start.z, start.valid, start.edges, start.masks)]
    empty = core.constraint_check_points(view, arrays[0], arrays[4][:0])
    rows = []
    for name, s in (("full", strip), ("empty", empty)):
        for _ in range(k):
            t0 = time.perf_counter()
            o = core.refine_strip(view, s, *arrays, tolerance=tolerance)
            rows.append({"strip": name, "points": s.size, "wall_s": time.perf_counter() - t0,
                         "scan_s": o.scan_seconds, "split_s": o.split_seconds,
                         "inserted": o.strip_inserted, "nodes_inserted": o.nodes_inserted})  # fmt: skip
            del o
    print("PROF " + json.dumps({"vertices": int(len(arrays[0])), "triangles": int(len(arrays[1])),
                                "edges": int(len(arrays[4])), "calls": rows}), file=sys.stderr)  # fmt: skip
    return out


edge_strip.run = run
sys.argv = ["rasputin", *sys.argv[split + 1 :]]
cli.app()
