"""Run `rasputin mesh` from a copy of tin_engine, optionally traced or patched.

The copy is `$PKG/tin_engine` (master f81b20b7's `src_python/tin_engine` from
`git archive`, with a Release `_core` built from the same commit). The editable
install's import finder is dropped so the copy loads; both `__file__`s are
printed to stderr on every run.

Usage:
    PKG=<dir> python launch.py [--trace OUT.jsonl] [--patch 1,2,3] [--pieces K] -- <mesh args>

--trace   records every GEOS `buffer` call made through shapely's
          `BaseGeometry.buffer` (caller file:line, vertices, distance, keyword
          arguments, seconds), and each call of `feature_input.source_region`,
          `feature_input.read_source`, `feature_input.query_features` (time to
          exhaust its rows) and `feature_input.read_json`. One JSON line each.
          Instrumentation is a Python wrapper per call, no profiler.
--patch   the candidate fixes of docs/increments/perf-audit.md (branch
          worktree-perf-audit, 1da4a144), as its run_patched.py made them:
          1  every buffer of the same polygon with the same arguments computed
             once (memo keyed on the WKB and the arguments);
          2  feature_input.source_region grown from the domain's convex hull
             instead of the convex hull of the grown domain;
          3  a mitred buffer of a polygon with no holes and 5,000 vertices or
             more grown in pieces of K edges (--pieces, default 1000), each
             overlapping the next by one edge, flat ends, united with the
             polygon. With 1 and 3 together the memo wraps the pieces.
"""

import json
import os
import sys
import time
from pathlib import Path

sys.meta_path[:] = [f for f in sys.meta_path if "skbc" not in type(f).__module__]
sys.path.insert(0, str(Path(os.environ["PKG"]).resolve()))

import shapely  # noqa: E402
from shapely.geometry import Polygon  # noqa: E402
from shapely.geometry.base import BaseGeometry  # noqa: E402

import tin_engine  # noqa: E402
import tin_engine._core as core  # noqa: E402
from tin_engine import feature_input  # noqa: E402

print("tin_engine from", tin_engine.__file__, file=sys.stderr)
print("_core from", core.__file__, file=sys.stderr)

args = sys.argv[1:]
trace_path, patches, pieces = None, "", 1000
while args and args[0] != "--":
    flag, value, args = args[0], args[1], args[2:]
    if flag == "--trace":
        trace_path = value
    elif flag == "--patch":
        patches = value
    elif flag == "--pieces":
        pieces = int(value)
    else:
        raise SystemExit(f"unknown flag {flag}")
args = args[1:]
print("patches", patches or "none", "pieces", pieces, file=sys.stderr)


if "3" in patches:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from launch_pieces import piecewise_buffer
    plain = BaseGeometry.buffer

    def piecewise(self, distance, *a, **k):  # type: ignore[no-untyped-def]
        mitre = k.get("join_style") == "mitre" and not a
        if not (mitre and isinstance(self, Polygon) and not self.interiors):
            return plain(self, distance, *a, **k)
        if len(self.exterior.coords) - 1 < 5000:
            return plain(self, distance, *a, **k)
        return piecewise_buffer(self, distance, pieces)

    BaseGeometry.buffer = piecewise  # type: ignore[method-assign]

if "1" in patches:
    inner = BaseGeometry.buffer
    memo: dict = {}

    def shared(self, distance, *a, **k):  # type: ignore[no-untyped-def]
        key = (shapely.to_wkb(self), distance, a, tuple(sorted(k.items())))
        if key not in memo:
            memo[key] = inner(self, distance, *a, **k)
        return memo[key]

    BaseGeometry.buffer = shared  # type: ignore[method-assign]

if "2" in patches:

    def hull_region(domain, dem_crs, source_crs):  # type: ignore[no-untyped-def]
        hull = shapely.convex_hull(domain.polygon)
        ring = shapely.segmentize(hull.buffer(feature_input.MARGIN).exterior, feature_input.DENSIFY)
        xy = shapely.get_coordinates(ring)
        if not feature_input.same_crs(source_crs, dem_crs):
            xy = feature_input.reprojector(dem_crs, source_crs)(xy)
        out = shapely.convex_hull(shapely.multipoints(xy))
        assert isinstance(out, Polygon)
        return out

    feature_input.source_region = hull_region

if trace_path:
    out = open(trace_path, "w")  # noqa: SIM115

    def emit(**row):  # type: ignore[no-untyped-def]
        out.write(json.dumps(row) + "\n")
        out.flush()

    traced_buffer = BaseGeometry.buffer

    def buffer(self, distance, *a, **k):  # type: ignore[no-untyped-def]
        f = sys._getframe(1)
        where = f"{Path(f.f_code.co_filename).name}:{f.f_lineno}"
        t0 = time.perf_counter()
        r = traced_buffer(self, distance, *a, **k)
        emit(call="buffer", where=where, vertices=int(shapely.get_num_coordinates(self)),
             distance=float(distance), args=list(a), kwargs=k, s=time.perf_counter() - t0)  # fmt: skip
        return r

    BaseGeometry.buffer = buffer  # type: ignore[method-assign]

    def wrap(name, consume=False):  # type: ignore[no-untyped-def]
        real = getattr(feature_input, name)

        def timed(*a, **k):  # type: ignore[no-untyped-def]
            f = sys._getframe(1)
            where = f"{Path(f.f_code.co_filename).name}:{f.f_lineno}"
            t0 = time.perf_counter()
            r = real(*a, **k)
            if consume:
                r = list(r)
            emit(call=name, where=where, s=time.perf_counter() - t0,
                 rows=len(r) if consume else None)  # fmt: skip
            return iter(r) if consume else r

        setattr(feature_input, name, timed)

    wrap("source_region")
    wrap("read_source")
    wrap("query_features", consume=True)
    wrap("read_json")

from tin_engine.cli import app  # noqa: E402

sys.argv = ["rasputin", *args]
t0 = time.perf_counter()
try:
    app()
finally:
    print(f"launch: app wall {time.perf_counter() - t0:.3f} s", file=sys.stderr)
