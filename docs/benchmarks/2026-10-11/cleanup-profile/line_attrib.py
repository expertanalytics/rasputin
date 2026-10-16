"""Run `rasputin <argv>`; charge each wrapped shapely call to the first frame
in tin_engine/*.py above it (file:line). Output: argv[1], a text table."""
import sys, time, functools, collections
out, sys.argv = sys.argv[1], ["rasputin", *sys.argv[2:]]
import shapely, shapely.strtree
import tin_engine.feature_input as fi
import tin_engine.border_simplify as bs
tot = collections.defaultdict(lambda: [0, 0.0])
def caller():
    f = sys._getframe(2)
    while f is not None:
        fn = f.f_code.co_filename
        if "/tin_engine/" in fn:
            return f"{fn.rsplit('/',1)[1]}:{f.f_lineno}"
        f = f.f_back
    return "?"
def wrap(mod, name):
    orig = getattr(mod, name)
    @functools.wraps(orig)
    def w(*a, **k):
        t = time.perf_counter()
        try: return orig(*a, **k)
        finally:
            d = tot[(name, caller())]; d[0] += 1; d[1] += time.perf_counter() - t
    setattr(mod, name, w)
for n in ("intersection", "union_all", "coverage_union_all", "coverage_clean", "make_valid",
          "symmetric_difference", "transform", "get_coordinates", "linestrings", "distance",
          "contains_properly", "is_simple", "get_rings", "get_parts", "intersects"):
    wrap(fi.shapely, n)
for n in ("covers", "is_valid", "intersects"):   # methods on the domain polygon
    orig = getattr(shapely.geometry.base.BaseGeometry, n)
    def mk(orig, n):
        def w(self, *a, **k):
            t = time.perf_counter()
            try: return orig(self, *a, **k)
            finally:
                d = tot[(n + "()", caller())]; d[0] += 1; d[1] += time.perf_counter() - t
        return w
    setattr(shapely.geometry.base.BaseGeometry, n, property(mk(orig.fget, n)) if isinstance(orig, property) else mk(orig, n))
oq = shapely.strtree.STRtree.query
def q(self, *a, **k):
    t = time.perf_counter()
    try: return oq(self, *a, **k)
    finally:
        d = tot[("STRtree.query", caller())]; d[0] += 1; d[1] += time.perf_counter() - t
shapely.strtree.STRtree.query = q
oi = shapely.strtree.STRtree.__init__
def i(self, *a, **k):
    t = time.perf_counter()
    try: return oi(self, *a, **k)
    finally:
        d = tot[("STRtree()", caller())]; d[0] += 1; d[1] += time.perf_counter() - t
shapely.strtree.STRtree.__init__ = i
# whole clean-up functions, by wall time
for name in ("_clean", "snap_to_outline", "_outline_lines", "_chains", "_loops", "pre_clip"):
    orig = getattr(fi, name)
    def mk2(orig, name):
        @functools.wraps(orig)
        def w(*a, **k):
            t = time.perf_counter()
            try: return orig(*a, **k)
            finally:
                d = tot[(name, "whole")]; d[0] += 1; d[1] += time.perf_counter() - t
        return w
    setattr(fi, name, mk2(orig, name))
o_ring = fi._Outline.ring
def ring(self, *a, **k):
    t = time.perf_counter()
    try: return o_ring(self, *a, **k)
    finally:
        d = tot[("_Outline.ring", "whole")]; d[0] += 1; d[1] += time.perf_counter() - t
fi._Outline.ring = ring
o_add = fi._Tally._add
def add(self, *a, **k):
    t = time.perf_counter()
    try: return o_add(self, *a, **k)
    finally:
        d = tot[("_Tally._add", "whole")]; d[0] += 1; d[1] += time.perf_counter() - t
fi._Tally._add = add
o_sb = bs.simplify_borders
def sb(*a, **k):
    t = time.perf_counter()
    try: return o_sb(*a, **k)
    finally:
        d = tot[("simplify_borders", "whole")]; d[0] += 1; d[1] += time.perf_counter() - t
fi.simplify_borders = sb
from tin_engine.cli import app
t0 = time.perf_counter()
try: app()
except SystemExit as exc: print("exit", exc.code)
wall = time.perf_counter() - t0
with open(out, "w") as f:
    f.write(f"# wall of app(): {wall:.3f} s\n# seconds  calls  function  charged-to (tin_engine file:line, or 'whole' = the function's own wall time)\n")
    for (n, c), (k, s) in sorted(tot.items(), key=lambda kv: -kv[1][1]):
        if s >= 0.005: f.write(f"{s:8.3f} {k:7d}  {n:20s} {c}\n")
