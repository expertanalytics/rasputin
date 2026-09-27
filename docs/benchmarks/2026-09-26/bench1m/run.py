"""run.py WT args...: run `rasputin` from worktree WT, bypassing the editable finder.

Wraps cli.refine and cli.refine_start_stride to print BENCH lines on stderr:
refine wall seconds (the _core.refine call alone) and the start stride; and the
whole app() body in seconds (imports and interpreter start-up excluded).
"""
import sys, time
wt = sys.argv[1]
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, wt + "/src_python")
import tin_engine
assert tin_engine.__file__.startswith(wt), tin_engine.__file__
import tin_engine.cli as cli
_r = cli.refine
def refine(*a, **k):
    t = time.perf_counter(); out = _r(*a, **k)
    print(f"BENCH refine_s {time.perf_counter() - t:.4f}", file=sys.stderr)
    return out
cli.refine = refine
_s = cli.refine_start_stride
def rss(*a, **k):
    s = _s(*a, **k); print(f"BENCH start_stride {s}", file=sys.stderr); return s
cli.refine_start_stride = rss
sys.argv = ["rasputin"] + sys.argv[2:]
t = time.perf_counter()
try:
    cli.app()
finally:
    print(f"BENCH app_s {time.perf_counter() - t:.4f}", file=sys.stderr)
