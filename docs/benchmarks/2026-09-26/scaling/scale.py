"""scale.py N: run the CLI at 1 m on the quarter circle with refine forced to N threads.
Prints BENCH refine_s for the _core.refine call alone."""
import sys, time
n = int(sys.argv[1])
import tin_engine.cli as cli
_r = cli.refine
def refine(*a, **k):
    k["threads"] = n
    t = time.perf_counter(); out = _r(*a, **k)
    print(f"BENCH refine_s {time.perf_counter() - t:.4f}", file=sys.stderr)
    return out
cli.refine = refine
sys.argv = ["rasputin"] + sys.argv[2:]
cli.app()
