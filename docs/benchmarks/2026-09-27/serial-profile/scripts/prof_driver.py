"""prof_driver.py --pkg DIR --threads N --repeat K [--pause S] -- <rasputin argv>

bench.py's child technique (pkg first on sys.path, wrap cli.refine), but the
wrapped refine is called K times with the same arguments, each timed, and the
RefineOutcome phase seconds printed per call as one JSON line (PHASES ...).
--pause S sleeps S seconds before the first call and prints the pid, so
`sample` can attach to a steady window of refine calls.
"""
import json
import os
import sys
import time

argv = sys.argv[1:]
split = argv.index("--")
own, rasputin = argv[:split], argv[split + 1 :]
pkg = own[own.index("--pkg") + 1]
threads = int(own[own.index("--threads") + 1])
repeat = int(own[own.index("--repeat") + 1])
pause = float(own[own.index("--pause") + 1]) if "--pause" in own else 0.0
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, pkg)
import tin_engine._core as core  # noqa: E402
import tin_engine.cli as cli  # noqa: E402

assert core.__file__.startswith(pkg), core.__file__
real = vars(cli)["refine"]


def refine(*args, **kwargs):
    kwargs["threads"] = threads
    if pause:
        print(f"PID {os.getpid()}", file=sys.stderr, flush=True)
        time.sleep(pause)
    out = None
    for i in range(repeat):
        t0 = time.perf_counter()
        out = real(*args, **kwargs)
        wall = time.perf_counter() - t0
        rec = dict(call=i, threads=threads, refine_s=wall, legalise_s=out.legalise_seconds,
                   quality_s=out.quality_seconds, scan_s=out.scan_seconds,
                   split_s=out.split_seconds, rounds=out.rounds, inserted=out.inserted,
                   flips=out.flips, max_error=out.max_error,
                   triangles=len(out.triangles), vertices=len(out.vertices))
        rec["rest_s"] = wall - rec["legalise_s"] - rec["quality_s"] - rec["scan_s"] - rec["split_s"]
        print("PHASES " + json.dumps(rec), file=sys.stderr, flush=True)
    return out


vars(cli)["refine"] = refine
cli.app(args=rasputin, prog_name="rasputin", standalone_mode=False)
