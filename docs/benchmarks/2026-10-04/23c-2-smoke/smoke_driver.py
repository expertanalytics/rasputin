"""Runs `rasputin mesh` from a package directory and records every `cli.refine`
call (one per piece in a cut run): start-quality inserts and skips, rounds,
inserted, max error, seconds. One JSON line per call to RECORD (@perf).
python smoke_driver.py PKG RECORD.jsonl -- <rasputin mesh argv>"""
import json
import sys
import time

split = sys.argv.index("--")
pkg, record = sys.argv[1], sys.argv[2]
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, pkg)
import tin_engine.cli as cli  # noqa: E402

assert cli.__file__.startswith(pkg), cli.__file__
real = cli.refine
fh = open(record, "w")


def refine(*a, **k):
    t0 = time.perf_counter()
    out = real(*a, **k)
    fh.write(json.dumps({"seconds": round(time.perf_counter() - t0, 4), "quality_inserted": out.quality_inserted,
                         "quality_skipped": out.quality_skipped, "rounds": out.rounds, "inserted": out.inserted,
                         "max_error": out.max_error, "status": str(getattr(out, "status", ""))}) + "\n")  # fmt: skip
    fh.flush()
    return out


cli.refine = refine
sys.argv = ["rasputin", *sys.argv[split + 1 :]]
cli.app()
