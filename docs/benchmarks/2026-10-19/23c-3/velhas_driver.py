"""Runs `rasputin mesh` from a package directory and records, one JSON line per
`cli.refine` call (one per piece in a cut run) to RECORD: start-quality inserts,
rounds, inserted, max error, seconds, and the peak physical footprint of the
interval before the call ("pre": setup, the previous piece's write and, on
23c-3, this piece's load) and inside it ("in"). A last line holds the run's
tail ("post"), the process's lifetime peak and every PhaseClock's phases
(--stats is refused on a cut run). Peaks come from macOS
`proc_pid_rusage` (RUSAGE_INFO_V4): `ri_interval_max_phys_footprint` resets on
every read, `ri_lifetime_max_phys_footprint` does not (@perf).
python velhas_driver.py PKG RECORD.jsonl -- <rasputin mesh argv>"""
import ctypes
import json
import os
import sys
import time

split = sys.argv.index("--")
pkg, record = sys.argv[1], sys.argv[2]
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, pkg)

_libc = ctypes.CDLL("libc.dylib")
_buf = (ctypes.c_uint64 * 64)()  # 16-byte uuid, then u64 fields; v4 is 35 fields


def peaks():
    """(interval max, lifetime max, current) phys footprint in bytes."""
    if _libc.proc_pid_rusage(os.getpid(), 4, ctypes.byref(_buf)) != 0:
        raise OSError("proc_pid_rusage")
    f = _buf[2:]  # skip the uuid's two words
    return f[33], f[28], f[7]


peaks()  # start the first interval here
import tin_engine.cli as cli  # noqa: E402

assert cli.__file__.startswith(pkg), cli.__file__
real = cli.refine
fh = open(record, "w")


def refine(*a, **k):
    pre = peaks()[0]
    t0 = time.perf_counter()
    out = real(*a, **k)
    sec = time.perf_counter() - t0
    inside, _, now = peaks()
    fh.write(json.dumps({"seconds": round(sec, 4), "quality_inserted": out.quality_inserted,
                         "quality_skipped": out.quality_skipped, "rounds": out.rounds, "inserted": out.inserted,
                         "max_error": out.max_error, "status": str(getattr(out, "status", "")),
                         "pre_peak": pre, "in_peak": inside, "after": now}) + "\n")  # fmt: skip
    fh.flush()
    return out


clocks = []


class KeptClock(cli.PhaseClock):  # bench.py's technique; a cut run may make several
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        clocks.append(self)


cli.refine = refine
cli.PhaseClock = KeptClock
sys.argv = ["rasputin", *sys.argv[split + 1 :]]
try:
    cli.app()
finally:
    post, life, now = peaks()
    phases = [dict(c.phases()) for c in clocks]
    fh.write(json.dumps({"post_peak": post, "lifetime_peak": life, "end": now, "phases": phases}) + "\n")
    fh.close()
