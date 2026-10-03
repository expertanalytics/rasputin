"""``rasputin`` with phase markers: the CLI, in this process, unchanged except
that a few calls are wrapped to append a line to a marker file when they begin
and end. A measurement script, not production code; nothing imports it.

    python phase_driver.py MARKERS.jsonl -- mesh ARGS...

A marker is ``{"t": <time.time()>, "ev": "B"|"E", "name": ..., "fp": ...,
"life_max": ...}``, the last two this process's footprint and lifetime peak
footprint at the marker, written and
flushed before (B) and after (E) the call, so `run_phases.py`'s memory samples
(same wall clock) can be cut by phase. Wrapped:

- every ``PhaseClock.phase`` block (the ``--stats`` rows timed in Python:
  ``decode``, ``start mesh: *``, ``check points: store`` (the freeze),
  ``trim``, ``write: *``);
- ``dem_input.assemble`` (the source mosaic), ``dem_input.resample`` (the
  target grid), and inside it the call that builds the returned ``DemTile``
  (marker ``resample: DemTile copy``, the name kept from 2026-10-02 so
  ``analyse.py`` reads both trees): before 15e the public constructor (a
  copy), from 15e on ``DemTile._adopt`` (no copy, fix 2);
- ``cli.refine`` (phase 1), ``final_check.run`` (store fill + phase 2) and
  ``final_check.refine_points`` (phase 2's refinement).

Wrappers read no result (an outcome's arrays are copies when read), except
``CheckPoints.size`` at phase 2's start, which is a count.
"""

from __future__ import annotations

import functools
import json
import sys
import ctypes
import os
import time
from contextlib import contextmanager
from typing import Any

from tin_engine import cli, dem_input, final_check, stats, target_grid

_out = open(sys.argv[1], "w")  # noqa: SIM115 - kept open for the run
_libc, _buf = ctypes.CDLL("/usr/lib/libSystem.B.dylib"), ctypes.create_string_buffer(512)


def _usage() -> dict[str, int]:
    """This process's ``phys_footprint`` now and its lifetime maximum
    (``struct rusage_info_v4``, offsets 72 and 240)."""
    _libc.proc_pid_rusage(os.getpid(), 4, _buf)
    q = lambda off: int.from_bytes(_buf.raw[off : off + 8], "little")  # noqa: E731
    return {"fp": q(72), "life_max": q(240)}


def mark(ev: str, name: str, **extra: Any) -> None:
    rec = {"t": time.time(), "ev": ev, "name": name, **_usage(), **extra}
    _out.write(json.dumps(rec) + "\n")
    _out.flush()


def wrap(module: Any, attr: str, name: str, extra: Any = None, after: Any = None) -> None:
    f = getattr(module, attr)

    @functools.wraps(f)
    def g(*a: Any, **k: Any) -> Any:
        mark("B", name, **(extra(a) if extra else {}))
        result: Any = None
        try:
            result = f(*a, **k)
            return result
        finally:
            mark("E", name, **(after(result) if after and result is not None else {}))

    setattr(module, attr, g)


_phase = stats.PhaseClock.phase


@contextmanager
def _marked(self: stats.PhaseClock, name: str) -> Any:
    mark("B", f"clock: {name}")
    try:
        with _phase(self, name):
            yield
    finally:
        mark("E", f"clock: {name}")


stats.PhaseClock.phase = _marked  # type: ignore[method-assign,assignment]
wrap(dem_input, "assemble", "assemble", after=lambda r: {"shape": list(r.tile.array.shape), "nbytes": r.tile.array.nbytes})
wrap(dem_input, "resample", "resample",
     lambda a: {"rows": a[0].rows, "cols": a[0].cols, "block_rows": _block_rows(a[0].cols)})
if hasattr(target_grid, "rows_per_block"):  # 15e: resample ends in DemTile._adopt

    class _AdoptOnly:  # target_grid uses the name DemTile only for this call at run time
        _adopt = staticmethod(target_grid.DemTile._adopt)

    target_grid.DemTile = _AdoptOnly  # type: ignore[misc,assignment]
    wrap(_AdoptOnly, "_adopt", "resample: DemTile copy")
    _block_rows = target_grid.rows_per_block
else:  # before 15e: the public constructor copies the canvas, 256-row blocks
    wrap(target_grid, "DemTile", "resample: DemTile copy")
    _block_rows = lambda cols: 256  # noqa: E731
wrap(cli, "refine", "refine phase 1")
wrap(final_check, "run", "final_check.run")
wrap(final_check, "refine_points", "refine_points", lambda a: {"store_size": a[0].size})

if __name__ == "__main__":
    sys.argv = ["rasputin", *sys.argv[sys.argv.index("--") + 1 :]]
    import tin_engine, tin_engine._core as _c
    print("tin_engine from", tin_engine.__file__, "_core from", _c.__file__, file=sys.stderr)
    mark("B", "process")
    try:
        cli.app()
    finally:
        mark("E", "process")
