"""Whether a call released the GIL: count a ticker thread's ticks while it runs.

`docs/increments/python-audit.md`, section 10 (R10). Test support only. Each
suite keeps its own thresholds and its own control test, because what a
held and a released call let through differs between the calls they measure.

Two measures. `ticks_during` counts ticks, which other processes' load
distorts: a starved thread hands the GIL over at the edges of a held call, and
under `pytest -n auto` a held call let up to 26 ticks through. `silence_during`
measures the longest stretch of the call with no tick, as a fraction of the
call. Load can only lengthen a silence, never fill one, because the ticker
cannot tick without the GIL: a held call stays near 1 and a released one near
0 however busy the machine is.
"""

from __future__ import annotations

import threading
import time
from collections.abc import Callable
from itertools import pairwise


class _Ticker:
    """A thread that ticks once per millisecond, and cannot tick without the GIL.

    `time.sleep` releases the GIL and must re-acquire it to continue, so a tick
    is direct evidence that the GIL was available. The tick rate is bounded by
    the sleep, so the ticker never starves the thread it is observing.
    """

    def __init__(self) -> None:
        self.ticks = 0
        self.stamps: list[float] = []
        self._stop = threading.Event()
        self._ready = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self) -> None:
        self._ready.set()
        while not self._stop.is_set():
            self.ticks += 1
            self.stamps.append(time.perf_counter())
            time.sleep(0.001)

    def __enter__(self) -> _Ticker:
        self._thread.start()
        assert self._ready.wait(timeout=10.0), "ticker never started"
        return self

    def __exit__(self, *exc: object) -> None:
        self._stop.set()
        self._thread.join(timeout=10.0)


def ticks_during[T](call: Callable[[], T]) -> tuple[T, int, float]:
    """Return the call's result, the ticks observed during it, and its duration."""
    with _Ticker() as ticker:
        time.sleep(0.02)  # let the ticker reach steady state
        before = ticker.ticks
        started = time.perf_counter()
        result = call()
        elapsed = time.perf_counter() - started
        after = ticker.ticks
    return result, after - before, elapsed


def silence_during[T](call: Callable[[], T]) -> tuple[T, float, float]:
    """Return the call's result, the longest stretch of it with no tick as a
    fraction of its duration, and its duration."""
    with _Ticker() as ticker:
        time.sleep(0.02)  # let the ticker reach steady state
        started = time.perf_counter()
        result = call()
        ended = time.perf_counter()
    marks = [started, *(t for t in ticker.stamps if started < t < ended), ended]
    longest = max(b - a for a, b in pairwise(marks))
    return result, longest / (ended - started), ended - started
