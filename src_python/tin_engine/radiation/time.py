"""Times for the sun's position, as UTC milliseconds since the Unix epoch.

Only timezone-aware times are taken: a time without a zone is refused rather
than guessed, since a guessed zone shifts every sun position by whole hours.
"""

from __future__ import annotations

import datetime as dt
from collections.abc import Iterable

import numpy as np
import numpy.typing as npt

_EPOCH = dt.datetime(1970, 1, 1, tzinfo=dt.UTC)
_MS = dt.timedelta(milliseconds=1)


def to_utc_ms(
    times: npt.NDArray[np.datetime64] | Iterable[dt.datetime | str],
) -> npt.NDArray[np.int64]:
    """Turn times into a 1-D int64 array of UTC milliseconds since the epoch.

    Takes a ``datetime64`` array, read as UTC, or an iterable of timezone-aware
    ``datetime``s and ISO 8601 strings with an offset (or ``Z``); finer than a
    millisecond is floored. A naive ``datetime`` or a string without an offset
    is a ``ValueError`` saying the time zone is missing.
    """
    if isinstance(times, np.ndarray) and times.dtype.kind == "M":
        return times.astype("datetime64[ms]").astype(np.int64).ravel()
    out: list[int] = []
    for t in times:
        d = dt.datetime.fromisoformat(t) if isinstance(t, str) else t
        if d.utcoffset() is None:
            raise ValueError(
                f"{t!r} has no time zone; give it one (an offset such as +00:00, or Z)"
            )
        out.append((d - _EPOCH) // _MS)
    return np.array(out, dtype=np.int64)
