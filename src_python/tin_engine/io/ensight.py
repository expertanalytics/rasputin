"""EnSight Gold binary encoders for an animation of values per triangle: bytes and text out.

The geometry is written once; each value at each frame is one small file. The
coordinates are float32, so x and y are stored as offsets from an origin the
caller gives (written in the geometry's second description line), which keeps
millimetres at UTM northings. The case file lists the time values one per
line: VTK's reader, which ParaView uses, reads at most 256 characters a line.
Pure: no path, no file.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt


def _text(line: str) -> bytes:
    return line.encode("ascii").ljust(80, b"\0")


def _int(value: int) -> bytes:
    return np.int32(value).tobytes()


def geometry_bytes(
    vertices: npt.ArrayLike, triangles: npt.ArrayLike, *, origin: tuple[float, float]
) -> bytes:
    """One part of `tria3` triangles (0-based corners in, 1-based written), the
    coordinates as float32 offsets from `origin`."""
    v = np.asarray(vertices, dtype=np.float64).reshape(-1, 3)
    t = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    header = ("C Binary", "rasputin snow", f"origin {origin[0]!r} {origin[1]!r}",
              "node id off", "element id off", "part")  # fmt: skip
    columns = (v[:, 0] - origin[0], v[:, 1] - origin[1], v[:, 2])
    return b"".join(
        [*map(_text, header), _int(1), _text("mesh"), _text("coordinates"), _int(len(v)),
         *(c.astype("<f4").tobytes() for c in columns),
         _text("tria3"), _int(len(t)), (t + 1).astype("<i4").tobytes()]
    )  # fmt: skip


def variable_bytes(name: str, values: npt.ArrayLike) -> bytes:
    """One float32 value per triangle of part 1, described as `name`."""
    data = np.asarray(values).astype("<f4").tobytes()
    return b"".join([_text(name), _text("part"), _int(1), _text("tria3"), data])


def frame_file(name: str, index: int) -> str:
    """The file of `name` at frame `index`, by the case file's pattern."""
    return f"{name}{index:04d}.var"


def case_text(
    *, geometry: str, timed: Sequence[str], static: Sequence[str], times_hours: Sequence[float]
) -> str:
    """The case file: `timed` variables on time set 1 (`<name>****.var`, frames
    from 0), `static` ones in `<name>.var`, the times one per line."""
    variables = [f"scalar per element: 1 {v} {v}****.var" for v in timed]
    variables += [f"scalar per element: {v} {v}.var" for v in static]
    lines = ["FORMAT", "type: ensight gold", "GEOMETRY", f"model: {geometry}", "VARIABLE",
             *variables, "TIME", "time set: 1", f"number of steps: {len(times_hours)}",
             "filename start number: 0", "filename increment: 1", "time values:",
             *map(repr, map(float, times_hours))]  # fmt: skip
    return "\n".join(lines) + "\n"
