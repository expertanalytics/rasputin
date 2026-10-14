"""An independent EnSight Gold reader. Test support only; nothing in `src_python` uses it.

It reads what `io/ensight.py` and `rasputin snow` write, the way VTK's
EnSight Gold binary reader (what ParaView opens a `.case` with) reads it, and
takes every count from the bytes it was handed, never from the writer. It knows
the subset the snow animation uses and refuses anything outside it:

- the case file: `FORMAT` (`type: ensight gold`), `GEOMETRY` (`model:` one
  file, no time set), `VARIABLE` (`scalar per element:` lines, with the time
  set `1` and a `*` pattern, or without a time set and a plain file name) and
  `TIME` (time set 1, its step count, start number, increment and the values,
  which VTK reads one line of at most 256 characters at a time);
- the geometry file: `C Binary`, two description lines, `node id off`,
  `element id off`, one `part` with `coordinates` (x, y, z as little-endian
  float32 blocks) and one `tria3` block of 1-based int32 corners;
- a per-element variable file: a description line, `part`, the part number,
  `tria3` and one float32 per triangle.

Every 80-byte text record is NUL- or space-padded ASCII. A count that does
not match the bytes, a corner outside the coordinates, or trailing bytes is
an `EnsightError`: a file VTK would half-read is a file this one raises on.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt

RECORD = 80
#: VTK's EnSight reader reads a case-file line into a 256-byte buffer.
VTK_CASE_LINE = 256


class EnsightError(ValueError):
    """A file outside the subset, or inconsistent with itself."""


@dataclass(frozen=True)
class Geometry:
    descriptions: tuple[str, str]
    part: int
    part_description: str
    points: npt.NDArray[np.float32]  # (n, 3)
    triangles: npt.NDArray[np.int64]  # (T, 3), 0-based


@dataclass(frozen=True)
class Variable:
    description: str
    part: int
    values: npt.NDArray[np.float32]


@dataclass(frozen=True)
class Case:
    geometry: str
    #: name -> file pattern with `*`s, for variables on time set 1.
    timed: dict[str, str]
    #: name -> file, for variables without a time set.
    static: dict[str, str]
    times: tuple[float, ...]
    start: int
    increment: int
    lines: tuple[str, ...]

    def frame_file(self, name: str, index: int) -> str:
        """The file of `name` at time step `index`, by the case's own pattern."""
        pattern = self.timed[name]
        stars = pattern.count("*")
        number = self.start + index * self.increment
        return pattern.replace("*" * stars, f"{number:0{stars}d}")


class _Cursor:
    def __init__(self, data: bytes, what: str) -> None:
        self.data, self.at, self.what = data, 0, what

    def text(self) -> str:
        chunk = self.take(RECORD)
        try:
            return chunk.rstrip(b"\0 ").decode("ascii")
        except UnicodeDecodeError as e:
            raise EnsightError(f"{self.what}: a text record is not ASCII") from e

    def int32(self) -> int:
        return int(np.frombuffer(self.take(4), "<i4")[0])

    def floats(self, count: int) -> npt.NDArray[np.float32]:
        return np.frombuffer(self.take(4 * count), "<f4").copy()

    def ints(self, count: int) -> npt.NDArray[np.int32]:
        return np.frombuffer(self.take(4 * count), "<i4").copy()

    def take(self, size: int) -> bytes:
        if size < 0 or self.at + size > len(self.data):
            raise EnsightError(
                f"{self.what}: ends at byte {len(self.data)}, needs {self.at + size}"
            )
        chunk = self.data[self.at : self.at + size]
        self.at += size
        return chunk

    def expect(self, word: str) -> None:
        got = self.text()
        if got != word:
            raise EnsightError(f"{self.what}: expected {word!r}, found {got!r}")

    def end(self) -> None:
        if self.at != len(self.data):
            raise EnsightError(f"{self.what}: {len(self.data) - self.at} trailing bytes")


def read_geometry(data: bytes) -> Geometry:
    c = _Cursor(data, "geometry")
    c.expect("C Binary")
    descriptions = (c.text(), c.text())
    c.expect("node id off")
    c.expect("element id off")
    c.expect("part")
    part = c.int32()
    part_description = c.text()
    c.expect("coordinates")
    n = c.int32()
    if n < 0:
        raise EnsightError(f"geometry: {n} nodes")
    xyz = np.column_stack([c.floats(n), c.floats(n), c.floats(n)])
    c.expect("tria3")
    count = c.int32()
    if count < 0:
        raise EnsightError(f"geometry: {count} triangles")
    corners = c.ints(3 * count).reshape(-1, 3).astype(np.int64)
    c.end()
    if corners.size and (corners.min() < 1 or corners.max() > n):
        raise EnsightError(f"geometry: a corner outside 1..{n}")
    return Geometry(descriptions, part, part_description, xyz, corners - 1)


def read_variable(data: bytes, triangles: int) -> Variable:
    c = _Cursor(data, "variable")
    description = c.text()
    c.expect("part")
    part = c.int32()
    c.expect("tria3")
    values = c.floats(triangles)
    c.end()
    return Variable(description, part, values)


_SCALAR = re.compile(r"scalar per element:\s+(?:(\d+)\s+)?(\S+)\s+(\S+)")


def read_case(text: str) -> Case:
    lines = tuple(text.splitlines())
    for line in lines:
        if len(line) >= VTK_CASE_LINE:
            raise EnsightError(f"case: a line of {len(line)} characters, VTK reads {VTK_CASE_LINE}")
    sections: dict[str, list[str]] = {}
    current = ""
    for line in lines:
        if not line.strip():
            continue
        if line.strip() in ("FORMAT", "GEOMETRY", "VARIABLE", "TIME"):
            current = line.strip()
            if current in sections:
                raise EnsightError(f"case: {current} twice")
            sections[current] = []
        elif not current:
            raise EnsightError(f"case: {line!r} before any section")
        else:
            sections[current].append(line.strip())
    if sections.get("FORMAT") != ["type: ensight gold"]:
        raise EnsightError(f"case: FORMAT is {sections.get('FORMAT')}")
    model = sections.get("GEOMETRY", [])
    if len(model) != 1 or not model[0].startswith("model:") or len(model[0].split()) != 2:
        raise EnsightError(f"case: GEOMETRY is {model}, not one model file without a time set")
    timed: dict[str, str] = {}
    static: dict[str, str] = {}
    for line in sections.get("VARIABLE", []):
        match = _SCALAR.fullmatch(line)
        if not match:
            raise EnsightError(f"case: variable line {line!r}")
        tset, name, file = match.groups()
        if tset is None:
            if "*" in file:
                raise EnsightError(f"case: {name} has a pattern and no time set")
            static[name] = file
        elif tset == "1" and re.fullmatch(r"[^*]*\*+[^*]*", file):
            timed[name] = file
        else:
            raise EnsightError(f"case: {line!r} is not on time set 1 with one run of '*'")
    time = sections.get("TIME", [])
    keys: dict[str, str] = {}
    values: list[float] = []
    in_values = False
    for line in time:
        if in_values:
            values += [float(word) for word in line.split()]
            continue
        key, _, rest = line.partition(":")
        keys[key.strip()] = rest.strip()
        if key.strip() == "time values":
            in_values = True
            values += [float(word) for word in rest.split()]
    if keys.get("time set") != "1":
        raise EnsightError(f"case: time set is {keys.get('time set')!r}")
    steps = int(keys["number of steps"])
    if steps != len(values):
        raise EnsightError(f"case: {steps} steps, {len(values)} time values")
    return Case(
        geometry=model[0].split()[1],
        timed=timed,
        static=static,
        times=tuple(values),
        start=int(keys["filename start number"]),
        increment=int(keys["filename increment"]),
        lines=lines,
    )


@dataclass(frozen=True)
class Animation:
    case: Case
    geometry: Geometry

    @classmethod
    def read(cls, case_path: Path) -> Animation:
        case = read_case(case_path.read_text(encoding="ascii"))
        geometry = read_geometry((case_path.parent / case.geometry).read_bytes())
        return cls(case, geometry)

    def frame(self, folder: Path, name: str, index: int) -> npt.NDArray[np.float32]:
        data = (folder / self.case.frame_file(name, index)).read_bytes()
        return read_variable(data, len(self.geometry.triangles)).values

    def static(self, folder: Path, name: str) -> npt.NDArray[np.float32]:
        data = (folder / self.case.static[name]).read_bytes()
        return read_variable(data, len(self.geometry.triangles)).values
