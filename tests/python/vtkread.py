"""An independent legacy-VTK reader. Test support only; nothing in `src_python` uses it.

Increment 13's pure suite goes through this parser for the same reason
increment 10's goes through `plyread`: it **reads the file it was handed**, and
takes every count, type and width from a declaration in the blob. It never asks
the writer where anything is (`PRINCIPLES.md` B1).

It is deliberately narrow and deliberately strict. It knows the subset of
legacy VTK 4.2 that `13-bundled-mesh.md` rules on -- `DATASET POLYDATA`, a
dataset `FIELD`, `POINTS`, `LINES`, `POLYGONS`, and a `CELL_DATA` holding
`SCALARS` and `FIELD` blocks -- and refuses any keyword outside it, any count
that does not match its declaration, and any trailing bytes. A file the real
reader would half-load is a file this one raises on.

What it takes from VTK rather than from the writer, and where that was checked:

- **ASCII strings** are one per line, with VTK's `%XX` escapes decoded. VTK's
  own writer emits `+proj=utm%20+zone=33`, and its reader decodes `%20` and
  `%25` (the design's finding 3). A raw space in an ASCII string is refused
  here, because ruling 7 says the writer never relies on it.
- **BINARY numbers** are big-endian, by the format's definition (ruling 3). A
  single newline may follow a binary block, as VTK's own writer emits one.
- **BINARY strings** are length-prefixed and not percent-encoded. The prefix's
  top two bits give the width of the length: `11` is six bits in one byte,
  `10` fourteen bits in two, `01` thirty bits in four, `00` sixty-two bits in
  eight. Measured by writing a `vtkStringArray` holding
  `+proj=utm +zone=33 a%b` with `vtkPolyDataWriter` 9.7.0 in binary mode: the
  value was written as `0xd6` (`0xc0 | 22`) followed by the 22 raw bytes. A
  64-character value -- the length of a sha256 fingerprint in hex -- was
  written with the two-byte prefix `0x8040`, and a 63-character one with
  `0xff`.

The parser records the order in which top-level sections appeared, because
ruling 4's cell order is a property of the file and not only of the arrays.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

#: Legacy type names onto big-endian numpy types. The ASCII path uses the same
#: table for the value type; only the binary path cares about the byte order.
_DTYPES = {
    "bit": "u1",
    "char": "i1",
    "unsigned_char": "u1",
    "short": ">i2",
    "unsigned_short": ">u2",
    "int": ">i4",
    "unsigned_int": ">u4",
    "long": ">i8",
    "unsigned_long": ">u8",
    "vtktypeint64": ">i8",
    "vtktypeuint64": ">u8",
    "float": ">f4",
    "double": ">f8",
}

_ESCAPE = re.compile(r"%([0-9A-Fa-f]{2})")


@dataclass(frozen=True)
class Array:
    """One data array as declared: its type, its component count, its values.

    `values` is a numpy array for numeric types and a tuple of decoded strings
    for `string`. `raw` holds the strings exactly as they stood in an ASCII
    file, before `%XX` decoding, so a test can assert on the encoding itself.
    """

    name: str
    type_name: str
    components: int
    values: np.ndarray | tuple[str, ...]
    raw: tuple[bytes, ...] = ()


@dataclass
class VtkFile:
    """Everything the file declared, in the shape a test asserts on."""

    version: str = ""
    title: str = ""
    encoding: str = ""
    dataset: str = ""
    sections: list[str] = field(default_factory=list)
    field_data_name: str | None = None
    field_data: dict[str, Array] = field(default_factory=dict)
    points_type: str = ""
    points: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    lines_header: tuple[int, int] | None = None
    lines: list[np.ndarray] = field(default_factory=list)
    polygons_header: tuple[int, int] | None = None
    polygons: list[np.ndarray] = field(default_factory=list)
    point_count: int | None = None
    point_scalars: dict[str, Array] = field(default_factory=dict)
    cell_count: int | None = None
    cell_sections: list[str] = field(default_factory=list)
    scalars: dict[str, Array] = field(default_factory=dict)
    lookup_tables: dict[str, str] = field(default_factory=dict)
    cell_fields: dict[str, dict[str, Array]] = field(default_factory=dict)

    @property
    def cells(self) -> int:
        return len(self.lines) + len(self.polygons)

    def cell_array(self, name: str) -> Array:
        """A cell array by name, from `SCALARS` or from any cell `FIELD` block."""
        if name in self.scalars:
            return self.scalars[name]
        for arrays in self.cell_fields.values():
            if name in arrays:
                return arrays[name]
        raise KeyError(f"no cell array {name!r}")


class _Cursor:
    """A byte position in the blob, and the three ways of consuming from it."""

    def __init__(self, blob: bytes) -> None:
        self.blob = blob
        self.at = 0
        self.binary = False

    def raw_line(self) -> bytes:
        stop = self.blob.find(b"\n", self.at)
        if stop < 0:
            raise ValueError(f"unterminated line at byte {self.at}")
        line = self.blob[self.at : stop]
        self.at = stop + 1
        return line

    def keyword_line(self) -> list[str] | None:
        """The next non-blank line, split, or None at the end of the file."""
        while self.at < len(self.blob):
            line = self.raw_line()
            if line.strip():
                return line.decode("ascii").split()
        return None

    def numbers(self, count: int, type_name: str) -> np.ndarray:
        dtype = np.dtype(_DTYPES[type_name])
        if self.binary:
            size = dtype.itemsize * count
            if self.at + size > len(self.blob):
                raise ValueError(f"{count} {type_name} values run past the end of the file")
            values = np.frombuffer(self.blob, dtype=dtype, count=count, offset=self.at)
            self.at += size
            self._optional_newline()
            return np.array(values)
        tokens: list[str] = []
        while len(tokens) < count:
            tokens += self.raw_line().decode("ascii").split()
        if len(tokens) != count:
            raise ValueError(f"expected {count} values, a line carried {len(tokens)}")
        convert = float if dtype.kind == "f" else int
        return np.array([convert(t) for t in tokens], dtype=dtype.newbyteorder("="))

    def strings(self, count: int) -> tuple[tuple[str, ...], tuple[bytes, ...]]:
        if self.binary:
            values = tuple(self._binary_string() for _ in range(count))
            self._optional_newline()
            return values, ()
        raw = tuple(self.raw_line() for _ in range(count))
        for line in raw:
            if b" " in line or b"\t" in line:
                raise ValueError(f"an ASCII string carries raw whitespace: {line!r}")
        decoded = tuple(
            _ESCAPE.sub(lambda m: chr(int(m.group(1), 16)), line.decode("ascii")) for line in raw
        )
        return decoded, raw

    def _binary_string(self) -> str:
        head = self.blob[self.at]
        width = {0b11: 1, 0b10: 2, 0b01: 4, 0b00: 8}[head >> 6]
        prefix = int.from_bytes(self.blob[self.at : self.at + width], "big")
        length = prefix & ((1 << (8 * width - 2)) - 1)
        self.at += width
        value = self.blob[self.at : self.at + length]
        if len(value) != length:
            raise ValueError("a binary string runs past the end of the file")
        self.at += length
        return value.decode("ascii")

    def _optional_newline(self) -> None:
        if self.blob[self.at : self.at + 1] == b"\n":
            self.at += 1

    def rest_is_blank(self) -> bool:
        return not self.blob[self.at :].strip()


def read_vtk(blob: bytes) -> VtkFile:
    """Parse a whole legacy VTK PolyData file, refusing anything outside the subset."""
    cur = _Cursor(blob)
    out = VtkFile()
    first = cur.raw_line().decode("ascii")
    prefix = "# vtk DataFile Version "
    if not first.startswith(prefix):
        raise ValueError(f"not a legacy VTK file: it begins {blob[:32]!r}")
    out.version = first[len(prefix) :]
    out.title = cur.raw_line().decode("ascii")
    out.encoding = cur.raw_line().decode("ascii").strip()
    if out.encoding not in ("ASCII", "BINARY"):
        raise ValueError(f"line 3 must be ASCII or BINARY; got {out.encoding!r}")
    cur.binary = out.encoding == "BINARY"
    words = cur.keyword_line()
    if not words or words[0] != "DATASET" or len(words) != 2:
        raise ValueError(f"line 4 must be DATASET <type>; got {words}")
    out.dataset = words[1]

    while (words := cur.keyword_line()) is not None:
        keyword = words[0]
        out.sections.append(keyword)
        if keyword == "FIELD" and out.cell_count is None:
            out.field_data_name = words[1]
            out.field_data = _field_arrays(cur, int(words[2]))
        elif keyword == "POINTS":
            count, out.points_type = int(words[1]), words[2]
            out.points = cur.numbers(3 * count, out.points_type).reshape(count, 3)
        elif keyword in ("LINES", "POLYGONS"):
            count, size = int(words[1]), int(words[2])
            cells = _cells(cur.numbers(size, "int"), count, keyword)
            if keyword == "LINES":
                out.lines_header, out.lines = (count, size), cells
            else:
                out.polygons_header, out.polygons = (count, size), cells
        elif keyword == "POINT_DATA":
            out.point_count = int(words[1])
            _point_scalars(cur, out)
        elif keyword == "CELL_DATA":
            out.cell_count = int(words[1])
            _cell_data(cur, out)
            break
        else:
            raise ValueError(f"unexpected section {keyword!r}")
    if not cur.rest_is_blank():
        raise ValueError(f"{len(blob) - cur.at} bytes left over after the declared sections")
    return out


def _cells(flat: np.ndarray, count: int, keyword: str) -> list[np.ndarray]:
    cells: list[np.ndarray] = []
    at = 0
    for _ in range(count):
        if at >= len(flat):
            raise ValueError(f"{keyword} declares {count} cells; the data ran out")
        n = int(flat[at])
        cells.append(np.array(flat[at + 1 : at + 1 + n]))
        at += 1 + n
    if at != len(flat):
        raise ValueError(f"{keyword} size is {len(flat)}; its {count} cells use {at}")
    return cells


def _field_arrays(cur: _Cursor, count: int) -> dict[str, Array]:
    arrays: dict[str, Array] = {}
    for _ in range(count):
        words = cur.keyword_line()
        if words is None or len(words) != 4:
            raise ValueError(f"a field array line is `name components tuples type`; got {words}")
        name, components, tuples, type_name = words[0], int(words[1]), int(words[2]), words[3]
        if name in arrays:
            raise ValueError(f"field array {name!r} declared twice")
        if type_name == "string":
            values, raw = cur.strings(components * tuples)
            arrays[name] = Array(name, type_name, components, values, raw)
        else:
            arrays[name] = Array(
                name, type_name, components, cur.numbers(components * tuples, type_name)
            )
    return arrays


def _point_scalars(cur: _Cursor, out: VtkFile) -> None:
    """One `SCALARS` block per point array; the writer puts CELL_DATA after."""
    assert out.point_count is not None
    words = cur.keyword_line()
    if not words or words[0] != "SCALARS":
        raise ValueError(f"POINT_DATA must hold SCALARS; got {words}")
    name, type_name = words[1], words[2]
    components = int(words[3]) if len(words) > 3 else 1
    table = cur.keyword_line()
    if not table or table[0] != "LOOKUP_TABLE":
        raise ValueError(f"SCALARS {name} must be followed by LOOKUP_TABLE; got {table}")
    values = cur.numbers(components * out.point_count, type_name)
    out.point_scalars[name] = Array(name, type_name, components, values)


def _cell_data(cur: _Cursor, out: VtkFile) -> None:
    assert out.cell_count is not None
    while (words := cur.keyword_line()) is not None:
        keyword = words[0]
        out.cell_sections.append(keyword)
        if keyword == "SCALARS":
            name, type_name = words[1], words[2]
            components = int(words[3]) if len(words) > 3 else 1
            table = cur.keyword_line()
            if not table or table[0] != "LOOKUP_TABLE":
                raise ValueError(f"SCALARS {name} must be followed by LOOKUP_TABLE; got {table}")
            out.lookup_tables[name] = table[1]
            values = cur.numbers(components * out.cell_count, type_name)
            out.scalars[name] = Array(name, type_name, components, values)
        elif keyword == "FIELD":
            arrays = _field_arrays(cur, int(words[2]))
            for array in arrays.values():
                tuples = len(array.values) // array.components
                if tuples != out.cell_count:
                    raise ValueError(
                        f"cell array {array.name} has {tuples} tuples for {out.cell_count} cells"
                    )
            out.cell_fields[words[1]] = arrays
        else:
            raise ValueError(f"unexpected CELL_DATA section {keyword!r}")


def lines_as_array(vtk: VtkFile) -> np.ndarray:
    """The `LINES` cells as an `(E, 2)` array, refusing a polyline of any other length."""
    if any(len(cell) != 2 for cell in vtk.lines):
        raise ValueError("a LINES cell with other than two points")
    return np.array(vtk.lines, dtype=np.int64).reshape(-1, 2)


def polygons_as_array(vtk: VtkFile) -> np.ndarray:
    """The `POLYGONS` cells as a `(T, 3)` array, refusing anything but triangles."""
    if any(len(cell) != 3 for cell in vtk.polygons):
        raise ValueError("a POLYGONS cell that is not a triangle")
    return np.array(vtk.polygons, dtype=np.int64).reshape(-1, 3)
