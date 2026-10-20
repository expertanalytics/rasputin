"""Reads back the binary legacy VTK that `rasputin mesh` writes: bytes in, arrays out.

Only rasputin's own files (second header line `rasputin mesh`), binary only:
the dataset strings `crs` and `land_cover_codes`, the points, the triangles
and the cell array `land_cover_code`; constraint lines and every other array
are skipped. Pure: no path, no file.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

#: Bytes per value of each legacy type the writer uses.
_SIZE = {"double": 8, "float": 4, "int": 4, "unsigned_int": 4, "unsigned_char": 1}


@dataclass(frozen=True)
class MeshArrays:
    """`vertices` `(n, 3)` float64, `triangles` `(T, 3)` int64, `crs` the text written;
    `land_cover_code` `(T,)` int32, one per triangle, or None when the file has
    none, and `land_cover_codes` the text saying what the codes are, or ""."""

    vertices: npt.NDArray[np.float64]
    triangles: npt.NDArray[np.int64]
    crs: str
    land_cover_code: npt.NDArray[np.int32] | None = None
    land_cover_codes: str = ""


class _Cursor:
    def __init__(self, data: bytes) -> None:
        self.data, self.at = data, 0

    def take(self, size: int) -> bytes:
        if self.at + size > len(self.data):
            raise ValueError("the mesh file is truncated")
        self.at += size
        return self.data[self.at - size : self.at]

    def line(self) -> list[str]:
        end = self.data.find(b"\n", self.at)
        if end < 0:
            raise ValueError("the mesh file is truncated")
        text = self.take(end + 1 - self.at)[:-1]
        return text.decode("ascii", errors="replace").split()

    def string(self) -> str:
        head = self.take(1)[0]
        width = {3: 1, 2: 2, 1: 4, 0: 8}[head >> 6]
        raw = bytes([head & 0x3F]) + self.take(width - 1)
        return self.take(int.from_bytes(raw, "big")).decode("ascii")


def read_mesh(data: bytes) -> MeshArrays:
    """The mesh in `data`. A ValueError for an ASCII file, another program's
    file, a file without a `crs` string, a truncated one, or a `land_cover_code`
    array that is not `int` or not one entry per cell."""
    head = data.split(b"\n", 3)
    if len(head) < 4 or not head[0].startswith(b"# vtk DataFile") or head[1] != b"rasputin mesh":
        raise ValueError("not a mesh written by rasputin (its second line is not 'rasputin mesh')")
    if head[2] != b"BINARY":
        raise ValueError("an ASCII mesh cannot be read here: write the mesh without --ascii")
    c = _Cursor(data)
    for _ in range(4):
        c.line()
    strings: dict[str, str] = {}
    for _ in range(int(c.line()[2])):
        name, _, count, kind = c.line()
        if kind == "string":
            values = [c.string() for _ in range(int(count))]
            if values:
                strings[name] = values[0]
        else:
            c.take(int(count) * _SIZE[kind])
        c.take(1)
    crs = strings.get("crs")
    if crs is None:
        raise ValueError("the mesh has no crs string; write it with rasputin mesh")
    _, n, _ = c.line()
    vertices = np.frombuffer(c.take(int(n) * 24), ">f8").reshape(-1, 3).astype(np.float64)
    c.take(1)
    keyword, cells, size = c.line()
    lines = 0
    if keyword == "LINES":
        lines = int(cells)
        c.take(int(size) * 4 + 1)
        keyword, cells, size = c.line()
    flat = np.frombuffer(c.take(int(size) * 4), ">i4").reshape(int(cells), -1)
    if keyword != "POLYGONS" or flat.shape[1] != 4 or (flat[:, 0] != 3).any():
        raise ValueError("the mesh's cells are not all triangles")
    c.take(1)
    codes, total = None, lines + len(flat)
    while c.at < len(c.data):
        words = c.line()
        if words[0] in ("POINT_DATA", "CELL_DATA"):
            values_n = int(words[1])
        elif words[0] == "SCALARS":
            c.line()
            c.take(values_n * _SIZE[words[2]] + 1)
        elif words[0] == "FIELD":
            for _ in range(int(words[2])):
                name, _, entries, kind = c.line()
                body = c.take(int(entries) * _SIZE[kind] + 1)
                if name != "land_cover_code":
                    continue
                if kind != "int":
                    raise ValueError(f"the mesh's land_cover_code array is {kind}, not int")
                if int(entries) != total:
                    msg = (
                        f"the mesh's land_cover_code array has {entries} entries for {total} cells"
                    )
                    raise ValueError(msg)
                codes = np.frombuffer(body[:-1], ">i4")[lines:].astype(np.int32)
    return MeshArrays(
        vertices, flat[:, 1:].astype(np.int64), crs, codes, strings.get("land_cover_codes", "")
    )
