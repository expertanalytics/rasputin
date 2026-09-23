"""An independent PLY reader. Test support only; nothing in `src_python` uses it.

This module exists so that increment 10's round trip goes through a parser that
**reads the header it was handed** instead of re-deriving the writer's byte
offsets (`PRINCIPLES.md` B1, and `10-mesh-output.md`'s "Testing" section). The
distinction is the whole value of the suite: a reader that computed "vertices
start at 3 * 8 * N bytes from the end of the header, because that is what
`write_ply` computes" would agree with the writer's arithmetic bug and report a
pass.

So every size here comes from a declaration in the blob:

- which elements exist, and how many of each, comes from the `element` lines;
- how wide each record is comes from summing the declared property types;
- where the body starts comes from locating `end_header`.

Nothing is read from `tin_engine`. The only thing assumed about the input is
that it is PLY, which is the thing under test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_END = b"end_header\n"

#: PLY type names mapped onto little-endian numpy types. Both the classic
#: names (`uchar`, `int`, `double`) and the explicit-width aliases PLY also
#: permits (`uint8`, `int32`, `float64`) are accepted, because which spelling
#: the writer chooses is not this suite's business -- the widths are.
_DTYPES = {
    "char": "i1",
    "int8": "i1",
    "uchar": "u1",
    "uint8": "u1",
    "short": "<i2",
    "int16": "<i2",
    "ushort": "<u2",
    "uint16": "<u2",
    "int": "<i4",
    "int32": "<i4",
    "uint": "<u4",
    "uint32": "<u4",
    "float": "<f4",
    "float32": "<f4",
    "double": "<f8",
    "float64": "<f8",
}


@dataclass(frozen=True)
class Property:
    """One declared property. `count_type` is None unless it is a list."""

    name: str
    type_name: str
    count_type: str | None = None

    @property
    def is_list(self) -> bool:
        return self.count_type is not None

    @property
    def dtype(self) -> np.dtype[np.generic]:
        return np.dtype(_DTYPES[self.type_name])


@dataclass(frozen=True)
class Element:
    """One declared element block: its name, its record count, its properties."""

    name: str
    count: int
    properties: tuple[Property, ...]

    @property
    def record_size(self) -> int:
        """Bytes per record. Only meaningful when no property is a list."""
        if any(p.is_list for p in self.properties):
            raise ValueError(f"element {self.name} has a list property; records vary in size")
        return sum(p.dtype.itemsize for p in self.properties)


@dataclass(frozen=True)
class Header:
    """What the header said, and where it stopped."""

    fmt: str
    version: str
    comments: tuple[str, ...]
    elements: tuple[Element, ...]
    body_offset: int

    def element(self, name: str) -> Element:
        for element in self.elements:
            if element.name == name:
                return element
        raise KeyError(f"no element {name!r}; the file declares {[e.name for e in self.elements]}")

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(element.name for element in self.elements)


def parse_header(blob: bytes) -> Header:
    """Parse the ASCII header of `blob`, refusing anything that is not PLY."""
    if not blob.startswith(b"ply\n"):
        raise ValueError(f"not a PLY file: it begins {blob[:16]!r}")
    stop = blob.find(_END)
    if stop < 0:
        raise ValueError("no end_header line")
    body_offset = stop + len(_END)
    lines = blob[: stop - 1].decode("ascii").split("\n")

    fmt = version = ""
    comments: list[str] = []
    elements: list[Element] = []
    properties: list[Property] = []
    for line in lines[1:]:
        word, _, rest = line.partition(" ")
        if word == "format":
            fmt, _, version = rest.partition(" ")
        elif word == "comment":
            comments.append(rest)
        elif word == "element":
            if elements:
                elements[-1] = _with_properties(elements[-1], properties)
            name, _, count = rest.partition(" ")
            elements.append(Element(name=name, count=int(count), properties=()))
            properties = []
        elif word == "property":
            properties.append(_parse_property(rest))
        else:
            raise ValueError(f"unexpected header line {line!r}")
    if elements:
        elements[-1] = _with_properties(elements[-1], properties)
    return Header(
        fmt=fmt,
        version=version,
        comments=tuple(comments),
        elements=tuple(elements),
        body_offset=body_offset,
    )


def _with_properties(element: Element, properties: list[Property]) -> Element:
    return Element(name=element.name, count=element.count, properties=tuple(properties))


def _parse_property(rest: str) -> Property:
    words = rest.split()
    if words[0] == "list":
        _, count_type, type_name, name = words
        return Property(name=name, type_name=type_name, count_type=count_type)
    type_name, name = words
    return Property(name=name, type_name=type_name)


def read_ply(blob: bytes) -> tuple[Header, dict[str, dict[str, np.ndarray]]]:
    """Parse a whole PLY file into `{element: {property: array}}`.

    Raises if the body does not end exactly where the declarations say it does.
    A file with trailing bytes, or one record short, is a defect the round trip
    must see rather than tolerate: both are what a wrong offset looks like from
    the outside.
    """
    header = parse_header(blob)
    if header.fmt == "binary_little_endian":
        data, consumed = _read_binary(header, blob)
    elif header.fmt == "ascii":
        data, consumed = _read_ascii(header, blob)
    else:
        raise ValueError(f"unsupported format {header.fmt!r}")
    if consumed != len(blob):
        raise ValueError(f"body is {len(blob) - header.body_offset} bytes, declarations want "
                         f"{consumed - header.body_offset}")
    return header, data


def _read_binary(header: Header, blob: bytes) -> tuple[dict[str, dict[str, np.ndarray]], int]:
    offset = header.body_offset
    data: dict[str, dict[str, np.ndarray]] = {}
    for element in header.elements:
        if any(p.is_list for p in element.properties):
            data[element.name], offset = _read_binary_lists(element, blob, offset)
            continue
        dtype = np.dtype([(p.name, p.dtype) for p in element.properties])
        needed = dtype.itemsize * element.count
        if offset + needed > len(blob):
            raise ValueError(f"element {element.name} runs past the end of the file")
        records = np.frombuffer(blob, dtype=dtype, count=element.count, offset=offset)
        data[element.name] = {p.name: np.array(records[p.name]) for p in element.properties}
        offset += needed
    return data, offset


def _read_binary_lists(
    element: Element, blob: bytes, offset: int
) -> tuple[dict[str, np.ndarray], int]:
    """Read a block record by record, because list properties vary in length."""
    rows: dict[str, list[np.ndarray]] = {p.name: [] for p in element.properties}
    for _ in range(element.count):
        for prop in element.properties:
            if prop.is_list:
                count_dtype = np.dtype(_DTYPES[str(prop.count_type)])
                (length,) = np.frombuffer(blob, dtype=count_dtype, count=1, offset=offset)
                offset += count_dtype.itemsize
                values = np.frombuffer(blob, dtype=prop.dtype, count=int(length), offset=offset)
                offset += prop.dtype.itemsize * int(length)
                rows[prop.name].append(np.array(values))
            else:
                (value,) = np.frombuffer(blob, dtype=prop.dtype, count=1, offset=offset)
                offset += prop.dtype.itemsize
                rows[prop.name].append(np.array(value))
    return {name: np.array(values) for name, values in rows.items()}, offset


def _read_ascii(header: Header, blob: bytes) -> tuple[dict[str, dict[str, np.ndarray]], int]:
    body = blob[header.body_offset :].decode("ascii")
    lines = body.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    data: dict[str, dict[str, np.ndarray]] = {}
    cursor = 0
    for element in header.elements:
        rows: dict[str, list[object]] = {p.name: [] for p in element.properties}
        for _ in range(element.count):
            tokens = lines[cursor].split()
            cursor += 1
            at = 0
            for prop in element.properties:
                if prop.is_list:
                    length = int(tokens[at])
                    at += 1
                    rows[prop.name].append([prop.dtype.type(t) for t in tokens[at : at + length]])
                    at += length
                else:
                    rows[prop.name].append(prop.dtype.type(tokens[at]))
                    at += 1
            if at != len(tokens):
                raise ValueError(f"element {element.name}: {len(tokens) - at} tokens left over")
        data[element.name] = {
            p.name: np.array(rows[p.name], dtype=p.dtype) for p in element.properties
        }
        cursor += 0
    if cursor != len(lines):
        raise ValueError(f"{len(lines) - cursor} lines left over after the declared elements")
    return data, len(blob)


def element_bytes(blob: bytes, name: str) -> bytes:
    """The raw bytes of one fixed-width element block, located from the header.

    Used for the byte-identity claim ruling 3 rests on: the surface file's
    vertex block and the constraint file's must be the same bytes, or vertex
    *i* does not mean the same point in both files.
    """
    header = parse_header(blob)
    offset = header.body_offset
    for element in header.elements:
        size = element.record_size * element.count
        if element.name == name:
            return blob[offset : offset + size]
        offset += size
    raise KeyError(f"no element {name!r} in {header.names}")


def vertex_array(data: dict[str, dict[str, np.ndarray]]) -> np.ndarray:
    """The vertex block as an (N, 3) array, from the x/y/z properties."""
    vertex = data["vertex"]
    return np.column_stack([vertex["x"], vertex["y"], vertex["z"]])
