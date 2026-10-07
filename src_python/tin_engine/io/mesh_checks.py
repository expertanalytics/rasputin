"""The mesh writers' shared checks: the text a file header may hold, and codes
that must fit 32 bits. One copy each, so `io/ply.py` and `io/vtk_legacy.py`
refuse the same values in the same words. Imports nothing first-party.

A PLY header and a legacy VTK string array are line-oriented ASCII. Any control
character, not just `\\n`, can forge a line: `\\r` terminates one for every
CRLF-tolerant reader. A non-ASCII character cannot be encoded at all, and
`--crs` is free text by ruling 5, so both are a documented `ValueError` a
caller meets, never a `UnicodeEncodeError` from inside a writer.
"""

from __future__ import annotations

from typing import Any

import numpy.typing as npt


def checked_ascii(value: str, what: str) -> bytes:
    """`value` as ASCII bytes; a control character (below U+0020, or U+007F)
    or a non-ASCII character is a ValueError naming `what` and the character."""
    bad = next((ch for ch in value if ch < " " or ch == "\x7f"), None)
    if bad is not None:
        raise ValueError(f"a {what} may not contain a control character; got {bad!r}")
    if not value.isascii():
        bad = next(ch for ch in value if not ch.isascii())
        raise ValueError(f"a {what} must be ASCII; got {bad!r} in {value!r}")
    return value.encode("ascii")


def check_int32(values: npt.NDArray[Any], name: str) -> None:
    """ValueError `<name> must fit int32` unless every value is in
    [-2**31, 2**31 - 1]; an empty array passes."""
    if len(values) and (values.min() < -(2**31) or values.max() >= 2**31):
        raise ValueError(f"{name} must fit int32")


__all__ = ["check_int32", "checked_ascii"]
