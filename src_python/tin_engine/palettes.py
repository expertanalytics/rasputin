"""A ParaView preset from a class map's palette.

The palettes themselves are data: each coded class-map file in
``tin_engine/data/class_maps/`` carries its classes' labels and colours.
CORINE's are natural colours after Patterson and Kelso (2004), not the
official CLC legend: a choice of taste, free to change.
"""

from __future__ import annotations

from collections.abc import Mapping

#: A code the table lacks draws magenta, so it is seen rather than blended in.
_NAN_COLOR = [1.0, 0.0, 1.0]


def paraview_preset(table: Mapping[int, tuple[str, str]], name: str) -> list[dict[str, object]]:
    """A categorical ParaView colour preset: `Annotations` (value, label, ...)
    and `IndexedColors` (r, g, b in 0 .. 1, one triple per value, in order)."""
    annotations: list[str] = []
    colours: list[float] = []
    for code, (label, colour) in sorted(table.items()):
        annotations += [str(code), f"{code} {label}"]
        colours += [int(colour[i : i + 2], 16) / 255 for i in (1, 3, 5)]
    return [
        {
            "Name": name,
            "Annotations": annotations,
            "IndexedColors": colours,
            "NanColor": _NAN_COLOR,
        }
    ]
