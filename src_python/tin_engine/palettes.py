"""Natural colours for land-cover codes, and a ParaView preset from them.

Increment 16c (`docs/increments/16c-landcover-labels.md`, R4). Data plus one
function; nothing first-party, nothing about meshes. The colours follow the
natural-colour convention of Patterson and Kelso (2004): the colour the ground
has from above, not the official CLC legend (which paints peat bogs blue and
the sea almost white). They are taste, and Ola's to change (Default D6).
"""

from __future__ import annotations

#: CORINE level-3 code -> (CLC `LABEL3`, `#rrggbb`); 0 is no polygon, and
#: every constraint line in a `.vtk`.
# fmt: off
CORINE_NATURAL: dict[int, tuple[str, str]] = {
    0: ("No polygon; constraint lines", "#3c3c3c"),
    111: ("Continuous urban fabric", "#8c5f5a"),
    112: ("Discontinuous urban fabric", "#a88d86"),
    121: ("Industrial or commercial units", "#8e8a96"),
    122: ("Road and rail networks and associated land", "#6e6e6e"),
    123: ("Port areas", "#7d8796"),
    124: ("Airports", "#a8a8a8"),
    131: ("Mineral extraction sites", "#b49b78"),
    132: ("Dump sites", "#8a7d64"),
    133: ("Construction sites", "#bcae98"),
    141: ("Green urban areas", "#86b86e"),
    142: ("Sport and leisure facilities", "#a3cf7e"),
    211: ("Non-irrigated arable land", "#e6d58c"),
    212: ("Permanently irrigated land", "#d9cc6e"),
    213: ("Rice fields", "#cdd89a"),
    221: ("Vineyards", "#9c7a44"),
    222: ("Fruit trees and berry plantations", "#a9ad5e"),
    223: ("Olive groves", "#8f9a52"),
    231: ("Pastures", "#b4d47a"),
    241: ("Annual crops associated with permanent crops", "#dccf94"),
    242: ("Complex cultivation patterns", "#d2c47c"),
    243: ("Land principally occupied by agriculture, with significant areas of natural vegetation", "#bcc47e"),  # noqa: E501
    244: ("Agro-forestry areas", "#a9b574"),
    311: ("Broad-leaved forest", "#4f8f3f"),
    312: ("Coniferous forest", "#1f5a2e"),
    313: ("Mixed forest", "#357438"),
    321: ("Natural grasslands", "#c2d68a"),
    322: ("Moors and heathland", "#a7c47f"),
    323: ("Sclerophyllous vegetation", "#8a9658"),
    324: ("Transitional woodland-shrub", "#7ea65a"),
    331: ("Beaches, dunes, sands", "#e9ddb2"),
    332: ("Bare rocks", "#8f8f8f"),
    333: ("Sparsely vegetated areas", "#b9b8a0"),
    334: ("Burnt areas", "#4b3f3a"),
    335: ("Glaciers and perpetual snow", "#eef6fb"),
    411: ("Inland marshes", "#6e9470"),
    412: ("Peat bogs", "#8a6642"),
    421: ("Salt marshes", "#7f9f8c"),
    422: ("Salines", "#d8d6cc"),
    423: ("Intertidal flats", "#b3bdb3"),
    511: ("Water courses", "#4c8ec4"),
    512: ("Water bodies", "#3e7bb6"),
    521: ("Coastal lagoons", "#5b93b3"),
    522: ("Estuaries", "#5188b4"),
    523: ("Sea and ocean", "#2b5d8e"),
}
# fmt: on

#: The palettes `rasputin palette` knows, by name, with their preset names.
PALETTES = {"corine": (CORINE_NATURAL, "rasputin CORINE natural")}

#: A code the table lacks draws magenta, so it is seen rather than blended in.
_NAN_COLOR = [1.0, 0.0, 1.0]


def paraview_preset(table: dict[int, tuple[str, str]], name: str) -> list[dict[str, object]]:
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
