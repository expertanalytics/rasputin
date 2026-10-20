"""GML 3.2 inputs: fixture A, cut from NIBIO's AR50 file, and builders for hand-made documents.

Fixture A, ``tests/fixtures/ar50/ar50_sample.gml``, is the real file's first
720 bytes (the declaration and the root element with its namespaces), an XML
comment carrying NIBIO's NLOD credit, six whole ``gml:featureMember``
elements copied byte for byte (the first five of the file and its feature
12,843, which has two holes), then the closing root tag. `AR50_FIGURES` is
what `docs/increments/73-probes/pick.py` measured on those six.

The builders write the small documents the decoder and converter suites
refuse or map: one string per case, a few hundred bytes each.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import NamedTuple

AR50_SAMPLE = Path(__file__).resolve().parents[1] / "fixtures" / "ar50" / "ar50_sample.gml"

#: The srsName on every AR50 polygon: ETRS89 / UTM 32N.
URI_25832 = "http://www.opengis.net/def/crs/EPSG/0/25832"
#: WGS 84 as an OGC URI: latitude first, as GML 3.2 writes it.
URI_4326 = "http://www.opengis.net/def/crs/EPSG/0/4326"

#: The leaf properties of an AR50 feature, in document order; the nested
#: ``identifikasjon``, ``kopidata`` and ``kvalitet`` are not among them.
AR50_VALUES = (
    "informasjon", "oppdateringsdato", "arealtype", "skogbonitet",
    "treslag", "jordbruk", "vegetasjonsdekke", "kartstandard",
)  # fmt: skip


class Ar50Figure(NamedTuple):
    gml_id: str
    holes: int
    exterior_positions: int  # the closing position included
    interior_positions: int
    area_m2: float  # rounded to 0.1 m2
    arealtype: str
    treslag: str
    skogbonitet: str


AR50_FIGURES = (
    Ar50Figure("id551e0bc0-b63d-45b7-a7ba-0ecdff9d7c95", 0, 37, 0, 15_315.0, "60", "32", "11"),
    Ar50Figure("idd58d2d03-6841-4a8e-8043-10e00987a82a", 0, 10, 0, 815.1, "60", "31", "13"),
    Ar50Figure("idfc6f8e2f-9c4b-445a-b840-1693a8483aec", 0, 15, 0, 4_899.6, "60", "31", "11"),
    Ar50Figure("id9dd0239a-c844-4fc2-a42a-1a6e70666a08", 0, 50, 0, 15_966.2, "60", "31", "11"),
    Ar50Figure("id40a10b70-c9f7-4ac4-b74b-30bd33ed551d", 0, 55, 0, 27_579.9, "60", "31", "11"),
    Ar50Figure("id0ed7f868-8d50-4dc3-bb29-53b1d9ecff7d", 2, 192, 60, 353_418.5, "30", "31", "18"),
)

GML32_NS = "http://www.opengis.net/gml/3.2"
HEAD = (
    '<?xml version="1.0" encoding="utf-8"?>\n'
    f'<gml:FeatureCollection xmlns:gml="{GML32_NS}" xmlns:app="urn:test"\n'
    '    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
)
TAIL = "</gml:FeatureCollection>\n"

SQUARE = "0 0 10 0 10 10 0 10 0 0"
SQUARE_2 = "20 0 30 0 30 10 20 10 20 0"
HOLE = "2 2 2 4 4 4 4 2 2 2"


def document(*features: str) -> io.BytesIO:
    """A GML 3.2 feature collection, one ``gml:featureMember`` per feature."""
    members = "".join(f"  <gml:featureMember>{f}</gml:featureMember>\n" for f in features)
    return io.BytesIO((HEAD + members + TAIL).encode())


def values(**pairs: str | None) -> str:
    """Leaf properties in the order given; None is ``xsi:nil``."""
    return "".join(
        f'<app:{k} xsi:nil="true"/>' if v is None else f"<app:{k}>{v}</app:{k}>"
        for k, v in pairs.items()
    )


def feature(gml_id: str, *geometries: str, props: str | None = None, type_name: str = "F") -> str:
    """An ``app:<type_name>`` with leaf ``props`` (default ``a`` = 60), then one
    ``app:g`` property per geometry."""
    body = values(a="60") if props is None else props
    geometry = "".join(f"<app:g>{g}</app:g>" for g in geometries)
    return f'<app:{type_name} gml:id="{gml_id}">{body}{geometry}</app:{type_name}>'


def ring(positions: str, dimension: int | None = None) -> str:
    attribute = "" if dimension is None else f' srsDimension="{dimension}"'
    return f"<gml:LinearRing><gml:posList{attribute}>{positions}</gml:posList></gml:LinearRing>"


def rings(exterior: str, *holes: str) -> str:
    inner = "".join(f"<gml:interior>{h}</gml:interior>" for h in holes)
    return f"<gml:exterior>{exterior}</gml:exterior>{inner}"


def srs(name: str | None) -> str:
    return "" if name is None else f' srsName="{name}"'


def polygon(exterior: str = SQUARE, *holes: str, srs_name: str | None = URI_25832) -> str:
    """A ``gml:Polygon`` of posList rings in two dimensions."""
    body = rings(ring(exterior), *(ring(h) for h in holes))
    return f"<gml:Polygon{srs(srs_name)}>{body}</gml:Polygon>"


def surface(*patches: str, srs_name: str | None = URI_25832) -> str:
    """A ``gml:Surface`` of one ``gml:PolygonPatch`` per exterior ring given."""
    body = "".join(f"<gml:PolygonPatch>{rings(ring(p))}</gml:PolygonPatch>" for p in patches)
    return f"<gml:Surface{srs(srs_name)}><gml:patches>{body}</gml:patches></gml:Surface>"
