"""River segments (NVE's ELVIS network or a user's own), read back (increment 29, PR 3).

`docs/increments/29-nve-reference-catchments.md`, "Placing the gauge" (the
model, the `kind` rule and step 4, exact copies) and "The station set". One
LineString per segment: any other geometry, a MultiLineString included, is
refused naming the segment's `objectid`, and so is an `objectid` that appears
twice. Exact copies within one `elvid` are dropped by :func:`drop_copies`,
so a user's file is cleaned as a fetched one is, and the count dropped is
returned for the commands to print.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict

from .station_set import _geometry_type, _unique, features_of, required

#: Two vertices closer than this, in metres, are the same vertex (step 4).
COPY_TOLERANCE = 0.01


class RiverSegment(BaseModel):
    """One mapped line, digitised downstream, in its file's CRS. `objekttype`
    is kept as served; `kind` is what "Placing the gauge" reads."""

    model_config = ConfigDict(frozen=True)

    objectid: int
    elvid: str | None
    vassdragsnr: str | None
    name: str | None
    objekttype: str | None
    kind: Literal["lake", "river"]
    line: tuple[tuple[float, float], ...]


def kind_of(objekttype: str | None, vatnlnr: Any) -> Literal["lake", "river"]:
    """Total: `lake` when the casefolded `objekttype` starts with `innsj`, or
    when it is null or blank and `vatnlnr` is set: not null and not 0, which
    is how NVE says "no lake"; `river` otherwise."""
    if objekttype is None or not objekttype.strip():
        return "lake" if vatnlnr is not None and vatnlnr != 0 else "river"
    return "lake" if objekttype.casefold().startswith("innsj") else "river"


def _segment(path: Path, feature: Mapping[str, Any]) -> RiverSegment:
    props = feature["properties"]
    try:
        line = tuple((float(x), float(y)) for x, y, *_ in feature["geometry"]["coordinates"])
    except (TypeError, ValueError) as exc:
        oid = props["objectid"]
        raise ValueError(f"{path.name}: segment {oid}: a vertex is not two numbers") from exc
    return RiverSegment(
        objectid=props["objectid"],
        elvid=props.get("elvid"),
        vassdragsnr=props.get("vassdragsnr"),
        name=props.get("elvenavn"),
        objekttype=props.get("objekttype"),
        kind=kind_of(props.get("objekttype"), props.get("vatnlnr")),
        line=line,
    )


def _same(a: RiverSegment, b: RiverSegment) -> bool:
    return len(a.line) == len(b.line) and all(
        math.dist(p, q) <= COPY_TOLERANCE for p, q in zip(a.line, b.line, strict=True)
    )


def drop_copies(
    segments: Iterable[RiverSegment],
) -> tuple[tuple[RiverSegment, ...], int]:
    """`segments` without exact copies (within one `elvid`, vertex lists equal
    to 1 cm), the smallest `objectid` kept, in input order; and how many went."""
    given = tuple(segments)
    kept: dict[str | None, list[RiverSegment]] = {}
    for s in sorted(given, key=lambda s: s.objectid):
        group = kept.setdefault(s.elvid, [])
        if not any(_same(s, k) for k in group):
            group.append(s)
    keep = {s.objectid for group in kept.values() for s in group}
    out = tuple(s for s in given if s.objectid in keep)
    return out, len(given) - len(out)


def read_segments(path: Path) -> tuple[tuple[RiverSegment, ...], str, int]:
    """The segments in `path` with exact copies dropped, the file's CRS, and
    the count of copies dropped."""
    features, crs = features_of(path)
    for f in features:
        _geometry_type(path, f, ("LineString",))
    _unique(path, [required(path, f, "objectid") for f in features], "objectid")
    segments, dropped = drop_copies(_segment(path, f) for f in features)
    return segments, crs, dropped


__all__ = ["COPY_TOLERANCE", "RiverSegment", "drop_copies", "kind_of", "read_segments"]
