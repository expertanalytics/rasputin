"""A cut run's index, `<out>.pieces/index.json`, and K4's check before it is written.

Increment 23c, `docs/increments/23-basin-scale.md`, "Output" and "Settled
after 23c's red step" (items 3, 5, 6 and 11, question B). `MeshIndex` is
frozen and refuses an unknown key: it is read by other programs (23d's
stitcher, a consumer), so a key it does not know is drift, not data.

`check_conformity` is K4: the records two pieces wrote of one seam edge are
equal bit for bit, or the run fails naming the seam and the two pieces. Bit
for bit means bytes, so a NaN z equals a NaN z and 0.0 differs from -0.0.
One record alone is a seam edge on the outline, which has one piece.

Like the rest of `io/`, this module imports nothing first-party but
`io/models.py`, and never `_core`.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Annotated, Self

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator

from .models import IndexWindow

PieceId = tuple[int, int, int]  # (cell row j, cell column i, component rank k)
_Rank = Annotated[int, Field(strict=True, ge=0)]  # no text, no float, no sign


class _Record(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


class PartitionRecord(_Record):
    """The rule's inputs and its cells, as `decompose.partition` gave them."""

    pieces: int = Field(ge=1)
    memory_budget: int = Field(ge=1)
    bytes_per_node: int = Field(ge=1)
    nx: int = Field(ge=1)
    ny: int = Field(ge=1)
    dx: int = Field(ge=1)
    dy: int = Field(ge=1)
    cols: int = Field(ge=1)
    rows: int = Field(ge=1)


class PieceWindow(IndexWindow):
    """A piece's window on the lattice; unknown keys refused, as in the index."""

    model_config = ConfigDict(frozen=True, extra="forbid")


class PieceCounts(_Record):
    """The piece's own counts, never the whole run's."""

    triangles: int = Field(ge=0)
    vertices: int = Field(ge=0)
    on_frozen: int = Field(ge=0)


class PieceEntry(_Record):
    """One piece. `file` and `sha256` are None for a piece all over NoData."""

    id: tuple[_Rank, _Rank, _Rank]
    file: str | None
    sha256: str | None
    counts: PieceCounts
    window: PieceWindow

    @model_validator(mode="after")
    def _file_and_sha256_together(self) -> Self:
        if (self.file is None) != (self.sha256 is None):
            raise ValueError("file and sha256 are null together or set together")
        return self


class MeshIndex(_Record):
    """`index.json`: the run, its partition and its pieces in piece-id order."""

    crs: str
    tolerance: float
    source: str  # the run's elevation_source sentence: identity and credit
    vocabulary: str  # the piece files' vocabulary fingerprint
    partition: PartitionRecord
    pieces: tuple[PieceEntry, ...]


@dataclass(frozen=True, slots=True)
class SeamRecord:
    """A piece's vertex sequence along one seam edge, from the edge's lower end.

    `edge` is the seam edge's index in the start triangulation's edge list.
    """

    piece: PieceId
    edge: int
    xyz: tuple[tuple[float, float, float], ...]


class ConformityError(ValueError):
    """Two pieces wrote one seam edge differently (K4)."""

    def __init__(self, edge: int, pieces: tuple[PieceId, PieceId]) -> None:
        super().__init__(
            f"seam edge {edge} differs between pieces {pieces[0]} and {pieces[1]}: "
            "the pieces do not conform"
        )
        self.edge = edge
        self.pieces = pieces


def check_conformity(records: Iterable[SeamRecord]) -> None:
    """Raise `ConformityError` on the first seam edge whose records differ in a byte."""
    first: dict[int, tuple[SeamRecord, bytes]] = {}
    for record in records:
        data = np.asarray(record.xyz, dtype=np.float64).tobytes()
        seen = first.setdefault(record.edge, (record, data))
        if seen[1] != data:
            raise ConformityError(record.edge, (seen[0].piece, record.piece))
