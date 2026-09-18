"""The synthetic gallery: eight shapes, each answering one question.

Declarative data. Nothing here is generated, because a fixture whose
coordinates are computed is a fixture nobody can check by reading -- every
number below was typed, and the only arithmetic is :data:`ORIGIN`, a constant
translation that puts the whole gallery at UTM 33N magnitudes. That translation
is the ordering ruling's one mitigation: the robustness behaviour a person is
being asked to judge is magnitude-sensitive, and a picture drawn in the unit
square silently exercises the easy case. Coordinates are written as local
metres so the shapes stay legible on the page.

Roles are the **strings** ``"outer"``, ``"hole"`` and ``"breakline"``. That is
forced rather than chosen: ``viz/`` never imports ``_core``, so a fixture cannot
name ``ChainRole``, and ``cli.py`` -- the composition root -- maps the strings
onto the enum before calling ``build_pslg``. Keeping the fixture itself a
``PslgLike`` is what lets a fixture the validator *rejects* still be drawn;
``degenerate`` is exactly that case and has no ``Pslg`` to hand over.

Ring winding is the validator's, not a preference: an ``outer`` ring must be
counter-clockwise and a ``hole`` clockwise, and neither stores its closing
index. **Three** of the eight are deliberate failures, and they fail at two
different depths: ``degenerate`` never reaches ``triangulate`` at all, because
the PSLG validator refuses it (``PslgError.DegenerateRing``), while
``not-noded`` and ``hole-in-hole`` are backend refusals of a valid PSLG
(``NotNoded`` and ``InvalidTopology``). All three are here because a failure
presentation nobody has looked at is a failure presentation that is wrong.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

#: UTM 33N-shaped easting and northing added to every local coordinate.
ORIGIN = np.array([430_000.0, 6_900_000.0], dtype=np.float64)


@dataclass(frozen=True, slots=True)
class FixtureChain:
    """One chain, structurally a ``ChainLike``: a run of ``chain_indices``."""

    begin: int
    count: int
    role: str
    properties: int
    """The chain's feature set, as a BARE mask, written as a number.

    ``viz/`` imports nothing first-party here, so this module cannot ask a
    vocabulary what bit ``river`` is and the number below has to be the number
    ``tin_engine.features.DEFAULT_VOCABULARY`` gives it. That agreement is the
    producer/consumer risk in its smallest form, and
    ``test_the_river_fixture_sets_the_river_bit_and_only_that_bit`` is where it
    is checked."""


@dataclass(frozen=True, slots=True)
class Fixture:
    """One gallery entry, structurally a ``PslgLike``, plus its own reason.

    ``description`` is what ``draw --list`` prints. A fixture without one is a
    fixture whose reason for existing has been lost, which is how a gallery
    becomes a list of eight opaque nouns.
    """

    name: str
    description: str
    vertices: npt.NDArray[np.float64]
    chains: tuple[FixtureChain, ...]
    chain_indices: npt.NDArray[np.uint32]

    def indices_of(self, c: int) -> npt.NDArray[np.uint32]:
        """Chain ``c``'s slice of the flat buffer, as ``Pslg`` would give it."""
        chain = self.chains[c]
        return self.chain_indices[chain.begin : chain.begin + chain.count]


def _fixture(
    name: str,
    description: str,
    points: Sequence[Sequence[float]],
    chains: Sequence[tuple[Sequence[int], str, int]],
) -> Fixture:
    """Flatten authored index runs into the ``(begin, count)`` form.

    The only thing this computes is where each chain starts, which is the one
    number a reader should not be asked to keep in their head.
    """
    flat: list[int] = []
    records: list[FixtureChain] = []
    for indices, role, properties in chains:
        records.append(FixtureChain(len(flat), len(indices), role, properties))
        flat.extend(int(i) for i in indices)
    return Fixture(
        name=name,
        description=description,
        vertices=np.array(points, dtype=np.float64) + ORIGIN,
        chains=tuple(records),
        chain_indices=np.array(flat, dtype=np.uint32),
    )


#: An open polyline across the interior, shared by `breakline-chain` and
#: `river`, which differ only in the property mask -- written once so that the
#: two pictures are the same picture with one stroke changed.
_CROSSING = [[60.0, 80.0], [350.0, 260.0], [640.0, 120.0]]

_BOX_700_500 = [[0.0, 0.0], [700.0, 0.0], [700.0, 500.0], [0.0, 500.0]]

CATCHMENT = _fixture(
    "catchment",
    "an outer ring with two lake holes and an interior breakline -- the shape the project is for",
    [
        # The outline, counter-clockwise from the south-west.
        [0.0, 0.0], [300.0, 40.0], [620.0, 0.0], [900.0, 180.0], [960.0, 520.0],
        [740.0, 760.0], [420.0, 900.0], [120.0, 760.0], [20.0, 480.0], [60.0, 220.0],
        # Lake one, clockwise.
        [200.0, 380.0], [300.0, 340.0], [260.0, 250.0], [160.0, 280.0],
        # Lake two, clockwise.
        [640.0, 640.0], [760.0, 600.0], [720.0, 500.0], [600.0, 520.0],
        # A breakline through the northern half, clear of both lakes.
        [120.0, 560.0], [340.0, 620.0], [520.0, 760.0],
    ],
    [
        (range(10), "outer", 0),
        ([10, 11, 12, 13], "hole", 0),
        ([14, 15, 16, 17], "hole", 0),
        ([18, 19, 20], "breakline", 0),
    ],
)

SLIVER_FAN = _fixture(
    "sliver-fan",
    "a fan of near-collinear constraints; what the triangulator does with extreme aspect ratios",
    [
        [0.0, 0.0], [1000.0, 0.0], [1000.0, 200.0], [0.0, 200.0],
        # The apex, shared by all three rays, and their far ends. The first ray
        # is a three-point polyline that is near-collinear without being
        # collinear -- 1000 m long and 1 m off straight.
        [10.0, 100.0], [500.0, 101.0], [990.0, 100.5], [990.0, 150.0], [990.0, 55.0],
    ],
    [
        (range(4), "outer", 0),
        ([4, 5, 6], "breakline", 0),
        ([4, 7], "breakline", 0),
        ([4, 8], "breakline", 0),
    ],
)

CORNER_HOLE = _fixture(
    "corner-hole",
    "a hole touching the outer ring at exactly one vertex -- InvalidTopology's neighbour",
    [
        [0.0, 0.0], [600.0, 0.0], [600.0, 600.0], [0.0, 600.0],
        [560.0, 400.0], [400.0, 560.0],
    ],
    [
        (range(4), "outer", 0),
        # Index 2 is the ring's north-east corner, reused: one shared vertex is
        # a topology worth looking at, two would be a slit and another fixture.
        ([2, 4, 5], "hole", 0),
    ],
)

HOLE_IN_HOLE = _fixture(
    "hole-in-hole",
    "a hole nested inside a second outline; what in-domain means, drawn",
    [
        [0.0, 0.0], [800.0, 0.0], [800.0, 800.0], [0.0, 800.0],
        [150.0, 650.0], [650.0, 650.0], [650.0, 150.0], [150.0, 150.0],
        [300.0, 500.0], [500.0, 500.0], [500.0, 300.0], [300.0, 300.0],
    ],
    [
        (range(4), "outer", 0),
        ([4, 5, 6, 7], "hole", 0),
        ([8, 9, 10, 11], "hole", 0),
    ],
)

BREAKLINE_CHAIN = _fixture(
    "breakline-chain",
    "an open breakline crossing the interior; where the Delaunay property visibly stops",
    _BOX_700_500 + _CROSSING,
    [(range(4), "outer", 0), ([4, 5, 6], "breakline", 0)],
)

RIVER = _fixture(
    "river",
    "the same open breakline carrying the river property, so that stroke is exercised early",
    _BOX_700_500 + _CROSSING,
    # Bit 0, spelled as the bare number a vocabulary would give for "river":
    # this module holds no vocabulary and may not import one.
    [(range(4), "outer", 0), ([4, 5, 6], "breakline", 1)],
)

NOT_NODED = _fixture(
    "not-noded",
    "two crossing constraints -- a deliberate failure, rendering the non-Ok presentation",
    [
        [0.0, 0.0], [700.0, 0.0], [700.0, 700.0], [0.0, 700.0],
        [100.0, 100.0], [600.0, 600.0], [100.0, 600.0], [600.0, 100.0],
    ],
    [
        (range(4), "outer", 0),
        ([4, 5], "breakline", 0),
        # Crosses the chain above at (350, 350), a point neither chain names:
        # that is what "not noded" means, and the noder (5b) is what fixes it.
        ([6, 7], "breakline", 0),
    ],
)

DEGENERATE = _fixture(
    "degenerate",
    "an all-collinear outer ring -- refused by the PSLG validator, drawn with its own words",
    [[0.0, 0.0], [200.0, 0.0], [400.0, 0.0], [600.0, 0.0]],
    [(range(4), "outer", 0)],
)

#: The gallery, by name. `draw --list` prints it; `draw NAME` looks a fixture up
#: here and hands it to the engine and to the renderer both.
GALLERY: Mapping[str, Fixture] = {
    fixture.name: fixture
    for fixture in (
        CATCHMENT,
        SLIVER_FAN,
        CORNER_HOLE,
        HOLE_IN_HOLE,
        BREAKLINE_CHAIN,
        RIVER,
        NOT_NODED,
        DEGENERATE,
    )
}
