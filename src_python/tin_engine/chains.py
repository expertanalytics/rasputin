"""The start mesh's chains: the domain's rings, then the features' lines (16b R2, R7).

Increment 16's R6 table as code. The domain's exterior is ``"outer"`` and each
hole ``"hole"``, mask 0 (16 U5; this half moved here from ``cli._domain_chains``);
then every feature line, in feature order, is a ``"breakline"`` with its
feature's mask. A closed line (an unclipped ring) repeats its first index
(increment 8, ruling 2). Shared edges go in twice: the noder merges them and
unions their bits (R7).

Roles leave as strings; ``cli.py`` maps them onto ``ChainRole``. No ``_core``.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import shapely

from tin_engine.domain import DomainPolygon
from tin_engine.features import EdgeVocabulary, TerrainFeature

Chain = tuple[list[int], str, int]


@dataclass(frozen=True, slots=True)
class StartChains:
    """``vertices`` ``(N, 2)`` in the DEM's CRS, and ``(indices, role, mask)``."""

    vertices: npt.NDArray[np.float64]
    chains: list[Chain]


def start_chains(
    domain: DomainPolygon, features: Iterable[TerrainFeature], vocabulary: EdgeVocabulary
) -> StartChains:
    """The domain's chains, then each feature's; every mask is checked against
    ``vocabulary``, which raises on a bit it does not name."""
    blocks: list[npt.NDArray[np.float64]] = []
    chains: list[Chain] = []
    count = 0

    def add(xy: npt.NDArray[np.float64], role: str, mask: int, closed: bool) -> None:
        nonlocal count
        indices = list(range(count, count + len(xy)))
        chains.append((indices + indices[:1] if closed else indices, role, mask))
        blocks.append(xy)
        count += len(xy)

    for k, ring in enumerate((domain.polygon.exterior, *domain.polygon.interiors)):
        add(shapely.get_coordinates(ring)[:-1], "outer" if k == 0 else "hole", 0, False)
    for feature in features:
        vocabulary.names(feature.mask)
        for line in feature.lines:
            xy = shapely.get_coordinates(line)
            closed = len(xy) > 2 and bool((xy[0] == xy[-1]).all())
            add(xy[:-1] if closed else xy, "breakline", feature.mask, closed)
    vertices = np.concatenate(blocks) if blocks else np.zeros((0, 2))
    return StartChains(vertices=vertices.astype(np.float64), chains=chains)
