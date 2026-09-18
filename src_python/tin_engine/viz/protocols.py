"""Structural descriptions of the two core types the renderer consumes.

These mirror ``_core.IndexedMesh2``, ``_core.Pslg`` and ``_core.Chain``, and
mirroring is all they do -- this module imports neither ``_core`` nor anything
that would drag the extension in. A member added here without a matching
accessor on the bound type is caught by the join in
``tests/python/test_viz_protocols.py``, which is the only place the two halves
meet before ``cli.py`` does.

Every member is a read-only property, so a mutating accessor cannot satisfy a
protocol by accident and the renderer has no way to write through to a core
buffer.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol

import numpy as np
import numpy.typing as npt


class ChainLike(Protocol):
    """One constraint chain: a run of a PSLG's flat index buffer."""

    @property
    def begin(self) -> int: ...
    @property
    def count(self) -> int: ...
    @property
    def role(self) -> object:
        """The chain's role.

        Typed as ``object`` because naming ``_core.ChainRole`` is exactly the
        import this module exists to avoid. Consumers compare roles for
        equality or use them as mapping keys; neither needs the concrete type.
        """

    @property
    def properties(self) -> int:
        """The chain's feature set, as a bare mask.

        An ``int`` because ``viz/`` may not hold a vocabulary: which bit means
        which feature is ``tin_engine.features``' answer and the composition
        root's to supply, exactly as ``role``'s concrete type is.
        """


class PslgLike(Protocol):
    """The input side of the picture: what was asked for."""

    @property
    def vertices(self) -> npt.NDArray[np.float64]: ...
    @property
    def chains(self) -> Sequence[ChainLike]: ...
    @property
    def chain_indices(self) -> npt.NDArray[np.uint32]: ...
    def indices_of(self, c: int) -> npt.NDArray[np.uint32]: ...


class MeshLike(Protocol):
    """The output side: what the triangulator produced."""

    @property
    def vertices(self) -> npt.NDArray[np.float64]: ...
    @property
    def triangles(self) -> npt.NDArray[np.uint32]: ...
    @property
    def constrained_edges(self) -> npt.NDArray[np.uint8]: ...
    @property
    def triangle_count(self) -> int: ...
    @property
    def empty(self) -> bool: ...
