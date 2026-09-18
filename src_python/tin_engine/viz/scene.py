"""The geometry mapping: ``(Pslg, IndexedMesh2)`` -> a drawable ``Scene``.

This module is increment 6b-i in full. It is pure computation over immutable
data, it imports the protocols and nothing else first-party, and it never sees a
file, a path or a CRS -- so it is testable against hand-built dataclasses with
no compiled extension in the process.

Two independent derivations of "this edge is a constraint" meet here: the mesh's
``constrained_edges`` bitmask and the input PSLG's chains. They can disagree,
and the scene **records the disagreement rather than reconciling it** -- see
``FindingKind``. That is only meaningful because the two are genuinely different
objects, which is why ``SceneEdge.constrained`` (the mask's verdict) and
``SceneEdge.role`` (the join's) stay separate fields.

Bit ``e`` of a triangle's mask is the edge ``(v[e], v[(e + 1) % 3])``, per
``indexed_mesh.hpp``. This is *not* CGAL's "edge ``e`` is opposite vertex ``e``"
convention, which is a rotation of it; a scene that rotates it puts every
constraint stroke on the wrong edge and still draws a plausible picture.
``_mesh_edges`` is the only place the convention is spelled out.
"""

from __future__ import annotations

import enum
import itertools
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from .protocols import MeshLike, PslgLike

#: Half-extent applied to a degenerate bbox axis, in world units (metres).
#: A zero-extent axis is a division by zero in the viewport transform, and a
#: picture that is an accidental point.
PAD = 0.5

Pair = tuple[int, int]


class SceneKind(enum.Enum):
    """What the picture is of, which decides how 6b-ii presents it."""

    MESH = "mesh"
    OK_BUT_EMPTY = "ok-but-empty"
    """A *successful* backend call that produced no triangles (04-cdt risk 5).
    Distinct from ``FAILED`` because it is the silent mode, and the viewer is
    its second guard."""
    FAILED = "failed"


class FindingKind(enum.Enum):
    """A disagreement between the mesh's mask and the input's chains."""

    MASKED_EDGE_WITHOUT_CHAIN = "masked-edge-without-chain"
    CHAIN_EDGE_WITHOUT_MASK = "chain-edge-without-mask"


@dataclass(frozen=True, slots=True)
class Finding:
    """One disagreement, on one undirected edge ``(a, b)`` with ``a < b``."""

    kind: FindingKind
    a: int
    b: int


@dataclass(frozen=True, slots=True)
class SceneEdge:
    """One undirected edge to draw, endpoints canonically ordered ``a < b``."""

    a: int
    b: int
    constrained: bool
    """The mesh mask's verdict. ``False`` whenever there is no mesh to ask."""
    role: object | None
    """The role of the first input chain containing this pair, or ``None``.
    Opaque: the scene compares roles and uses them as keys, never reads them."""
    is_river: bool
    """``True`` if *any* contributing chain was a river."""


@dataclass(frozen=True, slots=True)
class BBox:
    """World-space bounds of everything drawn, padded away from zero extent."""

    min_x: float
    min_y: float
    max_x: float
    max_y: float
    padded: bool


@dataclass(frozen=True, slots=True)
class Scene:
    """Everything 6b-ii needs, and nothing it does not."""

    kind: SceneKind
    vertices: npt.NDArray[np.float64]
    triangles: npt.NDArray[np.uint32]
    edges: tuple[SceneEdge, ...]
    findings: tuple[Finding, ...]
    bbox: BBox


def _key(a: int, b: int) -> Pair:
    """The canonical undirected form of an edge."""
    return (a, b) if a < b else (b, a)


def _mesh_edges(mesh: MeshLike) -> tuple[set[Pair], set[Pair]]:
    """Every undirected mesh edge, and the subset the mask calls constrained.

    Each interior edge is shared by two triangles; deduplicating here is what
    keeps it from being drawn at double weight, and what keeps a mask bit set on
    only one of the two incident triangles from being lost.
    """
    triangles = np.asarray(mesh.triangles, dtype=np.uint32)
    masks = np.asarray(mesh.constrained_edges, dtype=np.uint8)
    edges: set[Pair] = set()
    masked: set[Pair] = set()
    for triangle, mask in zip(triangles, masks, strict=True):
        corner = [int(v) for v in triangle]
        for e in range(3):
            pair = _key(corner[e], corner[(e + 1) % 3])
            edges.add(pair)
            if int(mask) >> e & 1:
                masked.add(pair)
    return edges, masked


def _chain_edges(pslg: PslgLike, closed_roles: Sequence[object]) -> dict[Pair, tuple[object, bool]]:
    """Join each undirected pair to ``(role, is_river)`` over the input chains.

    The earliest chain containing a pair wins the role -- roles are opaque and
    so cannot be ranked by value -- while ``is_river`` is the OR over every
    contributing chain, because the bit means "some chain here was a river".

    A ``Pslg`` never stores a ring's closing edge, and ``ChainLike.role`` is
    typed ``object`` so this module cannot name ``ChainRole``. ``closed_roles``
    is therefore supplied by ``cli.py``, the composition root that knows the
    enum; without it every ring's closure is a false finding. Membership uses
    only ``==``/``__hash__``.

    Fewer than three vertices is not a ring, so ``len(walk) > 2`` appends no
    closure. The invariant-critical half of that guard is the **one**-vertex
    chain: closing it emits ``_key(a, a)``, a self-loop, which violates
    ``SceneEdge``'s documented ``a < b`` ordering --
    ``test_a_single_vertex_ring_is_never_closed_into_a_self_loop`` is what
    defends it. Double emission is *not* the reason, whatever the design once
    said: ``joined`` is keyed on ``_key``, so closing a two-vertex chain would
    write the same key twice and the second write would only OR ``is_river``.
    """
    joined: dict[Pair, tuple[object, bool]] = {}
    for c, chain in enumerate(pslg.chains):
        walk = [int(i) for i in pslg.indices_of(c)]
        if len(walk) > 2 and any(chain.role == closed for closed in closed_roles):
            walk.append(walk[0])
        for a, b in itertools.pairwise(walk):
            pair = _key(a, b)
            previous = joined.get(pair)
            if previous is None:
                joined[pair] = (chain.role, chain.is_river)
            else:
                joined[pair] = (previous[0], previous[1] or chain.is_river)
    return joined


def _pad(lo: float, hi: float) -> tuple[float, float, bool]:
    """Widen a degenerate axis about its own centre, reporting whether it did."""
    if hi > lo:
        return lo, hi, False
    centre = (lo + hi) / 2.0
    return centre - PAD, centre + PAD, True


def _bbox(vertices: npt.NDArray[np.float64]) -> BBox:
    """Bounds of every vertex, or a unit box at the origin if there are none."""
    if vertices.size == 0:
        return BBox(-PAD, -PAD, PAD, PAD, True)
    lo = vertices.min(axis=0)
    hi = vertices.max(axis=0)
    min_x, max_x, padded_x = _pad(float(lo[0]), float(hi[0]))
    min_y, max_y, padded_y = _pad(float(lo[1]), float(hi[1]))
    return BBox(min_x, min_y, max_x, max_y, padded_x or padded_y)


def _findings(edges: Iterable[SceneEdge]) -> tuple[Finding, ...]:
    """The disagreements, in the edges' own (sorted) order."""
    found: list[Finding] = []
    for edge in edges:
        if edge.constrained and edge.role is None:
            found.append(Finding(FindingKind.MASKED_EDGE_WITHOUT_CHAIN, edge.a, edge.b))
        elif edge.role is not None and not edge.constrained:
            found.append(Finding(FindingKind.CHAIN_EDGE_WITHOUT_MASK, edge.a, edge.b))
    return tuple(found)


def build_scene(
    pslg: PslgLike,
    mesh: MeshLike | None = None,
    *,
    ok: bool = True,
    closed_roles: Sequence[object] = (),
) -> Scene:
    """Map a triangulation attempt onto the primitives a renderer draws.

    ``ok`` is the backend status reduced to the one bit the scene uses; the
    status name and message are pure passthrough to the header band and belong
    to ``svg.py``. A non-``ok`` run draws the input PSLG alone -- never a blank
    page -- and so does an ``ok`` run that produced no triangles.

    Raises:
        ValueError: if any drawn coordinate is not finite. The PSLG validator's
            stage 3 rejects these, so reaching here means something upstream is
            wrong; without the check a NaN poisons the bbox and the SVG renders
            as nothing, in silence.
    """
    drawable = mesh if ok and mesh is not None and not mesh.empty else None
    source = pslg.vertices if drawable is None else drawable.vertices
    vertices = np.asarray(source, dtype=np.float64)
    if not bool(np.isfinite(vertices).all()):
        raise ValueError("every scene coordinate must be finite; got NaN or infinity")

    if not ok:
        kind = SceneKind.FAILED
    elif drawable is None:
        kind = SceneKind.OK_BUT_EMPTY
    else:
        kind = SceneKind.MESH

    joined = _chain_edges(pslg, closed_roles)
    if drawable is None:
        # No mask to disagree with, so findings are suppressed: every chain edge
        # would otherwise be an alarm on a picture whose real alarm is the
        # status. The join itself still runs -- roles are kept, because the
        # PSLG-only picture is drawn in role colours.
        drawn, masked = set(joined), set[Pair]()
        triangles = np.zeros((0, 3), dtype=np.uint32)
    else:
        drawn, masked = _mesh_edges(drawable)
        drawn |= set(joined)
        triangles = np.asarray(drawable.triangles, dtype=np.uint32)

    edges = tuple(
        SceneEdge(a, b, (a, b) in masked, *joined.get((a, b), (None, False)))
        for a, b in sorted(drawn)
    )
    findings = _findings(edges) if drawable is not None else ()
    return Scene(kind, vertices, triangles, edges, findings, _bbox(vertices))
