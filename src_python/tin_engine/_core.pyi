"""Type stubs for the _core C++ extension.

Hand-written rather than generated: the pybind11 enums are declared as
``enum.Enum`` for the checker's benefit (see below), and no stub generator
infers that from the pybind11 signatures.
"""

from collections.abc import Iterable, Sequence
from enum import Enum
from typing import final, overload

import numpy as np
import numpy.typing as npt

#: The standard library's bounds checks compiled into this build: "libc++ fast",
#: "libc++ extensive", "libc++ debug", "libstdc++ assertions", or "none".
hardening: str

# The CDT surface (increment 6a). The three enums are spelled as ``enum.Enum``
# for the checker's benefit; at runtime they are pybind11 enum objects, which
# support the same ``name``, ``value``, ``__members__`` and call-by-value
# protocol but are not instances of ``enum.Enum``.
#
# IDENTITY IS THE PART THAT DOES NOT CARRY OVER. Compare with ``==``, never
# with ``is``. A real ``Enum`` member is a singleton, so ``is`` is the
# idiomatic comparison and no checker will flag it here -- but every read of an
# enum-valued attribute builds a fresh pybind11 object, so
# ``outcome.status is CdtStatus.Ok`` is False even when the status IS Ok, and
# so is ``outcome.status is outcome.status``. Only ``CdtStatus.Ok is
# CdtStatus.Ok``, comparing the class attribute with itself, holds. The trap
# therefore fires on exactly the values a caller tests: the ones that came out
# of C++. ``==``, ``!=``, ``in`` against a tuple, and use as a dict key all
# behave as expected, because they go through ``__eq__`` and ``__hash__``.

class ChainRole(Enum):
    """The role a constraint chain plays in the domain."""

    Outer = 0
    Hole = 1
    Breakline = 2

class PslgError(Enum):
    """Why the validator rejected a proposed PSLG."""

    NoOuterChain = 0
    ChainTooShort = 1
    IndexOutOfRange = 2
    NonFiniteVertex = 3
    StoredClosure = 4
    WrongWinding = 5
    DegenerateRing = 6
    VertexCountOverflow = 7

class CdtStatus(Enum):
    """What the triangulation backend did, grouped by what to do about it."""

    Ok = 0
    NotRun = 1
    NotNoded = 2
    DegenerateGeometry = 3
    InvalidTopology = 4
    MalformedInput = 5
    BackendFailure = 6

class NodeStatus(Enum):
    """What the noder did, grouped by what to do about it.

    ``Ok`` is 0 in this enumeration and in :class:`CdtStatus`, and the two are
    different types: compare a status against its own enumeration's member.
    """

    Ok = 0
    NotRun = 1
    InvalidSnapSpacing = 2
    CoordinateOutOfRange = 3
    RingCollapsed = 4
    RingDegenerateAfterSnap = 5
    NonSimpleRing = 6
    NotConverged = 7
    MalformedOutput = 8

@final
class Chain:
    """One constraint chain: a run of ``chain_indices``, its role, its set."""

    @property
    def begin(self) -> int: ...
    @property
    def count(self) -> int:
        """Number of DISTINCT vertices: a ring does not store its closure."""

    @property
    def role(self) -> ChainRole: ...
    @property
    def properties(self) -> int:
        """The feature property set, as a bare 32-bit mask.

        An ``int``: ``EdgeProperties`` is deliberately not bound, and the
        meaning of each bit lives in ``tin_engine.features.EdgeVocabulary``.
        """

@final
class PslgDiagnostic:
    """One reason a proposed PSLG was rejected."""

    @property
    def error(self) -> PslgError: ...
    @property
    def chain(self) -> int:
        """Index into ``chains``, or the ``kNoChain`` sentinel 2**32 - 1."""

    @property
    def vertex(self) -> int:
        """A VALID index into ``vertices``, or the ``kNoVertex`` sentinel
        2**32 - 1 -- never the offending out-of-range value. Indexing with it
        therefore cannot perform the out-of-bounds read that ``IndexOutOfRange``
        exists to prevent; the offending value is in ``message``."""

    @property
    def message(self) -> str: ...

@final
class PslgBuildResult:
    """A validated ``Pslg``, or the whole list of reasons why not."""

    @property
    def ok(self) -> bool: ...
    @property
    def pslg(self) -> Pslg | None: ...
    @property
    def diagnostics(self) -> list[PslgDiagnostic]: ...

@final
class Pslg:
    """A validated planar straight-line graph. Not constructible from Python."""

    @property
    def vertices(self) -> npt.NDArray[np.float64]:
        """Read-only ``(N, 2)`` view of the vertex buffer."""

    @property
    def chains(self) -> list[Chain]:
        """The chains, in order, as frozen ``Chain`` records.

        Unlike :attr:`vertices`, :attr:`chain_indices` and :attr:`indices_of`,
        which are zero-copy views, this **copies**: a fresh list is built on
        every read. Bind it once rather than re-reading it inside a loop.
        """
    @property
    def chain_indices(self) -> npt.NDArray[np.uint32]:
        """Read-only ``(M,)`` view of the flat index buffer."""

    def indices_of(self, c: int) -> npt.NDArray[np.uint32]:
        """Read-only view of chain ``c``'s indices. ``IndexError`` if out of range."""

@final
class NodedPslg:
    """A PSLG after snap rounding. Not constructible from Python.

    NOT a subclass of :class:`Pslg` and there is no conversion either way:
    :func:`triangulate` takes one of these, so un-noded input is
    unrepresentable at the entry point rather than diagnosed inside it. It has
    the same four shared accessors, because the renderer's scene join is
    structural.

    Node ids are **not** input vertex ids -- the node set is ordered by its
    lattice point -- so follow an input vertex across with
    :attr:`node_of_input_vertex`.
    """

    @property
    def vertices(self) -> npt.NDArray[np.float64]:
        """Read-only ``(N, 2)`` view of the node coordinates."""

    @property
    def chains(self) -> list[Chain]:
        """The chains, in order. Copies on every read, like ``Pslg.chains``."""

    @property
    def chain_indices(self) -> npt.NDArray[np.uint32]:
        """Read-only ``(M,)`` view of the flat index buffer."""

    def indices_of(self, c: int) -> npt.NDArray[np.uint32]:
        """Read-only view of chain ``c``'s indices. ``IndexError`` if out of range."""

    @property
    def grid_spacing(self) -> float:
        """The spacing this graph was noded at. The ``SnapGrid`` does not cross."""

    @property
    def edge_properties(self) -> npt.NDArray[np.uint32]:
        """Read-only ``(E,)`` view of the dense per-edge property masks.

        Index-aligned with the flat edge enumeration: chain ``c``'s edge ``k``
        sits at ``sum(edge_count(j) for j < c) + k``. Each entry is the union
        over every input chain that contributed geometry to that edge, and 0
        means *unclassified* rather than *wrong*.
        """

    @property
    def node_of_input_vertex(self) -> npt.NDArray[np.uint32]:
        """Read-only ``(N,)`` view mapping each input vertex onto its node id.

        Total over the input's vertex array, unreferenced vertices included.
        """

@final
class NodeOutcome:
    """A ``NodedPslg``, or a status and a message saying why not.

    ``ok()`` is a method, not a property -- ``if outcome.ok`` is truthy for a
    bound method and never sees a failure.
    """

    @property
    def status(self) -> NodeStatus: ...
    @property
    def message(self) -> str: ...
    @property
    def pslg(self) -> NodedPslg | None:
        """The noded graph, or None. Engaged iff ``status`` is ``Ok``."""

    def ok(self) -> bool: ...

@final
class IndexedMesh2:
    """A flat indexed triangle mesh, from ``triangulate`` or ``indexed_mesh``.

    Bit ``e`` of a mask is set iff the edge ``(v[e], v[(e + 1) % 3])`` is
    constrained -- not CGAL's "edge e is opposite vertex e", which is a
    rotation of this convention. The three arrays are read-only zero-copy views
    that keep the mesh alive.
    """

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

def indexed_mesh(
    vertices: npt.ArrayLike, triangles: npt.ArrayLike, constrained_edges: npt.ArrayLike
) -> IndexedMesh2:
    """An ``IndexedMesh2`` from arrays, copied. A shape, an index out of range,
    a mask above 7 or a non-finite coordinate is a ``ValueError``; orientation
    is ``refine``'s check, not this one's."""

@final
class CdtOutcome:
    """A status, a message and a mesh. ``ok()`` is a method, not a property."""

    @property
    def status(self) -> CdtStatus: ...
    @property
    def message(self) -> str: ...
    @property
    def mesh(self) -> IndexedMesh2: ...
    def ok(self) -> bool: ...

# Two ``@overload`` stubs and no union parameter, because pybind11 resolves
# this by exact type with conversions disabled, and both
# enumerations' ``Ok`` is integer 0, so a single union stub would type-check a
# call the runtime rejects. No implementation stub follows, because a stub file
# may not carry one: mypy rejects "an implementation for an overloaded function"
# in a ``.pyi``.
@overload
def describe(status: CdtStatus) -> str:
    """One sentence of prose for a ``CdtStatus``."""

@overload
def describe(status: NodeStatus) -> str:
    """One sentence of prose for a ``NodeStatus``."""

def build_pslg(
    vertices: npt.ArrayLike,
    chains: Iterable[tuple[Sequence[int], ChainRole, int]],
) -> PslgBuildResult:
    """Validate a constraint set. The coordinates are copied, so the result
    neither aliases nor keeps alive the array handed in. Invalid input is data,
    not an exception: a ``ValueError`` means the vertex array is not ``(N, 2)``
    or a chain's property mask is negative or carries a bit at or above 32, and
    a ``TypeError`` means a path or filename was passed where coordinates
    belong.
    """

def node(pslg: Pslg, spacing: float, max_rounds: int = ...) -> NodeOutcome:
    """Snap-round a validated PSLG, releasing the GIL for the duration.

    ``spacing`` has no default: the right value is a policy question about the
    data, and this layer is not the composition root. A refused noding is a
    status, not an exception; only a mis-shaped argument raises.
    """

def triangulate(pslg: NodedPslg, delaunay: bool = ...) -> CdtOutcome:
    """Triangulate a noded PSLG, releasing the GIL for the duration.

    Takes a ``NodedPslg`` and not a ``Pslg``: call :func:`node` first.
    """

@final
class RasterView:
    """A zero-copy view over a 2-D DEM array. Built by :func:`raster_view`;
    holds the array alive."""

def raster_view(
    array: npt.NDArray[np.float32] | npt.NDArray[np.float64],
    *,
    x_min: float,
    y_max: float,
    delta_x: float,
    delta_y: float,
    nodata: float | None = ...,
) -> RasterView:
    """View a C-contiguous float32 or float64 array without copying it. Any other
    dtype or layout is a ``TypeError``; rows and columns come from its shape."""

def sample(
    view: RasterView, points: npt.ArrayLike
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.bool_]]:
    """Bilinear ``(z, valid)`` at ``(N, 2)`` points. A point that is a DEM node
    bit for bit reads that node alone; any other point is invalid if one of its
    four corners is NoData or NaN. ``z`` is 0.0 where ``valid`` is False, never
    NaN. Releases the GIL."""

class RefineStatus(Enum):
    """Why :func:`refine` refused, or ``Ok``."""

    Ok = 0
    OutsideGrid = 1
    NotCounterClockwise = 2
    InvalidTolerance = 3

class RefineOutcome:
    """A status, a message, the refined mesh and four numbers. The arrays are
    read-only views that keep the outcome alive, and empty unless ``ok()``.
    Not final: :class:`PointRefineOutcome` extends it."""

    @property
    def status(self) -> RefineStatus: ...
    @property
    def message(self) -> str: ...
    def ok(self) -> bool: ...
    @property
    def vertices(self) -> npt.NDArray[np.float64]: ...
    @property
    def z(self) -> npt.NDArray[np.float64]: ...
    @property
    def valid(self) -> npt.NDArray[np.bool_]: ...
    @property
    def triangles(self) -> npt.NDArray[np.uint32]: ...
    @property
    def edges(self) -> npt.NDArray[np.uint32]: ...
    @property
    def masks(self) -> npt.NDArray[np.uint32]: ...
    @property
    def rounds(self) -> int: ...
    @property
    def inserted(self) -> int: ...
    @property
    def flips(self) -> int: ...
    @property
    def max_error(self) -> float: ...
    @property
    def uncovered(self) -> int: ...
    @property
    def carved(self) -> int: ...
    @property
    def legalise_seconds(self) -> float: ...
    @property
    def scan_seconds(self) -> float: ...
    @property
    def split_seconds(self) -> float: ...
    @property
    def quality_inserted(self) -> int: ...
    @property
    def quality_skipped(self) -> int: ...
    @property
    def quality_seconds(self) -> float: ...
    @property
    def feet(self) -> int: ...
    @property
    def feet_refused(self) -> int: ...

def refine(
    view: RasterView,
    mesh: IndexedMesh2,
    edges: npt.ArrayLike,
    masks: npt.ArrayLike,
    *,
    tolerance: float,
    threads: int = ...,
    min_angle_deg: float = ...,
    constraint_feet: bool = ...,
    frozen_mask: int = ...,
) -> RefineOutcome:
    """Refine a start mesh whose vertices lie in the DEM's node rectangle (off-node
    ones keep their position and get bilinear z) until every triangle is
    within ``tolerance`` of the DEM. Releases the GIL; the output does not
    depend on ``threads``. ``min_angle_deg`` > 0 first improves the start
    mesh's angles with DEM nodes; 0 is off. ``constraint_feet`` inserts a
    node's foot on a nearby constraint segment instead of the node. No vertex
    goes on an edge whose mask meets ``frozen_mask`` (unsigned; 0 is off)."""

@final
class CheckPoints:
    """The final check's store: points in the grid's frame with float32 z,
    filed by lattice cell. ``add`` blocks, ``freeze`` once, then pass it to
    :func:`refine_points`. ``add`` after ``freeze`` is a ``RuntimeError``."""

    def __init__(
        self, *, x_min: float, y_max: float, spacing: float, rows: int, cols: int
    ) -> None: ...
    def add(self, xy: npt.NDArray[np.float64], z: npt.NDArray[np.float32]) -> None:
        """File ``(N, 2)`` float64 points with ``(N,)`` float32 z; any other
        shape or dtype is a ``ValueError``."""
    def freeze(self) -> None: ...
    @property
    def size(self) -> int: ...
    @property
    def duplicates(self) -> int: ...
    @property
    def outside(self) -> int: ...

@final
class PointRefineOutcome(RefineOutcome):
    """What :func:`refine_points` returned: :class:`RefineOutcome`'s fields plus
    the check points that coincide with a start vertex."""

    @property
    def vertices(self) -> npt.NDArray[np.float64]:
        """The start vertices as given, then the inserted check points."""
    @property
    def max_error(self) -> float:
        """Largest ``|z - plane|`` over the check points, in valid triangles."""
    @property
    def coincident(self) -> int: ...
    @property
    def coincident_max_error(self) -> float: ...
    @property
    def on_frozen(self) -> int:
        """Check points on a frozen edge, never inserted, each counted once."""
    @property
    def on_frozen_max_error(self) -> float: ...
    @property
    def strip_points(self) -> int: ...
    @property
    def strip_inserted(self) -> int: ...
    @property
    def strip_max_error(self) -> float: ...
    @property
    def strip_refused(self) -> int: ...
    @property
    def strip_refused_max_error(self) -> float: ...
    @property
    def nodes_inserted(self) -> int:
        """``refine_strip`` only: DEM nodes its rescan inserted."""

@final
class ConstraintCheckPoints:
    """The edge strip's check points, filed by constraint edge; built only by
    :func:`constraint_check_points`, read-only."""

    @property
    def size(self) -> int: ...
    @property
    def no_data(self) -> int: ...
    @property
    def duplicates(self) -> int: ...
    @property
    def edge_count(self) -> int: ...

def constraint_check_points(
    view: RasterView, vertices: npt.ArrayLike, edges: npt.ArrayLike
) -> ConstraintCheckPoints:
    """Grid-line crossings of each constraint edge and the midpoints between
    neighbours, z from ``view``. A refused input is a ``ValueError``.
    Releases the GIL."""

def refine_strip(
    view: RasterView,
    strip: ConstraintCheckPoints,
    vertices: npt.ArrayLike,
    triangles: npt.ArrayLike,
    z: npt.ArrayLike,
    valid: npt.ArrayLike,
    edges: npt.ArrayLike,
    masks: npt.ArrayLike,
    *,
    tolerance: float,
    threads: int = ...,
    frozen_mask: int = ...,
) -> PointRefineOutcome:
    """The edge strip on the projected path: refine's output refined until every
    strip point is within ``tolerance``, the DEM's nodes rescanned in every
    triangle it writes. Nothing goes on an edge whose mask meets
    ``frozen_mask``; a strip point on one is a ``RuntimeError``. Releases the
    GIL."""

def refine_points(
    points: CheckPoints,
    vertices: npt.ArrayLike,
    triangles: npt.ArrayLike,
    z: npt.ArrayLike,
    valid: npt.ArrayLike,
    edges: npt.ArrayLike,
    masks: npt.ArrayLike,
    *,
    tolerance: float,
    threads: int = ...,
    strip: ConstraintCheckPoints | None = ...,
    frozen_mask: int = ...,
) -> PointRefineOutcome:
    """Refine phase 1's mesh until every check point in a frozen store, and
    every point of ``strip``, is within ``tolerance``. Releases the GIL; the
    output does not depend on ``threads``. A check point on an edge whose mask
    meets ``frozen_mask`` is counted in ``on_frozen``, not inserted."""

@final
class SeamOutcome:
    """What :func:`refine_seam` returned. ``a < b`` by ``(x, y)``; the arrays are
    read-only views, from ``a`` to ``b``, that keep the outcome alive."""

    @property
    def a(self) -> tuple[float, float]: ...
    @property
    def b(self) -> tuple[float, float]: ...
    @property
    def z_a(self) -> float | None: ...
    @property
    def z_b(self) -> float | None: ...
    @property
    def points(self) -> npt.NDArray[np.float64]:
        """``(K, 2)`` world points inserted on the seam."""
    @property
    def z(self) -> npt.NDArray[np.float64]: ...
    @property
    def s(self) -> npt.NDArray[np.float64]:
        """``(K,)`` parameters from ``a``, strictly increasing in ``(0, 1)``."""
    @property
    def check_points(self) -> int: ...
    @property
    def no_data(self) -> int: ...
    @property
    def max_error(self) -> float: ...

def refine_seam(
    view: RasterView, a: tuple[float, float], b: tuple[float, float], *, tolerance: float
) -> SeamOutcome:
    """The seam pass for the seam edge ``(a, b)``: check points inserted by a
    one-dimensional greedy until each is within ``tolerance``. The same output
    for either order of the ends. A bad tolerance, ``a == b`` or an end outside
    the node rectangle is a ``ValueError``. Releases the GIL."""

@final
class UpstreamOutcome:
    """What :func:`upstream` returned. Bounds are inclusive and meaningful
    when ``nodes_in`` > 0."""

    @property
    def mask(self) -> npt.NDArray[np.uint8]:
        """Read-only ``(rows, cols)``: 1 for a node in the catchment, 0 otherwise."""
    @property
    def nodes_in(self) -> int: ...
    @property
    def row_min(self) -> int: ...
    @property
    def row_max(self) -> int: ...
    @property
    def col_min(self) -> int: ...
    @property
    def col_max(self) -> int: ...
    @property
    def touches_edge(self) -> bool: ...
    @property
    def touches_nodata(self) -> bool: ...

def upstream(view: RasterView, seed: npt.ArrayLike) -> UpstreamOutcome:
    """Every node draining into a seed of the ``(rows, cols)`` mask
    (Priority-Flood); another shape is a ``ValueError``. Releases the GIL."""

@final
class AccumulateOutcome:
    """What :func:`accumulate` returned; every array is read-only ``(rows, cols)``."""

    @property
    def count(self) -> npt.NDArray[np.uint32]:
        """0 on NoData, else the nodes draining through the node, itself included."""
    @property
    def reach(self) -> npt.NDArray[np.uint8]:
        """Bit 0: the catchment touches the window's edge; bit 1: it touches NoData."""
    @property
    def flow_to(self) -> npt.NDArray[np.uint8]:
        """``3*(dr+1)+(dc+1)`` of the neighbour drained to; 255 outlet or NoData."""

def accumulate(view: RasterView) -> AccumulateOutcome:
    """Flow accumulation from :func:`upstream`'s flood; 2^32 nodes or more is
    a ``ValueError``. Releases the GIL."""

class ReduceStatus(Enum):
    """Why :func:`reduce_ring` refused, or ``Ok``."""

    Ok = 0
    InvalidTolerance = 1
    NotCounterClockwise = 2
    TooFewVertices = 3

@final
class ReduceOutcome:
    """What :func:`reduce_ring` returned: the ring, open, a status and counts."""

    @property
    def ring(self) -> npt.NDArray[np.float64]: ...
    @property
    def status(self) -> ReduceStatus: ...
    @property
    def collinear(self) -> int: ...
    @property
    def collapses(self) -> int: ...
    @property
    def rejected_crossing(self) -> int: ...
    @property
    def rejected_seed(self) -> int: ...
    @property
    def rejected_tolerance(self) -> int: ...

def reduce_ring(ring: npt.ArrayLike, tolerance: float, keep: npt.ArrayLike) -> ReduceOutcome:
    """Reduce an open counter-clockwise ``(N, 2)`` ring to ``tolerance``, keeping
    its area and the ``(K, 2)`` keep-points inside. Releases the GIL."""
