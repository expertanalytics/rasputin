"""Type stubs for the _core C++ extension.

Hand-written rather than generated: the dimensional asymmetry of ``cross`` --
scalar for 2D, vector for 3D -- is expressed with ``@overload``, and no stub
generator infers that from the pybind11 signatures.
"""

from collections.abc import Iterable, Sequence
from enum import Enum
from typing import final, overload

import numpy as np
import numpy.typing as npt

@final
class Point2:
    """A 2D point with double-precision coordinates. Immutable and hashable."""

    def __init__(self, x: float = ..., y: float = ...) -> None: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __add__(self, other: Point2) -> Point2: ...
    def __sub__(self, other: Point2) -> Point2: ...
    def __mul__(self, scalar: float) -> Point2: ...
    def __rmul__(self, scalar: float) -> Point2: ...
    def __truediv__(self, scalar: float) -> Point2: ...
    def __neg__(self) -> Point2: ...
    def __hash__(self) -> int: ...
    def __repr__(self) -> str: ...

@final
class Point3:
    """A 3D point with double-precision coordinates. Immutable and hashable."""

    def __init__(self, x: float = ..., y: float = ..., z: float = ...) -> None: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def z(self) -> float: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __add__(self, other: Point3) -> Point3: ...
    def __sub__(self, other: Point3) -> Point3: ...
    def __mul__(self, scalar: float) -> Point3: ...
    def __rmul__(self, scalar: float) -> Point3: ...
    def __truediv__(self, scalar: float) -> Point3: ...
    def __neg__(self) -> Point3: ...
    def __hash__(self) -> int: ...
    def __repr__(self) -> str: ...

@overload
def dot(a: Point2, b: Point2) -> float: ...
@overload
def dot(a: Point3, b: Point3) -> float: ...
@overload
def cross(a: Point2, b: Point2) -> float:
    """2D cross product: the scalar z-component of the 3D cross product."""

@overload
def cross(a: Point3, b: Point3) -> Point3:
    """3D cross product: a vector orthogonal to both operands."""

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
    """A flat indexed triangle mesh. Not constructible from Python.

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

# Two ``@overload`` stubs and no union parameter, the pattern ``cross`` above
# uses: pybind11 resolves this by exact type with conversions disabled, and both
# enumerations' ``Ok`` is integer 0, so a single union stub would type-check a
# call the runtime rejects. No implementation stub follows, because a stub file
# may not carry one: mypy rejects "an implementation for an overloaded function"
# in a ``.pyi``, and ``cross`` above declares two overloads and nothing else.
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
    """Bilinear ``(z, valid)`` at ``(N, 2)`` points. ``z`` is 0.0 where ``valid``
    is False, never NaN. Releases the GIL."""

class RefineStatus(Enum):
    """Why :func:`refine` refused, or ``Ok``."""

    Ok = 0
    OffLattice = 1
    NotCounterClockwise = 2
    InvalidTolerance = 3

@final
class RefineOutcome:
    """A status, a message, the refined mesh and four numbers. The arrays are
    read-only views that keep the outcome alive, and empty unless ``ok()``."""

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

def refine(
    view: RasterView,
    mesh: IndexedMesh2,
    edges: npt.ArrayLike,
    masks: npt.ArrayLike,
    *,
    tolerance: float,
    threads: int = ...,
) -> RefineOutcome:
    """Refine a start mesh whose vertices are DEM nodes until every triangle is
    within ``tolerance`` of the DEM. Releases the GIL; the output does not
    depend on ``threads``."""
