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

@final
class Chain:
    """One constraint chain: a run of ``chain_indices``, its role, its flag."""

    @property
    def begin(self) -> int: ...
    @property
    def count(self) -> int:
        """Number of DISTINCT vertices: a ring does not store its closure."""

    @property
    def role(self) -> ChainRole: ...
    @property
    def is_river(self) -> bool: ...

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

def describe(status: CdtStatus) -> str:
    """One sentence of prose for a ``CdtStatus``."""

def build_pslg(
    vertices: npt.ArrayLike,
    chains: Iterable[tuple[Sequence[int], ChainRole, bool]],
) -> PslgBuildResult:
    """Validate a constraint set. The coordinates are copied, so the result
    neither aliases nor keeps alive the array handed in. Invalid input is data,
    not an exception: a ``ValueError`` means the vertex array is not ``(N, 2)``,
    and a ``TypeError`` means a path or filename was passed where coordinates
    belong.
    """

def triangulate(pslg: Pslg, delaunay: bool = ...) -> CdtOutcome:
    """Triangulate an already-noded PSLG, releasing the GIL for the duration."""
