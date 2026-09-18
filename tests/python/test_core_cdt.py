"""The CDT binding surface: increment 6a.

Committed RED, before `bindings/core.cpp` exposes any of it -- see
`docs/increments/README.md`, "The loop". Every name reached here is looked up
on `_core` *inside* a test or fixture rather than imported at module scope, so
the failure mode is a named `AttributeError` on the symbol that is missing and
pytest still collects the file. A collection error would hide the next real
problem behind an import.

Scope is 6a only: the bindings, the stubs and `viz/protocols.py`. The renderer,
the fixtures and `rasputin draw` are 6b and appear nowhere here.

What the design (`docs/increments/06-cdt-viewer.md`, "The binding surface")
makes this suite responsible for, in descending order of what a defect would
cost:

* the zero-copy array contract, and above all that an array outlives a dropped
  mesh -- risk 4, the one genuinely dangerous line in the increment;
* that the prohibited surface stays absent -- risk 3, binding-surface gravity;
* that `build_pslg` returns the *whole* diagnostics vector;
* that `triangulate`'s `delaunay` flag is reachable both ways, and that the GIL
  is actually released (see `TestGil` for how that is measured, and for the
  control that proves the measurement can fail).
"""

from __future__ import annotations

import ast
import gc
import sys
import threading
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any, TypeVar

import numpy as np
import pytest

from tin_engine import _core

T = TypeVar("T")

REPO_ROOT = Path(__file__).resolve().parents[2]
STUB = REPO_ROOT / "src_python" / "tin_engine" / "_core.pyi"

# terrain::kNoVertex / kNoChain: uint32 max. Bound as-is per the design, which
# is why a consumer must be told that `vertex` is either this or a valid index.
NO_VERTEX = 2**32 - 1

# Realistic coordinate magnitudes, per the ordering ruling: a UTM 33N-shaped
# domain, not the unit square. Robustness behaviour is magnitude-sensitive and a
# fixture in [0,1]^2 silently exercises the easy case.
EAST = 430_000.0
NORTH = 6_900_000.0


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------


@pytest.fixture
def square() -> np.ndarray:
    """A counterclockwise outer ring at UTM33 magnitudes, closure not stored."""
    return np.array(
        [
            [EAST, NORTH],
            [EAST + 100.0, NORTH],
            [EAST + 100.0, NORTH + 80.0],
            [EAST, NORTH + 80.0],
        ],
        dtype=np.float64,
    )


@pytest.fixture
def square_pslg(square: np.ndarray) -> Any:
    result = _core.build_pslg(square, [([0, 1, 2, 3], _core.ChainRole.Outer, False)])
    assert result.ok, [d.message for d in result.diagnostics]
    return result.pslg


@pytest.fixture
def breakline_pslg(square: np.ndarray) -> Any:
    """Outer ring plus one interior breakline: a mesh with a constrained
    interior edge, which is the only configuration that distinguishes the mask
    convention from its CGAL rotation."""
    vertices = np.vstack([square, [[EAST + 20.0, NORTH + 40.0], [EAST + 80.0, NORTH + 40.0]]])
    result = _core.build_pslg(
        vertices,
        [
            ([0, 1, 2, 3], _core.ChainRole.Outer, False),
            ([4, 5], _core.ChainRole.Breakline, True),
        ],
    )
    assert result.ok, [d.message for d in result.diagnostics]
    return result.pslg


@pytest.fixture
def mesh(breakline_pslg: Any) -> Any:
    outcome = _core.triangulate(breakline_pslg)
    assert outcome.ok(), outcome.message
    return outcome.mesh


@pytest.fixture
def crossing_pslg(square: np.ndarray) -> Any:
    """Two breaklines that cross in the interior: valid as a PSLG (a Pslg
    promises no pairwise disjointness) and not noded, so the backend must fail
    rather than produce a mesh."""
    vertices = np.vstack(
        [
            square,
            [
                [EAST + 20.0, NORTH + 20.0],
                [EAST + 80.0, NORTH + 60.0],
                [EAST + 20.0, NORTH + 60.0],
                [EAST + 80.0, NORTH + 20.0],
            ],
        ]
    )
    result = _core.build_pslg(
        vertices,
        [
            ([0, 1, 2, 3], _core.ChainRole.Outer, False),
            ([4, 5], _core.ChainRole.Breakline, False),
            ([6, 7], _core.ChainRole.Breakline, False),
        ],
    )
    assert result.ok, [d.message for d in result.diagnostics]
    return result.pslg


@pytest.fixture
def broken_build() -> Any:
    """Input wrong in four independent ways at once.

    `pslg_builder.hpp` is explicit that bulk-wrong input reported one failure at
    a time turns one fix into N round trips. One breakage per validator stage,
    so no stage can mask another:

    * chain 0 is a Hole with two vertices          -- stage 1, ChainTooShort
    * chain 1 names vertex 99                      -- stage 2, IndexOutOfRange
    * vertex 4 is NaN                              -- stage 3, NonFiniteVertex
    * no chain has role Outer                      -- stage 6, NoOuterChain
    """
    vertices = np.array(
        [
            [EAST, NORTH],
            [EAST + 100.0, NORTH],
            [EAST + 100.0, NORTH + 80.0],
            [EAST, NORTH + 80.0],
            [float("nan"), NORTH + 10.0],
        ],
        dtype=np.float64,
    )
    chains = [
        ([0, 1], _core.ChainRole.Hole, False),
        ([2, 99], _core.ChainRole.Breakline, False),
    ]
    return _core.build_pslg(vertices, chains)


# --------------------------------------------------------------------------
# Enums and status round-tripping
# --------------------------------------------------------------------------


class TestChainRole:
    def test_has_exactly_the_three_roles(self) -> None:
        assert set(_core.ChainRole.__members__) == {"Outer", "Hole", "Breakline"}

    def test_values_match_the_cpp_enumerator_order(self) -> None:
        # core/pslg.hpp: Outer, Hole, Breakline. A renderer keying a colour off
        # the integer value depends on this, so it is pinned rather than left
        # to whatever pybind11 happened to emit.
        assert _core.ChainRole.Outer.value == 0
        assert _core.ChainRole.Hole.value == 1
        assert _core.ChainRole.Breakline.value == 2

    def test_round_trips_through_its_integer_value(self) -> None:
        for role in _core.ChainRole.__members__.values():
            assert _core.ChainRole(role.value) == role


class TestCdtStatus:
    NAMES = (
        "Ok",
        "NotRun",
        "NotNoded",
        "DegenerateGeometry",
        "InvalidTopology",
        "MalformedInput",
        "BackendFailure",
    )

    def test_has_exactly_the_seven_statuses(self) -> None:
        assert set(_core.CdtStatus.__members__) == set(self.NAMES)

    def test_ok_is_zero_and_every_status_round_trips(self) -> None:
        assert _core.CdtStatus.Ok.value == 0
        for status in _core.CdtStatus.__members__.values():
            assert _core.CdtStatus(status.value) == status

    def test_describe_reaches_python_for_every_status(self) -> None:
        # Including the non-Ok cases: 6b's header band presents exactly this
        # text, so a status whose description never crossed is a blank band.
        for name in self.NAMES:
            text = _core.describe(getattr(_core.CdtStatus, name))
            assert isinstance(text, str)
            assert text.strip(), f"{name} describes to nothing"

    def test_descriptions_are_distinct(self) -> None:
        # Wording is not pinned -- testing.md forbids pinning prose -- but a
        # switch arm that fell through to a neighbour's text is a real defect
        # and distinctness is what catches it.
        texts = [_core.describe(s) for s in _core.CdtStatus.__members__.values()]
        assert len(set(texts)) == len(texts)

    def test_describe_rejects_a_non_status(self) -> None:
        with pytest.raises(TypeError):
            _core.describe("Ok")


# --------------------------------------------------------------------------
# build_pslg
# --------------------------------------------------------------------------


class TestBuildPslg:
    def test_returns_a_pslg_and_no_diagnostics_on_valid_input(self, square: np.ndarray) -> None:
        result = _core.build_pslg(square, [([0, 1, 2, 3], _core.ChainRole.Outer, False)])
        assert result.ok is True
        assert result.pslg is not None
        assert list(result.diagnostics) == []

    def test_accepts_any_float64_convertible_array_like(self) -> None:
        # The design says "any (N, 2) float64-convertible array-like", so a
        # list of lists must work without the caller reaching for numpy.
        ring = [[EAST, NORTH], [EAST + 10.0, NORTH], [EAST + 10.0, NORTH + 10.0]]
        result = _core.build_pslg(ring, [([0, 1, 2], _core.ChainRole.Outer, False)])
        assert result.ok

    def test_rejects_a_vertex_array_that_is_not_n_by_2(self) -> None:
        bad = np.zeros((4, 3), dtype=np.float64)
        with pytest.raises(ValueError):
            _core.build_pslg(bad, [([0, 1, 2], _core.ChainRole.Outer, False)])

    def test_reports_failure_as_data_rather_than_raising(self, broken_build: Any) -> None:
        assert broken_build.ok is False
        assert broken_build.pslg is None

    def test_returns_the_whole_diagnostics_vector(self, broken_build: Any) -> None:
        # The load-bearing assertion of this class: four independent breakages
        # must yield at least four diagnostics with at least four distinct
        # error kinds. A binding returning only the first throws the builder's
        # entire bulk-reporting design away at the boundary.
        diagnostics = list(broken_build.diagnostics)
        assert len(diagnostics) >= 4, [d.message for d in diagnostics]
        assert len({d.error for d in diagnostics}) >= 4

    def test_every_diagnostic_carries_the_four_fields(self, broken_build: Any) -> None:
        for d in broken_build.diagnostics:
            assert d.error is not None
            assert isinstance(d.chain, int)
            assert isinstance(d.vertex, int)
            assert isinstance(d.message, str) and d.message.strip()

    def test_diagnostic_vertex_is_never_out_of_range(self, broken_build: Any) -> None:
        # pslg_builder.hpp anticipates precisely a Python binding dereferencing
        # this field: it is a valid index or the kNoVertex sentinel, never the
        # offending out-of-range value. IndexOutOfRange is in this fixture, so
        # the one diagnostic that could break the rule is present.
        for d in broken_build.diagnostics:
            assert d.vertex == NO_VERTEX or 0 <= d.vertex < 5

    def test_diagnostics_are_frozen(self, broken_build: Any) -> None:
        d = broken_build.diagnostics[0]
        with pytest.raises(AttributeError):
            d.message = "rewritten"

    def test_build_result_is_frozen(self, broken_build: Any) -> None:
        with pytest.raises(AttributeError):
            broken_build.ok = True


# --------------------------------------------------------------------------
# Pslg: opaque, read-only
# --------------------------------------------------------------------------


class TestPslg:
    def test_is_not_constructible_from_python(self) -> None:
        # There is exactly one legitimate producer and holding a Pslg is a
        # proof that validation ran. A Python constructor voids that proof.
        with pytest.raises(TypeError):
            _core.Pslg()

    def test_vertices_are_a_read_only_n_by_2_float64_view(
        self, square_pslg: Any, square: np.ndarray
    ) -> None:
        v = square_pslg.vertices
        assert v.shape == (4, 2)
        assert v.dtype == np.float64
        assert v.flags.writeable is False
        assert np.array_equal(v, square)

    def test_chain_indices_are_a_read_only_uint32_view(self, square_pslg: Any) -> None:
        idx = square_pslg.chain_indices
        assert idx.shape == (4,)
        assert idx.dtype == np.uint32
        assert idx.flags.writeable is False

    def test_chains_are_frozen_records(self, breakline_pslg: Any) -> None:
        chains = list(breakline_pslg.chains)
        assert len(chains) == 2
        outer, breakline = chains
        assert (outer.begin, outer.count) == (0, 4)
        assert outer.role == _core.ChainRole.Outer
        assert outer.is_river is False
        assert (breakline.begin, breakline.count) == (4, 2)
        assert breakline.role == _core.ChainRole.Breakline
        assert breakline.is_river is True
        with pytest.raises(AttributeError):
            outer.count = 99

    def test_indices_of_agrees_with_the_flat_buffer(self, breakline_pslg: Any) -> None:
        flat = breakline_pslg.chain_indices
        for c, chain in enumerate(breakline_pslg.chains):
            got = breakline_pslg.indices_of(c)
            assert np.array_equal(got, flat[chain.begin : chain.begin + chain.count])
            assert got.flags.writeable is False

    def test_indices_of_rejects_an_out_of_range_chain(self, square_pslg: Any) -> None:
        # The C++ precondition is a debug assert, which in a release build is
        # no precondition at all. Reached from Python it must be an exception.
        with pytest.raises(IndexError):
            square_pslg.indices_of(7)


# --------------------------------------------------------------------------
# IndexedMesh2: the array contract
# --------------------------------------------------------------------------


class TestMeshArrays:
    def test_vertices_shape_and_dtype(self, mesh: Any) -> None:
        v = mesh.vertices
        assert v.ndim == 2 and v.shape[1] == 2
        assert v.dtype == np.float64

    def test_triangles_shape_and_dtype(self, mesh: Any) -> None:
        t = mesh.triangles
        assert t.ndim == 2 and t.shape[1] == 3
        assert t.dtype == np.uint32

    def test_constrained_edges_shape_and_dtype(self, mesh: Any) -> None:
        m = mesh.constrained_edges
        assert m.shape == (mesh.triangles.shape[0],)
        assert m.dtype == np.uint8

    @pytest.mark.parametrize("name", ["vertices", "triangles", "constrained_edges"])
    def test_arrays_are_not_writeable(self, mesh: Any, name: str) -> None:
        arr = getattr(mesh, name)
        assert arr.flags.writeable is False
        with pytest.raises(ValueError):
            arr[0] = arr[0]

    @pytest.mark.parametrize("name", ["vertices", "triangles", "constrained_edges"])
    def test_arrays_are_zero_copy_views_onto_the_mesh(self, mesh: Any, name: str) -> None:
        arr = getattr(mesh, name)
        assert arr.flags.owndata is False, f"{name} copied the buffer"
        assert arr.base is not None, f"{name} has no base object to keep the owner alive"

    def test_triangle_count_and_empty_agree_with_the_arrays(self, mesh: Any) -> None:
        assert mesh.triangle_count == mesh.triangles.shape[0]
        assert mesh.empty is False
        assert isinstance(mesh.empty, bool)

    def test_masks_are_below_eight(self, mesh: Any) -> None:
        # indexed_mesh.hpp guarantee 3, restated at the boundary because it is
        # a debug assert in C++ and a release build never checks it.
        assert int(mesh.constrained_edges.max()) < 8

    def test_every_triangle_index_is_in_range(self, mesh: Any) -> None:
        assert int(mesh.triangles.max()) < mesh.vertices.shape[0]

    def test_mesh_vertices_begin_with_the_pslg_vertices(
        self, breakline_pslg: Any, mesh: Any
    ) -> None:
        # triangulate.hpp semantic obligation 1 and indexed_mesh.hpp guarantee
        # 4: index k means the same point on both sides. 6b's role join is
        # built entirely on this, and it is checkable only from here.
        n = breakline_pslg.vertices.shape[0]
        assert np.array_equal(mesh.vertices[:n], breakline_pslg.vertices)

    def test_mesh_is_not_constructible_from_python(self) -> None:
        with pytest.raises(TypeError):
            _core.IndexedMesh2()


class TestZeroCopyLifetime:
    """Risk 4: a py::array_t over an IndexedMesh2's vectors without a correct
    base object is a use-after-free whose symptom is a plausible-looking wrong
    picture. This is the single most valuable assertion in the increment."""

    def test_arrays_outlive_the_mesh_and_the_outcome(self, breakline_pslg: Any) -> None:
        outcome = _core.triangulate(breakline_pslg)
        assert outcome.ok()
        mesh = outcome.mesh
        vertices, triangles, mask = mesh.vertices, mesh.triangles, mesh.constrained_edges
        expected = (np.array(vertices), np.array(triangles), np.array(mask))

        del mesh, outcome
        gc.collect()
        # Churn the allocator: a freed buffer that is merely still mapped
        # compares equal and proves nothing, so make reuse likely.
        churn = [np.full(1 << 16, i, dtype=np.float64) for i in range(64)]
        assert len(churn) == 64

        assert np.array_equal(vertices, expected[0])
        assert np.array_equal(triangles, expected[1])
        assert np.array_equal(mask, expected[2])

    def test_each_array_holds_a_reference_to_its_owner(self, mesh: Any) -> None:
        # The mechanism behind the test above, asserted directly so that a
        # failure says "no keep-alive" rather than "the numbers changed".
        for name in ("vertices", "triangles", "constrained_edges"):
            arr = getattr(mesh, name)
            assert arr.base is not None
            assert sys.getrefcount(arr.base) > 1


# --------------------------------------------------------------------------
# triangulate
# --------------------------------------------------------------------------


class TestTriangulate:
    def test_succeeds_on_a_valid_pslg(self, breakline_pslg: Any) -> None:
        outcome = _core.triangulate(breakline_pslg)
        assert outcome.status == _core.CdtStatus.Ok
        assert outcome.ok() is True
        assert outcome.message == ""
        assert outcome.mesh.triangle_count > 0

    def test_delaunay_flag_is_reachable_in_both_settings(self, breakline_pslg: Any) -> None:
        # One bound bool, and the cheapest intuition in the increment: the user
        # gets to see both triangulations of the same input.
        on = _core.triangulate(breakline_pslg, delaunay=True)
        off = _core.triangulate(breakline_pslg, delaunay=False)
        assert on.ok() and off.ok()
        # The vertex set is fixed, so flipping cannot change the counts; the
        # two triangulations may legitimately differ edge for edge.
        assert on.mesh.vertices.shape == off.mesh.vertices.shape
        assert on.mesh.triangle_count == off.mesh.triangle_count

    def test_delaunay_defaults_to_true(self, breakline_pslg: Any) -> None:
        default = _core.triangulate(breakline_pslg)
        explicit = _core.triangulate(breakline_pslg, delaunay=True)
        assert np.array_equal(default.mesh.triangles, explicit.mesh.triangles)

    def test_reports_a_non_noded_input_as_a_failure_status(self, crossing_pslg: Any) -> None:
        # Which non-Ok status the backend picks is its own classification and
        # is not pinned here; that it fails, says something, and hands back no
        # mesh is the contract, and it is what 6b's failure presentation draws.
        outcome = _core.triangulate(crossing_pslg)
        assert outcome.status != _core.CdtStatus.Ok
        assert outcome.ok() is False
        assert outcome.message.strip()
        assert _core.describe(outcome.status).strip()

    @pytest.mark.parametrize("pslg_name", ["breakline_pslg", "crossing_pslg"])
    def test_outcome_never_pairs_ok_with_an_empty_mesh(
        self, request: pytest.FixtureRequest, pslg_name: str
    ) -> None:
        # testing.md, `cdt`: the silent mode is "forgot to add the outline",
        # which triangulates successfully to zero interior triangles. The
        # wrapper asserts this in debug only, so the boundary restates it.
        outcome = _core.triangulate(request.getfixturevalue(pslg_name))
        assert outcome.ok() == (outcome.mesh.triangle_count > 0)
        assert outcome.ok() == (outcome.message == "")

    def test_outcome_fields_are_read_only(self, breakline_pslg: Any) -> None:
        outcome = _core.triangulate(breakline_pslg)
        with pytest.raises(AttributeError):
            outcome.status = _core.CdtStatus.BackendFailure

    def test_rejects_anything_that_is_not_a_pslg(self, square: np.ndarray) -> None:
        with pytest.raises(TypeError):
            _core.triangulate(square)

    def test_concurrent_triangulation_of_one_pslg_agrees(self, breakline_pslg: Any) -> None:
        # A const Pslg is safe for concurrent read (testing.md, core geometry
        # -- PSLG). With the GIL released this stops being theoretical: two
        # threads genuinely execute the backend at once, and a binding that
        # caches or mutates through the const handle diverges here.
        results: list[Any] = []
        barrier = threading.Barrier(4)

        def run() -> None:
            barrier.wait(timeout=30.0)
            results.append(_core.triangulate(breakline_pslg).mesh.triangles)

        threads = [threading.Thread(target=run) for _ in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30.0)
            assert not t.is_alive()
        assert len(results) == 4
        for other in results[1:]:
            assert np.array_equal(results[0], other)


# --------------------------------------------------------------------------
# The GIL
# --------------------------------------------------------------------------


class _Ticker:
    """A thread that ticks once per millisecond, and cannot tick without the GIL.

    `time.sleep` releases the GIL and must re-acquire it to continue, so a tick
    is direct evidence that the GIL was available. The tick rate is bounded by
    the sleep, so the ticker never starves the thread it is observing.
    """

    def __init__(self) -> None:
        self.ticks = 0
        self._stop = threading.Event()
        self._ready = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self) -> None:
        self._ready.set()
        while not self._stop.is_set():
            self.ticks += 1
            time.sleep(0.001)

    def __enter__(self) -> _Ticker:
        self._thread.start()
        assert self._ready.wait(timeout=10.0), "ticker never started"
        return self

    def __exit__(self, *exc: object) -> None:
        self._stop.set()
        self._thread.join(timeout=10.0)


def ticks_during(call: Callable[[], T]) -> tuple[T, int, float]:
    """Return the call's result, the ticks observed during it, and its duration."""
    with _Ticker() as ticker:
        time.sleep(0.02)  # let the ticker reach steady state
        before = ticker.ticks
        started = time.perf_counter()
        result = call()
        elapsed = time.perf_counter() - started
        after = ticker.ticks
    return result, after - before, elapsed


# Thresholds from the control below, measured on CPython 3.14 (arm64): a
# GIL-holding C call of 0.67 s let through 1 tick; a releasing call of 0.30 s
# let through ~220. The gap is two orders of magnitude, so the thresholds do
# not need to be tight.
HELD_TICKS = 3
RELEASED_TICKS = 20


@pytest.fixture(scope="module")
def large_pslg() -> Any:
    """Enough work for the ticker to resolve. A jittered lattice strictly
    inside the outer ring: deterministic (a fixed seed), no duplicate points,
    and no point on a constraint edge."""
    rng = np.random.default_rng(20260918)
    side = 380
    step = np.linspace(1.0, 99.0, side)
    gx, gy = np.meshgrid(step, step * 0.8)
    jitter = rng.uniform(-0.2, 0.2, size=(2, side, side))
    interior = np.column_stack(
        [(gx + jitter[0]).ravel() + EAST, (gy + jitter[1]).ravel() + NORTH]
    )
    ring = np.array(
        [
            [EAST, NORTH],
            [EAST + 100.0, NORTH],
            [EAST + 100.0, NORTH + 80.0],
            [EAST, NORTH + 80.0],
        ]
    )
    result = _core.build_pslg(
        np.vstack([ring, interior]), [([0, 1, 2, 3], _core.ChainRole.Outer, False)]
    )
    assert result.ok, [d.message for d in result.diagnostics]
    return result.pslg


class TestGil:
    def test_the_ticker_distinguishes_a_held_gil_from_a_released_one(self) -> None:
        """The control. Without it a pass and a probe that never ran look the
        same, and the GIL test below would measure nothing.

        `sorted` on a shuffled list is a single C call that holds the GIL for
        its whole duration -- the broken X. `time.sleep` is one that releases.
        If this fails, the GIL test below is uninterpretable, whatever it says;
        a free-threaded build is the likeliest cause, and on such a build the
        whole question has to be re-asked rather than silently passed.
        """
        import random

        data = list(range(3_000_000))
        random.Random(12345).shuffle(data)

        _, held, held_elapsed = ticks_during(lambda: sorted(data))
        assert held_elapsed >= 0.05, f"control workload too fast ({held_elapsed:.3f}s)"
        assert held <= HELD_TICKS, f"a GIL-holding call let {held} ticks through"

        _, released, _ = ticks_during(lambda: time.sleep(0.3))
        assert released >= RELEASED_TICKS, f"a GIL-releasing call let only {released} ticks through"

    def test_triangulate_releases_the_gil(self, large_pslg: Any) -> None:
        outcome, ticks, elapsed = ticks_during(lambda: _core.triangulate(large_pslg))
        assert outcome.ok(), outcome.message
        # A probe able to fail: if the workload is too fast there is nothing to
        # observe, and this must say so rather than pass on an empty window.
        assert elapsed >= 0.05, (
            f"triangulation took only {elapsed:.3f}s -- too fast to measure GIL release; "
            "enlarge large_pslg rather than lowering this bound"
        )
        assert ticks >= RELEASED_TICKS, (
            f"only {ticks} ticks in {elapsed:.3f}s: triangulate appears to hold the GIL"
        )


# --------------------------------------------------------------------------
# The firewall: what must never cross
# --------------------------------------------------------------------------


class TestBoundary:
    """Risk 3, binding-surface gravity. The surface in the design is
    exhaustive; an addition is a design change, not a line in a PR. Asserting
    absence is cheap and it is what stops the boundary eroding by accretion."""

    PROHIBITED = (
        "PslgBuilder",
        "Segment2",
        "IndexedRing",
        "PointRing",
        "SnapGrid",
        "Box2",
    )

    # Substrings rather than names, because the thing being excluded is a
    # capability and not one spelling of it. Matched case-insensitively against
    # every public name the module exports.
    FORBIDDEN_SUBSTRINGS = (
        "detria",
        "kernel",
        "path",
        "file",
        "crs",
        "proj",
        "epsg",
        "adjacen",
        "neighbo",
    )

    def public_names(self) -> list[str]:
        return [n for n in dir(_core) if not n.startswith("_")]

    def test_no_prohibited_type_is_exported(self) -> None:
        exported = set(self.public_names())
        assert exported.isdisjoint(self.PROHIBITED), exported & set(self.PROHIBITED)

    def test_no_exported_name_suggests_io_crs_or_topology(self) -> None:
        offenders = [
            name
            for name in self.public_names()
            for bad in self.FORBIDDEN_SUBSTRINGS
            if bad in name.lower()
        ]
        assert offenders == []

    @pytest.mark.parametrize("obj_name", ["Pslg", "IndexedMesh2", "CdtOutcome"])
    def test_no_bound_type_exposes_io_crs_or_topology(self, obj_name: str) -> None:
        cls = getattr(_core, obj_name)
        offenders = [
            name
            for name in dir(cls)
            if not name.startswith("_")
            for bad in self.FORBIDDEN_SUBSTRINGS
            if bad in name.lower()
        ]
        assert offenders == []

    def test_build_pslg_does_not_accept_a_path(self, tmp_path: Path) -> None:
        # File decoding is Python's; the core never sees a path. A str or a
        # Path reaching build_pslg must be a TypeError and never an open().
        target = tmp_path / "outline.geojson"
        target.write_text("{}")
        for candidate in (str(target), target):
            with pytest.raises(TypeError):
                _core.build_pslg(candidate, [])

    def test_mesh_exposes_no_mutable_accessor(self, mesh: Any) -> None:
        # The three arrays are covered by their own writeable tests; this
        # catches a *fourth* accessor added later that forgets to be read-only.
        for name in dir(mesh):
            if name.startswith("_"):
                continue
            value = getattr(mesh, name)
            if isinstance(value, np.ndarray):
                assert value.flags.writeable is False, f"{name} is writeable"


# --------------------------------------------------------------------------
# Documentation and stubs
# --------------------------------------------------------------------------


NEW_SURFACE = (
    "ChainRole",
    "CdtStatus",
    "Pslg",
    "IndexedMesh2",
    "CdtOutcome",
    "PslgBuildResult",
    "build_pslg",
    "triangulate",
    "describe",
)


class TestDocumentation:
    """CLAUDE.md section 4 requires every pybind11-exposed surface to be
    documented; `tests/python/test_core.py` holds the same line for Point2."""

    @pytest.mark.parametrize("name", NEW_SURFACE)
    def test_every_bound_name_carries_a_docstring(self, name: str) -> None:
        assert getattr(_core, name).__doc__

    def test_the_diagnostic_vertex_sentinel_is_documented_where_python_sees_it(self) -> None:
        # The header's warning about dereferencing `vertex` was written for a
        # Python binding. It has to cross, or it protects nobody.
        text = STUB.read_text(encoding="utf-8")
        assert "kNoVertex" in text or "sentinel" in text


class TestStubs:
    """`_core.pyi` is 6a's own deliverable and mypy checks its consistency, not
    its completeness: a name absent from both the stub and the caller's code
    passes strict mode silently."""

    def test_stub_file_exists(self) -> None:
        assert STUB.is_file()

    @pytest.mark.parametrize("name", NEW_SURFACE)
    def test_stub_declares_the_new_surface(self, name: str) -> None:
        tree = ast.parse(STUB.read_text(encoding="utf-8"))
        declared = {
            node.name
            for node in tree.body
            if isinstance(node, ast.ClassDef | ast.FunctionDef)
        }
        assert name in declared
