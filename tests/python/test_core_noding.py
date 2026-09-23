"""The noder's binding surface: increment 5c.

Committed RED, before `bindings/core.cpp` exposes any of it -- see
`docs/increments/README.md`, "The loop". Every name reached here is looked up
on `_core` *inside* a test or a fixture rather than imported at module scope,
so the failure mode is a named `AttributeError` on the symbol that is missing
and pytest still collects the file. A collection error would hide the next real
problem behind an import.

**Rebuild before trusting a green run in this file.** `pytest` does not rebuild
the extension on this machine -- scikit-build-core's editable auto-rebuild needs
`ninja`, which is absent, so it silently no-ops and the run exercises the
previously installed `.so`. A green handback from a session that never rebuilt
is byte-identical to one that did (`.claude/REQUIRED-READING.md`, measured):

    cmake --build build-pyext -j --target _core
    cp build-pyext/_core.cpython-*-darwin.so .venv/lib/python3.*/site-packages/tin_engine/

This is 5c's **invariant-critical** suite (`docs/increments/05c-noder-wiring.md`,
"What is worth testing"), and the reason is that everything the increment adds
to the binding layer is a lifetime or a reinterpret:

* `edge_properties` is a `(E,)` uint32 view reinterpreted over a
  `std::vector<EdgeProperties>`. A wrong stride or a wrong length reads
  plausible garbage -- masks that are *almost* right are exactly the failure a
  picture cannot show, so the length is checked against an oracle derived from
  the chains' own roles rather than from the array itself;
* every view needs its base object, or dropping the `NodeOutcome` leaves Python
  holding freed memory (`06-cdt-viewer.md` risk 4);
* the GIL is released around `node()`, which is licensed by `node.hpp:5-11`'s
  purity claim and by nothing else, so the purity is measured here too.

What this file does **not** own: the CLI's presentation of a refusal
(`test_cli_draw.py`), the structural join onto `viz/protocols.py`
(`test_viz_protocols.py`), and `triangulate`'s own surface (`test_core_cdt.py`,
which pins that it no longer accepts a `Pslg` at all).
"""

from __future__ import annotations

import ast
import gc
import sys
import threading
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from tin_engine import _core
from tin_engine.features import DEFAULT_VOCABULARY

REPO_ROOT = Path(__file__).resolve().parents[2]
STUB = REPO_ROOT / "src_python" / "tin_engine" / "_core.pyi"

# The ordering ruling again: a UTM 33N-shaped domain, not the unit square.
# Snapping is magnitude-sensitive in a way nothing else in the project is --
# `kMaxGridIndex` is compared against `coordinate / spacing` -- so a noding
# suite in [0,1]^2 would exercise the one case that cannot overflow.
EAST = 430_000.0
NORTH = 6_900_000.0
ORIGIN = np.array([EAST, NORTH], dtype=np.float64)

# The composition root's default, restated as a local constant rather than
# imported from `cli.py`: this suite tests the binding, and taking the number
# from the module under test one layer up would make the two move together.
SPACING = 1e-3

NO_PROPERTIES = 0
RIVER = DEFAULT_VOCABULARY.mask("river")
ROAD = DEFAULT_VOCABULARY.mask("road")

NODE_STATUS_NAMES = (
    "Ok",
    "NotRun",
    "InvalidSnapSpacing",
    "CoordinateOutOfRange",
    "RingCollapsed",
    "RingDegenerateAfterSnap",
    "NonSimpleRing",
    "NotConverged",
    "MalformedOutput",
)

NEW_SURFACE = ("NodeStatus", "NodedPslg", "NodeOutcome", "node")


def build(vertices: np.ndarray, chains: list[Any]) -> Any:
    """A validated `Pslg`, or a failure that names why at the call site."""
    result = _core.build_pslg(vertices, chains)
    assert result.ok, [d.message for d in result.diagnostics]
    return result.pslg


def noded(pslg: Any, spacing: float = SPACING) -> Any:
    outcome = _core.node(pslg, spacing)
    assert outcome.ok(), f"{outcome.status}: {outcome.message}"
    return outcome.pslg


def expected_edge_count_of(graph: Any, c: int) -> int:
    chain = graph.chains[c]
    closed = (_core.ChainRole.Outer, _core.ChainRole.Hole)
    return chain.count if chain.role in closed else chain.count - 1


def expected_edge_count(graph: Any) -> int:
    """How many edges the graph has, from the chains' roles.

    The oracle for `edge_properties`' length, and it is deliberately not the
    array's own `size`: a ring carries a closing edge its index run does not
    store and an open breakline does not, so this is the one arithmetic in the
    increment that a reinterpret can get wrong without the picture showing it.
    """
    return sum(expected_edge_count_of(graph, c) for c in range(len(graph.chains)))


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------


@pytest.fixture
def square() -> np.ndarray:
    """A counterclockwise outer ring at UTM33 magnitudes, closure not stored."""
    return np.array(
        [[0.0, 0.0], [700.0, 0.0], [700.0, 700.0], [0.0, 700.0]], dtype=np.float64
    ) + ORIGIN


@pytest.fixture
def crossing_pslg(square: np.ndarray) -> Any:
    """The increment's sentence, as a PSLG: a road crossing a river.

    The same geometry as the `road-crosses-river` gallery fixture and as
    `test_noding_node.cpp:162-165`, at the same UTM magnitudes. The two
    breaklines meet at (350, 350) + ORIGIN, a point neither chain names, which
    is what "not noded" means -- and after `node()` it is a vertex of both.
    """
    vertices = np.vstack(
        [
            square,
            np.array([[100.0, 100.0], [600.0, 600.0], [100.0, 600.0], [600.0, 100.0]])
            + ORIGIN,
        ]
    )
    return build(
        vertices,
        [
            ([0, 1, 2, 3], _core.ChainRole.Outer, NO_PROPERTIES),
            ([4, 5], _core.ChainRole.Breakline, ROAD),
            ([6, 7], _core.ChainRole.Breakline, RIVER),
        ],
    )


@pytest.fixture
def disjoint_pslg(square: np.ndarray) -> Any:
    """An outer ring and one breakline that touches nothing: noding it must be
    a no-op on the topology, which is the control for every claim below that
    noding *changed* something."""
    vertices = np.vstack(
        [square, np.array([[100.0, 350.0], [600.0, 350.0]]) + ORIGIN]
    )
    return build(
        vertices,
        [
            ([0, 1, 2, 3], _core.ChainRole.Outer, NO_PROPERTIES),
            ([4, 5], _core.ChainRole.Breakline, RIVER),
        ],
    )


@pytest.fixture
def noded_crossing(crossing_pslg: Any) -> Any:
    return noded(crossing_pslg)


# --------------------------------------------------------------------------
# NodeStatus, and `describe` as an overload set
# --------------------------------------------------------------------------


class TestNodeStatus:
    def test_has_exactly_the_nine_statuses(self) -> None:
        assert set(_core.NodeStatus.__members__) == set(NODE_STATUS_NAMES)

    def test_ok_is_zero_and_every_status_round_trips(self) -> None:
        # Ok is 0 in BOTH enumerations, which is the whole hazard the overload
        # tests below exist for.
        assert _core.NodeStatus.Ok.value == 0
        for status in _core.NodeStatus.__members__.values():
            assert _core.NodeStatus(status.value) == status

    def test_describe_reaches_python_for_every_status(self) -> None:
        # Including the self-checks: `cli.py`'s band prints exactly this text,
        # so a status whose description never crossed is a blank band on the
        # one run where the user needs words.
        for name in NODE_STATUS_NAMES:
            text = _core.describe(getattr(_core.NodeStatus, name))
            assert isinstance(text, str)
            assert text.strip(), f"{name} describes to nothing"

    def test_descriptions_are_distinct(self) -> None:
        # Wording is not pinned -- testing.md forbids pinning prose -- but a
        # switch arm falling through to a neighbour's text is a real defect and
        # distinctness is what catches it.
        texts = [_core.describe(s) for s in _core.NodeStatus.__members__.values()]
        assert len(set(texts)) == len(texts)

    def test_a_status_is_not_equal_to_a_cdt_status_of_the_same_value(self) -> None:
        # The premise of the overload ruling, asserted rather than assumed:
        # `py::enum_` registers a distinct Python type per enumeration, so the
        # two zeros are different objects even though both are integer 0.
        assert _core.NodeStatus.Ok != _core.CdtStatus.Ok
        assert _core.NodeStatus.Ok.value == _core.CdtStatus.Ok.value


class TestDescribeOverloadSet:
    """`describe` becomes one name over two enumerations whose `Ok`s are both
    integer 0 (`05c-noder-wiring.md`, "describe becomes an overload set").

    The design's claim is that pybind11's first resolution pass runs every
    overload with conversions disabled, so each call lands on its own
    enumeration exactly. That claim is the thing pinned here: a binding that
    got it wrong returns the *other* enumeration's sentence, which is a
    plausible English sentence in the header band and wrong.
    """

    def test_each_enumeration_gets_its_own_sentence_for_zero(self) -> None:
        assert _core.describe(_core.NodeStatus.Ok) != _core.describe(_core.CdtStatus.Ok)

    def test_no_node_status_describes_to_a_cdt_status_sentence(self) -> None:
        # The stronger form: not merely that the two zeros differ, but that no
        # value of either enumeration is answered out of the other's table.
        cdt = {_core.describe(s) for s in _core.CdtStatus.__members__.values()}
        node = {_core.describe(s) for s in _core.NodeStatus.__members__.values()}
        assert cdt.isdisjoint(node)

    def test_the_cdt_overload_still_answers_every_cdt_status(self) -> None:
        # Adding an overload must not shadow the shipped one for any value.
        for status in _core.CdtStatus.__members__.values():
            assert _core.describe(status).strip()

    @pytest.mark.parametrize("bad", ["Ok", 0, 1.0, None])
    def test_describe_rejects_anything_that_is_not_a_status(self, bad: object) -> None:
        # `describe("Ok")` was a TypeError with one overload and must stay one
        # with two: a second overload is where an accidental implicit
        # conversion from `int` or `str` would first become reachable.
        with pytest.raises(TypeError):
            _core.describe(bad)

    def test_describe_rejects_a_third_enumeration(self) -> None:
        # `ChainRole.Outer` is also integer 0 and has no `describe`. Measured
        # on this tree before the second overload existed: a TypeError. If it
        # starts answering, the no-convert pass is not doing what the design
        # says it does.
        with pytest.raises(TypeError):
            _core.describe(_core.ChainRole.Outer)


# --------------------------------------------------------------------------
# node(): what it produces, and what it refuses
# --------------------------------------------------------------------------


class TestNode:
    def test_nodes_a_crossing_that_triangulate_used_to_refuse(
        self, crossing_pslg: Any
    ) -> None:
        outcome = _core.node(crossing_pslg, SPACING)
        assert outcome.status == _core.NodeStatus.Ok
        assert outcome.ok() is True
        assert outcome.message == ""
        assert outcome.pslg is not None

    def test_the_crossing_point_becomes_a_vertex_of_both_chains(
        self, crossing_pslg: Any, noded_crossing: Any
    ) -> None:
        # The increment's product, stated as geometry rather than as counts:
        # (350, 350) + ORIGIN is named by no input chain and must be a node of
        # both breaklines afterwards. The tolerance is the producer's own
        # relation -- a snapped node is within half a cell of the true
        # intersection, not on it.
        meeting = np.array([350.0, 350.0]) + ORIGIN
        vertices = np.asarray(noded_crossing.vertices)
        close = np.flatnonzero(np.abs(vertices - meeting).max(axis=1) <= SPACING)
        assert close.size == 1, "the crossing is not a single node"
        crossing_id = int(close[0])

        breaklines = [
            c
            for c, chain in enumerate(noded_crossing.chains)
            if chain.role == _core.ChainRole.Breakline
        ]
        assert len(breaklines) == 2
        for c in breaklines:
            assert crossing_id in set(np.asarray(noded_crossing.indices_of(c)).tolist())

        # Able to fail on its own: the input named no such vertex, so a `node`
        # that returned its argument unchanged would fail the assertion above
        # rather than passing it vacuously.
        original = np.asarray(crossing_pslg.vertices)
        assert np.abs(original - meeting).max(axis=1).min() > SPACING

    def test_each_breakline_gains_exactly_one_split(self, noded_crossing: Any) -> None:
        # Two segments crossing once become two segments of two edges each.
        breaklines = [
            c
            for c, chain in enumerate(noded_crossing.chains)
            if chain.role == _core.ChainRole.Breakline
        ]
        assert [noded_crossing.chains[c].count for c in breaklines] == [3, 3]

    def test_a_disjoint_input_keeps_its_chain_structure(self, disjoint_pslg: Any) -> None:
        # The control: noding input that needs no noding must not invent nodes.
        # Without it every assertion above is satisfied by a `node` that
        # splits everything at everything.
        out = noded(disjoint_pslg)
        assert [chain.count for chain in out.chains] == [
            chain.count for chain in disjoint_pslg.chains
        ]
        assert out.vertices.shape == disjoint_pslg.vertices.shape

    def test_max_rounds_defaults_to_four_and_is_reachable(self, crossing_pslg: Any) -> None:
        # The binding reproduces the C++ default rather than choosing one; that
        # it is spelled at all is what this checks, and that spelling it makes
        # no difference on input that converges in one round.
        default = _core.node(crossing_pslg, SPACING)
        explicit = _core.node(crossing_pslg, SPACING, max_rounds=4)
        assert default.status == explicit.status
        assert np.array_equal(default.pslg.vertices, explicit.pslg.vertices)

    def test_spacing_has_no_default(self, crossing_pslg: Any) -> None:
        # `NodeOptions::spacing` has none "EVER" (`node.hpp:83-89`): the right
        # value is a policy question and the binding is not the composition
        # root. `cli.py` defaults; this must not.
        with pytest.raises(TypeError):
            _core.node(crossing_pslg)

    def test_is_deterministic_across_repeated_calls(self, crossing_pslg: Any) -> None:
        first = noded(crossing_pslg)
        second = noded(crossing_pslg)
        assert np.array_equal(first.vertices, second.vertices)
        assert np.array_equal(first.chain_indices, second.chain_indices)
        assert np.array_equal(first.edge_properties, second.edge_properties)


class TestNodeRefusalsAreData:
    """A refused noding is a fact about the terrain and about the spacing, so
    it comes back as a status; only a mis-shaped argument raises. This is
    `build_pslg`'s rule, applied to the second engine entry point."""

    @pytest.mark.parametrize("spacing", [0.0, -1.0, float("nan"), float("inf")])
    def test_a_bad_spacing_is_a_status_and_not_an_exception(
        self, crossing_pslg: Any, spacing: float
    ) -> None:
        # The C++ half of the two refusal channels the design keeps apart: the
        # CLI refuses `--snap-spacing` as a usage error before drawing
        # anything, and the library refuses it as a precondition for callers
        # that are not the CLI. Neither is removable in favour of the other.
        outcome = _core.node(crossing_pslg, spacing)
        assert outcome.status == _core.NodeStatus.InvalidSnapSpacing
        assert outcome.ok() is False
        assert outcome.pslg is None
        assert outcome.message.strip()

    def test_a_spacing_too_fine_for_the_coordinates_is_a_status(
        self, crossing_pslg: Any
    ) -> None:
        # `kMaxGridIndex` is 2**51, so at spacing 1e-12 a UTM easting overflows
        # the lattice. This is the lower bound that makes `--snap-spacing` an
        # option rather than a hidden constant, and it is a diagnosis with a
        # named lever: coarsen, or re-project.
        outcome = _core.node(crossing_pslg, 1e-12)
        assert outcome.status == _core.NodeStatus.CoordinateOutOfRange
        assert outcome.pslg is None
        assert outcome.message.strip()

    def test_the_default_spacing_is_inside_the_representable_range(
        self, crossing_pslg: Any
    ) -> None:
        # A probe able to fail: the test above proves nothing about the default
        # unless the default is known to be on the other side of the boundary.
        assert _core.node(crossing_pslg, SPACING).status == _core.NodeStatus.Ok

    def test_never_returns_a_pslg_alongside_a_failure(self, crossing_pslg: Any) -> None:
        for spacing in (0.0, 1e-12, SPACING):
            outcome = _core.node(crossing_pslg, spacing)
            assert outcome.ok() == (outcome.pslg is not None)
            assert outcome.ok() == (outcome.message == "")

    def test_rejects_anything_that_is_not_a_pslg(self, square: np.ndarray) -> None:
        with pytest.raises(TypeError):
            _core.node(square, SPACING)

    def test_does_not_accept_a_path(self, tmp_path: Path) -> None:
        # File decoding is Python's; the core never sees a path. A `str` or a
        # `Path` reaching `node` must be a TypeError and never an `open()`.
        target = tmp_path / "outline.geojson"
        target.write_text("{}")
        for candidate in (str(target), target):
            with pytest.raises(TypeError):
                _core.node(candidate, SPACING)

    def test_rejects_a_noded_pslg(self, noded_crossing: Any) -> None:
        # Noding twice is a caller error, not a status: `node` takes the
        # un-noded graph, and the type system is where that is said.
        with pytest.raises(TypeError):
            _core.node(noded_crossing, SPACING)


# --------------------------------------------------------------------------
# NodedPslg: the four shared accessors, and the three it adds
# --------------------------------------------------------------------------


class TestNodedPslgShape:
    """`PslgLike` is a structural contract with three implementations --
    `Pslg`, `NodedPslg` and `viz.fixtures.Fixture`. The four shared accessors
    are bound by one function template precisely so that a member cannot be
    added to one and forgotten on the other."""

    @pytest.mark.parametrize("name", ["vertices", "chains", "chain_indices", "indices_of"])
    def test_carries_every_pslg_like_member(self, noded_crossing: Any, name: str) -> None:
        assert hasattr(noded_crossing, name)

    def test_vertices_are_an_n_by_two_float_array(self, noded_crossing: Any) -> None:
        v = noded_crossing.vertices
        assert v.ndim == 2 and v.shape[1] == 2
        assert v.dtype == np.float64

    def test_chain_indices_are_uint32(self, noded_crossing: Any) -> None:
        assert noded_crossing.chain_indices.dtype == np.uint32

    def test_indices_of_is_the_chains_own_run(self, noded_crossing: Any) -> None:
        flat = np.asarray(noded_crossing.chain_indices)
        for c, chain in enumerate(noded_crossing.chains):
            run = np.asarray(noded_crossing.indices_of(c))
            assert np.array_equal(run, flat[chain.begin : chain.begin + chain.count])

    def test_indices_of_guards_its_range(self, noded_crossing: Any) -> None:
        # The same `IndexError` guard `Pslg.indices_of` has, replacing a debug
        # assert that a release build never checks.
        with pytest.raises(IndexError):
            noded_crossing.indices_of(len(noded_crossing.chains))

    def test_every_index_is_in_range(self, noded_crossing: Any) -> None:
        assert int(np.asarray(noded_crossing.chain_indices).max()) < (
            noded_crossing.vertices.shape[0]
        )

    def test_roles_survive_noding(self, crossing_pslg: Any, noded_crossing: Any) -> None:
        # `cli.py`'s `closed_roles` is keyed on these, and they are `ChainRole`
        # values rather than the fixtures' strings -- which is the whole reason
        # the source and its vocabulary travel together.
        assert [chain.role for chain in noded_crossing.chains] == [
            chain.role for chain in crossing_pslg.chains
        ]

    def test_is_not_constructible_from_python(self) -> None:
        # The private-constructor proof (`core/noded_pslg.hpp:145-148`) must
        # not have a Python back door: constructing one here would be
        # constructing an unverified `NodedPslg`.
        with pytest.raises(TypeError):
            _core.NodedPslg()


class TestGridSpacing:
    def test_reports_the_spacing_it_was_handed(self, crossing_pslg: Any) -> None:
        for spacing in (SPACING, 0.1, 1.0):
            assert noded(crossing_pslg, spacing).grid_spacing == pytest.approx(spacing)

    def test_the_grid_itself_does_not_cross(self, noded_crossing: Any) -> None:
        # `SnapGrid` is not bound: a grid vocabulary in Python next to no
        # consumer. `grid_spacing` is the one value a caller needs and it is
        # the value it handed in.
        assert not hasattr(noded_crossing, "grid")


class TestEdgeProperties:
    def test_is_a_flat_uint32_array(self, noded_crossing: Any) -> None:
        # A view of bare masks, not of a bound `EdgeProperties`: Python already
        # has an integer with `|`, `&` and `bit_count()`, and `Chain.properties`
        # already crosses that way.
        props = noded_crossing.edge_properties
        assert props.ndim == 1
        assert props.dtype == np.uint32

    def test_holds_one_entry_per_output_edge(self, noded_crossing: Any) -> None:
        # The oracle is the roles' own arithmetic, not the array's size: a
        # ring carries a closing edge and an open breakline does not. A
        # reinterpret with the wrong length or stride fails here and nowhere
        # else, because the values it reads are plausible either way.
        assert noded_crossing.edge_properties.shape == (
            expected_edge_count(noded_crossing),
        )

    def test_carries_both_input_features_after_the_split(
        self, noded_crossing: Any
    ) -> None:
        # The picture the increment exists to draw: a road and a river, in two
        # colours, meeting at a constructed node. Both bits must survive, and
        # neither may be smeared over the other.
        masks = [int(m) for m in noded_crossing.edge_properties]
        assert sum(1 for m in masks if m & ROAD) >= 2
        assert sum(1 for m in masks if m & RIVER) >= 2
        assert all(m in (NO_PROPERTIES, ROAD, RIVER) for m in masks)

    def test_the_outer_rings_edges_are_unclassified(self, noded_crossing: Any) -> None:
        # The empty set is a legal value meaning *unclassified*, and the outer
        # ring was built with it. A mask smeared across chains shows up here.
        outer = next(
            c
            for c, chain in enumerate(noded_crossing.chains)
            if chain.role == _core.ChainRole.Outer
        )
        base = sum(
            expected_edge_count_of(noded_crossing, k) for k in range(outer)
        )
        count = expected_edge_count_of(noded_crossing, outer)
        masks = np.asarray(noded_crossing.edge_properties)[base : base + count]
        assert np.all(masks == NO_PROPERTIES)

    def test_two_features_within_one_cell_merge_into_one_edge(
        self, square: np.ndarray
    ) -> None:
        """A road running along a river: one output edge carrying both bits.

        The gallery cannot show this -- `Fixture` has no spacing field, and a
        per-fixture spacing would put a policy in the module that holds none --
        so the merge case is reachable only from here, by calling `node()`
        with a spacing coarse enough to bring the two within half a cell. That
        is the design's deferred fixture, exercised at the one layer where it
        costs nothing.
        """
        vertices = np.vstack(
            [
                square,
                np.array([[100.0, 350.0], [600.0, 350.0]]) + ORIGIN,
                np.array([[100.0, 350.004], [600.0, 350.004]]) + ORIGIN,
            ]
        )
        pslg = build(
            vertices,
            [
                ([0, 1, 2, 3], _core.ChainRole.Outer, NO_PROPERTIES),
                ([4, 5], _core.ChainRole.Breakline, RIVER),
                ([6, 7], _core.ChainRole.Breakline, ROAD),
            ],
        )
        merged = noded(pslg, 0.1)
        masks = {int(m) for m in merged.edge_properties}
        assert (RIVER | ROAD) in masks

        # Able to fail on its own: at the default spacing the two are four
        # millimetres apart and must stay two distinct edges, which is what
        # makes the default a default rather than a merge.
        apart = {int(m) for m in noded(pslg, SPACING).edge_properties}
        assert (RIVER | ROAD) not in apart

    def test_edge_base_does_not_cross(self, noded_crossing: Any) -> None:
        # Recorded as a decision rather than left to read as an oversight
        # (`05c-noder-wiring.md`, "NodedPslg's binding"): a Python caller holds
        # the dense array and `chains`, so `edge_base` is a prefix sum with no
        # consumer. It is one line to add when one turns up.
        assert not hasattr(noded_crossing, "edge_base")


class TestNodeOfInputVertex:
    def test_is_total_over_the_input_vertices(
        self, crossing_pslg: Any, noded_crossing: Any
    ) -> None:
        # Guarantee 9's replacement, and the only thing that answers "where did
        # vertex 4 go" without a coordinate search. Total over the INPUT array,
        # unreferenced vertices included.
        mapping = noded_crossing.node_of_input_vertex
        assert mapping.shape == (crossing_pslg.vertices.shape[0],)
        assert mapping.dtype == np.uint32

    def test_every_entry_is_a_node_of_the_noded_graph(self, noded_crossing: Any) -> None:
        mapping = np.asarray(noded_crossing.node_of_input_vertex)
        assert int(mapping.max()) < noded_crossing.vertices.shape[0]

    def test_each_input_vertex_maps_to_its_own_snapped_image(
        self, crossing_pslg: Any, noded_crossing: Any
    ) -> None:
        # The producer's relation, not the prose one: a snapped node is NEAR
        # its input vertex, never on it, so an exact-equality oracle would be
        # red on correct output. Half a cell per axis is the snap's own bound.
        original = np.asarray(crossing_pslg.vertices)
        nodes = np.asarray(noded_crossing.vertices)
        mapping = np.asarray(noded_crossing.node_of_input_vertex)
        offset = np.abs(nodes[mapping] - original).max()
        assert offset <= SPACING / 2.0 + SPACING * 1e-9


# --------------------------------------------------------------------------
# The arrays: writeability and lifetime
# --------------------------------------------------------------------------


NODED_ARRAYS = ("vertices", "chain_indices", "edge_properties", "node_of_input_vertex")


class TestArrayContract:
    @pytest.mark.parametrize("name", NODED_ARRAYS)
    def test_arrays_are_not_writeable(self, noded_crossing: Any, name: str) -> None:
        # A mutant that forgets this lets Python write through into a type
        # whose entire contract is that it was verified once and not touched
        # since.
        arr = getattr(noded_crossing, name)
        assert arr.flags.writeable is False
        with pytest.raises(ValueError):
            arr[0] = arr[0]

    @pytest.mark.parametrize("name", NODED_ARRAYS)
    def test_arrays_are_zero_copy_views(self, noded_crossing: Any, name: str) -> None:
        arr = getattr(noded_crossing, name)
        assert arr.flags.owndata is False, f"{name} copied the buffer"
        assert arr.base is not None, f"{name} has no base object to keep its owner alive"

    def test_no_accessor_returns_a_writeable_array(self, noded_crossing: Any) -> None:
        # Catches a FIFTH array added later that forgets to be read-only,
        # which the parametrized tests above cannot.
        for name in dir(noded_crossing):
            if name.startswith("_"):
                continue
            value = getattr(noded_crossing, name)
            if isinstance(value, np.ndarray):
                assert value.flags.writeable is False, f"{name} is writeable"


class TestZeroCopyLifetime:
    """Risk 4 again, for the two arrays 5c adds. A view without its owner is a
    use-after-free whose symptom is a plausible-looking wrong picture, and the
    outcome is the natural thing for a caller to drop."""

    def test_arrays_outlive_the_outcome_and_the_graph(self, crossing_pslg: Any) -> None:
        outcome = _core.node(crossing_pslg, SPACING)
        assert outcome.ok(), outcome.message
        graph = outcome.pslg
        views = tuple(getattr(graph, name) for name in NODED_ARRAYS)
        expected = tuple(np.array(view) for view in views)

        del graph, outcome
        gc.collect()
        # Churn the allocator: a freed buffer that is merely still mapped
        # compares equal and proves nothing, so make reuse likely.
        churn = [np.full(1 << 16, i, dtype=np.float64) for i in range(64)]
        assert len(churn) == 64

        for view, want in zip(views, expected, strict=True):
            # This pair is the check, in this order. Reading a freed buffer is
            # undefined rather than reliably wrong, so `array_equal` alone
            # could pass on allocator luck; `base` is a live `NodedPslg` only
            # if both links are held, and that holds or fails whatever the
            # bytes say.
            assert isinstance(view.base, _core.NodedPslg)
            assert sys.getrefcount(view.base) > 1
            assert np.array_equal(view, want)


class TestNodeOutcome:
    def test_mirrors_the_build_result_and_the_cdt_outcome(self, crossing_pslg: Any) -> None:
        outcome = _core.node(crossing_pslg, SPACING)
        assert isinstance(outcome.status, type(_core.NodeStatus.Ok))
        assert isinstance(outcome.message, str)
        assert outcome.ok() is True

    def test_ok_is_a_method_not_a_property(self, crossing_pslg: Any) -> None:
        # The same trap `CdtOutcome.ok()` carries: `if outcome.ok` is truthy
        # for a bound method, so a caller who forgets the parentheses never
        # sees a failure.
        outcome = _core.node(crossing_pslg, SPACING)
        assert callable(outcome.ok)

    def test_fields_are_read_only(self, crossing_pslg: Any) -> None:
        outcome = _core.node(crossing_pslg, SPACING)
        with pytest.raises(AttributeError):
            outcome.status = _core.NodeStatus.NotConverged

    def test_is_not_constructible_from_python(self) -> None:
        with pytest.raises(TypeError):
            _core.NodeOutcome()

    def test_not_run_never_crosses_from_a_real_call(self, crossing_pslg: Any) -> None:
        # `NotRun` is a default-constructed outcome and a self-check: no call
        # through the binding may produce it.
        for spacing in (SPACING, 0.0, 1e-12):
            assert _core.node(crossing_pslg, spacing).status != _core.NodeStatus.NotRun


# --------------------------------------------------------------------------
# End to end: the criterion a person checks
# --------------------------------------------------------------------------


class TestEndToEnd:
    def test_the_crossing_triangulates_after_noding(self, crossing_pslg: Any) -> None:
        # `05b-noder-driver.md:143-146`'s acceptance criterion with the call it
        # was waiting for. Before 5c this input reached `CdtStatus.NotNoded`
        # and zero triangles; that is the picture the increment changes.
        outcome = _core.node(crossing_pslg, SPACING)
        assert outcome.status == _core.NodeStatus.Ok
        result = _core.triangulate(outcome.pslg)
        assert result.status == _core.CdtStatus.Ok, result.message
        assert result.mesh.triangle_count > 0

    def test_the_mesh_contains_the_constructed_node_once(self, noded_crossing: Any) -> None:
        meeting = np.array([350.0, 350.0]) + ORIGIN
        mesh = _core.triangulate(noded_crossing).mesh
        vertices = np.asarray(mesh.vertices)
        hits = np.flatnonzero(np.abs(vertices - meeting).max(axis=1) <= SPACING)
        assert hits.size == 1

    def test_the_mesh_vertex_array_begins_with_the_noded_graphs(
        self, noded_crossing: Any
    ) -> None:
        # `triangulate.hpp` obligation 1, restated against the type that now
        # feeds it: index k means the same point on both sides, and `cli.py`'s
        # scene join is built entirely on that.
        mesh = _core.triangulate(noded_crossing).mesh
        n = noded_crossing.vertices.shape[0]
        assert np.array_equal(np.asarray(mesh.vertices)[:n], noded_crossing.vertices)


# --------------------------------------------------------------------------
# The GIL, and the purity that licenses releasing it
# --------------------------------------------------------------------------


class _Ticker:
    """A thread that ticks once per millisecond, and cannot tick without the
    GIL. `time.sleep` releases the GIL and must re-acquire it to continue, so a
    tick is direct evidence that the GIL was available."""

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


def ticks_during[T](call: Callable[[], T]) -> tuple[T, int, float]:
    """Return the call's result, the ticks observed during it, and its duration."""
    with _Ticker() as ticker:
        time.sleep(0.02)  # let the ticker reach steady state
        before = ticker.ticks
        started = time.perf_counter()
        result = call()
        elapsed = time.perf_counter() - started
        after = ticker.ticks
    return result, after - before, elapsed


RELEASED_TICKS = 20


@pytest.fixture(scope="module")
def ladder_pslg() -> Any:
    """Enough noding work for the ticker to resolve: a grid of crossing
    breaklines, every crossing a node the driver has to construct. Deterministic
    by construction -- no random numbers -- and the counts are chosen so that
    every crossing is strictly interior and no two coincide."""
    rungs = 48
    offsets = np.linspace(20.0, 680.0, rungs)
    ring = np.array([[0.0, 0.0], [700.0, 0.0], [700.0, 700.0], [0.0, 700.0]]) + ORIGIN
    points = [ring]
    chains: list[Any] = [([0, 1, 2, 3], _core.ChainRole.Outer, NO_PROPERTIES)]
    next_id = 4
    for y in offsets:
        points.append(np.array([[10.0, y], [690.0, y]]) + ORIGIN)
        chains.append(([next_id, next_id + 1], _core.ChainRole.Breakline, RIVER))
        next_id += 2
    for x in offsets + 0.5:
        points.append(np.array([[x, 10.0], [x, 690.0]]) + ORIGIN)
        chains.append(([next_id, next_id + 1], _core.ChainRole.Breakline, ROAD))
        next_id += 2
    return build(np.vstack(points), chains)


class TestGil:
    def test_the_ticker_distinguishes_a_held_gil_from_a_released_one(self) -> None:
        """The control. Without it a pass and a probe that never ran look the
        same, and the test below would measure nothing. `sorted` on a shuffled
        list holds the GIL for its whole duration; `time.sleep` releases it."""
        import random

        data = list(range(3_000_000))
        random.Random(12345).shuffle(data)

        _, held, held_elapsed = ticks_during(lambda: sorted(data))
        assert held_elapsed >= 0.05, f"control workload too fast ({held_elapsed:.3f}s)"
        assert held <= 3, f"a GIL-holding call let {held} ticks through"

        _, released, _ = ticks_during(lambda: time.sleep(0.3))
        assert released >= RELEASED_TICKS, f"a releasing call let only {released} ticks through"

    def test_node_releases_the_gil(self, ladder_pslg: Any) -> None:
        outcome, ticks, elapsed = ticks_during(lambda: _core.node(ladder_pslg, SPACING))
        assert outcome.ok(), outcome.message
        # A probe able to fail: if the workload is too fast there is nothing to
        # observe, and this must say so rather than pass on an empty window.
        assert elapsed >= 0.05, (
            f"noding took only {elapsed:.3f}s -- too fast to measure GIL release; "
            "enlarge ladder_pslg rather than lowering this bound"
        )
        assert ticks >= RELEASED_TICKS, (
            f"only {ticks} ticks in {elapsed:.3f}s: node appears to hold the GIL"
        )

    def test_concurrent_noding_of_one_pslg_is_bit_identical(self, crossing_pslg: Any) -> None:
        # The property `node.hpp:5-11` claims, and the ONLY thing that makes
        # releasing the GIL safe: `node<K>` is a pure function of (pslg,
        # options), so two threads on the same input produce identical output.
        # With the GIL released this stops being theoretical.
        results: list[Any] = []
        barrier = threading.Barrier(4)

        def run() -> None:
            barrier.wait(timeout=30.0)
            out = _core.node(crossing_pslg, SPACING)
            results.append((np.array(out.pslg.vertices), np.array(out.pslg.chain_indices)))

        threads = [threading.Thread(target=run) for _ in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30.0)
            assert not t.is_alive()
        assert len(results) == 4
        for other in results[1:]:
            assert np.array_equal(results[0][0], other[0])
            assert np.array_equal(results[0][1], other[1])


# --------------------------------------------------------------------------
# The firewall: what must never cross
# --------------------------------------------------------------------------


class TestBoundary:
    """The binding surface is exhaustive and an addition needs a reason in an
    increment file. `05c-noder-wiring.md`'s "what does not cross" table is the
    reason each of these is absent; asserting absence is what stops the
    boundary eroding by accretion."""

    PROHIBITED = (
        "SnapGrid",
        "GridPoint",
        "NodeOptions",
        "NodedPslgBuilder",
        "NodeSet",
        "EdgeProperties",
    )

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

    @pytest.mark.parametrize("obj_name", ["NodedPslg", "NodeOutcome"])
    def test_no_new_type_exposes_io_crs_or_topology(self, obj_name: str) -> None:
        cls = getattr(_core, obj_name)
        offenders = [
            name
            for name in dir(cls)
            if not name.startswith("_")
            for bad in self.FORBIDDEN_SUBSTRINGS
            if bad in name.lower()
        ]
        assert offenders == []

    def test_no_composite_entry_point_is_offered(self) -> None:
        # Composition is `cli.py`'s, per the blueprint's boundary 4. A
        # convenience `mesh(vertices, chains, spacing)` here would be a second
        # composition root, and the first thing it would swallow is which of
        # the three refusals the caller got.
        for name in ("mesh", "build_mesh", "make_mesh", "pipeline"):
            assert not hasattr(_core, name)


# --------------------------------------------------------------------------
# Documentation and stubs
# --------------------------------------------------------------------------


class TestDocumentation:
    @pytest.mark.parametrize("name", NEW_SURFACE)
    def test_every_bound_name_carries_a_docstring(self, name: str) -> None:
        assert getattr(_core, name).__doc__

    def test_describes_docstring_no_longer_claims_one_enumeration(self) -> None:
        # The shipped docstring says "Raises TypeError for anything that is not
        # a CdtStatus". That sentence becomes false on the day the second
        # overload lands, and a docstring that is false about the type it
        # rejects is worse than none.
        text = _core.describe.__doc__ or ""
        assert "NodeStatus" in text


class TestStubs:
    """`_core.pyi` is hand-written, and mypy checks its consistency rather than
    its completeness: a name absent from both the stub and the caller passes
    strict mode in silence. `cli.py` calls `node()`, so a missing stub is a
    mypy failure one layer up -- and a *wrong* stub is not."""

    def tree(self) -> ast.Module:
        return ast.parse(STUB.read_text(encoding="utf-8"))

    @pytest.mark.parametrize("name", NEW_SURFACE)
    def test_stub_declares_the_new_surface(self, name: str) -> None:
        declared = {
            node.name
            for node in self.tree().body
            if isinstance(node, ast.ClassDef | ast.FunctionDef)
        }
        assert name in declared

    def test_describe_is_declared_as_an_overload_set(self) -> None:
        # Two `@overload` stubs plus an implementation stub, the pattern
        # `cross` already uses in this file. One stub naming a union would
        # type-check a call that pybind11 resolves by exact type.
        overloads = [
            node
            for node in self.tree().body
            if isinstance(node, ast.FunctionDef)
            and node.name == "describe"
            and any(
                getattr(d, "id", getattr(d, "attr", None)) == "overload"
                for d in node.decorator_list
            )
        ]
        assert len(overloads) == 2

    def test_the_noded_pslg_stub_carries_the_three_extra_accessors(self) -> None:
        cls = next(
            node
            for node in ast.walk(self.tree())
            if isinstance(node, ast.ClassDef) and node.name == "NodedPslg"
        )
        members = {node.name for node in cls.body if isinstance(node, ast.FunctionDef)}
        assert {"grid_spacing", "edge_properties", "node_of_input_vertex"} <= members
        assert "edge_base" not in members

    def test_the_outcomes_pslg_is_optional_in_the_stub(self) -> None:
        # The `Optional` is the type system carrying "engaged iff status ==
        # Ok", and it is what makes `cli.py`'s `if outcome.pslg is None` a
        # narrowing rather than a defensive check mypy cannot see through.
        cls = next(
            node
            for node in ast.walk(self.tree())
            if isinstance(node, ast.ClassDef) and node.name == "NodeOutcome"
        )
        pslg = next(
            node
            for node in cls.body
            if isinstance(node, ast.FunctionDef) and node.name == "pslg"
        )
        assert ast.unparse(pslg.returns or ast.Constant(None)) == "NodedPslg | None"

    def test_triangulate_is_declared_over_the_noded_type(self) -> None:
        # The architectural product of 5b+5c: un-noded input is unrepresentable
        # at the entry point rather than diagnosed inside it. A stub still
        # saying `Pslg` would let `cli.py` pass one and fail at run time.
        fn = next(
            node
            for node in self.tree().body
            if isinstance(node, ast.FunctionDef) and node.name == "triangulate"
        )
        assert ast.unparse(fn.args.args[0].annotation or ast.Constant(None)) == "NodedPslg"
