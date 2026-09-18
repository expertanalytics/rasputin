"""`tin_engine.viz.scene`: increment 6b-i, and increment 6's only
invariant-critical suite.

Committed RED, before `viz/scene.py` exists. Every test imports the module
inside its body, so the intended failure is one named `ModuleNotFoundError` per
test rather than a collection error that takes the rest of `tests/python/` down
with it.

`06-cdt-viewer.md` names this file as the sole suite in increment 6 that carries
a mutation round, and names the mutant: **a renderer that reads the
constrained-edge mask under CGAL's opposite-vertex convention instead of
`indexed_mesh.hpp`'s.** That mutant draws a completely plausible picture with
the constraint strokes on the wrong edges, and nothing downstream looks wrong.
`TestMaskConvention` exists for it, and `MESH_TRIANGLES`/`MESH_MASKS` are chosen
so that the rotated read produces a constrained set **disjoint** from the
correct one -- see `EXPECTED_CONSTRAINED` and `ROTATED_CONSTRAINED`.

Three things about the oracle, per `.claude/REQUIRED-READING.md`'s
object-identity rule:

* The masks are hand-written integers and the expected edge sets are
  hand-written pairs. Neither is computed with the `(v[e], v[(e+1)%3])`
  arithmetic the module under test uses, because an oracle that borrows the
  producer's formula agrees with the producer's bug.
* The fakes are plain dataclasses, not `_core` objects. The suite runs with no
  compiled extension in the process, which is the whole reason `protocols.py`
  exists; `TestModuleIsolation` asserts the module keeps it that way.
* `TestFixtureSanity` runs first and is able to fail on its own: if the fixture
  ever stops containing a triangle with exactly one constrained bit, the
  mutation round silently stops measuring anything, because the two conventions
  agree on masks of 0 and 0b111.

Coordinates are UTM 33N-shaped (easting ~4.3e5, northing ~6.9e6) per the
design's ordering ruling: a picture drawn in the unit square silently exercises
the easy case.

## The API this suite pins, and where it departs from the design text

`build_scene(pslg, mesh=None, *, ok=True, closed_roles=())`.

* **`ok` is a parameter, not a status object.** `viz/` may not name
  `_core.CdtStatus`, and the status string and the backend message are pure
  passthrough to the header band, which is `svg.py`'s. The scene carries the
  classification (`SceneKind`) and nothing it does not use.
* **`SceneEdge.constrained` and `SceneEdge.role` are separate.** `constrained`
  is the mesh mask's verdict; `role` is the chain join's. Collapsing them loses
  the two disagreement states the design requires the renderer to draw, and it
  makes the mask-convention oracle vacuous -- see `masked_pairs`.
* **`closed_roles` is supplied by the caller.** The design has the scene close
  `Outer` and `Hole` chains -- a `Pslg` never stores a ring's closure -- but
  `ChainLike.role` is typed `object` precisely so `viz/` cannot name
  `ChainRole`, so the scene cannot tell a ring from a breakline on its own.
  Membership needs only `==`/`__hash__`, which `object` has, so the policy moves
  to the composition root that already knows the enum. Without it, every ring's
  closing edge is a false `MASKED_EDGE_WITHOUT_CHAIN` finding --
  `test_closure_policy_is_the_callers` pins exactly that.
"""

from __future__ import annotations

import ast
import enum
import importlib
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCENE_PATH = REPO_ROOT / "src_python" / "tin_engine" / "viz" / "scene.py"

EAST = 430_000.0
NORTH = 6_900_000.0


class Role(enum.Enum):
    """A stand-in for `_core.ChainRole`, which this suite may not import.

    The scene treats a role as an opaque token: it compares roles and uses them
    as mapping keys and never asks what they mean. A local enum is therefore a
    faithful fake, and using one keeps the extension out of the process.
    """

    Outer = "Outer"
    Hole = "Hole"
    Breakline = "Breakline"


CLOSED = (Role.Outer, Role.Hole)


@dataclass(frozen=True)
class FakeChain:
    begin: int
    count: int
    role: Role
    is_river: bool


@dataclass(frozen=True)
class FakePslg:
    vertices: npt.NDArray[np.float64]
    chains: tuple[FakeChain, ...]
    chain_indices: npt.NDArray[np.uint32]

    def indices_of(self, c: int) -> npt.NDArray[np.uint32]:
        chain = self.chains[c]
        return self.chain_indices[chain.begin : chain.begin + chain.count]


@dataclass(frozen=True)
class FakeMesh:
    vertices: npt.NDArray[np.float64]
    triangles: npt.NDArray[np.uint32]
    constrained_edges: npt.NDArray[np.uint8]

    @property
    def triangle_count(self) -> int:
        return int(self.triangles.shape[0])

    @property
    def empty(self) -> bool:
        return self.triangle_count == 0


# --- The fixture, drawn on paper first -------------------------------------
#
# A 100 x 80 m rectangle at UTM 33N magnitudes, with a horizontal breakline
# across its middle. Six vertices, six triangles, eleven undirected edges
# (Euler: (3*6 + 4)/2 = 11), five of them constrained.
#
#   3 +-------------------+ 2
#     | \               / |
#     |   4 ==========5   |        4=5 is the breakline (is_river)
#     | /               \ |
#   0 +-------------------+ 1
#
MESH_VERTICES = np.array(
    [
        [EAST + 0.0, NORTH + 0.0],
        [EAST + 100.0, NORTH + 0.0],
        [EAST + 100.0, NORTH + 80.0],
        [EAST + 0.0, NORTH + 80.0],
        [EAST + 20.0, NORTH + 40.0],
        [EAST + 80.0, NORTH + 40.0],
    ],
    dtype=np.float64,
)

# `indexed_mesh.hpp` guarantee 4: the mesh's vertex array *begins with* the
# PSLG's, it does not equal it -- the backend may append points of its own.
# Index 6 stands for one such point. It is placed outside the rectangle so that
# a scene built from the input vertices instead of the mesh's has a bbox that
# visibly omits it, and no triangle indexes it, so `MESH_TRIANGLES` and
# `MESH_MASKS` stay valid against this vertex array unchanged.
BACKEND_VERTEX = np.array([[EAST + 160.0, NORTH + 130.0]], dtype=np.float64)
MESH_VERTICES_PLUS_BACKEND = np.concatenate([MESH_VERTICES, BACKEND_VERTEX])

# Every triangle counterclockwise, which is what `triangulate` guarantees; the
# six of them tile the rectangle exactly (8000 m^2, checked by hand).
MESH_TRIANGLES = np.array(
    [
        [0, 1, 5],
        [0, 5, 4],
        [0, 4, 3],
        [1, 2, 5],
        [4, 5, 2],
        [4, 2, 3],
    ],
    dtype=np.uint32,
)

# Hand-written, NOT derived from the convention under test. Bit e is edge
# (v[e], v[(e+1)%3]):
#   (0,1,5): (0,1) constrained          -> 0b001
#   (0,5,4): (5,4) constrained          -> 0b010
#   (0,4,3): (3,0) constrained          -> 0b100
#   (1,2,5): (1,2) constrained          -> 0b001
#   (4,5,2): (4,5) constrained          -> 0b001
#   (4,2,3): (2,3) constrained          -> 0b010
# Every mask has exactly one bit set, which is the case that distinguishes the
# two conventions most sharply.
MESH_MASKS = np.array([0b001, 0b010, 0b100, 0b001, 0b001, 0b010], dtype=np.uint8)

# The four ring edges plus the breakline, written out as sorted index pairs.
EXPECTED_CONSTRAINED = frozenset({(0, 1), (1, 2), (2, 3), (0, 3), (4, 5)})

# What the CGAL-convention mutant would produce: bit e read as the edge
# opposite vertex e, (v[(e+1)%3], v[(e+2)%3]). Hand-derived per triangle:
# (1,5), (4,0), (0,4), (2,5), (5,2), (3,4). Every one of them is an interior
# edge, so the mutant's set is disjoint from the correct one.
ROTATED_CONSTRAINED = frozenset({(1, 5), (0, 4), (2, 5), (3, 4)})

ALL_EDGES = frozenset(
    {
        (0, 1),
        (0, 3),
        (0, 4),
        (0, 5),
        (1, 2),
        (1, 5),
        (2, 3),
        (2, 4),
        (2, 5),
        (3, 4),
        (4, 5),
    }
)

MESH_MEMBERS = ("vertices", "triangles", "constrained_edges", "triangle_count", "empty")
PSLG_MEMBERS = ("vertices", "chains", "chain_indices", "indices_of")
CHAIN_MEMBERS = ("begin", "count", "role", "is_river")


def scene_module() -> ModuleType:
    return importlib.import_module("tin_engine.viz.scene")


def pairs(edges: Any) -> frozenset[tuple[int, int]]:
    return frozenset((int(e.a), int(e.b)) for e in edges)


def masked_pairs(scene: Any) -> frozenset[tuple[int, int]]:
    """Edges the *mesh mask* called constrained. This is the object the mask
    convention is a claim about, and it is deliberately not `role is not None`:
    a masked edge that matches no chain has no role, so a role-based oracle goes
    vacuously green under exactly the rotation mutant this suite exists for.
    Measured -- with the role-based helper, two of the four `TestMaskConvention`
    assertions passed against the rotated mutant.
    """
    return pairs(e for e in scene.edges if e.constrained)


def roled_pairs(scene: Any) -> frozenset[tuple[int, int]]:
    """Edges the *input PSLG* called constraints, via the role join."""
    return pairs(e for e in scene.edges if e.role is not None)


def edge_at(scene: Any, a: int, b: int) -> Any:
    matches = [e for e in scene.edges if (int(e.a), int(e.b)) == (a, b)]
    assert len(matches) == 1, f"edge {(a, b)} appears {len(matches)} times"
    return matches[0]


def make_pslg(
    chains: list[tuple[list[int], Role, bool]],
    vertices: npt.NDArray[np.float64] | None = None,
) -> FakePslg:
    flat: list[int] = []
    records: list[FakeChain] = []
    for indices, role, is_river in chains:
        records.append(FakeChain(len(flat), len(indices), role, is_river))
        flat.extend(indices)
    return FakePslg(
        vertices=MESH_VERTICES if vertices is None else vertices,
        chains=tuple(records),
        chain_indices=np.array(flat, dtype=np.uint32),
    )


@pytest.fixture
def pslg() -> FakePslg:
    """The input side: a closed outer ring and an open river breakline.

    The ring does not store its closure -- `testing.md`'s PSLG catalog pins
    that -- so the (3, 0) edge exists only if the scene closes the chain.
    """
    return make_pslg(
        [
            ([0, 1, 2, 3], Role.Outer, False),
            ([4, 5], Role.Breakline, True),
        ]
    )


@pytest.fixture
def mesh() -> FakeMesh:
    return FakeMesh(MESH_VERTICES, MESH_TRIANGLES, MESH_MASKS)


@pytest.fixture
def mesh_with_backend_vertex() -> FakeMesh:
    """The same triangulation over a vertex array longer than the input's."""
    return FakeMesh(MESH_VERTICES_PLUS_BACKEND, MESH_TRIANGLES, MESH_MASKS)


@pytest.fixture
def empty_mesh() -> FakeMesh:
    """`Ok` with zero triangles: increment 4's risk 5, the silent mode."""
    return FakeMesh(
        MESH_VERTICES,
        np.zeros((0, 3), dtype=np.uint32),
        np.zeros((0,), dtype=np.uint8),
    )


@pytest.fixture
def scene(pslg: FakePslg, mesh: FakeMesh) -> Any:
    return scene_module().build_scene(pslg, mesh, closed_roles=CLOSED)


class TestFixtureSanity:
    """Probes that can fail on their own.

    If any of these stops holding, the mutation round below keeps passing while
    measuring nothing -- the two mask conventions agree on masks of 0 and 0b111,
    so a fixture that drifted to those would make `TestMaskConvention` vacuous.
    """

    @pytest.mark.parametrize(
        ("fake", "members"),
        [
            (FakeMesh(MESH_VERTICES, MESH_TRIANGLES, MESH_MASKS), MESH_MEMBERS),
            (make_pslg([([0, 1, 2, 3], Role.Outer, False)]), PSLG_MEMBERS),
            (FakeChain(0, 4, Role.Outer, False), CHAIN_MEMBERS),
        ],
    )
    def test_the_fakes_carry_every_protocol_member(
        self, fake: Any, members: tuple[str, ...]
    ) -> None:
        # Joined against the real bound types by test_viz_protocols.py; here it
        # only has to catch a typo in this file, so that a red test means the
        # scene is wrong and not that the fixture is.
        missing = [m for m in members if not hasattr(fake, m)]
        assert missing == []

    def test_every_triangle_has_exactly_one_constrained_bit(self) -> None:
        counts = [int(m).bit_count() for m in MESH_MASKS]
        assert counts == [1] * len(MESH_TRIANGLES)

    def test_the_rotated_convention_is_distinguishable(self) -> None:
        # The kill condition, stated as data. Disjoint is stronger than
        # unequal: under the mutant not one constraint stroke lands on a
        # correct edge.
        assert EXPECTED_CONSTRAINED.isdisjoint(ROTATED_CONSTRAINED)

    def test_the_rotated_convention_only_names_interior_edges(self) -> None:
        assert ROTATED_CONSTRAINED < ALL_EDGES

    def test_coordinates_are_at_utm_magnitudes(self) -> None:
        assert MESH_VERTICES[:, 0].min() >= 1.0e5
        assert MESH_VERTICES[:, 1].min() >= 1.0e6

    def test_indices_of_partitions_the_flat_buffer(self, pslg: FakePslg) -> None:
        recovered = np.concatenate([pslg.indices_of(c) for c in range(len(pslg.chains))])
        assert np.array_equal(recovered, pslg.chain_indices)


class TestModuleIsolation:
    """`scene.py` imports the protocols and nothing else first-party.

    That is what lets this entire suite run against hand-built dataclasses with
    no compiled extension in the process, and it is checked by reading the
    source: `tin_engine/__init__.py` imports `_core` itself, so a `sys.modules`
    assertion would pass for the wrong reason.
    """

    def source(self) -> ast.Module:
        assert SCENE_PATH.is_file(), f"{SCENE_PATH} does not exist"
        return ast.parse(SCENE_PATH.read_text(encoding="utf-8"))

    def imports(self) -> list[str]:
        names: list[str] = []
        for node in ast.walk(self.source()):
            if isinstance(node, ast.Import):
                names.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                names.append("." * node.level + (node.module or ""))
        return names

    def test_the_module_exists(self) -> None:
        assert SCENE_PATH.is_file()

    def test_it_does_not_import_the_extension(self) -> None:
        assert [name for name in self.imports() if "_core" in name] == []

    def test_it_imports_the_protocols(self) -> None:
        assert [name for name in self.imports() if name.endswith("protocols")] != []

    def test_it_imports_no_other_first_party_module(self) -> None:
        first_party = [
            name
            for name in self.imports()
            if name.startswith((".", "tin_engine")) and not name.endswith("protocols")
        ]
        assert first_party == []


class TestMaskConvention:
    """The mutation round.

    `indexed_mesh.hpp` defines bit `e` as edge `(v[e], v[(e+1)%3])` and says in
    so many words that this is not CGAL's opposite-vertex convention. A scene
    that rotates it is a wrong picture that does not look wrong.
    """

    def test_constrained_edges_are_the_ring_and_the_breakline(self, scene: Any) -> None:
        assert masked_pairs(scene) == EXPECTED_CONSTRAINED

    def test_the_rotated_read_is_not_what_the_scene_produced(self, scene: Any) -> None:
        # Redundant with the line above by construction, and kept because it
        # names the mutant at the point of failure rather than in a docstring.
        assert masked_pairs(scene).isdisjoint(ROTATED_CONSTRAINED)

    @pytest.mark.parametrize(
        ("triangle", "expected"),
        [(0, (0, 1)), (1, (4, 5)), (2, (0, 3)), (3, (1, 2)), (4, (4, 5)), (5, (2, 3))],
    )
    def test_each_single_bit_mask_names_its_own_edge(
        self, pslg: FakePslg, triangle: int, expected: tuple[int, int]
    ) -> None:
        # Per-triangle attribution, which is what the name claims: membership in
        # the whole scene's constrained set is already implied by the set
        # equality above and would pass for any triangle. The mesh is cut down
        # to this one triangle, so the masked set can only come from its bit.
        one = FakeMesh(
            MESH_VERTICES,
            MESH_TRIANGLES[triangle : triangle + 1],
            MESH_MASKS[triangle : triangle + 1],
        )
        scene = scene_module().build_scene(pslg, one, closed_roles=CLOSED)
        assert masked_pairs(scene) == frozenset({expected})

    def test_no_interior_edge_is_constrained(self, scene: Any) -> None:
        interior = ALL_EDGES - EXPECTED_CONSTRAINED
        assert masked_pairs(scene).isdisjoint(interior)


class TestEdgeDedup:
    def test_every_undirected_edge_appears_once(self, scene: Any) -> None:
        listed = [(int(e.a), int(e.b)) for e in scene.edges]
        assert sorted(listed) == sorted(ALL_EDGES)

    def test_an_interior_edge_shared_by_two_triangles_is_emitted_once(self, scene: Any) -> None:
        # (0, 5) is triangle 0's third edge and triangle 1's first. Drawn twice
        # it is drawn at double weight and the picture lies about density.
        assert edge_at(scene, 0, 5).role is None

    def test_a_constrained_edge_shared_by_two_triangles_is_not_double_classified(
        self, scene: Any
    ) -> None:
        # The breakline (4, 5) is bit 1 of triangle 1 and bit 0 of triangle 4.
        # `constrained` is the mask's verdict and so is the object
        # "double-classified" is a claim about; `role` is the chain join's, and
        # is asserted alongside only because the two must agree here.
        breakline = edge_at(scene, 4, 5)
        assert breakline.constrained is True
        assert breakline.role is Role.Breakline

    def test_endpoints_are_canonically_ordered(self, scene: Any) -> None:
        assert all(int(e.a) < int(e.b) for e in scene.edges)

    def test_edge_order_is_deterministic(self, pslg: FakePslg, mesh: FakeMesh) -> None:
        build = scene_module().build_scene
        first = build(pslg, mesh, closed_roles=CLOSED)
        second = build(pslg, mesh, closed_roles=CLOSED)
        assert [(e.a, e.b) for e in first.edges] == [(e.a, e.b) for e in second.edges]
        assert [(e.a, e.b) for e in first.edges] == sorted((e.a, e.b) for e in first.edges)

    def test_triangles_survive_unchanged(self, scene: Any) -> None:
        assert np.array_equal(np.asarray(scene.triangles), MESH_TRIANGLES)


class TestRoleJoin:
    """The constrained-edge mask and the PSLG's chain roles are two independent
    derivations of "this edge is a constraint". The join is sound because the
    mesh vertex array begins with the PSLG's and detria splits nothing; where
    the two disagree the scene **counts and surfaces** the disagreement rather
    than reconciling it (risk 2).
    """

    def test_ring_edges_carry_the_outer_role(self, scene: Any) -> None:
        for a, b in [(0, 1), (1, 2), (2, 3), (0, 3)]:
            assert edge_at(scene, a, b).role is Role.Outer

    def test_the_river_bit_reaches_the_edge(self, scene: Any) -> None:
        assert edge_at(scene, 4, 5).is_river is True

    def test_an_ordinary_constrained_edge_is_not_a_river(self, scene: Any) -> None:
        assert edge_at(scene, 0, 1).is_river is False

    def test_unconstrained_edges_have_no_role_and_no_river_bit(self, scene: Any) -> None:
        for a, b in sorted(ALL_EDGES - EXPECTED_CONSTRAINED):
            edge = edge_at(scene, a, b)
            assert edge.constrained is False
            assert edge.role is None
            assert edge.is_river is False

    def test_the_join_agrees_on_well_formed_input(self, scene: Any) -> None:
        assert list(scene.findings) == []

    def test_an_edge_in_two_chains_takes_the_first_chains_role(self, mesh: FakeMesh) -> None:
        # Stated precedence: the earliest chain containing the pair wins the
        # role, because `role` is opaque and cannot be ranked by value.
        pslg = make_pslg(
            [
                ([0, 1, 2, 3], Role.Outer, False),
                ([4, 5], Role.Breakline, True),
                ([0, 1], Role.Breakline, True),
            ]
        )
        scene = scene_module().build_scene(pslg, mesh, closed_roles=CLOSED)
        assert edge_at(scene, 0, 1).role is Role.Outer

    @pytest.mark.parametrize("river_first", [True, False])
    def test_the_river_bit_is_the_or_over_every_chain_on_the_edge(
        self, mesh: FakeMesh, river_first: bool
    ) -> None:
        # Mirrors `testing.md`'s noding invariant: the bit means "some
        # contributing chain was a river", so neither an earlier nor a later
        # chain may clear it. Both orders, because a single order is passed by
        # an implementation that simply takes the last chain's bit -- measured:
        # with only the river-last case, the overwrite mutant survived.
        bits = [True, False] if river_first else [False, True]
        pslg = make_pslg(
            [
                ([0, 1, 2, 3], Role.Outer, False),
                ([4, 5], Role.Breakline, True),
                ([0, 1], Role.Breakline, bits[0]),
                ([0, 1], Role.Breakline, bits[1]),
            ]
        )
        scene = scene_module().build_scene(pslg, mesh, closed_roles=CLOSED)
        assert edge_at(scene, 0, 1).is_river is True

    def test_a_masked_edge_matching_no_chain_is_a_finding(self, pslg: FakePslg) -> None:
        # Triangle 0's mask gains bit 2, edge (5, 0) -- an interior edge no
        # chain contains. This is the backend-mask half of the disagreement.
        masks = MESH_MASKS.copy()
        masks[0] = 0b101
        scene = scene_module().build_scene(
            pslg, FakeMesh(MESH_VERTICES, MESH_TRIANGLES, masks), closed_roles=CLOSED
        )
        kinds = scene_module().FindingKind
        assert [(f.kind, f.a, f.b) for f in scene.findings] == [
            (kinds.MASKED_EDGE_WITHOUT_CHAIN, 0, 5)
        ]

    def test_a_masked_edge_matching_no_chain_is_still_drawn(self, pslg: FakePslg) -> None:
        masks = MESH_MASKS.copy()
        masks[0] = 0b101
        scene = scene_module().build_scene(
            pslg, FakeMesh(MESH_VERTICES, MESH_TRIANGLES, masks), closed_roles=CLOSED
        )
        assert pairs(scene.edges) == ALL_EDGES
        # The mask verdict survives the failed join: alarm state, not silence.
        orphan = edge_at(scene, 0, 5)
        assert orphan.constrained is True
        assert orphan.role is None

    def test_a_chain_edge_in_no_mask_is_a_finding(self, mesh: FakeMesh) -> None:
        # (0, 2) is the rectangle's diagonal: a chain edge that is not an edge
        # of any triangle. This is the input half of the disagreement.
        pslg = make_pslg(
            [
                ([0, 1, 2, 3], Role.Outer, False),
                ([4, 5], Role.Breakline, True),
                ([0, 2], Role.Breakline, False),
            ]
        )
        scene = scene_module().build_scene(pslg, mesh, closed_roles=CLOSED)
        kinds = scene_module().FindingKind
        assert [(f.kind, f.a, f.b) for f in scene.findings] == [
            (kinds.CHAIN_EDGE_WITHOUT_MASK, 0, 2)
        ]

    def test_a_chain_edge_in_no_mask_is_still_drawn(self, mesh: FakeMesh) -> None:
        pslg = make_pslg(
            [
                ([0, 1, 2, 3], Role.Outer, False),
                ([4, 5], Role.Breakline, True),
                ([0, 2], Role.Breakline, False),
            ]
        )
        scene = scene_module().build_scene(pslg, mesh, closed_roles=CLOSED)
        drawn = edge_at(scene, 0, 2)
        assert drawn.role is Role.Breakline
        assert drawn.constrained is False

    def test_closure_policy_is_the_callers(self, pslg: FakePslg, mesh: FakeMesh) -> None:
        # With no `closed_roles`, the outer ring's (3, 0) closing edge belongs
        # to no chain and the disagreement is reported rather than guessed at.
        scene = scene_module().build_scene(pslg, mesh)
        kinds = scene_module().FindingKind
        assert [(f.kind, f.a, f.b) for f in scene.findings] == [
            (kinds.MASKED_EDGE_WITHOUT_CHAIN, 0, 3)
        ]

    def test_a_breakline_is_never_closed(self, mesh: FakeMesh) -> None:
        # What excludes a breakline is `Breakline not in closed_roles`, not the
        # chain's length. The fixture's two-vertex breakline cannot show that:
        # closing it re-emits `_key(4, 5)`, and `_chain_edges` accumulates into
        # a dict keyed on exactly that pair, so the second write only ORs
        # `is_river` and the edge still appears once under the mutant. A
        # three-vertex breakline can: closing [4, 5, 2] invents the constraint
        # (2, 4), which is a real mesh edge and so would silently be drawn as
        # one.
        pslg = make_pslg(
            [
                ([0, 1, 2, 3], Role.Outer, False),
                ([4, 5, 2], Role.Breakline, True),
            ]
        )
        scene = scene_module().build_scene(pslg, mesh, closed_roles=CLOSED)
        assert edge_at(scene, 2, 4).role is None
        assert len([e for e in scene.edges if (int(e.a), int(e.b)) == (4, 5)]) == 1

    def test_a_single_vertex_ring_is_never_closed_into_a_self_loop(
        self, mesh: FakeMesh
    ) -> None:
        # This is what the length guard defends, and the only thing it does:
        # a one-vertex chain closed emits `_key(a, a)`, a self-loop that
        # violates `SceneEdge`'s `a < b` invariant and that no viewport
        # transform can draw. A two-vertex chain closed is harmless -- see
        # above -- so the guard may not be relaxed below two.
        pslg = make_pslg(
            [
                ([0, 1, 2, 3], Role.Outer, False),
                ([4, 5], Role.Breakline, True),
                ([2], Role.Outer, False),
            ]
        )
        scene = scene_module().build_scene(pslg, mesh, closed_roles=CLOSED)
        assert [(int(e.a), int(e.b)) for e in scene.edges if int(e.a) >= int(e.b)] == []
        assert pairs(scene.edges) == ALL_EDGES


class TestClassification:
    """Increment 4's risk 5: `Ok` with zero triangles is a *successful* backend
    call. 6b-ii renders it differently from a failure, so the scene has to tell
    them apart rather than both reducing to "nothing to draw".
    """

    def test_a_populated_mesh_is_a_mesh_scene(self, scene: Any) -> None:
        assert scene.kind is scene_module().SceneKind.MESH

    def test_ok_with_zero_triangles_is_its_own_kind(
        self, pslg: FakePslg, empty_mesh: FakeMesh
    ) -> None:
        scene = scene_module().build_scene(pslg, empty_mesh, ok=True, closed_roles=CLOSED)
        assert scene.kind is scene_module().SceneKind.OK_BUT_EMPTY

    def test_ok_with_zero_triangles_is_not_a_failure(
        self, pslg: FakePslg, empty_mesh: FakeMesh
    ) -> None:
        scene = scene_module().build_scene(pslg, empty_mesh, ok=True, closed_roles=CLOSED)
        assert scene.kind is not scene_module().SceneKind.FAILED

    def test_a_missing_mesh_with_an_ok_status_is_the_same_kind(self, pslg: FakePslg) -> None:
        scene = scene_module().build_scene(pslg, None, ok=True, closed_roles=CLOSED)
        assert scene.kind is scene_module().SceneKind.OK_BUT_EMPTY

    def test_a_non_ok_status_is_a_failure_even_with_a_full_mesh(
        self, pslg: FakePslg, mesh: FakeMesh
    ) -> None:
        # The bindings guarantee no mesh accompanies a non-`Ok` status, so this
        # combination is a contradiction; the status wins, and no triangle from
        # a failed run is ever drawn.
        scene = scene_module().build_scene(pslg, mesh, ok=False, closed_roles=CLOSED)
        assert scene.kind is scene_module().SceneKind.FAILED

    @pytest.mark.parametrize("ok", [True, False])
    def test_a_scene_without_a_drawable_mesh_carries_no_triangles(
        self, pslg: FakePslg, mesh: FakeMesh, empty_mesh: FakeMesh, ok: bool
    ) -> None:
        supplied = empty_mesh if ok else mesh
        scene = scene_module().build_scene(pslg, supplied, ok=ok, closed_roles=CLOSED)
        assert np.asarray(scene.triangles).shape == (0, 3)

    @pytest.mark.parametrize("ok", [True, False])
    def test_a_scene_without_a_mesh_draws_the_input_pslg(
        self, pslg: FakePslg, mesh: FakeMesh, empty_mesh: FakeMesh, ok: bool
    ) -> None:
        # Never a blank page: the user sees exactly the geometry they submitted.
        supplied = empty_mesh if ok else mesh
        scene = scene_module().build_scene(pslg, supplied, ok=ok, closed_roles=CLOSED)
        assert roled_pairs(scene) == EXPECTED_CONSTRAINED
        assert pairs(scene.edges) == EXPECTED_CONSTRAINED

    @pytest.mark.parametrize("ok", [True, False])
    def test_a_scene_without_a_mesh_reports_no_findings(
        self, pslg: FakePslg, mesh: FakeMesh, empty_mesh: FakeMesh, ok: bool
    ) -> None:
        # There is no mask to disagree with, so every chain edge would be a
        # CHAIN_EDGE_WITHOUT_MASK finding -- five alarms on a picture whose real
        # alarm is the status. The join itself still runs and the roles it
        # attaches survive (`test_the_pslg_only_scene_keeps_the_roles`); what is
        # suppressed is the findings derived from comparing it to a mask that
        # was never produced.
        supplied = empty_mesh if ok else mesh
        scene = scene_module().build_scene(pslg, supplied, ok=ok, closed_roles=CLOSED)
        assert list(scene.findings) == []

    def test_the_pslg_only_scene_keeps_the_roles(self, pslg: FakePslg) -> None:
        scene = scene_module().build_scene(pslg, None, ok=False, closed_roles=CLOSED)
        assert edge_at(scene, 4, 5).role is Role.Breakline
        assert edge_at(scene, 0, 1).role is Role.Outer
        # No mesh, so no mask: nothing may claim a verdict that was never made.
        assert [e for e in scene.edges if e.constrained] == []


class TestBoundingBox:
    """A zero-extent bbox is a division by zero in 6b-ii's viewport transform,
    and a picture that is an accidental point.

    Every comparison here is absolute. `pytest.approx`'s default is relative,
    and at an easting of 4.3e5 that is a tolerance of +-0.43 m -- wide enough to
    pass a padding that puts the box off-centre by a quarter of a metre.
    Measured: with the default tolerance the off-centre mutant survived the
    whole suite.
    """

    def test_it_contains_every_vertex(self, scene: Any) -> None:
        box = scene.bbox
        assert box.min_x == pytest.approx(MESH_VERTICES[:, 0].min(), rel=0.0, abs=1e-9)
        assert box.max_x == pytest.approx(MESH_VERTICES[:, 0].max(), rel=0.0, abs=1e-9)
        assert box.min_y == pytest.approx(MESH_VERTICES[:, 1].min(), rel=0.0, abs=1e-9)
        assert box.max_y == pytest.approx(MESH_VERTICES[:, 1].max(), rel=0.0, abs=1e-9)

    def test_a_healthy_box_is_not_padded(self, scene: Any) -> None:
        assert scene.bbox.padded is False

    def collinear(self, points: list[tuple[float, float]]) -> Any:
        vertices = np.array(points, dtype=np.float64)
        indices = list(range(len(points)))
        pslg = make_pslg([(indices, Role.Breakline, False)], vertices=vertices)
        return scene_module().build_scene(pslg, None, ok=False)

    def test_a_zero_width_box_is_padded_about_its_centre(self) -> None:
        column = [(EAST, NORTH), (EAST, NORTH + 40.0), (EAST, NORTH + 80.0)]
        box = self.collinear(column).bbox
        assert box.padded is True
        assert box.max_x > box.min_x
        assert (box.min_x + box.max_x) / 2.0 == pytest.approx(EAST, rel=0.0, abs=1e-9)

    def test_padding_leaves_the_healthy_axis_alone(self) -> None:
        column = [(EAST, NORTH), (EAST, NORTH + 40.0), (EAST, NORTH + 80.0)]
        box = self.collinear(column).bbox
        assert box.min_y == pytest.approx(NORTH, rel=0.0, abs=1e-9)
        assert box.max_y == pytest.approx(NORTH + 80.0, rel=0.0, abs=1e-9)

    def test_a_zero_height_box_is_padded_too(self) -> None:
        row = [(EAST, NORTH), (EAST + 40.0, NORTH), (EAST + 80.0, NORTH)]
        box = self.collinear(row).bbox
        assert box.padded is True
        assert box.max_y > box.min_y
        assert (box.min_y + box.max_y) / 2.0 == pytest.approx(NORTH, rel=0.0, abs=1e-9)

    def test_coincident_vertices_pad_both_axes(self) -> None:
        box = self.collinear([(EAST, NORTH), (EAST, NORTH)]).bbox
        assert box.padded is True
        assert box.max_x > box.min_x
        assert box.max_y > box.min_y

    def test_an_empty_vertex_set_still_yields_a_usable_box(self) -> None:
        pslg = FakePslg(
            vertices=np.zeros((0, 2), dtype=np.float64),
            chains=(),
            chain_indices=np.zeros((0,), dtype=np.uint32),
        )
        box = scene_module().build_scene(pslg, None, ok=False).bbox
        assert box.padded is True
        assert box.max_x > box.min_x
        assert box.max_y > box.min_y
        assert all(np.isfinite([box.min_x, box.min_y, box.max_x, box.max_y]))

    @pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
    def test_a_non_finite_coordinate_is_refused(self, mesh: FakeMesh, bad: float) -> None:
        # The PSLG validator's stage 3 rejects these, so reaching here means
        # something upstream is wrong. Refusing costs a line; the failure mode
        # without it is an SVG that renders as nothing, in silence.
        vertices = MESH_VERTICES.copy()
        vertices[2, 1] = bad
        pslg = make_pslg([([0, 1, 2, 3], Role.Outer, False)], vertices=vertices)
        broken = FakeMesh(vertices, MESH_TRIANGLES, MESH_MASKS)
        with pytest.raises(ValueError, match="finite"):
            scene_module().build_scene(pslg, broken, closed_roles=CLOSED)


class TestBackendIntroducedVertices:
    """`indexed_mesh.hpp` guarantee 4 says the mesh's vertex array *begins
    with* the PSLG's -- not that it equals it. The backend may append points of
    its own, and the scene draws the mesh's array whenever there is a mesh.

    Every other fixture here builds the PSLG and the mesh from one array, so
    the two sources are indistinguishable: measured, `source = pslg.vertices`
    passed the whole suite. It is not a cosmetic difference. Under that mutant
    6b-ii's viewport transform indexes triangle corners into an array that is
    too short, the bbox omits every backend-introduced point, and the finiteness
    guard inspects a strictly narrower set than the one that gets drawn -- which
    is why the third test below is here and not in `TestBoundingBox`.
    """

    def test_the_scene_carries_the_meshs_vertices_not_the_inputs(
        self, pslg: FakePslg, mesh_with_backend_vertex: FakeMesh
    ) -> None:
        scene = scene_module().build_scene(pslg, mesh_with_backend_vertex, closed_roles=CLOSED)
        assert len(scene.vertices) == 7
        assert np.array_equal(np.asarray(scene.vertices), MESH_VERTICES_PLUS_BACKEND)

    def test_the_bbox_covers_a_backend_introduced_vertex(
        self, pslg: FakePslg, mesh_with_backend_vertex: FakeMesh
    ) -> None:
        # Absolute tolerance, for `TestBoundingBox`'s reason: at an easting of
        # 4.3e5 the default relative one is +-0.43 m.
        box = scene_module().build_scene(
            pslg, mesh_with_backend_vertex, closed_roles=CLOSED
        ).bbox
        assert box.max_x == pytest.approx(EAST + 160.0, rel=0.0, abs=1e-9)
        assert box.max_y == pytest.approx(NORTH + 130.0, rel=0.0, abs=1e-9)

    def test_the_finiteness_guard_inspects_the_meshs_vertices(self, pslg: FakePslg) -> None:
        vertices = MESH_VERTICES_PLUS_BACKEND.copy()
        vertices[6, 0] = np.nan
        broken = FakeMesh(vertices, MESH_TRIANGLES, MESH_MASKS)
        with pytest.raises(ValueError, match="finite"):
            scene_module().build_scene(pslg, broken, closed_roles=CLOSED)

    def test_a_scene_without_a_mesh_falls_back_to_the_input_vertices(
        self, pslg: FakePslg
    ) -> None:
        # The other direction: with no drawable mesh there is no longer array to
        # prefer, and the input's own vertices are what the picture is of.
        scene = scene_module().build_scene(pslg, None, ok=False, closed_roles=CLOSED)
        assert np.array_equal(np.asarray(scene.vertices), MESH_VERTICES)
