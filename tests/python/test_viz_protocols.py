"""`tin_engine.viz.protocols`: increment 6a's half of the renderer seam.

Committed RED with the binding suite. `protocols.py` is the only `viz/` module
in 6a; `scene.py`, `svg.py`, `style.py`, `fixtures.py` and `rasputin draw` are
6b and are not tested here.

Two things make this module worth a suite of its own, small as it is:

* **`viz/` never imports `_core`.** That is what lets the whole renderer be
  unit-tested against a twelve-line fake mesh with no compiled extension in the
  process, and what lets 6a and 6b be built in either order. It is checked by
  reading the source rather than by importing, because `tin_engine/__init__.py`
  imports `_core` itself -- so a `sys.modules` check would be asserting
  something about the package, not about `viz/`.
* **The protocols must actually describe the bound types.** A `MeshLike` whose
  member names drifted from `IndexedMesh2`'s accessors still type-checks, still
  passes every renderer test against a fake, and fails only at the composition
  root. Joining the two halves is exactly this file's job.

The module is imported inside the tests rather than at module scope, so a
missing `viz/` fails as a named `ModuleNotFoundError` in one test instead of an
error before collection.
"""

from __future__ import annotations

import ast
import importlib
import typing
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from tin_engine import _core
from tin_engine.features import DEFAULT_VOCABULARY

REPO_ROOT = Path(__file__).resolve().parents[2]
VIZ_DIR = REPO_ROOT / "src_python" / "tin_engine" / "viz"

EAST = 430_000.0
NORTH = 6_900_000.0

# The properties position of a chain spec is a mask of property bits, not a
# flag: `0` is the empty set -- an unclassified constraint -- and a named bit
# comes from a vocabulary rather than from a literal, so that the number at
# the call site carries its units. A `bool` here is the pre-migration
# `is_river` spelling and `build_pslg` refuses it.
NO_PROPERTIES = 0
RIVER = DEFAULT_VOCABULARY.mask("river")

MESH_MEMBERS = {"vertices", "triangles", "constrained_edges", "triangle_count", "empty"}
PSLG_MEMBERS = {"vertices", "chains", "chain_indices", "indices_of"}
CHAIN_MEMBERS = {"begin", "count", "role", "properties"}


def protocols() -> ModuleType:
    return importlib.import_module("tin_engine.viz.protocols")


def protocol_members(proto: type) -> set[str]:
    """The declared members of a Protocol, across the versions CI runs.

    `typing.get_protocol_members` is 3.13+; `__protocol_attrs__` is 3.12+ and
    covers the floor. A third branch on `typing._get_protocol_attrs` was
    carried for 3.11 and is gone with it -- private API kept for an
    interpreter the project no longer supports, which nothing would have
    flagged, since --cov does not measure the test tree.

    Written out rather than hand-listing the annotations, so that a member
    added to the protocol without a matching accessor on the bound type is
    caught here.
    """
    getter = getattr(typing, "get_protocol_members", None)
    if getter is not None:
        return set(getter(proto))  # type: ignore[arg-type]
    return set(proto.__protocol_attrs__)  # type: ignore[attr-defined]


@pytest.fixture
def mesh() -> Any:
    vertices = np.array(
        [
            [EAST, NORTH],
            [EAST + 100.0, NORTH],
            [EAST + 100.0, NORTH + 80.0],
            [EAST, NORTH + 80.0],
            [EAST + 20.0, NORTH + 40.0],
            [EAST + 80.0, NORTH + 40.0],
        ],
        dtype=np.float64,
    )
    result = _core.build_pslg(
        vertices,
        [
            ([0, 1, 2, 3], _core.ChainRole.Outer, NO_PROPERTIES),
            ([4, 5], _core.ChainRole.Breakline, RIVER),
        ],
    )
    assert result.ok, [d.message for d in result.diagnostics]
    outcome = _core.triangulate(result.pslg)
    assert outcome.ok(), outcome.message
    return outcome.mesh, result.pslg


class TestModuleShape:
    def test_viz_package_directory_exists(self) -> None:
        assert VIZ_DIR.is_dir()

    @pytest.mark.parametrize("name", ["MeshLike", "PslgLike", "ChainLike"])
    def test_declares_the_three_protocols(self, name: str) -> None:
        proto = getattr(protocols(), name)
        assert getattr(proto, "_is_protocol", False), f"{name} is not a typing.Protocol"

    def test_nothing_else_is_exported(self) -> None:
        # Risk 3 applied to the seam: protocols.py is ~25 lines and the place a
        # later increment will be tempted to park a helper. Anything public
        # beyond the three protocols is a design change.
        module = protocols()
        public = {
            name
            for name, value in vars(module).items()
            if not name.startswith("_") and getattr(value, "__module__", None) == module.__name__
        }
        assert public == {"MeshLike", "PslgLike", "ChainLike"}


class TestProtocolMembers:
    def test_mesh_like_members(self) -> None:
        assert protocol_members(protocols().MeshLike) == MESH_MEMBERS

    def test_pslg_like_members(self) -> None:
        assert protocol_members(protocols().PslgLike) == PSLG_MEMBERS

    def test_chain_like_members(self) -> None:
        assert protocol_members(protocols().ChainLike) == CHAIN_MEMBERS

    def test_no_protocol_admits_a_mutating_member(self) -> None:
        # The protocols describe read-only views. A `set_`/`clear`/`append`
        # member would let the renderer write through to a core buffer.
        forbidden = ("set_", "clear", "append", "resize", "insert")
        for proto in (protocols().MeshLike, protocols().PslgLike, protocols().ChainLike):
            for member in protocol_members(proto):
                assert not member.startswith(forbidden), f"{proto.__name__}.{member}"


class TestVizDoesNotImportCore:
    """Checked by reading the source: `tin_engine/__init__.py` imports `_core`,
    so importing anything under `tin_engine.viz` imports it transitively and a
    `sys.modules` assertion would pass for the wrong reason."""

    def viz_sources(self) -> list[Path]:
        return sorted(VIZ_DIR.glob("*.py"))

    def test_there_is_something_to_check(self) -> None:
        # A probe able to fail: an empty glob would let the next test pass
        # while measuring nothing.
        assert self.viz_sources(), f"no Python sources under {VIZ_DIR}"

    def test_no_viz_module_imports_the_extension(self) -> None:
        offenders: list[str] = []
        for path in self.viz_sources():
            for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
                if isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom):
                    names = [node.module or ""] + [alias.name for alias in node.names]
                else:
                    continue
                if any("_core" in name for name in names):
                    offenders.append(f"{path.name}:{node.lineno}")
        assert offenders == []


class TestBoundTypesSatisfyTheProtocols:
    """The join. Without this, `MeshLike` and `IndexedMesh2` can drift apart
    indefinitely and both halves of the increment stay green."""

    def test_indexed_mesh_satisfies_mesh_like(self, mesh: Any) -> None:
        real, _ = mesh
        missing = [m for m in protocol_members(protocols().MeshLike) if not hasattr(real, m)]
        assert missing == []

    def test_pslg_satisfies_pslg_like(self, mesh: Any) -> None:
        _, pslg = mesh
        missing = [m for m in protocol_members(protocols().PslgLike) if not hasattr(pslg, m)]
        assert missing == []

    def test_bound_chain_record_satisfies_chain_like(self, mesh: Any) -> None:
        _, pslg = mesh
        chain = pslg.chains[0]
        missing = [m for m in protocol_members(protocols().ChainLike) if not hasattr(chain, m)]
        assert missing == []
