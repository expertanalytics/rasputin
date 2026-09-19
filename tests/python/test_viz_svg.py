"""`tin_engine.viz.svg`, `viz.style` and `viz.fixtures`: increment 6b-ii.

Committed RED, before any of the three modules exists. The modules under test
are imported inside test bodies and inside fixtures, so the intended failure is
a named `ModuleNotFoundError` per test rather than a collection error that takes
the rest of `tests/python/` down with it. `viz/scene.py` *does* exist -- it is
6b-i's, merged -- so `build_scene` is imported at module scope and every scene
here is a real one.

`06-cdt-viewer.md` is explicit that this suite is **not** invariant-critical and
carries **no mutation round**: a renderer's failure mode is a wrong-looking
picture and a person catches that instantly. It is equally explicit that a
**golden-file SVG comparison is rejected outright**, and that **no assertion is
made on any style string**. Nothing below reads a colour, a stroke width, a dash
array or the stylesheet. What is asserted is structure and arithmetic:

* the document parses as XML in the SVG namespace;
* element counts equal the scene's primitive counts;
* the viewport transform flips y, preserves aspect, and centres the bbox;
* the failure presentations are non-blank and carry their status text;
* the gallery is the eight fixtures the design names, and each one is the shape
  its table row claims.

The `class` attribute *is* asserted, and the distinction matters. A class token
names **which stroke class an edge belongs to** -- `constrained`, `role-hole`,
`finding`, and at most one property token such as `river` -- which is
structure, and it is the only way to check that role colouring reaches the
right edges without naming a colour. The stylesheet that turns each token into
a colour gets no test at all, per the design.

**At most one property token**, chosen by `SvgStyle.property_strokes` in
descending priority: an edge carries a SET of properties and a polyline carries
one stroke, because two overlaid strokes on one line read as a rendering defect
rather than as two features. `TestPropertyStrokes` and
`TestStrokeClasses` hold both halves of that.

## Tolerances are absolute, everywhere a coordinate is compared

`pytest.approx` defaults to a **relative** tolerance. At an easting of 4.3e5
that is +-0.43 m, and a bbox pushed a quarter-metre off centre survived an
entire mutation round in 6b-i before `@tester` switched to `rel=0.0, abs=1e-9`;
`TestBoundingBox`'s docstring in `test_viz_scene.py` records the incident. The
viewport is where that bites next, because it is the one place UTM-magnitude
world coordinates and SVG user units meet: a relative tolerance on the world
side is 0.43 m of slack and on the SVG side is several pixels. Every comparison
in `TestViewport` and `TestGeometryElements` is therefore `rel=0.0, abs=...`.

The two predicates in `TestGallery` are the deliberate exception and are
*scaled* rather than absolute -- collinearity and near-collinearity are
statements about a normalised cross product, so an absolute bound on an
unnormalised one at UTM magnitudes would mean nothing. Both predicates have a
self-test that fails if the predicate stops discriminating.

## The API this suite pins

```
style.SvgStyle(width=..., height=..., margin=..., header_height=...,
               show_vertices=False)          # frozen Pydantic V2
svg.viewport(bbox, style) -> Viewport        # .scale, .point(x, y)
svg.render_svg(scene, style=None, *, title="", status="", message="",
               labels=False) -> str
fixtures.GALLERY: Mapping[str, Fixture]
```

`status` and `message` are **passthrough strings**: `viz/` may not name
`_core.CdtStatus`, so `cli.py` composes the status name and `describe(status)`
and hands the result over as text. Same reasoning as 6b-i's `ok` flag.

`Fixture` satisfies `PslgLike` structurally and its chain roles are the strings
`"outer"`, `"hole"` and `"breakline"`. That is forced, not chosen: `viz/` never
imports `_core`, so a fixture cannot name `ChainRole`, and `cli.py` -- the
composition root -- maps the strings to the enum when it calls `build_pslg`.
Keeping the fixture itself a `PslgLike` is what lets the `degenerate` fixture be
drawn at all; see `test_cli_draw.py::TestFailurePresentation`.
"""

from __future__ import annotations

import ast
import enum
import importlib
import math
import re
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest

from tin_engine.viz.scene import build_scene

REPO_ROOT = Path(__file__).resolve().parents[2]
VIZ = REPO_ROOT / "src_python" / "tin_engine" / "viz"

SVG_NS = "http://www.w3.org/2000/svg"

EAST = 430_000.0
NORTH = 6_900_000.0

GALLERY_NAMES = (
    "catchment",
    "sliver-fan",
    "corner-hole",
    "hole-in-hole",
    "breakline-chain",
    "river",
    "not-noded",
    "degenerate",
)

ROLE_NAMES = frozenset({"outer", "hole", "breakline"})

# Bare masks, as `cli.py` hands them down. `viz/` may not name a vocabulary --
# `fixtures.py` and `style.py` import nothing first-party at all, which
# `TestModuleIsolation` pins -- so the renderer knows bit positions, the
# stylesheet knows tokens, and only the composition root knows that bit 0 means
# "river".
#
# Two disjoint properties, because an edge that is both a road and a river is
# the case this whole widening exists for and the one a single bit cannot state.
RIVER_BIT = 0
ROAD_BIT = 1
NO_PROPERTIES = 0
RIVER = 1 << RIVER_BIT
ROAD = 1 << ROAD_BIT

#: The precedence the gallery draws under, in DESCENDING priority: first match
#: wins. Road over river, which is deliberately NOT bit order -- priority
#: belongs to the stylesheet, not to how a vocabulary chose to number its
#: features, and a renderer that sorted by bit would pass a list that happened
#: to agree with the numbering.
GALLERY_STROKES = ((ROAD_BIT, "road"), (RIVER_BIT, "river"))


class Role(enum.Enum):
    """A stand-in for `_core.ChainRole`, which this suite may not import.

    A role is opaque to `viz/`: the scene compares roles, and the renderer needs
    a *name* for the stroke class. `getattr(role, "name", str(role))` is what
    serves both an enum and the plain strings `fixtures.py` authors, and
    `test_a_string_role_names_its_own_stroke_class` pins the second half.
    """

    Outer = "Outer"
    Hole = "Hole"
    Breakline = "Breakline"


CLOSED = (Role.Outer, Role.Hole)


@dataclass(frozen=True)
class FakeChain:
    begin: int
    count: int
    role: object
    properties: int


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
# 6b-i's rectangle: 100 x 80 m at UTM 33N magnitudes, a horizontal river
# breakline across the middle, six triangles, eleven undirected edges, five of
# them constrained. Vertex 6 is a backend-introduced point -- `indexed_mesh.hpp`
# guarantee 4 says the mesh array *begins with* the PSLG's, it does not equal it
# -- placed strictly INSIDE the rectangle so that the bbox stays exactly
# 100 x 80 and every viewport number below can be computed by hand.
#
#   3 +-------------------+ 2
#     | \               / |
#     |   4 ==========5   |        4=5 is the breakline (properties: RIVER)
#     | /      x6       \ |
#   0 +-------------------+ 1
#
# The header-band counts are then all different -- 7 vertices, 6 triangles, 5
# constrained edges, 0 findings -- so a header that labels one count with
# another's number cannot pass `mentions_count`.
MESH_VERTICES = np.array(
    [
        [EAST + 0.0, NORTH + 0.0],
        [EAST + 100.0, NORTH + 0.0],
        [EAST + 100.0, NORTH + 80.0],
        [EAST + 0.0, NORTH + 80.0],
        [EAST + 20.0, NORTH + 40.0],
        [EAST + 80.0, NORTH + 40.0],
        [EAST + 50.0, NORTH + 20.0],
    ],
    dtype=np.float64,
)

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

# Hand-written, per `indexed_mesh.hpp`: bit e is the edge (v[e], v[(e+1)%3]).
MESH_MASKS = np.array([0b001, 0b010, 0b100, 0b001, 0b001, 0b010], dtype=np.uint8)

CONSTRAINED_PAIRS = frozenset({(0, 1), (1, 2), (2, 3), (0, 3), (4, 5)})

VERTEX_COUNT = 7
TRIANGLE_COUNT = 6
EDGE_COUNT = 11
CONSTRAINED_COUNT = 5

# --- The viewport, computed by hand ----------------------------------------
#
# With width 1000, height 800, margin 24 and a 72-unit header band, the map area
# is x in [24, 976] (952 wide) and y in [96, 776] (680 tall). The bbox is
# 100 x 80, so a uniform scale is min(952/100, 680/80) = min(9.52, 8.5) = 8.5:
# the drawing is 850 x 680, letterboxed horizontally and exactly filling the map
# area's height. Its centre is the map area's centre, (500, 436).
WIDTH = 1000
HEIGHT = 800
MARGIN = 24.0
HEADER = 72.0
SCALE = 8.5
CENTRE_X = 500.0
CENTRE_Y = 436.0
MAP_LEFT = 75.0  # 500 - 50 m * 8.5
MAP_RIGHT = 925.0
MAP_TOP = 96.0  # 436 - 40 m * 8.5, the world's MAXIMUM y
MAP_BOTTOM = 776.0


def viz_module(name: str) -> ModuleType:
    return importlib.import_module(f"tin_engine.viz.{name}")


def make_style(**overrides: Any) -> Any:
    """The style every hand-computed viewport number above was derived for."""
    kwargs: dict[str, Any] = {
        "width": WIDTH,
        "height": HEIGHT,
        "margin": MARGIN,
        "header_height": HEADER,
    }
    kwargs.update(overrides)
    return viz_module("style").SvgStyle(**kwargs)


def strokes(*pairs: tuple[int, str]) -> tuple[Any, ...]:
    """`PropertyStroke`s from `(bit, token)` pairs, in the order given.

    The order IS the draw priority -- first match wins -- so a helper that
    sorted would destroy the one thing these tests are about.
    """
    stroke = viz_module("style").PropertyStroke
    return tuple(stroke(bit=bit, token=token) for bit, token in pairs)


def make_pslg(
    chains: Sequence[tuple[Sequence[int], object, int]],
    vertices: npt.NDArray[np.float64] | None = None,
) -> FakePslg:
    flat: list[int] = []
    records: list[FakeChain] = []
    for indices, role, properties in chains:
        records.append(FakeChain(len(flat), len(indices), role, properties))
        flat.extend(int(i) for i in indices)
    return FakePslg(
        vertices=MESH_VERTICES if vertices is None else vertices,
        chains=tuple(records),
        chain_indices=np.array(flat, dtype=np.uint32),
    )


# --- XML helpers ------------------------------------------------------------


def tag(name: str) -> str:
    return f"{{{SVG_NS}}}{name}"


def parse(document: str) -> ET.Element:
    """Parse, with a failure message that shows the head of the document."""
    try:
        return ET.fromstring(document)
    except ET.ParseError as exc:  # pragma: no cover - only on a red suite
        pytest.fail(f"render_svg did not emit well-formed XML: {exc}\n{document[:400]}")


def group(root: ET.Element, gid: str) -> ET.Element:
    matches = [g for g in root.iter(tag("g")) if g.get("id") == gid]
    assert len(matches) == 1, f"expected exactly one <g id={gid!r}>, found {len(matches)}"
    return matches[0]


def group_or_none(root: ET.Element, gid: str) -> ET.Element | None:
    matches = [g for g in root.iter(tag("g")) if g.get("id") == gid]
    return matches[0] if matches else None


def children(root: ET.Element, gid: str, name: str) -> list[ET.Element]:
    return list(group(root, gid).iter(tag(name)))


def text_of(element: ET.Element) -> str:
    """All text under an element, whitespace collapsed."""
    return " ".join("".join(element.itertext()).split())


def classes(element: ET.Element) -> frozenset[str]:
    return frozenset((element.get("class") or "").split())


def property_tokens(element: ET.Element, vocabulary: Sequence[str]) -> list[str]:
    """Which of `vocabulary`'s tokens the element carries.

    A list rather than a set, because the assertion the property strokes need
    is a COUNT -- "at most one" -- and `frozenset` silently answers a different
    question about a duplicate.
    """
    return [t for t in (element.get("class") or "").split() if t in set(vocabulary)]


def edge_at(scene: Any, a: int, b: int) -> Any:
    matches = [e for e in scene.edges if (int(e.a), int(e.b)) == (a, b)]
    assert len(matches) == 1, f"edge {(a, b)} appears {len(matches)} times"
    return matches[0]


def header_line(document: ET.Element, needle: str) -> ET.Element:
    """The one `<text>` of the header band that carries `needle`.

    A class token is structure, not taste (see `svg.py`'s own convention), so
    `alarm` is assertable where a colour is not -- and it is the design's stated
    mechanism for "degeneracy and emptiness are the loudest thing on the page".
    Asserting it needs the individual line rather than the band's collapsed
    text, which is what this returns.
    """
    lines = [t for t in group(document, "header").iter(tag("text")) if needle in text_of(t)]
    assert len(lines) == 1, f"expected one header line saying {needle!r}, found {len(lines)}"
    return lines[0]


def numbers(text: str) -> list[float]:
    return [float(m) for m in re.findall(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", text)]


def mentions_count(text: str, word: str, n: int) -> bool:
    """`word` and `n` adjacent, in either order.

    Deliberately loose about the separator and the order -- "Triangles: 6" and
    "6 triangles" are both fine -- and deliberately strict about adjacency, so
    that a header labelling the vertex count with the triangle count fails. The
    fixture's four counts are pairwise distinct, which is what gives that teeth.
    """
    forward = re.search(rf"(?i)\b{word}\b\D{{0,3}}{n}\b", text)
    backward = re.search(rf"(?i)\b{n}\b\D{{0,3}}{word}\b", text)
    return bool(forward or backward)


def edge_element(root: ET.Element, a: int, b: int, scene: Any) -> ET.Element:
    """The drawn element for scene edge `(a, b)`.

    Located by position: the renderer emits one element per `Scene.edge`, in
    `Scene.edges` order. That is the only ordering assumption this suite makes,
    and it is made rather than inventing per-element ids, which would be a
    naming convention with no reader.
    """
    order = [(int(e.a), int(e.b)) for e in scene.edges]
    index = order.index((a, b))
    return children(root, "edges", "line")[index]


# --- Scenes -----------------------------------------------------------------


@pytest.fixture
def pslg() -> FakePslg:
    """A closed outer ring and an open river breakline, as 6b-i's suite uses."""
    return make_pslg([([0, 1, 2, 3], Role.Outer, NO_PROPERTIES), ([4, 5], Role.Breakline, RIVER)])


@pytest.fixture
def mesh() -> FakeMesh:
    return FakeMesh(MESH_VERTICES, MESH_TRIANGLES, MESH_MASKS)


@pytest.fixture
def scene(pslg: FakePslg, mesh: FakeMesh) -> Any:
    return build_scene(pslg, mesh, closed_roles=CLOSED)


@pytest.fixture
def empty_scene(pslg: FakePslg) -> Any:
    """`Ok` with zero triangles: increment 4's risk 5, the silent mode."""
    empty = FakeMesh(
        MESH_VERTICES, np.zeros((0, 3), dtype=np.uint32), np.zeros((0,), dtype=np.uint8)
    )
    return build_scene(pslg, empty, ok=True, closed_roles=CLOSED)


@pytest.fixture
def failed_scene(pslg: FakePslg) -> Any:
    return build_scene(pslg, None, ok=False, closed_roles=CLOSED)


@pytest.fixture
def finding_scene(pslg: FakePslg) -> Any:
    """Triangle 0 gains bit 2 -- edge (0, 5), which no chain contains."""
    masks = MESH_MASKS.copy()
    masks[0] = 0b101
    broken = FakeMesh(MESH_VERTICES, MESH_TRIANGLES, masks)
    return build_scene(pslg, broken, closed_roles=CLOSED)


@pytest.fixture
def document(scene: Any) -> ET.Element:
    # The gallery's precedence list, supplied here because `SvgStyle`'s default
    # is `()` and must stay `()`: a default naming `river` would be policy in
    # the one module that declares it holds none, exactly as `closed_roles`
    # is the composition root's in 6b-i.
    return parse(
        viz_module("svg").render_svg(
            scene, style=make_style(property_strokes=strokes(*GALLERY_STROKES))
        )
    )


class TestFixtureSanity:
    """Probes that can fail on their own.

    Every count assertion below is only as good as the scene it counts, and a
    header check whose four numbers were all `6` would pass a renderer that
    labelled them at random. These run first so that a red suite says the
    fixture drifted rather than that the renderer is wrong.
    """

    def test_the_scene_has_the_primitive_counts_the_assertions_assume(
        self, scene: Any
    ) -> None:
        assert len(scene.vertices) == VERTEX_COUNT
        assert len(scene.triangles) == TRIANGLE_COUNT
        assert len(scene.edges) == EDGE_COUNT
        constrained = frozenset((int(e.a), int(e.b)) for e in scene.edges if e.constrained)
        assert constrained == CONSTRAINED_PAIRS
        assert len(constrained) == CONSTRAINED_COUNT
        assert list(scene.findings) == []

    def test_the_four_header_counts_are_pairwise_distinct(self) -> None:
        counts = [VERTEX_COUNT, TRIANGLE_COUNT, CONSTRAINED_COUNT, 0]
        assert len(set(counts)) == len(counts)

    def test_the_bbox_is_the_one_the_viewport_numbers_were_computed_for(
        self, scene: Any
    ) -> None:
        box = scene.bbox
        assert box.padded is False
        assert box.max_x - box.min_x == pytest.approx(100.0, rel=0.0, abs=1e-9)
        assert box.max_y - box.min_y == pytest.approx(80.0, rel=0.0, abs=1e-9)

    def test_the_backend_vertex_does_not_widen_the_bbox(self, scene: Any) -> None:
        # If vertex 6 ever moves outside the rectangle, SCALE stops being 8.5
        # and every hand-computed number in TestViewport silently drifts.
        assert scene.bbox.max_x == pytest.approx(EAST + 100.0, rel=0.0, abs=1e-9)
        assert scene.bbox.max_y == pytest.approx(NORTH + 80.0, rel=0.0, abs=1e-9)

    def test_coordinates_are_at_utm_magnitudes(self) -> None:
        assert MESH_VERTICES[:, 0].min() >= 1.0e5
        assert MESH_VERTICES[:, 1].min() >= 1.0e6


class TestModuleIsolation:
    """`viz/` never imports `_core`, and `cli.py` is the sole composition root.

    `test_viz_scene.py::TestModuleIsolation` pins this for `scene.py`; the three
    new modules need the same, and for the same reason -- it is what lets this
    whole suite run against hand-built dataclasses. Checked by reading the
    source, because `tin_engine/__init__.py` imports `_core` itself and a
    `sys.modules` assertion would therefore pass for the wrong reason.
    """

    def imports(self, module: str) -> list[str]:
        path = VIZ / f"{module}.py"
        assert path.is_file(), f"{path} does not exist"
        names: list[str] = []
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                names.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                names.append("." * node.level + (node.module or ""))
        return names

    @pytest.mark.parametrize("module", ["svg", "style", "fixtures"])
    def test_the_module_exists(self, module: str) -> None:
        assert (VIZ / f"{module}.py").is_file()

    @pytest.mark.parametrize("module", ["svg", "style", "fixtures"])
    def test_it_does_not_import_the_extension(self, module: str) -> None:
        assert [name for name in self.imports(module) if "_core" in name] == []

    def test_the_renderer_imports_only_its_own_siblings(self) -> None:
        allowed = {".scene", ".style", ".protocols", "tin_engine.viz.scene"}
        first_party = [
            name
            for name in self.imports("svg")
            if name.startswith((".", "tin_engine")) and name not in allowed
        ]
        assert first_party == []

    def test_the_style_model_imports_no_first_party_module(self) -> None:
        # `SvgStyle` is the renderer's input and knows nothing about a scene.
        assert [n for n in self.imports("style") if n.startswith((".", "tin_engine"))] == []

    def test_the_gallery_imports_no_first_party_module(self) -> None:
        # Declarative data. A fixture that imported the scene builder would be
        # a fixture nobody can read as data.
        assert [n for n in self.imports("fixtures") if n.startswith((".", "tin_engine"))] == []

    def test_the_cli_names_the_enum_that_viz_may_not(self) -> None:
        # Named for what it actually checks: that `ChainRole` appears in
        # `cli.py` at all, which is the other half of the four `_core` checks
        # above -- they establish that no `viz/` module names it, this one that
        # the composition root does. *Exclusivity* is theirs; that the mapping
        # is total is `test_cli_draw.py`'s, which may import `_core` where this
        # suite may not.
        source = (REPO_ROOT / "src_python" / "tin_engine" / "cli.py").read_text(encoding="utf-8")
        assert "ChainRole" in source


class TestPackageExports:
    """`viz/__init__.py`'s re-exports, which `06-cdt-viewer.md` names."""

    def test_the_package_re_exports_the_three_public_names(self) -> None:
        package = importlib.import_module("tin_engine.viz")
        for name in ("render_svg", "Scene", "SvgStyle"):
            assert hasattr(package, name), f"tin_engine.viz does not export {name}"

    def test_the_exports_are_declared(self) -> None:
        package = importlib.import_module("tin_engine.viz")
        assert set(getattr(package, "__all__", ())) >= {"render_svg", "Scene", "SvgStyle"}


class TestSvgStyle:
    """Frozen Pydantic V2, and nothing about colour.

    The stylesheet is untested by design, so what is left to pin is that the
    model is immutable, that it rejects a field it does not know, and that the
    two numbers the viewport divides by cannot be zero or negative.
    """

    def test_a_default_style_constructs(self) -> None:
        style = viz_module("style").SvgStyle()
        assert style.width > 0
        assert style.height > 0

    def test_it_is_frozen(self) -> None:
        import pydantic

        style = viz_module("style").SvgStyle()
        with pytest.raises(pydantic.ValidationError):
            style.width = 42

    def test_an_unknown_field_is_refused(self) -> None:
        import pydantic

        # extra="forbid": a typo in a style override is a silent no-op
        # otherwise, and the symptom is a picture that looks almost right.
        with pytest.raises(pydantic.ValidationError):
            viz_module("style").SvgStyle(colour_of_everything="beige")

    @pytest.mark.parametrize(("field", "value"), [("width", 0), ("height", -1)])
    def test_a_non_positive_canvas_dimension_is_refused(self, field: str, value: int) -> None:
        import pydantic

        with pytest.raises(pydantic.ValidationError):
            viz_module("style").SvgStyle(**{field: value})

    def test_a_negative_margin_is_refused(self) -> None:
        import pydantic

        with pytest.raises(pydantic.ValidationError):
            viz_module("style").SvgStyle(margin=-1.0)

    def test_vertices_are_off_by_default(self) -> None:
        assert viz_module("style").SvgStyle().show_vertices is False


class TestPropertyStrokes:
    """`SvgStyle.property_strokes`: which property is drawn, and in what order.

    Structure, not taste, which is why it lives in `style.py` and not in the
    stylesheet: WHICH of an edge's properties gets a token is a decision about
    the document, and the colour that token resolves to is the CSS's. The list
    is in descending draw priority and the first match wins.

    `PropertyStroke.token`'s `^[a-z][a-z0-9_-]*$` is the pattern that is
    actually load-bearing, and this is the suite over it. `_edge_classes`
    (`svg.py:166`) joins TOKENS into one string and `_edges` (`svg.py:225`)
    interpolates that string into a `class="..."` attribute UNESCAPED -- only
    `_text` escapes -- and neither ever sees an `EdgeProperty.name`, because
    `viz/` may import no vocabulary at all. The one bridge is `cli.py:84`,
    which constructs a `PropertyStroke` from a feature name, so every name that
    can reach the attribute is re-validated here. A token containing a quote
    ends the attribute, and a stylesheet read from a configuration file is
    untrusted input.

    `features.EdgeProperty.name`'s narrower `^[a-z][a-z0-9_]*$` is defence in
    depth behind this pattern, not the other way round -- a quote-carrying
    feature name dies at this model whatever that one says.
    """

    def test_the_default_is_empty(self) -> None:
        # A default naming `river` would be policy in the module whose own
        # docstring says colours and vocabulary are somebody else's.
        assert viz_module("style").SvgStyle().property_strokes == ()

    def test_a_stroke_carries_a_bit_and_a_token(self) -> None:
        stroke = viz_module("style").PropertyStroke(bit=RIVER_BIT, token="river")
        assert stroke.bit == RIVER_BIT
        assert stroke.token == "river"

    def test_the_declared_order_is_preserved(self) -> None:
        # The order IS the priority. A model that sorted -- by bit, by token --
        # would silently replace the stylesheet's ruling with the vocabulary's
        # numbering, which is the thing the ruling exists to decouple from.
        style = make_style(property_strokes=strokes(*GALLERY_STROKES))
        assert [(s.bit, s.token) for s in style.property_strokes] == list(GALLERY_STROKES)

    @pytest.mark.parametrize("bit", [-1, 32, 33, 64])
    def test_a_bit_outside_the_ceiling_is_refused(self, bit: int) -> None:
        import pydantic

        # 32 is the C++ ceiling (`EdgeProperties::kMaxProperties`), restated at
        # the layer that meets untrusted input. A stroke on bit 32 is a stroke
        # no edge can ever carry, so it is a typo, caught at construction.
        with pytest.raises(pydantic.ValidationError):
            viz_module("style").PropertyStroke(bit=bit, token="river")

    @pytest.mark.parametrize("bit", [0, 31])
    def test_the_boundary_bits_are_accepted(self, bit: int) -> None:
        # Able to fail on its own: a range check written with the wrong
        # comparison refuses 31, and the refusals above would still pass.
        assert viz_module("style").PropertyStroke(bit=bit, token="river").bit == bit

    @pytest.mark.parametrize(
        "token",
        ['ri"ver', "ri ver", "River", "1river", "_river", "river>", "", "riv<er",
         "river\n", "riv\ner", "riv\u00e9r"],
        ids=["quote", "space", "capital", "leading-digit", "leading-underscore",
             "gt", "empty", "lt", "trailing-newline", "embedded-newline",
             "non-ascii"],
    )
    def test_a_token_that_would_escape_the_class_attribute_is_refused(
        self, token: str
    ) -> None:
        """The hostile-input set for the pattern that is the boundary.

        The two newline cases and the non-ASCII one mirror
        `test_features.py`'s, and belong here more than there: Python's `re`
        matches `$` BEFORE a trailing newline, so a validator hand-written as
        `re.match(r"^[a-z][a-z0-9_-]*$", token)` accepts `"river\\n"` while
        `Field(pattern=...)` -- pydantic's Rust engine, where `$` is
        end-of-haystack -- refuses it. `"river\\n"` alone carries that kill:
        without it, the hand-rolled mutant survives this class while dying in
        `test_features.py`, over a pattern no attribute depends on. The
        embedded-newline and non-ASCII cases are its companions, pinning the
        pattern's shape rather than discriminating against that mutant --
        probed against pydantic 2.13.5, both are refused by the real model and
        by the mutant alike.
        """
        import pydantic

        with pytest.raises(pydantic.ValidationError):
            viz_module("style").PropertyStroke(bit=RIVER_BIT, token=token)

    @pytest.mark.parametrize("token", ["river", "road-bridge", "river_bank", "r2"])
    def test_a_css_shaped_token_is_accepted(self, token: str) -> None:
        # The other half: a pattern that refused everything would pass every
        # refusal above. A hyphen is legal here and is NOT in
        # `EdgeProperty.name`'s pattern -- a CSS class may contain one and a
        # Python-side feature name may not.
        assert viz_module("style").PropertyStroke(bit=RIVER_BIT, token=token).token == token

    def test_a_stroke_is_frozen(self) -> None:
        import pydantic

        stroke = viz_module("style").PropertyStroke(bit=RIVER_BIT, token="river")
        with pytest.raises(pydantic.ValidationError):
            stroke.token = "road"

    def test_an_unknown_field_on_a_stroke_is_refused(self) -> None:
        import pydantic

        with pytest.raises(pydantic.ValidationError):
            viz_module("style").PropertyStroke(bit=RIVER_BIT, token="river", colour="blue")


class TestViewport:
    """World metres to SVG user units: one uniform scale, and a y flip.

    Every comparison is `rel=0.0`. The default relative tolerance is worthless
    on both sides of this transform -- +-0.43 m at an easting of 4.3e5 on the
    world side, several user units on the SVG side -- and 6b-i measured exactly
    that: an off-centre bbox survived a whole mutation round under it.
    """

    def view(self, scene: Any, **overrides: Any) -> Any:
        return viz_module("svg").viewport(scene.bbox, make_style(**overrides))

    def test_the_scale_is_the_smaller_of_the_two_fits(self, scene: Any) -> None:
        # min(952/100, 680/80): the height is the binding constraint, so a
        # transform that took the x fit would overflow the map area vertically.
        assert self.view(scene).scale == pytest.approx(SCALE, rel=0.0, abs=1e-9)

    def test_the_bottom_left_world_corner_maps_to_the_bottom_left(
        self, scene: Any
    ) -> None:
        x, y = self.view(scene).point(EAST, NORTH)
        assert x == pytest.approx(MAP_LEFT, rel=0.0, abs=1e-9)
        assert y == pytest.approx(MAP_BOTTOM, rel=0.0, abs=1e-9)

    def test_the_top_right_world_corner_maps_to_the_top_right(self, scene: Any) -> None:
        x, y = self.view(scene).point(EAST + 100.0, NORTH + 80.0)
        assert x == pytest.approx(MAP_RIGHT, rel=0.0, abs=1e-9)
        assert y == pytest.approx(MAP_TOP, rel=0.0, abs=1e-9)

    def test_the_y_axis_is_flipped(self, scene: Any) -> None:
        # One sign. Without it every picture is mirrored, every winding reads
        # backwards, and the user's intuition is trained on a reflection.
        view = self.view(scene)
        low = view.point(EAST, NORTH)[1]
        high = view.point(EAST, NORTH + 80.0)[1]
        assert high < low

    def test_the_x_axis_is_not_flipped(self, scene: Any) -> None:
        view = self.view(scene)
        assert view.point(EAST, NORTH)[0] < view.point(EAST + 100.0, NORTH)[0]

    def test_aspect_ratio_is_preserved(self, scene: Any) -> None:
        # Equal world lengths, one horizontal and one vertical, must come out
        # equal in user units: that is what "one uniform scale" means, and it is
        # what a picture is judged against when the eye compares a sliver's
        # width to its length.
        view = self.view(scene)
        across = view.point(EAST + 40.0, NORTH)[0] - view.point(EAST, NORTH)[0]
        up = view.point(EAST, NORTH)[1] - view.point(EAST, NORTH + 40.0)[1]
        assert across == pytest.approx(up, rel=0.0, abs=1e-9)
        assert across == pytest.approx(40.0 * SCALE, rel=0.0, abs=1e-9)

    def test_the_drawing_is_centred_in_the_map_area(self, scene: Any) -> None:
        # The bbox centre lands at the map area's centre. A transform that
        # anchors at a corner instead passes every corner test above taken one
        # at a time and puts the letterboxing all on one side.
        x, y = self.view(scene).point(EAST + 50.0, NORTH + 40.0)
        assert x == pytest.approx(CENTRE_X, rel=0.0, abs=1e-9)
        assert y == pytest.approx(CENTRE_Y, rel=0.0, abs=1e-9)

    def test_the_header_band_is_not_drawn_over(self, scene: Any) -> None:
        # The map area starts below the header band, so nothing in the drawing
        # can land in the top `HEADER` units. Checked on the world maximum,
        # which is the only point that can violate it.
        assert self.view(scene).point(EAST, NORTH + 80.0)[1] >= MARGIN + HEADER

    def test_every_scene_vertex_lands_inside_the_margins(self, scene: Any) -> None:
        view = self.view(scene)
        for x_world, y_world in np.asarray(scene.vertices):
            x, y = view.point(float(x_world), float(y_world))
            assert MARGIN - 1e-9 <= x <= WIDTH - MARGIN + 1e-9
            assert MARGIN + HEADER - 1e-9 <= y <= HEIGHT - MARGIN + 1e-9

    def test_a_padded_bbox_does_not_divide_by_zero(self) -> None:
        # 6b-i pads a degenerate axis rather than leaving a zero extent; this is
        # the other half of that contract, and the failure mode without it is a
        # ZeroDivisionError from a viewer, on a fixture a user typed by hand.
        column = np.array(
            [[EAST, NORTH], [EAST, NORTH + 40.0], [EAST, NORTH + 80.0]], dtype=np.float64
        )
        pslg = make_pslg([([0, 1, 2], Role.Breakline, NO_PROPERTIES)], vertices=column)
        degenerate = build_scene(pslg, None, ok=False)
        assert degenerate.bbox.padded is True
        view = viz_module("svg").viewport(degenerate.bbox, make_style())
        assert math.isfinite(view.scale)
        assert view.scale > 0.0

    def test_a_wider_style_rescales(self, scene: Any) -> None:
        # The transform reads the style rather than a constant: doubling the
        # canvas height doubles nothing else, but it does change the fit.
        assert self.view(scene, height=1480).scale > self.view(scene).scale


class TestDocumentShape:
    def test_it_parses_as_xml_in_the_svg_namespace(self, document: ET.Element) -> None:
        # Not a golden file: that the browser gets a well-formed SVG document is
        # the one thing about the text a person cannot check by looking at the
        # picture, because a malformed one shows them nothing at all.
        assert document.tag == tag("svg")

    def test_the_canvas_matches_the_style(self, document: ET.Element) -> None:
        assert numbers(document.get("width") or "") == [float(WIDTH)]
        assert numbers(document.get("height") or "") == [float(HEIGHT)]

    def test_it_declares_a_viewbox_covering_the_canvas(self, document: ET.Element) -> None:
        assert numbers(document.get("viewBox") or "") == [0.0, 0.0, float(WIDTH), float(HEIGHT)]

    def test_the_page_has_a_background_covering_the_canvas(self, document: ET.Element) -> None:
        # "A hole in the mesh is the absence of triangles, so it reads as a hole
        # only if the page reads as a page." The colour is the stylesheet's and
        # is not asserted; that a full-canvas rect exists at all is structure.
        rects = [
            r
            for r in document.iter(tag("rect"))
            if numbers(r.get("width") or "") == [float(WIDTH)]
            and numbers(r.get("height") or "") == [float(HEIGHT)]
        ]
        assert rects != []

    def test_rendering_is_deterministic(self, scene: Any) -> None:
        render = viz_module("svg").render_svg
        assert render(scene, style=make_style()) == render(scene, style=make_style())

    def test_the_style_argument_is_optional(self, scene: Any) -> None:
        assert parse(viz_module("svg").render_svg(scene)).tag == tag("svg")

    def test_nothing_is_written_to_disk(self, scene: Any, tmp_path: Path) -> None:
        # `render_svg` returns a `str`; no file is written below `cli.py`. That
        # is what makes it wrappable in `asyncio.to_thread` unchanged.
        before = set(tmp_path.iterdir())
        assert isinstance(viz_module("svg").render_svg(scene, style=make_style()), str)
        assert set(tmp_path.iterdir()) == before


class TestGeometryElements:
    """Element counts equal the scene's primitive counts, and the coordinates
    are the viewport's -- checked against numbers computed by hand, not against
    the module's own transform.
    """

    def test_one_polygon_per_triangle(self, document: ET.Element) -> None:
        assert len(children(document, "triangles", "polygon")) == TRIANGLE_COUNT

    def test_one_line_per_scene_edge(self, document: ET.Element) -> None:
        # Deduplication is 6b-i's, and this is where it becomes visible: an
        # interior edge emitted twice is drawn at double weight and the picture
        # lies about density.
        assert len(children(document, "edges", "line")) == EDGE_COUNT

    def test_a_triangles_corners_are_the_transformed_world_points(
        self, document: ET.Element
    ) -> None:
        # Triangle 0 is (0, 1, 5) = (0,0), (100,0), (80,40) in local metres.
        # By hand at scale 8.5 about the map centre: (75, 776), (925, 776),
        # (755, 436). Absolute tolerance, for this suite's stated reason.
        corners = numbers(children(document, "triangles", "polygon")[0].get("points") or "")
        expected = [75.0, 776.0, 925.0, 776.0, 755.0, 436.0]
        assert len(corners) == 6
        assert corners == pytest.approx(expected, rel=0.0, abs=1e-6)

    def test_an_edges_endpoints_are_the_transformed_world_points(
        self, document: ET.Element, scene: Any
    ) -> None:
        # Edge (0, 1) is the rectangle's bottom side: (75, 776) -> (925, 776).
        line = edge_element(document, 0, 1, scene)
        assert numbers(line.get("x1") or "") == pytest.approx([MAP_LEFT], rel=0.0, abs=1e-6)
        assert numbers(line.get("y1") or "") == pytest.approx([MAP_BOTTOM], rel=0.0, abs=1e-6)
        assert numbers(line.get("x2") or "") == pytest.approx([MAP_RIGHT], rel=0.0, abs=1e-6)
        assert numbers(line.get("y2") or "") == pytest.approx([MAP_BOTTOM], rel=0.0, abs=1e-6)

    def test_no_vertex_dots_by_default(self, document: ET.Element) -> None:
        dots = group_or_none(document, "vertices")
        assert dots is None or list(dots) == []

    def test_show_vertices_draws_one_dot_per_vertex(self, scene: Any) -> None:
        style = make_style(show_vertices=True)
        document = parse(viz_module("svg").render_svg(scene, style=style))
        assert len(children(document, "vertices", "circle")) == VERTEX_COUNT

    def test_no_labels_by_default(self, document: ET.Element) -> None:
        labels = group_or_none(document, "labels")
        assert labels is None or list(labels) == []

    def test_labels_draw_one_index_per_vertex(self, scene: Any) -> None:
        document = parse(
            viz_module("svg").render_svg(scene, style=make_style(), labels=True)
        )
        texts = children(document, "labels", "text")
        assert len(texts) == VERTEX_COUNT
        assert sorted(text_of(t) for t in texts) == sorted(str(i) for i in range(VERTEX_COUNT))


class TestStrokeClasses:
    """Which stroke class an edge belongs to -- structure, not colour.

    `06-cdt-viewer.md` requires four visually distinct constrained strokes plus
    an unconstrained one, coloured by the role of the input chain. The colours
    live in the stylesheet and are untested; that the *right edges* carry the
    right class is the part a wrong picture would get wrong while still looking
    plausible, and it is the one thing here worth an assertion.
    """

    def test_every_edge_element_is_classified(self, document: ET.Element) -> None:
        for line in children(document, "edges", "line"):
            assert classes(line), "an edge element carries no class attribute"

    def test_constrained_and_unconstrained_are_distinguished(
        self, document: ET.Element, scene: Any
    ) -> None:
        for edge, line in zip(scene.edges, children(document, "edges", "line"), strict=True):
            token = "constrained" if edge.constrained else "unconstrained"
            assert token in classes(line), f"edge {(edge.a, edge.b)} is missing {token!r}"

    def test_a_ring_edge_carries_its_roles_class(
        self, document: ET.Element, scene: Any
    ) -> None:
        assert "role-outer" in classes(edge_element(document, 0, 1, scene))

    def test_a_breakline_carries_its_own_role_class(
        self, document: ET.Element, scene: Any
    ) -> None:
        assert "role-breakline" in classes(edge_element(document, 4, 5, scene))

    def test_a_property_is_a_class_of_its_own(
        self, document: ET.Element, scene: Any
    ) -> None:
        # A river is a breakline plus a property, so it is a fourth stroke
        # rather than a fourth role.
        assert "river" in classes(edge_element(document, 4, 5, scene))

    def test_an_edge_with_no_properties_carries_no_property_token(
        self, document: ET.Element, scene: Any
    ) -> None:
        tokens = classes(edge_element(document, 0, 1, scene))
        assert "river" not in tokens
        assert "road" not in tokens

    def test_an_edge_with_two_properties_draws_exactly_one_stroke(
        self, mesh: FakeMesh
    ) -> None:
        # The finding this ruling rests on: at gallery scale two overlaid
        # strokes on one polyline read as a rendering defect, not as two
        # features. So the SET is carried on the edge, in full, and the
        # DOCUMENT names one token -- precedence is a drawing question, which
        # is why it is answered at the one layer a human is looking at.
        pslg = make_pslg(
            [
                ([0, 1, 2, 3], Role.Outer, NO_PROPERTIES),
                ([4, 5], Role.Breakline, RIVER | ROAD),
            ]
        )
        scene = build_scene(pslg, mesh, closed_roles=CLOSED)
        assert edge_at(scene, 4, 5).properties == RIVER | ROAD

        document = parse(
            viz_module("svg").render_svg(
                scene, style=make_style(property_strokes=strokes(*GALLERY_STROKES))
            )
        )
        line = edge_element(document, 4, 5, scene)
        # WHICH token wins is the next test's subject, not this one's, so the
        # expectation is the gallery list's own head rather than a literal:
        # `GALLERY_STROKES` is declared 800 lines up and an expectation spelled
        # `["river"]` here names a different object than the one the fixture
        # defines -- which is exactly how this test and the next came to
        # contradict each other. The order is pinned in place so that a reorder
        # reads as a fixture change and not as a renderer failure.
        assert GALLERY_STROKES[0] == (ROAD_BIT, "road"), "the gallery is road-first"
        assert property_tokens(line, ("river", "road")) == [GALLERY_STROKES[0][1]]

    def test_the_stroke_drawn_is_the_styles_first_match_not_the_lowest_bit(
        self, mesh: FakeMesh
    ) -> None:
        # Priority is explicit and independent of bit numbering: the same edge,
        # the same two properties, the two precedence lists. A renderer that
        # took the lowest set bit, or the vocabulary's order, draws `river`
        # under both and fails the second.
        pslg = make_pslg(
            [
                ([0, 1, 2, 3], Role.Outer, NO_PROPERTIES),
                ([4, 5], Role.Breakline, RIVER | ROAD),
            ]
        )
        scene = build_scene(pslg, mesh, closed_roles=CLOSED)

        def token_for(*pairs: tuple[int, str]) -> list[str]:
            document = parse(
                viz_module("svg").render_svg(
                    scene, style=make_style(property_strokes=strokes(*pairs))
                )
            )
            return property_tokens(edge_element(document, 4, 5, scene), ("river", "road"))

        assert token_for((RIVER_BIT, "river"), (ROAD_BIT, "road")) == ["river"]
        assert token_for((ROAD_BIT, "road"), (RIVER_BIT, "river")) == ["road"]

    def test_a_property_with_no_stroke_draws_no_token_and_loses_no_edge(
        self, mesh: FakeMesh
    ) -> None:
        # The one asymmetry, and it is deliberate. The DATA path may not lose a
        # property -- `EdgeVocabulary.names` raises on a bit nobody names --
        # while the DRAWING path may decline to draw one: the edge is still
        # drawn as the constrained breakline it is, in its role colour. A lost
        # property is a wrong answer; an undrawn one is a picture with less in
        # it, and those are different failures.
        pslg = make_pslg(
            [
                ([0, 1, 2, 3], Role.Outer, NO_PROPERTIES),
                ([4, 5], Role.Breakline, RIVER),
            ]
        )
        scene = build_scene(pslg, mesh, closed_roles=CLOSED)
        document = parse(
            viz_module("svg").render_svg(
                scene, style=make_style(property_strokes=strokes((ROAD_BIT, "road")))
            )
        )
        line = edge_element(document, 4, 5, scene)
        assert property_tokens(line, ("river", "road")) == []
        assert "role-breakline" in classes(line)
        assert "constrained" in classes(line)

    def test_the_default_style_names_no_property(self, scene: Any) -> None:
        # `SvgStyle.property_strokes` defaults to `()`, so a caller that
        # supplies no precedence gets no property stroke anywhere -- able to
        # fail on its own, and the assertion that keeps a vocabulary out of the
        # module whose docstring says it holds none.
        document = parse(viz_module("svg").render_svg(scene, style=make_style()))
        for line in children(document, "edges", "line"):
            assert "river" not in classes(line)
            assert "road" not in classes(line)

    def test_an_unconstrained_edge_carries_no_role(
        self, document: ET.Element, scene: Any
    ) -> None:
        tokens = classes(edge_element(document, 0, 5, scene))
        assert [t for t in tokens if t.startswith("role-")] == []

    def test_a_string_role_names_its_own_stroke_class(self, mesh: FakeMesh) -> None:
        # `fixtures.py` authors roles as the strings "outer"/"hole"/"breakline",
        # because `viz/` may not name `ChainRole`. The renderer must name a
        # stroke class from a role it cannot introspect:
        # `getattr(role, "name", str(role))` serves the enum and the string
        # alike, and this is the half that only the string exercises.
        pslg = make_pslg([([0, 1, 2, 3], "outer", NO_PROPERTIES), ([4, 5], "breakline", RIVER)])
        scene = build_scene(pslg, mesh, closed_roles=("outer", "hole"))
        document = parse(viz_module("svg").render_svg(scene, style=make_style()))
        assert "role-outer" in classes(edge_element(document, 0, 1, scene))

    def test_a_hole_role_is_a_class_of_its_own(self, mesh: FakeMesh) -> None:
        # Drawn PSLG-only so the hole chain, which no mask agrees with, is not
        # also a finding -- what is under test here is the role, not the join.
        pslg = make_pslg(
            [([0, 1, 2, 3], Role.Outer, NO_PROPERTIES), ([4, 5, 6], Role.Hole, NO_PROPERTIES)],
        )
        scene = build_scene(pslg, None, ok=False, closed_roles=CLOSED)
        document = parse(viz_module("svg").render_svg(scene, style=make_style()))
        assert "role-hole" in classes(edge_element(document, 4, 5, scene))


class TestFindings:
    """Disagreement is drawn, not reconciled (risk 2).

    `Scene.findings` already carries both kinds; the renderer's obligation is
    that the offending edge is visibly marked and that the count reaches the
    header band, so that a non-zero count is a finding about the backend's mask
    at a glance.
    """

    def test_the_fixture_really_disagrees(self, finding_scene: Any) -> None:
        # Able to fail on its own: if the scene ever stops producing a finding,
        # every assertion below goes vacuously green.
        assert len(finding_scene.findings) == 1

    def test_a_finding_edge_is_marked(self, finding_scene: Any) -> None:
        document = parse(viz_module("svg").render_svg(finding_scene, style=make_style()))
        assert "finding" in classes(edge_element(document, 0, 5, finding_scene))

    def test_an_agreeing_edge_is_not_marked(self, finding_scene: Any) -> None:
        document = parse(viz_module("svg").render_svg(finding_scene, style=make_style()))
        assert "finding" not in classes(edge_element(document, 0, 1, finding_scene))

    def test_the_finding_count_reaches_the_header(self, finding_scene: Any) -> None:
        document = parse(viz_module("svg").render_svg(finding_scene, style=make_style()))
        assert mentions_count(text_of(group(document, "header")), "findings", 1)

    def test_a_clean_scene_reports_zero_findings(self, document: ET.Element) -> None:
        # Reported rather than omitted: a header that hides a zero makes the
        # absence of the word indistinguishable from a renderer that forgot it.
        assert mentions_count(text_of(group(document, "header")), "findings", 0)


class TestHeaderBand:
    def test_it_carries_every_count(self, document: ET.Element) -> None:
        header = text_of(group(document, "header"))
        assert mentions_count(header, "triangles", TRIANGLE_COUNT)
        assert mentions_count(header, "vertices", VERTEX_COUNT)
        assert mentions_count(header, "constrained", CONSTRAINED_COUNT)

    def test_it_carries_the_bbox(self, document: ET.Element) -> None:
        header = text_of(group(document, "header"))
        found = numbers(header)
        for bound in (EAST, NORTH, EAST + 100.0, NORTH + 80.0):
            assert any(abs(n - bound) < 1.0 for n in found), f"{bound} is not in the header"

    def test_the_title_is_passed_through(self, scene: Any) -> None:
        document = parse(
            viz_module("svg").render_svg(scene, style=make_style(), title="Gaula, UTM 33N")
        )
        assert "Gaula, UTM 33N" in text_of(group(document, "header"))

    def test_the_status_is_passed_through(self, scene: Any) -> None:
        document = parse(
            viz_module("svg").render_svg(scene, style=make_style(), status="Ok")
        )
        assert "Ok" in text_of(group(document, "header"))

    def test_a_title_with_markup_in_it_cannot_break_the_document(self, scene: Any) -> None:
        # The title is a free string from the CLI boundary. Unescaped, `<` ends
        # the document as far as a browser is concerned and the user gets a
        # blank page from a file that was written successfully.
        hostile = '</svg><script>alert("x")</script> & <'
        document = parse(
            viz_module("svg").render_svg(scene, style=make_style(), title=hostile)
        )
        assert document.tag == tag("svg")
        assert hostile in text_of(group(document, "header"))
        assert list(document.iter(tag("script"))) == []


class TestLegendAndScaleBar:
    def test_the_legend_names_every_role_stroke(self, document: ET.Element) -> None:
        legend = text_of(group(document, "legend")).lower()
        for stroke in ("outer", "hole", "breakline", "unconstrained"):
            assert stroke in legend, f"the legend does not name {stroke!r}"

    def test_the_legend_names_every_property_stroke_the_style_declares(
        self, document: ET.Element
    ) -> None:
        # The legend is derived from `style.property_strokes`, not written into
        # `svg.py`: a hard-coded "river breakline" row is a vocabulary in the
        # module that may not hold one, and it names a stroke the default style
        # can never emit. Both tokens, because a legend that named only the
        # first would pass a one-property gallery.
        legend = text_of(group(document, "legend")).lower()
        for _bit, token in GALLERY_STROKES:
            assert token in legend, f"the legend does not name {token!r}"

    def test_the_legend_names_no_property_the_style_does_not_declare(
        self, scene: Any
    ) -> None:
        # Able to fail on its own, and the assertion that catches the
        # hard-coded row: under the default style there is no property stroke,
        # so there is nothing to put a legend row next to.
        document = parse(viz_module("svg").render_svg(scene, style=make_style()))
        legend = text_of(group(document, "legend")).lower()
        assert "river" not in legend
        assert "road" not in legend

    def test_the_scale_bar_is_labelled_in_world_units(self, document: ET.Element) -> None:
        assert re.search(r"\d+\s*m\b", text_of(group(document, "scale-bar")))

    def test_the_scale_bar_is_as_long_as_it_says_it_is(self, document: ET.Element) -> None:
        # The one assertion that makes the scale bar mean anything: its drawn
        # length must be its stated distance under the same transform the
        # drawing used. 8.5 is hand-computed, not read back from the module, so
        # a bar drawn under some other scale cannot agree with it.
        bar = group(document, "scale-bar")
        label = re.search(r"(\d+(?:\.\d+)?)\s*m\b", text_of(bar))
        assert label is not None
        line = next(iter(bar.iter(tag("line"))))
        drawn = abs(numbers(line.get("x2") or "")[0] - numbers(line.get("x1") or "")[0])
        assert drawn == pytest.approx(float(label.group(1)) * SCALE, rel=0.0, abs=1e-6)

    def test_the_scale_bar_fits_on_the_page(self, document: ET.Element) -> None:
        line = next(iter(group(document, "scale-bar").iter(tag("line"))))
        drawn = abs(numbers(line.get("x2") or "")[0] - numbers(line.get("x1") or "")[0])
        assert 0.0 < drawn <= WIDTH - 2 * MARGIN


class TestFailurePresentation:
    """The renderer never emits a blank page.

    Increment 4's risk 5 is the `Ok`-with-zero-triangles mode: "forgot the
    outline" triangulates successfully to nothing, and in a release build with
    the assert compiled out one test is the only guard. The viewer is the second
    guard, and it only works if the failure is conspicuous rather than blank --
    which is why three of the eight gallery fixtures exist to put this
    presentation in front of a person rather than have it assumed --
    `not-noded`, `degenerate` and, since `4482649`, `hole-in-hole`.
    """

    @pytest.fixture(params=["empty_scene", "failed_scene"])
    def undrawable(self, request: pytest.FixtureRequest) -> Any:
        return request.getfixturevalue(request.param)

    def test_no_triangle_is_drawn(self, undrawable: Any) -> None:
        document = parse(viz_module("svg").render_svg(undrawable, style=make_style()))
        triangles = group_or_none(document, "triangles")
        assert triangles is None or list(triangles) == []

    def test_the_input_pslg_is_drawn_alone(self, undrawable: Any) -> None:
        # Never blank: the user sees exactly the geometry they submitted.
        document = parse(viz_module("svg").render_svg(undrawable, style=make_style()))
        assert len(children(document, "edges", "line")) == len(undrawable.edges)
        assert len(undrawable.edges) == CONSTRAINED_COUNT

    def test_the_roles_still_colour_the_input(self, undrawable: Any) -> None:
        document = parse(viz_module("svg").render_svg(undrawable, style=make_style()))
        assert "role-outer" in classes(edge_element(document, 0, 1, undrawable))

    def test_the_status_and_message_are_shown(self, failed_scene: Any) -> None:
        document = parse(
            viz_module("svg").render_svg(
                failed_scene,
                style=make_style(),
                status="NotNoded",
                message="Point 4 is exactly on a constrained edge",
            )
        )
        header = text_of(group(document, "header"))
        assert "NotNoded" in header
        assert "Point 4 is exactly on a constrained edge" in header
        assert "alarm" in classes(header_line(document, "NotNoded")), "the status is not loud"

    def test_ok_but_empty_says_so_in_those_words(self, empty_scene: Any) -> None:
        # The design's own wording, and the whole point of the mode having its
        # own `SceneKind`: a successful call that produced nothing must not read
        # like a failure and must not read like a success.
        document = parse(
            viz_module("svg").render_svg(empty_scene, style=make_style(), status="Ok")
        )
        assert "Ok BUT EMPTY" in text_of(group(document, "header"))

    def test_a_failure_does_not_claim_to_be_empty_but_ok(self, failed_scene: Any) -> None:
        document = parse(
            viz_module("svg").render_svg(failed_scene, style=make_style(), status="NotNoded")
        )
        assert "Ok BUT EMPTY" not in text_of(group(document, "header"))

    def test_a_padded_bbox_is_declared(self) -> None:
        # "The header band says so, rather than the picture being an accidental
        # point": a padded box means the extent on screen is not the extent of
        # the data, and only the header can say that.
        column = np.array([[EAST, NORTH], [EAST, NORTH + 80.0]], dtype=np.float64)
        pslg = make_pslg([([0, 1], Role.Breakline, NO_PROPERTIES)], vertices=column)
        scene = build_scene(pslg, None, ok=False)
        assert scene.bbox.padded is True
        document = parse(viz_module("svg").render_svg(scene, style=make_style()))
        assert "padded" in text_of(group(document, "header")).lower()
        assert "alarm" in classes(header_line(document, "padded")), "the padding is not loud"


# --- The gallery ------------------------------------------------------------


def cross(o: Sequence[float], a: Sequence[float], b: Sequence[float]) -> float:
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def sine_of(o: Sequence[float], a: Sequence[float], b: Sequence[float]) -> float:
    """|sin| of the angle at `o`, which is the cross product normalised.

    Scaled rather than absolute, and that is deliberate: collinearity is a
    statement about direction, so at UTM magnitudes an absolute bound on an
    unnormalised cross product would call two segments a kilometre apart
    parallel. This is the one place in the suite where a relative tolerance is
    the correct instrument rather than the trap.
    """
    u = (a[0] - o[0], a[1] - o[1])
    v = (b[0] - o[0], b[1] - o[1])
    lengths = math.hypot(*u) * math.hypot(*v)
    return abs(cross(o, a, b)) / lengths if lengths else 0.0


def properly_cross(
    p: Sequence[float], q: Sequence[float], r: Sequence[float], s: Sequence[float]
) -> bool:
    """Do segments `pq` and `rs` cross at an interior point of both?

    Strict orientation signs on both sides, so a shared endpoint or a T-junction
    is not a crossing -- the `not-noded` fixture has to be a *crossing*, since a
    T-junction is a different unnoded defect and triangulates differently.
    """
    def sign(value: float) -> int:
        return (value > 0.0) - (value < 0.0)

    # Strict on both sides, and `sign` rather than `> 0` is what makes it
    # strict: with a boolean test a zero determinant reads as "negative", so a
    # shared endpoint -- determinant exactly zero -- is reported as a crossing.
    # Measured: the self-test below failed under `(d1 > 0) != (d2 > 0)`, which
    # would have let `corner-hole`'s touching rings satisfy `not-noded`'s test.
    d1, d2 = sign(cross(r, s, p)), sign(cross(r, s, q))
    d3, d4 = sign(cross(p, q, r)), sign(cross(p, q, s))
    return d1 * d2 < 0 and d3 * d4 < 0


class TestGalleryPredicates:
    """The two geometric oracles, self-tested.

    Neither is borrowed from the code under test -- `fixtures.py` is data and
    computes nothing -- but a predicate that cannot fail measures nothing, so
    each is shown to separate a positive case from a negative one before it is
    used on a fixture.
    """

    def test_the_collinearity_predicate_separates_the_two_cases(self) -> None:
        on = sine_of((EAST, NORTH), (EAST + 50.0, NORTH), (EAST + 90.0, NORTH))
        off = sine_of((EAST, NORTH), (EAST + 50.0, NORTH), (EAST + 90.0, NORTH + 40.0))
        assert on < 1.0e-9
        assert off > 1.0e-2

    def test_the_crossing_predicate_separates_the_two_cases(self) -> None:
        a, b = (EAST, NORTH), (EAST + 100.0, NORTH + 100.0)
        c, d = (EAST, NORTH + 100.0), (EAST + 100.0, NORTH)
        assert properly_cross(a, b, c, d)
        assert not properly_cross(a, b, (EAST + 200.0, NORTH), (EAST + 200.0, NORTH + 100.0))

    def test_a_shared_endpoint_is_not_a_crossing(self) -> None:
        # What separates `not-noded` from `corner-hole`: one is two constraints
        # through each other, the other is two rings meeting at a vertex.
        a, b = (EAST, NORTH), (EAST + 100.0, NORTH)
        assert not properly_cross(a, b, b, (EAST + 100.0, NORTH + 100.0))


def gallery() -> Any:
    return viz_module("fixtures").GALLERY


def chain_points(fixture: Any, c: int) -> list[tuple[float, float]]:
    vertices = np.asarray(fixture.vertices, dtype=np.float64)
    return [(float(vertices[i][0]), float(vertices[i][1])) for i in fixture.indices_of(c)]


def chains_with(fixture: Any, role: str) -> list[int]:
    return [c for c, chain in enumerate(fixture.chains) if chain.role == role]


def bounds(points: Sequence[Sequence[float]]) -> tuple[float, float, float, float]:
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return min(xs), min(ys), max(xs), max(ys)


class TestGallery:
    """Eight fixtures, each answering one question a person can ask of a
    picture. Declarative data -- no procedural generation, because a fixture
    whose coordinates are computed is a fixture nobody can check by reading.
    """

    def test_it_holds_exactly_the_eight_named_fixtures(self) -> None:
        assert sorted(gallery()) == sorted(GALLERY_NAMES)

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_each_fixture_describes_itself(self, name: str) -> None:
        # `--list` prints these, and a fixture whose description is empty is a
        # fixture whose reason for existing has been lost.
        fixture = gallery()[name]
        assert fixture.name == name
        assert fixture.description.strip()

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_each_fixture_is_a_pslg_like(self, name: str) -> None:
        # Structural, because `cli.py` hands the fixture itself to `build_scene`
        # -- which is what lets a fixture the PSLG validator REJECTS still be
        # drawn. See `test_cli_draw.py::TestFailurePresentation`.
        fixture = gallery()[name]
        for member in ("vertices", "chains", "chain_indices", "indices_of"):
            assert hasattr(fixture, member), f"{name} has no {member}"
        for chain in fixture.chains:
            for member in ("begin", "count", "role", "properties"):
                assert hasattr(chain, member), f"{name}'s chain has no {member}"

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_each_vertex_array_is_n_by_two_float64(self, name: str) -> None:
        vertices = np.asarray(gallery()[name].vertices)
        assert vertices.dtype == np.float64
        assert vertices.ndim == 2
        assert vertices.shape[1] == 2
        assert vertices.shape[0] >= 3

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_every_coordinate_is_finite(self, name: str) -> None:
        assert np.isfinite(np.asarray(gallery()[name].vertices)).all()

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_coordinates_are_at_utm_33n_magnitudes(self, name: str) -> None:
        # The ordering ruling's one mitigation: the robustness behaviour a
        # person is asked to judge is magnitude-sensitive, and a picture drawn
        # in the unit square silently exercises the easy case.
        vertices = np.asarray(gallery()[name].vertices)
        assert vertices[:, 0].min() >= 1.0e5
        assert vertices[:, 1].min() >= 1.0e6

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_the_domain_is_small_enough_to_hold_in_the_head(self, name: str) -> None:
        # A fixture is a shape a person compares the picture against, so it is
        # catchment-sized: hundreds of metres, not hundreds of kilometres.
        vertices = np.asarray(gallery()[name].vertices)
        assert (vertices.max(axis=0) - vertices.min(axis=0)).max() <= 1.0e4

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_every_chain_index_is_in_range(self, name: str) -> None:
        fixture = gallery()[name]
        count = len(np.asarray(fixture.vertices))
        assert fixture.chains
        for c in range(len(fixture.chains)):
            indices = [int(i) for i in fixture.indices_of(c)]
            assert indices, f"{name}'s chain {c} is empty"
            assert all(0 <= i < count for i in indices), f"{name}'s chain {c} is out of range"

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_indices_of_partitions_the_flat_buffer(self, name: str) -> None:
        fixture = gallery()[name]
        recovered = np.concatenate(
            [np.asarray(fixture.indices_of(c)) for c in range(len(fixture.chains))]
        )
        assert np.array_equal(recovered, np.asarray(fixture.chain_indices))

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_every_role_is_one_the_cli_can_map(self, name: str) -> None:
        # `cli.py` maps these strings onto `_core.ChainRole`; a role outside the
        # vocabulary is a KeyError at the composition root, on a fixture.
        roles = {chain.role for chain in gallery()[name].chains}
        assert roles <= ROLE_NAMES, f"{name} carries roles outside {sorted(ROLE_NAMES)}"

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_every_fixture_declares_an_outer_ring(self, name: str) -> None:
        # No outer chain is `PslgError.NoOuterChain`, which is a rejection at
        # validation and not the triangulation failure any fixture is there to
        # show -- `degenerate` included.
        assert chains_with(gallery()[name], "outer"), f"{name} has no outer chain"

    def test_the_catchment_is_the_shape_the_project_is_for(self) -> None:
        fixture = gallery()["catchment"]
        assert len(chains_with(fixture, "hole")) == 2, "two lake holes"
        assert chains_with(fixture, "breakline"), "a braided breakline"

    @pytest.mark.parametrize("name", ["catchment", "hole-in-hole"])
    def test_each_hole_sits_strictly_inside_an_outer_ring(self, name: str) -> None:
        # Bounding-box containment rather than point-in-polygon: weaker, but it
        # is checkable by reading the coordinates, which is the same standard
        # the fixtures themselves are held to. `corner-hole` is deliberately
        # excluded -- its hole touches the ring, which is its whole point.
        fixture = gallery()[name]
        outer = bounds(chain_points(fixture, chains_with(fixture, "outer")[0]))
        for c in chains_with(fixture, "hole"):
            hole = bounds(chain_points(fixture, c))
            assert outer[0] < hole[0] and outer[1] < hole[1]
            assert hole[2] < outer[2] and hole[3] < outer[3]

    def test_hole_in_hole_really_nests_twice(self) -> None:
        # "A hole nested inside a second outline": two holes, the inner one
        # strictly inside the outer one, which is a different picture from the
        # catchment's two disjoint lakes.
        fixture = gallery()["hole-in-hole"]
        holes = chains_with(fixture, "hole")
        assert len(holes) >= 2
        boxes = sorted(
            (bounds(chain_points(fixture, c)) for c in holes),
            key=lambda b: (b[2] - b[0]) * (b[3] - b[1]),
        )
        inner, outer = boxes[0], boxes[-1]
        assert outer[0] < inner[0] and outer[1] < inner[1]
        assert inner[2] < outer[2] and inner[3] < outer[3]

    def test_the_corner_hole_touches_the_ring_at_exactly_one_vertex(self) -> None:
        # `InvalidTopology`'s neighbour: one shared vertex is a topology worth
        # looking at, two would be a slit and a different fixture.
        fixture = gallery()["corner-hole"]
        outer = {int(i) for i in fixture.indices_of(chains_with(fixture, "outer")[0])}
        hole = {int(i) for i in fixture.indices_of(chains_with(fixture, "hole")[0])}
        assert len(outer & hole) == 1

    def test_the_breakline_chain_is_open(self) -> None:
        # An open chain crossing the interior: if its ends coincided it would be
        # a ring, and "where the Delaunay property visibly stops" would not be
        # what the picture showed.
        fixture = gallery()["breakline-chain"]
        chain = chains_with(fixture, "breakline")[0]
        points = chain_points(fixture, chain)
        assert len(points) >= 3
        assert points[0] != points[-1]

    def test_the_river_fixture_sets_the_river_bit_and_only_that_bit(self) -> None:
        # The property stroke is exercised here so that it is known to work
        # before 5b depends on it.
        #
        # `== RIVER`, not "is non-empty": the gallery holds a BARE mask --
        # `fixtures.py` imports nothing first-party, so it cannot ask a
        # vocabulary -- and the number it writes down has to be the number
        # `DEFAULT_VOCABULARY` gives bit 0 to `river`. That agreement is risk
        # 20 in its smallest form, and this is the only place it is checkable.
        masks = [int(chain.properties) for chain in gallery()["river"].chains]
        assert RIVER in masks
        assert set(masks) <= {NO_PROPERTIES, RIVER}

    def test_only_the_river_fixture_carries_a_property(self) -> None:
        # Able to fail on its own: if every fixture set a bit, the test above
        # would pass without the `river` fixture existing at all.
        classified = {
            name
            for name in GALLERY_NAMES
            if any(chain.properties for chain in gallery()[name].chains)
        }
        assert classified == {"river"}

    def test_the_sliver_fan_is_near_collinear_without_being_collinear(self) -> None:
        # The aspect-ratio question: near-collinear is the hard case for the
        # triangulator, exactly collinear is a different fixture (`degenerate`).
        fixture = gallery()["sliver-fan"]
        sines = [
            sine_of(*chain_points(fixture, c)[i : i + 3])
            for c in range(len(fixture.chains))
            for i in range(max(0, len(chain_points(fixture, c)) - 2))
        ]
        tight = [s for s in sines if 0.0 < s < 1.0e-2]
        assert tight, "no near-collinear triple: the fan has no sliver in it"

    def test_the_degenerate_fixture_is_exactly_collinear(self) -> None:
        points = [
            (float(x), float(y)) for x, y in np.asarray(gallery()["degenerate"].vertices)
        ]
        assert len(points) >= 3
        assert all(sine_of(points[0], points[1], p) < 1.0e-9 for p in points[2:])

    def test_the_not_noded_fixture_really_crosses(self) -> None:
        # A deliberate failure fixture. If its constraints stopped crossing it
        # would triangulate cleanly and the non-`Ok` presentation nobody has
        # looked at would go back to being assumed rather than seen.
        fixture = gallery()["not-noded"]
        segments = [
            (points[i], points[i + 1])
            for c in range(len(fixture.chains))
            for points in [chain_points(fixture, c)]
            for i in range(len(points) - 1)
        ]
        crossings = [
            (i, j)
            for i, (p, q) in enumerate(segments)
            for j, (r, s) in enumerate(segments)
            if i < j and properly_cross(p, q, r, s)
        ]
        assert crossings, "no two constraints cross"

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_every_fixture_renders_without_a_mesh(self, name: str) -> None:
        # The PSLG-only picture, which is what every fixture falls back to when
        # the backend refuses it, and which must therefore never blank.
        fixture = gallery()[name]
        scene = build_scene(fixture, None, ok=False, closed_roles=("outer", "hole"))
        document = parse(viz_module("svg").render_svg(scene, style=make_style()))
        assert document.tag == tag("svg")
        assert children(document, "edges", "line") != []
