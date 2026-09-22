"""``(Scene, SvgStyle) -> str``: the picture, as text.

Hand-written SVG, per ``06-cdt-viewer.md``'s dependency ruling -- no plotting
library, nothing beyond the standard library and NumPy. The module is a pure
function of its arguments: it opens no file, resolves no path and returns a
``str``, which is what lets a future GUI wrap it in ``asyncio.to_thread``
unchanged and what keeps every path in this project on ``cli.py``'s side.

Two conventions are worth stating before reading the emitters.

**Class tokens are structure; colour is not.** Every edge element carries the
stroke class it belongs to -- ``constrained``/``unconstrained``, ``role-outer``
and friends, ``finding``, and at most one property token -- and
:data:`STYLESHEET` turns those tokens into strokes. A property token this module
cannot name arrives from ``SvgStyle.property_strokes``, which is how the
precedence stays the composition root's decision and this module stays free of a
vocabulary; a token with no CSS rule simply resolves to the stroke it already
had. The tokens are what the suite asserts, because "the right edges
carry the right class" is the part a wrong picture gets wrong while still
looking plausible; the colours are taste and are deliberately untested.

**Roles are opaque.** ``viz/`` may not name ``_core.ChainRole``, so a role
arrives either as that enum (from ``cli.py``) or as one of the strings
``fixtures.py`` authors. ``getattr(role, "name", role)`` names both, and that
one expression is the entire concession this module makes to not knowing the
type.

Edges are emitted in ``Scene.edges`` order, which is sorted and therefore
interleaves constrained and unconstrained strokes rather than painting the
constrained ones last. They are distinguished by width and opacity in the
stylesheet instead of by paint order.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from xml.sax.saxutils import escape

import numpy as np

from .scene import BBox, Scene, SceneEdge, SceneKind
from .style import PropertyStroke, SvgStyle

SVG_NS = "http://www.w3.org/2000/svg"

#: Untested by design (`06-cdt-viewer.md`: a golden-file comparison is rejected
#: outright and no assertion is made on any style string). The page is a light
#: neutral that is *not* white, which is load-bearing for comprehension rather
#: than taste: a hole in the mesh is the absence of triangles, so it reads as a
#: hole only if the page reads as a page.
STYLESHEET = """
.page { fill: #eceae4; }
.triangle { fill: #f7f5ef; fill-opacity: 0.75; stroke: none; }
line.unconstrained { stroke: #9aa0a6; stroke-width: 0.6; }
line.constrained { stroke: #444444; stroke-width: 2.0; }
line.role-outer { stroke: #1b5e20; }
line.role-hole { stroke: #b06000; }
line.role-breakline { stroke: #1a4f8a; }
line.river { stroke: #0b8fb0; stroke-width: 2.6; }
line.road { stroke: #8a5a00; stroke-width: 2.2; }
line.finding { stroke: #c0132b; stroke-width: 3.2; stroke-dasharray: 6 3; }
circle.vertex { fill: #444444; }
text { font-family: sans-serif; font-size: 11px; fill: #23282d; }
text.status { font-size: 18px; font-weight: bold; }
text.alarm { fill: #c0132b; }
"""

LINE_HEIGHT = 13.0
VERTEX_RADIUS = 2.0


@dataclass(frozen=True, slots=True)
class Viewport:
    """World metres to SVG user units: one uniform scale, and a y flip.

    A single ``scale`` for both axes is what makes the picture worth judging by
    eye -- the eye compares a sliver's width to its length, and two scales would
    make that comparison a lie. The flip is one sign: world y increases upward
    and SVG y increases downward, so without it every picture is mirrored and
    every winding reads backwards.
    """

    scale: float
    world_x: float
    world_y: float
    svg_x: float
    svg_y: float

    def point(self, x: float, y: float) -> tuple[float, float]:
        """Map one world point into user units, about the two centres."""
        return (
            self.svg_x + (x - self.world_x) * self.scale,
            self.svg_y - (y - self.world_y) * self.scale,
        )


def viewport(bbox: BBox, style: SvgStyle) -> Viewport:
    """Fit ``bbox`` into the map area, centred, at one scale.

    The map area is the canvas less the margins and the header band, so nothing
    drawn can land on the band. The scale is the *smaller* of the two fits: the
    larger one overflows the area on the other axis.

    ``bbox`` is a precondition, not an input to check. ``build_scene`` has
    already refused any non-finite coordinate and padded a zero-extent axis, so
    neither division here can fail -- and re-checking would be the viewport
    asserting about a box that is already meaningless if it were wrong.
    """
    span_x = bbox.max_x - bbox.min_x
    span_y = bbox.max_y - bbox.min_y
    area_width = style.width - 2.0 * style.margin
    area_height = style.height - 2.0 * style.margin - style.header_height
    return Viewport(
        scale=min(area_width / span_x, area_height / span_y),
        world_x=(bbox.min_x + bbox.max_x) / 2.0,
        world_y=(bbox.min_y + bbox.max_y) / 2.0,
        svg_x=style.width / 2.0,
        svg_y=style.margin + style.header_height + area_height / 2.0,
    )


def _num(value: float) -> str:
    """A coordinate, to the micro-unit, with no trailing noise.

    Fixed point rather than ``repr``: it keeps the document diffable and the
    same scene renders byte-identically twice, which is the only thing
    ``test_rendering_is_deterministic`` can check without pinning a style.
    """
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return "0" if text in {"", "-", "-0"} else text


def _role_name(role: object) -> str:
    """The stroke-class stem of an opaque role: an enum's name, or the string."""
    return str(getattr(role, "name", role)).lower()


def _group(gid: str, rows: list[str]) -> str:
    """One ``<g id=...>``, always emitted -- empty when it has nothing to draw.

    Emitting the empty group rather than omitting it keeps "this picture has no
    triangles" and "this renderer forgot triangles" the same shape in the
    document, which is what the failure presentation is judged on.
    """
    return "\n".join([f'<g id="{gid}">', *rows, "</g>"])


def _property_token(properties: int, strokes: tuple[PropertyStroke, ...]) -> str | None:
    """The one token an edge's property set is drawn with, or ``None``.

    First match over ``strokes``, which is in descending draw priority: an edge
    carries a SET and a polyline carries one stroke, because two overlaid
    strokes on one line read as a rendering defect rather than as two features.
    First match rather than the lowest set bit, because priority is the
    stylesheet's ruling and a bit position is a vocabulary's numbering.

    A property with no stroke contributes nothing and the edge is still drawn as
    the constraint it is. That is the deliberate asymmetry: the DATA path may not
    lose a property -- ``features.EdgeVocabulary.names`` raises on a bit nobody
    names -- while the DRAWING path may decline to draw one. A lost property is a
    wrong answer; an undrawn one is a picture with less in it.
    """
    return next((s.token for s in strokes if properties >> s.bit & 1), None)


def _edge_classes(
    edge: SceneEdge, findings: frozenset[tuple[int, int]], style: SvgStyle
) -> str:
    """Every stroke class this edge belongs to, space separated.

    ``constrained`` is the mesh mask's verdict and ``role-*`` the input chains';
    they are separate tokens because they are separate derivations, and
    ``finding`` marks the edges where the two disagree (risk 2). A property is a
    role's stroke plus one token, so it is a fourth class rather than a fourth
    role.
    """
    tokens = ["constrained" if edge.constrained else "unconstrained"]
    if edge.role is not None:
        tokens.append(f"role-{_role_name(edge.role)}")
    token = _property_token(edge.properties, style.property_strokes)
    if token is not None:
        tokens.append(token)
    if (edge.a, edge.b) in findings:
        tokens.append("finding")
    return " ".join(tokens)


def _text(x: float, y: float, content: str, css: str = "") -> str:
    """One ``<text>``, with its content escaped.

    The title and the backend's message are free strings from outside. Left
    unescaped, a ``<`` in either ends the document as far as a browser is
    concerned and the user gets a blank page from a file that was written
    successfully -- a failure mode indistinguishable from a renderer bug.
    """
    css_attr = f' class="{css}"' if css else ""
    return f'<text x="{_num(x)}" y="{_num(y)}"{css_attr}>{escape(content)}</text>'


def _triangles(scene: Scene, view: Viewport) -> str:
    """One ``<polygon>`` per triangle, in the scene's order."""
    vertices = np.asarray(scene.vertices, dtype=np.float64)
    rows: list[str] = []
    for triangle in np.asarray(scene.triangles, dtype=np.uint32):
        corners = [view.point(float(vertices[i][0]), float(vertices[i][1])) for i in triangle]
        points = " ".join(f"{_num(x)},{_num(y)}" for x, y in corners)
        rows.append(f'<polygon class="triangle" points="{points}"/>')
    return _group("triangles", rows)


def _edges(scene: Scene, view: Viewport, style: SvgStyle) -> str:
    """One ``<line>`` per scene edge -- one, because 6b-i deduplicated them.

    An interior edge emitted twice is drawn at double weight and the picture
    lies about the mesh's density, which is precisely the kind of wrong picture
    that does not look wrong.
    """
    vertices = np.asarray(scene.vertices, dtype=np.float64)
    findings = frozenset((f.a, f.b) for f in scene.findings)
    rows: list[str] = []
    for edge in scene.edges:
        x1, y1 = view.point(float(vertices[edge.a][0]), float(vertices[edge.a][1]))
        x2, y2 = view.point(float(vertices[edge.b][0]), float(vertices[edge.b][1]))
        rows.append(
            f'<line class="{_edge_classes(edge, findings, style)}"'
            f' x1="{_num(x1)}" y1="{_num(y1)}"'
            f' x2="{_num(x2)}" y2="{_num(y2)}"/>'
        )
    return _group("edges", rows)


def _vertices(scene: Scene, view: Viewport, style: SvgStyle) -> str:
    """One dot per vertex, off unless the style asks."""
    rows: list[str] = []
    if style.show_vertices:
        for x_world, y_world in np.asarray(scene.vertices, dtype=np.float64):
            x, y = view.point(float(x_world), float(y_world))
            rows.append(
                f'<circle class="vertex" cx="{_num(x)}" cy="{_num(y)}" r="{VERTEX_RADIUS}"/>'
            )
    return _group("vertices", rows)


def _labels(scene: Scene, view: Viewport, labels: bool) -> str:
    """One index per vertex. Bounded at the CLI, not here: the refusal needs a
    message and an exit code, and neither belongs to a pure function."""
    rows: list[str] = []
    if labels:
        for index, (x_world, y_world) in enumerate(np.asarray(scene.vertices, dtype=np.float64)):
            x, y = view.point(float(x_world), float(y_world))
            rows.append(_text(x + 2.0 * VERTEX_RADIUS, y - VERTEX_RADIUS, str(index)))
    return _group("labels", rows)


def _status_line(scene: Scene, status: str) -> str:
    """What the band says happened, in the design's own words.

    ``Ok BUT EMPTY`` is a successful call that produced nothing -- increment 4's
    risk 5, the silent mode -- and it must read neither like a failure nor like
    a success, which is why it gets a wording of its own rather than a status
    name alone.
    """
    if scene.kind is SceneKind.OK_BUT_EMPTY:
        return f"{status or 'Ok'} BUT EMPTY"
    return status


def _header(scene: Scene, style: SvgStyle, title: str, status: str, message: str) -> str:
    """Status, message, the four counts, the bbox and the caller's title.

    The counts are the picture's own arithmetic written down, so that a person
    can check what they are looking at against what the engine says it is. A
    zero finding count is *printed* rather than omitted: an absent word cannot
    be told apart from a renderer that forgot it.
    """
    box = scene.bbox
    constrained = sum(1 for edge in scene.edges if edge.constrained)
    alarming = scene.kind is not SceneKind.MESH
    lines = [
        (_status_line(scene, status), "status alarm" if alarming else "status"),
        (message, ""),
        (title, ""),
        (
            f"Vertices: {len(scene.vertices)}  Triangles: {len(scene.triangles)}  "
            f"Constrained: {constrained}  Findings: {len(scene.findings)}",
            "",
        ),
        (
            f"bbox [{_num(box.min_x)}, {_num(box.min_y)}] to "
            f"[{_num(box.max_x)}, {_num(box.max_y)}]" + (" (padded)" if box.padded else ""),
            "alarm" if box.padded else "",
        ),
    ]
    rows = [
        _text(style.margin, style.margin + LINE_HEIGHT * (row + 1), content, css)
        for row, (content, css) in enumerate(line for line in lines if line[0])
    ]
    return _group("header", rows)


def _legend(style: SvgStyle) -> str:
    """The stroke classes, named. Without it the colours mean nothing.

    The property rows are DERIVED from ``style.property_strokes`` rather than
    written down here: a hard-coded ``river`` row would be a vocabulary in the
    module that holds none, and it would name a stroke the default style can
    never emit -- a legend entry for a class no edge in the document carries.
    """
    entries = [
        ("role-outer constrained", "outer ring"),
        ("role-hole constrained", "hole ring"),
        ("role-breakline constrained", "breakline"),
        *(
            (f"constrained {stroke.token}", stroke.token)
            for stroke in style.property_strokes
        ),
        ("unconstrained", "unconstrained edge"),
    ]
    rows: list[str] = []
    y = style.height - style.margin - LINE_HEIGHT * len(entries)
    for css, label in entries:
        rows.append(
            f'<line class="{css}" x1="{_num(style.margin)}" y1="{_num(y)}"'
            f' x2="{_num(style.margin + 24.0)}" y2="{_num(y)}"/>'
        )
        rows.append(_text(style.margin + 30.0, y + 4.0, label))
        y += LINE_HEIGHT
    return _group("legend", rows)


def _nice_length(span: float) -> tuple[float, str]:
    """A round world distance of about a quarter of the drawing, and its label.

    1, 2 or 5 times a power of ten, so the bar reads as a measurement rather
    than as an artefact of the bbox. The label is formatted at the decade's own
    precision rather than with ``%g``, because ``%g`` reaches for an exponent on
    a small enough domain and ``1e-05 m`` is a label no reader parses the way
    the bar was drawn.
    """
    target = span / 4.0
    decade = math.floor(math.log10(target))
    base = 10.0**decade
    length = next((step * base for step in (5.0, 2.0, 1.0) if step * base <= target), base)
    decimals = max(0, -math.floor(math.log10(length)))
    return length, f"{length:.{decimals}f}"


def _scale_bar(scene: Scene, view: Viewport, style: SvgStyle) -> str:
    """A bar that is as long as it says it is, under the drawing's own scale.

    That equality is the only thing that makes the bar mean anything: a bar
    drawn at some other scale is a ruler that lies, and a picture the user is
    asked to form an opinion about with a lying ruler on it is worse than one
    with no ruler at all.
    """
    length, label = _nice_length(scene.bbox.max_x - scene.bbox.min_x)
    right = style.width - style.margin
    left = right - length * view.scale
    y = style.height - style.margin - LINE_HEIGHT
    return _group(
        "scale-bar",
        [
            f'<line class="constrained" x1="{_num(left)}" y1="{_num(y)}"'
            f' x2="{_num(right)}" y2="{_num(y)}"/>',
            _text(left, y + LINE_HEIGHT, f"{label} m"),
        ],
    )


def render_svg(
    scene: Scene,
    style: SvgStyle | None = None,
    *,
    title: str = "",
    status: str = "",
    message: str = "",
    labels: bool = False,
) -> str:
    """Draw a scene. Returns the document; writes nothing.

    ``status`` and ``message`` are passthrough text. ``viz/`` may not name
    ``_core.CdtStatus``, so ``cli.py`` composes the status name and its prose
    and hands the result over as strings -- the same reasoning that makes
    ``build_scene``'s verdict a bare ``ok`` flag.

    The page is never blank. When there is no drawable mesh the input PSLG is
    drawn alone in its role colours, with the engine's own words above it, so
    that a failure is conspicuous rather than absent.
    """
    style = SvgStyle() if style is None else style
    view = viewport(scene.bbox, style)
    return (
        "\n".join(
            [
                f'<svg xmlns="{SVG_NS}" width="{style.width}" height="{style.height}"'
                f' viewBox="0 0 {style.width} {style.height}">',
                f"<style>{STYLESHEET}</style>",
                f'<rect class="page" x="0" y="0" width="{style.width}" height="{style.height}"/>',
                _triangles(scene, view),
                _edges(scene, view, style),
                _vertices(scene, view, style),
                _labels(scene, view, labels),
                _header(scene, style, title, status, message),
                _legend(style),
                _scale_bar(scene, view, style),
                "</svg>",
            ]
        )
        + "\n"
    )
