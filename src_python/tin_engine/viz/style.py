"""The renderer's one input that is not geometry: the page it draws on.

``SvgStyle`` is deliberately thin. Colours, stroke widths and dash arrays live
in ``svg.py``'s stylesheet as CSS, because a class token on an element is
structure -- which edge belongs to which stroke class -- while the colour that
token resolves to is taste, and ``06-cdt-viewer.md`` rules that the stylesheet
is untested. What is modelled here is the arithmetic the viewport divides by,
which is the part a wrong value silently corrupts.

Frozen and ``extra="forbid"``: a style is passed down a call chain and shared,
so mutating one would change a picture already described, and a typo in an
override is a silent no-op otherwise -- the symptom being a picture that looks
almost right.

This module imports nothing first-party. ``SvgStyle`` is the renderer's input
and knows nothing about a scene, which is why ``06-cdt-viewer.md`` pairs it with
6b-ii rather than with the scene builder.
"""

from __future__ import annotations

from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field


class PropertyStroke(BaseModel):
    """One property bit, and the class token an edge carrying it is drawn with.

    Structure rather than taste, which is why it lives here and not in the
    stylesheet: WHICH of an edge's properties gets a token is a decision about
    the document, and the colour that token resolves to is the CSS's.

    ``bit`` is a position, never a name -- this module holds no vocabulary, so
    the caller that knows bit 0 means "river" is the composition root. 32 is the
    C++ ceiling (``EdgeProperties::kMaxProperties``), restated at the layer that
    meets untrusted input: a stroke on bit 32 is a stroke no edge can carry, so
    it is a typo and is caught at construction.

    ``token``'s pattern is a security boundary, not a style preference.
    ``svg.py``'s ``_edge_classes`` interpolates it into a ``class="..."``
    attribute unescaped -- only ``_text`` escapes -- so a token containing a
    quote read from a configuration file would end the attribute. A hyphen is
    legal here and deliberately is not in ``features.EdgeProperty.name``: a CSS
    class may contain one and a Python-side feature name may not.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    bit: Annotated[int, Field(ge=0, lt=32)]
    token: Annotated[str, Field(pattern=r"^[a-z][a-z0-9_-]*$")]


class SvgStyle(BaseModel):
    """Canvas geometry, in SVG user units, plus the two optional overlays."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    width: Annotated[int, Field(gt=0)] = 1000
    """Canvas width. Positive, because the map area's width divides a scale."""

    height: Annotated[int, Field(gt=0)] = 800
    """Canvas height. Positive, for the same reason."""

    margin: Annotated[float, Field(ge=0.0)] = 24.0
    """Blank border on all four sides. Zero is legal; negative is not."""

    header_height: Annotated[float, Field(ge=0.0)] = 72.0
    """Band reserved at the top for status, counts and title. The map area
    starts below it, so nothing in the drawing can land on top of it."""

    show_vertices: bool = False
    """Draw a dot per vertex. Off by default: on a dense mesh the dots merge
    into a grey wash and hide the triangulation they sit on."""

    property_strokes: tuple[PropertyStroke, ...] = ()
    """Which property an edge carrying several is drawn as, in DESCENDING draw
    priority: first match wins, and a polyline carries at most one of them,
    because two overlaid strokes on one line read as a rendering defect rather
    than as two features.

    The order is the priority, and the model does not sort -- sorting by bit
    would silently replace this ruling with however a vocabulary chose to number
    its features. Empty by default: a default naming ``river`` would be policy
    in the module whose own docstring says the vocabulary is somebody else's, so
    the list is ``cli.py``'s to supply, exactly as ``closed_roles`` is."""
