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
