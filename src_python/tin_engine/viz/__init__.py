"""Rendering a CDT as a picture a person can form an opinion about.

``protocols`` describes the two core types drawn, ``scene`` maps a mesh and
its input onto what is drawn, ``svg`` writes the picture, ``style`` sizes the
page, and ``fixtures`` is a gallery of shapes to draw. Nothing in this package imports
``tin_engine._core``: the renderer is written against the protocols below, so
it is testable against a hand-built fake with no compiled extension in the
process, and ``cli.py`` is the single composition root that joins the two.
"""

from .scene import Scene, build_scene
from .style import SvgStyle
from .svg import render_svg

__all__ = ["Scene", "SvgStyle", "build_scene", "render_svg"]
