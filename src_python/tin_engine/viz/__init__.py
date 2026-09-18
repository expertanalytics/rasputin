"""Rendering a CDT as a picture a person can form an opinion about.

Increment 6a ships only ``protocols``; the scene builder, the SVG writer, the
style model and the fixture gallery are 6b. Nothing in this package imports
``tin_engine._core``: the renderer is written against the protocols below, so
it is testable against a hand-built fake with no compiled extension in the
process, and ``cli.py`` is the single composition root that joins the two.
"""
