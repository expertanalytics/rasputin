"""Rasputin terrain meshing engine."""

import importlib.metadata

from tin_engine._core import Point2, Point3, cross, dot

__all__ = ["Point2", "Point3", "cross", "dot", "installed_version"]


def installed_version() -> str:
    """The installed rasputin version, or a placeholder in a source checkout.

    Without `pip install` there is no rasputin metadata; `rasputin version`
    prints this placeholder and a fetch records it in the cache manifest.
    """
    try:
        return importlib.metadata.version("rasputin")
    except importlib.metadata.PackageNotFoundError:
        return "unknown (rasputin is not installed)"
