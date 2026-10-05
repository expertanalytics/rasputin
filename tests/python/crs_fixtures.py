"""CRS spellings and a no-point-moved guard, shared by the CRS suites.

Audit PR B (`docs/increments/python-audit.md`, section 9): a CRS spelt as a
PROJ string of an EPSG code is that code, so a site that meets one must move
no point. `refuse_point_moves` makes any point-moving `Transformer` method fail
the test. Building a transformer stays allowed: `crs.same_crs` builds one to
compare two CRSs and moves nothing with it.
"""

from __future__ import annotations

import warnings
from typing import Any

import pytest
from pyproj import CRS, Transformer

POINT_MOVING = ("transform", "itransform", "transform_bounds")
#: UTM 33 on GRS80 with longitudes counted from Paris: PROJ identifies it as
#: EPSG:25833 at confidence 70, yet its points lie 185 to 215 km off
#: (`test_crs.py`'s not-the-same pairs, `test_domain.py`'s refused flag).
UTM33_PARIS = "+proj=utm +zone=33 +ellps=GRS80 +units=m +pm=paris +no_defs"


def proj4_of(epsg: int) -> str:
    """pyproj's PROJ string of `epsg`, without its lossy-conversion warning."""
    with warnings.catch_warnings(action="ignore", category=UserWarning):
        return CRS.from_epsg(epsg).to_proj4()


def refuse_point_moves(monkeypatch: pytest.MonkeyPatch) -> None:
    """From here on, each method in `POINT_MOVING` raises `AssertionError`
    naming itself."""

    def refusing(name: str) -> Any:
        def refuse(*args: Any, **kwargs: Any) -> Any:
            raise AssertionError(f"Transformer.{name} was called")

        return refuse

    for name in POINT_MOVING:
        monkeypatch.setattr(Transformer, name, refusing(name))
