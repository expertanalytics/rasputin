"""CRS text and the one reprojection site (increment 15b).

`docs/increments/15-dem-mosaic.md` R2, R9 and I8. `_transformer` holds the only
`Transformer.from_crs` in `src_python/`, always with `always_xy=True`: x is
easting or longitude and y northing or latitude, whatever the CRS's own axis
order, which is the convention a GeoTIFF's model space and GeoJSON both use.

Every refusal is a `ValueError`. pyproj's `CRSError` is a `RuntimeError`, and
would escape every `except ValueError` that turns a refusal into a usage error.

Pure: pyproj and numpy. No paths, no `_core`; CRS never crosses into the core.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt
from pyproj import CRS, Transformer
from pyproj.exceptions import CRSError

Xy = npt.NDArray[np.float64]


def parse_crs(text: str | CRS) -> CRS:
    """Anything `CRS.from_user_input` accepts, or a `ValueError` naming it."""
    if isinstance(text, CRS):
        return text
    try:
        return CRS.from_user_input(text)
    except CRSError as exc:
        raise ValueError(f"cannot read the CRS {text!r}: {exc}") from exc


def reprojector(src: str | CRS, dst: str | CRS) -> Callable[[Any], Xy]:
    """A map from an `(N, 2)` array-like of `(x, y)` in `src` to an `(N, 2)`
    float64 array in `dst`: pyproj's `always_xy` transform, bit for bit.
    A point with no image comes back as `inf`; the caller decides."""
    transformer = _transformer(src, dst)

    def apply(xy: Any) -> Xy:
        points = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        x, y = transformer.transform(points[:, 0], points[:, 1])
        return np.column_stack([np.asarray(x, np.float64), np.asarray(y, np.float64)])

    return apply


def transform_description(src: str | CRS, dst: str | CRS) -> str:
    """What PROJ picked for `src` to `dst` (a datum shift included), for the record."""
    return str(_transformer(src, dst).description)


def crs_label(crs: str | CRS) -> str:
    """`EPSG:n` when pyproj finds an exact EPSG code (so OGC's CRS84 is not
    `EPSG:4326`), otherwise ASCII text pyproj parses back to the same CRS."""
    parsed = parse_crs(crs)
    code = parsed.to_epsg(min_confidence=100)
    if code is not None:
        return f"EPSG:{code}"
    return parsed.to_string().encode("ascii", "backslashreplace").decode("ascii")


def _transformer(src: str | CRS, dst: str | CRS) -> Transformer:
    return Transformer.from_crs(parse_crs(src), parse_crs(dst), always_xy=True)
