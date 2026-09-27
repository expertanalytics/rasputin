"""`tin_engine.io.geotiff.read_meta`: the header phase alone (increment 15a, R3).

G1 and G2 of the parked design, carried by `docs/increments/15-dem-mosaic.md`:
`read_meta` is `decode_dem` stopped before the pixels. It must agree with
`decode_dem(...).meta` on every valid file, refuse every header defect with
the same message, and **not** read the pixels, so a file whose strip is
damaged still gives its header.

HOW THIS FILE GOES RED: `read_meta` is looked up in a fixture, so each test
fails on its own with `AttributeError` while it is missing.
"""

from __future__ import annotations

import io
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest

from geotiff_fixtures import (
    GT_RASTER_TYPE,
    HAS_CODECS,
    KARTVERKET,
    METRE,
    PIXEL_IS_AREA,
    REFUSALS,
    VERTICAL_UNITS,
    Refusal,
    corrupt_deflate,
    elevations,
    micro_tiff,
    needs_codecs,
    truncated_strip,
    with_keys,
)
from tin_engine.io import geotiff
from tin_engine.io.models import GeoTiffError

VALID: dict[str, Callable[[], io.BytesIO]] = {
    "baseline": micro_tiff,
    "area-registered": lambda: micro_tiff(geokeys=with_keys({GT_RASTER_TYPE: PIXEL_IS_AREA})),
    "nodata-tag": lambda: micro_tiff(nodata="-32767"),
    "vertical-metres": lambda: micro_tiff(geokeys=with_keys({VERTICAL_UNITS: METRE})),
    "int16-promoted": lambda: micro_tiff(elevations(np.int16)),
    "float64": lambda: micro_tiff(elevations(np.float64)),
    "deflate": lambda: micro_tiff(compression="deflate"),
    "larger": lambda: micro_tiff(elevations(rows=7, cols=9)),
}


@pytest.fixture
def read_meta() -> Callable[..., Any]:
    return geotiff.read_meta  # type: ignore[attr-defined, no-any-return]


@pytest.mark.parametrize("build", VALID.values(), ids=VALID.keys())
def test_g1_read_meta_equals_the_decoded_meta(
    read_meta: Callable[..., Any], build: Callable[[], io.BytesIO]
) -> None:
    assert read_meta(build()) == geotiff.decode_dem(build()).meta


def test_g1_with_a_caller_nodata(read_meta: Callable[..., Any]) -> None:
    assert (
        read_meta(micro_tiff(), nodata=-9999.0)
        == geotiff.decode_dem(micro_tiff(), nodata=-9999.0).meta
    )


@needs_codecs
def test_g1_on_the_real_tile(read_meta: Callable[..., Any]) -> None:
    with KARTVERKET.open("rb") as stream:
        meta = read_meta(stream)
    with KARTVERKET.open("rb") as stream:
        assert meta == geotiff.decode_dem(stream).meta


def _header_refusals() -> list[Any]:
    # With the codecs extra present, the LZW-relabelled file is not refused by
    # its header: it fails later, decoding pixels, which read_meta never does.
    return [
        pytest.param(
            r, id=r.name, marks=pytest.mark.skipif(HAS_CODECS, reason="LZW decodes with the extra")
        )
        if r.name == "refuses_missing_codec"
        else pytest.param(r, id=r.name)
        for r in REFUSALS
    ]


@pytest.mark.parametrize("refusal", _header_refusals())
def test_g2_every_header_refusal_fires_with_the_same_message(
    read_meta: Callable[..., Any], refusal: Refusal
) -> None:
    with pytest.raises(GeoTiffError) as from_decode:
        geotiff.decode_dem(refusal.build(), **refusal.decode_kwargs)
    with pytest.raises(GeoTiffError) as from_header:
        read_meta(refusal.build(), **refusal.decode_kwargs)
    assert str(from_header.value) == str(from_decode.value)


def test_a_bool_nodata_is_a_type_error(read_meta: Callable[..., Any]) -> None:
    with pytest.raises(TypeError):
        read_meta(micro_tiff(), nodata=True)


@pytest.mark.parametrize(
    "build", [truncated_strip, corrupt_deflate], ids=["truncated-strip", "corrupt-deflate"]
)
def test_the_pixels_are_not_read(
    read_meta: Callable[..., Any], build: Callable[[], io.BytesIO]
) -> None:
    """The fixture's header is the baseline's; only its strip is damaged."""
    with pytest.raises(GeoTiffError, match="pixel data"):
        geotiff.decode_dem(build())
    meta = read_meta(build())
    assert (meta.rows, meta.cols) == (3, 4)
    assert meta == geotiff.decode_dem(micro_tiff(compression="deflate")).meta
