"""The micro-GeoTIFF fixtures carry the defect they are named for, and no other.

Green before `tin_engine.io.geotiff` exists, and it has to be: this file is
the evidence that `test_io_geotiff.py`'s red is red for the right reason. It
reads every fixture with tifffile alone. No production code is imported.

Two properties per refusal fixture:

1. its `witness` holds, so the named defect is really in the bytes;
2. the valid baseline does *not* satisfy that witness, so the witness can
   fail (`PRINCIPLES.md` A3) and the defect is not shared with the baseline.
"""

from __future__ import annotations

import numpy as np
import pytest
import tifffile

from geotiff_fixtures import (
    BASE_KEYS,
    COLS,
    DELTA_X,
    DELTA_Y,
    FLOATING_POINT_PREDICTOR,
    PACKBITS,
    REFUSALS,
    ROWS,
    TIE_X,
    TIE_Y,
    Refusal,
    elevations,
    floating_point_predictor_tiff,
    micro_tiff,
    packbits_tiff,
)

BY_NAME = pytest.mark.parametrize("refusal", REFUSALS, ids=[r.name for r in REFUSALS])


def test_catalogue_names_all_eighteen_refusals_once() -> None:
    """Design §5 and §13 count eighteen; the catalogue must too, with no repeat."""
    names = [r.name for r in REFUSALS]
    assert len(names) == 18
    assert len(set(names)) == 18


def test_baseline_is_a_valid_point_registered_projected_tile() -> None:
    tif = tifffile.TiffFile(micro_tiff())
    page = tif.pages.first
    geo = tif.geotiff_metadata

    assert len(tif.pages) == 1
    assert page.shape == (ROWS, COLS)
    assert page.dtype == np.float32
    assert page.samplesperpixel == 1
    assert page.compression == 1
    assert geo is not None
    assert geo["ModelTiepoint"] == [0.0, 0.0, 0.0, TIE_X, TIE_Y, 0.0]
    assert geo["ModelPixelScale"] == [DELTA_X, DELTA_Y, 0.0]
    assert {k: int(v) for k, v in geo.items() if k.endswith("GeoKey")} == {
        "GTModelTypeGeoKey": BASE_KEYS[1024],
        "GTRasterTypeGeoKey": BASE_KEYS[1025],
        "ProjectedCSTypeGeoKey": BASE_KEYS[3072],
    }
    assert 34264 not in page.tags
    assert 42113 not in page.tags
    np.testing.assert_array_equal(page.asarray(), elevations())


@BY_NAME
def test_fixture_carries_its_defect(refusal: Refusal) -> None:
    assert refusal.witness(tifffile.TiffFile(refusal.build()))


@BY_NAME
def test_baseline_does_not_carry_the_defect(refusal: Refusal) -> None:
    if refusal.decode_kwargs:
        # The contradiction is between the tag and the caller; the baseline has
        # no tag, so the witness is checked against a tag-bearing twin instead.
        assert not refusal.witness(tifffile.TiffFile(micro_tiff(nodata="-9999")))
    else:
        assert not refusal.witness(tifffile.TiffFile(micro_tiff()))


def test_floating_point_predictor_fixture_is_a_float_tile_declaring_predictor_3() -> None:
    """§8: GDAL writes `PREDICTOR=3` for float DEMs. The baseline declares no predictor."""
    page = tifffile.TiffFile(floating_point_predictor_tiff()).pages.first
    assert page.predictor == FLOATING_POINT_PREDICTOR
    assert page.dtype == np.float32
    assert page.compression in (8, 32946)  # Deflate, either code
    assert page.shape == (ROWS, COLS)
    assert tifffile.TiffFile(micro_tiff()).pages.first.predictor == 1


def test_packbits_fixture_is_real_packbits_that_tifffile_decodes() -> None:
    """The strip really is PackBits, and tifffile alone decodes it to the baseline.

    That is what makes the PackBits test in `test_io_geotiff.py` a statement
    about the reader and not about the fixture.
    """
    page = tifffile.TiffFile(packbits_tiff()).pages.first
    assert page.compression == PACKBITS
    assert page.tags[279].value == (ROWS * COLS * 4 + 1,)
    np.testing.assert_array_equal(page.asarray(), elevations())
