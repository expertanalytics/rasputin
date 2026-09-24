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

import functools
import importlib.util
import math
from typing import Any

import numpy as np
import pytest
import tifffile

from geotiff_fixtures import (
    BASE_KEYS,
    COLS,
    DELTA_X,
    DELTA_Y,
    EPSG_COMPOUND,
    EPSG_COMPOUND_GEOGRAPHIC,
    EPSG_GEOCENTRIC,
    EPSG_PROJECTED_3D,
    EPSG_UTM33,
    FLOATING_POINT_PREDICTOR,
    GEOGRAPHIC_TYPE,
    PACKBITS,
    PROJECTED_CS_TYPE,
    REFUSALS,
    ROWS,
    SCALE,
    THUNDERSCAN,
    TIE_X,
    TIE_Y,
    TIEPOINT,
    UNDECODABLE,
    UNDECODABLE_WITH_CODECS,
    UNKNOWN_COMPRESSION,
    Refusal,
    elevations,
    floating_point_predictor_tiff,
    micro_tiff,
    packbits_tiff,
    with_compression_tag,
    with_keys,
)

BY_NAME = pytest.mark.parametrize("refusal", REFUSALS, ids=[r.name for r in REFUSALS])


def test_catalogue_names_all_nineteen_refusals_once() -> None:
    """Design §5 names eighteen, and round 3 adds refusal 13a
    (`refuses_compound_crs`). The catalogue has each once."""
    names = [r.name for r in REFUSALS]
    assert len(names) == 19
    assert len(set(names)) == 19
    assert "refuses_compound_crs" in names


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


@pytest.mark.parametrize("code", [EPSG_UTM33, EPSG_GEOCENTRIC], ids=["projected", "geocentric"])
def test_a_crs_code_in_the_geographic_key_alone_survives_the_write(code: int) -> None:
    """§5 refusal 13, round 2: the `absent_geographic_*` cases of
    `test_refuses_missing_crs` are about 2048 carrying a non-geographic code
    with no 3072. The bytes must say exactly that, or the case tests nothing."""
    tif = tifffile.TiffFile(
        micro_tiff(geokeys=with_keys({PROJECTED_CS_TYPE: None, GEOGRAPHIC_TYPE: code}))
    )
    geo = tif.geotiff_metadata
    assert geo is not None
    assert "ProjectedCSTypeGeoKey" not in geo
    assert int(geo["GeographicTypeGeoKey"]) == code


def test_non_finite_tie_point_k_and_z_survive_the_write() -> None:
    """§4, round 2, reading (a): the tie point's NaN K and infinite Z are in the
    file, so a reader that ignores them is ignoring something real."""
    tie = (0.0, 0.0, math.nan, TIE_X, TIE_Y, math.inf)
    written = tifffile.TiffFile(micro_tiff(tiepoint=tie)).pages.first.tags[33922].value
    assert math.isnan(written[2])
    assert written[5] == math.inf
    assert tuple(written[i] for i in (0, 1, 3, 4)) == (0.0, 0.0, TIE_X, TIE_Y)


# ---------------------------------------------------------------------------
# Round 3 fixtures (§5 refusals 1-3, 8, 11 and 13a; §14 choice C)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("code", [EPSG_COMPOUND, EPSG_PROJECTED_3D, EPSG_COMPOUND_GEOGRAPHIC])
def test_a_three_axis_crs_code_survives_the_write(code: int) -> None:
    """Refusal 13a is about 3072 carrying a 3-axis CRS. The code reads back."""
    geo = tifffile.TiffFile(
        micro_tiff(geokeys=with_keys({PROJECTED_CS_TYPE: code}))
    ).geotiff_metadata
    assert geo is not None
    assert int(geo["ProjectedCSTypeGeoKey"]) == code


@pytest.mark.parametrize(
    ("tag", "kwargs", "count"),
    [
        (33922, {"tiepoint": TIEPOINT[:3]}, 3),
        (33922, {"tiepoint": TIEPOINT[:5]}, 5),
        (33922, {"tiepoint": (*TIEPOINT, 1.0)}, 7),
        (33922, {"tiepoint": (*TIEPOINT, 1.0, 2.0, 3.0)}, 9),
        (33550, {"scale": SCALE[:2]}, 2),
        (34264, {"tiepoint": None, "scale": None, "transformation": (1.0, 2.0, 3.0)}, 3),
        (
            34264,
            {"tiepoint": None, "scale": None, "transformation": tuple(map(float, range(12)))},
            12,
        ),
    ],
    ids=[
        "tiepoint_3",
        "tiepoint_5",
        "tiepoint_7",
        "tiepoint_9",
        "scale_2",
        "transform_3",
        "transform_12",
    ],
)
def test_a_wrong_length_tag_survives_the_write(
    tag: int, kwargs: dict[str, Any], count: int
) -> None:
    """§5, after refusal 5 (round 3): the tag has exactly the count the case
    names. The transformation cases carry no tie point and no scale."""
    tags = tifffile.TiffFile(micro_tiff(**kwargs)).pages.first.tags
    assert len(tags[tag].value) == count
    if tag == 34264:
        assert 33922 not in tags and 33550 not in tags


@pytest.mark.parametrize(
    ("tag", "kwargs"), [(33922, {"tiepoint": (0.0,)}), (33550, {"scale": (DELTA_X,)})]
)
def test_a_one_value_double_tag_reads_back_as_a_bare_float(
    tag: int, kwargs: dict[str, Any]
) -> None:
    """§5 (round 3), measured: "a DOUBLE tag of count 1 reads back as a bare
    `float`, not a tuple". That is the scalar path the `*_1_value` cases take."""
    value = tifffile.TiffFile(micro_tiff(**kwargs)).pages.first.tags[tag].value
    assert type(value) is float


def test_useframes_turns_a_full_resolution_second_page_into_a_frame() -> None:
    """§5 refusal 8 (round 3), `unparsed_frame`: with tifffile's private
    `_useframes=True`, page 1 comes back as a `TiffFrame`. Parsed normally, the
    same page is full resolution, so the frame hides a second image."""
    stream = micro_tiff(extra_pages=[(elevations(), 0)])
    parsed = tifffile.TiffFile(stream).pages[1]
    assert isinstance(parsed, tifffile.TiffPage) and not parsed.is_reduced
    stream.seek(0)
    framed = functools.partial(tifffile.TiffFile, _useframes=True)(stream)
    assert len(framed.pages) == 2
    assert isinstance(framed.pages[1], tifffile.TiffFrame)
    assert not isinstance(framed.pages[1], tifffile.TiffPage)


def test_unknown_and_undecodable_compression_codes() -> None:
    """§5 refusal 11 (round 3): 60000 is not in tifffile's enum; 32809 is
    (`THUNDERSCAN`), and tifffile has no decoder for it either way."""
    assert UNKNOWN_COMPRESSION not in tifffile.COMPRESSION
    assert tifffile.COMPRESSION(THUNDERSCAN).name == "THUNDERSCAN"
    assert THUNDERSCAN not in tifffile.TIFF.DECOMPRESSORS
    for code in (UNKNOWN_COMPRESSION, THUNDERSCAN):
        page = tifffile.TiffFile(with_compression_tag(micro_tiff(), code)).pages.first
        assert int(page.compression) == code


def _tifffile_alone(stream: Any) -> None:
    with tifffile.TiffFile(stream) as tif:
        tif.pages.first.asarray()


@pytest.mark.parametrize("name", list(UNDECODABLE))
def test_undecodable_fixture_defeats_tifffile_alone(name: str) -> None:
    """§14, choice C: each stream makes tifffile (or its codec) raise with no
    production code involved, and the baseline does not. So the wrapping test
    in `test_io_geotiff.py` is red because the wrap is missing."""
    build, _ = UNDECODABLE[name]
    with pytest.raises(Exception):  # noqa: B017 - the type is tifffile's and varies (§14)
        _tifffile_alone(build())
    _tifffile_alone(micro_tiff())


@pytest.mark.skipif(
    importlib.util.find_spec("imagecodecs") is None, reason="writing LZW needs imagecodecs"
)
def test_corrupt_lzw_fixture_defeats_tifffile_alone() -> None:
    build, _ = UNDECODABLE_WITH_CODECS["corrupt_lzw"]
    assert tifffile.TiffFile(build()).pages.first.compression == 5
    with pytest.raises(Exception):  # noqa: B017 - imagecodecs' error class (§14)
        _tifffile_alone(build())
