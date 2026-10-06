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

from geographic_fixtures import GEO_ASCII_PARAMS, GEO_DOUBLE_PARAMS
from geotiff_fixtures import (
    AUSTRIA_KEYS,
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
    GEOG_TOWGS84,
    GEOGRAPHIC_TYPE,
    GRADS_WITH_PARIS_MERIDIAN,
    NOT_ONE_NUMBER_DEFECTS,
    PACKBITS,
    PARAMETRIC_DEFECTS,
    PROJECTED_CS_TYPE,
    REFUSALS,
    ROWS,
    SCALE,
    THUNDERSCAN,
    TIE_X,
    TIE_Y,
    TIEPOINT,
    TOWGS84_SEVEN,
    TOWGS84_THREE,
    UNDECODABLE,
    UNDECODABLE_WITH_CODECS,
    UNKNOWN_COMPRESSION,
    ParametricDefect,
    Refusal,
    austria_keys,
    austria_tiff,
    elevations,
    floating_point_predictor_tiff,
    lambert_on,
    micro_tiff,
    packbits_tiff,
    transverse_mercator_on,
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
    """What `decode_dem` asks of tifffile: open, walk every page, decode page 0."""
    with tifffile.TiffFile(stream) as tif:
        tuple(tif.pages)
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


def test_truncated_second_ifd_opens_and_fails_in_the_page_walk() -> None:
    """The stream opens, page 0 decodes, and only reaching page 1 fails. So a
    wrap around `TiffFile(source)` alone does not catch it: the page walk must
    be inside the "TIFF structure" stage too (reviewer's finding on round 3)."""
    build, stage = UNDECODABLE["truncated_second_ifd"]
    assert stage == "TIFF structure"
    with tifffile.TiffFile(build()) as tif:
        assert tif.pages.first.shape == (ROWS, COLS)
        tif.pages.first.asarray()
        with pytest.raises(tifffile.TiffFileError, match="corrupted IFD structure"):
            tif.pages[1]
    with (
        tifffile.TiffFile(build()) as tif,
        pytest.raises(tifffile.TiffFileError, match="corrupted IFD structure"),
    ):
        tuple(tif.pages)


@pytest.mark.skipif(
    importlib.util.find_spec("imagecodecs") is None, reason="writing LZW needs imagecodecs"
)
def test_corrupt_lzw_fixture_defeats_tifffile_alone() -> None:
    build, _ = UNDECODABLE_WITH_CODECS["corrupt_lzw"]
    assert tifffile.TiffFile(build()).pages.first.compression == 5
    with pytest.raises(Exception):  # noqa: B017 - imagecodecs' error class (§14)
        _tifffile_alone(build())


# ---------------------------------------------------------------------------
# A CRS given by parameters (`docs/increments/geotiff-crs-by-parameters.md`,
# section 8): the Austrian key set, its variants, and one defect per check.
# ---------------------------------------------------------------------------

#: Section 1's table as tifffile names and decodes it (enums compare as ints).
AUSTRIA_AS_READ = {
    "GTModelTypeGeoKey": 1,
    "GTRasterTypeGeoKey": 1,
    "GTCitationGeoKey": "MGI_Austria_Lambert",
    "GeographicTypeGeoKey": 4312,
    "GeogCitationGeoKey": "MGI",
    "GeogAngularUnitsGeoKey": 9102,
    "GeogSemiMajorAxisGeoKey": 6377397.155,
    "GeogInvFlatteningGeoKey": 299.1528128000033,
    "ProjectedCSTypeGeoKey": 32767,
    "ProjectionGeoKey": 32767,
    "ProjCoordTransGeoKey": 8,
    "ProjLinearUnitsGeoKey": 9001,
    "ProjStdParallel1GeoKey": 46.0,
    "ProjStdParallel2GeoKey": 49.0,
    "ProjFalseOriginLongGeoKey": 13.33333333300013,
    "ProjFalseOriginLatGeoKey": 47.5,
    "ProjFalseOriginEastingGeoKey": 400000.0,
    "ProjFalseOriginNorthingGeoKey": 400000.0,
}


def _geokeys(stream: Any) -> dict[str, Any]:
    geo = tifffile.TiffFile(stream).geotiff_metadata
    assert geo is not None
    return {k: v for k, v in geo.items() if k.endswith("GeoKey")}


def test_austrian_fixture_reads_back_as_the_files_eighteen_keys() -> None:
    """Exact equality, values included: the fixture is the real file's key set
    (section 1, read from the file with tifffile 2026.9.20), on the baseline grid."""
    assert len(AUSTRIA_KEYS) == 18
    assert _geokeys(austria_tiff()) == AUSTRIA_AS_READ
    np.testing.assert_array_equal(tifffile.TiffFile(austria_tiff()).asarray(), elevations())


def test_parameter_keys_are_doubles_and_citations_ascii() -> None:
    """The on-disk form a GDAL-written file has: 34736 and 34737 are present
    and the directory points into them (location != 0) for those keys."""
    page = tifffile.TiffFile(austria_tiff()).pages.first
    assert GEO_DOUBLE_PARAMS in page.tags and GEO_ASCII_PARAMS in page.tags
    directory = page.tags[34735].value
    locations = {directory[i]: directory[i + 1] for i in range(4, len(directory), 4)}
    assert locations[3084] == GEO_DOUBLE_PARAMS
    assert locations[1026] == GEO_ASCII_PARAMS
    assert locations[3072] == 0


def test_lambert_on_another_datum_drops_only_the_ellipsoid_keys() -> None:
    expected = {
        k: v
        for k, v in AUSTRIA_AS_READ.items()
        if k not in ("GeogSemiMajorAxisGeoKey", "GeogInvFlatteningGeoKey")
    }
    assert _geokeys(austria_tiff(lambert_on(4326))) == {**expected, "GeographicTypeGeoKey": 4326}


def test_transverse_mercator_variant_carries_the_natural_origin_keys() -> None:
    geo = _geokeys(austria_tiff(transverse_mercator_on(4258, 15.0)))
    assert geo["ProjCoordTransGeoKey"] == 1
    assert geo["GeographicTypeGeoKey"] == 4258
    assert (
        geo["ProjNatOriginLatGeoKey"], geo["ProjNatOriginLongGeoKey"],
        geo["ProjScaleAtNatOriginGeoKey"], geo["ProjFalseEastingGeoKey"],
        geo["ProjFalseNorthingGeoKey"],
    ) == (0.0, 15.0, 0.9996, 500000.0, 0.0)  # fmt: skip
    assert not any(k.startswith(("ProjStdParallel", "ProjFalseOrigin")) for k in geo)
    assert "GeogSemiMajorAxisGeoKey" not in geo and "GeogInvFlatteningGeoKey" not in geo


#: Red test 8's defects and red test 13's, each one change to `AUSTRIA_KEYS`.
_ONE_CHANGE = (*PARAMETRIC_DEFECTS, *NOT_ONE_NUMBER_DEFECTS)
BY_DEFECT = pytest.mark.parametrize("defect", _ONE_CHANGE, ids=[d.name for d in _ONE_CHANGE])


@BY_DEFECT
def test_parametric_fixture_carries_its_defect(defect: ParametricDefect) -> None:
    assert defect.witness(tifffile.TiffFile(defect.build()).geotiff_metadata or {})


@BY_DEFECT
def test_austrian_fixture_does_not_carry_the_defect(defect: ParametricDefect) -> None:
    assert not defect.witness(tifffile.TiffFile(austria_tiff()).geotiff_metadata or {})


@BY_DEFECT
def test_parametric_fixture_differs_from_the_austrian_set_in_one_key(
    defect: ParametricDefect,
) -> None:
    """Exactly one change: the one key the defect names, added, removed or changed."""
    before, after = _geokeys(austria_tiff()), _geokeys(defect.build())
    changed = {k for k in before.keys() | after.keys() if before.get(k) != after.get(k)}
    assert len(changed) == 1, changed


@pytest.mark.parametrize("shift", [TOWGS84_SEVEN, TOWGS84_THREE], ids=["seven", "three"])
def test_towgs84_reads_back_as_all_its_doubles(shift: tuple[float, ...]) -> None:
    """2062 with seven doubles and with three: one key, its whole tuple."""
    assert _geokeys(austria_tiff(austria_keys({GEOG_TOWGS84: shift})))["GeogTOWGS84GeoKey"] == shift


@pytest.mark.parametrize(
    "defect", NOT_ONE_NUMBER_DEFECTS, ids=[d.name for d in NOT_ONE_NUMBER_DEFECTS]
)
def test_not_one_number_survives_the_write(defect: ParametricDefect) -> None:
    """Red test 13's witness, stated without the witness function: the key
    reads back as a float NaN, a float +inf, or a tuple of two finite floats."""
    (written,) = defect.changes.values()
    _, name, _ = defect.must_name
    value = _geokeys(defect.build())[name]
    if isinstance(written, tuple):
        assert isinstance(value, tuple) and len(value) == 2
        assert all(isinstance(v, float) and math.isfinite(v) for v in value)
    elif isinstance(written, float) and math.isnan(written):
        assert isinstance(value, float) and math.isnan(value)
    else:
        assert value == math.inf


def test_grads_fixture_carries_the_unit_and_the_paris_meridian() -> None:
    """Red test 14's fixture: 2054 = 9105 (grad) and 2061 = 2.5969213 (Paris in
    grads), and otherwise the Austrian set; the Austrian set has no 2061."""
    before, after = _geokeys(austria_tiff()), _geokeys(austria_tiff(GRADS_WITH_PARIS_MERIDIAN))
    assert "GeogPrimeMeridianLongGeoKey" not in before
    assert after["GeogAngularUnitsGeoKey"] == 9105
    assert after["GeogPrimeMeridianLongGeoKey"] == 2.5969213
    changed = {k for k in before.keys() | after.keys() if before.get(k) != after.get(k)}
    assert changed == {"GeogAngularUnitsGeoKey", "GeogPrimeMeridianLongGeoKey"}
