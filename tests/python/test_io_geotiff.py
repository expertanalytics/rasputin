"""`tin_engine.io.geotiff.decode_dem`: increment 11, a GeoTIFF stream to a `DemTile`.

The design is `docs/increments/11-raster-ingestion.md`; section numbers below
are its. The reader's job is refusal (ruling 3), so most of this file is the
eighteen named refusals of §5, one crafted micro-TIFF each, built by
`geotiff_fixtures.py`. `test_geotiff_fixtures.py` shows, with tifffile alone,
that every fixture carries exactly the defect its refusal names.

HOW THIS FILE GOES RED. Production modules are imported inside the `decode`
and `geotiff_error` fixtures rather than at module scope. A module-scope import
of a missing module is one collection error that aborts the whole session;
here each test fails on its own, at fixture setup, with `ModuleNotFoundError`,
and the rest of `tests/python` still runs and reports.

WHAT A REFUSAL MUST SAY. §5, "The exception type": one `GeoTiffError`, a
`ValueError`, whose message names the tag or GeoKey by number and by name, and
the file's value. Matching is case-insensitive substring, so the suite pins
*what* is named and not the sentence around it. A bare number such as a page
index is matched as a whole number (`_names_numbers`), so `3` is not found
inside `317`. Every refusal names a tag (amended, problem 4): degenerate
shape, ambiguous pages, unsupported dtype and missing codec name the baseline
TIFF tag their decision reads, as in §5's table.
"""

from __future__ import annotations

import ast
import importlib.util
import math
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import tifffile

from geotiff_fixtures import (
    COLS,
    DELTA_X,
    DELTA_Y,
    EPSG_FEET,
    EPSG_GEOCENTRIC,
    EPSG_UNRESOLVABLE,
    EPSG_UTM33,
    EPSG_WGS84,
    FLOATING_POINT_PREDICTOR,
    FOOT,
    GEOG_LINEAR_UNITS,
    GEOGRAPHIC_TYPE,
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    LZW,
    METRE,
    PIXEL_IS_AREA,
    PROJ_LINEAR_UNITS,
    PROJECTED_CS_TYPE,
    REFUSALS,
    ROWS,
    TIE_X,
    TIE_Y,
    USER_DEFINED,
    VERTICAL_UNITS,
    elevations,
    floating_point_predictor_tiff,
    micro_tiff,
    packbits_tiff,
    with_compression_tag,
    with_keys,
)

HERE = Path(__file__).resolve().parent
KARTVERKET = HERE.parent / "fixtures" / "dem_archive" / "7908_3_10m_z33.tif"
HAS_CODECS = importlib.util.find_spec("imagecodecs") is not None
needs_codecs = pytest.mark.skipif(not HAS_CODECS, reason="LZW needs the `codecs` extra (§8)")
without_codecs = pytest.mark.skipif(HAS_CODECS, reason="tests the refusal when the extra is absent")

REFUSAL = {r.name: r for r in REFUSALS}

Decode = Callable[..., Any]


@pytest.fixture
def decode() -> Decode:
    from tin_engine.io.geotiff import decode_dem

    return decode_dem


@pytest.fixture
def geotiff_error() -> type[Exception]:
    from tin_engine.io.models import GeoTiffError

    return GeoTiffError


@pytest.fixture
def refused(decode: Decode, geotiff_error: type[Exception]) -> Callable[..., str]:
    """Decode, require `GeoTiffError`, require each `names` token in the message."""

    def check(stream: Any, *names: str, **decode_kwargs: Any) -> str:
        with pytest.raises(geotiff_error) as info:
            decode(stream, **decode_kwargs)
        message = str(info.value)
        for token in names:
            assert token.lower() in message.lower(), f"{token!r} not named in {message!r}"
        return message

    return check


def refuse_named(refused: Callable[..., str], name: str) -> str:
    entry = REFUSAL[name]
    return refused(entry.build(), *entry.must_name, **entry.decode_kwargs)


def _names_numbers(message: str, *numbers: int) -> None:
    """Each number appears in `message` as a whole number, not as part of another."""
    for number in numbers:
        pattern = rf"(?<![\d.]){number}(?!\d)"
        assert re.search(pattern, message), f"{number} not named in {message!r}"


# ---------------------------------------------------------------------------
# The exception type (§5)
# ---------------------------------------------------------------------------


def test_geotiff_error_is_a_value_error(geotiff_error: type[Exception]) -> None:
    assert issubclass(geotiff_error, ValueError)


# ---------------------------------------------------------------------------
# The happy path, on micro-TIFFs
# ---------------------------------------------------------------------------


class TestPlacement:
    """§4: the node grid, from the tie point, scale and registration."""

    def test_point_registered_nodes_sit_on_the_tie_point(self, decode: Decode) -> None:
        meta = decode(micro_tiff()).meta
        assert (meta.x_min, meta.y_max) == (TIE_X, TIE_Y)
        assert (meta.delta_x, meta.delta_y) == (DELTA_X, DELTA_Y)
        assert meta.pixel_is_area is False

    def test_area_registered_nodes_are_shifted_inward_half_a_cell(self, decode: Decode) -> None:
        """delta_x != delta_y, so a swapped half-cell cannot pass."""
        stream = micro_tiff(geokeys=with_keys({GT_RASTER_TYPE: PIXEL_IS_AREA}))
        meta = decode(stream).meta
        assert meta.x_min == TIE_X + DELTA_X / 2
        assert meta.y_max == TIE_Y - DELTA_Y / 2
        assert (meta.delta_x, meta.delta_y) == (DELTA_X, DELTA_Y)
        assert meta.pixel_is_area is True

    @pytest.mark.parametrize(
        ("i", "j", "area", "x_min", "y_max"),
        [
            (2.0, 1.0, False, TIE_X - 2 * DELTA_X, TIE_Y + 1 * DELTA_Y),
            (2.0, 1.0, True, TIE_X - 2 * DELTA_X + DELTA_X / 2, TIE_Y + 1 * DELTA_Y - DELTA_Y / 2),
            (0.5, 0.5, True, TIE_X, TIE_Y),
        ],
        ids=["point_registered", "area_registered", "area_registered_at_a_cell_centre"],
    )
    def test_tie_point_off_pixel_zero_is_converted(
        self, decode: Decode, i: float, j: float, area: bool, x_min: float, y_max: float
    ) -> None:
        """§4, amended (problem 7): raster point (I, J) maps to model point (X, Y).

        I is the column and J the row, and I != J with delta_x != delta_y, so
        swapping either the axes or the spacings cannot pass. Every expected
        value is exact in binary floating point, so equality is exact.
        """
        keys = with_keys({GT_RASTER_TYPE: PIXEL_IS_AREA}) if area else None
        meta = decode(micro_tiff(tiepoint=(i, j, 0.0, TIE_X, TIE_Y, 0.0), geokeys=keys)).meta
        assert (meta.x_min, meta.y_max) == (x_min, y_max)
        assert (meta.delta_x, meta.delta_y) == (DELTA_X, DELTA_Y)
        assert meta.pixel_is_area is area

    @pytest.mark.parametrize(
        ("k", "z"), [(7.0, 100.0), (math.nan, math.inf)], ids=["finite", "non_finite"]
    )
    def test_tie_point_k_and_z_are_ignored(self, decode: Decode, k: float, z: float) -> None:
        """§4: with ScaleZ zero there is no vertical mapping for K and Z to offset.

        Amended (round 2), reading (a): "ignored" covers non-finite K and Z too.
        Refusal 1's finiteness rule is for I, J, X and Y only.
        `test_geotiff_fixtures.py` shows the NaN and inf survive the write.
        """
        meta = decode(micro_tiff(tiepoint=(0.0, 0.0, k, TIE_X, TIE_Y, z))).meta
        assert (meta.x_min, meta.y_max) == (TIE_X, TIE_Y)

    def test_dimensions_come_from_the_array_shape_as_python_ints(self, decode: Decode) -> None:
        """Prior art §5.7: `(12.0, 16.0) == (12, 16)` passed. Types are checked too."""
        tile = decode(micro_tiff())
        assert (tile.meta.rows, tile.meta.cols) == tile.array.shape == (ROWS, COLS)
        assert type(tile.meta.rows) is int
        assert type(tile.meta.cols) is int

    def test_epsg_is_recorded(self, decode: Decode) -> None:
        assert decode(micro_tiff()).meta.epsg == EPSG_UTM33


class TestArray:
    """§7: numpy 2-D, C-contiguous, read-only, values intact."""

    def test_values_survive_decode(self, decode: Decode) -> None:
        np.testing.assert_array_equal(decode(micro_tiff()).array, elevations())

    def test_array_is_two_dimensional_and_c_contiguous(self, decode: Decode) -> None:
        array = decode(micro_tiff()).array
        assert isinstance(array, np.ndarray)
        assert array.ndim == 2
        assert array.flags.c_contiguous

    def test_array_is_read_only(self, decode: Decode) -> None:
        array = decode(micro_tiff()).array
        assert not array.flags.writeable
        with pytest.raises(ValueError, match="read-only"):
            array[0, 0] = 1.0

    def test_deflate_tile_decodes_without_optional_dependencies(self, decode: Decode) -> None:
        tile = decode(micro_tiff(compression="deflate"))
        np.testing.assert_array_equal(tile.array, elevations())

    def test_packbits_tile_decodes_without_optional_dependencies(self, decode: Decode) -> None:
        """§8, amended (problem 1): PackBits reads with tifffile's built-in decoder.

        A reader that kept a list of scheme names, instead of asking tifffile,
        would refuse this file. The probe must not.
        """
        np.testing.assert_array_equal(decode(packbits_tiff()).array, elevations())


class TestFrozenModels:
    def test_meta_is_frozen(self, decode: Decode) -> None:
        from pydantic import ValidationError

        meta = decode(micro_tiff()).meta
        with pytest.raises(ValidationError, match="frozen"):
            meta.x_min = 0.0

    def test_tile_is_frozen(self, decode: Decode) -> None:
        from pydantic import ValidationError

        tile = decode(micro_tiff())
        with pytest.raises(ValidationError, match="frozen"):
            tile.array = elevations()


class TestShapeAgreement:
    """§7, amended (problem 6), B1: `RasterMeta` keeps `rows` and `cols` as
    `StrictInt`, and `DemTile` refuses a tile where they disagree with
    `array.shape`. The error is Pydantic's `ValidationError`, not
    `GeoTiffError`: it is a programming error, not a fact about a file.
    """

    @pytest.fixture
    def models(self) -> Any:
        from tin_engine.io import models

        return models

    def _meta(self, decode: Decode, models: Any, **changes: Any) -> Any:
        fields = {**decode(micro_tiff()).meta.model_dump(), **changes}
        return models.RasterMeta.model_validate(fields)

    def test_matching_shape_is_accepted(self, decode: Decode, models: Any) -> None:
        """The control: without it the refusals below could pass by refusing everything."""
        tile = models.DemTile(meta=self._meta(decode, models), array=elevations())
        assert (tile.meta.rows, tile.meta.cols) == tile.array.shape == (ROWS, COLS)

    @pytest.mark.parametrize(
        ("rows", "cols"),
        [(ROWS + 1, COLS), (ROWS, COLS - 1), (COLS, ROWS)],
        ids=["rows_differ", "cols_differ", "transposed"],
    )
    def test_tile_refuses_meta_disagreeing_with_array_shape(
        self,
        decode: Decode,
        models: Any,
        geotiff_error: type[Exception],
        rows: int,
        cols: int,
    ) -> None:
        from pydantic import ValidationError

        meta = self._meta(decode, models, rows=rows, cols=cols)
        with pytest.raises(ValidationError) as info:
            models.DemTile(meta=meta, array=elevations())
        assert not isinstance(info.value, geotiff_error)

    @pytest.mark.parametrize("field", ["rows", "cols"])
    def test_dimensions_are_strict_ints(self, decode: Decode, models: Any, field: str) -> None:
        """Prior art §5.7: the legacy's `(12.0, 16.0)` would be rejected at the type."""
        from pydantic import ValidationError

        with pytest.raises(ValidationError, match=field):
            self._meta(decode, models, **{field: float(ROWS if field == "rows" else COLS)})

    @pytest.mark.parametrize(
        ("nodata", "source"),
        [(math.nan, "tag"), (math.inf, "caller"), (-math.inf, "caller")],
        ids=["nan", "inf", "negative_inf"],
    )
    def test_meta_refuses_non_finite_nodata(
        self,
        decode: Decode,
        models: Any,
        geotiff_error: type[Exception],
        nodata: float,
        source: str,
    ) -> None:
        """§6, amended (round 2): `nodata` is `allow_inf_nan=False`, so a
        `RasterMeta` holding a NaN or infinite sentinel cannot be built. As in
        B1 this is a programming error, so Pydantic's error, not `GeoTiffError`.
        The source is never "absent", so the absent-means-None validator cannot
        be what refuses it."""
        from pydantic import ValidationError

        with pytest.raises(ValidationError, match="nodata") as info:
            self._meta(decode, models, nodata=nodata, nodata_source=source)
        assert not isinstance(info.value, geotiff_error)
        assert any(e["type"] == "finite_number" for e in info.value.errors())

    def test_meta_accepts_finite_nodata(self, decode: Decode, models: Any) -> None:
        """The control for the refusal above: a finite sentinel still builds."""
        meta = self._meta(decode, models, nodata=-9999.0, nodata_source="caller")
        assert meta.nodata == -9999.0


PROMOTION = [
    # (file dtype, array dtype, an extreme value the promotion must keep exact)
    ("int8", "float32", -128),
    ("uint8", "float32", 255),
    ("int16", "float32", -32768),
    ("uint16", "float32", 65535),
    ("int32", "float64", -(2**31)),
    ("uint32", "float64", 2**32 - 1),
    ("float32", "float32", np.float32(1234.5678)),
    ("float64", "float64", 1e300),
]


@pytest.mark.parametrize(
    ("file_dtype", "array_dtype", "extreme"), PROMOTION, ids=[p[0] for p in PROMOTION]
)
def test_promotion_table(decode: Decode, file_dtype: str, array_dtype: str, extreme: Any) -> None:
    """§7's table, one row per case, asserting the exact dtype and exact values."""
    source = elevations(file_dtype)
    source[-1, -1] = extreme
    array = decode(micro_tiff(source)).array
    assert array.dtype == np.dtype(array_dtype)
    np.testing.assert_array_equal(array, source.astype(array_dtype))
    assert array[-1, -1] == source[-1, -1]


class TestNotRefused:
    """§5, "What is explicitly *not* refused"."""

    def test_unknown_geokeys_are_ignored(self, decode: Decode) -> None:
        """Prior art §5.10: the legacy raised on any key outside its own enum."""
        stream = micro_tiff(geokeys=with_keys({GEOG_LINEAR_UNITS: METRE, 5000: 7}))
        assert decode(stream).meta.epsg == EPSG_UTM33

    def test_geographic_type_key_beside_a_projected_crs_is_accepted(self, decode: Decode) -> None:
        """Ruling 7: the reference DEM carries 2048 = 4258 *and* 3072 = 25833."""
        stream = micro_tiff(geokeys=with_keys({GEOGRAPHIC_TYPE: 4258}))
        assert decode(stream).meta.epsg == EPSG_UTM33

    def test_gdal_metadata_is_ignored(self, decode: Decode) -> None:
        stream = micro_tiff(gdal_metadata="<GDALMetadata><Item>x</Item></GDALMetadata>")
        assert decode(stream).meta.epsg == EPSG_UTM33

    def test_explicit_metre_linear_unit_is_accepted(self, decode: Decode) -> None:
        stream = micro_tiff(geokeys=with_keys({PROJ_LINEAR_UNITS: METRE}))
        assert decode(stream).meta.epsg == EPSG_UTM33

    def test_reduced_resolution_pages_are_accepted_and_page_zero_is_read(
        self, decode: Decode
    ) -> None:
        stream = micro_tiff(extra_pages=[(elevations(rows=2, cols=2) + 100, 1)])
        np.testing.assert_array_equal(decode(stream).array, elevations())


class TestVerticalUnit:
    """§5 refusal 14: absent means metres, and the assumption is recorded."""

    def test_absent_vertical_unit_is_assumed(self, decode: Decode) -> None:
        assert decode(micro_tiff()).meta.vertical_unit_assumed is True

    def test_declared_metre_vertical_unit_is_not_assumed(self, decode: Decode) -> None:
        stream = micro_tiff(geokeys=with_keys({VERTICAL_UNITS: METRE}))
        assert decode(stream).meta.vertical_unit_assumed is False


class TestNoData:
    """§6, ruling 9."""

    def test_absent_tag_means_no_sentinel(self, decode: Decode) -> None:
        meta = decode(micro_tiff()).meta
        assert meta.nodata is None
        assert meta.nodata_source == "absent"

    def test_absent_tag_is_never_guessed_from_the_array(self, decode: Decode) -> None:
        """The array holds both classic sentinels; neither may be sniffed."""
        source = elevations()
        source[0, 0], source[0, 1] = -9999.0, -32767.0
        meta = decode(micro_tiff(source)).meta
        assert meta.nodata is None
        assert meta.nodata_source == "absent"

    @pytest.mark.parametrize("file_dtype", ["float32", "int16"])
    def test_tag_is_parsed_and_recorded(self, decode: Decode, file_dtype: str) -> None:
        meta = decode(micro_tiff(elevations(file_dtype), nodata="-32767")).meta
        assert meta.nodata == -32767.0
        assert meta.nodata_source == "tag"

    @pytest.mark.parametrize(
        ("file_dtype", "text", "value"),
        [
            ("uint8", "0", 0.0),
            ("uint8", "255", 255.0),
            ("int16", "-32768", -32768.0),
            ("uint32", "4294967295", 4294967295.0),
            ("float64", "0.1", 0.1),
        ],
        ids=["uint8_min", "uint8_max", "int16_min", "uint32_max", "float64_fraction"],
    )
    def test_sentinel_representable_in_the_file_dtype_is_accepted(
        self, decode: Decode, file_dtype: str, text: str, value: float
    ) -> None:
        """§6, amended (problem 5): the check is against the *file* dtype, both
        bounds included. `0.1` is refused on float32 below and accepted here."""
        meta = decode(micro_tiff(elevations(file_dtype), nodata=text)).meta
        assert meta.nodata == value
        assert meta.nodata_source == "tag"

    @pytest.mark.parametrize("text", ["nan", "NaN"])
    def test_nan_tag_yields_no_sentinel(self, decode: Decode, text: str) -> None:
        """§6: a NaN sentinel matches nothing under `==`; `is_nodata` catches NaN.
        The file still made a declaration, so the source is the tag (5a)."""
        meta = decode(micro_tiff(nodata=text)).meta
        assert meta.nodata is None
        assert meta.nodata_source == "tag"

    def test_caller_sentinel_without_tag(self, decode: Decode) -> None:
        meta = decode(micro_tiff(), nodata=-9999.0).meta
        assert meta.nodata == -9999.0
        assert meta.nodata_source == "caller"

    def test_caller_nan_without_tag_yields_no_sentinel(self, decode: Decode) -> None:
        """§6, amended (round 2): a NaN sentinel becomes None whoever declared it.
        The caller made the declaration and the file did not, so "caller"."""
        meta = decode(micro_tiff(), nodata=math.nan).meta
        assert meta.nodata is None
        assert meta.nodata_source == "caller"

    def test_caller_sentinel_agreeing_with_tag_is_accepted(self, decode: Decode) -> None:
        """§6 table: the file would have given the same answer alone, so "tag"."""
        meta = decode(micro_tiff(nodata="-32767"), nodata=-32767.0).meta
        assert meta.nodata == -32767.0
        assert meta.nodata_source == "tag"

    def test_caller_nan_agreeing_with_nan_tag_is_accepted(self, decode: Decode) -> None:
        """§6: "equal" treats NaN as equal to NaN, so this is not a contradiction."""
        meta = decode(micro_tiff(nodata="nan"), nodata=math.nan).meta
        assert meta.nodata is None
        assert meta.nodata_source == "tag"


# ---------------------------------------------------------------------------
# The eighteen refusals (§5). Test names are the design's names.
# ---------------------------------------------------------------------------


def _tiepoint(**changes: float) -> tuple[float, ...]:
    """The baseline tie point with some of I, J, X, Y replaced."""
    values = {"i": 0.0, "j": 0.0, "x": TIE_X, "y": TIE_Y, **changes}
    return (values["i"], values["j"], 0.0, values["x"], values["y"], 0.0)


_TIEPOINT_NAMES = ("33922", "ModelTiepointTag")


@pytest.mark.parametrize(
    ("changes", "names"),
    [
        ({"tiepoint": None}, _TIEPOINT_NAMES),
        ({"scale": None}, ("33550", "ModelPixelScaleTag")),
        ({"tiepoint": _tiepoint(i=math.nan)}, (*_TIEPOINT_NAMES, "nan")),
        ({"tiepoint": _tiepoint(j=math.inf)}, (*_TIEPOINT_NAMES, "inf")),
        ({"tiepoint": _tiepoint(x=math.nan)}, (*_TIEPOINT_NAMES, "nan")),
        ({"tiepoint": _tiepoint(y=-math.inf)}, (*_TIEPOINT_NAMES, "-inf")),
    ],
    ids=["tiepoint_absent", "scale_absent", "i_nan", "j_inf", "x_nan", "y_negative_inf"],
)
def test_refuses_missing_georeferencing(
    refused: Callable[..., str], changes: dict[str, Any], names: tuple[str, ...]
) -> None:
    """§5 refusal 1, amended (problem 7): a non-finite origin is as unusable as none."""
    refused(micro_tiff(**changes), *names)


def test_refuses_model_transformation(refused: Callable[..., str]) -> None:
    """Refused even though a tie point and scale are also present."""
    refuse_named(refused, "refuses_model_transformation")


def test_refuses_multiple_tiepoints(refused: Callable[..., str]) -> None:
    refuse_named(refused, "refuses_multiple_tiepoints")


def test_refuses_nonzero_pixel_scale_z(refused: Callable[..., str]) -> None:
    refuse_named(refused, "refuses_nonzero_pixel_scale_z")


@pytest.mark.parametrize(
    ("scale", "value"),
    [
        ((0.0, DELTA_Y, 0.0), "0"),
        ((-DELTA_X, DELTA_Y, 0.0), "-10"),
        ((DELTA_X, -DELTA_Y, 0.0), "-5"),
        ((math.nan, DELTA_Y, 0.0), "nan"),
        ((DELTA_X, math.inf, 0.0), "inf"),
    ],
    ids=["x_zero", "x_negative", "y_negative", "x_nan", "y_inf"],
)
def test_refuses_nonpositive_pixel_scale(
    refused: Callable[..., str], scale: tuple[float, ...], value: str
) -> None:
    refused(micro_tiff(scale=scale), "33550", "ModelPixelScaleTag", value)


@pytest.mark.parametrize(
    ("raster_type", "names"),
    [
        (None, ("1025", "GTRasterTypeGeoKey")),
        (USER_DEFINED, ("1025", "GTRasterTypeGeoKey", "32767")),
    ],
    ids=["absent", "user_defined"],
)
def test_refuses_unknown_raster_type(
    refused: Callable[..., str], raster_type: int | None, names: tuple[str, ...]
) -> None:
    refused(micro_tiff(geokeys=with_keys({GT_RASTER_TYPE: raster_type})), *names)


_LENGTH = ("257", "ImageLength")
_WIDTH = ("256", "ImageWidth")


@pytest.mark.parametrize(
    ("rows", "cols", "names"),
    [(1, COLS, _LENGTH), (ROWS, 1, _WIDTH), (1, 1, (*_LENGTH, *_WIDTH))],
    ids=["one_row", "one_column", "one_cell"],
)
def test_refuses_degenerate_shape(
    refused: Callable[..., str], rows: int, cols: int, names: tuple[str, ...]
) -> None:
    """§5 table: the offending dimension's tag with its value, and "need at least 2"."""
    message = refused(micro_tiff(elevations(rows=rows, cols=cols)), *names, "at least 2")
    _names_numbers(message, 1)


@pytest.mark.parametrize(
    ("extra_pages", "subfiletype", "index", "count"),
    [
        ([(elevations(), 0)], 0, 1, 2),
        ([(elevations(rows=2, cols=2), 1), (elevations(), 0)], 0, 2, 3),
        ([(elevations(), 2)], 2, 1, 2),
    ],
    ids=["second_page_full", "third_page_full_after_a_reduced_one", "second_page_is_a_page"],
)
def test_refuses_ambiguous_pages(
    refused: Callable[..., str],
    extra_pages: list[tuple[np.ndarray, int]],
    subfiletype: int,
    index: int,
    count: int,
) -> None:
    """§5 table: `NewSubfileType` (254) of the first offending page with its value
    (0 when absent), that page's index and the page count. Indices are 0-based:
    §5 refusal 8 calls the page that is read "page 0"."""
    message = refused(micro_tiff(extra_pages=extra_pages), "254", "NewSubfileType")
    _names_numbers(message, subfiletype, index, count)


def test_refuses_multi_sample(refused: Callable[..., str]) -> None:
    refuse_named(refused, "refuses_multi_sample")


@pytest.mark.parametrize(
    ("file_dtype", "sample_format", "bits"),
    [
        ("int64", 2, 64),
        ("uint64", 1, 64),
        ("float16", 3, 16),
        ("complex64", 6, 64),
        ("bool", 1, 1),
    ],
)
def test_refuses_unsupported_dtype(
    refused: Callable[..., str], file_dtype: str, sample_format: int, bits: int
) -> None:
    """§5 table: `SampleFormat` (339) and `BitsPerSample` (258) with their values,
    and the numpy dtype. The values are tifffile's, measured in §5."""
    names = ("339", "SampleFormat", "258", "BitsPerSample", file_dtype)
    message = refused(micro_tiff(elevations(file_dtype)), *names)
    _names_numbers(message, sample_format, bits)


@without_codecs
@pytest.mark.parametrize(
    ("build", "names", "value"),
    [
        (lambda: with_compression_tag(micro_tiff(), LZW), ("259", "Compression", "LZW"), LZW),
        (
            floating_point_predictor_tiff,
            ("317", "Predictor", "floating"),
            FLOATING_POINT_PREDICTOR,
        ),
    ],
    ids=["lzw_compression", "floating_point_predictor"],
)
def test_refuses_missing_codec(
    refused: Callable[..., str], build: Callable[[], Any], names: tuple[str, ...], value: int
) -> None:
    """§8: not tifffile's bare error. Names the tag, its value, the scheme and the
    `codecs` extra. "floating" is the scheme, however it is spelled (tifffile's
    enum says `FLOATINGPOINT`). Amended (problem 1): the predictor counts too."""
    message = refused(build(), *names)
    _names_numbers(message, value)
    assert re.search(r"(?<!image)codecs", message), f"extra not named in {message!r}"


@needs_codecs
def test_missing_codec_direction_reads_when_extra_present(decode: Decode) -> None:
    """The other direction of §8: with the extra installed, LZW decodes."""
    tile = decode(micro_tiff(compression="lzw"))
    np.testing.assert_array_equal(tile.array, elevations())


_PROJECTED = ("3072", "ProjectedCSTypeGeoKey")
_GEOGRAPHIC = ("2048", "GeographicTypeGeoKey")


@pytest.mark.parametrize(
    ("changes", "names"),
    [
        ({PROJECTED_CS_TYPE: None}, _PROJECTED),
        ({PROJECTED_CS_TYPE: USER_DEFINED}, (*_PROJECTED, "32767")),
        ({PROJECTED_CS_TYPE: EPSG_UNRESOLVABLE}, (*_PROJECTED, str(EPSG_UNRESOLVABLE))),
        (
            {PROJECTED_CS_TYPE: USER_DEFINED, GEOGRAPHIC_TYPE: EPSG_WGS84},
            (*_PROJECTED, "32767"),
        ),
        (
            {PROJECTED_CS_TYPE: None, GEOGRAPHIC_TYPE: EPSG_UNRESOLVABLE},
            (*_PROJECTED, *_GEOGRAPHIC, str(EPSG_UNRESOLVABLE)),
        ),
        (
            {PROJECTED_CS_TYPE: None, GEOGRAPHIC_TYPE: USER_DEFINED},
            (*_PROJECTED, *_GEOGRAPHIC, "32767"),
        ),
        (
            {PROJECTED_CS_TYPE: None, GEOGRAPHIC_TYPE: EPSG_UTM33},
            (*_PROJECTED, *_GEOGRAPHIC, str(EPSG_UTM33), "Projected CRS"),
        ),
        (
            {PROJECTED_CS_TYPE: None, GEOGRAPHIC_TYPE: EPSG_GEOCENTRIC},
            (*_PROJECTED, *_GEOGRAPHIC, str(EPSG_GEOCENTRIC), "Geocentric CRS"),
        ),
    ],
    ids=[
        "absent",
        "user_defined",
        "unresolvable",
        "user_defined_beside_geographic",
        "absent_geographic_unresolvable",
        "absent_geographic_user_defined",
        "absent_geographic_projected",
        "absent_geographic_geocentric",
    ],
)
def test_refuses_missing_crs(
    refused: Callable[..., str], changes: dict[int, int | None], names: tuple[str, ...]
) -> None:
    """§5 refusal 12, amended (problem 3): neither key yields a resolvable code.

    2048 is consulted only when 3072 is *absent*. So a user-defined 3072 beside
    a geographic 2048 is "no CRS", naming 3072's value, and not a geographic
    refusal; and when 2048 was consulted and failed, the message names it too.

    Amended (round 2): a code reached through 2048 is never accepted. One that
    resolves to a non-geographic CRS is this refusal, and the message names the
    CRS's `type_name`. `absent_geographic_projected` was accepted before round
    2; `test_geotiff_fixtures.py` shows its 2048 really reads back as 25833.
    """
    refused(micro_tiff(geokeys=with_keys(changes)), *names)


def test_refuses_geographic_crs(refused: Callable[..., str]) -> None:
    """The constructed CRS is geographic, so it is refused. Greenfield (§3)."""
    refuse_named(refused, "refuses_geographic_crs")


@pytest.mark.parametrize(
    ("code", "type_name"),
    [(EPSG_WGS84, "Geographic 2D CRS"), (EPSG_GEOCENTRIC, "Geocentric CRS")],
    ids=["geographic", "geocentric"],
)
def test_refuses_geographic_crs_names_the_crs_type(
    refused: Callable[..., str], code: int, type_name: str
) -> None:
    """§5 refusal 13, amended (round 2): the refusal fires whenever the 3072 CRS
    is not projected, which includes geocentric. So the message names the
    constructed CRS's `type_name` (pyproj 3.8.0's spelling) and the code, rather
    than asserting "geographic" of a CRS that is not."""
    keys = with_keys({PROJECTED_CS_TYPE: code})
    refused(micro_tiff(geokeys=keys), *_PROJECTED, str(code), type_name)


def test_a_realistically_encoded_geographic_file_is_refused(refused: Callable[..., str]) -> None:
    """§5 refusal 13, amended (problem 3): ModelTypeGeographic, 2048 = 4326, no
    3072. Refused *as geographic*, naming the key that carried the code."""
    keys = with_keys({GT_MODEL_TYPE: 2, GEOGRAPHIC_TYPE: EPSG_WGS84, PROJECTED_CS_TYPE: None})
    refused(micro_tiff(geokeys=keys), *_GEOGRAPHIC, str(EPSG_WGS84))


@pytest.mark.parametrize(
    ("projected", "linear", "names"),
    [
        (EPSG_UTM33, FOOT, ("3076", "ProjLinearUnitsGeoKey", "9002")),
        (EPSG_FEET, None, (str(EPSG_FEET),)),
        (EPSG_FEET, METRE, (str(EPSG_FEET),)),
    ],
    ids=["key_says_feet", "crs_axes_in_feet", "key_and_crs_disagree"],
)
def test_refuses_non_metre_linear_unit(
    refused: Callable[..., str], projected: int, linear: int | None, names: tuple[str, ...]
) -> None:
    keys = with_keys({PROJECTED_CS_TYPE: projected, PROJ_LINEAR_UNITS: linear})
    refused(micro_tiff(geokeys=keys), *names)


def test_refuses_non_metre_vertical_unit(refused: Callable[..., str]) -> None:
    refuse_named(refused, "refuses_non_metre_vertical_unit")


@pytest.mark.parametrize("text", ["void", "12abc", "1_000"])
def test_refuses_unparseable_nodata(refused: Callable[..., str], text: str) -> None:
    """Measured: tifffile's `page.nodata` logs a warning and returns 0 for 'void'.
    The tag's text is what must be parsed, and refused. `float()` accepts
    '1_000'; no TIFF writer means that, so §6 refuses text containing `_`."""
    refused(micro_tiff(nodata=text), "42113", "GDAL_NODATA", text)


@pytest.mark.parametrize(
    ("file_dtype", "text"),
    [
        ("float32", "0.1"),
        ("float32", "16777217"),
        ("float32", "1e40"),
        ("float32", "1e400"),
        ("float32", "inf"),
        ("float32", "-inf"),
        ("uint8", "0.5"),
        ("uint8", "-9999"),
        ("uint8", "256"),
        ("uint8", "nan"),
    ],
    ids=[
        "float_inexact_fraction",
        "float_beyond_float32_mantissa",
        "float_overflows_float32",
        "float_overflows_every_float",
        "float_infinity",
        "float_negative_infinity",
        "int_fraction",
        "int_below_range",
        "int_above_range",
        "int_nan",
    ],
)
def test_refuses_nodata_not_representable(
    refused: Callable[..., str], file_dtype: str, text: str
) -> None:
    """§6, amended (problem 5): no cell of the *file* dtype can equal the value, so
    `v == *nodata_` would disable NoData silently. On uint8, `0.5` and `-9999`
    pass a float32 round trip and still match no cell (5c). `1e400` parses to
    inf, and infinities are refused (5d). The message names the file dtype."""
    source = elevations(file_dtype)
    refused(micro_tiff(source, nodata=text), "42113", "GDAL_NODATA", text, file_dtype)


@pytest.mark.parametrize(
    ("file_dtype", "caller", "value"),
    [
        ("uint8", 0.5, "0.5"),
        ("uint8", -9999.0, "-9999"),
        ("uint8", math.nan, "nan"),
        ("float32", 0.1, "0.1"),
        ("float32", math.inf, "inf"),
    ],
    ids=["int_fraction", "int_below_range", "int_nan", "float_inexact", "float_infinity"],
)
def test_refuses_caller_nodata_not_representable(
    refused: Callable[..., str], file_dtype: str, caller: float, value: str
) -> None:
    """§6 (5b): the caller's sentinel passes the same check, under the same name.
    The message names the `nodata=` argument, not tag 42113."""
    stream = micro_tiff(elevations(file_dtype))
    message = refused(stream, "nodata", value, file_dtype, nodata=caller)
    assert "42113" not in message, message


def test_representability_is_checked_before_contradiction(refused: Callable[..., str]) -> None:
    """§6: the contradiction check runs after both values pass the representability
    check. A valid tag and an unrepresentable caller value is therefore the
    caller's refusal, which does not name tag 42113, and not a contradiction."""
    message = refused(micro_tiff(nodata="-32767"), "nodata", "0.1", "float32", nodata=0.1)
    assert "42113" not in message, message


@pytest.mark.parametrize(
    ("tag", "caller", "names"),
    [
        ("-32767", -9999.0, REFUSAL["refuses_contradictory_nodata_override"].must_name),
        ("-32767", math.nan, ("42113", "GDAL_NODATA", "-32767", "nan")),
        ("nan", -32767.0, ("42113", "GDAL_NODATA", "nan", "-32767")),
    ],
    ids=["finite_both", "caller_nan_tag_finite", "caller_finite_tag_nan"],
)
def test_refuses_contradictory_nodata_override(
    refused: Callable[..., str], tag: str, caller: float, names: tuple[str, ...]
) -> None:
    """§6: "equal" treats NaN as equal to NaN and to nothing else (round 2). With
    `test_caller_nan_agreeing_with_nan_tag_is_accepted` this pins both sides.
    The first case is the catalogue's fixture, with the catalogue's names."""
    refused(micro_tiff(nodata=tag), *names, nodata=caller)


# ---------------------------------------------------------------------------
# Refusals are real exceptions (§5; prior art §5.9)
# ---------------------------------------------------------------------------

_UNDER_O = """
import sys
from geotiff_fixtures import REFUSALS
from tin_engine.io.geotiff import decode_dem
from tin_engine.io.models import GeoTiffError
passed = []
for r in REFUSALS:
    try:
        decode_dem(r.build(), **r.decode_kwargs)
    except GeoTiffError:
        passed.append(r.name)
expected = {r.name for r in REFUSALS} - set(sys.argv[1:])
missed = sorted(expected - set(passed))
print("MISSED", missed)
sys.exit(1 if missed else 0)
"""


def test_every_refusal_survives_python_dash_o() -> None:
    """Under `-O` a bare `assert` vanishes. Every catalogue refusal must still raise.

    The codec refusal is excluded when the extra is installed: its fixture then
    reaches a real LZW decoder and is not a refusal case at all.
    """
    excluded = ["refuses_missing_codec"] if HAS_CODECS else []
    result = subprocess.run(
        [sys.executable, "-O", "-c", _UNDER_O, *excluded],
        capture_output=True,
        text=True,
        cwd=HERE,
        env={"PYTHONPATH": str(HERE), "PATH": ""},
        check=False,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "MISSED []" in result.stdout


# ---------------------------------------------------------------------------
# Ruling 1: decode only, never `_core`
# ---------------------------------------------------------------------------


def _imported_modules(source: Path) -> set[str]:
    tree = ast.parse(source.read_text(encoding="utf-8"))
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            prefix = "." * node.level + (node.module or "")
            found.add(prefix)
            found.update(f"{prefix}.{alias.name}".replace("..", ".") for alias in node.names)
    return found


@pytest.mark.parametrize("module", ["geotiff.py", "models.py"])
def test_module_never_imports_core(module: str) -> None:
    """§2: tifffile, pyproj, numpy, pydantic and stdlib; first-party only `io.models`."""
    source = HERE.parents[1] / "src_python" / "tin_engine" / "io" / module
    imported = _imported_modules(source)
    assert not any("_core" in name for name in imported), sorted(imported)
    first_party = {n for n in imported if n.startswith(("tin_engine", "."))}
    allowed = {".models", ".models.", "tin_engine.io.models"}
    assert all(any(n.startswith(a) for a in allowed) for n in first_party), sorted(first_party)


# ---------------------------------------------------------------------------
# The real Kartverket fixture (§12). LZW, so both need the `codecs` extra.
# ---------------------------------------------------------------------------


@needs_codecs
def test_kartverket_fixture_decodes(decode: Decode) -> None:
    with KARTVERKET.open("rb") as stream:
        tile = decode(stream)
    assert tile.array.shape == (5051, 5051)
    assert tile.array.dtype == np.float32
    assert not tile.array.flags.writeable
    assert (tile.meta.rows, tile.meta.cols) == (5051, 5051)
    assert (tile.meta.delta_x, tile.meta.delta_y) == (10.0, 10.0)
    assert tile.meta.epsg == 25833
    assert tile.meta.nodata == -32767.0
    assert tile.meta.nodata_source == "tag"
    assert tile.meta.pixel_is_area is True
    assert tile.meta.vertical_unit_assumed is False


@needs_codecs
def test_kartverket_fixture_is_shifted_half_a_cell(decode: Decode) -> None:
    """§4: tie point (799745, 7950255) at 10 m, area-registered. The only check
    that ruling 4 is implemented and not merely documented."""
    geo = tifffile.TiffFile(KARTVERKET).geotiff_metadata
    assert geo is not None
    assert geo["ModelTiepoint"][3:5] == [799745.0, 7950255.0]
    with KARTVERKET.open("rb") as stream:
        meta = decode(stream).meta
    assert meta.x_min == 799750.0
    assert meta.y_max == 7950250.0
