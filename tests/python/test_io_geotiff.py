"""`tin_engine.io.geotiff.decode_dem`: increment 11, a GeoTIFF stream to a `DemTile`.

The design is `docs/increments/11-raster-ingestion.md`; section numbers below
are its. The reader's job is refusal (ruling 3), so most of this file is the
nineteen named refusals of §5 (eighteen, and round 3's 13a), one crafted micro-TIFF each, built by
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
import functools
import importlib.util
import math
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pyproj
import pytest
import tifffile

from geotiff_fixtures import (
    COLS,
    DELTA_X,
    DELTA_Y,
    EPSG_COMPOUND,
    EPSG_COMPOUND_GEOGRAPHIC,
    EPSG_COMPOUND_HORIZONTAL,
    EPSG_FEET,
    EPSG_GEOCENTRIC,
    EPSG_PROJECTED_3D,
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
    SCALE,
    THUNDERSCAN,
    TIE_X,
    TIE_Y,
    TIEPOINT,
    UNDECODABLE,
    UNDECODABLE_WITH_CODECS,
    UNKNOWN_COMPRESSION,
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

    def test_tile_does_not_share_the_callers_buffer(self, decode: Decode) -> None:
        """§7, amended (round 3): `DemTile` copies. A view of an already
        C-contiguous caller array shared its buffer, so a write by the caller
        changed the tile that every refinement thread samples."""
        from tin_engine.io.models import DemTile

        mine = elevations(np.float32)
        assert mine.flags.writeable and mine.flags.c_contiguous
        tile = DemTile(meta=decode(micro_tiff()).meta, array=mine)
        assert not np.shares_memory(mine, tile.array)
        mine[0, 0] = 42.0
        np.testing.assert_array_equal(tile.array, elevations(np.float32))

    @pytest.mark.parametrize("origin", ["decoded", "caller_built"])
    def test_tile_array_cannot_be_made_writeable(self, decode: Decode, origin: str) -> None:
        """§7, amended (round 3): the tile never hands out a writeable handle.
        A read-only view whose base is writeable can have its flag set back;
        a view of a read-only owned copy cannot (numpy: "cannot set WRITEABLE
        flag to True of this array")."""
        from tin_engine.io.models import DemTile

        tile = decode(micro_tiff())
        if origin == "caller_built":
            tile = DemTile(meta=tile.meta, array=elevations(np.float32))
        with pytest.raises(ValueError, match="WRITEABLE"):
            tile.array.flags.writeable = True
        assert not tile.array.flags.writeable


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

    def test_meta_refuses_absent_source_with_a_sentinel(
        self, decode: Decode, models: Any, geotiff_error: type[Exception]
    ) -> None:
        """§6: "`nodata_source == "absent"` means `nodata is None`". A finite
        sentinel is used so `allow_inf_nan` cannot be what refuses it, and the
        error must come from the model validator (an empty `loc`), not a field."""
        from pydantic import ValidationError

        with pytest.raises(ValidationError, match="absent") as info:
            self._meta(decode, models, nodata=-9999.0, nodata_source="absent")
        assert not isinstance(info.value, geotiff_error)
        assert [e["loc"] for e in info.value.errors()] == [()]

    @pytest.mark.parametrize("source", ["absent", "tag", "caller"])
    def test_meta_accepts_no_sentinel_from_any_source(
        self, decode: Decode, models: Any, source: str
    ) -> None:
        """§6: the implication runs one way only. `None` with "tag" or "caller"
        is the NaN rows of the table, and "absent" with `None` is the first row."""
        meta = self._meta(decode, models, nodata=None, nodata_source=source)
        assert (meta.nodata, meta.nodata_source) == (None, source)

    @pytest.mark.parametrize(
        "array",
        [
            elevations(np.float32).reshape(ROWS, COLS, 1),
            elevations(np.float64)[np.newaxis, :, :],
            elevations(np.int16),
            elevations(np.float16),
            elevations(np.bool_),
            elevations(np.complex128),
        ],
        ids=["3d_trailing", "3d_leading", "int16", "float16", "bool", "complex128"],
    )
    def test_tile_refuses_array_not_2d_float32_or_float64(
        self, decode: Decode, models: Any, geotiff_error: type[Exception], array: Any
    ) -> None:
        """§7: `array` is "numpy 2-D, C-contiguous, float32 or float64". The
        refusal must come from the `array` field: a 3-D array would also fail
        the shape agreement, so checking the `loc` is what makes the ndim cases
        test this check rather than that one."""
        from pydantic import ValidationError

        with pytest.raises(ValidationError, match="2-D float32 or float64") as info:
            models.DemTile(meta=self._meta(decode, models), array=array)
        assert not isinstance(info.value, geotiff_error)
        assert [e["loc"] for e in info.value.errors()] == [("array",)]

    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    def test_tile_accepts_2d_float32_and_float64(
        self, decode: Decode, models: Any, dtype: Any
    ) -> None:
        """The control for the refusal above: both permitted dtypes build."""
        tile = models.DemTile(meta=self._meta(decode, models), array=elevations(dtype))
        assert tile.array.dtype == dtype


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

    @pytest.mark.parametrize(
        "flag", [True, False, np.True_, np.False_], ids=["true", "false", "np_true", "np_false"]
    )
    def test_caller_bool_is_refused(self, decode: Decode, flag: Any) -> None:
        """§6, amended (round 3): a `bool` is not a sentinel. `False` became 0.0
        and deleted every sea-level cell. It is `TypeError`, not `GeoTiffError`:
        a programming error in the caller, not a fact about a file (§7, B1).
        The message names the `nodata=` argument and the type."""
        with pytest.raises(TypeError) as info:
            decode(micro_tiff(), nodata=flag)
        message = str(info.value)
        assert "nodata" in message, message
        assert "bool" in message.lower(), message


# ---------------------------------------------------------------------------
# The nineteen refusals (§5, with round 3's 13a). Test names are the design's names.
# ---------------------------------------------------------------------------


def _tiepoint(**changes: float) -> tuple[float, ...]:
    """The baseline tie point with some of I, J, X, Y replaced."""
    values = {"i": 0.0, "j": 0.0, "x": TIE_X, "y": TIE_Y, **changes}
    return (values["i"], values["j"], 0.0, values["x"], values["y"], 0.0)


_TIEPOINT_NAMES = ("33922", "ModelTiepointTag")
_SCALE_NAMES = ("33550", "ModelPixelScaleTag")


@pytest.mark.parametrize(
    ("changes", "names", "count"),
    [
        ({"tiepoint": None}, _TIEPOINT_NAMES, None),
        ({"scale": None}, _SCALE_NAMES, None),
        ({"tiepoint": _tiepoint(i=math.nan)}, (*_TIEPOINT_NAMES, "nan"), None),
        ({"tiepoint": _tiepoint(j=math.inf)}, (*_TIEPOINT_NAMES, "inf"), None),
        ({"tiepoint": _tiepoint(x=math.nan)}, (*_TIEPOINT_NAMES, "nan"), None),
        ({"tiepoint": _tiepoint(y=-math.inf)}, (*_TIEPOINT_NAMES, "-inf"), None),
        ({"tiepoint": TIEPOINT[:1]}, _TIEPOINT_NAMES, 1),
        ({"tiepoint": TIEPOINT[:3]}, _TIEPOINT_NAMES, 3),
        ({"tiepoint": TIEPOINT[:5]}, _TIEPOINT_NAMES, 5),
        ({"scale": SCALE[:1]}, _SCALE_NAMES, 1),
        ({"scale": SCALE[:2]}, _SCALE_NAMES, 2),
    ],
    ids=[
        "tiepoint_absent",
        "scale_absent",
        "i_nan",
        "j_inf",
        "x_nan",
        "y_negative_inf",
        "tiepoint_1_value",
        "tiepoint_3_values",
        "tiepoint_5_values",
        "scale_1_value",
        "scale_2_values",
    ],
)
def test_refuses_missing_georeferencing(
    refused: Callable[..., str],
    changes: dict[str, Any],
    names: tuple[str, ...],
    count: int | None,
) -> None:
    """§5 refusal 1, amended (problem 7): a non-finite origin is as unusable as none.

    Amended (round 3): a short tie point, or a scale with any count but 3, is
    refusal 1 too, and the count is the diagnosis, so the message names it.
    These are checked from `page.tags` before `geotiff_metadata`, which would
    otherwise raise tifffile's bare `ValueError` on the reshape. A count of 1
    is the scalar path: a one-value DOUBLE tag reads back as a bare `float`
    (`test_geotiff_fixtures.py` shows it), which the green code iterated.
    """
    message = refused(micro_tiff(**changes), *names)
    if count is not None:
        _names_numbers(message, count)


_TRANSFORMATION_NAMES = ("34264", "ModelTransformationTag")


@pytest.mark.parametrize(
    "build",
    [
        REFUSAL["refuses_model_transformation"].build,
        lambda: micro_tiff(tiepoint=None, scale=None, transformation=(1.0, 0.0, 0.0)),
        lambda: micro_tiff(tiepoint=None, scale=None, transformation=tuple(map(float, range(12)))),
    ],
    ids=["16_values_beside_tiepoint_and_scale", "3_values_alone", "12_values_alone"],
)
def test_refuses_model_transformation(
    refused: Callable[..., str], build: Callable[[], Any]
) -> None:
    """§5 refusal 2: refused even though a tie point and scale are also present.

    Amended (round 3): with any length. The 3- and 12-value cases carry no tie
    point and no scale, so a refusal naming 34264 shows that 34264 is checked
    before either is looked at (and before `geotiff_metadata`, whose 4x4
    reshape used to raise first)."""
    refused(build(), *_TRANSFORMATION_NAMES)


@pytest.mark.parametrize(
    ("build", "count"),
    [
        (lambda: micro_tiff(tiepoint=(*TIEPOINT, 1.0)), 7),
        (lambda: micro_tiff(tiepoint=(*TIEPOINT, 1.0, 1.0, 0.0)), 9),
        (REFUSAL["refuses_multiple_tiepoints"].build, 12),
    ],
    ids=["7_values", "9_values", "12_values"],
)
def test_refuses_multiple_tiepoints(
    refused: Callable[..., str], build: Callable[[], Any], count: int
) -> None:
    """§5 refusal 3, amended (round 3): more than 6 values, whatever the count.
    7 and 9 are not a GCP list, so the message gives the count instead. 12 is
    the catalogue's fixture."""
    message = refused(build(), *_TIEPOINT_NAMES)
    _names_numbers(message, count)


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
    ("extra_pages", "subfiletype", "index", "count", "frames"),
    [
        ([(elevations(), 0)], 0, 1, 2, False),
        ([(elevations(rows=2, cols=2), 1), (elevations(), 0)], 0, 2, 3, False),
        ([(elevations(), 2)], 2, 1, 2, False),
        ([(elevations(), 0)], None, 1, 2, True),
    ],
    ids=[
        "second_page_full",
        "third_page_full_after_a_reduced_one",
        "second_page_is_a_page",
        "unparsed_frame",
    ],
)
def test_refuses_ambiguous_pages(
    refused: Callable[..., str],
    monkeypatch: pytest.MonkeyPatch,
    extra_pages: list[tuple[np.ndarray, int]],
    subfiletype: int | None,
    index: int,
    count: int,
    frames: bool,
) -> None:
    """§5 table: `NewSubfileType` (254) of the first offending page with its value
    (0 when absent), that page's index and the page count. Indices are 0-based:
    §5 refusal 8 calls the page that is read "page 0".

    Amended (round 3), `unparsed_frame`: an extra page tifffile returns as a
    `TiffFrame` is refused, because a frame has no `is_reduced` to consult.
    The message names 254, the index and the count, and no value (there is
    none to read). The page is really full resolution, so the green code's
    skip dropped a second image without a word. The test reaches for
    tifffile's private `_useframes=True` because no public route yields a
    frame: tifffile makes them only for LSM, NDPI and ScanImage files, never
    a DEM. `test_geotiff_fixtures.py` shows the patch does produce a frame.
    """
    if frames:
        framed = functools.partial(tifffile.TiffFile, _useframes=True)
        monkeypatch.setattr(tifffile, "TiffFile", framed)
    message = refused(micro_tiff(extra_pages=extra_pages), "254", "NewSubfileType")
    _names_numbers(message, *(n for n in (subfiletype, index, count) if n is not None))


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
        (
            lambda: with_compression_tag(micro_tiff(), UNKNOWN_COMPRESSION),
            ("259", "Compression"),
            UNKNOWN_COMPRESSION,
        ),
    ],
    ids=["lzw_compression", "floating_point_predictor", "unknown_compression"],
)
def test_refuses_missing_codec(
    refused: Callable[..., str], build: Callable[[], Any], names: tuple[str, ...], value: int
) -> None:
    """§8: not tifffile's bare error. Names the tag, its value, the scheme and the
    `codecs` extra. "floating" is the scheme, however it is spelled (tifffile's
    enum says `FLOATINGPOINT`). Amended (problem 1): the predictor counts too.

    Amended (round 3), `unknown_compression`: a number tifffile's enum does not
    know has no scheme name, and the message prints `Compression (259) = 60000`
    with no "(unknown)" after it."""
    message = refused(build(), *names)
    _names_numbers(message, value)
    assert re.search(r"(?<!image)codecs", message), f"extra not named in {message!r}"
    assert "(unknown)" not in message.lower(), message


@needs_codecs
@pytest.mark.parametrize(
    ("code", "names"),
    [(THUNDERSCAN, ("THUNDERSCAN",)), (UNKNOWN_COMPRESSION, ())],
    ids=["thunderscan", "unknown_compression"],
)
def test_undecodable_scheme_with_extra_present_gives_no_install_advice(
    refused: Callable[..., str], code: int, names: tuple[str, ...]
) -> None:
    """§5 refusal 11, amended (round 3): with imagecodecs installed, THUNDERSCAN
    (32809) is still not decodable, so advising the user to install the extra
    would be false. The message names the tag, value and scheme, and no
    `pip install`. Runs in CI's second pytest step, which installs the extra (§12).

    `unknown_compression` is the same rule for a number with no scheme name,
    and carries the no-"(unknown)" check into the branch where the extra is
    present; `test_refuses_missing_codec` only runs where it is absent."""
    stream = with_compression_tag(micro_tiff(), code)
    message = refused(stream, "259", "Compression", *names)
    _names_numbers(message, code)
    assert "pip install" not in message.lower(), message
    assert "(unknown)" not in message.lower(), message


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
    [
        (EPSG_WGS84, "Geographic 2D CRS"),
        (EPSG_GEOCENTRIC, "Geocentric CRS"),
        (EPSG_COMPOUND_GEOGRAPHIC, "Compound CRS"),
    ],
    ids=["geographic", "geocentric", "compound_geographic"],
)
def test_refuses_geographic_crs_names_the_crs_type(
    refused: Callable[..., str], code: int, type_name: str
) -> None:
    """§5 refusal 13, amended (round 2): the refusal fires whenever the 3072 CRS
    is not projected, which includes geocentric. So the message names the
    constructed CRS's `type_name` (pyproj 3.8.0's spelling) and the code, rather
    than asserting "geographic" of a CRS that is not.

    Amended (round 3), `compound_geographic`: 9707 (WGS 84 + EGM96 height) is
    compound, but not projected, so it stays refusal 13. With the geocentric
    case (3 axes too), this pins that refusal 13a runs after the
    `is_projected` test. Its message says "not a projected CRS", which 13a's
    does not; that is what tells the two apart."""
    keys = with_keys({PROJECTED_CS_TYPE: code})
    message = refused(micro_tiff(geokeys=keys), *_PROJECTED, str(code), type_name)
    assert "4096" not in message, f"refusal 13a fired before refusal 13: {message!r}"


@pytest.mark.parametrize(
    ("code", "names"),
    [
        (
            EPSG_COMPOUND,
            ("Compound CRS", str(EPSG_COMPOUND_HORIZONTAL), "4096"),
        ),
        (EPSG_PROJECTED_3D, ("Projected CRS", "4096")),
    ],
    ids=["compound_projected", "projected_3d"],
)
def test_refuses_compound_crs(
    refused: Callable[..., str], code: int, names: tuple[str, ...]
) -> None:
    """§5 refusal 13a, added in round 3 by user ruling: the 3072 CRS passes
    `is_projected` but is not a 2-D horizontal CRS. Both used to be accepted.

    5972 is compound, and pyproj reports it as projected because its
    horizontal part is. 9895 is a "Projected CRS" with an ellipsoidal-height
    axis: not compound, so only the axis count catches it. The message names
    3072 and the code, the CRS's `type_name`, the axis count (3), the
    horizontal sub-CRS's code when there is one (11022 for 5972, measured, not
    25832), and `VerticalGeoKey (4096)` as where the vertical part belongs."""
    keys = with_keys({PROJECTED_CS_TYPE: code})
    message = refused(micro_tiff(geokeys=keys), *_PROJECTED, str(code), *names)
    _names_numbers(message, 3)


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
# Errors from inside tifffile or a codec (§14, choice C)
# ---------------------------------------------------------------------------

_UNDECODABLE_CASES = [
    *(pytest.param(build, stage, id=name) for name, (build, stage) in UNDECODABLE.items()),
    *(
        pytest.param(build, stage, id=name, marks=needs_codecs)
        for name, (build, stage) in UNDECODABLE_WITH_CODECS.items()
    ),
]


@pytest.mark.parametrize(("build", "stage"), _UNDECODABLE_CASES)
def test_undecodable_input_is_a_geotiff_error(
    decode: Decode, geotiff_error: type[Exception], build: Callable[[], Any], stage: str
) -> None:
    """§14, choice C: the three calls into tifffile are wrapped. What escaped
    before was `TiffFileError`, `struct.error`, a bare `ValueError`, `zlib.error`
    or an imagecodecs `RuntimeError`, and the type changed with the extra.

    Now it is `GeoTiffError`, chained `from` the original, and the message
    names the stage ("TIFF structure" or "pixel data"), the original type and
    the original text. `test_geotiff_fixtures.py` shows tifffile alone fails
    on every stream here."""
    with pytest.raises(geotiff_error) as info:
        decode(build())
    cause = info.value.__cause__
    assert cause is not None, "the original error is not chained"
    assert not isinstance(cause, geotiff_error), cause
    message = str(info.value)
    assert stage.lower() in message.lower(), f"stage {stage!r} not named in {message!r}"
    assert type(cause).__name__ in message, f"{type(cause).__name__} not named in {message!r}"
    assert str(cause) in message, f"{str(cause)!r} not in {message!r}"


def test_geokey_directory_failure_is_a_geotiff_error(
    decode: Decode, geotiff_error: type[Exception], monkeypatch: pytest.MonkeyPatch
) -> None:
    """§14, choice C: the third wrapped call. No crafted file makes
    `geotiff_metadata` raise once round 3's tag checks run first (§5: a
    corrupt 34735 is logged, not raised), so the failure is planted."""

    def broken(self: Any) -> Any:
        raise RuntimeError("planted GeoKey failure")

    monkeypatch.setattr(tifffile.TiffFile, "geotiff_metadata", property(broken))
    with pytest.raises(geotiff_error) as info:
        decode(micro_tiff())
    assert isinstance(info.value.__cause__, RuntimeError)
    message = str(info.value)
    for token in ("GeoKey directory", "RuntimeError", "planted GeoKey failure"):
        assert token.lower() in message.lower(), f"{token!r} not named in {message!r}"


def test_memory_error_is_not_wrapped(decode: Decode, monkeypatch: pytest.MonkeyPatch) -> None:
    """§14, choice C: `MemoryError` passes through untouched. It says nothing
    about the file. It guards the pixel-data stage from swallowing it, which a
    bare `except Exception` would not (`MemoryError` is an `Exception`)."""

    def exhausted(self: Any, *args: Any, **kwargs: Any) -> Any:
        raise MemoryError("planted")

    monkeypatch.setattr(tifffile.TiffPage, "asarray", exhausted)
    with pytest.raises(MemoryError, match="planted"):
        decode(micro_tiff())


def test_reader_bug_is_not_wrapped(decode: Decode, monkeypatch: pytest.MonkeyPatch) -> None:
    """§14, choice C, and why not A: only the calls into tifffile are wrapped,
    so the reader's own defects surface as themselves. The round-3 `TypeError`
    for a one-value 33550 is the example: wrapped, it would have read as a file
    refusal and never been filed as a bug.

    The planted bug is in CRS resolution, which runs inside the `with` block
    between `geotiff_metadata` and `asarray`, and outside every wrapped call.
    It is planted on `pyproj.CRS.from_epsg` rather than on a private helper,
    because §5 splits `_placement` in round 3 and the successor's name is not
    fixed; §5 does fix that the code is resolved with `from_epsg`. It is what
    fails under option A, one wrap around the whole body."""

    def buggy(code: Any) -> Any:
        raise TypeError("planted reader bug")

    monkeypatch.setattr(pyproj.CRS, "from_epsg", staticmethod(buggy))
    with pytest.raises(TypeError, match="planted reader bug") as info:
        decode(micro_tiff())
    assert info.value.__cause__ is None


def test_page_refusal_bug_is_not_wrapped(decode: Decode, monkeypatch: pytest.MonkeyPatch) -> None:
    """§14, choice C: the "TIFF structure" stage covers opening the file and
    walking its pages, and not the reader's refusal logic over those pages
    (§5 refusal 8). A bug in that logic must surface as itself.

    The bug is planted on tifffile's public `TiffPage.is_reduced`, which the
    refusal consults for each extra page (§5: a reduced page is allowed).
    Collecting the pages does not read it, so the plant fires only in the
    reader's own check, whatever that check's helper is called. The file is
    valid, with a reduced second page, so nothing else can raise.

    Red while the refusal loop runs inside the stage: the `TypeError` comes out
    as `GeoTiffError: TIFF structure: ... TypeError`."""

    def buggy(self: Any) -> Any:
        raise TypeError("planted refusal bug")

    stream = micro_tiff(extra_pages=[(elevations(rows=2, cols=2), 1)])
    monkeypatch.setattr(tifffile.TiffPage, "is_reduced", property(buggy))
    with pytest.raises(TypeError, match="planted refusal bug") as info:
        decode(stream)
    assert info.value.__cause__ is None


# ---------------------------------------------------------------------------
# Refusals are real exceptions (§5; prior art §5.9)
# ---------------------------------------------------------------------------

_UNDER_O = """
import sys
from geotiff_fixtures import REFUSALS
from tin_engine.io.geotiff import decode_dem
from tin_engine.io.models import GeoTiffError
excluded = set(sys.argv[1:])
unknown = excluded - {r.name for r in REFUSALS}
if unknown:
    sys.exit(f"UNKNOWN EXCLUSION {sorted(unknown)}")
passed = []
expected = [r for r in REFUSALS if r.name not in excluded]
for r in expected:
    try:
        decode_dem(r.build(), **r.decode_kwargs)
    except GeoTiffError:
        passed.append(r.name)
missed = sorted({r.name for r in expected} - set(passed))
print("CHECKED", len(expected))
print("MISSED", missed)
sys.exit(1 if missed else 0)
"""


def test_every_refusal_survives_python_dash_o() -> None:
    """Under `-O` a bare `assert` vanishes. Every catalogue refusal must still raise.

    The codec refusal is excluded when the extra is installed: its fixture then
    reaches a real LZW decoder and is not a refusal case at all. Excluded means
    not decoded: its fake LZW strip would make the real decoder raise its own
    error, which is not a `GeoTiffError` and would crash the script. An
    exclusion naming no catalogue entry fails, so a rename cannot silently
    widen what is skipped.
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
    assert f"CHECKED {len(REFUSALS) - len(excluded)}" in result.stdout, result.stdout


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
    """§4: tie point (799745, 7950255) at 10 m, area-registered. A local check
    that the micro-TIFF result holds on the one real product; ruling 4 itself is
    checked by `TestPlacement.test_area_registered_nodes_are_shifted_inward_half_a_cell`,
    which runs everywhere (§12, round 3)."""
    geo = tifffile.TiffFile(KARTVERKET).geotiff_metadata
    assert geo is not None
    assert geo["ModelTiepoint"][3:5] == [799745.0, 7950255.0]
    with KARTVERKET.open("rb") as stream:
        meta = decode(stream).meta
    assert meta.x_min == 799750.0
    assert meta.y_max == 7950250.0
