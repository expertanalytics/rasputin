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
*what* is named and not the sentence around it. Four refusals have no tag to
name (degenerate shape, ambiguous pages, unsupported dtype, missing codec);
for those the suite asserts only what §5 says each one names.
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
    EPSG_UNRESOLVABLE,
    EPSG_UTM33,
    EPSG_WGS84,
    FOOT,
    GEOG_LINEAR_UNITS,
    GEOGRAPHIC_TYPE,
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
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
    micro_tiff,
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

    @pytest.mark.parametrize("text", ["nan", "NaN"])
    def test_nan_tag_yields_no_sentinel(self, decode: Decode, text: str) -> None:
        """§6: a NaN sentinel matches nothing under `==`; `is_nodata` catches NaN."""
        assert decode(micro_tiff(nodata=text)).meta.nodata is None

    def test_caller_sentinel_without_tag(self, decode: Decode) -> None:
        meta = decode(micro_tiff(), nodata=-9999.0).meta
        assert meta.nodata == -9999.0
        assert meta.nodata_source == "caller"

    def test_caller_sentinel_agreeing_with_tag_is_accepted(self, decode: Decode) -> None:
        assert decode(micro_tiff(nodata="-32767"), nodata=-32767.0).meta.nodata == -32767.0


# ---------------------------------------------------------------------------
# The eighteen refusals (§5). Test names are the design's names.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("omit", "names"),
    [("tiepoint", ("33922", "ModelTiepointTag")), ("scale", ("33550", "ModelPixelScaleTag"))],
)
def test_refuses_missing_georeferencing(
    refused: Callable[..., str], omit: str, names: tuple[str, ...]
) -> None:
    refused(micro_tiff(**{omit: None}), *names)


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


@pytest.mark.parametrize(("rows", "cols"), [(1, COLS), (ROWS, 1), (1, 1)])
def test_refuses_degenerate_shape(refused: Callable[..., str], rows: int, cols: int) -> None:
    refused(micro_tiff(elevations(rows=rows, cols=cols)))


def test_refuses_ambiguous_pages(refused: Callable[..., str]) -> None:
    refuse_named(refused, "refuses_ambiguous_pages")


def test_refuses_multi_sample(refused: Callable[..., str]) -> None:
    refuse_named(refused, "refuses_multi_sample")


@pytest.mark.parametrize("file_dtype", ["int64", "uint64", "float16", "complex64", "bool"])
def test_refuses_unsupported_dtype(refused: Callable[..., str], file_dtype: str) -> None:
    refused(micro_tiff(elevations(file_dtype)), file_dtype)


@without_codecs
def test_refuses_missing_codec(refused: Callable[..., str]) -> None:
    """§8: not tifffile's bare error. Names the scheme and the `codecs` extra."""
    message = refuse_named(refused, "refuses_missing_codec")
    assert re.search(r"(?<!image)codecs", message), f"extra not named in {message!r}"


@needs_codecs
def test_missing_codec_direction_reads_when_extra_present(decode: Decode) -> None:
    """The other direction of §8: with the extra installed, LZW decodes."""
    tile = decode(micro_tiff(compression="lzw"))
    np.testing.assert_array_equal(tile.array, elevations())


@pytest.mark.parametrize(
    ("projected", "names"),
    [
        (None, ("3072", "ProjectedCSTypeGeoKey")),
        (USER_DEFINED, ("3072", "ProjectedCSTypeGeoKey", "32767")),
        (EPSG_UNRESOLVABLE, ("3072", "ProjectedCSTypeGeoKey", str(EPSG_UNRESOLVABLE))),
    ],
    ids=["absent", "user_defined", "unresolvable"],
)
def test_refuses_missing_crs(
    refused: Callable[..., str], projected: int | None, names: tuple[str, ...]
) -> None:
    refused(micro_tiff(geokeys=with_keys({PROJECTED_CS_TYPE: projected})), *names)


def test_refuses_geographic_crs(refused: Callable[..., str]) -> None:
    """The constructed CRS is geographic, so it is refused. Greenfield (§3)."""
    refuse_named(refused, "refuses_geographic_crs")


def test_a_realistically_encoded_geographic_file_is_refused(refused: Callable[..., str]) -> None:
    """ModelTypeGeographic, 2048 = 4326, no 3072: whichever refusal fires, it must be one."""
    keys = with_keys({GT_MODEL_TYPE: 2, GEOGRAPHIC_TYPE: EPSG_WGS84, PROJECTED_CS_TYPE: None})
    refused(micro_tiff(geokeys=keys))


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


@pytest.mark.parametrize("text", ["void", "12abc"])
def test_refuses_unparseable_nodata(refused: Callable[..., str], text: str) -> None:
    """Measured: tifffile's `page.nodata` logs a warning and returns 0 for 'void'.
    The tag's text is what must be parsed, and refused."""
    refused(micro_tiff(nodata=text), "42113", "GDAL_NODATA", text)


@pytest.mark.parametrize(
    "text",
    ["0.1", "16777217", "1e40"],
    ids=["inexact_fraction", "beyond_float32_mantissa", "overflows_float32"],
)
def test_refuses_nodata_not_representable(refused: Callable[..., str], text: str) -> None:
    """Float32 tile: `v == *nodata_` would match no cell, disabling NoData silently."""
    refused(micro_tiff(nodata=text), "42113", "GDAL_NODATA", text)


def test_refuses_contradictory_nodata_override(refused: Callable[..., str]) -> None:
    refuse_named(refused, "refuses_contradictory_nodata_override")


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
