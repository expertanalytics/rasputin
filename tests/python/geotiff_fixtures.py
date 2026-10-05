"""Micro-GeoTIFF builders for increment 11's suite (`11-raster-ingestion.md` §12).

Every stream here is an `io.BytesIO` written by `tifffile`, 4x4 or smaller,
uncompressed or Deflate. No optional dependency is needed to build one. Three
builders patch the written bytes afterwards, because tifffile cannot *write*
what they need without `imagecodecs`: an LZW `Compression`, a floating-point
`Predictor`, and a PackBits strip (design §8).

The valid baseline is deliberately asymmetric so a swap cannot pass: 3 rows by
4 columns, `delta_x` 10 and `delta_y` 5, every cell a distinct value. Each
refusal fixture is that baseline with **exactly one** defect, so a refusal can
only fire for the reason its name gives. `test_geotiff_fixtures.py` checks the
defect is present using tifffile alone, with no production code involved.

GeoKeys are written as a raw `GeoKeyDirectoryTag` (34735), by
`geographic_fixtures.geokey_tags`: an `int` value is an inline SHORT, a
`float` or a tuple of floats goes to `GeoDoubleParamsTag` (34736) and a `str`
to `GeoAsciiParamsTag` (34737). That is the on-disk form tifffile's
`geotiff_metadata` decodes, so the reader under test sees what a real file
would give it. A directory of SHORTs alone writes neither params tag, so such
a file is byte for byte what it was when only SHORTs could be written.
"""

from __future__ import annotations

import importlib.util
import io
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import tifffile

from geographic_fixtures import geokey_tags

#: The one real DEM in the tree, and the markers for what decoding it needs.
#: Shared here since increment 12, whose CLI suite is the second user.
KARTVERKET = Path(__file__).resolve().parents[1] / "fixtures" / "dem_archive" / "7908_3_10m_z33.tif"
HAS_CODECS = importlib.util.find_spec("imagecodecs") is not None
needs_codecs = pytest.mark.skipif(not HAS_CODECS, reason="LZW needs the `codecs` extra (§8)")
without_codecs = pytest.mark.skipif(HAS_CODECS, reason="tests the refusal when the extra is absent")

# Tag numbers, named the way tifffile and the design name them.
MODEL_PIXEL_SCALE = 33550
MODEL_TIEPOINT = 33922
MODEL_TRANSFORMATION = 34264
GEOKEY_DIRECTORY = 34735
GDAL_METADATA = 42112
GDAL_NODATA = 42113
COMPRESSION = 259
STRIP_OFFSETS = 273
STRIP_BYTE_COUNTS = 279
PREDICTOR = 317
SAMPLE_FORMAT = 339

# GeoKey ids.
GT_MODEL_TYPE = 1024
GT_RASTER_TYPE = 1025
GEOGRAPHIC_TYPE = 2048
GEOG_LINEAR_UNITS = 2052
PROJECTED_CS_TYPE = 3072
PROJ_LINEAR_UNITS = 3076
VERTICAL_UNITS = 4099

PIXEL_IS_AREA = 1
PIXEL_IS_POINT = 2
USER_DEFINED = 32767
METRE = 9001
FOOT = 9002

#: ETRS89 / UTM 33N, the real fixture's projected CRS.
EPSG_UTM33 = 25833
#: NAD83 / New York Long Island (ftUS): projected, but its axes are US feet.
EPSG_FEET = 2263
#: WGS 84, geographic.
EPSG_WGS84 = 4326
#: WGS 84 geocentric: neither projected nor geographic (pyproj 3.8.0,
#: `type_name` "Geocentric CRS").
EPSG_GEOCENTRIC = 4978
#: Not an EPSG code pyproj can resolve (measured: `CRSError`).
EPSG_UNRESOLVABLE = 9999
#: ETRS89-NOR / UTM 32N + NN2000 height: a compound CRS that pyproj reports as
#: projected, with 3 axes (§5 refusal 13a). Its horizontal part is 11022.
EPSG_COMPOUND = 5972
EPSG_COMPOUND_HORIZONTAL = 11022
#: LUREF / Luxembourg TM (3D): "Projected CRS" with an ellipsoidal-height axis,
#: not compound, so only the axis count catches it (§5 refusal 13a).
EPSG_PROJECTED_3D = 9895
#: WGS 84 + EGM96 height: compound with a geographic horizontal part, so
#: `is_projected` is false and it stays refusal 13 (§5 refusal 13a, "Order").
EPSG_COMPOUND_GEOGRAPHIC = 9707

TIE_X = 500_000.0
TIE_Y = 6_600_000.0
DELTA_X = 10.0
DELTA_Y = 5.0
ROWS = 3
COLS = 4

TIEPOINT: tuple[float, ...] = (0.0, 0.0, 0.0, TIE_X, TIE_Y, 0.0)
SCALE: tuple[float, ...] = (DELTA_X, DELTA_Y, 0.0)
BASE_KEYS: Mapping[int, int] = {
    GT_MODEL_TYPE: 1,  # ModelTypeProjected
    GT_RASTER_TYPE: PIXEL_IS_POINT,
    PROJECTED_CS_TYPE: EPSG_UTM33,
}


def elevations(dtype: Any = np.float32, rows: int = ROWS, cols: int = COLS) -> np.ndarray:
    """A `rows x cols` grid of distinct values: cell (r, c) holds `10 r + c`."""
    r, c = np.indices((rows, cols))
    return np.asarray(10 * r + c).astype(dtype)


#: A GeoKey value as written: inline SHORT, double(s) in 34736, or ASCII in 34737.
GeoKeyValue = int | float | tuple[float, ...] | str


def with_keys(
    changes: Mapping[int, GeoKeyValue | None] | None = None,
    base: Mapping[int, GeoKeyValue] = BASE_KEYS,
) -> dict[int, GeoKeyValue]:
    """`base` (default `BASE_KEYS`) with `changes` applied; `None` deletes that key."""
    merged: dict[int, GeoKeyValue | None] = {**base, **(changes or {})}
    return {k: v for k, v in merged.items() if v is not None}


def _geokey_tags(geokeys: Mapping[int, GeoKeyValue]) -> list[tuple[int, str, int, Any, bool]]:
    """Split by value type into SHORTs, doubles and texts, for `geokey_tags`."""
    shorts = {k: v for k, v in geokeys.items() if isinstance(v, int)}
    doubles = {k: v for k, v in geokeys.items() if isinstance(v, float | tuple)}
    texts = {k: v for k, v in geokeys.items() if isinstance(v, str)}
    return geokey_tags(shorts, doubles, texts)


def _georeference_tags(
    *,
    tiepoint: Sequence[float] | None,
    scale: Sequence[float] | None,
    geokeys: Mapping[int, GeoKeyValue] | None,
    transformation: Sequence[float] | None,
    nodata: str | None,
    gdal_metadata: str | None,
) -> list[tuple[int, str, int, Any, bool]]:
    tags: list[tuple[int, str, int, Any, bool]] = []
    if tiepoint is not None:
        tags.append((MODEL_TIEPOINT, "d", len(tiepoint), tuple(tiepoint), True))
    if scale is not None:
        tags.append((MODEL_PIXEL_SCALE, "d", len(scale), tuple(scale), True))
    if transformation is not None:
        tags.append((MODEL_TRANSFORMATION, "d", len(transformation), tuple(transformation), True))
    if geokeys is not None:
        tags += _geokey_tags(geokeys)
    if nodata is not None:
        tags.append((GDAL_NODATA, "s", 0, nodata, True))
    if gdal_metadata is not None:
        tags.append((GDAL_METADATA, "s", 0, gdal_metadata, True))
    return tags


def micro_tiff(
    array: np.ndarray | None = None,
    *,
    tiepoint: Sequence[float] | None = TIEPOINT,
    scale: Sequence[float] | None = SCALE,
    geokeys: Mapping[int, GeoKeyValue] | None = None,
    transformation: Sequence[float] | None = None,
    nodata: str | None = None,
    gdal_metadata: str | None = None,
    compression: str | None = None,
    extra_pages: Sequence[tuple[np.ndarray, int]] = (),
    **write_kwargs: Any,
) -> io.BytesIO:
    """A GeoTIFF in memory, positioned at 0.

    `geokeys=None` means `BASE_KEYS`; pass `{}` for an empty directory. A
    `tiepoint` or `scale` of `None` omits that tag. `extra_pages` are
    `(array, subfiletype)` pairs written after page 0; subfiletype 1 is
    FILETYPE_REDUCEDIMAGE.
    """
    data = elevations() if array is None else array
    tags = _georeference_tags(
        tiepoint=tiepoint,
        scale=scale,
        geokeys=dict(BASE_KEYS) if geokeys is None else geokeys,
        transformation=transformation,
        nodata=nodata,
        gdal_metadata=gdal_metadata,
    )
    stream = io.BytesIO()
    with tifffile.TiffWriter(stream) as writer:
        writer.write(data, extratags=tags, compression=compression, **write_kwargs)
        for page, subfiletype in extra_pages:
            writer.write(page, subfiletype=subfiletype)
    stream.seek(0)
    return stream


def with_short_tag(stream: io.BytesIO, tag: int, old: int, new: int) -> io.BytesIO:
    """Rewrite the IFD entry `tag`, a little-endian SHORT of count 1, from `old` to `new`.

    The assertion guards against the entry being found zero times or twice, so
    a patch can never land on some other bytes that happen to match.
    """
    buffer = bytearray(stream.getvalue())
    entry = tag.to_bytes(2, "little") + b"\x03\x00\x01\x00\x00\x00" + old.to_bytes(2, "little")
    assert buffer.count(entry) == 1, f"expected exactly one tag {tag} entry with value {old}"
    at = buffer.index(entry) + 8
    buffer[at : at + 2] = new.to_bytes(2, "little")
    return io.BytesIO(bytes(buffer))


def with_compression_tag(stream: io.BytesIO, code: int) -> io.BytesIO:
    """Rewrite an uncompressed file's `Compression` (259) value to `code`.

    tifffile cannot *write* LZW without `imagecodecs`, but the refusal under
    test must fire before any decode, so the strip's bytes never need to be
    valid LZW.
    """
    return with_short_tag(stream, COMPRESSION, 1, code)


LZW = 5
PACKBITS = 32773
#: In tifffile's `COMPRESSION` enum, but not decodable even with imagecodecs
#: (§5 refusal 11, round 3, measured with imagecodecs 2026.8.16).
THUNDERSCAN = 32809
#: Not a member of tifffile's `COMPRESSION` enum at all.
UNKNOWN_COMPRESSION = 60000
FLOATING_POINT_PREDICTOR = 3


def floating_point_predictor_tiff() -> io.BytesIO:
    """A Deflate float32 tile declaring `Predictor` (317) = 3, as GDAL writes float DEMs.

    tifffile refuses to write predictor 3 without `imagecodecs`, and refuses
    predictor 2 for floats. So an int32 tile is written with predictor 2 and
    then relabelled: `SampleFormat` 2 -> 3 (float) and `Predictor` 2 -> 3. The
    strip is not valid floating-point-predicted data. It never has to be: the
    refusal is a capability probe that runs before any decode (§8).
    """
    stream = micro_tiff(elevations(np.int32), compression="deflate", predictor=2)
    stream = with_short_tag(stream, SAMPLE_FORMAT, 2, 3)
    return with_short_tag(stream, PREDICTOR, 2, FLOATING_POINT_PREDICTOR)


def packbits_tiff() -> io.BytesIO:
    """The baseline tile with its one strip re-encoded as real PackBits.

    tifffile cannot write PackBits without `imagecodecs` but, measured, reads
    it with its built-in decoder (§8, problem 1). The strip is appended to the
    file as one PackBits literal run (header byte `n - 1`, then `n` bytes), and
    `StripOffsets`, `StripByteCounts` and `Compression` are repointed at it.
    """
    stream = micro_tiff()
    page = tifffile.TiffFile(stream).pages.first
    raw = elevations().tobytes()
    assert len(raw) <= 128, "one literal run holds at most 128 bytes"
    strip = bytes([len(raw) - 1]) + raw
    buffer = bytearray(stream.getvalue())
    offset = len(buffer)
    buffer += strip
    for tag, value in ((STRIP_OFFSETS, offset), (STRIP_BYTE_COUNTS, len(strip))):
        entry = page.tags[tag]
        assert entry.count == 1 and entry.dtype == 4, "expected one inline LONG"
        buffer[entry.valueoffset : entry.valueoffset + 4] = value.to_bytes(4, "little")
    return with_compression_tag(io.BytesIO(bytes(buffer)), PACKBITS)


# --------------------------------------------------------------------------
# Streams tifffile itself cannot read (§14, choice C). Each names the stage of
# `decode_dem` whose call into tifffile raises: "TIFF structure" for
# `TiffFile(source)` and the page walk, "pixel data" for `asarray`.
# --------------------------------------------------------------------------


def _strip(stream: io.BytesIO) -> tuple[bytearray, int, int]:
    """The file's bytes, and the offset and byte count of its one strip."""
    page = tifffile.TiffFile(stream).pages.first
    (offset,), (count,) = page.dataoffsets, page.databytecounts
    return bytearray(stream.getvalue()), offset, count


def truncated_header() -> io.BytesIO:
    """Cut inside the 8-byte header, after the byte order and magic number."""
    return io.BytesIO(micro_tiff().getvalue()[:6])


def truncated_ifd() -> io.BytesIO:
    """Cut ten bytes into the first IFD, which the header points at."""
    data = micro_tiff().getvalue()
    first_ifd = int.from_bytes(data[4:8], "little")
    return io.BytesIO(data[: first_ifd + 10])


def truncated_second_ifd() -> io.BytesIO:
    """A two-page file cut six bytes into the second IFD.

    Page 0 and its strip are whole, so `TiffFile(source)` opens it; only the
    page walk reaches IFD 1 and fails there. This is the case that needs the
    walk inside the "TIFF structure" stage, not just the open.
    """
    data = micro_tiff(extra_pages=[(elevations(), 1)]).getvalue()
    second_ifd = tifffile.TiffFile(io.BytesIO(data)).pages[1].offset
    return io.BytesIO(data[: second_ifd + 6])


def truncated_strip() -> io.BytesIO:
    """Cut halfway through the strip. The strip is the file's last bytes, so
    the IFD and every tag value survive and only the pixel read can fail."""
    data, offset, count = _strip(micro_tiff())
    assert offset + count == len(data), "the strip must be the last thing in the file"
    return io.BytesIO(bytes(data[: offset + count // 2]))


def _corrupt_strip(stream: io.BytesIO) -> io.BytesIO:
    data, offset, count = _strip(stream)
    data[offset : offset + count] = b"\xff" * count
    return io.BytesIO(bytes(data))


def corrupt_deflate() -> io.BytesIO:
    """A Deflate tile whose strip is overwritten with 0xFF: not a zlib stream."""
    return _corrupt_strip(micro_tiff(compression="deflate"))


def corrupt_lzw() -> io.BytesIO:
    """An LZW tile whose strip is overwritten with 0xFF. Writing LZW needs
    imagecodecs, so only a `needs_codecs` test may build this."""
    return _corrupt_strip(micro_tiff(compression="lzw"))


UNDECODABLE: Mapping[str, tuple[Callable[[], io.BytesIO], str]] = {
    "empty": (lambda: io.BytesIO(b""), "TIFF structure"),
    "not_tiff": (lambda: io.BytesIO(b"this is not a TIFF file\n" * 4), "TIFF structure"),
    "truncated_header": (truncated_header, "TIFF structure"),
    "truncated_ifd": (truncated_ifd, "TIFF structure"),
    "truncated_second_ifd": (truncated_second_ifd, "TIFF structure"),
    "truncated_strip": (truncated_strip, "pixel data"),
    "corrupt_deflate": (corrupt_deflate, "pixel data"),
}
UNDECODABLE_WITH_CODECS: Mapping[str, tuple[Callable[[], io.BytesIO], str]] = {
    "corrupt_lzw": (corrupt_lzw, "pixel data"),
}


# --------------------------------------------------------------------------
# The catalogue: one representative fixture per named refusal (design §5).
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Refusal:
    """One named refusal: how to build it, and what its message must name.

    `witness` inspects the fixture with tifffile alone and returns True when
    the defect is present; it is what lets the suite say a red test is red
    for the right reason.
    """

    name: str
    build: Callable[[], io.BytesIO]
    witness: Callable[[tifffile.TiffFile], bool]
    must_name: tuple[str, ...]
    decode_kwargs: Mapping[str, Any] = field(default_factory=dict)


def _geo(tif: tifffile.TiffFile) -> dict[Any, Any]:
    return tif.geotiff_metadata or {}


def _tag(tif: tifffile.TiffFile, code: int) -> Any:
    tag = tif.pages.first.tags.get(code)
    return None if tag is None else tag.value


_IDENTITY_TRANSFORM = (
    DELTA_X, 0.0, 0.0, TIE_X,
    0.0, -DELTA_Y, 0.0, TIE_Y,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0,
)  # fmt: skip

REFUSALS: tuple[Refusal, ...] = (
    Refusal(
        "refuses_missing_georeferencing",
        lambda: micro_tiff(tiepoint=None),
        lambda t: _tag(t, MODEL_TIEPOINT) is None and _tag(t, MODEL_PIXEL_SCALE) is not None,
        ("33922", "ModelTiepointTag"),
    ),
    Refusal(
        "refuses_model_transformation",
        lambda: micro_tiff(transformation=_IDENTITY_TRANSFORM),
        lambda t: _tag(t, MODEL_TRANSFORMATION) is not None,
        ("34264", "ModelTransformationTag"),
    ),
    Refusal(
        "refuses_multiple_tiepoints",
        lambda: micro_tiff(tiepoint=(*TIEPOINT, 3.0, 2.0, 0.0, TIE_X + 30, TIE_Y - 10, 0.0)),
        lambda t: len(_tag(t, MODEL_TIEPOINT)) == 12,
        ("33922", "ModelTiepointTag"),
    ),
    Refusal(
        "refuses_nonzero_pixel_scale_z",
        lambda: micro_tiff(scale=(DELTA_X, DELTA_Y, 2.5)),
        lambda t: _tag(t, MODEL_PIXEL_SCALE)[2] == 2.5,
        ("33550", "ModelPixelScaleTag", "2.5"),
    ),
    Refusal(
        "refuses_nonpositive_pixel_scale",
        lambda: micro_tiff(scale=(DELTA_X, -DELTA_Y, 0.0)),
        lambda t: _tag(t, MODEL_PIXEL_SCALE)[1] < 0,
        ("33550", "ModelPixelScaleTag", "-5"),
    ),
    Refusal(
        "refuses_unknown_raster_type",
        lambda: micro_tiff(geokeys=with_keys({GT_RASTER_TYPE: None})),
        lambda t: "GTRasterTypeGeoKey" not in _geo(t),
        ("1025", "GTRasterTypeGeoKey"),
    ),
    Refusal(
        "refuses_degenerate_shape",
        lambda: micro_tiff(elevations(rows=1)),
        lambda t: t.pages.first.shape == (1, COLS),
        ("257", "ImageLength", "at least 2"),
    ),
    Refusal(
        "refuses_ambiguous_pages",
        lambda: micro_tiff(extra_pages=[(elevations(), 0)]),
        lambda t: len(t.pages) == 2 and not t.pages[1].is_reduced,
        ("254", "NewSubfileType"),
    ),
    Refusal(
        "refuses_multi_sample",
        lambda: micro_tiff(
            elevations()[..., None].repeat(2, axis=2),
            photometric="minisblack",
            planarconfig="contig",
        ),
        lambda t: t.pages.first.samplesperpixel == 2,
        ("SamplesPerPixel", "2"),
    ),
    Refusal(
        "refuses_unsupported_dtype",
        lambda: micro_tiff(elevations(np.int64)),
        lambda t: t.pages.first.dtype == np.int64,
        ("339", "SampleFormat", "258", "BitsPerSample", "int64"),
    ),
    Refusal(
        "refuses_missing_codec",
        lambda: with_compression_tag(micro_tiff(), LZW),
        lambda t: t.pages.first.compression == LZW,
        ("259", "Compression", "LZW"),
    ),
    Refusal(
        "refuses_missing_crs",
        lambda: micro_tiff(geokeys=with_keys({PROJECTED_CS_TYPE: None})),
        lambda t: "ProjectedCSTypeGeoKey" not in _geo(t),
        ("3072", "ProjectedCSTypeGeoKey"),
    ),
    Refusal(
        "refuses_geographic_crs",
        lambda: micro_tiff(geokeys=with_keys({PROJECTED_CS_TYPE: EPSG_WGS84})),
        lambda t: _geo(t).get("ProjectedCSTypeGeoKey") == EPSG_WGS84,
        ("4326",),
    ),
    Refusal(
        "refuses_compound_crs",
        lambda: micro_tiff(geokeys=with_keys({PROJECTED_CS_TYPE: EPSG_COMPOUND})),
        lambda t: _geo(t).get("ProjectedCSTypeGeoKey") == EPSG_COMPOUND,
        ("3072", "ProjectedCSTypeGeoKey", str(EPSG_COMPOUND), "Compound CRS"),
    ),
    Refusal(
        "refuses_non_metre_linear_unit",
        lambda: micro_tiff(geokeys=with_keys({PROJ_LINEAR_UNITS: FOOT})),
        lambda t: _geo(t).get("ProjLinearUnitsGeoKey") == FOOT,
        ("3076", "ProjLinearUnitsGeoKey", "9002"),
    ),
    Refusal(
        "refuses_non_metre_vertical_unit",
        lambda: micro_tiff(geokeys=with_keys({VERTICAL_UNITS: FOOT})),
        lambda t: _geo(t).get("VerticalUnitsGeoKey") == FOOT,
        ("4099", "VerticalUnitsGeoKey", "9002"),
    ),
    Refusal(
        "refuses_unparseable_nodata",
        lambda: micro_tiff(nodata="void"),
        lambda t: _tag(t, GDAL_NODATA) == "void",
        ("42113", "GDAL_NODATA", "void"),
    ),
    Refusal(
        "refuses_nodata_not_representable",
        lambda: micro_tiff(nodata="0.1"),
        lambda t: _tag(t, GDAL_NODATA) == "0.1" and t.pages.first.dtype == np.float32,
        ("42113", "GDAL_NODATA", "0.1", "float32"),
    ),
    Refusal(
        "refuses_contradictory_nodata_override",
        lambda: micro_tiff(nodata="-32767"),
        lambda t: _tag(t, GDAL_NODATA) == "-32767",
        ("42113", "GDAL_NODATA", "-32767", "-9999"),
        decode_kwargs={"nodata": -9999.0},
    ),
)


# --------------------------------------------------------------------------
# A projected CRS given by parameters (3072 = 32767), read as the EPSG code it
# matches: `docs/increments/geotiff-crs-by-parameters.md`, section 8.
# --------------------------------------------------------------------------

# GeoKey ids beyond the ones above (GeoTIFF 1.0 names, as tifffile spells them).
GT_CITATION = 1026
GEOG_CITATION = 2049
GEOG_GEODETIC_DATUM = 2050
GEOG_PRIME_MERIDIAN = 2051
GEOG_ANGULAR_UNITS = 2054
GEOG_ELLIPSOID = 2056
GEOG_SEMI_MAJOR_AXIS = 2057
GEOG_SEMI_MINOR_AXIS = 2058
GEOG_INV_FLATTENING = 2059
GEOG_PRIME_MERIDIAN_LONG = 2061
GEOG_TOWGS84 = 2062
PROJECTION = 3074
PROJ_COORD_TRANS = 3075
PROJ_STD_PARALLEL_1 = 3078
PROJ_STD_PARALLEL_2 = 3079
PROJ_NAT_ORIGIN_LONG = 3080
PROJ_NAT_ORIGIN_LAT = 3081
PROJ_FALSE_EASTING = 3082
PROJ_FALSE_NORTHING = 3083
PROJ_FALSE_ORIGIN_LONG = 3084
PROJ_FALSE_ORIGIN_LAT = 3085
PROJ_FALSE_ORIGIN_EASTING = 3086
PROJ_FALSE_ORIGIN_NORTHING = 3087
PROJ_SCALE_AT_NAT_ORIGIN = 3092

#: ProjCoordTransGeoKey (3075) codes.
CT_TRANSVERSE_MERCATOR = 1
CT_LAMBERT_CONF_CONIC_2SP = 8
CT_ALBERS_EQUAL_AREA = 11

DEGREE = 9102
GRAD = 9105
#: MGI, ETRS89, Xian 1980: geographic 2D CRSs, by EPSG code.
EPSG_MGI = 4312
EPSG_ETRS89 = 4258
EPSG_XIAN_1980 = 4610
#: MGI's geodetic datum (a datum code, not a CRS code), for 2050.
EPSG_MGI_DATUM = 6312
#: MGI / Austria Lambert: what the Austrian file is, by its parameters.
EPSG_AUSTRIA_LAMBERT = 31287
#: The file's longitude of origin; EPSG:31287 has 13.3333333333333.
AUSTRIA_LON = 13.33333333300013
#: Bessel 1841, as the Austrian file's 2057 and 2059 carry it.
BESSEL_A = 6377397.155
BESSEL_RF = 299.1528128000033

#: The eighteen GeoKeys of `../rasputin_data/austria_dgm10/dhm_at_lamb_10m_2018.tif`,
#: section 1's table, written with the types that file uses (citations ASCII,
#: parameters and ellipsoid doubles, the rest SHORT).
AUSTRIA_KEYS: Mapping[int, GeoKeyValue] = {
    GT_MODEL_TYPE: 1,
    GT_RASTER_TYPE: PIXEL_IS_AREA,
    GT_CITATION: "MGI_Austria_Lambert",
    GEOGRAPHIC_TYPE: EPSG_MGI,
    GEOG_CITATION: "MGI",
    GEOG_ANGULAR_UNITS: DEGREE,
    GEOG_SEMI_MAJOR_AXIS: BESSEL_A,
    GEOG_INV_FLATTENING: BESSEL_RF,
    PROJECTED_CS_TYPE: USER_DEFINED,
    PROJECTION: USER_DEFINED,
    PROJ_COORD_TRANS: CT_LAMBERT_CONF_CONIC_2SP,
    PROJ_LINEAR_UNITS: METRE,
    PROJ_STD_PARALLEL_1: 46.0,
    PROJ_STD_PARALLEL_2: 49.0,
    PROJ_FALSE_ORIGIN_LONG: AUSTRIA_LON,
    PROJ_FALSE_ORIGIN_LAT: 47.5,
    PROJ_FALSE_ORIGIN_EASTING: 400000.0,
    PROJ_FALSE_ORIGIN_NORTHING: 400000.0,
}

#: The ellipsoid keys, which only fit the MGI base; a variant on another datum
#: drops them (their check, section 4 row 6, has its own red test).
_BESSEL_KEYS: Mapping[int, None] = {GEOG_SEMI_MAJOR_AXIS: None, GEOG_INV_FLATTENING: None}
_LAMBERT_KEYS = (
    PROJ_STD_PARALLEL_1, PROJ_STD_PARALLEL_2, PROJ_FALSE_ORIGIN_LONG,
    PROJ_FALSE_ORIGIN_LAT, PROJ_FALSE_ORIGIN_EASTING, PROJ_FALSE_ORIGIN_NORTHING,
)  # fmt: skip


def austria_keys(changes: Mapping[int, GeoKeyValue | None] | None = None) -> dict[int, GeoKeyValue]:
    """`AUSTRIA_KEYS` with `changes` applied; `None` deletes that key."""
    return with_keys(changes, base=AUSTRIA_KEYS)


def lambert_on(geographic: int, **changes: float) -> dict[int, GeoKeyValue]:
    """The Austrian Lambert parameters on another base geographic CRS, without
    the Bessel ellipsoid keys. `changes` replaces parameters by key name, e.g.
    `PROJ_FALSE_ORIGIN_EASTING=400000.001`."""
    named: dict[int, GeoKeyValue | None] = {globals()[k]: v for k, v in changes.items()}
    return austria_keys({**_BESSEL_KEYS, GEOGRAPHIC_TYPE: geographic, **named})


def transverse_mercator_on(
    geographic: int,
    longitude: float,
    *,
    scale: float = 0.9996,
    false_easting: float = 500000.0,
    false_northing: float = 0.0,
    latitude: float = 0.0,
) -> dict[int, GeoKeyValue]:
    """The Austrian key set with its method swapped for transverse Mercator (3075 = 1)
    on `geographic`: the natural-origin keys in, the Lambert keys and the
    Bessel ellipsoid keys out. The citations stay; no citation is read."""
    changes: dict[int, GeoKeyValue | None] = {k: None for k in (*_LAMBERT_KEYS, *_BESSEL_KEYS)}
    changes |= {
        GEOGRAPHIC_TYPE: geographic,
        PROJ_COORD_TRANS: CT_TRANSVERSE_MERCATOR,
        PROJ_NAT_ORIGIN_LAT: latitude,
        PROJ_NAT_ORIGIN_LONG: longitude,
        PROJ_SCALE_AT_NAT_ORIGIN: scale,
        PROJ_FALSE_EASTING: false_easting,
        PROJ_FALSE_NORTHING: false_northing,
    }
    return austria_keys(changes)


def austria_tiff(geokeys: Mapping[int, GeoKeyValue] | None = None) -> io.BytesIO:
    """The usual 3 x 4 baseline carrying `geokeys` (default `AUSTRIA_KEYS`)."""
    return micro_tiff(geokeys=AUSTRIA_KEYS if geokeys is None else geokeys)


#: GeogTOWGS84GeoKey (2062) as a seven-parameter (Helmert) and a
#: three-parameter (translation) shift. The values are the MGI-to-WGS 84 shift
#: of EPSG:1618's family; any values would do, since the key is refused present.
TOWGS84_SEVEN = (577.326, 90.129, 463.919, 5.137, 1.474, 5.297, 2.4232)
TOWGS84_THREE = (577.326, 90.129, 463.919)


@dataclass(frozen=True)
class ParametricDefect:
    """One check of section 4 (rows 1 to 11) on its one defect.

    `changes` applied to `AUSTRIA_KEYS` is the fixture; `witness` reads the
    tifffile `geotiff_metadata` dict and holds when the defect is there;
    `must_name` is what the message must hold (case-insensitive), beyond `P`.
    """

    name: str
    row: int
    changes: Mapping[int, GeoKeyValue | None]
    witness: Callable[[dict[Any, Any]], bool]
    must_name: tuple[str, ...]

    def build(self) -> io.BytesIO:
        return austria_tiff(austria_keys(self.changes))


def _has(name: str) -> Callable[[dict[Any, Any]], bool]:
    return lambda g: name in g


def _lacks(name: str) -> Callable[[dict[Any, Any]], bool]:
    return lambda g: name not in g


def _is(name: str, value: Any) -> Callable[[dict[Any, Any]], bool]:
    return lambda g: name in g and g[name] == value


PARAMETRIC_DEFECTS: tuple[ParametricDefect, ...] = (
    ParametricDefect(
        "projection_by_code", 1, {PROJECTION: 16033},
        _is("ProjectionGeoKey", 16033), ("3074", "ProjectionGeoKey", "16033"),
    ),
    ParametricDefect(
        "method_absent", 2, {PROJ_COORD_TRANS: None},
        _lacks("ProjCoordTransGeoKey"), ("3075", "ProjCoordTransGeoKey"),
    ),
    ParametricDefect(
        "method_albers", 3, {PROJ_COORD_TRANS: CT_ALBERS_EQUAL_AREA},
        _is("ProjCoordTransGeoKey", CT_ALBERS_EQUAL_AREA), ("3075", "ProjCoordTransGeoKey", "11"),
    ),
    ParametricDefect(
        "datum_absent", 4, {GEOGRAPHIC_TYPE: None},
        _lacks("GeographicTypeGeoKey"), ("2048", "GeographicTypeGeoKey", "absent"),
    ),
    ParametricDefect(
        "datum_user_defined", 4, {GEOGRAPHIC_TYPE: USER_DEFINED},
        _is("GeographicTypeGeoKey", USER_DEFINED), ("2048", "GeographicTypeGeoKey"),
    ),
    ParametricDefect(
        "datum_unresolvable", 4, {GEOGRAPHIC_TYPE: EPSG_UNRESOLVABLE},
        _is("GeographicTypeGeoKey", EPSG_UNRESOLVABLE), ("2048", "GeographicTypeGeoKey", "9999"),
    ),
    ParametricDefect(
        "datum_projected", 4, {GEOGRAPHIC_TYPE: EPSG_UTM33},
        _is("GeographicTypeGeoKey", EPSG_UTM33),
        ("2048", "GeographicTypeGeoKey", "25833", "Projected CRS"),
    ),
    ParametricDefect(
        "datum_key_present", 5, {GEOG_GEODETIC_DATUM: EPSG_MGI_DATUM},
        _is("GeogGeodeticDatumGeoKey", EPSG_MGI_DATUM),
        ("2050", "GeogGeodeticDatumGeoKey", "2048"),
    ),
    ParametricDefect(
        "semi_major_disagrees", 6, {GEOG_SEMI_MAJOR_AXIS: 6378137.0},
        _is("GeogSemiMajorAxisGeoKey", 6378137.0),
        ("2057", "GeogSemiMajorAxisGeoKey", "6378137", "EPSG:4312"),
    ),
    ParametricDefect(
        "inverse_flattening_disagrees", 6, {GEOG_INV_FLATTENING: 299.0},
        _is("GeogInvFlatteningGeoKey", 299.0),
        ("2059", "GeogInvFlatteningGeoKey", "EPSG:4312"),
    ),
    ParametricDefect(
        "towgs84_seven", 7, {GEOG_TOWGS84: TOWGS84_SEVEN},
        _is("GeogTOWGS84GeoKey", TOWGS84_SEVEN),
        ("2062", "GeogTOWGS84GeoKey", "577.326", "2.4232"),
    ),
    ParametricDefect(
        "towgs84_three", 7, {GEOG_TOWGS84: TOWGS84_THREE},
        _is("GeogTOWGS84GeoKey", TOWGS84_THREE),
        ("2062", "GeogTOWGS84GeoKey", "577.326", "463.919"),
    ),
    ParametricDefect(
        "angular_unit_grad", 8, {GEOG_ANGULAR_UNITS: GRAD},
        _is("GeogAngularUnitsGeoKey", GRAD), ("2054", "GeogAngularUnitsGeoKey", "9105"),
    ),
    ParametricDefect(
        "linear_unit_absent", 9, {PROJ_LINEAR_UNITS: None},
        _lacks("ProjLinearUnitsGeoKey"), ("3076", "ProjLinearUnitsGeoKey", "absent"),
    ),
    ParametricDefect(
        "linear_unit_foot", 10, {PROJ_LINEAR_UNITS: FOOT},
        _is("ProjLinearUnitsGeoKey", FOOT), ("3076", "ProjLinearUnitsGeoKey", "9002"),
    ),
    ParametricDefect(
        "parameter_absent", 11, {PROJ_FALSE_ORIGIN_LONG: None},
        _lacks("ProjFalseOriginLongGeoKey"), ("3084", "ProjFalseOriginLongGeoKey", "absent"),
    ),
)  # fmt: skip
