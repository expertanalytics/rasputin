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

GeoKeys are written as a raw `GeoKeyDirectoryTag` (34735) of inline SHORTs.
That is the on-disk form tifffile's `geotiff_metadata` decodes, so the reader
under test sees what a real file would give it.
"""

from __future__ import annotations

import io
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import tifffile

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
#: Not an EPSG code pyproj can resolve (measured: `CRSError`).
EPSG_UNRESOLVABLE = 9999

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


def with_keys(changes: Mapping[int, int | None] | None = None) -> dict[int, int]:
    """`BASE_KEYS` with `changes` applied; a `None` value deletes that key."""
    merged: dict[int, int | None] = {**BASE_KEYS, **(changes or {})}
    return {k: v for k, v in merged.items() if v is not None}


def _geokey_directory(geokeys: Mapping[int, int]) -> tuple[int, str, int, list[int], bool]:
    shorts = [1, 1, 0, len(geokeys)]
    for key_id, value in sorted(geokeys.items()):
        shorts += [key_id, 0, 1, value]  # location 0: value is inline
    return (GEOKEY_DIRECTORY, "H", len(shorts), shorts, True)


def _georeference_tags(
    *,
    tiepoint: Sequence[float] | None,
    scale: Sequence[float] | None,
    geokeys: Mapping[int, int] | None,
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
        tags.append(_geokey_directory(geokeys))
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
    geokeys: Mapping[int, int] | None = None,
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
