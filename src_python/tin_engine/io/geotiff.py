"""A GeoTIFF stream to a `DemTile`, or a `GeoTiffError` that says why not.

Increment 11, `docs/increments/11-raster-ingestion.md`. The reader's job is
refusal (ruling 3): every refusal in §5 raises `GeoTiffError` naming the tag or
GeoKey by number and name and the file's value. None is a bare `assert`,
because `python -O` strips those (prior art §5.9).

The module takes a binary stream and never a path (§2), and imports neither
`_core` nor anything first-party except `io.models` (ruling 1). GeoKeys come
from tifffile's own `geotiff_metadata` (ruling 5); the CRS is resolved only
through `pyproj.CRS.from_epsg` and tested on the constructed CRS (rulings 6, 7).
No transform happens here, so there is no `Transformer` (ruling 8).
"""

from __future__ import annotations

import importlib.util
import math
from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any, BinaryIO, Literal, NoReturn

import numpy as np
import pyproj
import tifffile
from pyproj.exceptions import CRSError

from .models import DemTile, GeoTiffError, RasterMeta

USER_DEFINED = 32767
METRE = 9001
PIXEL_IS_AREA = 1
PIXEL_IS_POINT = 2

#: §7: file dtype to array dtype. Every row is exact; anything else is refused.
PROMOTION: dict[np.dtype[Any], np.dtype[Any]] = {
    np.dtype(src): np.dtype(dst)
    for dst, sources in {
        "float32": ("int8", "uint8", "int16", "uint16", "float32"),
        "float64": ("int32", "uint32", "float64"),
    }.items()
    for src in sources
}

NodataSource = Literal["tag", "caller", "absent"]


def decode_dem(source: BinaryIO, *, nodata: float | None = None) -> DemTile:
    """Decode page 0 of the GeoTIFF in `source` into a `DemTile`.

    `nodata` is the caller asserting a sentinel (§6). It must be representable
    in the file's dtype, and must agree with tag 42113 when that is present. A
    `bool` is refused with `TypeError`: it is a caller's mistake, not a file's.

    Raises `GeoTiffError` for every §5 refusal, and (§14, choice C) when
    tifffile or a codec fails on the file; that message names the stage, the
    original type and text, and chains the original. `MemoryError` and the
    reader's own bugs are not wrapped.

    The whole page is decoded, so memory is O(rows x cols), briefly doubled
    while `DemTile` takes its read-only copy (§7). `source` is read but not
    closed; closing it is the caller's.
    """
    if isinstance(nodata, bool | np.bool_):
        raise TypeError(f"nodata= must be a number or None, not {type(nodata).__name__}")
    with _stage("TIFF structure"):
        tif = tifffile.TiffFile(source)
    with tif:
        with _stage("TIFF structure"):
            page = _single_page(tif)
        dtype = _check_page(page)
        tie, scale = _georeferencing(page)
        with _stage("GeoKey directory"):
            geokeys: dict[str, Any] = tif.geotiff_metadata or {}
        x_min, y_max, delta_x, delta_y, area = _placement(tie, scale, geokeys)
        epsg = _projected_epsg(geokeys)
        vertical = geokeys.get("VerticalUnitsGeoKey")
        if vertical is not None and int(vertical) != METRE:
            raise GeoTiffError(
                f"VerticalUnitsGeoKey (4099) = {int(vertical)}; only metres ({METRE}) are read"
            )
        sentinel, source_of = _nodata(page, dtype, nodata)
        with _stage("pixel data"):
            array = page.asarray().astype(PROMOTION[dtype], copy=False)
    rows, cols = array.shape
    meta = RasterMeta(
        x_min=x_min,
        y_max=y_max,
        delta_x=delta_x,
        delta_y=delta_y,
        cols=cols,
        rows=rows,
        epsg=epsg,
        nodata=sentinel,
        nodata_source=source_of,
        pixel_is_area=area,
        vertical_unit_assumed=vertical is None,
    )
    return DemTile(meta=meta, array=array)


@contextmanager
def _stage(stage: str) -> Iterator[None]:
    """§14, choice C: a failure inside tifffile or a codec becomes `GeoTiffError`.

    Wraps call sites, not types, so the reader's own code stays outside and its
    bugs surface as bugs. A refusal raised inside passes through unchanged.
    """
    try:
        yield
    except (MemoryError, GeoTiffError):
        raise
    except Exception as error:
        raise GeoTiffError(
            f"{stage}: tifffile could not decode the file: {type(error).__name__}: {error}"
        ) from error


def _single_page(tif: tifffile.TiffFile) -> tifffile.TiffPage:
    """Page 0, refusing a second full-resolution page (§5 refusal 8, 0-based)."""
    pages = tif.pages
    for index in range(1, len(pages)):
        extra = pages[index]
        where = f"on page {index} of {len(pages)}: only page 0 may be full resolution"
        if not isinstance(extra, tifffile.TiffPage):
            # A TiffFrame has no subfiletype to consult (round 3).
            raise GeoTiffError(f"NewSubfileType (254) could not be read (a frame) {where}")
        if not extra.is_reduced:
            raise GeoTiffError(f"NewSubfileType (254) = {int(extra.subfiletype)} {where}")
    return pages.first


def _check_page(page: tifffile.TiffPage) -> np.dtype[Any]:
    """Samples, shape, dtype and codec (§5 refusals 7, 9, 10, 11). Returns the file dtype."""
    if page.samplesperpixel != 1:
        raise GeoTiffError(f"SamplesPerPixel (277) = {page.samplesperpixel}; a DEM has one band")
    short = [
        f"{name} ({code}) = {value}"
        for name, code, value in (
            ("ImageLength", 257, page.imagelength),
            ("ImageWidth", 256, page.imagewidth),
        )
        if value < 2
    ]
    if short:
        raise GeoTiffError(f"{' and '.join(short)}: need at least 2 rows and 2 columns")
    dtype = np.dtype(page.dtype or "V")  # tifffile gives None for a layout numpy lacks
    if dtype not in PROMOTION:
        raise GeoTiffError(
            f"SampleFormat (339) = {int(page.sampleformat)}, BitsPerSample (258) = "
            f"{page.bitspersample} (numpy {dtype}): not in the promotion table"
        )
    # §8, problem 1: ask tifffile, never a list of scheme names. Membership
    # resolves the codec and answers False instead of raising.
    tiff = tifffile.TIFF
    for name, code, value, known, enum in (
        ("Compression", 259, int(page.compression), tiff.DECOMPRESSORS, tifffile.COMPRESSION),
        ("Predictor", 317, int(page.predictor), tiff.PREDICTORS, tifffile.PREDICTOR),
    ):
        if value not in known:
            scheme = f" ({enum(value).name})" if value in enum else ""
            advice = (
                "cannot be decoded by tifffile, even with imagecodecs installed"
                if importlib.util.find_spec("imagecodecs") is not None
                else "cannot be decoded as installed; the `codecs` extra (imagecodecs) "
                "may decode it: pip install 'rasputin[codecs]'"
            )
            raise GeoTiffError(f"{name} ({code}) = {value}{scheme} {advice}")
    return dtype


def _georeferencing(page: tifffile.TiffPage) -> tuple[list[float], list[float]]:
    """Tie point and scale from `page.tags` alone (§5 refusals 1-5, round 3).

    Runs before `geotiff_metadata`, whose reshape of a wrong-length 33922 or
    34264 would otherwise raise first.
    """
    tags = page.tags
    if 34264 in tags:
        raise GeoTiffError("ModelTransformationTag (34264) is present; only north-up is read")
    for code, name in ((33922, "ModelTiepointTag"), (33550, "ModelPixelScaleTag")):
        if code not in tags:
            raise GeoTiffError(f"{name} ({code}) is absent: the file is not georeferenced")
    # A one-value DOUBLE tag reads back as a bare float.
    tie = [float(v) for v in np.atleast_1d(tags[33922].value)]
    scale = [float(v) for v in np.atleast_1d(tags[33550].value)]
    if len(tie) > 6:
        raise GeoTiffError(
            f"ModelTiepointTag (33922) has {len(tie)} values; exactly 6 (one tie point) are read"
        )
    if len(tie) < 6:
        raise GeoTiffError(f"ModelTiepointTag (33922) has {len(tie)} values; need 6")
    if len(scale) != 3:
        raise GeoTiffError(f"ModelPixelScaleTag (33550) has {len(scale)} values; need 3")
    if not all(math.isfinite(v) for v in (*tie[:2], *tie[3:5])):
        raise GeoTiffError(f"ModelTiepointTag (33922) = {tie}: need finite I, J, X, Y")
    if not all(math.isfinite(v) and v > 0 for v in scale[:2]):
        raise GeoTiffError(f"ModelPixelScaleTag (33550) = {scale}: need finite positive X, Y")
    if scale[2] != 0.0:
        raise GeoTiffError(f"ModelPixelScaleTag (33550) ScaleZ = {scale[2]}; must be 0")
    return tie, scale


def _placement(
    tie: list[float], scale: list[float], geokeys: dict[str, Any]
) -> tuple[float, float, float, float, bool]:
    """The node grid's origin and spacing (§4, §5 refusal 6)."""
    raster_type = geokeys.get("GTRasterTypeGeoKey")
    if raster_type is None or int(raster_type) not in (PIXEL_IS_AREA, PIXEL_IS_POINT):
        shown = "absent" if raster_type is None else f"= {int(raster_type)}"
        raise GeoTiffError(f"GTRasterTypeGeoKey (1025) {shown}: registration unknown")
    area = int(raster_type) == PIXEL_IS_AREA
    (i, j, _, x, y, _), (dx, dy, _) = tie, scale
    half = 0.5 if area else 0.0
    return x - i * dx + half * dx, y + j * dy - half * dy, dx, dy, area


def _projected_epsg(geokeys: dict[str, Any]) -> int:
    """Rulings 6 and 7, §5 refusals 12-14: a projected, metre CRS from 3072 alone."""
    projected = geokeys.get("ProjectedCSTypeGeoKey")
    if projected is None:
        _refuse_through_2048(geokeys.get("GeographicTypeGeoKey"))
    code = int(projected)
    crs = _resolve(code)
    if crs is None:
        raise GeoTiffError(f"ProjectedCSTypeGeoKey (3072) = {code} is not a resolvable EPSG code")
    if not crs.is_projected:
        raise GeoTiffError(
            f"ProjectedCSTypeGeoKey (3072) = {code} is a {crs.type_name}, not a projected CRS"
        )
    if crs.is_compound or len(crs.axis_info) != 2:  # §5 refusal 13a
        subs = crs.sub_crs_list
        horizontal_code = subs[0].to_epsg() if subs else None
        part = "" if horizontal_code is None else f"; horizontal part EPSG {horizontal_code}"
        raise GeoTiffError(
            f"ProjectedCSTypeGeoKey (3072) = {code} is a {crs.type_name} with "
            f"{len(crs.axis_info)} axes{part}. Put the horizontal CRS in 3072 and the "
            "vertical CRS in VerticalGeoKey (4096)"
        )
    linear = geokeys.get("ProjLinearUnitsGeoKey")
    if linear is not None and int(linear) != METRE:
        raise GeoTiffError(f"ProjLinearUnitsGeoKey (3076) = {int(linear)}; only metres ({METRE})")
    horizontal = [a for a in crs.axis_info if a.direction not in ("up", "down")]
    if any(a.unit_conversion_factor != 1.0 for a in horizontal):
        units = sorted({a.unit_name for a in horizontal})
        raise GeoTiffError(f"ProjectedCSTypeGeoKey (3072) = {code} has axes in {units}, not metres")
    return code


def _refuse_through_2048(geographic: Any) -> NoReturn:
    """3072 is absent. A 2048 code is never accepted; pick the refusal that is true."""
    if geographic is None:
        raise GeoTiffError("ProjectedCSTypeGeoKey (3072) is absent: no CRS")
    code = int(geographic)
    crs = _resolve(code)
    consulted = f"ProjectedCSTypeGeoKey (3072) is absent; GeographicTypeGeoKey (2048) = {code}"
    if crs is None:
        raise GeoTiffError(f"{consulted} is not a resolvable EPSG code either: no CRS")
    if crs.is_geographic:
        raise GeoTiffError(f"{consulted} is a {crs.type_name}; a projected CRS is required")
    raise GeoTiffError(f"{consulted} is a {crs.type_name}; a projected CRS belongs in 3072")


def _resolve(code: int) -> pyproj.CRS | None:
    if code == USER_DEFINED:
        return None
    try:
        return pyproj.CRS.from_epsg(code)
    except CRSError:
        return None


def _nodata(
    page: tifffile.TiffPage, dtype: np.dtype[Any], caller: float | None
) -> tuple[float | None, NodataSource]:
    """§6: parse tag 42113's text, check both sentinels, then compare them."""
    tag: float | None = None
    if 42113 in page.tags:
        text = str(page.tags[42113].value)
        where = f"GDAL_NODATA (42113) = {text!r}"
        try:
            tag = None if "_" in text else float(text)  # float() accepts "1_0"
        except ValueError:
            tag = None
        if tag is None:
            raise GeoTiffError(f"{where} is not a number")
        _representable(tag, dtype, where)
    if caller is not None:
        _representable(float(caller), dtype, f"nodata= argument {float(caller)!r}")
        if tag is not None and not _same(tag, float(caller)):
            raise GeoTiffError(f"{where} contradicts the nodata= argument {float(caller)!r}")
    origin: NodataSource = "tag" if tag is not None else "caller"
    value = tag if tag is not None else caller
    if value is None:
        return None, "absent"
    return (None if math.isnan(value) else float(value)), origin


def _representable(value: float, dtype: np.dtype[Any], where: str) -> None:
    """Refuse a sentinel no cell of the *file* dtype can equal (§6, problem 5)."""
    if dtype.kind in "iu":
        info = np.iinfo(dtype)
        ok = math.isfinite(value) and value.is_integer() and info.min <= value <= info.max
    else:
        with np.errstate(over="ignore"):
            ok = math.isnan(value) or (math.isfinite(value) and float(dtype.type(value)) == value)
    if not ok:
        raise GeoTiffError(f"{where} cannot be held by a {dtype} cell")


def _same(a: float, b: float) -> bool:
    return a == b or (math.isnan(a) and math.isnan(b))


__all__ = ["decode_dem"]
