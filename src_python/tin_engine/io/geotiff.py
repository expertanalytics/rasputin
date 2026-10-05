"""A GeoTIFF stream to a `DemTile`, or a `GeoTiffError` that says why not.

Increment 11, `docs/increments/11-raster-ingestion.md`. The reader's job is
refusal (ruling 3): every refusal in §5 raises `GeoTiffError` naming the tag or
GeoKey by number and name and the file's value. None is a bare `assert`,
because `python -O` strips those (prior art §5.9).

The module takes a binary stream and never a path (§2), and imports neither
`_core` nor anything first-party except `io.models` and `crs` (ruling 1).
GeoKeys come from tifffile's own `geotiff_metadata` (ruling 5); the CRS is
resolved through `pyproj.CRS.from_epsg` and tested on the constructed CRS
(rulings 6, 7), or, for 3072 = 32767, built from its GeoKeys and named by the
EPSG code it matches (`docs/increments/geotiff-crs-by-parameters.md`).
No transform happens here, so there is no `Transformer` (ruling 8).
"""

from __future__ import annotations

import functools
import importlib.util
import math
from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any, BinaryIO, Literal

import numpy as np
import pyproj
import tifffile
from pyproj.exceptions import CRSError

from .. import crs as crs_rules
from .models import DemTile, GeoTiffError, RasterMeta

USER_DEFINED = 32767
METRE = 9001
DEGREE = 9102
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

#: The words every refusal of a CRS given by parameters begins with.
P = f"ProjectedCSTypeGeoKey (3072) = {USER_DEFINED} (a CRS given by parameters)"

#: 3075 to the EPSG method and its parameters, each (EPSG name, EPSG code, GeoKey,
#: unit), in the EPSG method's own order: PROJ's candidate search compares
#: parameters by position. One parameter per line, packed so the table reads as one.
# fmt: off
METHODS: dict[int, tuple[str, int, tuple[tuple[str, int, str, str], ...]]] = {
    1: ("Transverse Mercator", 9807, (
        ("Latitude of natural origin", 8801, "ProjNatOriginLatGeoKey", "degree"),
        ("Longitude of natural origin", 8802, "ProjNatOriginLongGeoKey", "degree"),
        ("Scale factor at natural origin", 8805, "ProjScaleAtNatOriginGeoKey", "unity"),
        ("False easting", 8806, "ProjFalseEastingGeoKey", "metre"),
        ("False northing", 8807, "ProjFalseNorthingGeoKey", "metre"))),
    8: ("Lambert Conic Conformal (2SP)", 9802, (
        ("Latitude of false origin", 8821, "ProjFalseOriginLatGeoKey", "degree"),
        ("Longitude of false origin", 8822, "ProjFalseOriginLongGeoKey", "degree"),
        ("Latitude of 1st standard parallel", 8823, "ProjStdParallel1GeoKey", "degree"),
        ("Latitude of 2nd standard parallel", 8824, "ProjStdParallel2GeoKey", "degree"),
        ("Easting at false origin", 8826, "ProjFalseOriginEastingGeoKey", "metre"),
        ("Northing at false origin", 8827, "ProjFalseOriginNorthingGeoKey", "metre"))),
}
# fmt: on
#: GeoKeys that would name the datum other than through 2048 (row 5).
DATUM_KEYS = ("GeogGeodeticDatumGeoKey", "GeogPrimeMeridianGeoKey", "GeogEllipsoidGeoKey")


def read_meta(source: BinaryIO, *, nodata: float | None = None) -> RasterMeta:
    """The header phase of `decode_dem` alone (increment 15a, R3): no pixel is read.

    Every refusal `decode_dem` makes from the header fires here with the same
    message, because both run the same `_header`. A file whose pixels are
    damaged still gives its `RasterMeta`; `decode_dem` then refuses it at the
    "pixel data" stage. `source` is read but not closed.
    """
    return read_header(source, nodata=nodata)[0]


def read_header(
    source: BinaryIO, *, nodata: float | None = None
) -> tuple[RasterMeta, np.dtype[Any]]:
    """`read_meta`, and the dtype `decode_dem` will return: the file's through
    `PROMOTION` (15a S2, so a mosaic's cap can count the canvas it will build)."""
    with _tiff(source, nodata) as tif:
        meta, dtype = _header(tif, nodata)
    return meta, PROMOTION[dtype]


def read_page(
    source: BinaryIO, *, nodata: float | None = None, geographic: bool = False
) -> tuple[RasterMeta, np.dtype[Any], tifffile.TiffPage]:
    """`read_header`, and the full-resolution page itself, for `io.cog` to
    decode block by block (23a-1). The same `_header`, so every refusal is the
    same. A prefix of the file that ends before its first block is enough;
    the page decodes blocks from their bytes alone after the file is closed.
    `geographic` is 23a-2's flag, kept for its callers: since 15c-2 every
    entry point reads a geographic 2D CRS (D6), so it changes nothing."""
    with _tiff(source, nodata) as tif:
        meta, dtype = _header(tif, nodata)
        page = tif.pages.first
    return meta, PROMOTION[dtype], page


def decode_dem(source: BinaryIO, *, nodata: float | None = None) -> DemTile:
    """Decode page 0 of the GeoTIFF in `source` into a `DemTile`.

    `nodata` is the caller asserting a sentinel (§6). It must be representable
    in the file's dtype, and must agree with tag 42113 when that is present. A
    `bool` is refused with `TypeError`: it is a caller's mistake, not a file's.

    Raises `GeoTiffError` for every §5 refusal, and (§14, choice C) when
    tifffile or a codec fails on the file; that message names the stage, the
    original type and text, and chains the original. Only the calls into
    tifffile are wrapped: `MemoryError` passes through, and the reader's own
    code, refusals included, runs outside the wrapped regions, so its bugs
    surface as themselves.

    The whole page is decoded, so memory is O(rows x cols), briefly doubled
    while `DemTile` takes its read-only copy (§7). `source` is read but not
    closed; closing it is the caller's.
    """
    with _tiff(source, nodata) as tif:
        meta, dtype = _header(tif, nodata)
        with _stage("pixel data"):
            array = tif.pages.first.asarray().astype(PROMOTION[dtype], copy=False)
    return DemTile(meta=meta, array=array)


@contextmanager
def _tiff(source: BinaryIO, nodata: float | None) -> Iterator[tifffile.TiffFile]:
    """The caller-argument check and the TIFF structure, shared by both phases."""
    if isinstance(nodata, bool | np.bool_):
        raise TypeError(f"nodata= must be a number or None, not {type(nodata).__name__}")
    with _stage("TIFF structure"):
        tif = tifffile.TiffFile(source)
    with tif:
        yield tif


def _header(tif: tifffile.TiffFile, nodata: float | None) -> tuple[RasterMeta, np.dtype[Any]]:
    """Everything before the pixels (R3): the §5 header refusals, then the meta.

    Returns the file dtype too, which `decode_dem` promotes by. `rows` and
    `cols` come from ImageLength and ImageWidth, which `_check_page` has
    bounded; `DemTile` checks them against the decoded array's shape.
    """
    with _stage("TIFF structure"):
        page, pages = tif.pages.first, tuple(tif.pages)  # parses every IFD
    _single_page(pages)
    dtype = _check_page(page)
    tie, scale = _georeferencing(page)
    with _stage("GeoKey directory"):
        geokeys: dict[str, Any] = tif.geotiff_metadata or {}
    x_min, y_max, delta_x, delta_y, area = _placement(tie, scale, geokeys)
    # 15c-2, D6: with 3072 absent, a geographic 2D CRS in degrees from 2048.
    projected = geokeys.get("ProjectedCSTypeGeoKey")
    geographic = projected is None
    if projected is None:
        epsg = _geographic_epsg(geokeys)
    elif int(projected) == USER_DEFINED:
        epsg = _parametric_epsg(geokeys)
    else:
        epsg = _projected_epsg(geokeys)
    vertical = geokeys.get("VerticalUnitsGeoKey")
    if vertical is not None and int(vertical) != METRE:
        raise GeoTiffError(
            f"VerticalUnitsGeoKey (4099) = {int(vertical)}; only metres ({METRE}) are read"
        )
    sentinel, source_of = _nodata(page, dtype, nodata)
    meta = RasterMeta(
        x_min=x_min,
        y_max=y_max,
        delta_x=delta_x,
        delta_y=delta_y,
        cols=int(page.imagewidth),
        rows=int(page.imagelength),
        epsg=epsg,
        geographic=geographic,
        nodata=sentinel,
        nodata_source=source_of,
        pixel_is_area=area,
        vertical_unit_assumed=vertical is None,
    )
    return meta, dtype


@contextmanager
def _stage(stage: str) -> Iterator[None]:
    """§14, choice C: a failure inside tifffile or a codec becomes `GeoTiffError`.

    Wraps call sites, not types, so the reader's own code stays outside and its
    bugs surface as bugs. No refusal runs inside a stage. The type is named
    with its top-level module, so `struct.error` does not read as "error".
    """
    try:
        yield
    except MemoryError:
        raise
    except Exception as error:
        kind = type(error)
        module = kind.__module__.partition(".")[0]
        name = kind.__qualname__ if module == "builtins" else f"{module}.{kind.__qualname__}"
        raise GeoTiffError(
            f"{stage}: tifffile could not decode the file: {name}: {error}"
        ) from error


def _single_page(pages: tuple[tifffile.TiffPage | tifffile.TiffFrame, ...]) -> None:
    """Refuse a second full-resolution page (§5 refusal 8, 0-based)."""
    for index, extra in enumerate(pages[1:], start=1):
        where = f"on page {index} of {len(pages)}: only page 0 may be full resolution"
        if not isinstance(extra, tifffile.TiffPage):
            # A TiffFrame has no subfiletype to consult (round 3).
            raise GeoTiffError(f"NewSubfileType (254) could not be read (a frame) {where}")
        if not extra.is_reduced:
            raise GeoTiffError(f"NewSubfileType (254) = {int(extra.subfiletype)} {where}")


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
            f"ModelTiepointTag (33922) has {_count(tie)}; exactly 6 (one tie point) are read"
        )
    if len(tie) < 6:
        raise GeoTiffError(f"ModelTiepointTag (33922) has {_count(tie)}; need 6")
    if len(scale) != 3:
        raise GeoTiffError(f"ModelPixelScaleTag (33550) has {_count(scale)}; need 3")
    if not all(math.isfinite(v) for v in (*tie[:2], *tie[3:5])):
        raise GeoTiffError(f"ModelTiepointTag (33922) = {tie}: need finite I, J, X, Y")
    if not all(math.isfinite(v) and v > 0 for v in scale[:2]):
        raise GeoTiffError(f"ModelPixelScaleTag (33550) = {scale}: need finite positive X, Y")
    if scale[2] != 0.0:
        raise GeoTiffError(f"ModelPixelScaleTag (33550) ScaleZ = {scale[2]}; must be 0")
    return tie, scale


def _count(values: list[float]) -> str:
    return "1 value" if len(values) == 1 else f"{len(values)} values"


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
    code = int(geokeys["ProjectedCSTypeGeoKey"])
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


def _parametric_epsg(geokeys: dict[str, Any]) -> int:
    """3072 = 32767: the CRS built from the GeoKeys on the datum 2048 names, read
    as the EPSG code it matches, or the refusal of section 4's row that fires."""
    projection = geokeys.get("ProjectionGeoKey")
    if projection is not None and int(projection) != USER_DEFINED:
        raise GeoTiffError(
            f"{P}, but ProjectionGeoKey (3074) = {int(projection)} names its projection by "
            "EPSG code, which is not read; only a projection given by its parameters is"
        )
    transform = geokeys.get("ProjCoordTransGeoKey")
    if transform is None:
        raise GeoTiffError(
            f"{P} has no ProjCoordTransGeoKey (3075), so its projection method is unknown"
        )
    if int(transform) not in METHODS:
        raise GeoTiffError(
            f"{P} uses ProjCoordTransGeoKey (3075) = {int(transform)} "
            f"({getattr(transform, 'name', 'unknown')}); only transverse Mercator (1) and "
            "Lambert conic conformal with two standard parallels (8) are read"
        )
    method, _, parameters = METHODS[int(transform)]
    g, base = _parametric_base(geokeys)
    linear = geokeys.get("ProjLinearUnitsGeoKey")
    if linear is None:
        raise GeoTiffError(
            f"{P}: ProjLinearUnitsGeoKey (3076) is absent, so the unit of its false easting "
            "and northing is unknown"
        )
    if int(linear) != METRE:
        raise GeoTiffError(
            f"{P}: ProjLinearUnitsGeoKey (3076) = {int(linear)}; only metres ({METRE})"
        )
    for _, _, key, _ in parameters:
        if geokeys.get(key) is None:
            number = int(tifffile.TIFF.GEO_KEYS[key])
            raise GeoTiffError(f"{P}: {method} needs {key} ({number}), which is absent")
    values = tuple(float(geokeys[key]) for _, _, key, _ in parameters)
    matches = _epsg_matches_cached(g, int(transform), values)
    if not matches:
        raise GeoTiffError(
            f"{P}: it is {method} on EPSG:{g} ({base.name}), and no EPSG projected CRS on that "
            "datum has these parameters, so it cannot be named. Reproject the file to a CRS "
            "with an EPSG code first"
        )
    if not all(crs_rules.same_crs(f"EPSG:{matches[0]}", f"EPSG:{n}") for n in matches[1:]):
        codes = ", ".join(f"EPSG:{n}" for n in matches)
        raise GeoTiffError(
            f"{P}: its parameters match {codes}, which are not the same CRS, so it cannot be named"
        )
    return matches[0]


def _parametric_base(geokeys: dict[str, Any]) -> tuple[int, pyproj.CRS]:
    """Rows 4 to 8: the datum from 2048 alone, and no key that says otherwise."""
    geographic = geokeys.get("GeographicTypeGeoKey")
    g = None if geographic is None else int(geographic)
    base = None if g is None else _resolve(g)
    if g is None or base is None or not base.is_geographic or len(base.axis_info) != 2:
        shown = (
            "is absent" if g is None
            else f"= {g}, which is not an EPSG code" if base is None
            else f"= {g}, a {base.type_name}"
        )  # fmt: skip
        raise GeoTiffError(
            f"{P}: GeographicTypeGeoKey (2048) {shown}; its datum must be given there as an "
            "EPSG geographic CRS code"
        )
    for key in DATUM_KEYS:
        if key in geokeys:
            raise GeoTiffError(
                f"{P}: {key} ({int(tifffile.TIFF.GEO_KEYS[key])}) is present; the datum is "
                "read only from GeographicTypeGeoKey (2048)"
            )
    ellipsoid, meridian = base.ellipsoid, base.prime_meridian
    assert ellipsoid is not None and meridian is not None  # a geographic CRS has both
    greenwich_east = math.degrees(meridian.longitude * float(meridian.unit_conversion_factor))
    for key, expected, what in (
        ("GeogSemiMajorAxisGeoKey", ellipsoid.semi_major_metre, ellipsoid.name),
        ("GeogSemiMinorAxisGeoKey", ellipsoid.semi_minor_metre, ellipsoid.name),
        ("GeogInvFlatteningGeoKey", ellipsoid.inverse_flattening, ellipsoid.name),
        ("GeogPrimeMeridianLongGeoKey", greenwich_east, "prime meridian"),
    ):
        value = geokeys.get(key)
        # 1e-10 relative is PROJ's; the absolute 1e-10 matters only at a prime meridian of 0.
        if value is not None and not math.isclose(
            float(value), expected, rel_tol=1e-10, abs_tol=1e-10
        ):
            raise GeoTiffError(
                f"{P}: {key} ({int(tifffile.TIFF.GEO_KEYS[key])}) = {float(value)!r} disagrees "
                f"with EPSG:{g}'s {what} ({expected!r})"
            )
    towgs84 = geokeys.get("GeogTOWGS84GeoKey")
    if towgs84 is not None:
        raise GeoTiffError(
            f"{P}: GeogTOWGS84GeoKey (2062) = {towgs84} gives the file's own datum shift to "
            "WGS 84, which is not read; the datum is read only from GeographicTypeGeoKey (2048)"
        )
    units = geokeys.get("GeogAngularUnitsGeoKey")
    if units is not None and int(units) != DEGREE:
        raise GeoTiffError(
            f"{P}: GeogAngularUnitsGeoKey (2054) = {int(units)}; only degrees ({DEGREE})"
        )
    return g, base


@functools.cache
def _epsg_matches_cached(
    geographic: int, transform: int, values: tuple[float, ...]
) -> tuple[int, ...]:
    """`crs.epsg_matches` of the CRS built from these GeoKeys, once per distinct key set
    (a mosaic reads every tile's header)."""
    return crs_rules.epsg_matches(_projected_crs(geographic, transform, values))


def _projected_crs(geographic: int, transform: int, values: tuple[float, ...]) -> pyproj.CRS:
    """PROJJSON: the base by its EPSG code, the method and parameters with their EPSG
    names and codes, and axes east then north in metres (GeoTIFF 1.1's fixed order)."""
    method, method_code, parameters = METHODS[transform]

    def epsg(name: str, code: int) -> dict[str, Any]:
        return {"name": name, "id": {"authority": "EPSG", "code": code}}

    return pyproj.CRS.from_json_dict({
        "type": "ProjectedCRS", "name": "unknown",
        "base_crs": pyproj.CRS.from_epsg(geographic).to_json_dict(),
        "conversion": {"name": "unknown", "method": epsg(method, method_code), "parameters": [
            {**epsg(name, code), "value": value, "unit": unit}
            for (name, code, _, unit), value in zip(parameters, values, strict=True)
        ]},
        "coordinate_system": {"subtype": "Cartesian", "axis": [
            {"name": axis, "abbreviation": axis[0], "direction": direction, "unit": "metre"}
            for axis, direction in (("Easting", "east"), ("Northing", "north"))
        ]},
    })  # fmt: skip


def _geographic_epsg(geokeys: dict[str, Any]) -> int:
    """3072 is absent: a geographic 2D CRS in degrees from 2048 (15c-2, D6), or
    the refusal that is true of what 2048 holds."""
    geographic = geokeys.get("GeographicTypeGeoKey")
    if geographic is None:
        raise GeoTiffError("ProjectedCSTypeGeoKey (3072) is absent: no CRS")
    code = int(geographic)
    crs = _resolve(code)
    consulted = f"ProjectedCSTypeGeoKey (3072) is absent; GeographicTypeGeoKey (2048) = {code}"
    if crs is None:
        raise GeoTiffError(f"{consulted} is not a resolvable EPSG code either: no CRS")
    if not crs.is_geographic:
        raise GeoTiffError(f"{consulted} is a {crs.type_name}; a projected CRS belongs in 3072")
    if len(crs.axis_info) != 2:
        raise GeoTiffError(f"{consulted} is a {crs.type_name}; only a geographic 2D CRS is read")
    units = geokeys.get("GeogAngularUnitsGeoKey")
    if units is not None and int(units) != DEGREE:
        raise GeoTiffError(f"GeogAngularUnitsGeoKey (2054) = {int(units)}; only degrees ({DEGREE})")
    if any(a.unit_name != "degree" for a in crs.axis_info):
        found = sorted({a.unit_name for a in crs.axis_info})
        raise GeoTiffError(f"{consulted} has axes in {found}, not degrees")
    return code


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


__all__ = ["decode_dem", "read_header", "read_meta", "read_page"]
