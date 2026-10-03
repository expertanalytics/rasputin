"""Geographic micro-GeoTIFFs and the independent oracles for increment 15c-2.

`docs/increments/15c-geographic-dem.md`, "Tests for @tester" (G1-G10). Test
support only; nothing in `src_python` uses it.

The builders write GeoKeys the way a real file carries them: SHORT values
inline in the GeoKeyDirectoryTag (34735), DOUBLE values in
GeoDoubleParamsTag (34736) and ASCII values, `|`-terminated, in
GeoAsciiParamsTag (34737). `geotiff_fixtures.micro_tiff` writes inline SHORTs
only, which is enough for a projected file but not for the OpenTopography
COG's key set ("Measured by @architect": 2049, 2057 and 2059 beside 2048).

The oracles share no code with `tin_engine`: positions come from pyproj's
own `Transformer` (not `crs.reprojector`), membership from a barycentric test
over every output triangle, and the plane from the output file's vertices.
"""

from __future__ import annotations

import io
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt
import shapely
import tifffile
from pyproj import Transformer
from shapely.geometry import Polygon

#: ANADEM's spacing in degrees, both axes (15-dem-mosaic.md B2; the COG's header).
ANADEM_STEP = 0.00026949458523585647
#: Copernicus GLO-30: one arc-second.
GLO30_STEP = 1.0 / 3600.0
NODATA = -9999.0
NODATA_TEXT = "-9999"

#: Upper-left node of the synthetic Velhas-area tiles (lon, lat), inside UTM 23S.
LON0, LAT0 = -44.0, -19.0

# GeoKey ids and values (GeoTIFF 1.1, OGC 19-008r4).
GT_MODEL_TYPE = 1024
GT_RASTER_TYPE = 1025
GEOGRAPHIC_TYPE = 2048
GEOG_CITATION = 2049
GEOG_ANGULAR_UNITS = 2054
GEOG_SEMI_MAJOR_AXIS = 2057
GEOG_INV_FLATTENING = 2059
PROJECTED_CS_TYPE = 3072
MODEL_PROJECTED, MODEL_GEOGRAPHIC = 1, 2
PIXEL_IS_AREA, PIXEL_IS_POINT = 1, 2
RADIAN, DEGREE, GRAD = 9101, 9102, 9105
GEOKEY_DIRECTORY, GEO_DOUBLE_PARAMS, GEO_ASCII_PARAMS = 34735, 34736, 34737
MODEL_PIXEL_SCALE, MODEL_TIEPOINT, GDAL_NODATA = 33550, 33922, 42113

#: GRS80, as the COG's 2057 and 2059 carry it.
GRS80_A, GRS80_RF = 6378137.0, 298.257222101


def geokey_tags(
    shorts: Mapping[int, int],
    doubles: Mapping[int, float] | None = None,
    texts: Mapping[int, str] | None = None,
) -> list[tuple[int, str, int, Any, bool]]:
    """The 34735/34736/34737 extratags for these keys, sorted by key id."""
    entries: list[tuple[int, int, int, int]] = []
    double_values: list[float] = []
    ascii_text = ""
    for key, value in (shorts or {}).items():
        entries.append((key, 0, 1, int(value)))
    for key, number in (doubles or {}).items():
        entries.append((key, GEO_DOUBLE_PARAMS, 1, len(double_values)))
        double_values.append(float(number))
    for key, text in (texts or {}).items():
        entries.append((key, GEO_ASCII_PARAMS, len(text) + 1, len(ascii_text)))
        ascii_text += text + "|"
    directory = [1, 1, 0, len(entries)]
    for entry in sorted(entries):
        directory += list(entry)
    tags: list[tuple[int, str, int, Any, bool]] = [
        (GEOKEY_DIRECTORY, "H", len(directory), directory, True)
    ]
    if double_values:
        tags.append((GEO_DOUBLE_PARAMS, "d", len(double_values), tuple(double_values), True))
    if ascii_text:
        tags.append((GEO_ASCII_PARAMS, "s", 0, ascii_text, True))
    return tags


def geographic_keys(
    epsg: int = 4326, *, area: bool = True, **changes: int | None
) -> dict[int, int]:
    """ModelTypeGeographic, the registration, 2048 = `epsg` and 2054 = degrees;
    `changes` maps a key's name in this module (e.g. `GEOG_ANGULAR_UNITS`) to a
    new value, `None` deleting it."""
    keys: dict[int, int | None] = {
        GT_MODEL_TYPE: MODEL_GEOGRAPHIC,
        GT_RASTER_TYPE: PIXEL_IS_AREA if area else PIXEL_IS_POINT,
        GEOGRAPHIC_TYPE: epsg,
        GEOG_ANGULAR_UNITS: DEGREE,
    }
    for name, value in changes.items():
        keys[globals()[name]] = value
    return {k: v for k, v in keys.items() if v is not None}


def tiff(
    array: npt.NDArray[Any],
    *,
    x0: float,
    y0: float,
    step_x: float,
    step_y: float | None = None,
    shorts: Mapping[int, int],
    doubles: Mapping[int, float] | None = None,
    texts: Mapping[int, str] | None = None,
    nodata: str | None = NODATA_TEXT,
    overview: bool = False,
    **write_kwargs: Any,
) -> io.BytesIO:
    """A GeoTIFF in memory whose tie point (0, 0) sits at model `(x0, y0)`."""
    dy = step_x if step_y is None else step_y
    tags = [
        (MODEL_TIEPOINT, "d", 6, (0.0, 0.0, 0.0, x0, y0, 0.0), True),
        (MODEL_PIXEL_SCALE, "d", 3, (step_x, dy, 0.0), True),
        *geokey_tags(shorts, doubles, texts),
    ]
    if nodata is not None:
        tags.append((GDAL_NODATA, "s", 0, nodata, True))
    stream = io.BytesIO()
    with tifffile.TiffWriter(stream) as writer:
        writer.write(array, extratags=tags, **write_kwargs)
        if overview:
            writer.write(np.ascontiguousarray(array[::2, ::2]), subfiletype=1, **write_kwargs)
    stream.seek(0)
    return stream


def rough(rows: int, cols: int, *, seed: int = 15) -> npt.NDArray[np.float32]:
    """Terrain with relief at two cells' wavelength, so a resampled grid is
    metres off its source between nodes. float32, as ANADEM and GLO-30 are."""
    r, c = np.indices((rows, cols), dtype=np.float64)
    noise = np.random.default_rng(seed).normal(0.0, 1.5, (rows, cols))
    z = 700 + 25 * np.sin(c / 3.1) * np.cos(r / 2.7) + 9 * np.sin((r + 2 * c) / 1.3) + noise
    return z.astype(np.float32)


def geographic_tile_tiff(
    array: npt.NDArray[Any],
    *,
    lon0: float = LON0,
    lat0: float = LAT0,
    step: float = ANADEM_STEP,
    epsg: int = 4326,
    **write_kwargs: Any,
) -> io.BytesIO:
    """A point-registered geographic tile whose node (0, 0) is `(lon0, lat0)`."""
    return tiff(
        array,
        x0=lon0,
        y0=lat0,
        step_x=step,
        shorts=geographic_keys(epsg, area=False),
        **write_kwargs,
    )


def projected_tiff(
    array: npt.NDArray[Any], *, x0: float, y0: float, step: float, epsg: int
) -> io.BytesIO:
    """A point-registered projected tile in metres, node (0, 0) at `(x0, y0)`."""
    shorts = {GT_MODEL_TYPE: MODEL_PROJECTED, GT_RASTER_TYPE: PIXEL_IS_POINT}
    shorts[PROJECTED_CS_TYPE] = epsg
    return tiff(array, x0=x0, y0=y0, step_x=step, shorts=shorts)


def nodes(
    x0: float, y0: float, dx: float, dy: float, shape: tuple[int, int]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Every node's model `(x, y)`, row-major, for a north-up grid."""
    r, c = np.indices(shape, dtype=np.float64)
    return (x0 + c * dx).ravel(), (y0 - r * dy).ravel()


def project(
    src: str, dst: str, x: npt.NDArray[np.float64], y: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """pyproj's own `always_xy` transform, as an `(N, 2)` array: the oracle's."""
    t = Transformer.from_crs(src, dst, always_xy=True)
    px, py = t.transform(x, y)
    return np.column_stack([np.asarray(px, np.float64), np.asarray(py, np.float64)])


def project_ring(
    src: str, dst: str, ring: Sequence[tuple[float, float]]
) -> list[tuple[float, float]]:
    xy = project(src, dst, np.array([p[0] for p in ring]), np.array([p[1] for p in ring]))
    return [(float(a), float(b)) for a, b in xy]


@dataclass(frozen=True)
class SourceCheck:
    """What the independent final check found at the source nodes."""

    nodes: int
    outside_mesh: int
    over: int
    over_interior: int
    over_strip: int
    worst: float


def source_check(
    points: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.integer[Any]],
    xy: npt.NDArray[np.float64],
    z: npt.NDArray[np.float64],
    domain: Polygon,
    tolerance: float,
    strip: float,
    slack: float = 1e-4,
) -> SourceCheck:
    """J2 at the source nodes, by brute force, in the mesh's own CRS.

    `xy`, `z`: valid source nodes already projected; only those inside
    `domain` are checked. Every (node, triangle) pair whose closed triangle
    holds the node (barycentric, 1e-12 slack, so a node on an edge is tested
    in both triangles) gives `|plane - z|`, the plane from the mesh's own
    vertices. A node is over when that exceeds `tolerance + slack` (the store
    keeps a node's position to 2 um, D4). `strip` splits the count: nodes
    within that distance of the domain's boundary are the edge strip
    (Surprise 3). A node inside the domain but in no triangle is counted in
    `outside_mesh`; within 1e-2 m of the boundary that is the noder's snap.
    """
    inside = shapely.contains_xy(domain, xy[:, 0], xy[:, 1])
    px, py, pz = xy[inside, 0], xy[inside, 1], z[inside]
    in_strip = shapely.contains_xy(domain.boundary.buffer(strip), px, py)
    origin = points[:, :2].min(axis=0)
    v = points[:, :2] - origin
    a, b, c = v[triangles[:, 0]], v[triangles[:, 1]], v[triangles[:, 2]]
    za, zb, zc = points[triangles[:, 0], 2], points[triangles[:, 1], 2], points[triangles[:, 2], 2]
    two_a = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    assert (np.abs(two_a) > 0).all(), "a zero-area output triangle"
    worst_error = np.full(px.size, -np.inf)
    for lo in range(0, px.size, 128):
        x = px[lo : lo + 128, None] - origin[0]
        y = py[lo : lo + 128, None] - origin[1]

        def weight(p: Any, q: Any, x: Any = x, y: Any = y) -> Any:
            return (
                (q[:, 0] - p[:, 0]) * (y - p[:, 1]) - (q[:, 1] - p[:, 1]) * (x - p[:, 0])
            ) / two_a

        wa, wb, wc = weight(b, c), weight(c, a), weight(a, b)
        holds = (wa >= -1e-12) & (wb >= -1e-12) & (wc >= -1e-12)
        error = np.abs(wa * za + wb * zb + wc * zc - pz[lo : lo + 128, None])
        worst_error[lo : lo + 128] = np.where(holds, error, -np.inf).max(axis=1)
    located = np.isfinite(worst_error)
    near_edge = shapely.contains_xy(domain.boundary.buffer(1e-2), px, py)
    assert (located | near_edge).all(), "a source node well inside the domain is in no triangle"
    over = located & (worst_error > tolerance + slack)
    return SourceCheck(
        nodes=int(px.size),
        outside_mesh=int((~located).sum()),
        over=int(over.sum()),
        over_interior=int((over & ~in_strip).sum()),
        over_strip=int((over & in_strip).sum()),
        worst=float(worst_error[located].max()) if located.any() else 0.0,
    )
