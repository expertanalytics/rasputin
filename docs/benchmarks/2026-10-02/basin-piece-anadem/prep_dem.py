"""Fetch the DEM window under a domain, then project and resample it onto a square grid.

A measurement script, not production code (ROADMAP basin item 2.1). It does by
hand, for one piece, what increment 15c is to build: a geographic DEM projected
onto the target TIN CRS and resampled onto a square Cartesian grid there (Q6,
ruled 2026-09-30). Nothing imports it.

    python prep_dem.py fetch    SOURCE OUTLINE DATA_DIR   # source window, the source's own CRS
    python prep_dem.py resample SOURCE OUTLINE DATA_DIR [--spacing 30] [--epsg 31983]
    python prep_dem.py check    SOURCE OUTLINE DATA_DIR   # B6-style: resampled vs source

SOURCE is ``anadem`` (ANADEM v1, OpenTopography's single continental COG,
EPSG:4674, by HTTP range reads of its first 8 MiB, which hold every IFD, and
then of only the 512 x 512 blocks under the window) or ``glo30`` (Copernicus GLO-30 DSM, 1-degree
COGs on AWS, by range reads of the 1024 x 1024 blocks under the window). The
window is the domain's projected box plus a margin, mapped back to longitude
and latitude and grown by two source cells.

The target grid is node-registered (PixelIsPoint) in EPSG:31983 (SIRGAS 2000 /
UTM 23S) by default, at 30 m, with node coordinates on multiples of the
spacing. z is bilinear in the source's own (longitude, latitude) index space;
a target node whose four source neighbours are not all valid is NoData.
ANADEM is SIRGAS 2000 geographic (EPSG:4674), the target's own datum; GLO-30
is WGS 84 (EPSG:4326), which PROJ takes to SIRGAS 2000 by a null
transformation. The source CRS comes from the ``SOURCE_EPSG`` table, which
``fetch`` checks against the COG's GeoKeys for ANADEM. ``fetch`` writes the
source window in that CRS, and ``resample`` and ``check`` take it from the
same table, not from the window file.
"""

from __future__ import annotations

import io
import json
import math
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pyproj
import shapely
import shapely.ops
import tifffile
from shapely.geometry import shape

NODATA = -9999.0
# ANADEM's own host (metadados.snirh.gov.br) answers 403; OpenTopography's copy.
ANADEM_URL = "https://opentopography.s3.sdsc.edu/raster/ANADEM/ANADEM_be/anadem_v1_compressed_COG.tif"
ANADEM_PREFIX = 8 * 2**20  # every IFD lies in it (docs/increments/23-probes/anadem_cog.py)
SOURCE_EPSG = {"anadem": 4674, "glo30": 4326}
GLO30_URL = (
    "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_{ns}{lat:02d}_00_"
    "{ew}{lon:03d}_00_DEM/Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif"
)
MARGIN_CELLS = 4  # target cells beyond the domain's projected box


def _range(url: str, start: int, length: int, tries: int = 5) -> bytes:
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{start + length - 1}"})
    for k in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            if len(data) != length:
                raise OSError(f"{url}: asked {length} bytes at {start}, got {len(data)}")
            return data
        except Exception:  # noqa: BLE001 - public servers; retry, then raise
            if k == tries - 1:
                raise
            time.sleep(3 * (k + 1))
    raise AssertionError


class _Head(io.RawIOBase):
    """Enough of a seekable file for tifffile to parse the header: 64 KiB blocks."""

    def __init__(self, url: str) -> None:
        self.url, self.pos, self.cache = url, 0, {}
        head = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60)
        self.size = int(head.headers["Content-Length"])

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = 0) -> int:
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def readinto(self, buf) -> int:  # type: ignore[no-untyped-def]
        n, out = min(len(buf), self.size - self.pos), bytearray()
        while len(out) < n:
            b, o = divmod(self.pos + len(out), 1 << 16)
            if b not in self.cache:
                lo = b << 16
                self.cache[b] = _range(self.url, lo, min(1 << 16, self.size - lo))
            out += self.cache[b][o : o + n - len(out)]
        buf[:n] = out
        self.pos += n
        return n


def _window_blocks(url: str, r0: int, r1: int, c0: int, c1: int) -> tuple[np.ndarray, dict]:
    """Rows [r0, r1) x cols [c0, c1) of the file's first page, reading only its blocks."""
    tif = tifffile.TiffFile(io.BufferedReader(_Head(url), buffer_size=1 << 16))
    page = tif.pages.first
    th, tw = page.tilelength, page.tilewidth
    across = -(-page.imagewidth // tw)
    out = np.full((r1 - r0, c1 - c0), NODATA, dtype=np.float32)
    jobs = [(br, bc) for br in range(r0 // th, (r1 - 1) // th + 1)
            for bc in range(c0 // tw, (c1 - 1) // tw + 1)]  # fmt: skip

    def one(job: tuple[int, int]) -> tuple[int, int, np.ndarray]:
        br, bc = job
        i = br * across + bc
        raw = _range(url, page.dataoffsets[i], page.databytecounts[i])
        return br, bc, np.asarray(page.decode(raw, i)[0]).reshape(th, tw)

    with ThreadPoolExecutor(8) as pool:
        for br, bc, blk in pool.map(one, jobs):
            rr0, cc0 = br * th, bc * tw
            a0, a1 = max(r0, rr0), min(r1, rr0 + th, page.imagelength)
            b0, b1 = max(c0, cc0), min(c1, cc0 + tw, page.imagewidth)
            out[a0 - r0 : a1 - r0, b0 - c0 : b1 - c0] = blk[a0 - rr0 : a1 - rr0, b0 - cc0 : b1 - cc0]
    nodata = page.tags[42113].value if 42113 in page.tags else None
    info = {"tie": list(page.tags[33922].value), "scale": list(page.tags[33550].value),
            "raster_type": int((tif.geotiff_metadata or {}).get("GTRasterTypeGeoKey", 0)),
            "shape": list(page.shape), "blocks_read": len(jobs), "nodata_tag": nodata}  # fmt: skip
    if nodata is not None and float(nodata) != NODATA:
        out[out == np.float32(float(nodata))] = NODATA
    return out, info


_ANADEM_PAGE: list = []


def _anadem_page():  # type: ignore[no-untyped-def]
    """The full-resolution page of the COG, parsed from one 8 MiB prefix read, once per process.

    tifffile does not raise on a short prefix; it returns pages without
    offsets, so the offsets are counted against the blocks (as the probe does).
    """
    if not _ANADEM_PAGE:
        tif = tifffile.TiffFile(io.BytesIO(_range(ANADEM_URL, 0, ANADEM_PREFIX)))
        page = tif.pages.first
        blocks = -(-page.imagelength // page.tilelength) * -(-page.imagewidth // page.tilewidth)
        if len(page.dataoffsets) != blocks:
            raise SystemExit(f"ANADEM header not within {ANADEM_PREFIX} bytes")
        geo = tif.geotiff_metadata or {}
        if int(geo.get("GTRasterTypeGeoKey", 0)) != 1 or int(geo.get("GeographicTypeGeoKey", 0)) != 4674:
            raise SystemExit(f"ANADEM: not area-registered EPSG:4674: {geo}")
        _ANADEM_PAGE.append(page)
    return _ANADEM_PAGE[0]


def fetch_anadem(lon0: float, lon1: float, lat0: float, lat1: float):
    """Node (r, c) of the result is at (lon_first + c d, lat_first - r d)."""
    page = _anadem_page()
    x, y = page.tags[33922].value[3:5]
    d = page.tags[33550].value[0]
    # Area-registered (B2): node c sits at the centre of cell c.
    c0, c1 = math.floor((lon0 - x) / d - 0.5), math.ceil((lon1 - x) / d - 0.5) + 1
    r0, r1 = math.floor((y - lat1) / d - 0.5), math.ceil((y - lat0) / d - 0.5) + 1
    if c0 < 0 or r0 < 0 or c1 > page.imagewidth or r1 > page.imagelength:
        raise SystemExit("window not inside ANADEM")
    th, tw = page.tilelength, page.tilewidth
    across = -(-page.imagewidth // tw)
    z = np.full((r1 - r0, c1 - c0), NODATA, dtype=np.float32)
    jobs = [(br, bc) for br in range(r0 // th, (r1 - 1) // th + 1)
            for bc in range(c0 // tw, (c1 - 1) // tw + 1)]  # fmt: skip
    nbytes = sum(page.databytecounts[br * across + bc] for br, bc in jobs)

    def one(job: tuple[int, int]) -> tuple[int, int, np.ndarray]:
        br, bc = job
        i = br * across + bc
        raw = _range(ANADEM_URL, page.dataoffsets[i], page.databytecounts[i])
        return br, bc, np.asarray(page.decode(raw, i)[0]).reshape(th, tw)

    with ThreadPoolExecutor(8) as pool:
        for br, bc, blk in pool.map(one, jobs):
            rr0, cc0 = br * th, bc * tw
            a0, a1 = max(r0, rr0), min(r1, rr0 + th, page.imagelength)
            b0, b1 = max(c0, cc0), min(c1, cc0 + tw, page.imagewidth)
            z[a0 - r0 : a1 - r0, b0 - c0 : b1 - c0] = blk[a0 - rr0 : a1 - rr0, b0 - cc0 : b1 - cc0]
    nodata = float(page.tags[42113].value)
    z[(z == np.float32(nodata)) | ~np.isfinite(z)] = NODATA
    info = {"source": ANADEM_URL, "shape": list(page.shape), "tie": [x, y], "spacing_deg": d,
            "raster_type": 1, "blocks_read": len(jobs), "bytes_read": int(nbytes),
            "nodata_tag": nodata}  # fmt: skip
    return z, x + (c0 + 0.5) * d, y - (r0 + 0.5) * d, d, info


def fetch_glo30(lon0: float, lon1: float, lat0: float, lat1: float):
    """1-degree tiles, 3600 x 3600 below 50 degrees, PixelIsPoint, tie at the tile's
    north-west integer corner; tiles do not overlap, so they abut on one lattice."""
    d = 1.0 / 3600.0
    west, north = math.floor(lon0), math.ceil(lat1)
    c0, c1 = math.floor((lon0 - west) / d), math.ceil((lon1 - west) / d) + 1
    r0, r1 = math.floor((north - lat1) / d), math.ceil((north - lat0) / d) + 1
    z = np.full((r1 - r0, c1 - c0), NODATA, dtype=np.float32)
    infos = []
    for tn in range(north, math.floor(lat0), -1):  # tile whose top row is latitude tn
        for tw in range(west, math.ceil(lon1)):
            tr0, tc0 = (north - tn) * 3600, (tw - west) * 3600
            a0, a1 = max(r0, tr0), min(r1, tr0 + 3600)
            b0, b1 = max(c0, tc0), min(c1, tc0 + 3600)
            if a0 >= a1 or b0 >= b1:
                continue
            south = tn - 1  # GLO-30 names a tile by its south-west corner
            url = GLO30_URL.format(ns="S" if south < 0 else "N", lat=abs(south),
                                   ew="W" if tw < 0 else "E", lon=abs(tw))  # fmt: skip
            blk, info = _window_blocks(url, a0 - tr0, a1 - tr0, b0 - tc0, b1 - tc0)
            if info["shape"] != [3600, 3600] or info["raster_type"] != 2:
                raise SystemExit(f"{url}: not a 3600 x 3600 PixelIsPoint tile: {info}")
            z[a0 - r0 : a1 - r0, b0 - c0 : b1 - c0] = blk
            infos.append({"source": url, **info})
    return z, west + c0 * d, north - r0 * d, d, {"tiles": infos}


def _outline(path: Path) -> tuple[shapely.Polygon, str]:
    doc = json.loads(path.read_text())
    crs = doc.get("crs", {}).get("properties", {}).get("name", "EPSG:4326")
    return shape(doc["features"][0]["geometry"]), crs


def _target_box(outline: Path, epsg: int, h: float) -> tuple[float, float, float, float]:
    poly, crs = _outline(outline)
    fwd = pyproj.Transformer.from_crs(crs, f"EPSG:{epsg}", always_xy=True)
    x, y = fwd.transform(*np.asarray(poly.exterior.coords).T)
    m = MARGIN_CELLS * h
    return (math.floor((x.min() - m) / h) * h, math.floor((y.min() - m) / h) * h,
            math.ceil((x.max() + m) / h) * h, math.ceil((y.max() + m) / h) * h)  # fmt: skip


def _source_box(box: tuple[float, float, float, float], epsg: int, grow: float, src_epsg: int):
    inv = pyproj.Transformer.from_crs(f"EPSG:{epsg}", f"EPSG:{src_epsg}", always_xy=True)
    t = np.linspace(0.0, 1.0, 201)
    x0, y0, x1, y1 = box
    xs = np.concatenate([x0 + (x1 - x0) * t, np.full_like(t, x1), x1 - (x1 - x0) * t, np.full_like(t, x0)])
    ys = np.concatenate([np.full_like(t, y0), y0 + (y1 - y0) * t, np.full_like(t, y1), y1 - (y1 - y0) * t])
    lon, lat = inv.transform(xs, ys)
    return lon.min() - grow, lon.max() + grow, lat.min() - grow, lat.max() + grow


def write_geotiff(path: Path, z: np.ndarray, x0: float, y0: float, dx: float, dy: float,
                  epsg: int, geographic: bool) -> None:  # fmt: skip
    """Node-registered GeoTIFF: node (0, 0) at (x0, y0); NoData -9999 in tag 42113."""
    keys = {1024: 2 if geographic else 1, 1025: 2, (2048 if geographic else 3072): epsg}
    shorts = [1, 1, 0, len(keys)]
    for k, v in sorted(keys.items()):
        shorts += [k, 0, 1, v]
    tags = [(33922, "d", 6, (0.0, 0.0, 0.0, x0, y0, 0.0), True),
            (33550, "d", 3, (dx, dy, 0.0), True),
            (34735, "H", len(shorts), shorts, True),
            (42113, "s", 0, "-9999", True)]  # fmt: skip
    tifffile.imwrite(path, z.astype(np.float32), extratags=tags, compression="deflate",
                     predictor=3, tile=(512, 512))  # fmt: skip


def read_geotiff(path: Path) -> tuple[np.ndarray, float, float, float, float]:
    with tifffile.TiffFile(path) as tif:
        p = tif.pages.first
        tie, sc = p.tags[33922].value, p.tags[33550].value
        return p.asarray(), tie[3], tie[4], sc[0], sc[1]


def _names(source: str, outline: Path, data: Path, epsg: int, h: float) -> tuple[Path, Path]:
    stem = outline.name.split("_outline")[0]
    return (data / f"{stem}_{source}_window_epsg{SOURCE_EPSG[source]}.tif",
            data / "derived" / f"{stem}_{source}_epsg{epsg}_{h:g}m.tif")  # fmt: skip


def fetch(source: str, outline: Path, data: Path, epsg: int = 31983, h: float = 30.0) -> None:
    box = _target_box(outline, epsg, h)
    grow = 2.0 / 3600.0
    lon0, lon1, lat0, lat1 = _source_box(box, epsg, grow, SOURCE_EPSG[source])
    t0 = time.perf_counter()
    z, lonf, latf, d, info = (fetch_anadem if source == "anadem" else fetch_glo30)(lon0, lon1, lat0, lat1)
    secs = time.perf_counter() - t0
    win, _ = _names(source, outline, data, epsg, h)
    write_geotiff(win, z, lonf, latf, d, d, SOURCE_EPSG[source], geographic=True)
    meta = {"source": source, "window_lonlat": [lon0, lon1, lat0, lat1], "shape": list(z.shape),
            "first_node": [lonf, latf], "spacing_deg": d, "nodata_nodes": int((z == NODATA).sum()),
            "z_range": [float(z[z != NODATA].min()), float(z[z != NODATA].max())],
            "fetch_s": round(secs, 1), "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "info": info}  # fmt: skip
    win.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    print(json.dumps({k: v for k, v in meta.items() if k != "info"}))
    print(win)


def bilinear(z: np.ndarray, fx: np.ndarray, fy: np.ndarray) -> np.ndarray:
    """z at fractional (col, row); NoData where any of the four neighbours is."""
    c = np.floor(fx).astype(np.int64)
    r = np.floor(fy).astype(np.int64)
    ok = (c >= 0) & (r >= 0) & (c < z.shape[1] - 1) & (r < z.shape[0] - 1)
    c, r = np.where(ok, c, 0), np.where(ok, r, 0)
    tx, ty = fx - c, fy - r
    q = [z[r, c], z[r, c + 1], z[r + 1, c], z[r + 1, c + 1]]
    ok &= np.logical_and.reduce([v != NODATA for v in q])
    out = (q[0] * (1 - tx) * (1 - ty) + q[1] * tx * (1 - ty) + q[2] * (1 - tx) * ty + q[3] * tx * ty)
    return np.where(ok, out, NODATA)


def resample(source: str, outline: Path, data: Path, epsg: int = 31983, h: float = 30.0) -> None:
    win, out = _names(source, outline, data, epsg, h)
    src, lonf, latf, d, _ = read_geotiff(win)
    src = src.astype(np.float64)
    x0, y0, x1, y1 = _target_box(outline, epsg, h)
    xs = np.arange(x0, x1 + h / 2, h)
    ys = np.arange(y1, y0 - h / 2, -h)  # north to south: row 0 is the top
    inv = pyproj.Transformer.from_crs(f"EPSG:{epsg}", f"EPSG:{SOURCE_EPSG[source]}", always_xy=True)
    t0 = time.perf_counter()
    grid = np.empty((ys.size, xs.size), dtype=np.float32)
    step = 256

    def rows(i: int) -> None:
        gx, gy = np.meshgrid(xs, ys[i : i + step])
        lon, lat = inv.transform(gx, gy)
        grid[i : i + step] = bilinear(src, (lon - lonf) / d, (latf - lat) / d)

    with ThreadPoolExecutor(8) as pool:  # pyproj and numpy release the GIL in the bulk work
        list(pool.map(rows, range(0, ys.size, step)))
    secs = time.perf_counter() - t0
    out.parent.mkdir(parents=True, exist_ok=True)
    write_geotiff(out, grid, float(xs[0]), float(ys[0]), h, h, epsg, geographic=False)
    meta = {"source_window": win.name, "epsg": epsg, "spacing_m": h, "shape": list(grid.shape),
            "first_node": [float(xs[0]), float(ys[0])], "nodata_nodes": int((grid == NODATA).sum()),
            "resample_s": round(secs, 2)}  # fmt: skip
    out.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta))
    print(out)


def check(source: str, outline: Path, data: Path, epsg: int = 31983, h: float = 30.0) -> None:
    """B6 on this piece: the resampled grid, read bilinearly at each source node
    inside the domain, against the source value there."""
    win, out = _names(source, outline, data, epsg, h)
    src, lonf, latf, d, _ = read_geotiff(win)
    grid, gx0, gy0, _, _ = read_geotiff(out)
    poly, crs = _outline(outline)
    fwd = pyproj.Transformer.from_crs(f"EPSG:{SOURCE_EPSG[source]}", f"EPSG:{epsg}", always_xy=True)
    to_ll = pyproj.Transformer.from_crs(crs, f"EPSG:{SOURCE_EPSG[source]}", always_xy=True)
    poly_ll = shapely.ops.transform(to_ll.transform, poly)
    errs = []
    for i in range(0, src.shape[0], 512):
        rr, cc = np.mgrid[i : min(i + 512, src.shape[0]), 0 : src.shape[1]]
        lon, lat = lonf + cc * d, latf - rr * d
        inside = shapely.contains_xy(poly_ll, lon, lat) & (src[rr, cc] != NODATA)
        px, py = fwd.transform(lon[inside], lat[inside])
        back = bilinear(grid.astype(np.float64), (px - gx0) / h, (gy0 - py) / h)
        ok = back != NODATA
        errs.append(np.abs(back[ok] - src[rr[inside][ok], cc[inside][ok]]))
    e = np.concatenate(errs)
    res = {"source_nodes_in_domain": int(e.size), "max": float(e.max()),
           "p999": float(np.quantile(e, 0.999)), "p99": float(np.quantile(e, 0.99)),
           "median": float(np.median(e)),
           "share_over": {str(t): float((e > t).mean()) for t in (1, 2, 5, 10, 20)}}  # fmt: skip
    (out.parent / (out.stem + "_check.json")).write_text(json.dumps(res, indent=1))
    print(json.dumps(res))


if __name__ == "__main__":
    cmd, source, outline, data = sys.argv[1], sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
    rest = sys.argv[5:]
    kw = {}
    if "--spacing" in rest:
        kw["h"] = float(rest[rest.index("--spacing") + 1])
    if "--epsg" in rest:
        kw["epsg"] = int(rest[rest.index("--epsg") + 1])
    {"fetch": fetch, "resample": resample, "check": check}[cmd](source, outline, data, **kw)
