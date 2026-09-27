"""Increment 15's measurements on ANADEM v1, by HTTP range reads (no full download).

    python docs/increments/15-probes/anadem_probe.py headers   # B2, B3: CRS, lattice, overlaps
    python docs/increments/15-probes/anadem_probe.py seams     # B4: do overlaps agree?
    python docs/increments/15-probes/anadem_probe.py resample  # B6: resampling error
    python docs/increments/15-probes/anadem_probe.py frame     # B5, B7: sizes, edge bending

Needs network access to metadados.snirh.gov.br, plus tifffile, numpy and pyproj.
It is a measurement script, not production code, and nothing imports it.
"""

from __future__ import annotations

import io
import math
import sys
import urllib.request

import numpy as np
import pyproj
import tifffile

URL = "https://metadados.snirh.gov.br/files/anadem_v1_tiles/anadem_v1_{}.tif"
BASIN_TILES = ("22L", "23K", "23L", "23M", "24K", "24L", "24M", "25L", "25M")
NODATA = -9999.0
LCC = "+proj=lcc +lat_1=-10 +lat_2=-18.5 +lat_0=-14 +lon_0=-42 +ellps=WGS84 +units=m +no_defs"


class RangeFile(io.RawIOBase):
    """A seekable read-only file over HTTP range requests, cached in 64 KiB blocks."""

    BLOCK = 1 << 16

    def __init__(self, url: str) -> None:
        self.url, self.pos, self.cache, self.fetched = url, 0, {}, 0
        head = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=30)
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

    def _block(self, b: int) -> bytes:
        if b not in self.cache:
            start, end = b * self.BLOCK, min((b + 1) * self.BLOCK, self.size) - 1
            req = urllib.request.Request(self.url, headers={"Range": f"bytes={start}-{end}"})
            self.cache[b] = urllib.request.urlopen(req, timeout=60).read()
            self.fetched += end - start + 1
        return self.cache[b]

    def readinto(self, buf) -> int:  # type: ignore[no-untyped-def]
        n, out = min(len(buf), self.size - self.pos), bytearray()
        while len(out) < n:
            b, o = divmod(self.pos + len(out), self.BLOCK)
            out += self._block(b)[o : o + n - len(out)]
        buf[:n] = out
        self.pos += n
        return n


def open_tile(name: str) -> tuple[io.BufferedReader, tifffile.TiffFile]:
    f = io.BufferedReader(RangeFile(URL.format(name)), buffer_size=1 << 16)
    return f, tifffile.TiffFile(f)


def header(name: str) -> tuple[float, float, int, int, float]:
    """(lon, lat) of the upper-left corner, rows, cols, spacing."""
    _, tif = open_tile(name)
    page = tif.pages.first
    tie, scale = page.tags[33922].value, page.tags[33550].value
    return tie[3], tie[4], page.imagelength, page.imagewidth, scale[0]


def block(name: str, block_row: int, block_col: int) -> np.ndarray:
    """One 512 x 512 TIFF block, decoded, as float64."""
    f, tif = open_tile(name)
    page = tif.pages.first
    index = block_row * -(-page.imagewidth // 512) + block_col
    f.seek(page.dataoffsets[index])
    data = f.read(page.databytecounts[index])
    return np.asarray(page.decode(data, index)[0], dtype=np.float64).reshape(512, 512)


def headers() -> None:
    rows = {}
    for name in BASIN_TILES:
        f, tif = open_tile(name)
        page, keys = tif.pages.first, tif.geotiff_metadata or {}
        rows[name] = header(name)
        print(
            name, page.shape, page.dtype, page.compression.name, "tiled", page.tilewidth,
            "pages", len(tif.pages), "reduced", [p.is_reduced for p in tif.pages[1:]],
            "nodata", page.tags[42113].value if 42113 in page.tags else None,
            "model", keys.get("GTModelTypeGeoKey"), "raster", keys.get("GTRasterTypeGeoKey"),
            "geographic", keys.get("GeographicTypeGeoKey"), "bytes read", f.raw.fetched,
        )
    x0, y0, _, _, d = rows["23L"]
    print(f"spacing {d!r} deg = {d * 3600:.6f} arc-seconds")
    boxes = {}
    for name, (x, y, r, c, dd) in rows.items():
        kx, ky = (x - x0) / d, (y0 - y) / d
        boxes[name] = (round(ky), round(ky) + r, round(kx), round(kx) + c)
        print(
            f"{name}: same spacing {dd == d}; offset cols {kx:.6f} rows {ky:.6f} "
            f"(off-integer {abs(kx - round(kx)):.1e}, {abs(ky - round(ky)):.1e}); "
            f"lon [{x:.4f}, {x + c * d:.4f}] lat [{y - r * d:.4f}, {y:.4f}]; "
            f"{r * c / 1e6:.0f} M cells, {r * c * 4 / 2**30:.2f} GiB float32"
        )
    for a, b in (("23L", "24L"), ("23L", "23K"), ("23L", "23M"), ("24L", "24K"), ("24L", "24M"),
                 ("23K", "24K"), ("22L", "23L"), ("24L", "25L"), ("23M", "24M")):
        ra, rb = boxes[a], boxes[b]
        print(f"{a}/{b}: overlap rows {min(ra[1], rb[1]) - max(ra[0], rb[0])}, "
              f"cols {min(ra[3], rb[3]) - max(ra[2], rb[2])}")


def compare(label: str, a: np.ndarray, b: np.ndarray) -> None:
    va, vb = a != NODATA, b != NODATA
    both = va & vb
    diff = np.abs(a[both] - b[both])
    worst = diff.max() if diff.size else None
    print(f"{label}: {a.size} cells, both valid {int(both.sum())}, "
          f"one valid {int((va ^ vb).sum())}, equal {int((a[both] == b[both]).sum())}, "
          f"max |diff| {worst}")


def seams() -> None:
    # 23L | 24L at 42 W: 24L starts 22264 columns east of 23L, same rows. Row block 29.
    a, b = block("23L", 29, 43), block("24L", 29, 0)
    c0 = 22264 - 43 * 512
    compare("23L|24L aligned          ", a[:, c0 : c0 + 8], b[:, 0:8])
    compare("23L|24L shifted 1 column ", a[:, c0 + 1 : c0 + 8], b[:, 0:7])
    # 23L | 23K: 23K starts 30163 rows south and 1 column west of 23L.
    tr = 30163 // 512
    r0 = 30163 - tr * 512
    a, b = block("23L", tr, 20), block("23K", 0, 20)
    compare("23L|23K aligned          ", a[r0 : r0 + 8, 0:511], b[0:8, 1:512])
    compare("23L|23K shifted 1 row    ", a[r0 + 1 : r0 + 8, 0:511], b[0:7, 1:512])
    # 24L | 24K: 87-row overlap, same columns.
    a = np.vstack([block("24L", tr, 5), block("24L", tr + 1, 5)])[r0 : r0 + 87]
    compare("24L|24K 87 rows          ", a, block("24K", 0, 5)[0:87])


def bilinear(z: np.ndarray, fx: np.ndarray, fy: np.ndarray) -> np.ndarray:
    c = np.clip(np.floor(fx).astype(int), 0, z.shape[1] - 2)
    r = np.clip(np.floor(fy).astype(int), 0, z.shape[0] - 2)
    tx, ty = fx - c, fy - r
    return (z[r, c] * (1 - tx) * (1 - ty) + z[r, c + 1] * tx * (1 - ty)
            + z[r + 1, c] * (1 - tx) * ty + z[r + 1, c + 1] * tx * ty)


def resample() -> None:
    """3 x 3 blocks of 23K over the Serra do Espinhaco, resampled onto square LCC grids."""
    r_block, c_block = 18, 30
    s = np.vstack([np.hstack([block("23K", r_block + i, c_block + j) for j in range(3)])
                   for i in range(3)])
    x, y, _, _, d = header("23K")
    lon0 = x + d / 2 + c_block * 512 * d  # area-registered: nodes at cell centres
    lat0 = y - d / 2 - r_block * 512 * d
    print(f"window lon [{lon0:.4f}, {lon0 + 1535 * d:.4f}] "
          f"lat [{lat0 - 1535 * d:.4f}, {lat0:.4f}], "
          f"NoData cells {int((s == NODATA).sum())}, z [{s.min():.1f}, {s.max():.1f}]")
    fwd = pyproj.Transformer.from_crs("EPSG:4326", LCC, always_xy=True)
    inv = pyproj.Transformer.from_crs(LCC, "EPSG:4326", always_xy=True)
    rr, cc = np.mgrid[0:1536, 0:1536]
    lon, lat = lon0 + cc * d, lat0 - rr * d
    px, py = fwd.transform(lon, lat)
    m = 40
    inner = (slice(2 * m, -2 * m), slice(2 * m, -2 * m))
    for h in (30.0, 20.0, 10.0):
        xs = np.arange(px[m:-m, m:-m].min(), px[m:-m, m:-m].max(), h)
        ys = np.arange(py[m:-m, m:-m].max(), py[m:-m, m:-m].min(), -h)
        gx, gy = np.meshgrid(xs, ys)
        glon, glat = inv.transform(gx, gy)
        resampled = bilinear(s, (glon - lon0) / d, (lat0 - glat) / d)
        back = bilinear(resampled, (px[inner] - xs[0]) / h, (ys[0] - py[inner]) / h)
        e = np.abs(back - s[inner])
        print(f"h = {h:4.0f} m: {resampled.size / 1e6:.2f} M target nodes; |resampled - source| at "
              f"source nodes: max {e.max():.2f} m, p99.9 {np.quantile(e, 0.999):.2f}, "
              f"p99 {np.quantile(e, 0.99):.2f}, median {np.median(e):.3f}")
    ident = bilinear(s, (lon[100:-100, 100:-100] - lon0) / d, (lat0 - lat[100:-100, 100:-100]) / d)
    control = np.abs(ident - s[100:-100, 100:-100]).max()
    print(f"control, resampled onto its own nodes: max {control:.1e}")


def frame() -> None:
    d = 0.00026949458523585647
    geod = pyproj.Geod(ellps="WGS84")
    for lat in (-7, -14, -21):
        ew = geod.inv(-42, lat, -42 + d, lat)[2]
        ns = geod.inv(-42, lat, -42, lat + d)[2]
        print(f"lat {lat}: cell E-W {ew:.3f} m, N-S {ns:.3f} m, ratio {ns / ew:.4f}")
    cols, rows = 12 / d, 14 / d
    print(f"basin bbox 48-36 W, 21-7 S: {cols:.0f} x {rows:.0f} = {cols * rows / 1e9:.2f} G nodes, "
          f"{cols * rows * 4 / 2**30:.1f} GiB float32")
    cell = geod.inv(-42, -14, -42 + d, -14)[2] * geod.inv(-42, -14, -42, -14 + d)[2]
    print(f"basin 636 920 km2 at the 14 S cell area: {636920e6 / cell / 1e6:.0f} M nodes")
    print(f"frame spacing d * pi * a / 180 = {d * math.pi * 6378137.0 / 180!r}")
    fwd = pyproj.Transformer.from_crs("EPSG:4326", LCC, always_xy=True)
    for lat in (-7, -14, -21):
        for km in (1, 5, 20, 50):
            worst = 0.0
            for deg in range(0, 180, 15):
                a, half = math.radians(deg), 0.5 * km * 1000 / 111320
                lo = [-42 - half * math.cos(a), -42 + half * math.cos(a), -42.0]
                la = [lat - half * math.sin(a), lat + half * math.sin(a), float(lat)]
                (x0, x1, xm), (y0, y1, ym) = fwd.transform(lo, la)
                worst = max(worst, math.hypot(xm - (x0 + x1) / 2, ym - (y0 + y1) / 2))
            print(f"lat {lat}, edge {km:2d} km: straight lattice edge vs straight LCC edge, "
                  f"midpoints {worst:.4f} m apart")


if __name__ == "__main__":
    {"headers": headers, "seams": seams, "resample": resample, "frame": frame}[sys.argv[1]]()
