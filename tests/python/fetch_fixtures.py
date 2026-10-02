"""A local HTTP range server and micro COGs for 23a-2's suites (F1).

`docs/increments/23-basin-scale.md`, "The fetch step and the tile cache" and
its test list F1-F11. No test touches the network: every URL is
`http://127.0.0.1:<port>/...`, served by `RangeServer` from bytes held in
memory, and a catalogue entry naming it is built here and passed to `fetch`
(or put into `SOURCES` by `monkeypatch` for the CLI).

The server is `http.server` in a thread with its own handler. It answers a
single `Range: bytes=a-b` with 206 and `Content-Range: bytes a-b/total`,
clipping `b` to the file as RFC 9110 §14.1.2 and S3 do, logs every request's
range, and can be told to misbehave: fail with 500 after k responses, answer
every request with a fixed status (500, 404), answer 200 ignoring `Range`,
send a short body, send a wrong `Content-Range`, or report another total or
`Last-Modified` on block responses (any range not starting at 0) while the
prefix stays as it was. Its own suite (`test_fetch_server.py`) shows each
fault reaches a plain `http.client` caller, so a refusal the fetch step must
make is one the fixture can provoke.

Files. Every served file is padded past its last byte with zeros (legal
after a TIFF's last referenced byte) to at least `PADDED` bytes, so the
fetch step's first 1 MiB prefix request is satisfiable in full and a short
prefix is never the file's own end. The projected source is 23a-1's tiled
Deflate fixture (`cog_fixtures.build`, EPSG:25833, 16 x 16 tiles, one
overview); `wide` is an uncompressed one whose block rows lie more than
64 KiB apart, so coalescing has gaps to respect; `geographic` is EPSG:4326;
`glo30_tile` is one GLO-30-style 1° tile named as the bucket names them.

PINNED BY THESE SUITES (the record names the parts, not their spelling):
- `tin_engine.fetch.http`: `RangeClient(*, delays=(1, 2, 4), timeout=...)`,
  `.get(url, start, stop) -> RangeResponse(data, total, last_modified)`,
  `.get_text(url) -> str`; `FetchError(Exception)`.
- `tin_engine.fetch.plan`: `FetchRequest(source, domain=None, box=None,
  out_crs=None, margin=4, connections=8, dry_run=False, refresh=False)`,
  exactly one of `domain` and `box` (else `ValueError`);
  `source_box(request, meta) -> Bounds` in the source CRS (`meta.epsg`);
  `plan_object(object_id, url, page, meta, box, present) -> ObjectPlan`,
  pure, `present` the block indices already cached, a box meeting no block a
  `FetchError`; `ObjectPlan(object_id, url, blocks, ranges, bytes)` with
  `blocks` every block meeting the box, ascending, `ranges` the coalesced
  `(start, stop)` byte ranges of the missing non-sparse ones and `bytes` the
  sum of the ranges' lengths; `coalesce(spans) -> ((start, stop), ...)` over
  `(start, stop)` block spans; `parse_prefix(prefix, *, nodata) ->
  (RasterMeta, dtype, TiffPage) | None`, `None` while the prefix lacks an
  offset or byte count of any block, geographic CRSs allowed.
- `tin_engine.fetch.run`: `async fetch(request, source, client, writer) ->
  FetchReport(objects, needed, present, fetched, empty, bytes, requests,
  no_tile, seconds)`; `fetch` enters `writer` itself and a dry run never
  does. A dry run's report holds the counts the real run then reports.
- `tin_engine.io.repository`: `CacheWriter(root, source)`; entering it takes
  the lock, and a held lock is a `CacheError`. `CachedRequest(region_sha256,
  date)` (decided 9 renames `domain_sha256`).
- `tin_engine.sources`: `RemoteSource.cite: tuple[str, ...] = ()`,
  `notice(source) -> str`; GLO-30 tiles are named
  `Copernicus_DSM_COG_10_<N|S><lat:02>_00_<E|W><lon:03>_00_DEM`, one per
  line of the tile list, and `url_template` takes `{name}`.
- Refusals: a changed remote, a non-206 answer, a bad `Content-Range` or
  length, and a header CRS other than the catalogue's are `FetchError`; an
  `OSError` writing the cache propagates as itself.
"""

from __future__ import annotations

import io
import socketserver
import threading
import time
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any

import numpy as np
import tifffile

from cog_fixtures import VARIANTS, build
from geotiff_fixtures import (
    EPSG_WGS84,
    GEOGRAPHIC_TYPE,
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    PIXEL_IS_POINT,
    micro_tiff,
)

MIB = 1 << 20
GAP = 64 * 1024
MAX_RANGE = 8 * MIB
PADDED = 3 * MIB // 2
LAST_MODIFIED = "Thu, 01 Oct 2026 00:00:00 GMT"
CHANGED = "Fri, 02 Oct 2026 12:00:00 GMT"
PROJECTED_CRS = "EPSG:25833"

GEOGRAPHIC_KEYS: Mapping[int, int] = {
    GT_MODEL_TYPE: 2,  # ModelTypeGeographic
    GT_RASTER_TYPE: PIXEL_IS_POINT,
    GEOGRAPHIC_TYPE: EPSG_WGS84,
}


def padded(data: bytes, size: int = PADDED) -> bytes:
    """`data` followed by zeros up to `size` bytes (no shorter than it was)."""
    return data + bytes(max(0, size - len(data)))


# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------


def projected() -> bytes:
    """23a-1's tiled Deflate fixture: 50 x 70 nodes, 10 x 5 m, 4 x 5 blocks."""
    from test_cli_mesh_mosaic import terrain

    return padded(build(VARIANTS[0], array=terrain(50, 70)))


def wide() -> bytes:
    """48 x 1600 uncompressed float32 in 16 x 16 tiles: 3 x 100 blocks of
    1 KiB, so a block row is 100 KiB long and a box a few blocks wide has
    gaps over 64 KiB between its block rows."""
    rng = np.random.default_rng(232)
    array = rng.uniform(100, 200, (48, 1600)).astype(np.float32)
    stream = micro_tiff(array, tile=(16, 16), extra_pages=[(array[::2, ::2], 1)])
    return padded(stream.getvalue())


def with_long_header(size: int) -> bytes:
    """A tiled file whose ImageDescription (tag 270, written before the
    offsets) is `size` bytes, so the block offsets lie past `size`."""
    array = np.arange(48 * 64, dtype=np.float32).reshape(48, 64)
    stream = micro_tiff(
        array,
        tile=(16, 16),
        compression="deflate",
        extra_pages=[(array[::2, ::2], 1)],
        description="x" * size,
        metadata=None,
    )
    return padded(stream.getvalue())


def geographic(lon: float = 15.0, lat: float = 59.5, step: float = 0.001) -> bytes:
    """EPSG:4326, 40 x 48 nodes from (lon, lat) eastwards and southwards."""
    from test_cli_mesh_mosaic import terrain

    stream = micro_tiff(
        terrain(40, 48),
        tiepoint=(0.0, 0.0, 0.0, lon, lat, 0.0),
        scale=(step, step, 0.0),
        geokeys=GEOGRAPHIC_KEYS,
        tile=(16, 16),
        compression="deflate",
        extra_pages=[(terrain(20, 24), 1)],
    )
    return padded(stream.getvalue())


def glo30_name(lat: int, lon: int) -> str:
    """The bucket's name of the 1° tile whose south-west corner is (lat, lon)."""
    ns, ew = ("N" if lat >= 0 else "S"), ("E" if lon >= 0 else "W")
    return f"Copernicus_DSM_COG_10_{ns}{abs(lat):02d}_00_{ew}{abs(lon):03d}_00_DEM"


def glo30_tile(lat: int, lon: int) -> bytes:
    """A 1° tile as a 40 x 40 node COG, top-left node at (lon, lat + 1)."""
    from test_cli_mesh_mosaic import terrain

    stream = micro_tiff(
        terrain(40, 40),
        tiepoint=(0.0, 0.0, 0.0, float(lon), float(lat + 1), 0.0),
        scale=(1 / 40, 1 / 40, 0.0),
        geokeys=GEOGRAPHIC_KEYS,
        tile=(16, 16),
        compression="deflate",
        extra_pages=[(terrain(20, 20), 1)],
    )
    return padded(stream.getvalue())


def page_of(data: bytes) -> tifffile.TiffPage:
    page = tifffile.TiffFile(io.BytesIO(data)).pages.first
    assert isinstance(page, tifffile.TiffPage)
    return page


def block_span(page: tifffile.TiffPage, index: int) -> tuple[int, int]:
    offset = int(page.dataoffsets[index])
    return offset, offset + int(page.databytecounts[index])


# --------------------------------------------------------------------------
# The server
# --------------------------------------------------------------------------


@dataclass
class Logged:
    """One request: its path, the range asked (`stop` exclusive, `None` with
    no `Range`), the status answered, the body bytes the socket accepted, and
    whether the handler has finished."""

    path: str
    start: int | None
    stop: int | None
    status: int = 0
    sent: int = 0
    done: bool = False


@dataclass
class Served:
    data: bytes
    last_modified: str = LAST_MODIFIED


@dataclass
class Faults:
    """What the server does wrong; all off by default."""

    fail_after: int | None = None  # 500 for every response after this many
    status: int | None = None  # every response this status
    ignore_range: bool = False  # 200 with the whole body
    short: bool = False  # one byte fewer than Content-Length, then close
    wrong_content_range: bool = False  # the reported start one byte late
    block_total: int | None = None  # added to the total on block responses
    block_last_modified: str | None = None  # Last-Modified on block responses


class _Quick(ThreadingHTTPServer):
    """No `socket.getfqdn` in `server_bind`: a reverse lookup of 127.0.0.1
    can take seconds on macOS, and nothing here needs the server's name."""

    def server_bind(self) -> None:
        socketserver.TCPServer.server_bind(self)
        self.server_name, self.server_port = "127.0.0.1", int(self.server_address[1])


class RangeServer:
    """`http://127.0.0.1:<port>/<path>` for each of `files`, in a thread."""

    def __init__(self) -> None:
        self.files: dict[str, Served] = {}
        self.faults = Faults()
        self.log: list[Logged] = []
        self._lock = threading.Lock()
        self._answered = 0
        self._httpd = _Quick(("127.0.0.1", 0), _handler(self))
        self._httpd.daemon_threads = True
        self._thread = threading.Thread(
            target=self._httpd.serve_forever, kwargs={"poll_interval": 0.02}, daemon=True
        )

    def __enter__(self) -> RangeServer:
        self._thread.start()
        return self

    def __exit__(self, *exc: object) -> None:
        self._httpd.shutdown()
        self._httpd.server_close()

    def url(self, path: str) -> str:
        return f"http://127.0.0.1:{self._httpd.server_address[1]}/{path}"

    def put(self, path: str, data: bytes, last_modified: str = LAST_MODIFIED) -> str:
        self.files[path] = Served(data, last_modified)
        return self.url(path)

    def reset(self) -> None:
        """Faults off, log and response count cleared."""
        with self._lock:
            self.faults, self.log, self._answered = Faults(), [], 0

    def ranges(self, path: str) -> list[tuple[int, int]]:
        """Every `(start, stop)` asked of `path`, in order of arrival."""
        return [
            (e.start, e.stop)
            for e in self.log
            if e.path == path and e.start is not None and e.stop is not None
        ]

    def block_ranges(self, path: str) -> list[tuple[int, int]]:
        """The ranges asked of `path` that are not a header prefix (start 0)."""
        return [(a, b) for a, b in self.ranges(path) if a != 0]

    def wait_done(self, timeout: float = 10.0) -> None:
        """Until every logged request's handler has finished."""
        deadline = time.monotonic() + timeout
        while not all(e.done for e in self.log):
            assert time.monotonic() < deadline, "a handler did not finish"
            time.sleep(0.01)

    def _next(self, entry: Logged) -> int:
        with self._lock:
            self.log.append(entry)
            self._answered += 1
            return self._answered


def _handler(server: RangeServer) -> type[BaseHTTPRequestHandler]:
    class Handler(BaseHTTPRequestHandler):
        timeout = 10  # a write to a client that stopped reading ends

        def log_message(self, format: str, *args: Any) -> None:
            pass

        def do_GET(self) -> None:
            path = self.path.lstrip("/")
            start, stop = _parse_range(self.headers.get("Range"))
            entry = Logged(path, start, stop)
            count = server._next(entry)
            try:
                self._answer(entry, count)
            finally:
                entry.done = True

        def _answer(self, entry: Logged, count: int) -> None:
            faults = server.faults
            served = server.files.get(entry.path)
            if faults.fail_after is not None and count > faults.fail_after:
                return self._status(entry, 500)
            if faults.status is not None:
                return self._status(entry, faults.status)
            if served is None:
                return self._status(entry, 404)
            data, total = served.data, len(served.data)
            if entry.start is None or faults.ignore_range:
                entry.status = 200
                self.send_response(200)
                self.send_header("Content-Length", str(total))
                self.send_header("Last-Modified", served.last_modified)
                self.end_headers()
                return self._write(entry, data)
            if entry.start >= total:
                return self._status(entry, 416)
            last = min(entry.stop if entry.stop is not None else total, total) - 1
            body = data[entry.start : last + 1]
            block = entry.start != 0
            reported_start = entry.start + (1 if faults.wrong_content_range else 0)
            reported_total = total + ((faults.block_total or 0) if block else 0)
            modified = (block and faults.block_last_modified) or served.last_modified
            entry.status = 206
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {reported_start}-{last}/{reported_total}")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Last-Modified", modified)
            self.end_headers()
            if faults.short:
                self._write(entry, body[:-1])
                self.close_connection = True
                return None
            return self._write(entry, body)

        def _status(self, entry: Logged, status: int) -> None:
            entry.status = status
            self.send_response(status)
            self.send_header("Content-Length", "0")
            self.end_headers()

        def _write(self, entry: Logged, body: bytes) -> None:
            view = memoryview(body)
            try:
                for at in range(0, len(view), 64 * 1024):
                    self.wfile.write(view[at : at + 64 * 1024])
                    entry.sent += len(view[at : at + 64 * 1024])
                self.wfile.flush()
            except (BrokenPipeError, ConnectionResetError, TimeoutError):
                self.close_connection = True

    return Handler


def _parse_range(header: str | None) -> tuple[int | None, int | None]:
    """`bytes=a-b` as `(a, b + 1)`; `bytes=a-` as `(a, None)`; else `(None, None)`."""
    if header is None or not header.startswith("bytes="):
        return None, None
    first, _, last = header.removeprefix("bytes=").partition("-")
    if "," in last or not first:
        return None, None
    return int(first), (int(last) + 1 if last else None)


def snapshot(root: Any) -> dict[str, bytes]:
    """`{relative path: bytes}` for every file under `root` (directories as
    `<path>/` with no bytes, so an empty directory shows)."""
    from pathlib import Path

    base = Path(root)
    out: dict[str, bytes] = {}
    for path in sorted(base.rglob("*")):
        key = str(path.relative_to(base))
        out[key + "/" if path.is_dir() else key] = b"" if path.is_dir() else path.read_bytes()
    return out


def block_files(object_dir: Any) -> dict[str, bytes]:
    """`{"<row>/<col>.bin": bytes}` under `<object>/blocks`, `.part` files included."""
    from pathlib import Path

    blocks = Path(object_dir) / "blocks"
    if not blocks.is_dir():
        return {}
    return {
        str(p.relative_to(blocks)): p.read_bytes() for p in sorted(blocks.rglob("*")) if p.is_file()
    }


def serve() -> Iterator[RangeServer]:
    with RangeServer() as server:
        yield server
