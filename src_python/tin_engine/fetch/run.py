"""The fetch: headers, identity, plan, bounded async download (increment 23a-2).

`docs/increments/23-basin-scale.md`, "The cache" and "Downloading". `client`
and `writer` are parameters, so the run is tested against a local server and
a cache in `tmp_path`, and an API worker awaits it in its own loop.

Order of writes (a crash at any point leaves a cache 23a-1 reads): headers and
the manifest first, so a new object is listed before any of its blocks; then
blocks, each request's only after its whole body arrived and every length
matched (decided 4); then the request appended to the manifest. A dry run
reads headers and the cache, and writes nothing.
"""

from __future__ import annotations

import asyncio
import datetime
import hashlib
import importlib.metadata
import math
import time
from collections.abc import Callable
from typing import Any

import tifffile
from pydantic import BaseModel, ConfigDict
from pyproj import CRS

from tin_engine.crs import parse_crs
from tin_engine.fetch.http import FetchError, RangeClient
from tin_engine.fetch.plan import FetchRequest, ObjectPlan, parse_prefix, plan_object, source_box
from tin_engine.io.models import RasterMeta
from tin_engine.io.repository import CachedObject, CachedRequest, CacheManifest, CacheWriter
from tin_engine.mosaic import Bounds
from tin_engine.sources import RemoteSource, notice

MIB = 1 << 20
MAX_HEADER = 64 * MIB
#: GLO-30's tiles are chosen from a box grown as if the spacing were this many
#: degrees: more than any tile's below 80° (3", 0.00083°), so the box each tile's
#: blocks are planned from contains the one its own spacing would give.
TILE_SPACING = 0.001

Header = tuple[bytes, int, str, RasterMeta, tifffile.TiffPage]


class FetchReport(BaseModel):
    """What a fetch did, or, for a dry run, what the real run will do."""

    model_config = ConfigDict(frozen=True)

    objects: int
    needed: int
    present: int
    fetched: int
    empty: int
    bytes: int
    requests: int
    no_tile: tuple[str, ...]
    seconds: float
    plans: tuple[ObjectPlan, ...] = ()


async def fetch(
    request: FetchRequest,
    source: RemoteSource,
    client: RangeClient,
    writer: CacheWriter,
    *,
    progress: Callable[[int, int], None] | None = None,
) -> FetchReport:
    """Copy into `writer`'s cache what a mesh of `request`'s region reads.
    `progress(done, total)` is called with bytes after each request."""
    if request.dry_run:
        return await _run(request, source, client, writer, progress)
    with writer:
        if request.refresh:
            writer.discard()
        return await _run(request, source, client, writer, progress)


async def _run(
    request: FetchRequest,
    source: RemoteSource,
    client: RangeClient,
    writer: CacheWriter,
    progress: Callable[[int, int], None] | None,
) -> FetchReport:
    t0 = time.perf_counter()
    known = None if request.refresh else writer.manifest()
    entries = dict(known.objects) if known else {}
    gate = asyncio.Semaphore(request.connections)

    async def header(object_id: str, url: str) -> Header:
        async with gate:
            return await asyncio.to_thread(_header, client, source, object_id, url, entries)

    no_tile: tuple[str, ...] = ()
    if source.kind == "one-cog":
        assert source.url is not None
        object_id = source.url.rsplit("/", 1)[-1].split("?")[0].rsplit(".", 1)[0]
        headers = {object_id: (source.url, await header(object_id, source.url))}
        box = source_box(request, headers[object_id][1][3])
    else:
        box, listed, unlisted = await _tiles(request, source, client)
        no_tile = tuple(unlisted)
        assert source.url_template is not None
        urls = {n: source.url_template.format(name=n) for n in listed}
        got = await asyncio.gather(*(header(n, u) for n, u in urls.items()))
        headers = dict(zip(urls, zip(urls.values(), got, strict=True), strict=True))
    plans: list[tuple[ObjectPlan, tifffile.TiffPage, set[int]]] = []
    for object_id, (url, (_, _, _, meta, page)) in sorted(headers.items()):
        blocks = plan_object(object_id, url, page, meta, box, ()).blocks
        present = writer.present(object_id, page, blocks) if object_id in entries else set()
        plans.append((plan_object(object_id, url, page, meta, box, present), page, present))
    empty = [
        (p.object_id, page, i)
        for p, page, present in plans
        for i in p.blocks
        if i not in present and page.databytecounts[i] == 0
    ]
    needed = sum(len(p.blocks) for p, _, _ in plans)
    have = sum(len(present) for _, _, present in plans)
    report = FetchReport(
        objects=len(plans),
        needed=needed,
        present=have,
        fetched=needed - have - len(empty),
        empty=len(empty),
        bytes=sum(p.bytes for p, _, _ in plans),
        requests=sum(len(p.ranges) for p, _, _ in plans),
        no_tile=no_tile,
        seconds=0.0,
        plans=tuple(p for p, _, _ in plans),
    )
    if request.dry_run:
        return report.model_copy(update={"seconds": time.perf_counter() - t0})

    new = {i: h for i, (_, h) in headers.items() if i not in entries}
    for object_id, (prefix, total, modified, _, page) in new.items():
        writer.put_header(object_id, prefix)
        entries[object_id] = CachedObject(
            url=headers[object_id][0],
            content_length=total,
            last_modified=modified,
            header_sha256=hashlib.sha256(prefix).hexdigest(),
            header_bytes=len(prefix),
            block_shape=tuple(int(n) for n in page.chunks[:2]),
        )
    epsg = next(iter(headers.values()))[1][3].epsg
    manifest = CacheManifest(
        source=source.id,
        crs=f"EPSG:{epsg}",
        rasputin_version=importlib.metadata.version("rasputin"),
        objects=dict(sorted(entries.items())),
        requests=known.requests if known else (),
    )
    if new or known is None:
        writer.put_manifest(manifest)
    writer.put_notice(notice(source))
    for object_id, page, index in empty:
        writer.put_block(object_id, page, index, b"")

    done = 0

    def take(plan: ObjectPlan, page: tifffile.TiffPage, start: int, stop: int) -> None:
        nonlocal done
        got = client.get(plan.url, start, stop)
        entry = entries[plan.object_id]
        if (got.total, got.last_modified) != (entry.content_length, entry.last_modified):
            raise FetchError(_changed(plan.url, entry, got.total, got.last_modified))
        for i in plan.blocks:
            offset, count = int(page.dataoffsets[i]), int(page.databytecounts[i])
            if count and start <= offset and offset + count <= stop:
                writer.put_block(plan.object_id, page, i, got.data[offset - start :][:count])
        done += stop - start
        if progress is not None:
            progress(done, report.bytes)

    async def one(plan: ObjectPlan, page: tifffile.TiffPage, start: int, stop: int) -> None:
        async with gate:
            await asyncio.to_thread(take, plan, page, start, stop)

    jobs = [one(p, page, a, b) for p, page, _ in plans for a, b in p.ranges]
    failed = [r for r in await asyncio.gather(*jobs, return_exceptions=True) if r is not None]
    if failed:
        raise failed[0]
    asked = CachedRequest(region_sha256=_region(request, source), date=datetime.date.today())
    if asked not in manifest.requests:
        writer.put_manifest(manifest.model_copy(update={"requests": (*manifest.requests, asked)}))
    return report.model_copy(update={"seconds": time.perf_counter() - t0})


def _header(
    client: RangeClient,
    source: RemoteSource,
    object_id: str,
    url: str,
    entries: dict[str, CachedObject],
) -> Header:
    """The object's header prefix, parsed. A known object's is re-read at its
    recorded length and must be the recorded one; a new object's prefix starts
    at 1 MiB and doubles until the full page's block offsets are all in it."""
    entry = entries.get(object_id)
    size = entry.header_bytes if entry else MIB
    while True:
        got = client.get(url, 0, size)
        seen = (hashlib.sha256(got.data).hexdigest(), got.total, got.last_modified)
        if entry is not None and seen != (
            entry.header_sha256,
            entry.content_length,
            entry.last_modified,
        ):
            raise FetchError(_changed(url, entry, got.total, got.last_modified))
        parsed = parse_prefix(got.data, nodata=source.nodata)
        if parsed is not None:
            break
        if entry is not None or size >= MAX_HEADER or len(got.data) < size:
            raise FetchError(f"{url}: no complete header in the first {MAX_HEADER // MIB} MiB")
        size *= 2
    meta, _, page = parsed
    if CRS.from_epsg(meta.epsg) != parse_crs(source.crs):
        raise FetchError(
            f"{object_id}: the header's CRS is EPSG:{meta.epsg}; the catalogue's {source.id} "
            f"is {source.crs}"
        )
    return got.data, got.total, got.last_modified, meta, page


def _changed(url: str, entry: CachedObject, total: int, modified: str) -> str:
    return (
        f"{url}: the remote copy changed (length {total}, Last-Modified {modified!r}; the cache "
        f"has {entry.content_length}, {entry.last_modified!r}); --refresh re-fetches it"
    )


async def _tiles(
    request: FetchRequest, source: RemoteSource, client: RangeClient
) -> tuple[Bounds, list[str], list[str]]:
    """GLO-30: the box, and the 1° tiles meeting it, listed and not."""
    assert source.tile_list_url is not None
    text = await asyncio.to_thread(client.get_text, source.tile_list_url)
    listed = set(text.split())
    epsg = parse_crs(source.crs).to_epsg()
    assert epsg is not None
    nominal = RasterMeta(
        x_min=0.0, y_max=0.0, delta_x=TILE_SPACING, delta_y=TILE_SPACING, cols=1, rows=1,
        epsg=epsg, nodata=None, nodata_source="absent", pixel_is_area=False,
        vertical_unit_assumed=True,
    )  # fmt: skip
    box = source_box(request, nominal)
    names = [
        _tile_name(lat, lon)
        for lat in range(math.floor(box.y_min), math.ceil(box.y_max))
        for lon in range(math.floor(box.x_min), math.ceil(box.x_max))
    ]
    return box, [n for n in names if n in listed], sorted(n for n in names if n not in listed)


def _tile_name(lat: int, lon: int) -> str:
    ns, ew = ("N" if lat >= 0 else "S"), ("E" if lon >= 0 else "W")
    return f"Copernicus_DSM_COG_10_{ns}{abs(lat):02d}_00_{ew}{abs(lon):03d}_00_DEM"


def _region(request: FetchRequest, source: RemoteSource) -> str:
    """The sha256 of the region as given: the domain's WKB and CRS, or the
    box's four numbers and its CRS (decided 9)."""
    if request.domain is not None:
        data: Any = request.domain.polygon.wkb + request.domain.crs.encode()
    else:
        assert request.box is not None
        b = request.box
        data = repr((b.x_min, b.y_min, b.x_max, b.y_max, request.out_crs or source.crs)).encode()
    return hashlib.sha256(data).hexdigest()


__all__ = ["FetchReport", "fetch"]
