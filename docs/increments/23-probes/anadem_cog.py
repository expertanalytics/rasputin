"""Increment 23's probe of OpenTopography's ANADEM COG as a block cache source.

    python docs/increments/23-probes/anadem_cog.py DATA

DATA is `../rasputin_data/sao_francisco_piece` (for the BHO outlines). Reads
the COG's first 8 MiB with one range request, parses every IFD from that
prefix alone (a strict reader refuses any read past it), prints the
full-resolution page's block layout, decodes one block fetched by its own
range request, and counts the blocks (and their compressed bytes) meeting the
basin's and the Velhas piece's outlines. Needs network access. A measurement
script, not production code; nothing imports it.
"""

from __future__ import annotations

import io
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import shapely
import tifffile
from shapely.geometry import shape
from shapely.prepared import prep

URL = "https://opentopography.s3.sdsc.edu/raster/ANADEM/ANADEM_be/anadem_v1_compressed_COG.tif"
TIE = (-82.51654705336689, 14.079475111062084)  # 15c-geographic-dem.md, "Measured by @architect"
STEP = 0.00026949458523585647
PREFIX = 8 * 2**20


class Strict(io.BytesIO):
    """A prefix that refuses to be read past its end, so a parse that succeeds used it alone."""

    def read(self, size: int | None = -1) -> bytes:
        if size is not None and size > 0 and self.tell() + size > len(self.getbuffer()):
            raise EOFError(f"read past the {len(self.getbuffer())}-byte prefix")
        return super().read(size)


def ranged(start: int, length: int) -> bytes:
    req = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{start + length - 1}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        assert r.status == 206, r.status
        return bytes(r.read())


def main(data: Path) -> None:
    with urllib.request.urlopen(urllib.request.Request(URL, method="HEAD"), timeout=60) as r:
        keys = ("Content-Length", "ETag", "Last-Modified", "Accept-Ranges")
        print({k: r.headers[k] for k in keys})
    page = tifffile.TiffFile(Strict(ranged(0, PREFIX))).pages[0]
    offs = np.asarray(page.dataoffsets)
    counts = np.asarray(page.databytecounts)
    tile = page.tilewidth
    print("shape", page.shape, "tile", tile, "compression", page.compression, "dtype", page.dtype)
    print("blocks", len(offs), "first byte", int(offs.min()), "median", int(np.median(counts)))
    per_row = -(-page.shape[1] // tile)
    col, row = int((-43.9 - TIE[0]) / STEP), int((TIE[1] + 19.9) / STEP)
    k = (row // tile) * per_row + col // tile
    block, _, _ = page.decode(ranged(int(offs[k]), int(counts[k])), int(k))
    print("decoded block", k, block.shape, "heights", float(block.min()), "-", float(block.max()))
    for name, file in (
        ("basin", "bho2017_level2_76_raw.geojson"),
        ("velhas", "bho2017_5k_76949_outline_epsg4674.geojson"),
    ):
        outline = shape(json.loads((data / file).read_text())["features"][0]["geometry"])
        grown = prep(outline.buffer(3 * STEP))
        x0, y0, x1, y1 = outline.buffer(3 * STEP).bounds
        hits = [
            r * per_row + c
            for r in range(int((TIE[1] - y1) / STEP) // tile, int((TIE[1] - y0) / STEP) // tile + 1)
            for c in range(int((x0 - TIE[0]) / STEP) // tile, int((x1 - TIE[0]) / STEP) // tile + 1)
            if grown.intersects(
                shapely.box(
                    TIE[0] + c * tile * STEP,
                    TIE[1] - (r + 1) * tile * STEP,
                    TIE[0] + (c + 1) * tile * STEP,
                    TIE[1] - r * tile * STEP,
                )
            )
        ]
        print(f"{name}: {len(hits)} blocks, {counts[hits].sum() / 2**30:.2f} GiB compressed")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
