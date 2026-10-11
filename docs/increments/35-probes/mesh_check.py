"""The design's four changes, applied as monkeypatches to the installed
`tin_engine`, against master in the same process: `rasputin mesh` run twice
per case (the quick check's arguments plus `--out <scratch>.vtk --binary`),
the `.vtk`'s SHA-256 and the `features:` stderr line compared, the
`features clip` and `clean-up` times from `--stats` printed. Run as

    .venv/bin/python docs/increments/35-probes/mesh_check.py numedalslagen|lagan [scratch dir]

The patches are the probe's own copies of what @developer will write (section
3 of docs/increments/35-cleanup-speed.md): threaded pyproj transform (3.1),
`_untouched` by DE-9IM pattern (3.2), no `covers` where `_untouched` held
(3.3), one outline STRtree per clean-up (3.4)."""

from __future__ import annotations

import hashlib
import io
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stderr
from pathlib import Path

import numpy as np
import shapely
from pyproj import Transformer
from shapely.geometry import LineString

import tin_engine.crs as crs
import tin_engine.feature_input as fi
from tin_engine.features import TerrainFeature

DATA = "/Users/skavhaug/projects/rasputin_data"
ARGS = {
    "numedalslagen": [
        "mesh", "--dem", f"{DATA}/DTM10_UTM33_20260925",
        "--domain", f"{DATA}/numedalslagen_outline_nve.geojson",
        "--features", f"{DATA}/corine2018_dtm10_utm33.gpkg", "--features-layer", "corine2018",
        "--features-map", "corine", "--tolerance", "10",
    ],
    "lagan": [
        "mesh", "--dem", "glo30", "--cache", f"{DATA}/sweden_glo30_cache", "--out-crs", "EPSG:3006",
        "--domain", f"{DATA}/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson",
        "--features",
        f"{DATA}/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg",
        "--features-layer", "U2018_CLC2018_V2020_20u1", "--features-crs", "EPSG:3035",
        "--features-map", "corine", "--tolerance", "10",
    ],
}

# --- 3.1: the transform of a large array in chunks on threads, one Transformer each ---
POINTS_PER_CHUNK, WORKERS = 250_000, max(1, min(8, os.cpu_count() or 1))


def reprojector(src, dst):
    source, target = crs.parse_crs(src), crs.parse_crs(dst)
    local = threading.local()

    def transformer() -> Transformer:
        if not hasattr(local, "t"):
            local.t = Transformer.from_crs(source, target, always_xy=True)
        return local.t

    def part(points):
        x, y = transformer().transform(points[:, 0], points[:, 1])
        return np.column_stack([np.asarray(x, np.float64), np.asarray(y, np.float64)])

    def apply(xy):
        points = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        chunks = min(WORKERS, len(points) // POINTS_PER_CHUNK)
        if chunks < 2:
            return part(points)
        with ThreadPoolExecutor(chunks) as pool:
            return np.concatenate(list(pool.map(part, np.array_split(points, chunks))))

    return apply


# --- 3.2: a line whose interior lies in the domain's interior, ends in its closure ---
def _untouched(line, domain):
    xy = shapely.get_coordinates(line)
    return (
        len(xy) >= 2
        and not (xy[1:] == xy[:-1]).all(axis=1).any()
        and bool(shapely.relate_pattern(domain, line, "T**F**FF*"))
        and bool(shapely.is_simple(line))
    )


# --- 3.3: covers only where _untouched did not hold ---
def _add(self, fid, mask, kept, code, polygon):
    coded = {"code": code, "polygon": polygon}
    dom = self.domain.polygon
    whole = [fi._untouched(line, dom) for line in kept]
    lines = tuple(
        piece
        for line, w in zip(kept, whole, strict=True)
        for piece in ((line,) if w else shapely.get_parts(shapely.intersection(line, dom)))
        if isinstance(piece, LineString) and piece.length > 0
    )
    if lines:
        self.features.append(TerrainFeature(fid=fid, mask=mask, lines=lines, **coded))
        self.clipped += not all(w or dom.covers(g) for g, w in zip(kept, whole, strict=True))
    elif polygon is not None and polygon.intersects(self.inner):
        self.features.append(TerrainFeature(fid=fid, mask=mask, lines=(), **coded))
    else:
        self.outside += 1


# --- 3.4: the outline's segment tree built once per outline ---
_trees: dict[int, tuple] = {}


def _outline_lines(polygon, outline):
    if id(outline) not in _trees:
        rims = [shapely.get_coordinates(r) for r in shapely.get_rings(outline)]
        segments = shapely.linestrings(
            np.concatenate([np.stack([r[:-1], r[1:]], axis=1) for r in rims])
        )
        _trees[id(outline)] = (segments, shapely.STRtree(segments))
    segments, tree = _trees[id(outline)]
    out = []
    for ring in shapely.get_rings(shapely.get_parts(polygon)):
        xy = shapely.get_coordinates(ring)
        edges = shapely.linestrings(np.stack([xy[:-1], xy[1:]], axis=1))
        ei, si = tree.query(edges, predicate="dwithin", distance=fi.IN_LINE)
        near = [
            shapely.distance(shapely.points(xy[ei + k]), segments[si]) <= fi.IN_LINE for k in (0, 1)
        ]
        on = np.zeros(len(xy) - 1, dtype=bool)
        on[ei[near[0] & near[1]]] = True
        out += fi._chains(xy[:-1], on.tolist())
    return tuple(out)


# --- 3.1 too: the cover polygons moved in one shapely.transform call, after the loop ---
def _take(self, source, own, rows):
    name, cmap = source.path.name, source.class_map
    src = fi.parse_crs(own)
    move = None if fi.same_crs(src, self.dem) else fi.reprojector(src, self.dem)
    bound, region = None, fi.source_region(self.domain, self.dem, self.dem)
    shapely.prepare(region)
    if move is not None and src.is_geographic:
        steps = fi.transform_definition(src, self.dem)
        if "gridshift" not in steps and "deformation" not in steps:
            bound = fi._bound(self.dem, fi.source_region(self.domain, self.dem, src))
    cover = []
    for fid, geometry, value in rows:
        t0 = time.perf_counter()
        if geometry is None or geometry.is_empty:
            self.empty += 1
            continue
        geometry = shapely.force_2d(geometry)
        if geometry.geom_type not in fi.ACCEPTED:
            raise fi.FeatureError(
                f"{name}: feature {fid} is a {geometry.geom_type}; input is polygons and lines"
            )
        mask = fi._mask(cmap, self.vocabulary, value, f"{name}: feature {fid}")
        if mask is None:
            continue
        code = fi._code(value, f"{name}: feature {fid}") if cmap.codes else None
        polygonal = code is not None and geometry.geom_type in ("Polygon", "MultiPolygon")
        moved = polygon = None
        if polygonal and code is not None and self.cleanup is not None:
            cover.append((fid, mask, code, geometry))  # moved below, all at once
            self.seconds += time.perf_counter() - t0
            continue
        if move is not None and bound is not None and bound[1](geometry):
            chains = fi.pre_clip(geometry, bound[2], bound[0])
            kept = [shapely.transform(c, move) for c in chains]
            polygon = shapely.transform(geometry, move) if polygonal else None
        else:
            moved = shapely.transform(geometry, move) if move is not None else geometry
            kept = [moved]
            polygon = moved if polygonal else None
        if not all(fi._finite(g) for g in kept):
            raise fi.FeatureError(f"{name}: feature {fid} has a vertex with no image in the DEM")
        if moved is not None:
            kept = []
            if region.intersects(fi._linework(moved)):
                kept = list(fi.pre_clip(moved, region))
        self._add(fid, mask, kept, code, polygon)
        self.seconds += time.perf_counter() - t0
    if cover and self.cleanup is not None:
        t0 = time.perf_counter()
        if move is not None:
            geoms = shapely.transform(np.array([c[3] for c in cover], dtype=object), move)
            for (fid, *_), g in zip(cover, geoms, strict=True):
                if not fi._finite(g):
                    raise fi.FeatureError(
                        f"{name}: feature {fid} has a vertex with no image in the DEM"
                    )
            cover = [(*c[:3], g) for c, g in zip(cover, geoms, strict=True)]
        self.seconds += time.perf_counter() - t0
        self._clean(cover, region, self.cleanup)


MASTER = (
    crs.reprojector, fi.reprojector, fi._untouched, fi._Tally._add, fi._outline_lines,
    fi._Tally._take,
)


def patch(on: bool) -> None:
    if on:
        crs.reprojector = fi.reprojector = reprojector
        fi._untouched, fi._Tally._add, fi._outline_lines = _untouched, _add, _outline_lines
        fi._Tally._take = _take
    else:
        (crs.reprojector, fi.reprojector, fi._untouched, fi._Tally._add, fi._outline_lines,
         fi._Tally._take) = MASTER


def run(case: str, scratch: Path, label: str) -> None:
    from tin_engine.cli import app

    out, stats = scratch / f"{case}-{label}.vtk", scratch / f"{case}-{label}-stats.md"
    sys.argv = ["rasputin", *ARGS[case], "--out", str(out), "--binary", "--stats", str(stats)]
    err = io.StringIO()
    t = time.perf_counter()
    with redirect_stderr(err):
        try:
            app()
        except SystemExit as exc:
            if exc.code not in (0, None):
                print(err.getvalue())
                raise
    took = time.perf_counter() - t
    digest = hashlib.sha256(out.read_bytes()).hexdigest()
    features = next(line for line in err.getvalue().splitlines() if line.startswith("features:"))
    pattern = r"\|\s*(features clip(?:: clean-up)?)\s*\|\s*([\d.]+)"
    rows = {m[1]: float(m[2]) for m in re.finditer(pattern, stats.read_text())}
    print(f"{case} {label}: {took:.2f} s wall, {rows}, {features}\n  sha256 {digest}")


def main(case: str, scratch: str) -> None:
    sd = Path(scratch)
    sd.mkdir(parents=True, exist_ok=True)
    patch(False)
    run(case, sd, "master")
    patch(True)
    run(case, sd, "patched")
    a, b = (
        hashlib.sha256((sd / f"{case}-{k}.vtk").read_bytes()).hexdigest()
        for k in ("master", "patched")
    )
    print(f"{case}: meshes byte-identical: {a == b}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "/tmp/35-probes")
