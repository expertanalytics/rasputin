"""Time a prototype of increment 107's overlay on the quick check's two land-cover
catchments, and count the leftover pieces whose neighbour border is found
exactly and within a band.

Run (one run each, no mesh build):
  .venv/bin/python -I docs/increments/107-probes/prefer_timing.py numedalslagen RASPUTIN_DATA TMP
  .venv/bin/python -I docs/increments/107-probes/prefer_timing.py lagan RASPUTIN_DATA TMP
TMP receives NVE's lakes for Numedalslågen's box as GeoJSON (fetched once, paged).
The overlay is the design's (section 4): only lakes meeting the domain, each
base polygon cut by the union of the lakes it meets, unclipped; leftover pieces
meeting the domain; their neighbour found among the kept polygons subdivided
into pieces of at most 64 coordinates, border lengths summed per class.
"""

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import shapely

from tin_engine.feature_input import read_class_polygons, read_polygons
from tin_engine.legends import load
from tin_engine.polygon_fractions import subdivide

NVE_LAKES = "https://kart.nve.no/enterprise/rest/services/Innsjodatabase2/MapServer/5/query"
REPLACED = 512
BANDS = (0.0, 0.001, 0.01, 1.0)


def say(text: str) -> None:
    print(text, flush=True)


def nve_lakes(bounds: tuple[float, ...], out: Path) -> Path:
    if out.exists():
        return out
    features, offset = [], 0
    while True:
        params = {
            "geometry": ",".join(map(str, bounds)),
            "geometryType": "esriGeometryEnvelope",
            "inSR": "25833",
            "outSR": "25833",
            "spatialRel": "esriSpatialRelIntersects",
            "where": "1=1",
            "outFields": "vatnlnr",
            "resultOffset": str(offset),
            "f": "geojson",
        }
        url = f"{NVE_LAKES}?{urllib.parse.urlencode(params)}"
        with urllib.request.urlopen(url, timeout=300) as r:
            page = json.load(r)
        features += page["features"]
        offset += len(page["features"])
        if not page["features"] or not page.get("exceededTransferLimit"):
            break
    crs = {"type": "name", "properties": {"name": "EPSG:25833"}}
    out.write_text(json.dumps({"type": "FeatureCollection", "crs": crs, "features": features}))
    return out


def case(name: str, data: Path, tmp: Path) -> tuple[Path, str | None, str, Path, Path]:
    if name == "numedalslagen":
        domain = data / "numedalslagen_outline_nve.geojson"
        outline = shapely.from_geojson(domain.read_text())
        lakes = nve_lakes(outline.bounds, tmp / "numedalslagen_nve_lakes.geojson")
        return data / "corine2018_dtm10_utm33.gpkg", "corine2018", "EPSG:25833", domain, lakes
    clc = data / "corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg"
    domain = data / "sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson"
    lakes = data / "sweden_smhi_svar/wfs-20261021/vattenytor_2016_3006.geojson"
    return clc, "U2018_CLC2018_V2020_20u1", "EPSG:3006", domain, lakes


def neighbours(left: np.ndarray, kept: np.ndarray, codes: np.ndarray, band: float) -> int:
    """How many leftover pieces share a border of positive length with a kept
    polygon. Kept polygons near a piece are subdivided as `assign` does; the
    border length is summed per class over every (piece, subdivided piece) pair."""
    rings = shapely.boundary(left)
    if band > 0.0:
        rings = shapely.buffer(rings, band)
    near = np.unique(shapely.STRtree(kept).query(rings)[1])
    parts, part_codes = [], []
    for k in near:
        cut = subdivide(kept[k])
        parts += cut
        part_codes += [codes[k]] * len(cut)
    tree = np.array(parts, dtype=object)
    i, j = shapely.STRtree(tree).query(rings)
    lengths = shapely.length(shapely.intersection(rings[i], shapely.boundary(tree[j])))
    pair_code = np.array(part_codes)[j]
    keys, inverse = np.unique(np.column_stack([i, pair_code]), axis=0, return_inverse=True)
    summed = np.bincount(inverse.reshape(-1), weights=lengths)
    found = np.unique(keys[summed > 0.0, 0])
    return len(found)


def main(name: str, data: Path, tmp: Path) -> None:
    clc, layer, crs, domain, lake_file = case(name, data, tmp)
    outline = shapely.union_all(shapely.get_parts(shapely.from_geojson(domain.read_text())))
    xy = np.vstack([np.asarray(p.exterior.coords) for p in shapely.get_parts(outline)])
    t = time.perf_counter()
    base, _ = read_class_polygons(clc, layer, load("corine"), crs, xy)
    say(f"{name}: base {len(base)} polygons read, {time.perf_counter() - t:.2f} s")
    t = time.perf_counter()
    lakes, _ = read_polygons(lake_file, None, crs, xy)
    say(f"  preferred {len(lakes)} lakes read, {time.perf_counter() - t:.2f} s")
    t = time.perf_counter()
    geoms = np.array([g for _, g in base], dtype=object)
    codes = np.array([c for c, _ in base])
    pref = np.array(lakes, dtype=object)
    pref = pref[shapely.STRtree(pref).query(outline, predicate="intersects")]
    b, p = shapely.STRtree(pref).query(geoms, predicate="intersects")
    cut = geoms.copy()
    for i in np.unique(b):
        cut[i] = shapely.difference(geoms[i], shapely.union_all(pref[p[b == i]]))
    keep = (codes != REPLACED) & ~shapely.is_empty(cut)
    left = shapely.get_parts(cut[codes == REPLACED])
    left = left[shapely.intersects(left, outline)]
    say(f"  cut {len(np.unique(b))} base polygons by {len(pref)} lakes: "
        f"{len(left)} leftover pieces, {time.perf_counter() - t:.2f} s")
    for band in BANDS:
        t = time.perf_counter()
        n = neighbours(left, cut[keep], codes[keep], band)
        say(f"  band {band} m: {n} of {len(left)} pieces found a neighbour, "
            f"{time.perf_counter() - t:.2f} s")
    t = time.perf_counter()
    before = sum(len(subdivide(g)) for g in geoms)
    after = sum(len(subdivide(g)) for g in cut[keep]) + sum(len(subdivide(g)) for g in pref)
    after += sum(len(subdivide(g)) for g in left)
    say(f"  pieces for assign (64 coordinates): {before} before, {after} after, "
        f"{time.perf_counter() - t:.1f} s")


if __name__ == "__main__":
    main(sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]))
