"""Emulate increment 110 on Sula at 20 m without changing code.

The fix puts the preferred lakes into the land cover's --features polygons
(the same `preferred.overlay` the fractions use) before the clean-up, and the
preferred source is no --features source of its own. Emulated here by writing
CORINE overlaid with NVE's lakes to a GeoJSON in EPSG:25833 and meshing with
it as the only --features source. Prints, per run: triangles, triangles under
300 m2 in the shore band and elsewhere, area and lake share.

The band: points within 40 m of both a CORINE polygon border and an NVE
shoreline inside the domain (as `docs/research/sula-direct-lake-shore-density.md`).

Run: .venv/bin/python -I docs/increments/110-probes/direct_prefer_probe.py <work dir>
(needs the clc2018 and nve-lakes caches and the Sula files under
~/rasputin-output/mcp/ola/work/).
"""

import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import numpy as np
import shapely
from tin_engine.meshing_rules import load

from tin_engine.feature_input import read_class_polygons, read_polygons
from tin_engine.preferred import overlay

HOME = Path.home()
WORK = HOME / "rasputin-output/mcp/ola/work"
DOMAIN = WORK / "sula-73-27-0/sula-official-catchment.geojson"
CLC = HOME / ".cache/rasputin/clc2018/features.gpkg"
NVE = HOME / ".cache/rasputin/nve-lakes/features.gpkg"
CLC_RULES = WORK / "uploads/sula-corine-constraints.toml"
NVE_RULES = WORK / "uploads/sula-nve-lake-constraints.toml"
CRS = "EPSG:25833"
SMALL = 300.0
BAND = 40.0
RASPUTIN = str(Path(sys.executable).parent / "rasputin")


def common(out: Path) -> list[str]:
    return [
        RASPUTIN, "mesh", "--out", str(out), "--dem", "dtm10-norway", "--tolerance", "20",
        "--domain", str(DOMAIN), "--out-crs", CRS, "--features-repair", "0.1",
        "--features-tolerance", "20", "--features-outline-snap", "0",
        "--land-cover", "clc2018", "--prefer", "lake=nve-lakes", "--threads", "4",
        "--memory-budget", "4G", "--z-unit", "m",
    ]  # fmt: skip


def triangles(path: Path) -> tuple[np.ndarray, float, float]:
    """Triangle polygons, the area and the lake share of the cells."""
    con = sqlite3.connect(path)
    blobs = [b for (b,) in con.execute("select geom from triangles")]
    cells = con.execute("select area_m2, lake from cells").fetchall()
    con.close()
    geoms = []
    for b in blobs:
        envelope = (0, 32, 48, 48, 64)[(b[3] >> 1) & 7]
        geoms.append(shapely.from_wkb(bytes(b[8 + envelope :])))
    area = sum(a for a, _ in cells)
    lake = sum(a * f for a, f in cells) / area
    return np.array(geoms, dtype=object), area, lake


def main(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    doc = json.loads(DOMAIN.read_text())
    domain = shapely.geometry.shape(doc["features"][0]["geometry"])
    xy = shapely.get_coordinates(domain.exterior)
    found, _ = read_class_polygons(CLC, None, load("corine"), CRS, xy)
    lakes, _ = read_polygons(NVE, None, CRS, xy)
    x0, y0, x1, y1 = domain.buffer(500).bounds
    cover = np.array([[[x0, y0], [x1, y0], [x1, y1]], [[x0, y0], [x1, y1], [x0, y1]]])
    over = overlay(found, [(512, lakes)], cover)
    features = [
        {"type": "Feature", "properties": {"Code_18": str(c)}, "geometry": g.__geo_interface__}
        for c, g in over.polygons
    ]
    merged = work / "clc-with-nve-lakes.geojson"
    named = {"type": "name", "properties": {"name": CRS}}
    merged.write_text(json.dumps({"type": "FeatureCollection", "crs": named, "features": features}))
    clc_lines = shapely.intersection(
        shapely.union_all([g.boundary for _, g in found]), domain
    )
    nve_lines = shapely.intersection(shapely.union_all([g.boundary for g in lakes]), domain)
    band = shapely.intersection(clc_lines.buffer(BAND), nve_lines.buffer(BAND))
    print(f"band {band.area / 1e6:.3f} km2")
    runs = {
        "direct, as today": [
            "--features", str(CLC), "--features", str(NVE),
            "--features-rules", str(CLC_RULES), "--features-rules", str(NVE_RULES),
        ],
        "direct, fix emulated": [
            "--features", str(merged), "--features-rules", str(CLC_RULES),
        ],
    }  # fmt: skip
    shapely.prepare(band)
    for name, extra in runs.items():
        out = work / (name.split(",")[-1].strip().replace(" ", "-") + ".gpkg")
        subprocess.run([*common(out), *extra], check=True, capture_output=True)
        tri, area, lake = triangles(out)
        small = shapely.area(tri) < SMALL
        inside = shapely.contains_xy(band, *shapely.get_coordinates(shapely.centroid(tri)).T)
        print(
            f"{name}: {len(tri)} triangles, small in band {int((small & inside).sum())}, "
            f"small elsewhere {int((small & ~inside).sum())}, area {area / 1e6:.4f} km2, "
            f"lake {lake:.3%}"
        )


if __name__ == "__main__":
    main(Path(sys.argv[1]))
