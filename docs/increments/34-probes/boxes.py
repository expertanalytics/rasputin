"""Increment 34 design probe: rasputin master's uniform meshes of the two
simulation windows, from the outline of the window (README item 3).

Usage: python boxes.py PYTHON RASPUTIN_PY PKG OUT_DIR [tolerances...]
where RASPUTIN_PY runs tin_engine.cli's app and PKG holds tin_engine with a
built _core. Writes the box domains as GeoJSON in OUT_DIR.
"""

import json
import subprocess
import sys

PY, RASPUTIN, PKG, OUT = sys.argv[1:5]
D = "/Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925"
BOXES = {
    "romsdalen": (118000, 6938000, 128000, 6948000),
    "geilo-al": (127000, 6735000, 137000, 6745000),
}
for name, (x0, y0, x1, y1) in BOXES.items():
    path = f"{OUT}/{name}_box.geojson"
    ring = [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]
    with open(path, "w") as out:
        json.dump(
            {
                "type": "FeatureCollection",
                "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::25833"}},
                "features": [
                    {
                        "type": "Feature",
                        "properties": {},
                        "geometry": {"type": "Polygon", "coordinates": [ring]},
                    }
                ],
            },
            out,
        )
    for t in sys.argv[5:] or ["10", "2"]:
        r = subprocess.run(
            [
                PY,
                RASPUTIN,
                "mesh",
                "--dem",
                D,
                "--domain",
                path,
                "--tolerance",
                t,
                "--out",
                f"{OUT}/{name}_{t}.vtk",
            ],
            env={"PYTHONPATH": PKG},
            capture_output=True,
            text=True,
        )
        print(
            name,
            t,
            [
                line
                for line in (r.stdout + r.stderr).splitlines()
                if "triangles." in line or "rror" in line
            ],
        )
