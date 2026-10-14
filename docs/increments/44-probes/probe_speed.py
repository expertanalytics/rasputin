# ruff: noqa: E501
"""Probe for increment 44, P8: what every mesh run pays for the catalogue as
TOML, and what `assign` pays for reading the CRS in `row_areas`.

Run as <venv python> probe_speed.py <tree>. Each figure is the median (and
the largest) of 15 fresh interpreters, so the first-call costs are in it.
"""
import json
import os
import statistics
import subprocess
import sys

os.environ["PROJ_NETWORK"] = "OFF"
tree = sys.argv[1]
py = sys.executable


def fresh(code: str) -> list[float]:
    """Run `code` in 15 fresh interpreters; each prints one number of ms."""
    out = []
    for _ in range(15):
        got = subprocess.run([py, "-c", code], capture_output=True, text=True, check=True)
        out.append(float(got.stdout.split()[-1]))
    return out


def show(label: str, ms: list[float]) -> None:
    print(f"{label}: median {statistics.median(ms):.3f} ms, largest {max(ms):.3f} ms")


# The draft catalogue: today's four entries written as TOML tables keyed by id,
# with GLO-30's two new fields, in the shape of the design's 4.2.
DRAFT = r'''
import tin_engine
assert tin_engine.__file__.startswith(TREE), tin_engine.__file__
from tin_engine.sources import SOURCES, STATION_SOURCES
def val(v):
    if isinstance(v, str): return json.dumps(v)
    if isinstance(v, list): return "[" + ", ".join(val(x) for x in v) + "]"
    return repr(v)
lines = []
for table, entries in (("source", SOURCES), ("station", STATION_SOURCES)):
    for key, e in entries.items():
        lines.append(f"[{table}.{key}]")
        d = e.model_dump(mode="json"); d.pop("id")
        if key == "glo30":
            d["tile_name"] = "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM"
            d["tile_spacing"] = [0.0027778, 0.00027778]
        lines += [f"{k} = {val(v)}" for k, v in d.items() if v is not None]
        lines.append("")
text = "\n".join(lines)
'''
setup = f"import json, time; TREE={tree!r}\n" + DRAFT
text = {}
exec(setup, text)
print(f"draft catalogue: {len(text['text'].encode())} bytes, {text['text'].count(chr(10))} lines")

# Load = parse the text, then validate every entry, as `catalogue` would; the
# models are today's: the two new fields are popped, and an absent `nodata`
# is given as None, the default the design adds.
LOAD = setup + r'''
from tin_engine.sources import RemoteSource, StationSource
import tomllib
t0 = time.perf_counter()
data = tomllib.loads(text)
for key, d in data["source"].items():
    d = {"nodata": None} | {k: v for k, v in d.items() if k not in ("tile_name", "tile_spacing")}
    RemoteSource.model_validate({"id": key, **d})
for key, d in data["station"].items():
    StationSource.model_validate({"id": key, **d})
print((time.perf_counter() - t0) * 1e3)
'''
show("P8a parse and validate the catalogue, first call", fresh(LOAD))

READ = "import time; from importlib import resources; import tin_engine.sources; t0 = time.perf_counter(); (resources.files('tin_engine') / 'data' / 'nve_hrd_2025.csv').read_text(encoding='utf-8'); print((time.perf_counter() - t0) * 1e3)"
show("P8a2 read a shipped 4.8 kB data file through importlib.resources", fresh(READ))

IMPORT_CLI = "import time; t0 = time.perf_counter(); import tin_engine.cli; print((time.perf_counter() - t0) * 1e3)"
show("P8b import tin_engine.cli today (every mesh run pays it)", fresh(IMPORT_CLI))

TOMLLIB = "import sys, time; import tin_engine.cli; print(0.0 if 'tomllib' in sys.modules else 1.0)"
print("P8c tomllib already imported by tin_engine.cli today:", fresh(TOMLLIB)[0] == 0.0)

# The CRS read the design adds to `row_areas`, once per chunk of 2**18 triangles.
CRS_READ = r'''
import time
from pyproj import Geod
from tin_engine.crs import parse_crs
t0 = time.perf_counter()
crs = parse_crs(CODE)
geographic = crs.is_geographic
e = crs.ellipsoid
geod = Geod(a=e.semi_major_metre, rf=e.inverse_flattening) if geographic else None
print((time.perf_counter() - t0) * 1e3)
'''
for code in ("EPSG:4326", "EPSG:3035"):
    show(f"P8d CRS read and Geod for {code}, first call", fresh(f"CODE={code!r}\n" + CRS_READ))
REPEAT = CRS_READ.replace("t0 = time.perf_counter()\n", "parse_crs(CODE).ellipsoid\nt0 = time.perf_counter()\n", 1)
show("P8e the same for EPSG:4326, second call in one process", fresh("CODE='EPSG:4326'\n" + REPEAT))

# Against the per-row loop `row_areas` runs today on a window of 1,000 rows.
ROWS = r'''
import time
from tin_engine.io.models import ClassWindow
import numpy as np
from tin_engine.land_cover_fractions import row_areas
w = ClassWindow(array=np.zeros((1000, 4), np.uint8), x_min=-40.0, y_max=-10.0,
                delta_x=0.00026949458523585647, delta_y=0.00026949458523585647)
t0 = time.perf_counter()
row_areas(w)
print((time.perf_counter() - t0) * 1e3)
'''
try:
    show("P8f row_areas today on a 1,000-row MapBiomas window", fresh(ROWS))
except subprocess.CalledProcessError as exc:
    print("P8f could not build the window:", exc.stderr.strip().splitlines()[-1])
print(json.dumps({"python": sys.version.split()[0]}))
