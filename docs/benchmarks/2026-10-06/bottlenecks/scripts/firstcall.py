import sys, time, pathlib, importlib
import tin_engine.dem_input as di
from tin_engine.domain import read_domain
pre = sys.argv[2] == "warm"
if pre:
    t = time.perf_counter()
    for m in ("imagecodecs", "tifffile"): importlib.import_module(m)
    print("imports", round(time.perf_counter() - t, 3))
c = sys.argv[1]
dom = read_domain(pathlib.Path(f"/Users/skavhaug/projects/rasputin_scratch/norway/{c}/{c}_outline_nve.geojson"))
for i in range(2):
    t = time.perf_counter(); di.open_dem(di.DemRequest(sources=[pathlib.Path("/Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925")], domain=dom)); print("open_dem run", i + 1, round(time.perf_counter() - t, 3))
mods = sorted(m for m in sys.modules if m.startswith(("imagecodecs", "tifffile")))
print("modules", len(mods))
