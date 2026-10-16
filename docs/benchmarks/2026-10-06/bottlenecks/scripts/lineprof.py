"""Line-profile the three phases of `rasputin mesh` in-process.

Wraps functions at runtime with line_profiler; production code is unchanged.
usage: lineprof.py <catchment> <out.txt>
"""
import sys
from line_profiler import LineProfiler
import tin_engine.cli as cli
import tin_engine.feature_input as fi
import tin_engine.landcover as lc
import tin_engine.dem_input as di
import tin_engine.mosaic as mo
import tin_engine.io.cog as cog

c, out = sys.argv[1], sys.argv[2]
lp = LineProfiler()
for mod, name in [(fi, "pre_clip"), (fi, "_runs"), (lc, "regions"), (lc, "_lookup"),
                  (di, "open_dem"), (mo, "assemble"), (mo, "_decide"), (cog, "decode_window")]:
    w = lp(getattr(mod, name)); setattr(mod, name, w)
fi._Tally._take = lp(fi._Tally._take)
w = lp(lc.label_triangles); lc.label_triangles = w; cli.label_triangles = w
di.assemble = mo.assemble
if hasattr(cli, "open_dem"): cli.open_dem = di.open_dem
D = "/Users/skavhaug/projects/rasputin_data"
args = ["mesh", "--dem", f"{D}/DTM10_UTM33_20260925",
        "--domain", f"/Users/skavhaug/projects/rasputin_scratch/norway/{c}/{c}_outline_nve.geojson",
        "--features", f"{D}/corine2018_dtm10_utm33.gpkg", "--features-layer", "corine2018",
        "--features-map", "corine", "--tolerance", "10", "--binary",
        "--out", out.replace(".txt", ".vtk")]
try:
    cli.app(args, standalone_mode=False)
finally:
    with open(out, "w") as f:
        lp.print_stats(stream=f, output_unit=1e-3, stripzeros=True)
