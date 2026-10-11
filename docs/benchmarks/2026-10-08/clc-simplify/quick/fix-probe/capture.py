import pickle, sys
import tin_engine.feature_input as fi
from tin_engine import cli
real = fi.snap_to_outline
def cap(polys, outline, d):
    pickle.dump((polys, outline, d), open(sys.argv[1], "wb"))
    raise SystemExit(0)
fi.snap_to_outline = cap
D = "/Users/skavhaug/projects/rasputin_data"
args = ["mesh", "--dem", f"{D}/DTM10_UTM33_20260925", "--domain", f"{D}/numedalslagen_outline_nve.geojson",
        "--features", f"{D}/corine2018_dtm10_utm33.gpkg", "--features-layer", "corine2018",
        "--features-map", "corine", "--tolerance", "10", "--out", "x.vtk"] + sys.argv[2:]
cli.app(args)
