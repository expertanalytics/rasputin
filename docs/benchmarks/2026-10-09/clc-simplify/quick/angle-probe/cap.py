import pickle
import tin_engine.feature_input as fi
from tin_engine import cli
real = fi.simplify_borders
def wrap(polys, band, *a, **kw):
    out = real(polys, band, *a, **kw); kw = dict(kw, args=a)
    pickle.dump((list(polys), out.polygons, band, kw, fi_domain[0]), open("simp.pkl", "wb"))
    return out
fi_domain = [None]
orig_clean = fi._Reader._clean if hasattr(fi, "_Reader") else None
fi.simplify_borders = wrap
D = "/Users/skavhaug/projects/rasputin_data"
import shapely, json
dom = shapely.geometry.shape(json.load(open(f"{D}/numedalslagen_outline_nve.geojson"))["features"][0]["geometry"])
fi_domain[0] = dom
try:
    cli.app(["mesh", "--dem", f"{D}/DTM10_UTM33_20260925", "--domain", f"{D}/numedalslagen_outline_nve.geojson",
         "--features", f"{D}/corine2018_dtm10_utm33.gpkg", "--features-layer", "corine2018",
         "--features-map", "corine", "--tolerance", "10", "--ascii", "--out", "fix2.vtk"])
except SystemExit: pass
