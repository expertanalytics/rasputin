"""Run `rasputin mesh` in-process and save label_triangles' inputs (npz + WKB)."""
import sys, pickle, numpy as np, shapely
import tin_engine.cli as cli
import tin_engine.landcover as lc
c = sys.argv[1]; S = sys.argv[2]
orig = lc.label_triangles
def grab(vertices, triangles, edges, *, polygons, margin):
    np.savez(f"{S}/lc_{c}.npz", vertices=np.asarray(vertices), triangles=np.asarray(triangles), edges=np.asarray(edges), margin=margin)
    with open(f"{S}/lc_{c}_polys.pkl", "wb") as f:
        pickle.dump([(shapely.to_wkb(p), code) for p, code in polygons], f)
    out = orig(vertices, triangles, edges, polygons=polygons, margin=margin)
    np.save(f"{S}/lc_{c}_codes.npy", out.codes)
    return out
cli.label_triangles = grab
D = "/Users/skavhaug/projects/rasputin_data"
cli.app(["mesh", "--dem", f"{D}/DTM10_UTM33_20260925", "--domain",
  f"/Users/skavhaug/projects/rasputin_scratch/norway/{c}/{c}_outline_nve.geojson", "--features",
  f"{D}/corine2018_dtm10_utm33.gpkg", "--features-layer", "corine2018", "--features-map", "corine",
  "--tolerance", "10", "--binary", "--out", f"{S}/cap_{c}.vtk"], standalone_mode=False)
