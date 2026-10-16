import sys, pathlib, shapely, numpy as np
import tin_engine.feature_input as fi
from tin_engine.domain import read_domain
from tin_engine.crs import parse_crs
c = sys.argv[1]
dom = read_domain(pathlib.Path(f"/Users/skavhaug/projects/rasputin_scratch/norway/{c}/{c}_outline_nve.geojson")); dem = parse_crs("EPSG:25833")
region = fi.source_region(dom, dem, dem); shapely.prepare(region)
rows = fi.read_source(pathlib.Path("/Users/skavhaug/projects/rasputin_data/corine2018_dtm10_utm33.gpkg"), "corine2018", "Code_18", lambda own: region.bounds).rows
g = np.array([r[1] for r in rows])
box = shapely.box(*region.bounds)
nv = shapely.get_num_coordinates(g)
hit_box = shapely.intersects(g, box); hit_reg = shapely.intersects(region, g); hit_dom = shapely.intersects(dom.polygon, g)
print(c, "read", len(g), "| envelope meets hull bbox", int(shapely.intersects(shapely.envelope(g), box).sum()),
      "| geometry meets hull bbox", int(hit_box.sum()), "| meets hull", int(hit_reg.sum()), "| meets domain", int(hit_dom.sum()),
      "| vertices read", int(nv.sum()), "| vertices of features meeting hull", int(nv[hit_reg].sum()),
      "| domain area km2", round(dom.polygon.area/1e6), "| hull area km2", round(region.area/1e6), "| hull bbox area km2", round(box.area/1e6))
