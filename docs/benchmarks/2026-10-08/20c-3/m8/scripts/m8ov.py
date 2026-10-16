"""M8: is the input cover overlapping? sum of areas minus area of the union, before and after."""
import sys, numpy as np, shapely
for f in sys.argv[1:]:
    a, b = np.load(f, allow_pickle=True)
    for tag, G in (("before", shapely.from_wkb(list(a))), ("after", shapely.from_wkb(list(b)))):
        s = shapely.area(G).sum(); u = shapely.union_all(G).area
        print(f"{f.split('/')[-1]} {tag}: sum of areas - union = overlap {s - u:.2f} m2; valid coverage {shapely.coverage_is_valid(G)}")
