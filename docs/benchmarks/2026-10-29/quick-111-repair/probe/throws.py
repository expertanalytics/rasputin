"""Which invalid shapes make GEOS's intersection with a region throw."""
import shapely
from shapely import MultiPolygon, Polygon, box

reg = box(-10, -10, 10, 10).buffer(1)
cases = {
    "bowtie": Polygon([(0, 0), (4, 4), (4, 0), (0, 4)]),
    "hole crosses shell": Polygon([(0, 0), (4, 0), (4, 4), (0, 4)], [[(1, 1), (5, 1), (5, 2), (1, 2)]]),
    "overlapping parts": MultiPolygon([box(0, 0, 2, 2), box(1, 1, 3, 3)]),
    "hole outside shell": Polygon([(0, 0), (4, 0), (4, 4), (0, 4)], [[(6, 6), (7, 6), (7, 7), (6, 7)]]),
    "spike (zero-width)": Polygon([(0, 0), (4, 0), (4, 4), (2, 4), (2, 6), (2, 4), (0, 4)]),
}
for k, g in cases.items():
    try:
        r = shapely.intersection(g, reg)
        out = f"no throw, area {r.area:.3f} (input area {g.area:.3f})"
    except shapely.errors.GEOSException as e:
        out = f"THROWS {str(e)[:60]}"
    print(f"{k:22s} valid={g.is_valid}  {out}")
