# ruff: noqa: E501, B905
"""Increment 44 probe 4: a Geod from the raster CRS's ellipsoid against today's Geod(ellps="WGS84")."""
from pyproj import CRS, Geod

ref = Geod(ellps="WGS84")
dy, y = 0.00026949458523585647, 5.422904791408522  # the MapBiomas Collection 11 header

def rows(g):
    return [abs(g.polygon_area_perimeter([0, dy, dy, 0], [y - r * dy] * 2 + [y - r * dy - dy] * 2)[0])
            for r in range(0, 146501, 7)]

for label, g in [("CRS(EPSG:4326).get_geod()", CRS("EPSG:4326").get_geod()),
                 ("Geod(a=, rf=) from CRS(EPSG:4326).ellipsoid",
                  Geod(a=CRS("EPSG:4326").ellipsoid.semi_major_metre, rf=CRS("EPSG:4326").ellipsoid.inverse_flattening))]:
    a, b = rows(ref), rows(g)
    print(f"{label}: f {g.f!r} (ellps=WGS84: {ref.f!r}); MapBiomas rows sampled {len(a)}, differing {sum(p != q for p, q in zip(a, b))}")
s = CRS("+proj=longlat +R=6371000 +no_defs").ellipsoid
print("a sphere's ellipsoid: inverse_flattening", s.inverse_flattening)
