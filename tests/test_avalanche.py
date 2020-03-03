from pyproj import CRS
from datetime import datetime, timedelta
from rasputin.geometry import GeoPolygon
from shapely.geometry import Polygon
from rasputin.avalanche import build_regions, ARegion, BRegion, get_avalanche_warnings


def test_read_avalanche_regions():
    crs = CRS.from_epsg(32633)
    assert len(build_regions(region_type=ARegion, crs=crs)) > 20
    assert len(build_regions(region_type=BRegion, crs=crs)) > 10


def test_read_avalanche_regions_with_domain():

    poly = Polygon(((8.95, 61.50),
                    (8.95, 61.40),
                    (8.90, 61.40),
                    (8.90, 61.50),
                    (8.95, 61.50)))
    start_time = datetime.now()
    end_time = start_time + timedelta(days=2)


    crs = CRS.from_string("+init=EPSG:4326")
    domain = GeoPolygon(polygon=poly, crs=crs)
    regions = build_regions(region_type=ARegion, crs=crs, domain=domain)
    print(len(regions))
    warnings = get_avalanche_warnings(start_date=start_time, end_date=end_time, regions=regions)
    print(len(warnings))
