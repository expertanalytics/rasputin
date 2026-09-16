import pytest
from pathlib import Path
import os
import pyproj
from shapely.geometry import Polygon
from rasputin.geometry import GeoPolygon
from rasputin.reader import RasterRepository

assert "RASPUTIN_DATA_DIR" in os.environ, "You need to set RASPUTIN_DATA_DIR prior to running this test"
data_dir = os.environ["RASPUTIN_DATA_DIR"]


def test_get():
    x0 = 8.54758671814368
    y0 = 60.898468
    repo = RasterRepository(directory=Path(data_dir) / "dem_archive")
    input_crs = pyproj.CRS.from_string("+init=EPSG:4326")
    polygon = Polygon.from_bounds(xmin=x0 - 0.1, xmax=x0 + 0.1, ymin=y0 - 0.1, ymax=y0 + 0.1)
    geo_polygon = GeoPolygon(crs=input_crs, polygon=polygon)
    result = repo.read(domain=geo_polygon)
    assert result
