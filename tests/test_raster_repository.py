import pytest
from pathlib import Path
import os
import pyproj
from rasputin.reader import RasterRepository

if "RASPUTIN_DATA_DIR" in os.environ:
    data_dir = Path(os.environ["RASPUTIN_DATA_DIR"])
else:
    data_dir = Path(os.environ["HOME"]) /"projects" / "rasputin_data" / "dem_archive"


def test_get():
    x0 = 8.530918
    y0 = 60.898468
    repo = RasterRepository(directory=data_dir)
    input_coordinate_system = pyproj.Proj(init="EPSG:4326").definition_string()
    target_coordinate_system = pyproj.Proj(init="EPSG:32633").definition_string()
    result = repo.read(x=x0,
                       y=y0,
                       dx=5000,
                       dy=5000,
                       input_coordinate_system=input_coordinate_system,
                       target_coordinate_system=target_coordinate_system)
    assert result
