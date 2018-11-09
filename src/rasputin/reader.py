from PIL import Image
import json
import numpy as np
from pathlib import Path
from logging import getLogger

from rasputin import triangulate_dem


def read_raster_file(*,
                     filepath: Path,
                     start_x: int=0,
                     stop_x: int=-1,
                     start_y: int=0,
                     stop_y: int=-1) -> triangulate_dem.PointVector:
    logger = getLogger()
    logger.debug(f"Reading raster file {filepath}")
    assert filepath.exists()
    with Image.open(filepath) as image:
        dx, dy, _ = image.tag.named().get("ModelPixelScaleTag", (1, 1, 0))
        d = np.array(image.getdata()).reshape(image.size[0], image.size[1])[start_x:stop_x, start_y:stop_y]
        model_tie_point = image.tag.named().get("ModelTiepointTag", None)
    x0 = y0 = 0.0
    if model_tie_point:
        x0 = model_tie_point[3]
        y0 = model_tie_point[4]
    raster_coordinates = triangulate_dem.PointVector()
    for (i, j), h in np.ndenumerate(d):
        raster_coordinates.append([x0 + (start_x + i)*dx, y0 + (start_y + j)*dy, h])
    logger.debug("Done")
    return raster_coordinates


def read_sun_posisions(*, filepath: Path) -> triangulate_dem.ShadowVector:
    assert filepath.exists()
    with filepath.open("r") as ff:
        pass
    return triangulate_dem.ShadowVector()
