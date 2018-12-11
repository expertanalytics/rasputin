from typing import Tuple, Dict, Any, Optional
from enum import Enum
from PIL import Image
import json
from itertools import islice
import numpy as np
from pathlib import Path
from logging import getLogger
from rasputin import triangulate_dem


class AsciiGeoKeys(Enum):
    GTCitationGeoKey = 1026
    GeogCitationGeoKey = 2049


class CodedGeoKeys(Enum):
    GTModelTypeGeoKey = 1024
    GTRasterTypeGeoKey = 1025
    GeogAngularUnitsGeoKey = 2054
    ProjectedCSTypeGeoKey = 3072
    ProjLinearUnitsGeoKey = 3076


def extract_geo_keys(*, named_tags: dict) -> Dict[str, Any]:
    geo_key_directory = np.array(named_tags["GeoKeyDirectoryTag"]).reshape(-1, 4)
    geo_keys = geo_key_directory[:, 0]
    info = {}
    for k in AsciiGeoKeys:
        if k.value in geo_keys:
            row = geo_key_directory[np.where(k.value == geo_keys)].flatten()
            if row[1] == 34737:
                info[k.name] = named_tags["GeoAsciiParamsTag"][0][row[3]: row[3] + row[2]]
    for k in CodedGeoKeys:
        if k.value in geo_keys:
            row = geo_key_directory[np.where(k.value == geo_keys)].flatten()
            info[k.name] = row[3]
    return info


def img_slice(img, start_i, stop_i, start_j, stop_j):
    myiter = iter(img.getdata())
    m, n = img.size
    if start_i > 0:
        try:
            next(islice(myiter, start_i*n, start_i*n))
        except StopIteration:
            pass
    for i in range(stop_i - start_i):
        if stop_i <= i < start_i:
            try:
                next(islice(myiter, n, n))
            except StopIteration:
                pass
        else:
            for v in islice(myiter, start_j, stop_j):
                yield v
            try:
                next(islice(myiter, n - stop_j, n - stop_j))
            except StopIteration:
                pass



def read_raster_file(*,
                     filepath: Path,
                     start_x: int,
                     stop_x: Optional[int],
                     start_y: int,
                     stop_y: Optional[int]) -> Tuple[triangulate_dem.PointVector, Dict[str, str]]:
    logger = getLogger()
    logger.debug(f"Reading raster file {filepath}")
    assert filepath.exists()
    with Image.open(filepath) as image:
        if None in [stop_x, stop_y]:
            m, n = image.size
            if stop_x is None:
                stop_x = n
            if stop_y is None:
                stop_y = m
        count = (stop_x - start_x)*(stop_y - start_y)
        named_tags = image.tag.named()
        dx, dy, _ = named_tags.get("ModelPixelScaleTag", (1, 1, 0))
        info = extract_geo_keys(named_tags=named_tags)
        model_tie_point = named_tags.get("ModelTiepointTag", None)
        d = np.fromiter(img_slice(image, start_y, stop_y, start_x, stop_x), dtype="float", count=count).reshape(stop_y - start_y, -1)
        #d = np.array(image.getdata()).reshape(image.size[0], image.size[1])[start_x:stop_x, start_y:stop_y]
    x0 = y0 = 0.0
    if model_tie_point:
        x0 = model_tie_point[3]
        y0 = model_tie_point[4]
    raster_coordinates = triangulate_dem.PointVector()
    for (i, j), h in np.ndenumerate(d):
        raster_coordinates.append([x0 + (start_x + j)*dx, y0 - (start_y + i)*dy, h])
    logger.debug("Done")
    return raster_coordinates, info


def read_sun_posisions(*, filepath: Path) -> triangulate_dem.ShadowVector:
    assert filepath.exists()
    with filepath.open("r") as ff:
        pass
    return triangulate_dem.ShadowVector()
