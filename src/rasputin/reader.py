from typing import Tuple, Dict, Any, Optional
from enum import Enum
from PIL import Image
from PIL.TiffImagePlugin import TiffImageFile
import json
from itertools import islice
import numpy as np
from pathlib import Path
from logging import getLogger
from rasputin import triangulate_dem


class KeyValueTags(Enum):
    # Parameter value types:
    #    'short'  values are stored in the geokey directory
    #    'double'                   in tag 34736 at given offset
    #    'ascii'                    in tag 34737 at offset:offset+count
    GeoShortParamsTag = 0
    GeoDoubleParamsTag = 34736
    GeoAsciiParamsTag = 34737


class GeoKeys(Enum):
    # Configuration keys
    GTModelTypeGeoKey = 1024
    GTRasterTypeGeoKey = 1025
    GTCitationGeoKey = 1026

    # Geographic CS parameter keys
    GeographicTypeGeoKey = 2048
    GeogCitationGeoKey = 2049
    GeogGeodeticDatumGeoKey = 2050
    GeogPrimeMeridianGeoKey = 2051
    GeogPrimeMeridianLongGeoKey = 2061
    GeogLinearUnitSizeGeoKey = 2053
    GeogAngularUnitsGeoKey = 2054
    GeogAngularUnitSizeGeoKey = 2055
    GeogEllipsoidGeoKey = 2056
    GeogSemiMajorAxisGeoKey = 2057
    GeogSemiMinorAxisGeoKey = 2058
    GeogInvFlatteningGeoKey = 2059
    GeogAzimuthUnitsGeoKey = 2060

    # Projected CS parameters keys
    ProjectedCSTypeGeoKey = 3072
    PCSCitationGeoKey = 3073

    # Projection definition keys
    ProjectionGeoKey = 3074
    ProjCoordTransGeoKey = 3075
    ProjLinearUnitsGeoKey = 3076
    ProjLinearUnitSizeGeoKey = 3077
    ProjStdParallel1GeoKey = 3078
    ProjStdParallel2GeoKey = 3079
    ProjNatOriginLongGeoKey = 3080
    ProjNatOriginLatGeoKey = 3081
    ProjFalseEastingGeoKey = 3082
    ProjFalseNorthingGeoKey = 3083
    ProjFalseOriginLongGeoKey = 3084
    ProjFalseOriginLatGeoKey = 3085
    ProjFalseOriginEastingGeoKey = 3086
    ProjFalseOriginNorthingGeoKey = 3087
    ProjCenterLongGeoKey = 3088
    ProjCenterLatGeoKey = 3089
    ProjCenterEastingGeoKey = 3090
    ProjCenterNorthingGeoKey = 3091 # Incorrectly named in document
    ProjScaleAtNatOriginGeoKey = 3092
    ProjScaleAtCenterGeoKey = 3093
    ProjAzimuthAngleGeoKey = 3094
    ProjStraightVertPoleLongGeoKey = 3095


def extract_geo_keys(*, image: TiffImageFile) -> Dict[str, Any]:
    """ Extract GeoKeys from image and return as Python dictionary. """

    # Using .tag_v2 instead of legacy .tag
    image_tags = image.tag_v2
    try:
        GeoKeyDirectory = np.asarray(image_tags[34735], dtype="ushort").reshape((-1,4))
    except KeyError:
        raise RuntimeError("Image is missing GeoKeyDirectory required by GeoTiff v1.0.")

    Header = GeoKeyDirectory[0,:]
    KeyDirectoryVersion = Header[0]
    KeyRevision = (Header[1], Header[2])
    NumberOfKeys = Header[3]

    assert KeyDirectoryVersion == 1 # Only existing version
    assert len(GeoKeyDirectory) == NumberOfKeys + 1

    # Read all the geokey fields
    geo_keys = {}
    for (key_id, location, count, value_offset) in GeoKeyDirectory[1:]:
        key_name = GeoKeys(key_id).name

        if KeyValueTags(location) == KeyValueTags.GeoShortParamsTag:
            # Short values are stored in value_offset field in the directory
            key_value = value_offset

            # Translate special integer values for readability:
            if key_value == 0:
                key_value = "undefined"

            elif key_value == 32767:
                key_value = "user-defined"

        if KeyValueTags(location) == KeyValueTags.GeoDoubleParamsTag:
            tiff_tag = KeyValueTags.GeoDoubleParamsTag.value
            offset = value_offset
            key_value = image_tags[tiff_tag][offset]

        if KeyValueTags(location) == KeyValueTags.GeoAsciiParamsTag:
            tiff_tag = KeyValueTags.GeoAsciiParamsTag.value
            offset = value_offset
            ascii_val = image_tags[tiff_tag][offset:offset + count]
            key_value = ascii_val.replace("|", "\n").strip()

        geo_keys[key_name] = key_value

    return geo_keys


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
