from typing import Tuple, Dict, Any, Optional
from enum import Enum
from PIL import Image
from PIL.TiffImagePlugin import TiffImageFile
import json
import re
from itertools import islice
import numpy as np
from pathlib import Path
from logging import getLogger
from rasputin import triangulate_dem

class GeoTiffTags(Enum):
    ModelTiePointTag = 33922
    ModelPixelScaleTag = 33550
    GeoKeyDirectoryTag = 34735

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


def _isinteger(obj: Any) -> bool:
    return np.issubdtype(type(obj), np.integer)


def extract_geo_keys(*, image: TiffImageFile) -> Dict[str, Any]:
    """ Extract GeoKeys from image and return as Python dictionary. """

    # Using .tag_v2 instead of legacy .tag
    image_tags = image.tag_v2
    geo_key_tag = GeoTiffTags.GeoKeyDirectoryTag.value
    try:
        GeoKeyDirectory = np.asarray(image_tags[geo_key_tag], dtype="ushort").reshape((-1,4))
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


class GeoKeysParser(object):
    """
    This class provides functionality for converting a dict of GeoTIFF keys
    into a useable PROJ.4 projection specification.

    The string output can be used with pyproj to do coordinate transformations:

        import pyproj
        proj_str = GeoKeyParser(geokeys).to_proj4()
        proj = pyproj.Proj(proj_str)


    NOTE:
    This class is incomplete and many GeoKeys are not handled. However, extending the class
    with new handlers or expanding existing ones should be mostly straight forward, but some
    care needs to be taken in order to support geographical systems not having degrees as unit.
    """
    def __init__(self, geokeys):
        self.geokeys = geokeys
        self.dict = dict()
        self.parse()


    def parse(self):
        # Try to extract a proj4 keyword argument for each (key, value) pair in the geokeys
        unhandled = []
        for (geokey_name, geokey_value) in self.geokeys.items():
            # Get handler
            try:
                handler = getattr(self, f"_{geokey_name}")
            except AttributeError:
                # Skip the GeoKey since no handler defined
                unhandled.append(geokey_name)
                continue

            # Get a dict update from handler
            update = handler(geokey_value)
            if update is None:
                continue

            # Raise error on conflicts
            for (name, val) in update.items():
                if name in self.dict:
                    old_val = self.dict[name]
                    if old_val != val:
                        raise ValueError(f"Conflicting values for {name}: {old_val} and {val}")

                else:
                    self.dict[name] = val


    def to_proj4(self) -> str:
        epsg = self.dict.get("EPSG")
        if epsg is not None:
            # The EPSG code completely specifies the projection. No more parameters needed
            proj4_str = f"+init=EPSG:{epsg}"

        else:
            parts = [f"+{name}={value}" for (name, value) in self.dict.items()]
            proj4_str = " ".join(parts)

        # Add no_defs to proj4 string to avoid use of default values.
        # Pyproj should raise an error for incomplete specification (e.g. missing ellps)
        # Howerer, inconsistentencise ignored (e.g. when 'b' and 'f' are both provided)
        proj4_str += " +no_defs"

        return proj4_str


    @staticmethod
    def _ProjectedCSTypeGeoKey(value):
        if _isinteger(value) and 20000 <= value <= 32760:
            return {"EPSG": value}


    @staticmethod
    def _GeoGraphicTypeGeo(value):
        if _isinteger(value) and 4000 <= value < 5000:
            return {"EPSG": value}


    @staticmethod
    def _GeogInvFlatteningGeoKey(value):
        if np.isscalar(value):
            return {"rf": float(value)}


    @staticmethod
    def _GeogPrimeMeridianLongGeoKey(value):
        if np.isscalar(value):
            return {"pm" : float(value)}


    @staticmethod
    def _GeogSemiMajorAxisGeoKey(value):
        if np.isscalar(value):
            return {"a" : float(value)}


    @staticmethod
    def _ProjectionGeoKey(value):
        if _isinteger(value) and 16000 <= value < 16200:
            # UTM projections are denoted
            #    160zz (nothern hemisphere) or
            #    161zz
            # where zz in zone number

            zone = value % 16000
            south = bool(zone // 100)
            if south:
                zone %= 100

            return dict(proj="utm", zone=zone, south=south)


    @staticmethod
    def _ProjLinearUnitsGeoKey(value):
        # More units could be supported here
        unit_name = {9001: "m"}[value]
        return {"units" : unit_name}


    @staticmethod
    def _GeogCitationGeoKey(value):
        update = {}

        # Parse ellipsoid. GRS and WGS is probably enough
        gcs_re = {"GRS":    ("GRS[ ,_]{0,1}(19|)\d\d", "\d\d$"),
                  "WGS":    ("WGS[ ,_]{0,1}(19|)\d\d", "\d\d$"),
                  "sphere": ("sphere", None)}

        for (name, (pattern, postfix_pattern)) in gcs_re.items():
            s = re.search(pattern, value, flags=2)
            if s:
                # Add postfix if applicable
                if postfix_pattern:
                    postfix = re.search(postfix_pattern, s.group()).group()
                    name = name + postfix

                update["ellps"] = name
                break

        # TODO: Could parse for datum and prime meridian here as well
        return update


def get_projection(image: TiffImageFile) -> str:
    geokeys = extract_geokeys(image)
    return GeoKeysParser(geokeys).to_proj4()


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
                     stop_y: Optional[int]) -> Tuple[triangulate_dem.PointVector, Dict[str, Any]]:
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

        model_pixel_scale = image.tag_v2.get(GeoTiffTags.ModelPixelScaleTag.value)
        model_tie_point = image.tag_v2.get(GeoTiffTags.ModelTiePointTag.value)
        info = extract_geo_keys(image=image)

        d = np.fromiter(img_slice(image, start_y, stop_y, start_x, stop_x), dtype="float", count=count).reshape(stop_y - start_y, -1)
        #d = np.array(image.getdata()).reshape(image.size[0], image.size[1])[start_x:stop_x, start_y:stop_y]

    dx = dy = 1.0
    if model_pixel_scale:
        dx, dy, _ = model_pixel_scale

    x0 = y0 = 0.0
    I = J = 0
    if model_tie_point:
        J, I, _, x0, y0, _ = model_tie_point

    raster_coordinates = triangulate_dem.PointVector()
    for (i, j), h in np.ndenumerate(d):
        raster_coordinates.append([x0 + (start_x + j - J)*dx, y0 - (start_y + i - I)*dy, h])
    logger.debug("Done")
    return raster_coordinates, info


def read_sun_posisions(*, filepath: Path) -> triangulate_dem.ShadowVector:
    assert filepath.exists()
    with filepath.open("r") as ff:
        pass
    return triangulate_dem.ShadowVector()
