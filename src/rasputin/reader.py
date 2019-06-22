from typing import Dict, Any, Optional,  Tuple, List
from dataclasses import dataclass
from enum import Enum
from PIL import Image
from PIL.TiffImagePlugin import TiffImageFile
from rasputin.geometry import GeoPolygon
import re

import numpy as np
from pathlib import Path
from logging import getLogger
from rasputin import triangulate_dem
from shapely.geometry import Polygon


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
    ProjCenterNorthingGeoKey = 3091  # Incorrectly named in document
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
        GeoKeyDirectory = np.asarray(image_tags[geo_key_tag],
                                     dtype="ushort").reshape((-1, 4))
    except KeyError:
        raise RuntimeError("Image is missing GeoKeyDirectory required by GeoTiff v1.0.")

    Header = GeoKeyDirectory[0, :]
    KeyDirectoryVersion = Header[0]
    KeyRevision = (Header[1], Header[2])
    NumberOfKeys = Header[3]

    assert KeyDirectoryVersion == 1  # Only existing version
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

def get_image_bounds(image: TiffImageFile) -> tuple:
    m, n = image.size

    scale_idx = GeoTiffTags.ModelPixelScaleTag.value
    delta_x, delta_y, _ = image.tag_v2.get(scale_idx, (1.0, 1.0, 0.0))

    tiepoint_idx = GeoTiffTags.ModelTiePointTag.value
    j_tag, i_tag, _, x_tag, y_tag, _ = image.tag_v2.get(tiepoint_idx, (0, 0, 0, 0, 0, 0))

    x_min = x_tag - delta_x * j_tag
    y_max = y_tag + delta_y * i_tag
    x_max = x_tag + delta_x * (n - 1 - j_tag)
    y_min = y_tag - delta_y * (m - 1 - i_tag)

    return (x_min, y_min, x_max, y_max)

class GeoKeysInterpreter(object):
    """
    This class provides functionality for converting a dict of GeoTIFF keys
    into a useable PROJ.4 projection specification.

    The string output can be used with pyproj to do coordinate transformations:

        import pyproj
        proj_str = GeoKeyInterpreter(geokeys).to_proj4()
        proj = pyproj.Proj(proj_str)

    NOTE:
    This class is incomplete and many GeoKeys are not handled. Extending the class with new handlers
    or expanding existing ones should be mostly straight forward, but some care needs to be taken in
    order to support angular units that are not degrees.
    """
    def __init__(self, geokeys):
        self.geokeys = geokeys
        self.dict = dict()
        self.flags = set()
        self.interpret()

    def interpret(self):
        # Try to extract a proj4 keyword argument for each (key, value) pair in the geokeys
        ignored_keys = []
        for (geokey_name, geokey_value) in self.geokeys.items():
            # Get handler
            try:
                handler = getattr(self, f"_{geokey_name}")
            except AttributeError:
                # Skip the GeoKey when no handler defined
                ignored_keys.append(geokey_name)
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

        if "south" in self.dict:
            if self.dict["south"]:
                self.flags.add("south")
            del self.dict["south"]
        logger = getLogger()
        logger.debug(f"Ignored GeoKeys: {ignored_keys}")

    def to_proj4(self) -> str:
        """
        Returns a string to be used with pyproj.
        """
        epsg = self.dict.get("EPSG")
        if epsg is not None:
            # The EPSG code completely specifies the projection. No more parameters needed
            proj4_str = f"+init=EPSG:{epsg}"

        else:
            parts = [f"+{name}={value}" for (name, value) in self.dict.items()]
            parts.extend([f"+{name}" for name in self.flags])
            proj4_str = " ".join(parts)

        # Add no_defs to proj4 string to avoid use of default values.
        # Pyproj should raise an error for incomplete specification (e.g. missing ellps)
        # Howerer, inconsistentencies may be ignored (e.g. when 'b' and 'rf' are both provided)
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
            return {"pm": float(value)}

    @staticmethod
    def _GeogSemiMajorAxisGeoKey(value):
        if np.isscalar(value):
            return {"a": float(value)}

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
        return {"units": unit_name}

    @staticmethod
    def _GeogCitationGeoKey(value):
        update = {}

        # Parse ellipsoid. GRS and WGS is probably enough
        gcs_re = {"GRS":    ("GRS[ ,_]{0,1}(19|)\\d\\d", "\\d\\d$"),
                  "WGS":    ("WGS[ ,_]{0,1}(19|)\\d\\d", "\\d\\d$"),
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


def identify_projection(*, image: TiffImageFile) -> str:
    geokeys = extract_geo_keys(image=image)
    return GeoKeysInterpreter(geokeys).to_proj4()


@dataclass
class Rasterdata:

    array: np.ndarray
    x_min: float
    y_max: float
    delta_x: float
    delta_y: float
    coordinate_system: str
    info: Dict[str, Any]

    @property
    def box(self):
        x_max = self.x_min + (self.array.shape[1] - 1) * self.delta_x
        y_min = self.y_max - (self.array.shape[0] - 1) * self.delta_y
        return (self.x_min, y_min, x_max, self.y_max)

    @property
    def to_cpp(self) -> triangulate_dem.raster_data_float:
        return triangulate_dem.raster_data_float(self.array,
                                                 self.x_min,
                                                 self.y_max,
                                                 self.delta_x,
                                                 self.delta_y)


def read_raster_file(*,
                     filepath: Path,
                     polygon: Optional[Polygon] = None) -> Rasterdata:
    logger = getLogger()
    logger.debug(f"Reading raster file {filepath}")
    assert filepath.exists()


    with Image.open(filepath) as image:
    # Determine image extents
        m, n = image.size

        # Use pixel coordinates if image does not define an affine transformation
        # This should only happen if the image is not a GeoTIFF
        model_pixel_scale = image.tag_v2.get(GeoTiffTags.ModelPixelScaleTag.value, (1.0, 1.0, 0.0))
        model_tie_point = image.tag_v2.get(GeoTiffTags.ModelTiePointTag.value, (0, 0, 0, 0, 0, 0))

        info = extract_geo_keys(image=image)
        coordinate_system = identify_projection(image=image)

        delta_x, delta_y, _ = model_pixel_scale
        j_tag, i_tag, _, x_tag, y_tag, _ = model_tie_point

        x_min = x_tag - delta_x * j_tag
        y_max = y_tag + delta_y * i_tag

        x_max = x_tag + delta_x * (n - 1 - j_tag)
        y_min = y_tag - delta_y * (m - 1 - i_tag)

        # If polygon is provided we crop the raster image to the polygon
        if polygon:
            # Get maximal extents of the polygon within the raster domain
            raster_domain = (Polygon([(x_min, y_min), (x_max, y_min),
                                      (x_max, y_max), (x_min, y_max)])
                             # buffer to avoid potential issues from CG roundoff errors
                             .buffer(0.1 * max(delta_x, delta_y)))
            #polygon = polygon.intersection(raster_domain)

            # Find indices for the box
            x_min_p, y_min_p, x_max_p, y_max_p = polygon.bounds
            box = (np.clip(np.floor((x_min_p - x_min)/delta_x), 0, n-1),
                   np.clip(np.floor((y_max - y_max_p)/delta_y), 0, n-1),
                   np.clip(np.ceil((x_max_p - x_min)/delta_x) + 1, 1, n),
                   np.clip(np.ceil((y_max - y_min_p)/delta_y) + 1, 1, n))

            # NOTE: Cropping does not include last indices
            image = image.crop(box=box)

            # Find extents of the cropped image
            x_min = x_tag + (box[0] - j_tag) * delta_x
            y_max = y_tag - (box[1] - i_tag) * delta_y

        # NOTE: Storing the array ensures the pointer to the data stays alive
        image_array = np.asarray(image)

    logger.debug("Done")

    return Rasterdata(array=image_array, x_min=x_min, y_max=y_max,
                      delta_x=delta_x, delta_y=delta_y, info=info,
                      coordinate_system=coordinate_system)


class RasterRepository:

    def __init__(self, *, directory: Path) -> None:
        self.directory = directory

    def get_intersections(self,
                          *,
                          target_polygon: GeoPolygon) -> List[Rasterdata]:

        parts = []

        raster_files = self.directory.glob("*.tif")
        for filepath in raster_files:
            geo_polygon = GeoPolygon.from_raster_file(filepath=filepath)

            if target_polygon.intersects(geo_polygon):
                polygon = target_polygon.transform(target_projection=geo_polygon.projection).polygon
                part = read_raster_file(filepath=filepath,
                                        polygon=polygon)
                parts.append(part)

                target_polygon = target_polygon.difference(geo_polygon)

                if target_polygon.polygon.area < 1e-10:
                    break
        return parts

    def coordinate_system(self, domain: GeoPolygon) -> str:
        raster_files = self.directory.glob("*.tif")
        for filepath in raster_files:
            geo_polygon = GeoPolygon.from_raster_file(filepath=filepath)
            if domain.intersects(geo_polygon):
                return geo_polygon.projection.definition_string()
        return ""

    def read(self,
             *,
             domain: GeoPolygon) -> Tuple[triangulate_dem.raster_list, triangulate_dem.simple_polygon]:

        data = triangulate_dem.raster_list()
        for part in self.get_intersections(target_polygon=domain):
            data.add_raster(part.to_cpp)
        cgal_polygon = domain.to_cpp
        return data, cgal_polygon


def read_sun_posisions(*, filepath: Path) -> triangulate_dem.shadow_vector:
    assert filepath.exists()
    with filepath.open("r") as ff:
        pass
    return triangulate_dem.ShadowVector()
