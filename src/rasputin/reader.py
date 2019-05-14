from typing import Tuple, Dict, Any, Optional, List
from dataclasses import make_dataclass, dataclass
from enum import Enum
from PIL import Image, TiffImagePlugin
from PIL.TiffImagePlugin import TiffImageFile
from shapely import geometry
from shapely.ops import snap
import pyproj
import re

from itertools import islice
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
    def _cpp(self):
        return triangulate_dem.RasterData_float(self.array,
                                                self.x_min, self.y_max,
                                                self.delta_x, self.delta_y)


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

        # If polygon is provided we crop the raster image to the polygon
        if polygon:
            # Get maximal extents of the polygon
            x_min_p, y_min_p, x_max_p, y_max_p = polygon.bounds

            # Find indices for the box
            box = (min(max(0, int((x_min_p - x_min)/delta_x)), n-1),
                   min(max(0, int((y_max - y_max_p)/delta_y)), m-1),
                   min(max(1, int((x_max_p - x_min)/delta_x) + 2), n),
                   min(max(1, int((y_max - y_min_p)/delta_y) + 2), m))


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


class GeoPolygon:

    def __init__(self, *, polygon: geometry.Polygon, proj: pyproj.Proj):
        self.polygon = polygon
        self.proj = proj


class RasterRepository:

    def __init__(self, *, directory: Path):
        self.directory = directory
        self.shapes = {}  # Used for debugging

    def read(self, *,
             x: float,
             y: float,
             dx: float,
             dy: float,
             input_coordinate_system: str,
             target_coordinate_system: str) -> triangulate_dem.PointVector:
        input_proj = pyproj.Proj(input_coordinate_system)
        target_proj = pyproj.Proj(target_coordinate_system)

        target_x, target_y = pyproj.transform(input_proj, target_proj, x, y)
        target_bbox = geometry.Polygon.from_bounds(round(target_x - dx),
                                                   round(target_y - dy),
                                                   round(target_x + dx),
                                                   round(target_y + dy))
        remainding_bbox = geometry.Polygon(target_bbox)
        self.shapes["target_bbox"] = target_bbox

        files = self._extract_files(x=x, y=y, dx=dx, dy=dy, coordinate_system=input_proj)
        if not files:
            raise RuntimeError("No raster files found for given center point.")
        pts = triangulate_dem.PointVector()
        for i, file in enumerate(files):
            with Image.open(file) as img:
                img_geo_bbox = self._get_bounding_box(image=img)
                self.shapes[f"{file}_bbox"] = img_geo_bbox
                if target_proj == img_geo_bbox.proj:
                    mapped_remainder = remainding_bbox
                else:
                    x = np.round(remainding_bbox.bounds[0::2])
                    y = np.round(remainding_bbox.bounds[1::2])
                    xm, ym = pyproj.transform(target_proj, img_geo_bbox.proj, x, y)
                    xm = np.round(xm)
                    ym = np.round(ym)
                    mapped_remainder = geometry.Polygon.from_bounds(xm[0], ym[0], xm[1], ym[1])
                if not mapped_remainder.intersects(img_geo_bbox.polygon):
                    continue
                img_bbox = snap(img_geo_bbox.polygon, mapped_remainder, 2)
                self.shapes[f"{file}_mapped_remainder"] = mapped_remainder
                mapped_intersection = mapped_remainder.intersection(img_bbox)
                self.shapes[f"{file}_mapped_intersection"] = mapped_intersection
                pts.extend(self._extract_coords(image=img, bounds=mapped_intersection.bounds))
                x0, y0, x1, y1 = np.round(mapped_intersection.bounds)
                x_m, y_m = pyproj.transform(img_geo_bbox.proj, target_proj, [x0, x1], [y0, y1])
                x_m = np.round(x_m)
                y_m = np.round(y_m)
                intersection = snap(geometry.Polygon.from_bounds(x_m[0],
                                                                 y_m[0],
                                                                 x_m[1],
                                                                 y_m[1]), remainding_bbox, 2)
                self.shapes[f"{file}_intersection"] = intersection
                remainding_bbox = remainding_bbox.difference(intersection)
                if not remainding_bbox:
                    break
                self.shapes[f"{file}_remainding_bbox"] = remainding_bbox
        return pts

    def _extract_coords(self, *,
                        image: TiffImagePlugin.TiffImageFile,
                        bounds: Tuple[float, float, float, float]) -> triangulate_dem.PointVector:
        m, n = image.size
        x_min, y_min, x_max, y_max = bounds
        dx, dy, _ = image.tag_v2.get(GeoTiffTags.ModelPixelScaleTag.value)
        X0, _, _, Y1 = self._get_bounding_box(image=image).polygon.bounds

        # Here we actually need to snap coordinates back to bounds...
        j0 = int(np.floor((x_min - X0)/dx))
        j1 = int(np.ceil((x_max - X0)/dx))
        i0 = int(np.floor((Y1 - y_max)/dy))
        i1 = int(np.ceil((Y1 - y_min)/dy))

        if j1 <= 0 or i1 <= 0 or i0 >= n or j0 >= m:
            # Empty view window, raise an error since PIL does not
            raise ValueError("Selected view window is outside of image bounds")
        d = image.crop(box=(j0, i0, j1, i1))
        x0d, x1d = X0 + j0*dx, X0 + j1*dx
        y0d, y1d = Y1 - i1*dy, Y1 - i0*dy
        pv = triangulate_dem.rasterdata_to_pointvector(d, x0d, y0d, x1d, y1d, dx, dy)
        return pv

    def _extract_files(self, *,
                       x: float, y: float, dx: float, dy: float,
                       coordinate_system: pyproj.Proj) -> List[Path]:

        files = self.directory.glob("*.tif")
        path_list = []
        for file in files:
            with Image.open(file) as img:
                file_coordinate_system = GeoKeysInterpreter(extract_geo_keys(image=img)).to_proj4()
                l_proj = pyproj.Proj(file_coordinate_system)
                loc_x, loc_y = pyproj.transform(coordinate_system, l_proj, x, y)
                total_raster_bbox = geometry.Polygon.from_bounds(loc_x - dx,
                                                                 loc_y - dy,
                                                                 loc_x + dx,
                                                                 loc_y + dy)

                bounding_box = self._get_bounding_box(image=img)
                if bounding_box.polygon.intersects(total_raster_bbox):
                    path_list.append(file)
        return path_list

    def _get_bounding_box(self, *, image: TiffImagePlugin.TiffImageFile) -> GeoPolygon:
        m, n = image.size
        model_pixel_scale = image.tag_v2.get(GeoTiffTags.ModelPixelScaleTag.value)
        model_tie_point = image.tag_v2.get(GeoTiffTags.ModelTiePointTag.value)
        dx, dy, _ = model_pixel_scale
        jt, it, _, xt, yt, _ = model_tie_point
        # Determine image coordinates, assuming here that the tie point is associated with pixel center.
        # Otherwise, coordinates should be shiftet half pixel size in both directions
        X0, X1 = xt - jt * dx, xt + (n - 1 - jt) * dx
        Y1, Y0 = yt + it * dy, yt - (m - 1 - it) * dy
        polygon = geometry.Polygon.from_bounds(X0, Y0, X1, Y1)
        proj = pyproj.Proj(GeoKeysInterpreter(extract_geo_keys(image=image)).to_proj4())
        return GeoPolygon(polygon=polygon, proj=proj)


def read_sun_posisions(*, filepath: Path) -> triangulate_dem.ShadowVector:
    assert filepath.exists()
    with filepath.open("r") as ff:
        pass
    return triangulate_dem.ShadowVector()
