from typing import List, Optional, Tuple, Dict, Any
from pathlib import Path
import numpy as np
from pyproj import CRS, Transformer
from PIL import Image
from rasputin.reader import extract_geo_keys, GeoKeysInterpreter, GeoTiffTags
from rasputin.land_cover_repository import LandCoverBaseType, LandCoverMetaInfoBase, LandCoverRepository
import rasputin.triangulate_dem as td
from rasputin.geometry import GeoPoints, GeoPolygon
from rasputin.material_specification import lake_material, terrain_material

Image.MAX_IMAGE_PIXELS = None


class LandCoverType(LandCoverBaseType):
    crop_type_1 = 11
    crop_type_2 = 14
    crop_type_3 = 20
    crop_type_4 = 30
    forest_type_1 = 40
    forest_type_2 = 50
    forest_type_3 = 60
    forest_type_4 = 70
    forest_type_5 = 90
    forest_type_6 = 100
    shrub_type_1 = 110
    shrub_type_2 = 120
    shrub_type_3 = 130
    vegetation_type_1 = 140
    vegetation_type_2 = 150
    flood_type_1 = 160
    flood_type_2 = 170
    flood_type_3 = 180
    artificial = 190
    bare = 200
    water = 210
    snow_and_ice = 220
    no_data = 230


class LandCoverMetaInfo(LandCoverMetaInfoBase):

    @classmethod
    def describe(cls, *, land_cover_type=LandCoverType) -> str:
        description = {
            LandCoverType.crop_type_1: "Post-flooding or irrigated croplands (or aquatic)",
            LandCoverType.crop_type_2: "Rainfed croplands",
            LandCoverType.crop_type_3: "Mosaic cropland (50-70%) / vegetation (grassland/shrubland/forest) (20-50%)",
            LandCoverType.crop_type_4: "Mosaic vegetation (grassland/shrubland/forest) (50-70%) / cropland (20-50%)",
            LandCoverType.forest_type_1: "Closed to open (>15%) broadleaved evergreen or semi-deciduous forest (>5m)",
            LandCoverType.forest_type_2: "Closed (>40%) broadleaved deciduous forest (>5m)",
            LandCoverType.forest_type_3: "Open (15-40%) broadleaved deciduous forest/woodland (>5m)",
            LandCoverType.forest_type_4: "Closed (>40%) needleleaved evergreen forest (>5m)",
            LandCoverType.forest_type_5: "Open (15-40%) needleleaved deciduous or evergreen forest (>5m)",
            LandCoverType.forest_type_6: "Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)",
            LandCoverType.shrub_type_1: "Mosaic forest or shrubland (50-70%) / grassland (20-50%)",
            LandCoverType.shrub_type_2: "Mosaic grassland (50-70%) / forest or shrubland (20-50%)",
            LandCoverType.shrub_type_3: "Closed to open (>15%) (broadleaved or needleleaved, evergreen or deciduous) shrubland (<5m)",
            LandCoverType.vegetation_type_1: "Closed to open (>15%) herbaceous vegetation (grassland, savannas or lichens/mosses)",
            LandCoverType.vegetation_type_2: "Sparse (<15%) vegetation",
            LandCoverType.flood_type_1: "Closed to open (>15%) broadleaved forest regularly flooded (semi-permanently or temporarily) - Fresh or brackish water",
            LandCoverType.flood_type_2: "Closed (>40%) broadleaved forest or shrubland permanently flooded - Saline or brackish water",
            LandCoverType.flood_type_3: "Closed to open (>15%) grassland or woody vegetation on regularly flooded or waterlogged soil - Fresh, brackish or saline water",
            LandCoverType.artificial: "Artificial surfaces and associated areas (Urban areas >50%)",
            LandCoverType.bare: "Bare areas",
            LandCoverType.water: "Water bodies",
            LandCoverType.snow_and_ice: "Permanent snow and ice",
            LandCoverType.no_data: "No data (burnt areas, clouds,…)"}
        return description[land_cover_type]

    @classmethod
    def material(cls, *, land_cover_type: LandCoverBaseType) -> Dict[str, Any]:
        if land_cover_type == LandCoverType.water:
            return lake_material
        return terrain_material

    @classmethod
    def color(cls, *, land_cover_type: LandCoverType) -> Tuple[int, int, int]:
        """Default colors for types in dataset."""
        colors = {
            LandCoverType.crop_type_1: (170, 240, 240),
            LandCoverType.crop_type_2: (255, 255, 100),
            LandCoverType.crop_type_3: (220, 240, 100),
            LandCoverType.crop_type_4: (205, 205, 102),
            LandCoverType.forest_type_1: (0, 100, 0),
            LandCoverType.forest_type_2: (0, 160, 0),
            LandCoverType.forest_type_3: (170, 200, 0),
            LandCoverType.forest_type_4: (0, 60, 0),
            LandCoverType.forest_type_5: (40, 100, 0),
            LandCoverType.forest_type_6: (120, 130, 0),
            LandCoverType.shrub_type_1: (140, 160, 0),
            LandCoverType.shrub_type_2: (190, 150, 0),
            LandCoverType.shrub_type_3: (150, 100, 0),
            LandCoverType.vegetation_type_1: (255, 180, 50),
            LandCoverType.vegetation_type_2: (255, 235, 175),
            LandCoverType.flood_type_1: (0, 120, 90),
            LandCoverType.flood_type_2: (0, 150, 120),
            LandCoverType.flood_type_3: (0, 220, 130),
            LandCoverType.artificial: (195, 20, 0),
            LandCoverType.bare: (255, 245, 215),
            LandCoverType.water: (0, 70, 200),
            LandCoverType.snow_and_ice: (255, 255, 255),
            LandCoverType.no_data: (0, 0, 0)
        }
        return colors[land_cover_type]


class GlobCovRepository(LandCoverRepository):

    land_cover_type = LandCoverType
    land_cover_meta_info_type = LandCoverMetaInfo

    def __init__(self, *, path: Path) -> None:
        super().__init__()
        self.path = path / "GLOBCOVER_L4_200901_200912_V2.3.tif"
        assert self.path.is_file()
        with Image.open(self.path) as image:
            self.geo_keys = extract_geo_keys(image=image)
            self.gk_interpreter = GeoKeysInterpreter(self.geo_keys)
            self.data_crs = CRS.from_string("+init=EPSG:4326")
            self.model_tie_point = image.tag_v2.get(GeoTiffTags.ModelTiePointTag.value)
            self.model_pixel_scale = image.tag_v2.get(GeoTiffTags.ModelPixelScaleTag.value)
            self.M, self.N = image.size

    def read(self,
             *,
             land_type: LandCoverType,
             geo_points: GeoPoints,
             domain: GeoPolygon):
        return self.land_cover(land_types=[land_type], geo_points=geo_points, domain=domain)

    def constraints(self, *, domain: GeoPolygon) -> List[GeoPolygon]:
        return []

    def land_cover(self,
                   *,
                   land_types: Optional[List[LandCoverType]],
                   geo_points: GeoPoints,
                   domain: GeoPolygon) -> np.ndarray:
        pts_crs = geo_points.crs
        if pts_crs != self.data_crs:
            proj = Transformer.from_crs(pts_crs, self.data_crs)
            xy = np.dstack(proj.transform(geo_points.xy[:, 0], geo_points.xy[:, 1])).reshape(-1, 2)
        else:
            xy = geo_points.xy
        with Image.open(self.path) as image:
            return self._extract_land_types(image=image,
                                            land_types=land_types,
                                            indices=self._raster_indices(xy=xy))

    def _raster_indices(self, *, xy: np.ndarray) -> td.index_vector:
        pts = td.point2_vector(xy.tolist())
        dx, dy, _ = self.model_pixel_scale
        jt, it, _, xt, yt, _ = self.model_tie_point
        X0 = xt - jt*dx
        Y1 = yt + it * dy
        return td.coordinates_to_indices(X0, Y1, dx, dy, self.M, self.N, pts)

    def _extract_land_types(self,
                            *,
                            image,
                            land_types: Optional[List[LandCoverType]],
                            indices: td.index_vector) -> np.ndarray:
        # TODO: This does not work, perhaps figure out why?
        #all_land_types = np.asarray(td.extract_uint8_buffer_values(indices, image))
        all_land_types = np.asarray([image.getpixel(tuple(idx)) for idx in indices])
        if land_types is None:
            return all_land_types
        lcts = [lct.value for lct in land_types]
        func = np.vectorize(lambda t: t not in lcts)
        all_land_types[func(all_land_types)] = 0
        return all_land_types

