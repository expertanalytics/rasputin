import numpy as np
from pyproj import Proj
from pathlib import Path
from enum import Enum


class GeoPoints:

    def __init__(self, *, xy: np.ndarray, projection: Proj) -> None:
        self.xy = xy
        self.projection = projection


class LandCoverType(Enum):
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

    @staticmethod
    def describe(*, type="LandCoverType") -> str:
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
        return description[type]


class GlobCovRepository:

    def __init__(self, *, path: Path) -> None:
        self.path = path
        assert self.path.is_dir()

    def read(self, *, land_type: LandCoverType, geo_points: GeoPoints) -> np.ndarray:
        pass
