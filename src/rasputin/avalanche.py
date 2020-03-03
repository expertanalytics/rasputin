from typing import List, Dict, Any, Union, Type, Optional, Tuple
from enum import Enum
import numpy as np
import requests
import json
import pyproj
import datetime
from shapely.geometry import Polygon
from .geometry import GeoPolygon

base_avalanche_url = "https://api01.nve.no/hydrology/forecast/avalanche/v4.0.2/api"


def get_forecasts(*, x: float, y: float, crs: pyproj.CRS) -> Dict[str, Any]:
    target_crs = pyproj.CRS.from_epsg(4326)
    proj = pyproj.Transformer.from_crs(crs, target_crs)
    lat, lon = proj.transform(x, y)
    lang = 1
    now = datetime.datetime.now()
    start = now.strftime("%Y-%m-%d")
    end = (now + datetime.timedelta(days=1)).strftime("%Y-%m-%d")

    query = f"/AvalancheWarningByCoordinates/Detail/{lat}/{lon}/{lang}/{start}/{end}"
    full = base_avalanche_url + query
    response = requests.get(full)
    if response.ok:
        res = json.loads(response.text)
        return res
    return {}


class Region:

    def __init__(self, *,
                 region_id: int,
                 region_name: str,
                 geo_polygon: GeoPolygon) -> None:
        self.region_id = region_id
        self.region_name = region_name
        self.geo_polygon = geo_polygon

    def __hash__(self) -> int:
        return hash((self.region_id, self.region_name, self.geo_polygon))

    def __eq__(self, other) -> bool:
        if not isinstance(other, self.__class__):
            return False
        if self.region_id != other.region_id:
            return False
        if self.region_name != other.region_name:
            return False
        if self.geo_polygon != other.geo_polygon:
            return False
        return True


class ARegion(Region):
    API_ID = 10


class BRegion(Region):
    API_ID = 20


class AvalancheType(Enum):
    not_given = 0
    slab = 10
    loose_snow = 20


class AvalancheProblemType(Enum):
    not_given = 0
    loose_dry = 3
    loose_wet = 5
    new_snow_slab = 7
    wind_slab = 10
    new_snow = 20
    persistent_slab = 30
    persistent_deep_slab = 37
    wet_snow = 40
    wet_slab = 45
    glide = 50


class AvalancheExt(Enum):
    cornice = 40
    slush = 30
    glide = 27
    wet_slab = 25
    dry_slab = 20
    loose_wet = 15
    loose_dry = 10
    not_given = 0


class AvalancheCause(Enum):
    loose_snow = 24
    water_pooling = 22
    wet_snow_or_melt_near_ground = 20
    buried_weak_layer_of_faceted_snow_beneath_crust = 19
    buried_weak_layer_of_faceted_snow_above_crust = 18
    buried_weak_layer_of_faceted_snow_near_ground = 16
    poor_bonding_between_layers_in_wind_deposited_snow = 15
    poor_bonding_between_crust_and_overlying_snow = 14
    buried_weak_layer_of_faceted_snow_near_surface = 13
    buried_weak_layer_of_surface_hoar = 11
    buried_weak_layer_of_new_snow = 10
    not_given = 0


class AvalancheProbability(Enum):
    very_likely = 7
    likely = 5
    possible = 3
    unlikely = 2
    not_given = 0


class AvalancheTrigger(Enum):
    very_easy = 60
    easy = 50
    difficult = 40
    very_difficult = 30
    spontaneous = 22
    low_additional_load = 21
    high_additional_load = 10
    not_given = 0


class AvalancheSize(Enum):
    not_given = 0
    small = 1
    medium = 2
    large = 3
    very_large = 4
    extreme = 5


class AvalanchePropagation(Enum):
    also_in_moderately_steep_terrain = 5
    most_steep_slopes = 4
    widespread_steep_sloped = 3
    specific_steep_slopes = 2
    isolated_steep_slopes = 1
    not_given = 0


class AvalancheDangerLevel(Enum):
    not_given = 0
    low = 1
    moderate = 2
    considerable = 3
    high = 4
    extereme = 5


class AvalancheDanger:

    varsom_angles = [( -22.5,  22.5),
                     (  22.5,  67.5),
                     (  67.5, 112.5),
                     ( 112.5, 157.5),
                     ( 157.5,-157.5),
                     (-157.5,-112.5),
                     (-112.5, -67.5),
                     ( -67.5, -22.5)]

    def __init__(self, *,
                 expositions: str,
                 height1: int,
                 height2: int,
                 height_fill: int,
                 avy_type: AvalancheType,
                 avy_problem: AvalancheProblemType,
                 avy_ext: AvalancheExt,
                 avy_cause: AvalancheCause,
                 avy_prob: AvalancheProbability,
                 avy_trigger: AvalancheTrigger,
                 avy_size: AvalancheSize,
                 avy_prop: AvalanchePropagation,
                 avy_danger: AvalancheDangerLevel) -> None:
        expositions = [bool(int(s)) for s in expositions]
        angles = np.asarray(self.varsom_angles)/180*np.pi
        exposed_angles = angles[expositions].tolist()
        self.exposed_angles = exposed_angles
        self.h1 = height1
        self.h2 = height2
        self.hf = height_fill
        self.avy_type = avy_type
        self.avy_problem = avy_problem
        self.avy_ext = avy_ext
        self.avy_cause = avy_cause
        self.prob = avy_prob
        self.avy_trigger = avy_trigger
        self.avy_size = avy_size
        self.avy_prop = avy_prop
        self.avy_danger = avy_danger

    @property
    def danger_intervals(self) -> List[Tuple[float, float]]:
        if self.hf == 1:
            return [(self.h1, 20000)]
        elif self.hf == 2:
            return [(0, self.h1)]
        elif self.hf == 3:
            return [(0, self.h1), (self.h2, 20000)]
        elif self.hf == 4:
            return [(self.h1, self.h2)]


def build_regions(*,
                  region_type: Union[Type[ARegion], Type[BRegion]],
                  crs: pyproj.CRS,
                  domain: Optional[GeoPolygon] = None
                  ) -> Union[List[ARegion], List[BRegion]]:
    """Build avalanche region lists.

    :param region_type: Either ARegion or BRegion, corresponding to the different avalanche area types
                        in Norway. See  http://api.nve.no/doc/snoeskredvarsel/#regionsummary
                        Short story: A Regions are published daily from Dec 1 to May 31, and B Regions only
                        when avalanche dangers are either 4 or 5.
    :param crs:         Coordinate system for the region polygons.
    :param domain:      If given, only return the avalanche regions limited by this domain.
    :return:            List of regions, possibly intersected by the domain polygon.
    """
    query =  f"/Region/{region_type.API_ID}"
    full = base_avalanche_url + query
    response = requests.get(full)
    if not response.ok:
        raise requests.HTTPError("Could not read regions from Varsom.no")
    regions_info = json.loads(response.text)
    data_crs = pyproj.CRS.from_epsg(4326)
    region_list = []
    for region_info in regions_info:
        fpoly = region_info["Polygon"]
        assert len(fpoly) == 1
        fpoly = fpoly[0]
        polygon = Polygon([tuple(float(x) for x in coord.split(",")) for coord in fpoly.split()])
        geo_polygon = GeoPolygon(polygon=polygon, crs=data_crs)
        if domain is not None:
            if domain.intersects(geo_polygon):
                region_list.append(region_type(region_id=region_info["Id"],
                                               region_name=region_info["Name"],
                                               geo_polygon=geo_polygon.transform(target_crs=crs).intersection(domain)))
        else:
            region_list.append(region_type(region_id=region_info["Id"],
                                           region_name=region_info["Name"],
                                           geo_polygon=geo_polygon.transform(target_crs=crs)))
    return region_list

def get_avalanche_warnings(*,
                           start_date: datetime.datetime,
                           end_date: datetime.datetime,
                           regions: List[Region]) -> Dict[Region, List[AvalancheDanger]]:
    result = {}
    start = start_date.strftime("%Y-%m-%d")
    end = end_date.strftime("%Y-%m-%d")
    query = "/AvalancheWarningByRegion/Detail/{RegionId}/{LanguageKey}/" + f"{start}/{end}"
    for region in regions:
        full = base_avalanche_url + query.format(RegionId=region.region_id, LanguageKey=2)
        response = requests.get(full)
        if not response.ok:
            raise requests.HTTPError(f"Could not fetch avalanche warning for region {region.region_name}")
        region_warnings = json.loads(response.text)
        for warning in region_warnings:
            if warning["AvalancheProblems"] is None:
                continue
            for problem in warning["AvalancheProblems"]:
                ad = AvalancheDanger(
                    expositions=problem["ValidExpositions"],
                    height1=problem["ExposedHeight1"],
                    height2=problem["ExposedHeight2"],
                    height_fill=problem["ExposedHeightFill"],
                    avy_type=AvalancheType(problem["AvalancheTypeId"]),
                    avy_problem=AvalancheProblemType(problem["AvalancheProblemTypeId"]),
                    avy_ext=AvalancheExt(problem["AvalancheExtId"]),
                    avy_cause=AvalancheCause(problem["AvalCauseId"]),
                    avy_prob=AvalancheProbability(problem["AvalProbabilityId"]),
                    avy_trigger=AvalancheTrigger(problem["AvalTriggerSimpleId"]),
                    avy_size=AvalancheSize(problem["DestructiveSizeExtId"]),
                    avy_prop=AvalanchePropagation(problem["AvalPropagationId"]),
                    avy_danger=AvalancheDangerLevel(int(warning["DangerLevel"]))
                )
                if region not in result:
                    result[region] = []
                result[region].append(ad)
    return result
