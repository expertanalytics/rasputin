from typing import List, Optional, Tuple, Dict, Any
from abc import ABC, abstractmethod
from enum import Enum
import numpy as np
from rasputin.geometry import GeoPoints, GeoPolygon


class LandCoverBaseType(Enum):
    pass


class LandCoverMetaInfoBase(ABC):

    @classmethod
    @abstractmethod
    def color(cls, *, land_cover_type: LandCoverBaseType) -> Tuple[int, int, int]:
        pass

    @classmethod
    @abstractmethod
    def describe(cls, *, land_cover_type: LandCoverBaseType) -> str:
        pass

    @classmethod
    @abstractmethod
    def material(cls, *, land_cover_type: LandCoverBaseType) -> Dict[str, Any]:
        pass

class LandCoverRepository(ABC):

    @property
    @abstractmethod
    def land_cover_type(self) -> LandCoverBaseType:
        pass

    @property
    @abstractmethod
    def land_cover_meta_info_type(self) -> LandCoverMetaInfoBase:
        pass

    @abstractmethod
    def constraints(self, *, domain: GeoPolygon) -> List[GeoPolygon]:
        pass

    @abstractmethod
    def land_cover(self,
                   *,
                   land_cover_types: Optional[List[LandCoverBaseType]],
                   geo_points: GeoPoints,
                   domain: GeoPolygon) -> np.ndarray:
        pass

