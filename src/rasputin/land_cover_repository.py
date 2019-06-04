from typing import List, Optional
from abc import ABC, abstractmethod
from enum import Enum
import numpy as np
from rasputin.reader import GeoPolygon
from rasputin.geometry import GeoPoints


class LandCoverTypeBase(ABC, Enum):
    pass


class LandCoverRepository(ABC):

    @abstractmethod
    def constraints(self, *, domain: GeoPolygon) -> List[GeoPolygon]:
        pass

    @abstractmethod
    def land_cover(self, *,
                   land_cover_types: Optional[List[LandCoverTypeBase]],
                   geo_points: GeoPoints,
                   domain: GeoPolygon) -> np.ndarray:
        pass