from typing import Optional, Dict, List
from owslib.wfs import WebFeatureService
from owslib.fes import PropertyIsLike
from owslib.etree import etree


class GeoPolygon:
    pass


class LandCoverType(Enum):
    pass


class InspireLandtypeRepository:

    def __init__(self):
        pass

    @property
    def land_cover_types(self, *, constraint: Optional[GeoPolygon] = None) -> List[LandCoverType]:
        pass

    def read(self, *, land_type: LandCoverType, constraint: GeoPolygon) -> List[GeoPolygon]:
        return self.read_all(land_types=[land_type], constraint=constraint)[land_type]

    def read_all(self, *, land_types: List[LandCoverType], constraint: GeoPolygon) -> Dict[LandCoverType, List[GeoPolygon]:
        pass

# Some test code using owslib. Have not yet figured out how to use the WebFeatureService to extract a subset of land surface types as
# polygons. Below is a simple test summing up my understanding up to now. See also:
# * https://www.linz.govt.nz/data/linz-data-service/guides-and-documentation/wfs-filtering-by-attribute-or-feature
# * https://geopython.github.io/OWSLib/
# * https://kartkatalog.geonorge.no/metadata/kartverket/inspire-landcovervector-wfs/c1f5f3de-ce2d-4104-bee3-2560c5a5a948

def test()
    wfs = WebFeatureService("https://wfs.geonorge.no/skwms1/wfs.inspire-lcv", version="2.0.0")
    filter = PropertyIsLike(propertyname="lcv:LandCoverObservation", literal="21")
    filterxml = etree.tostring(filter.toXML()).decode("utf-8")
    response = wfs.getfeature(typename="lcv:LandCoverUnit", filter=filterxml)
