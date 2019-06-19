from typing import Dict, List, Tuple, Optional
from pathlib import Path
from lxml import etree
import numpy as np
from pyproj import transform
from pyproj import Proj
from shapely.geometry import Polygon, Point
from rasputin.geometry import GeoPoints, GeoPolygon
from rasputin.land_cover_repository import LandCoverBaseType, LandCoverRepository, LandCoverMetaInfoBase
from rasputin.material import lake_material, terrain_material


class LandCoverType(LandCoverBaseType):
    # Artificial
    urban_fabric_cont = 111
    urban_fabrid_discont = 112
    industrial_unit = 121
    road_and_rail = 122
    port = 123
    airport = 124
    mineral_extraction = 131
    dump_site = 132
    constrution_site = 133
    urban_green = 141
    sport_and_leisure = 142

    # Agricultural
    arable_land_non_irr = 211
    permanent_irr = 212
    rice_field = 213
    vinyard = 221
    fruit_and_berry = 222
    olive_grove = 223
    pasture = 231
    mix_annual_permament_crop = 241
    complex_cultivation = 242
    mix_agri_natural = 243
    agro_forestry = 244

    # Forest
    broad_leaved = 311
    coniferous = 312
    mixed_forest = 313
    natural_grass = 321
    moors_and_heath = 322
    sclerophyllous = 323
    transitional_woodland_shrub = 324
    beach_dune_sand = 331
    bare_rock = 332
    sparse_veg = 333
    burnt = 334
    glacier_and_snow = 335

    # Wetland
    inland_march = 411
    peat_bog = 412
    salt_march = 421
    saline = 422
    intertidal_flat = 423

    # Water bodies
    water_course = 511
    water_body = 512
    coastal_lagoon = 521
    estuary = 522
    sea_and_ocean = 523


class LandCoverMetaInfo(LandCoverMetaInfoBase):

    @classmethod
    def color(cls, *, land_cover_type: LandCoverType) -> Tuple[int, int, int]:

        return {111: (230,   0,  77),
                112: (255,   0,   0),
                121: (204,  77, 242),
                122: (204,   0,   0),
                123: (230, 204, 204),
                124: (230, 204, 230),
                131: (166,   0, 204),
                132: (166,  77,   0),
                133: (255,  77, 255),
                141: (255, 166, 255),
                142: (255, 230, 255),
                211: (255, 255, 168),
                212: (255, 255,   0),
                213: (230, 230,   0),
                221: (230, 128,   0),
                222: (242, 166,  77),
                223: (230, 166,   0),
                231: (230, 230,  77),
                241: (255, 230, 166),
                242: (255, 230,  77),
                243: (230, 204,  77),
                244: (242, 204, 166),
                311: (128, 255,   0),
                312: (  0, 166,   0),
                313: ( 77, 255,   0),
                321: (204, 242,  77),
                322: (166, 255, 128),
                323: (166, 230,  77),
                324: (166, 242,   0),
                331: (230, 230, 230),
                332: (204, 204, 204),
                333: (204, 255, 204),
                334: (  0,   0,   0),
                335: (166, 230, 204),
                411: (166, 166, 255),
                412: ( 77,  77, 255),
                421: (204, 204, 255),
                422: (230, 230, 255),
                423: (166, 166, 230),
                511: (  0, 204, 242),
                512: (128, 242, 230),
                521: (  0, 255, 166),
                522: (166, 255, 230),
                523: (230, 242, 255)}[land_cover_type.value]

    @classmethod
    def describe(cls, *, land_cover_type=LandCoverType) -> str:
        return ""

    @classmethod
    def material(cls, *, land_cover_type=LandCoverType) -> str:
        if land_cover_type.value > 500:
            return lake_material
        return terrain_material


class GMLRepository(LandCoverRepository):

    land_cover_type = LandCoverType
    land_cover_meta_info_type = LandCoverMetaInfo

    def __init__(self, path: Path) -> None:
        self.path = path
        files = list(path.glob("*.gml"))
        assert len(files) == 1
        self.fn = files[0]
        self.parser = etree.XMLParser(encoding="utf-8", recover=True, huge_tree=True)
        self.source_proj = None

    def _parse_polygon(self, polygon: etree._Element, nsmap: Dict[str, str]) -> Polygon:
        gml = f"{{{nsmap['gml']}}}"
        ogr = f"{{{nsmap['ogr']}}}"
        outer_boundary = next(polygon.iter(f"{gml}outerBoundaryIs"))
        citer = (float(x) for x in next(outer_boundary.iter(f"{gml}coordinates")).text.replace(",", " ").split())
        shell = np.fromiter(citer, dtype='d').reshape(-1, 2)
        holes = []
        for inner_boundary in polygon.iter(f"{gml}innerBoundaryIs"):
            citer = (float(x) for x in next(inner_boundary.iter(f"{gml}coordinates")).text.replace(",", " ").split())
            holes.append(np.fromiter(citer, dtype='d').reshape(-1, 2))
        return Polygon(shell=shell, holes=holes)

    def _assign_source_proj(self, root: etree._Element) -> None:
        gml = f"{{{root.nsmap['gml']}}}"
        elm = next(root.iter(f"{gml}featureMember"))
        self.source_proj = Proj(next(elm.iter(f"{gml}Polygon")).attrib["srsName"])

    def read(self, domain: GeoPolygon) -> Dict[LandCoverType, List[Polygon]]:
        tree = etree.parse(str(self.fn), self.parser)
        root = tree.getroot()
        self._assign_source_proj(root)
        assert "gml" in root.nsmap, "Not a GML file!"
        assert "ogr" in root.nsmap, "Can not find land cover types!"
        gml = f"{{{root.nsmap['gml']}}}"
        ogr = f"{{{root.nsmap['ogr']}}}"
        result = {}
        for elm in root.iter(f"{gml}featureMember"):
            code = LandCoverType(int(next(elm.iter(f"{ogr}clc18_kode")).text))
            polygon = self._parse_polygon(next(elm.iter(f"{gml}Polygon")), root.nsmap)
            if polygon.intersects(domain.polygon):
                if code not in result:
                    result[code] = []
                result[code].append(polygon.intersection(domain.polygon))
        return result

    def constraints(self, *, domain: GeoPolygon) -> List[GeoPolygon]:
        return []

    def land_cover(self,
                   *,
                   land_types: Optional[List[LandCoverType]],
                   geo_points: GeoPoints,
                   domain: GeoPolygon) -> np.ndarray:
        tree = etree.parse(str(self.fn), self.parser)
        root = tree.getroot()
        self._assign_source_proj(root)
        target_proj = geo_points.projection
        if target_proj.definition_string() != self.source_proj.definition_string():
            xy = np.dstack(transform(target_proj, self.source_proj, geo_points.xy[:, 0], geo_points.xy[:, 1]))[0]
        else:
            xy = geo_points.xy
        domain = domain.transform(target_projection=self.source_proj).buffer(50)

        # At this point, we have a consistent coordinate system for domain and xy all in the self.source_proj, so
        # omit checks for speed.

        land_cover_types = self.read(domain)

        # TODO: Move to C++ for speed!
        faces = np.zeros(len(xy), dtype="int")
        for i, pt in enumerate(xy):
            pt = Point(pt)
            found = False
            for key, p_list in land_cover_types.items():
                for p in p_list:
                    if p.intersects(pt):
                        faces[i] = key.value
                        found = True
                        break
                if found:
                    break
        return faces

