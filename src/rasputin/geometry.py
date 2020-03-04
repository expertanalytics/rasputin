from typing import Tuple, List, Optional, Union
import io
import shutil
from shapely import ops, wkt, wkb
from pathlib import Path
from dataclasses import dataclass
from shapely import geometry
import numpy as np
from pkg_resources import resource_filename
import pyproj
from rasputin import triangulate_dem as td
from rasputin.py2js import point_vector_to_lines, face_and_point_vector_to_lines
from rasputin.mesh_utils import vertex_field_to_vertex_values, face_field_to_vertex_values
from rasputin.material import terrain_material


class Geometry:

    def __init__(self, *,
                 mesh: "Mesh",
                 crs: pyproj.CRS,
                 base_color: Optional[Tuple[float, float, float]] = None,
                 material: Optional[str] = None):
        self.mesh = mesh
        self.crs = crs
        if base_color is None:
            self.base_color = (1.0, 1.0, 1.0)
        else:
            self.base_color = base_color
        self._surface_normals = None
        self._point_normals = None
        self._aspects = None
        self._slopes = None
        self._colors = None
        self._uvs = None
        self.texture = np.empty((1024,1024))
        self.scale_x = 1
        self.scale_y = 1
        self._material = material or terrain_material

    def extract_faces(self, faces: np.ndarray) -> "Geometry":
        return self.__class__(mesh=self.mesh.copy(),
                              crs=self.crs,
                              base_color=self.base_color,
                              material=self.material)

    def split_by_colors(self) -> List["Geometry"]:
        print(len(self.colors))
        print(len(self._colors))
        print(len(self.faces))
        sub_geometries = []
        if self._colors is None:
            raise RuntimeError("No colors to split by")
        for color in set([tuple(c) for c in self._colors]):
            color = np.array(color)
            sub_geom = self.extract_faces(faces=self.faces[np.all(self._colors==color, axis=1)])
            sub_geom.base_color = color
            sub_geometries.append(sub_geom)
        return sub_geometries


    @property
    def points(self) -> np.ndarray:
        return self.mesh.points

    @property
    def faces(self) -> np.ndarray:
        return self.mesh.faces

    @property
    def point_normals(self) -> np.ndarray:
        return self.mesh.point_normals

    @property
    def face_normals(self) -> np.ndarray:
        return self.mesh.face_normals

    @property
    def aspects(self) -> np.ndarray:
        if self._aspects is None:
            self._aspects = np.asarray(td.compute_aspects(self.face_normals))
        return self._aspects

    @property
    def slopes(self) -> np.ndarray:
        if self._slopes is None:
            self._slopes = np.asarray(td.compute_slopes(self.face_normals))
        return self._slopes

    @property
    def colors(self) -> np.ndarray:
        if self._colors is None:
            self._colors = np.empty((len(self.mesh.faces), len(self.base_color)))
            self._colors[:, 0] = self.base_color[0]
            self._colors[:, 1] = self.base_color[1]
            self._colors[:, 2] = self.base_color[2]
        return self._colors

    @colors.setter
    def colors(self, colors: np.ndarray) -> None:
        assert colors.shape == (len(self.mesh.faces), 3)
        self._colors = colors

    @property
    def color_lines(self):
        return face_field_to_vertex_values(face_field=self.colors, faces=self.faces)

    @property
    def material(self) -> str:
        return self._material

    @material.setter
    def material(self, material: str):
        self._material = material

    def as_javascript(self, *, file_handle: io.TextIOWrapper) -> True:
        if not len(self.faces):
            return False
        file_handle.write("{")
        # write vertices
        file_handle.write(f"vertices: ")
        for line in face_and_point_vector_to_lines(name=None,
                                                   face_vector=self.mesh.faces,
                                                   point_vector=self.mesh.points):
            file_handle.write(f"{line}\n")
        file_handle.write(",\n")

        file_handle.write("normals: ")
        point_normals = vertex_field_to_vertex_values(vertex_field=np.asarray(self.point_normals),
                                                      faces=self.faces,
                                                      points=self.points)
        for line in point_vector_to_lines(name=None, point_vector=point_normals):
            file_handle.write(f"{line}\n")
        file_handle.write(",\n")

        # write colors
        file_handle.write("colors: ")
        for line in point_vector_to_lines(name=None, point_vector=self.color_lines):
            file_handle.write(f"{line}\n")
        file_handle.write(",\n")

        file_handle.write(f"uvs: ")
        for line in face_and_point_vector_to_lines(name=None,
                                                   face_vector=self.mesh.faces,
                                                   point_vector=self.uvs):
            file_handle.write(f"{line}\n")
        file_handle.write(",\n")

        # write material
        file_handle.write(f"material_constructor: {self.material}")
        file_handle.write("}")
        return True
    
    @property
    def uvs(self) -> np.ndarray:
        w, h = self.texture.shape
        x = self.points[:, 0]
        y = self.points[:, 1]
        
        u = (x*self.scale_x)/float(w)
        v = (y*self.scale_y)/float(h)
        return np.array([u, v]).T

def write_scene(*, geometries: List[Geometry], output: Path):
    if output.exists():
        shutil.rmtree(output)
    templates_path = Path(resource_filename(__name__, "web"))
    shutil.copytree(templates_path, output)
    data_js = output / "data.js"
    with data_js.open("w") as tf:
        tf.write("const geometries = [\n")
        for geometry in geometries:
            if geometry.as_javascript(file_handle=tf):
                tf.write(",\n")
        tf.write("];\n")
        tf.write("const data = {geometries};\n")


class GeoPoints:

    def __init__(self, *, xy: np.ndarray, crs: pyproj.CRS) -> None:
        self.xy = xy
        assert len(self.xy.shape) == 2 and self.xy.shape[-1] == 2
        self.crs = crs


@dataclass
class GeoPoint:

    point: geometry.Point
    crs: pyproj.CRS


class GeoPolygon:

    def __init__(self, *, polygon: geometry.Polygon, crs: pyproj.CRS) -> None:
        self.polygon = polygon
        self.crs = crs

    @classmethod
    def from_raster_file(cls, *, filepath: Path) -> "GeoPolygon":
        from rasputin.reader import identify_projection, get_image_extents
        from PIL import Image
        with Image.open(filepath) as image:
            projection_str = identify_projection(image=image)
            image_extents = get_image_extents(image=image)

        return GeoPolygon(polygon=geometry.Polygon.from_bounds(*image_extents.box),
                          crs=pyproj.CRS.from_proj4(projection_str))

    @classmethod
    def from_polygon_file(cls, *, filepath: Path, crs: pyproj.CRS) -> "GeoPolygon":
        if filepath.suffix.lower() == ".wkb":
            with filepath.open("rb") as pfile:
                polygon = wkb.loads(pfile.read())
        elif filepath.suffix.lower() == ".wkt":
            with filepath.open("r") as pfile:
                polygon = wkt.loads(pfile.read())
        else:
            raise ValueError("Cannot determine file type.")

        return GeoPolygon(polygon=polygon, crs=crs)

    def transform(self, *, target_crs: pyproj.CRS) -> "GeoPolygon":
        if target_crs.to_authority() != self.crs.to_authority():
            proj = pyproj.Transformer.from_crs(self.crs, target_crs)
            polygon = ops.transform(proj.transform, self.polygon)
        else:
            polygon = self.polygon
        return GeoPolygon(polygon=polygon, crs=target_crs)

    def intersection(self, other: "GeoPolygon") -> "GeoPolygon":
        other = other.transform(target_crs=self.crs)
        return GeoPolygon(polygon=self.polygon.intersection(other.polygon),
                          crs=self.crs)

    def intersects(self, other: "GeoPolygon") -> bool:
        op = other.transform(target_crs=self.crs).polygon
        sp = self.polygon
        return sp.intersects(op) and not sp.touches(op)

    def difference(self, other: "GeoPolygon") -> "GeoPolygon":
        other = other.transform(target_crs=self.crs)
        return GeoPolygon(polygon=self.polygon.difference(other.polygon),
                          crs=self.crs)

    def buffer(self, value: float) -> "GeoPolygon":
        return GeoPolygon(polygon=self.polygon.buffer(value), crs=self.crs)

    def to_cpp(self) -> Union[td.simple_polygon, td.polygon]:
        # Note: Shapely polygons have one vertex repeated and orientation dos not matter
        #       CGAL polygons have no vertex repeated and orientation mattters
        from shapely.geometry.polygon import orient

        # Get point sequences as arrays
        exterior = np.asarray(orient(self.polygon).exterior)[:-1]
        interiors = [np.asarray(hole)[:-1]
                     for hole in orient(self.polygon, -1).interiors]

        cgal_polygon = td.simple_polygon(exterior)
        for interior in interiors:
            cgal_polygon = cgal_polygon.difference(td.simple_polygon(interior))
            # Intersection result is MultiPolygon
            if cgal_polygon.num_parts() == 1:
                cgal_polygon = cgal_polygon[0]
        return cgal_polygon
