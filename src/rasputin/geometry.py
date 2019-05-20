from typing import Tuple, List, Optional
import io
import shutil
from pathlib import Path
import numpy as np
from pkg_resources import resource_filename
from rasputin import triangulate_dem as td
from rasputin.py2js import point_vector_to_lines, face_and_point_vector_to_lines
from rasputin.mesh_utils import vertex_field_to_vertex_values


lake_material = """\
THREE.MeshPhongMaterial( {
    specular: 0xffffff,
    shininess: 25,
    color: 0x006994,
    reflectivity: 0.3,
} )
"""

avalanche_material = """\
THREE.MeshPhongMaterial( {
    specular: 0xffffff,
    shininess: 75,
    reflectivity: 0.3,
    vertexColors: THREE.FaceColors
} )
"""

terrain_material = """\
THREE.MeshPhysicalMaterial( {
    metalness: 0.0,
    roughness: 0.5,
    reflectivity: 0.7,
    clearCoat: 0.0,
    vertexColors: THREE.FaceColors
} )
"""


class Geometry:

    def __init__(self, *,
                 points: td.point3_vector,
                 faces: td.face_vector,
                 projection: str,
                 base_color: Optional[Tuple[float, float, float]],
                 material: Optional[str]):
        self.points = points
        self.faces = faces
        self.projection = projection
        if base_color is None:
            self.base_color = (1.0, 1.0, 1.0)
        else:
            self.base_color = base_color
        self._surface_normals = None
        self._point_normals = None
        self._aspects = None
        self._slopes = None
        self._colors = None
        self._material = material or terrain_material

    def consolidate(self) -> "Geometry":
        consolidated_points, consolidated_faces = td.consolidate(self.points, self.faces)
        return Geometry(points=consolidated_points,
                        faces=consolidated_faces,
                        projection=self.projection,
                        base_color=self.base_color,
                        material=self.material)

    @property
    def point_normals(self) -> td.point3_vector:
        if self._point_normals is None:
            self._point_normals = td.point_normals(self.points, self.faces)
        return self._point_normals

    @property
    def surface_normals(self) -> td.point3_vector:
        if self._surface_normals is None:
            self._surface_normals = td.surface_normals(self.points, self.faces)
        return self._surface_normals

    @property
    def aspects(self) -> np.ndarray:
        if self._aspects is None:
            self._aspects = np.asarray(td.compute_aspect(self.surface_normals))
        return self._aspects

    @property
    def slopes(self) -> np.ndarray:
        if self._slopes is None:
            self._slopes = np.asarray(td.compute_slopes(self.surface_normals))
        return self._slopes

    @property
    def colors(self) -> np.ndarray:
        if self._colors is None:
            self._colors = np.empty((3*len(self.faces), len(self.base_color)))
            self._colors[:, 0] = self.base_color[0]
            self._colors[:, 1] = self.base_color[1]
            self._colors[:, 2] = self.base_color[2]
        return self._colors

    @property
    def material(self) -> str:
        return f"new {self._material}"

    @material.setter
    def material(self, material: str):
        self._material = material

    def as_javascript(self, *, file_handle: io.TextIOWrapper) -> True:
        if not self.faces:
            return False
        file_handle.write("{")
        # write vertices
        file_handle.write(f"vertices: ")
        for line in face_and_point_vector_to_lines(name=None,
                                                   face_vector=self.faces,
                                                   point_vector=self.points):
            file_handle.write(f"{line}\n")
        file_handle.write(",\n")

        #write normals
        file_handle.write("normals: ")
        point_normals = vertex_field_to_vertex_values(vertex_field=np.asarray(self.point_normals),
                                                      faces=self.faces,
                                                      points=self.points)
        for line in point_vector_to_lines(name=None, point_vector=point_normals):
            file_handle.write(f"{line}\n")
        file_handle.write(",\n")

        # write colors
        file_handle.write("colors: ")
        for line in point_vector_to_lines(name=None, point_vector=self.colors):
            file_handle.write(f"{line}\n")
        file_handle.write(",\n")

        # write material
        file_handle.write(f"material: {self.material}\n")
        file_handle.write("}")
        return True


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


    pass
