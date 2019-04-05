from typing import Tuple
import numpy as np
from rasputin import triangulate_dem as td


class Geometry:

    def __init__(self, *,
                 points: td.PointVector,
                 faces: td.FaceVector,
                 base_color: Tuple[float, float, float]):
        self.points = points
        self.faces = faces
        self._base_color = base_color
        self._surface_normals = None
        self._point_normals = None
        self._aspects = None
        self._slopes = None
        self._colors = None
        self._material = None

    @property
    def point_normals(self) -> td.PointVector:
        if self._point_normals is None:
            self._point_normals = td.point_normals(self.points, self.faces)
        return self._point_normals

    @property
    def surface_normals(self) -> td.PointVector:
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
            self._colors = np.empty((len(self.points), len(self._base_color)))
            self.colors[:, 0] = self._base_color[0]
            self.colors[:, 1] = self._base_color[1]
            self.colors[:, 2] = self._base_color[2]
        return self.colors

    @property
    def material(self) -> str:
        """TODO: Make a system for generating different materials """
        if self._material is None:
            self._material = """\
new THREE.MeshPhysicalMaterial( {
        metalness: 0.0,
        roughness: 0.5,
        reflectivity: 0.7,
        clearCoat: 0.0,
        vertexColors: THREE.FaceColors
    }  
"""
        return self._material

    @material.setter
    def material(self, material: str):
        self._material = material


    def as_javascript(self) -> str:
        from rasputin.html_writer import vertex_field_to_vertex_values, point_vector_to_lines, face_and_point_vector_to_lines
        res = "{"
        # write vertices
        vertices = "\n".join(face_and_point_vector_to_lines(name=None,
                                                            face_vector=self.faces,
                                                            point_vector=self.points))
        res += f"vertices: {vertices},\n"

        #write normals
        normals = vertex_field_to_vertex_values(vertex_field=np.asarray(self.point_normals),
                                                faces=self.faces, points=self.points)
        normals = "\n".join(point_vector_to_lines(name=None, point_vector=normals))
        res += f"normals: {normals},\n"

        # write colors
        colors = "\n".join(point_vector_to_lines(name=None, point_vector=self.colors))
        res += f"colors: {colors},\n"

        # write material

        res += f"material: {self.material}"

        res += "}"

