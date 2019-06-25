import numpy as np
import typing as tp

from . import triangulate_dem
from .reader import GeoPolygon, Rasterdata
import meshio

class Mesh(object):
    def __init__(self, cpp_mesh: triangulate_dem.Mesh):
        self._cpp = cpp_mesh

    @classmethod
    def from_raster(cls,
                    data: tp.Union[tp.List[Rasterdata], Rasterdata],
                    domain: tp.Optional[GeoPolygon] = None) -> "Mesh":
        # Extract cpp objects
        if isinstance(data, list):
            if data[0].array.dtype == np.float64:
                rasterdata_cpp = triangulate_dem.raster_list_double()
            else:
                rasterdata_cpp = triangulate_dem.raster_list_float()

            for raster in data:
                rasterdata_cpp.add_raster(raster.to_cpp())

        else:
            rasterdata_cpp = data.to_cpp()

        if domain:
            mesh = cls(triangulate_dem.make_mesh(rasterdata_cpp, domain.to_cpp()))
        else:
            mesh = cls(triangulate_dem.make_mesh(rasterdata_cpp))

        return mesh

    @property
    def num_points(self) -> int:
        return self._cpp.num_vertices

    @property
    def num_edges(self) -> int:
        return self._cpp.num_edges

    @property
    def num_faces(self) -> int:
        return self._cpp.num_faces

    @property
    def characteristic(self) -> int:
        return self.num_points - self.num_edges + self.num_faces

    @property
    def points(self) -> triangulate_dem.point3_vector:
        array = np.asarray(self._cpp.points)
        array.flags.writeable = False
        return array

    @property
    def faces(self) -> triangulate_dem.face_vector:
        array = np.asarray(self._cpp.faces)
        array.flags.writeable = False
        return array

    @property
    def face_normals(self) -> triangulate_dem.point3_vector:
        # Only compute normals when needed, and only once
        if not hasattr(self, "_face_normals"):
            self._face_normals = triangulate_dem.surface_normals(self._cpp.points,
                                                                 self._cpp.faces)

        array = np.asarray(self._face_normals)
        array.flags.writeable = False
        return array

    def simplify(self,
                 *,
                 ratio: tp.Optional[float] = None,
                 max_size: tp.Optional[int] = None) -> "Mesh":

        """
        Simplify mesh by edge collapse, using the Lindstrom-Turk cost functional
        for minimising the deformations.

        :ratio:    Target ratio for edges in result mesh to edges in initial mesh
        :max_size: Maximum number of edges in the result mesh
        :returns:  Mesh

        """
        result = self
        if ratio is not None:
            result = self.__class__(self._cpp.lindstrom_turk_by_ratio(ratio))

        if max_size is not None:
            result = self.__class__(self._cpp.lindstrom_turk_by_size(max_size))

        # If no valid criteria return a copy (consistent with e.g. ratio > 1)
        if ratio is None and max_size is None:
            result = self.copy()

        return result

    def copy(self) -> "Mesh":
        return self.__class__(self._cpp.copy())

    def write(self, filename: str):
        """
        Write mesh to file using meshio.
        """
        meshio.write_points_cells(filename,
                                  self.points,
                                  dict(triangle=self.faces))
