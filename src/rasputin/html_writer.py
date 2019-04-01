from pathlib import Path
import shutil
import numpy as np
import seaborn as sns
from typing import List, Union, Optional, Tuple
from pkg_resources import resource_filename

from rasputin import triangulate_dem
from rasputin.avalanche import varsom_angles


def face_vector_to_lines(*,
                         name: str,
                         face_vector: triangulate_dem.FaceVector)-> List[str]:
    lines = [f'const {name} = [\n']
    for v in face_vector:
        lines.append(f"    {', '.join([str(s) for s in v])},\n")
    lines.append('];\n\n')
    return lines


def point_vector_to_lines(*,
                          name: Optional[str],
                          point_vector: Union[np.ndarray, triangulate_dem.PointVector])-> List[str]:
    if name is not None:
        lines = [f'const {name} = new Float32Array( [\n']
    else:
        lines = ["new Float32Array( [\n"]
    for v in point_vector:
        lines.append(f"    {', '.join([str(s) for s in v])},\n")
    if name is not None:
        lines.append(']);\n\n')
    else:
        lines.append("])\n")
    return lines


def face_and_point_vector_to_lines(*,
                                   name: Optional[str],
                                   face_vector: triangulate_dem.FaceVector,
                                   point_vector: triangulate_dem.PointVector) -> List[str]:
    if name is not None:
        lines = [f'const {name} = new Float32Array( [\n']
    else:
        lines = [f"new Float32Array([\n"]
    for f in face_vector:
        for i in f:
            lines.append(f"    {', '.join([str(s) for s in point_vector[i]])},\n")
    if name is not None:
        lines.append(']);\n\n')
    else:
        lines.append("])")
    return lines


def face_field_to_vertex_values(*,
                                face_field: np.ndarray,
                                faces: triangulate_dem.FaceVector,
                                points: triangulate_dem.PointVector) -> np.ndarray:
    """

    :param field: scalar field over faces
    :param faces: face index vector
    :param points: points
    :return: field values for points
    """
    assert len(face_field) == len(faces)
    vertex_field = np.empty((3*len(faces), face_field.shape[1]))
    for i, face in enumerate(faces):
        vertex_field[i*3:i*3 + 3] = face_field[i]
    return vertex_field


def vertex_field_to_vertex_values(*,
                                  vertex_field: np.ndarray,
                                  faces: triangulate_dem.FaceVector,
                                  points: triangulate_dem.PointVector) -> np.ndarray:
    """

    :param field: color (vector) field over points
    :param faces: face index vector
    :param points: points
    :return: field values for points
    """
    assert len(vertex_field) == len(points)
    res = np.empty((3*len(faces), vertex_field.shape[1]))
    for i, face in enumerate(faces):
        res[i*3:i*3 + 3] = [vertex_field[f] for f in face]
    return res


def color_field_by_height(*, points: triangulate_dem.PointVector) -> np.ndarray:
    vertex_heights = np.asarray(points)[:, 2].copy()
    factor = 2
    vertex_heights -= min(vertex_heights)
    vertex_heights /= factor*max(vertex_heights)
    vertex_heights += 1 - 1/factor
    colors = np.empty((len(vertex_heights), 3))
    colors[:, 0] = vertex_heights
    colors[:, 1] = vertex_heights
    colors[:, 2] = vertex_heights
    return colors

def add_slope_colors(*,
                     normals: triangulate_dem.PointVector,
                     colors: np.ndarray) -> np.ndarray:
    slopes = np.asanyarray(triangulate_dem.compute_slopes(normals))
    colors[slopes < 5.0e-2] = [0, 0, 1]
    #colors[np.logical_and(slopes < 55/180*np.pi, slopes > 30/180*np.pi)] = [0.5, 0, 0]
    colors[slopes >= 55/180*np.pi] = [0.2, 0.2, 0.2]
    return colors


def color_field_by_slope(*,
                         normals: triangulate_dem.PointVector) -> np.ndarray:
    slopes = np.asanyarray(triangulate_dem.compute_slopes(normals))
    colors = np.empty((len(slopes), 3))
    colors[slopes < 5.0e-2] = [0, 0, 1]
    colors[slopes > 5.0e-2] = [1, 1, 1]
    colors[slopes > 30/180*np.pi] = [1, 0, 0]
    return colors


def mesh_with_avalanche_danger(*,
                               normals: triangulate_dem.PointVector,
                               points: triangulate_dem.PointVector,
                               faces: triangulate_dem.FaceVector,
                               avalanche_problems: list) -> List[Tuple[triangulate_dem.PointVector, triangulate_dem.FaceVector]]:
    aspects = np.asarray(triangulate_dem.compute_aspect(normals))
    slopes = np.asarray(triangulate_dem.compute_slopes(normals))
    colors = np.ones((len(aspects), 3))
    angles = np.asarray(varsom_angles)/180*np.pi
    mesh_list = []
    for problem in avalanche_problems:
        counter = 0
        point_map = {}
        av_points = triangulate_dem.PointVector()
        av_faces = triangulate_dem.FaceVector()
        expositions = [int(s) for s in problem["expositions"]]
        exposed_angles = angles[expositions]
        level = problem["level"]
        heights = problem["heights"]
        for i, face in enumerate(faces):
            avg_h = np.mean([points[f][2] for f in face])
            risk = False
            for hmin, hmax in heights:
                if hmin < avg_h < hmax:
                    aspect = aspects[i]
                    for a_min, a_max in exposed_angles:
                        if a_min < aspect < a_max:
                            if slopes[i] > 30/180*np.pi:
                                risk = True
                                break
                if risk:
                    break
            if risk:
                for idx in face:
                    if idx not in point_map:
                        point_map[idx] = counter
                        counter += 1
                av_faces.append([point_map[idx] for idx in face])
        rpm = {v: k for (k, v) in point_map.items()}
        for k in sorted(rpm):
            pt = points[rpm[k]]
            pt[2] += 1
            av_points.append(pt)
        mesh_list.append((av_points, av_faces))
    return mesh_list

def color_field_by_avalanche_danger(*,
                                    normals: triangulate_dem.PointVector,
                                    points: triangulate_dem.PointVector,
                                    faces: triangulate_dem.FaceVector,
                                    avalanche_problems: list,
                                    colors: Optional[np.ndarray]) -> np.ndarray:
    aspects = np.asarray(triangulate_dem.compute_aspect(normals))
    slopes = np.asarray(triangulate_dem.compute_slopes(normals))
    if colors is None:
        colors = np.ones((len(aspects), 3))
    else:
        colors = colors.copy()
    angles = np.asarray(varsom_angles)/180*np.pi
    for problem in avalanche_problems:
        expositions = [bool(int(s)) for s in problem["expositions"]]
        exposed_angles = angles[expositions]
        level = problem["level"]
        heights = problem["heights"]
        print(level, heights)
        for i, face in enumerate(faces):
            avg_h = np.mean([points[f][2] for f in face])
            risk = False
            for hmin, hmax in heights:
                if hmin < avg_h < hmax:
                    aspect = aspects[i]
                    for a_min, a_max in exposed_angles:
                        if a_min < aspect < a_max:
                            if slopes[i] > 30/180*np.pi:
                                risk = True
                                break
                        elif a_min > a_max:
                            if a_min < aspect or aspect < a_max:
                                if slopes[i] > 30/180*np.pi:
                                    risk = True
                                    break
                if risk:
                    break
            if risk:
                if level == 1:
                    colors[i] = [0, 1, 0]
                elif level == 2:
                    colors[i] = [1, 1, 0]
                elif level == 3:
                    colors[i] = [1, 0, 0]
                elif level == 4:
                    colors[i] = [0, 0, 0]
    return colors


def color_field_by_aspect(*, normals: triangulate_dem.PointVector)->np.ndarray:

    palette = sns.color_palette("pastel", n_colors=len(varsom_angles))
    aspects = np.asanyarray(triangulate_dem.compute_aspect(normals))
    colors = np.empty((len(aspects), 3))
    for angle, col in zip(varsom_angles, palette):
        range = np.logical_and(aspects >= angle[0]/180*np.pi, aspects < angle[1]/180*np.pi)
        colors[range] = col
    return colors


class Geometry:

    def __init__(self, *,
                 points: triangulate_dem.PointVector,
                 faces: triangulate_dem.FaceVector,
                 base_color: Tuple[float, float, float]):
        self.points = points
        self.faces = faces
        self._base_color = base_color
        self._surface_normals = None
        self._point_normals = None
        self._aspects = None
        self._slopes = None
        self._colors = None

    @property
    def point_normals(self) -> triangulate_dem.PointVector:
        if self._point_normals is None:
            self._point_normals = triangulate_dem.point_normals(self.points, self.faces)
        return self._point_normals

    @property
    def surface_normals(self) -> triangulate_dem.PointVector:
        if self._surface_normals is None:
            self._surface_normals = triangulate_dem.surface_normals(self.points, self.faces)
        return self._surface_normals

    @property
    def aspects(self) -> np.ndarray:
        if self._aspects is None:
            self._aspects = np.asarray(triangulate_dem.compute_aspect(self.surface_normals))
        return self._aspects

    @property
    def slopes(self) -> np.ndarray:
        if self._slopes is None:
            self._slopes = np.asarray(triangulate_dem.compute_slopes(self.surface_normals))
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
        return """\
new THREE.MeshPhysicalMaterial( {
        metalness: 0.0,
        roughness: 0.5,
        reflectivity: 0.7,
        clearCoat: 0.0,
        vertexColors: THREE.FaceColors
    }  
"""

    def as_javascript(self) -> str:
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


def write_geometries(*,
                  geometries: List[Geometry], output_dir: Path) -> None:
    if output_dir.exists():
        shutil.rmtree(str(output_dir))

    templates_path = Path(resource_filename(__name__, "web"))

    # copy all template files
    shutil.copytree(str(templates_path), str(output_dir))

    data_js = output_dir / "data.js"

    with data_js.open("w") as tf:
        tf.write("const geometries = [\n")
        for geometry in geometries:
            tf.write(geometry.as_javascript() + ",\n")

        tf.write("];\n")
        tf.write("const data = {geometries};\n")


def write_mesh(*,
               pts: triangulate_dem.PointVector,
               faces: triangulate_dem.FaceVector,
               normals: triangulate_dem.PointVector,
               features: list,
               vertex_field: Optional[np.ndarray]=None,
               face_field: Optional[np.ndarray]=None,
               output_dir: Path) -> None:
    """
    This function writes a output html web-site code for the give mesh to output dir.
    """
    if output_dir.exists():
        shutil.rmtree(str(output_dir))

    templates_path = Path(resource_filename(__name__, "web"))

    # copy all template files
    shutil.copytree(str(templates_path), str(output_dir))

    data_js = output_dir / "data.js"

    lines = []
    lines.extend(face_and_point_vector_to_lines(name='vertices', face_vector=faces, point_vector=pts))
    lines.append("const normals = new Float32Array( [\n")
    for x,y,z in vertex_field_to_vertex_values(vertex_field=np.asarray(normals), faces=faces, points=pts):
        lines.append(f"{x}, {y}, {z},\n")
    lines.append("\n] );\n\n")

    lines.append("const face_field = new Float32Array( [\n")
    if face_field is not None:
        for r,g,b in face_field_to_vertex_values(face_field=face_field, faces=faces, points=pts):
            lines.append(f"{r}, {g}, {b},\n")
    lines.append("\n] );\n\n")

    lines.append("const vertex_field = new Float32Array( [\n")
    if vertex_field is not None:
        for r,g,b in vertex_field_to_vertex_values(vertex_field=vertex_field, faces=faces, points=pts):
            lines.append(f"{r}, {g}, {b},\n")
    lines.append("\n] );\n\n")
    lines.append("const features = [ \n")
    for faces in features:
        lines.extend(face_and_point_vector_to_lines(name=None, face_vector=faces, point_vector=pts))
        lines.append(",\n")
    if features:
        del lines[-1]

    lines.append("];")

    lines.append("\nconst data = {vertices, normals, face_field, vertex_field, features};\n")

    with data_js.open(mode="w") as f:
        f.writelines(lines)
