from pathlib import Path
import shutil
import numpy as np
import seaborn as sns
from typing import List, Union, Optional
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
                          name: str,
                          point_vector: triangulate_dem.PointVector)-> List[str]:
    lines = [f'const {name} = new Float32Array( [\n']
    for v in point_vector:
        lines.append(f"    {', '.join([str(s) for s in v])},\n")
    lines.append(']);\n\n')
    return lines


def face_and_point_vector_to_lines(*,
                                   name: str,
                                   face_vector: triangulate_dem.FaceVector,
                                   point_vector: triangulate_dem.PointVector) -> List[str]:
    lines = [f'const {name} = new Float32Array( [\n']
    for f in face_vector:
        for i in f:
            lines.append(f"    {', '.join([str(s) for s in point_vector[i]])},\n")
    lines.append(']);\n\n')
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


def color_field_by_avalanche_danger(*,
                                    normals: triangulate_dem.PointVector,
                                    points: triangulate_dem.PointVector,
                                    faces: triangulate_dem.FaceVector,
                                    avalanche_problems: list) -> np.ndarray:
    aspects = np.asarray(triangulate_dem.compute_aspect(normals))
    slopes = np.asarray(triangulate_dem.compute_slopes(normals))
    colors = np.ones((len(aspects), 3))
    angles = np.asarray(varsom_angles)/180*np.pi
    for problem in avalanche_problems:
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


def write_mesh(*,
               pts: triangulate_dem.PointVector,
               faces: triangulate_dem.FaceVector,
               normals: triangulate_dem.PointVector,
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
    lines.append("\nconst data = {vertices, normals, face_field, vertex_field};\n")

    with data_js.open(mode="w") as f:
        f.writelines(lines)
