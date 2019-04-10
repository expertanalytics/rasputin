from pathlib import Path
import shutil
import numpy as np
import seaborn as sns
from typing import List, Union, Optional, Tuple, Generator
from pkg_resources import resource_filename

from rasputin import triangulate_dem
from rasputin.avalanche import varsom_angles
from rasputin.geometry import Geometry


def face_vector_to_lines(*,
                         name: str,
                         face_vector: triangulate_dem.FaceVector)-> Generator[str, None, None]:
    yield f'const {name} = ['
    for v in face_vector:
        yield f"    {', '.join([str(s) for s in v])},"
    yield "];\n"


def point_vector_to_lines(*,
                          name: Optional[str],
                          point_vector: Union[np.ndarray, triangulate_dem.PointVector])-> Generator[str, None, None]:
    if name is not None:
        yield f'const {name} = new Float32Array( ['
    else:
        yield "new Float32Array( ["
    for v in point_vector:
        yield f"    {', '.join([str(s) for s in v])},"
    if name is not None:
        yield '] );'
    else:
        yield "] )"


def face_and_point_vector_to_lines(*,
                                   name: Optional[str],
                                   face_vector: triangulate_dem.FaceVector,
                                   point_vector: triangulate_dem.PointVector) -> Generator[str, None, None]:
    if name is not None:
        yield f'const {name} = new Float32Array( ['
    else:
        yield f"new Float32Array( ["
    for f in face_vector:
        for i in f:
            yield f"    {', '.join([str(s) for s in point_vector[i]])},"
    if name is not None:
        yield "] );"
    else:
        yield "] )"


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
    #colors[slopes < 5.0e-2] = [0, 0, 1]
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
    aspects = np.asarray(triangulate_dem.compute_aspects(normals))
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
    aspects = np.asarray(triangulate_dem.compute_aspects(normals))
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
    aspects = np.asanyarray(triangulate_dem.compute_aspects(normals))
    colors = np.empty((len(aspects), 3))
    for angle, col in zip(varsom_angles, palette):
        range = np.logical_and(aspects >= angle[0]/180*np.pi, aspects < angle[1]/180*np.pi)
        colors[range] = col
    return colors



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
        shutil.rmtree(output_dir)

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
