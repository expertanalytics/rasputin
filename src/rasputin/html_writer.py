from pathlib import Path
import shutil
import numpy as np
from typing import List, Union, Optional
from pkg_resources import resource_filename

from rasputin import triangulate_dem


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
    vertex_field = np.empty(3*len(faces))
    for i, face in enumerate(faces):
        vertex_field[i*3:i*3 + 3] = face_field[i]
    return vertex_field


def vertex_field_to_vertex_values(*,
                                  vertex_field: np.ndarray,
                                  faces: triangulate_dem.FaceVector,
                                  points: triangulate_dem.PointVector) -> np.ndarray:
    """

    :param field: scalar field over points
    :param faces: face index vector
    :param points: points
    :return: field values for points
    """
    assert len(vertex_field) == len(points)
    res = np.empty(3*len(faces))
    for i, face in enumerate(faces):
        res[i*3:i*3 + 3] = [vertex_field[f] for f in face]
    return res


def color_field_by_height(*, points: triangulate_dem.PointVector) -> np.ndarray:
    vertex_heights = np.asarray(points)[:, 2].copy()
    factor = 2
    vertex_heights -= min(vertex_heights)
    vertex_heights /= factor*max(vertex_heights)
    vertex_heights += 1 - 1/factor
    return vertex_heights


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
    lines.extend(point_vector_to_lines(name='normals', point_vector=normals))

    lines.append("const face_field = new Float32Array( [\n")
    if face_field is not None:
        for v in face_field_to_vertex_values(face_field=face_field, faces=faces, points=pts):
            lines.append(f"{v}, {v}, {v},\n")
    lines.append("\n] );\n\n")

    lines.append("const vertex_field = new Float32Array( [\n")
    if vertex_field is not None:
        for v in vertex_field_to_vertex_values(vertex_field=vertex_field, faces=faces, points=pts):
            lines.append(f"{v}, {v}, {v},\n")
    lines.append("\n] );\n\n")
    lines.append("\nconst data = {vertices, normals, face_field, vertex_field};\n")

    with data_js.open(mode="w") as f:
        f.writelines(lines)
