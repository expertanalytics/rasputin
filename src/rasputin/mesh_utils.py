import numpy as np
from rasputin import triangulate_dem as td

def face_field_to_vertex_values(*,
                                face_field: np.ndarray,
                                faces: td.FaceVector,
                                points: td.PointVector) -> np.ndarray:
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
                                  faces: td.FaceVector,
                                  points: td.PointVector) -> np.ndarray:
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

