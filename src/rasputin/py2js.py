from typing import Generator, Optional
import numpy as np


def face_vector_to_lines(*,
                         name: str,
                         face_vector: np.ndarray)-> Generator[str, None, None]:
    yield f'const {name} = ['
    for v in face_vector:
        yield f"    {', '.join([str(s) for s in v])},"
    yield "];\n"


def point_vector_to_lines(*,
                          name: Optional[str],
                          point_vector: np.ndarray)-> Generator[str, None, None]:
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
                                   face_vector: np.ndarray,
                                   point_vector: np.ndarray) -> Generator[str, None, None]:
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



def construct_material(spec: dict) -> str:
    params = spec['material_params']
    params_str = ",\n".join([ f"        {k}: {v}" for k,v in params.items()])
    constructor = f"""function (geometry){{
"""
    if "texture_file_name" in spec:
        constructor += f"""\
    const texture_loader = new THREE.TextureLoader();
    var texture = texture_loader.load("{spec["texture_file_name"]}");
    texture.name = "{spec["land_cover"]}";
    texture.wrapS = THREE.RepeatWrapping;
    texture.wrapT = THREE.RepeatWrapping;
"""
    constructor += f"""    var material = new {spec['material_type']}({{
{params_str} 
    }});
"""
    if "texture_file_name" in spec:
        constructor += "    return [new THREE.Mesh(geometry, material), texture];"
    else:
        constructor += "    return [new THREE.Mesh(geometry, material), null];"
    constructor += """
}
"""
    return constructor

