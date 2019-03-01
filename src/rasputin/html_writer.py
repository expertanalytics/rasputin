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


def write_mesh(*,
               pts: triangulate_dem.PointVector,
               faces: triangulate_dem.FaceVector,
               normals: triangulate_dem.PointVector,
               output_dir: Path) -> None:
    """
    This function writes a output html web-site code for the give mesh to output dir.
    """
    if output_dir.exists():
        shutil.rmtree(str(output_dir))

    templates_path = Path(resource_filename(__name__, "web"))
    #templates_path = Path(__file__).parent.parent.parent / "web"

    # copy all template files
    shutil.copytree(str(templates_path), str(output_dir))

    data_js = output_dir / "data.js"

    lines = []
    lines.extend(face_and_point_vector_to_lines(name='vertices', face_vector=faces, point_vector=pts))
    lines.extend(point_vector_to_lines(name='normals', point_vector=normals))


    face_scalar = np.zeros(3*len(faces))
    for i, face in enumerate(faces):
        face_scalar[i*3:i*3 + 3] = [pts[f][2] for f in face]
        if i == 100:
            print([pts[f][2] for f in face])

    factor = 4
    face_scalar -= min(face_scalar)
    face_scalar /= factor*max(face_scalar)
    face_scalar += 1 - 1/factor
    print(face_scalar[100*3:101*3])
    lines.append("const face_colors = new Float32Array( [\n")
    for h in face_scalar:
        r, g, b = np.random.random(3)
        #lines.append(f"{r}, {g}, {b},\n")
        lines.append(f"{h}, {h}, {h},\n")
    lines.append("\n] );\n\n")

    lines.append("const vertex_colors = new Float32Array( [\n")

    heights = np.array([p[2] for p in pts])
    heights -= min(heights)
    heights /= 2*max(heights)
    heights += 0.5
    for h in heights:
        lines.append(f"{h}, {h}, {h},\n")
    lines.append("\n] );\n\n")
    lines.append("\nconst data = {vertices, normals, face_colors, vertex_colors};\n")

    with data_js.open(mode="w") as f:
        f.writelines(lines)
