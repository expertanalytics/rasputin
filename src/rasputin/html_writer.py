from pathlib import Path
import shutil
from typing import List, Union, Optional

from . import triangulate_dem


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

    templates_path = Path(__file__).parent.parent.parent / "web"

    # copy all template files
    shutil.copytree(str(templates_path), str(output_dir))

    data_js = output_dir / "data.js"

    lines = []
    lines.extend(face_vector_to_lines(name='indices', face_vector=faces))
    lines.extend(point_vector_to_lines(name='vertices', point_vector=pts))
    lines.extend(point_vector_to_lines(name='normals', point_vector=normals))

    lines.append("const colors = new Float32Array( [\n")
    for p in pts:
        lines.append("0.0,  1.0,  0.0,\n")
    lines.append("\n] );\n\n")

    lines.append("\nconst data = {indices, vertices, normals, colors};\n")

    with data_js.open(mode="w") as f:
        f.writelines(lines)
