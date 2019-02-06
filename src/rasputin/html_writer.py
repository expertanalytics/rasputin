from pathlib import Path
import shutil
from typing import List, Union, Optional

from . import triangulate_dem


def convert_variable_to_lines(*,
                              name: str,
                              variable: Union[triangulate_dem.PointVector, triangulate_dem.FaceVector],
                              type: Optional[str]=None)-> List[str]:
    if type:
        lines = [f'const {name} = new {type}( [\n']
    else:
        lines = [f'const {name} = [\n']
    for v in variable:
        lines.append(f"    {', '.join([str(s) for s in v])},\n")

    lines.append(']')
    if type:
        lines.append(')')
    lines.append(';\n\n')
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
    lines.extend(convert_variable_to_lines(name='indices', variable=faces))
    lines.extend(convert_variable_to_lines(name='vertices', type='Float32Array', variable=pts))

    lines.append("const colors = new Float32Array( [\n")
    for p in pts:
        lines.append("0.0,  1.0,  0.0,\n")
    lines.append("\n] );\n\n")

    lines.append("\nconst data = {indices, vertices, colors};\n")

    with data_js.open(mode="w") as f:
        f.writelines(lines)
