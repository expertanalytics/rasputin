from pathlib import Path
import shutil

from rasputin import triangulate_dem




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
    with data_js.open(mode="w") as f:
        f.write("const indices = [\n")
        for face in faces:
            f.write(f"{', '.join([str(f) for f in face])},\n")
        f.write("\n];\n\n")
        f.write("\nconst vertices = new Float32Array( [\n")
        for p in pts:
            f.write(f"{', '.join([str(s) for s in p])},\n")
        f.write("\n] );\n\n")
        f.write("const colors = new Float32Array( [\n")
        for face in faces:
            f.write("0.0,  1.0,  0.0,\n")
            f.write("1.0,  0.0,  0.0,\n")
            f.write("0.0,  0.0,  1.0,\n\n")
        f.write("\n] );\n\n")

        f.write("\nconst data = {indices, vertices, colors};")
