import os
import sys
from pathlib import Path
import numpy as np
import pyproj
import argparse

from rasputin import triangulate_dem
from rasputin.html_writer import (write_mesh, add_slope_colors, color_field_by_avalanche_danger)
from rasputin.reader import RasterRepository
from rasputin.triangulate_dem import lindstrom_turk_by_ratio
from rasputin.geometry import Geometry
from rasputin import avalanche
from rasputin.avalanche import varsom_angles


def web_visualize():
    """
    Avalance Forecast Visualization Example.

    Items to develop:
     * Construct several color maps for each avalanche danger, and make the web app offer selections
     * Alternative approach is to partition the whole area into disjoint topologies and switch between these
     * Use different textures for different terrain types (under development)
     * Use more of the information from varsom.no (and perhaps alpha blending) to better display avalanche dangers

    """
    if "RASPUTIN_DATA_DIR" in os.environ:
        data_dir = Path(os.environ["RASPUTIN_DATA_DIR"])
    else:
        #  data_dir = Path(os.environ["HOME"]) /"projects" / "rasputin_data" / "dem_archive"
        data_dir = Path(".") / "dem_archive"
        print(f"WARNING: No raster archive directory specified, assuming {data_dir.absolute()}")
    files = data_dir.glob("*.tif")
    if not files:
        raise RuntimeError(f"No GeoTIFF files found in {data_dir.absolute()}, giving up.")

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-lat", type=float, default=60.898468, help="Latitude of center coordinate")
    arg_parser.add_argument("-lon", type=float, default=8.530918, help="Longitude of center coordinate")
    arg_parser.add_argument("-dx", type=float, default=5000, help="Distance in meters")
    arg_parser.add_argument("-dy", type=float, default=5000, help="Distance in meters")
    arg_parser.add_argument("-ratio", type=float, default=0.4, help="Mesh coarsening factor in [0, 1]")
    arg_parser.add_argument("-a", type=bool, default=False, help="Add avalanche forecast")
    arg_parser.add_argument("-output", type=str, default="web", help="Output directory")
    res = arg_parser.parse_args(sys.argv[1:])

    # Define area of interest
    x0 = res.lon
    y0 = res.lat
    dx = res.dx
    dy = res.dy

    input_coordinate_system = pyproj.Proj(init="EPSG:4326").definition_string()
    target_coordinate_system = pyproj.Proj(init="EPSG:32633").definition_string()

    avalanche_problems = []
    if res.a:
        avalanche_details, avalanche_meta = avalanche.get_forecasts(x=x0, y=y0, proj=input_coordinate_system)
        if avalanche_details:
            if "AvalancheProblems" in avalanche_details:
                level = int(avalanche_details["DangerLevel"])
                for p in avalanche_details["AvalancheProblems"]:
                    expositions = p["ValidExpositions"]
                    hf = p["ExposedHeightFill"]
                    danger_interval = []
                    if hf == 1:
                        min_h = p["ExposedHeight1"]
                        max_h = 20000
                        danger_interval = [(min_h, max_h)]
                    elif hf == 2:
                        min_h = 0
                        max_h = p["ExposedHeight1"]
                        danger_interval = [(min_h, max_h)]
                    elif hf == 3:
                        min_h1 = 0
                        max_h1 = p["ExposedHeight1"]
                        min_h2 = p["ExposedHeight2"]
                        max_h2 = 20000
                        danger_interval = [(min_h1, max_h1), (min_h2, max_h2)]
                    elif hf == 4:
                        min_h = p["ExposedHeight1"]
                        max_h = p["ExposedHeight2"]
                        danger_interval = [(min_h, max_h)]
                    avalanche_problems.append({"expositions": expositions,
                                               "level": level,
                                               "heights": danger_interval})
    raster_coords = RasterRepository(directory=data_dir).read(x=x0,
                                                              y=y0,
                                                              dx=dx,
                                                              dy=dy,
                                                              input_coordinate_system=input_coordinate_system,
                                                              target_coordinate_system=target_coordinate_system)

    points, faces = lindstrom_turk_by_ratio(raster_coords, res.ratio)
    lakes, terrain = triangulate_dem.extract_lakes(points, faces)
    geometries = []
    lake_geometry = Geometry(points=points, faces=lakes, base_color=(0, 0, 1))
    terrain_geometry = Geometry(points=points, faces=terrain, base_color=(1, 1, 1))
    geometries.append(lake_geometry)

    if res.a:
        problem = avalanche_problems[0]
        expositions = [bool(int(s)) for s in problem["expositions"]]
        angles = np.asarray(varsom_angles)/180*np.pi
        exposed_angles = angles[expositions]
        avalanche_risk, safe_terrain = triangulate_dem.extract_avalanche_expositions(terrain_geometry.points,
                                                                                     terrain_geometry.faces,
                                                                                     exposed_angles,
                                                                                     np.asarray(problem["heights"]))
        geometries.append(Geometry(points=points, faces=avalanche_risk, base_color=(1, 0, 0)))
        geometries.append(Geometry(points=points, faces=safe_terrain, base_color=(1, 1, 1)))

        #colors = color_field_by_avalanche_danger(normals=normals,
        #                                        points=points,
        #                                         faces=terrain,
        #                                         avalanche_problems=avalanche_problems,
        #                                         colors=colors)
    output = Path(res.output).absolute()
    #write_geometries(geometries=(lakes_geometry, terrain_geometry))
    #write_mesh(pts=points,
    #           faces=terrain,
    #           normals=point_normals,
    #           features=[lakes],
    #           output_dir=output,
    #           face_field=colors)
    print(f"""Successfully generated a web_gl based TIN visualizer in {output}.
To see it, please run:
cd {output}
python -m http.server 8080
Then visit http://localhost:8080
""")


if __name__ == "__main__":
    web_visualize()

