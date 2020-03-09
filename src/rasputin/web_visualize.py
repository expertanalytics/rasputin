import os
import sys
from pathlib import Path
import numpy as np
import pyproj
import argparse
import logging
import yaml

from rasputin import rasputin_data_dir
from rasputin.tin_repository import TinRepository
from rasputin import triangulate_dem
from rasputin.reader import RasterRepository
from rasputin.geometry import Geometry, write_scene
from rasputin.material import avalanche_material, lake_material, terrain_material
from rasputin import avalanche
from rasputin.avalanche import varsom_angles
from rasputin.py2js import construct_material


def depr_web_visualize():
    return

    """
    Avalance Forecast Visualization Example.

    Items to develop:
     * Construct several color maps for each avalanche danger, and make the web app offer selections
     * Alternative approach is to partition the whole area into disjoint topologies and switch between these
     * Use different textures for different terrain types (under development)
     * Use more of the information from varsom.no (and perhaps alpha blending) to better display avalanche dangers

    """
    data_dir = rasputin_data_dir / "dem_archive"
    files = data_dir.glob("*.tif")
    if not files:
        raise RuntimeError(f"No GeoTIFF files found in {data_dir.absolute()}, giving up.")

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-lat", type=float, default=60.898468, help="Latitude of center coordinate")
    arg_parser.add_argument("-lon", type=float, default=8.530918, help="Longitude of center coordinate")
    arg_parser.add_argument("-dx", type=float, default=5000, help="Distance in meters")
    arg_parser.add_argument("-dy", type=float, default=5000, help="Distance in meters")
    arg_parser.add_argument("-ratio", type=float, default=0.4, help="Mesh coarsening factor in [0, 1]")
    arg_parser.add_argument("-a", default=False, action="store_true", help="Add avalanche forecast")
    arg_parser.add_argument("-output", type=str, default="web", help="Output directory")
    res = arg_parser.parse_args(sys.argv[1:])

    # Define area of interest
    x0 = res.lon
    y0 = res.lat
    dx = res.dx
    dy = res.dy

    input_crs = pyproj.CRS.from_string("+init=EPSG:4326")
    target_crs = pyproj.CRS.from_epsg(32633)

    avalanche_problems = []
    if res.a:
        avalanche_details, avalanche_meta = avalanche.get_forecasts(x=x0, y=y0, crs=input_crs)
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
                        danger_interval = [[min_h, max_h]]
                    elif hf == 2:
                        min_h = 0
                        max_h = p["ExposedHeight1"]
                        danger_interval = [[min_h, max_h]]
                    elif hf == 3:
                        min_h1 = 0
                        max_h1 = p["ExposedHeight1"]
                        min_h2 = p["ExposedHeight2"]
                        max_h2 = 20000
                        danger_interval = [[min_h1, max_h1], [min_h2, max_h2]]
                    elif hf == 4:
                        min_h = p["ExposedHeight1"]
                        max_h = p["ExposedHeight2"]
                        danger_interval = [[min_h, max_h]]
                    avalanche_problems.append({"expositions": expositions,
                                               "level": level,
                                               "heights": danger_interval})
    raise DeprecationWarning("This code is deprecated!")
    raster_coords = RasterRepository(directory=data_dir).read(x=x0,
                                                              y=y0,
                                                              dx=dx,
                                                              dy=dy,
                                                              input_coordinate_system=input_coordinate_system,
                                                              target_coordinate_system=target_coordinate_system)

    points, faces = lindstrom_turk_by_ratio(raster_coords, res.ratio)
    lakes, terrain = triangulate_dem.extract_lakes(points, faces)
    geometries = []
    lake_geometry = Geometry(points=points, faces=lakes, base_color=(0, 0, 1), material=lake_material)
    terrain_geometry = Geometry(points=points, faces=terrain, base_color=(1, 1, 1), material=terrain_material)
    geometries.append(lake_geometry)

    if res.a:
        problem = avalanche_problems[0]
        expositions = [bool(int(s)) for s in problem["expositions"]]
        angles = np.asarray(varsom_angles)/180*np.pi
        exposed_angles = triangulate_dem.point2_vector(angles[expositions].tolist())
        heights = triangulate_dem.point2_vector(problem["heights"])

        avalanche_risk, safe_terrain = triangulate_dem.extract_avalanche_expositions(terrain_geometry.points,
                                                                                     terrain_geometry.faces,
                                                                                     exposed_angles,
                                                                                     heights)
        geometries.append(Geometry(points=points,
                                   faces=avalanche_risk,
                                   base_color=(1, 0, 0),
                                   material=avalanche_material))
        geometries.append(Geometry(points=points,
                                   faces=safe_terrain,
                                   base_color=(1, 1, 1),
                                   material=terrain_material))

    else:
        geometries.append(terrain_geometry)

    output = Path(res.output).absolute()
    write_scene(geometries=geometries, output=output)
    print(f"""Successfully generated a web_gl based TIN visualizer in {output}.
To see it, please run:
cd {output}
{Path(sys.executable).name} -m http.server 8080
Then visit http://localhost:8080
""")


def visualize_tin():
    logging.basicConfig(level=logging.CRITICAL, format='Rasputin[%(levelname)s]: %(message)s')
    logger = logging.getLogger()

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("uid", type=str, help="Tin repository uid for mesh to visualize.")
    arg_parser.add_argument("-output", type=str, default="web_viz", help="Directory for web content.")
    arg_parser.add_argument("-material", type=str, default="", help="Optional texture/material spec in yaml")
    arg_parser.add_argument("-silent", action="store_true", help="Run in silent mode")
    res = arg_parser.parse_args(sys.argv[1:])
    if not res.silent:
        logger.setLevel(logging.INFO)

    tin_archive = rasputin_data_dir / "tin_archive"
    tin_repo = TinRepository(path=tin_archive)
    info = tin_repo.info(uid=res.uid)

    material_by_id = {}
    if res.material:
        materials = yaml.load(open(res.material))
        for material in materials:
            material_ids = material.pop("land_type_ids")
            for material_id in material_ids:
                material_by_id[material_id] = material

    geometries = []
    for land_cover in info["info"]["land_covers"]:
        id = land_cover[0]
        geometry = tin_repo.extract(uid=res.uid, face_id=id)
        if id in material_by_id:
            geometry.material_spec = material_by_id[id]
        geometry.material_spec["land_cover"] = land_cover[1].replace('_', ' ')
        geometries.append(geometry)
    output = Path(res.output).absolute()
    write_scene(geometries=geometries, output=output)
    print(f"""Successfully generated a web_gl based TIN visualizer in {output}.
To see it, please run:
cd {output}
{Path(sys.executable).name} -m http.server 8080
Then visit http://localhost:8080
""")

