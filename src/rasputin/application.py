import os
import sys
from pathlib import Path
import pprint
import pyproj
import argparse
from rasputin.reader import RasterRepository
from rasputin.tin_repository import TinRepository
from rasputin.triangulate_dem import lindstrom_turk_by_ratio


def store_tin():
    """
    Avalance Forecast Visualization Example.

    Items to develop:
     * Construct several color maps for each avalanche danger, and make the web app offer selections
     * Alternative approach is to partition the whole area into disjoint topologies and switch between these
     * Use different textures for different terrain types (under development)
     * Use more of the information from varsom.no (and perhaps alpha blending) to better display avalanche dangers

    """
    if "RASPUTIN_DATA_DIR" in os.environ:
        dem_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "dem_archive"
        tin_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "tin_archive"
    else:
        #  data_dir = Path(os.environ["HOME"]) /"projects" / "rasputin_data" / "dem_archive"
        dem_archive = Path(".") / "dem_archive"
        tin_archive = Path(".") / "tin_archive"
        print(f"WARNING: No data directory specified, assuming dem_archive {dem_archive.absolute()}")
        print(f"WARNING: No data directory specified, assuming tin_archive {tin_archive.absolute()}")
    try:
        next(dem_archive.glob("*.tif"))
    except StopIteration as si:
        raise RuntimeError(f"No GeoTIFF files found in {dem_archive.absolute()}, giving up.")
    if not tin_archive.exists():
        tin_archive.mkdir(parents=True)
    elif tin_archive.exists() and not tin_archive.is_dir():
        raise RuntimeError(f"{tin_archive} exists and is not a directory, giving up.")

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-lat", type=float, default=60.898468, help="Latitude of center coordinate")
    arg_parser.add_argument("-lon", type=float, default=8.530918, help="Longitude of center coordinate")
    arg_parser.add_argument("-dx", type=float, default=5000, help="Distance in meters")
    arg_parser.add_argument("-dy", type=float, default=5000, help="Distance in meters")
    arg_parser.add_argument("-ratio", type=float, default=0.4, help="Mesh coarsening factor in [0, 1]")
    arg_parser.add_argument("-override", action="store_true", help="Replace existing archive entry")
    arg_parser.add_argument("uid", type=str, help="Unique ID for the result TIN")
    res = arg_parser.parse_args(sys.argv[1:])

    tr = TinRepository(path=tin_archive)
    if res.uid in tr.content:
        if not res.override:
            raise RuntimeError(f"Tin archive {tin_archive.absolute()} already contains uid {res.uid}.")
        else:
            tr.delete(res.uid)
    # Define area of interest
    x0 = res.lon
    y0 = res.lat
    dx = res.dx
    dy = res.dy

    input_coordinate_system = pyproj.Proj(init="EPSG:4326").definition_string()
    target_coordinate_system = pyproj.Proj(init="EPSG:32633").definition_string()

    raster_coords = RasterRepository(directory=dem_archive).read(x=x0,
                                                                 y=y0,
                                                                 dx=dx,
                                                                 dy=dy,
                                                                 input_coordinate_system=input_coordinate_system,
                                                                 target_coordinate_system=target_coordinate_system)
    points, faces = lindstrom_turk_by_ratio(raster_coords, res.ratio)
    tr.save(uid=res.uid, points=points, faces=faces)
    meta = tr.content[res.uid]
    print(f"Successfully added uid='{res.uid}' to the tin archive {tin_archive.absolute()}, with meta info:")
    pprint.PrettyPrinter(indent=4).pprint(meta)
