import sys
from pathlib import Path
import argparse
from logging import getLogger

from shapely.geometry import Polygon

from rasputin.writer import write
from rasputin.reader import read_raster_file
from rasputin.calculate import compute_shade
from rasputin.triangulate_dem import lindstrom_turk_by_ratio
from rasputin.triangulate_dem import lindstrom_turk_by_size
from rasputin.triangulate_dem import surface_normals, orient_tin, compute_slopes


def geo_tiff_reader():
    logger = getLogger()
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("input", type=str, metavar="FILENAME")
    arg_parser.add_argument("-output", type=str, default="output.off", help="Surface mesh file name")

    coords = arg_parser.add_argument_group("Window coordinates",
                                           "Coordinates of axially aligned view-window to the raster data")
    coords.add_argument("-x0", type=float, default=None, help="Start coordinate in x-direction")
    coords.add_argument("-x1", type=float, default=None, help="Stop coordinate in x-direction")
    coords.add_argument("-y0", type=float, default=None, help="Start coordinate in y-direction")
    coords.add_argument("-y1", type=float, default=None, help="Stop coordinate in y-direction")

    arg_parser.add_argument("-sun_x", type=float, help="Sun ray x component")
    arg_parser.add_argument("-sun_y", type=float, help="Sun ray y component")
    arg_parser.add_argument("-sun_z", type=float, help="Sun ray z component")
    arg_parser.add_argument("-n", action="store_true", help="Compute surface normals")
    arg_parser.add_argument("-slope", action="store_true", help="Compute slope")
    arg_parser.add_argument("-loglevel", type=int, default=1, help="Verbosity")
    group = arg_parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-ratio", type=float, help="Edge ratio between original meshed raster and result")
    group.add_argument("-size", type=int, help="Number of edges in the generated surface mesh")

    res = arg_parser.parse_args(sys.argv[1:])
    logger.setLevel(res.loglevel)

    # Determine region of interest
    coordinates = (res.x0, res.y0, res.x1, res.y1)
    if all(coord is None for coord in coordinates):
        polygon = None

    elif all(coord is not None for coord in coordinates):
        polygon = Polygon([(res.x0, res.y0),
                           (res.x1, res.y0),
                           (res.x1, res.y1),
                           (res.x0, res.y1)])

    else:
        raise ValueError("Either all or none of the window coordinates must be provided")

    # Read from raster
    rasterdata = read_raster_file(filepath=Path(res.input),
                                   polygon=polygon)
    m, n = rasterdata.array.shape
    logger.debug(f"Original: {m * n}")
    logger.critical(rasterdata.info)
    if res.ratio:
        pts, faces = lindstrom_turk_by_ratio(rasterdata._cpp, res.ratio)
    else:
        pts, faces = lindstrom_turk_by_size(rasterdata._cpp, res.size)
    logger.debug(f"Result: {len(pts)}")

    # Compute normals and re-orient before shadow computation
    normals = orient_tin(pts, faces)

    # Compute additional fields
    fields = {}

    if res.sun_x is not None or res.sun_y is not None or res.sun_z is not None:
        assert res.sun_x is not None
        assert res.sun_y is not None
        assert res.sun_z is not None
        fields["shade"] = compute_shade(pts=pts, faces=faces, sun_ray=[res.sun_x, res.sun_y, res.sun_z])

    if res.n:
        fields["surface_normal"] = normals

    if res.slope:
        slopes = compute_slopes(normals)
        fields["slope"] = slopes

    output_path = Path(res.output)
    write(pts=pts, faces=faces, filepath=output_path, fields=fields)
