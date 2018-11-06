import sys
from pathlib import Path
import argparse
from logging import getLogger
from rasputin.reader import read_raster_file
from rasputin.writer import write_ascii_off
from rasputin.writer import write_vtk
from rasputin.calculate import compute_shade
from rasputin.triangulate_dem import lindstrom_turk_by_ratio
from rasputin.triangulate_dem import lindstrom_turk_by_size


def geo_tiff_reader():
    logger = getLogger()
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("input", type=str, metavar="FILENAME")
    arg_parser.add_argument("-output", type=str, default="output.off", help="Surface mesh file name")
    arg_parser.add_argument("-start_x", type=int, default=0, help="Start index in x-direction")
    arg_parser.add_argument("-stop_x", type=int, default=-1, help="Stop index in x-direction")
    arg_parser.add_argument("-start_y", type=int, default=0, help="Start index in y-direction")
    arg_parser.add_argument("-stop_y", type=int, default=-1, help="Stop index in y-direction")
    arg_parser.add_argument("-sun_x", type=float, help="Sun ray x component")
    arg_parser.add_argument("-sun_y", type=float, help="Sun ray y component")
    arg_parser.add_argument("-sun_z", type=float, help="Sun ray z component")
    arg_parser.add_argument("-loglevel", type=int, default=1, help="Verbosity")
    group = arg_parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-ratio", type=float, help="Edge ratio between original meshed raster and result")
    group.add_argument("-size", type=float, help="Number of edges in the generated surface mesh")
    res = arg_parser.parse_args(sys.argv[1:])
    logger.setLevel(res.loglevel)
    raster_coords = read_raster_file(filepath=Path(res.input),
                                     start_x=res.start_x,
                                     stop_x=res.stop_x,
                                     start_y=res.start_y,
                                     stop_y=res.stop_y)
    logger.debug(f"Original: {len(raster_coords)}")
    if res.ratio:
        pts, faces = lindstrom_turk_by_ratio(raster_coords, res.ratio)
    else:
        pts, faces = lindstrom_turk_by_size(raster_coords, res.size)
    logger.debug(f"Result: {len(raster_coords)}")
    if res.sun_x is not None or res.sun_y is not None or res.sun_z is not None:
        assert res.sun_x is not None
        assert res.sun_y is not None
        assert res.sun_z is not None
        shade_field = compute_shade(pts, faces, [res.sun_x, res.sun_y, res.sun_z])
    else:
        shade_field = None
    output_path = Path(res.output)
    if output_path.suffix == ".vtk":
        write_vtk(pts, faces, shade_field, output_path)
    elif output_path.suffix == ".off":
        write_ascii_off(pts, faces, output_path)
    else:
        raise RuntimeError(f"Unknown output file format {output_path.suffix}")
