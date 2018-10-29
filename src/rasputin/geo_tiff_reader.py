from typing import Callable, Optional, Union
import sys
import argparse
import logging
import numpy as np
from PIL import Image
from . import triangulate_dem


logger = logging.getLogger()


def read_raster_file(*,
                     filename: str,
                     start_x: int=0,
                     stop_x: int=-1,
                     start_y: int=0,
                     stop_y: int=-1) -> triangulate_dem.PointVector:
    logger.debug(f"Reading raster file {filename}")
    image = Image.open(filename)
    dx, dy, _ = image.tag.named().get("ModelPixelScaleTag", (1,1,0))
    model_tie_point = image.tag.named().get("ModelTiepointTag", None)
    x0 = y0 = 0.0
    if model_tie_point:
        x0 = model_tie_point[3]
        y0 = model_tie_point[4]
    d = np.array(image.getdata()).reshape(image.size[0], image.size[1])[start_x:stop_x, start_y:stop_y]
    raster_coordinates = triangulate_dem.PointVector()
    for (i, j), h in np.ndenumerate(d):
        raster_coordinates.append([x0 + (start_x + i)*dx, y0 + (start_y + j)*dy, h])
    logger.debug("Done")
    return raster_coordinates


def save_to_ascii_off(pts: triangulate_dem.PointVector,
                      faces: triangulate_dem.FaceVector,
                      name: str) -> None:
    with open(name, "w") as tf:
        tf.write("OFF\n")
        tf.write(f"{len(pts)} {len(faces)} 0\n")
        for pt in pts:
            tf.write(f"{pt[0]} {pt[1]} {pt[2]}\n")
        for fc in faces:
            tf.write(f"3 {fc[0]} {fc[1]} {fc[2]}\n")

def geo_tiff_reader():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("input", type=str, metavar="FILENAME")
    arg_parser.add_argument("-output", type=str, default="output.off", help="Surface mesh file name")
    arg_parser.add_argument("-start_x", type=int, default=0, help="Start index in x-direction")
    arg_parser.add_argument("-stop_x", type=int, default=-1, help="Stop index in x-direction")
    arg_parser.add_argument("-start_y", type=int, default=0, help="Start index in y-direction")
    arg_parser.add_argument("-stop_y", type=int, default=-1, help="Stop index in y-direction")
    arg_parser.add_argument("-loglevel", type=int, default=1, help="Verbosity")
    group = arg_parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-ratio", type=float, help="Edge ratio between original meshed raster and result")
    group.add_argument("-size", type=float, help="Number of edges in the generated surface mesh")
    res = arg_parser.parse_args(sys.argv[1:])
    logger.setLevel(res.loglevel)
    raster_coords = read_raster_file(filename=res.input,
                                     start_x=res.start_x,
                                     stop_x=res.stop_x,
                                     start_y=res.start_y,
                                     stop_y=res.stop_y)
    logger.debug(f"Original: {len(raster_coords)}")
    if res.ratio:
        pts, faces = triangulate_dem.lindstrom_turk_by_ratio(raster_coords, res.ratio)
    else:
        pts, faces = triangulate_dem.lindstrom_turk_by_ratio(raster_coords, res.size)
    logger.debug(f"Result: {len(raster_coords)}")
    save_to_ascii_off(pts, faces, res.output)
    logger.info(f"TIN written to {res.output}")
