from pathlib import Path
import os
import numpy as np
import pyproj
from shapely.geometry import Polygon
from rasputin.gml_repository import GMLRepository, LandCoverType
from descartes import PolygonPatch
import matplotlib.pyplot as plt


def test_gml_repository():
    if "RASPUTIN_DATA_DIR" not in os.environ:
        raise RuntimeError("Please set RASPUTIN_DATA_DIR")
    path = Path(os.environ["RASPUTIN_DATA_DIR"]) / "corine"
    input_coordinate_system = pyproj.Proj(init="EPSG:4326")
    target_coordinate_system = pyproj.Proj(init="EPSG:32633")

    x = np.array([8.5, 8.55, 8.55, 8.5])
    y = np.array([60.55, 60.55, 60.45, 60.45])

    x, y = pyproj.transform(input_coordinate_system, target_coordinate_system, x, y)
    domain = Polygon(shell=list(zip(x, y)))

    repos = GMLRepository(path=path)
    response = repos.read(domain=domain)
    fig = plt.figure()
    ax = fig.gca()
    for key in response:
        color = [c/255 for c in LandCoverType.color(land_cover_type=key)]
        for p in response[key]:
            ax.add_patch(PolygonPatch(p, fc=color, alpha=0.5))
    ax.set_xbound(min(x), max(x))
    ax.set_ybound(min(y), max(y))
    plt.show()
