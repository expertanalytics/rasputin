import pytest
from pathlib import Path
import os
import pyproj
from rasputin.reader import RasterRepository

if "RASPUTIN_DATA_DIR" in os.environ:
    data_dir = Path(os.environ["RASPUTIN_DATA_DIR"])
else:
    data_dir = Path(os.environ["HOME"]) /"projects" / "rasputin_data" / "dem_archive"


@pytest.mark.skip
def test_get():
    #x0 = 8.530918
    x0 = 8.54758671814368
    y0 = 60.898468
    repo = RasterRepository(directory=data_dir)
    input_coordinate_system = pyproj.Proj(init="EPSG:4326").definition_string()
    target_coordinate_system = pyproj.Proj(init="EPSG:32633").definition_string()
    result = repo.read(x=x0,
                       y=y0,
                       dx=210,
                       dy=210,
                       input_coordinate_system=input_coordinate_system,
                       target_coordinate_system=target_coordinate_system)
    assert result



def test_run():
    from descartes import PolygonPatch
    import matplotlib.pyplot as plt
    from pathlib import Path
    import os
    import pyproj
    from rasputin.reader import RasterRepository
    if "RASPUTIN_DATA_DIR" in os.environ:
        data_dir = Path(os.environ["RASPUTIN_DATA_DIR"])
    else:
        data_dir = Path(os.environ["HOME"]) /"projects" / "rasputin_data" / "dem_archive"

    x0 = 8.54550
    y0 = 60.9
    repo = RasterRepository(directory=data_dir)
    input_coordinate_system = pyproj.Proj(init="EPSG:4326").definition_string()
    target_coordinate_system = pyproj.Proj(init="EPSG:32633").definition_string()
    result = repo.read(x=x0,
                       y=y0,
                       dx=50,
                       dy=5.0,
                       input_coordinate_system=input_coordinate_system,
                       target_coordinate_system=target_coordinate_system)

    print([r for r in result])

    xbound = repo.shapes["target_bbox"].bounds[0::2]
    ybound = repo.shapes["target_bbox"].bounds[1::2]

    f1 = "/Users/skavhaug/projects/rasputin_data/dem_archive/6701_1_10m_z33.tif"
    f2 = "/Users/skavhaug/projects/rasputin_data/dem_archive/6701_4_10m_z33.tif"
    i1 = repo.shapes.get(f"{f1}_intersection", None)
    i2 = repo.shapes.get(f"{f2}_intersection", None)

    r1 = repo.shapes.get(f"{f1}_remainding_bbox")
    r2 = repo.shapes.get(f"{f2}_remainding_bbox")

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if i1:
        ax.add_patch(PolygonPatch(i1, fc="red", alpha=0.5))
    if i2:
        ax.add_patch(PolygonPatch(i2, fc="blue", alpha=0.5))
    ax.set_xbound(*xbound)
    ax.set_ybound(*ybound)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if r1:
        ax.add_patch(PolygonPatch(r1, fc="green", alpha=0.5))
    if r2:
        ax.add_patch(PolygonPatch(r2, fc="yellow", alpha=0.5))
    ax.set_xbound(*xbound)
    ax.set_ybound(*ybound)

    plt.show()
