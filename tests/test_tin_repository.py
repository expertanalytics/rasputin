import pytest
from tempfile import TemporaryDirectory
from pathlib import Path
import numpy as np
from rasputin.tin_repository import TinRepository
from rasputin.mesh import Mesh
from rasputin.geometry import Geometry
from pyproj import CRS

@pytest.fixture
def tin():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype='d')
    faces = np.array([[0, 1, 2]], dtype='i')
    return pts, faces


def test_store_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    crs = CRS.from_epsg(32633)
    geom_tag = "test_geom"
    mesh = Mesh.from_points_and_faces(points=pts, faces=faces)
    geom = {geom_tag: Geometry(mesh=mesh,
                               crs=crs,
                               base_color=(0, 0, 0),
                               material="dummy2")}
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(uid=uid, geometries=geom)
        assert uid in tr.content
        print(tr.content[uid])
        geom_info = tr.content[uid]["tins"][geom_tag]
        assert geom_info["num_points"] == len(pts)
        assert geom_info["num_faces"] == len(faces)
        assert "timestamp" in tr.content[uid]
        assert "projection" in geom_info


def test_store_and_load_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    crs = CRS(32633)
    mesh = Mesh.from_points_and_faces(points=pts, faces=faces)
    geom = {"test_geom": Geometry(mesh=mesh,
                                  crs=crs,
                                  base_color=(0,0,0),
                                  material="dummy2")}
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(uid=uid, geometries=geom)
        geom = tr.read(uid=uid)["test_geom"]
        assert np.allclose(geom.points, pts)
        assert np.allclose(geom.faces, faces)


def test_store_and_delete_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    crs = CRS.from_epsg(32633)
    mesh = Mesh.from_points_and_faces(points=pts, faces=faces)
    geom = {"test_geom": Geometry(mesh=mesh,
                                  crs=crs,
                                  base_color=(0,0,0),
                                  material="dummy2")}
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(uid=uid, geometries=geom)
        tr.delete(uid=uid)
        assert uid not in tr.content
        assert not (archive / f"{uid}.h5").exists()

