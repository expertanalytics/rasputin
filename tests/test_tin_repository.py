import pytest
from tempfile import TemporaryDirectory
from pathlib import Path
from rasputin.triangulate_dem import point3_vector, face_vector
from rasputin.tin_repository import TinRepository
from rasputin.geometry import Geometry
from pyproj import Proj

@pytest.fixture
def tin():
    pts = point3_vector([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    faces = face_vector([[0, 1, 2]])
    return pts, faces


def test_store_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    proj = Proj(init="EPSG:32633")
    geom_tag = "test_geom"
    geom = {geom_tag: Geometry(points=pts,
                                  faces=faces,
                                  projection=proj,
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
    proj = Proj(init="EPSG:32633")
    geom = {"test_geom": Geometry(points=pts,
                                  faces=faces,
                                  projection=proj,
                                  base_color=(0,0,0),
                                  material="dummy2")}
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(uid=uid, geometries=geom)
        geom = tr.read(uid=uid)["test_geom"]
        assert geom.points == pts
        assert geom.faces == faces

def test_store_and_delete_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    proj = Proj(init="EPSG:32633")
    geom = {"test_geom": Geometry(points=pts,
                                  faces=faces,
                                  projection=proj,
                                  base_color=(0,0,0),
                                  material="dummy2")}
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(uid=uid, geometries=geom)
        tr.delete(uid=uid)
        assert uid not in tr.content
        assert not (archive / f"{uid}.h5").exists()

