import pytest
from tempfile import TemporaryDirectory
from pathlib import Path
from rasputin.triangulate_dem import PointVector, FaceVector
from rasputin.tin_repository import TinRepository

@pytest.fixture
def tin():
    pts = PointVector([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    faces = FaceVector([[0, 1, 2]])
    return pts, faces


def test_store_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(points=pts, faces=faces, uid=uid)
        assert uid in tr.content
        assert tr.content[uid]["num_points"] == len(pts)
        assert tr.content[uid]["num_faces"] == len(faces)
        assert "timestamp" in tr.content[uid]
        assert "projection" in tr.content[uid]


def test_store_and_load_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(points=pts, faces=faces, uid=uid)
        pts2, faces2 = tr.read(uid=uid)
        assert pts2 == pts
        assert faces2 == faces

def test_store_and_delete_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(points=pts, faces=faces, uid=uid)
        tr.delete(uid=uid)
        assert uid not in tr.content
        assert not (archive / f"{uid}.h5").exists()

