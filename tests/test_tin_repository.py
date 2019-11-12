import pytest
from tempfile import TemporaryDirectory
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
from rasputin.tin_repository import TinRepository, ShadeRepository
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
    mesh = Mesh.from_points_and_faces(points=pts, faces=faces, proj4_str=crs.to_proj4())
    geom = Geometry(mesh=mesh,
                    crs=crs,
                    base_color=(0, 0, 0),
                    material="dummy2")
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(uid=uid, geometry=geom)
        assert uid in tr.content
        geom_info = tr.info(uid=uid)
        assert geom_info["tin"]["num_points"] == len(pts)
        assert geom_info["tin"]["num_faces"] == len(faces)
        assert "timestamp" in tr.content[uid]
        assert "projection" in geom_info["tin"]


def test_store_and_load_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    crs = CRS(32633)
    mesh = Mesh.from_points_and_faces(points=pts, faces=faces, proj4_str=crs.to_proj4())
    geom = Geometry(mesh=mesh,
                    crs=crs,
                    base_color=(0,0,0),
                    material="dummy2")
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(uid=uid, geometry=geom)
        geom = tr.read(uid=uid)
        assert np.allclose(geom.points, pts)
        assert np.allclose(geom.faces, faces)


def test_store_and_delete_tin(tin):
    pts, faces = tin
    uid = "test_tin"
    crs = CRS.from_epsg(32633)
    mesh = Mesh.from_points_and_faces(points=pts, faces=faces, proj4_str=crs.to_proj4())
    geom = Geometry(mesh=mesh,
                    crs=crs,
                    base_color=(0,0,0),
                    material="dummy2")
    with TemporaryDirectory() as directory:
        archive = Path(directory)
        tr = TinRepository(path=archive)
        tr.save(uid=uid, geometry=geom)
        tr.delete(uid=uid)
        assert uid not in tr.content
        assert not (archive / f"{uid}.h5").exists()


def test_store_shade(tin):
    pts, faces = tin
    uid = "test_tin"
    shade_uid = "test_tin_shade"
    crs = CRS(32633)
    mesh = Mesh.from_points_and_faces(points=pts, faces=faces, proj4_str=crs.to_proj4())
    geom = Geometry(mesh=mesh,
                    crs=crs,
                    base_color=(0,0,0),
                    material="dummy2")

    with TemporaryDirectory() as directory:
        tin_archive = Path(directory) / "tin_archive"
        shade_archive_path = Path(directory) / "shade_archive"
        tin_archive.mkdir()
        shade_archive_path.mkdir()
        tr = TinRepository(path=tin_archive)
        tr.save(uid=uid, geometry=geom)
        mesh = tr.read(uid=uid).mesh
        start_year = 2019
        start_month = 10
        start_day = 25
        end_year = 2019
        end_month = 10
        end_day = 26
        start_time = datetime(start_year, start_month, start_day)
        end_time = datetime(end_year, end_month, end_day)
        frequency = 3600
        dt = timedelta(seconds=frequency)
        t = start_time
        shade_repo = ShadeRepository(path=shade_archive_path)
        with shade_repo.open(tin_repo=tr,
                             tin_uid=uid,
                             shade_uid=shade_uid,
                             overwrite=True) as shade_writer:
            while t <= end_time:
                shade = mesh.shade(t.timestamp())
                shade_writer.save(t.timestamp(), shade)
                t += dt
        
        info = shade_repo.info(uid=shade_uid)
        assert info
        assert info["tin_uid"] == uid
        assert info.get("timestamps", False)
        assert isinstance(info["timestamps"], list)

