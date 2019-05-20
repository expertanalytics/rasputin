from typing import Dict, Any
import numpy as np
from datetime import datetime
from h5py import File
from pathlib import Path
from rasputin.triangulate_dem import point3_vector, face_vector, consolidate
from rasputin.geometry import Geometry


class TinRepository:

    def __init__(self, *, path: Path) -> None:
        self.path = path

    def info(self, *, uid: str) -> Dict[str, Any]:
        info = dict(tins=dict())
        with File(self.path / f"{uid}.h5", "r") as archive:
            info["timestamp"] = archive.attrs["timestamp"]
            for name in archive["tins"].keys():
                grp = archive["tins"][name]
                info["tins"][name] = dict(projection=grp["points"].attrs["projection"],
                                          num_points=grp["points"].shape[0],
                                          num_faces=grp["faces"].shape[0],
                                          color=grp["faces"].attrs["color"])
            return info

    def read(self, *, uid: str) -> Dict[str, Geometry]:
        filename = self.path / f"{uid}.h5"
        if not filename.exists():
            raise FileNotFoundError(f"File {filename.absolute()} not found.")
        geometries = dict()
        with File(filename, "r") as archive:
            root_group = archive["tins"]
            for name in root_group.keys():
                group = root_group[name]
                pts = group["points"][:]
                projection = group["points"].attrs["projection"]
                faces = group["faces"][:]
                color = group["faces"].attrs["color"]
                geometries[name] = Geometry(points=point3_vector(pts.tolist()),
                                            faces=face_vector(faces.tolist()),
                                            projection=projection,
                                            base_color=color,
                                            material=None)
        return geometries

    @property
    def content(self) -> Dict[str, Dict[str, Any]]:
        files = self.path.glob("*.h5")
        meta_info = dict()
        for f in files:
            meta_info[f.stem] = self.info(uid=f.stem)
        return meta_info

    def save(self, *, uid: str, geometries: Dict[str, Geometry], consolidate_mesh: bool=False) -> None:
        filename = self.path / f"{uid}.h5"
        if filename.exists():
            raise FileExistsError(f"Archive already has a data set with uid {uid}.")
        with File(filename, "w") as archive:
            archive.attrs["timestamp"] = datetime.utcnow().timestamp()
            root_grp = archive.create_group("tins")
            for name, geom in geometries.items():
                grp = root_grp.create_group(name)
                if consolidate_mesh:
                    points, faces = consolidate(geom.points, geom.faces)
                else:
                    points, faces = geom.points, geom.faces
                h5_points = grp.create_dataset(name="points", data=np.asarray(points), dtype="d")
                h5_points.attrs["projection"] = geom.projection
                h5_faces = grp.create_dataset(name="faces", data=np.asarray(faces), dtype="i")
                h5_faces.attrs["color"] = geom.base_color

    def delete(self, uid: str) -> None:
        if uid in self.content:
            (self.path / f"{uid}.h5").unlink()
