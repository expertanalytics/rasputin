from typing import Tuple, Dict, Any
import numpy as np
from datetime import datetime
from h5py import File
from pathlib import Path
from rasputin.triangulate_dem import PointVector, FaceVector


class TinRepository:

    def __init__(self, *, path: Path) -> None:
        self.path = path

    def read(self, *, uid: str) -> Tuple[PointVector, FaceVector]:
        filename = self.path / f"{uid}.h5"
        if not filename.exists():
            raise FileNotFoundError(f"File {filename.absolute()} not found.")
        with File(filename, "r") as archive:
            pts = archive["tin"]["points"][:]
            faces = archive["tin"]["faces"][:]
        return PointVector(pts.tolist()), FaceVector(faces.tolist())

    @property
    def content(self) -> Dict[str, Dict[str, Any]]:
        files = self.path.glob("*.h5")
        meta_info = dict()
        for f in files:
            with File(f, "r") as archive:
                meta_info[f.stem] = dict(timestamp=archive.attrs["timestamp"],
                                         projection=archive["tin"]["points"].attrs["projection"],
                                         num_points=archive["tin"]["points"].shape[0],
                                         num_faces=archive["tin"]["faces"].shape[0])
        return meta_info

    def save(self, *, uid: str, points: PointVector, faces: FaceVector):
        filename = self.path / f"{uid}.h5"
        if filename.exists():
            raise FileExistsError(f"Archive already has a data set with uid {uid}.")
        with File(filename, "w") as archive:
            grp = archive.create_group("tin")
            obj = grp.create_dataset(name="points", data=np.asarray(points), dtype="d")
            obj.attrs["projection"] = "NotSet"
            grp.create_dataset(name="faces", data=np.asarray(faces), dtype="i")
            archive.attrs["timestamp"] = datetime.utcnow().timestamp()

    def delete(self, uid: str) -> None:
        if uid in self.content:
            (self.path / f"{uid}.h5").unlink()
