from typing import Tuple, Dict, Any
import numpy as np
from datetime import datetime
from pathlib import Path
from json import loads, dumps
from rasputin.triangulate_dem import PointVector, FaceVector


class TinRepository:

    def __init__(self, *, path: Path) -> None:
        self.path = path

    def read(self, *, uid: str) -> Tuple[PointVector, FaceVector]:
        data = np.load(self.path / f"{uid}.npz")
        pts = data["points"]
        faces = data["faces"]
        return PointVector(pts.tolist()), FaceVector(faces.tolist())

    @property
    def content(self) -> Dict[str, Dict[str, Any]]:
        files = self.path.glob("*.npz")
        meta_info = dict()
        for f in files:
            with f.with_suffix(".meta").open("r") as meta:
                meta_info[f.stem] = loads(meta.read())
        return meta_info

    def save(self, *, uid: str, points: PointVector, faces: FaceVector):
        if (self.path / f"{uid}.npz").exists():
            raise RuntimeError(f"Archive already has a data set with uid {uid}.")
        if (self.path / f"{uid}.meta").exists():
            raise RuntimeError(f"Archive already has a data set with uid {uid}.")
        np.savez(self.path / f"{uid}.npz", points=points, faces=faces)
        with (self.path / f"{uid}.meta").open("w") as meta:
            meta.write(dumps(dict(num_points=len(points),
                                  num_faces=len(faces),
                                  timestamp=datetime.utcnow().timestamp())))

    def delete(self, uid: str) -> None:
        if uid in self.content:
            (self.path / f"{uid}.npz").unlink()
            (self.path / f"{uid}.meta").unlink()
