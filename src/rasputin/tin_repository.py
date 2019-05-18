from typing import Tuple, Dict, Any
import numpy as np
from datetime import datetime
from pathlib import Path
from json import loads, dumps
from rasputin.triangulate_dem import point3_vector, face_vector


class TinRepository:

    def __init__(self, *, path: Path) -> None:
        self.path = path

    def read(self, *, uid: str) -> Tuple[point3_vector, face_vector]:
        data = np.load(self.path / f"{uid}.npz")
        pts = data["points"]
        faces = data["faces"]
        return point3_vector(pts.tolist()), face_vector(faces.tolist())

    def info(self, *, uid: str) -> Dict[str, Any]:
        with (self.path / f"{uid}.meta").open("r") as meta:
            return loads(meta.read())

    @property
    def content(self) -> Dict[str, Dict[str, Any]]:
        files = self.path.glob("*.npz")
        meta_info = dict()
        for f in files:
            meta_info[f.stem] = self.info(uid=f.stem)
        return meta_info

    def save(self,
             *,
             uid: str,
             points: point3_vector,
             faces: face_vector,
             projection: str) -> None:
        if (self.path / f"{uid}.npz").exists():
            raise RuntimeError(f"Archive already has a data set with uid {uid}.")
        if (self.path / f"{uid}.meta").exists():
            raise RuntimeError(f"Archive already has a data set with uid {uid}.")
        np.savez(self.path / f"{uid}.npz", points=points, faces=faces)
        with (self.path / f"{uid}.meta").open("w") as meta:
            meta.write(dumps(dict(num_points=len(points),
                                  num_faces=len(faces),
                                  timestamp=datetime.utcnow().timestamp(),
                                  projection=projection)))

    def delete(self, uid: str) -> None:
        if uid in self.content:
            (self.path / f"{uid}.npz").unlink()
            (self.path / f"{uid}.meta").unlink()
