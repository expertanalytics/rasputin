from typing import Dict, Any, Optional
import numpy as np
from datetime import datetime
from h5py import File
from lxml import etree
from pathlib import Path
from pyproj import CRS
from rasputin.geometry import Geometry
from rasputin.land_cover_repository import LandCoverRepository
from rasputin.mesh import Mesh


class TinRepository:

    def __init__(self, *, path: Path) -> None:
        self.path = path

    def info(self, *, uid: str) -> Dict[str, Any]:
        info = dict()
        with File(self.path / f"{uid}.h5", "r") as archive:
            info["timestamp"] = archive.attrs["timestamp"]
            terrain_grp = archive["tin"]
            info_grp = archive.get("information", None)
            info["tin"] = dict(projection=terrain_grp["points"].attrs["projection"],
                               num_points=terrain_grp["points"].shape[0],
                               num_faces=terrain_grp["faces"].shape[0])
            if "face_fields" in terrain_grp:
                face_grp = terrain_grp["face_fields"]
                ffs = dict()
                for f_name in face_grp:
                    ffs[f_name] = face_grp[f_name].shape
                info["tin"]["face_fields"] = ffs
            if info_grp:
                info["info"] = dict()
                info["info"]["land_cover_type"] = info_grp.attrs["land_cover_type"]
                info["info"]["land_covers"] = info_grp["land_covers"].shape
            return info

    def read(self, *, uid: str) -> Geometry:
        filename = self.path / f"{uid}.h5"
        if not filename.exists():
            raise FileNotFoundError(f"File {filename.absolute()} not found.")
        with File(filename, "r") as archive:
            tin_grp_name = "tin"
            tin_group = archive[tin_grp_name]
            pts = tin_group["points"][:]
            projection = tin_group["points"].attrs["projection"]
            faces = tin_group["faces"][:]
            mesh = Mesh.from_points_and_faces(points=pts, faces=faces, proj4_str=projection)
            geometry = Geometry(mesh=mesh,
                                crs=CRS.from_proj4(projection))
            if "face_fields" in tin_group and "cover_color" in tin_group["face_fields"]:
                geometry.colors = tin_group["face_fields"]["cover_color"][:]
            return geometry

    def extract(self, * uid: str, face_id: int) -> Geometry:
        filename = self.path / f"{uid}.h5"
        if not filename.exists():
            raise FileNotFoundError(f"File {filename.absolute()} not found.")
        with File(filename, "r") as archive:
            tin_grp_name = "tin"
            tin_group = archive[tin_grp_name]
            if "face_fields" not in tin_group:
                raise IOError("No face fields in dataset.")
            face_group = tin_group["face_fields"]
            if "cover_type" not in face_group:
                raise IOError("No cover type information in dataset.")
            if "cover_color" not in face_group:
                raise IOError("No face color information in dataset.")
            pts = tin_group["points"][:]
            projection = tin_group["points"].attrs["projection"]
            faces = tin_group["faces"][:]
            mesh = Mesh.from_points_and_faces(points=pts, faces=faces, proj4_str=projection)
            indices = [i for (i, f) in enumerate(tin_group["face_fields"]["cover_type"]) if f == face_id]
            if not indices:
                raise IOError(f"face_id {face_id} not found in dataset")
            base_color = tin_group["face_fields"]["cover_color"][indices[0]]
            mesh = mesh.extract_sub_mesh(np.asarray(indices))
            geometry = Geometry(mesh=mesh,
                                crs=CRS.from_proj4(projection),
                                base_color=base_color)
            return geometry

    @property
    def content(self) -> Dict[str, Dict[str, Any]]:
        files = self.path.glob("*.h5")
        meta_info = dict()
        for f in files:
            meta_info[f.stem] = self.info(uid=f.stem)
        return meta_info

    def save(self, *,
             uid: str,
             geometry: Geometry,
             land_cover_repository: Optional[LandCoverRepository] = None,
             face_fields: Optional[Dict[str, np.ndarray]] = None) -> None:
        xdmf_filename = self.path / f"{uid}.xdmf"
        h5_base = f"{uid}.h5"
        h5_filename = self.path / h5_base
        if h5_filename.exists():
            raise FileExistsError(f"Archive already has a data set with uid {uid}.")
        root = etree.XML('''\
<?xml version="1.0" ?> 
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []> 
<Xdmf Version="3.0"> 
</Xdmf>''')
        tree = etree.ElementTree(root)
        domain = etree.SubElement(root, "Domain")

        with File(h5_filename, "w") as archive:
            timestamp = datetime.utcnow().timestamp()
            archive.attrs["timestamp"] = timestamp
            etree.SubElement(root,
                             "Information",
                             Name="timestamp",
                             Value=str(timestamp))
            tin_grp_name = "tin"
            info_grp_name = "information"
            tin_grp = archive.create_group(tin_grp_name)
            points, faces = geometry.points, geometry.faces
            grid = etree.SubElement(domain, "Grid", Name=tin_grp_name)
            x_geom = etree.SubElement(grid,
                                    "Geometry",
                                    GeometryType="XYZ")
            pts_elm = etree.SubElement(x_geom,
                                       "DataItem",
                                       Format="HDF",
                                       DataType="Float",
                                       Precision="8",
                                       Dimensions=f"{points.shape[0]} {points.shape[1]}")
            pts_elm.text = f"{h5_base}:/{tin_grp_name}/points"
            topo = etree.SubElement(grid,
                                    "Topology",
                                    NumberOfElements=str(faces.shape[0]),
                                    TopologyType="Triangle")
            f_elm = etree.SubElement(topo,
                                     "DataItem",
                                     Format="HDF",
                                     Precision="4",
                                     DataType="Int",
                                     Dimensions=f"{faces.shape[0]} {faces.shape[1]}")
            attr = etree.SubElement(grid,
                                    "Attribute",
                                    Name="Color",
                                    Center="Cell",
                                    AttributeType="Vector")
            attr_elm = etree.SubElement(attr,
                                        "DataItem",
                                        Format="HDF",
                                        Precision="8",
                                        DataType="Float",
                                        Dimensions=f"{faces.shape[0]} 3")
            attr_elm.text = f"{h5_base}:/{tin_grp_name}/face_color"
            etree.SubElement(pts_elm,
                             "Information",
                             Name="projection",
                             Value=geometry.crs.to_proj4())
            f_elm.text = f"{h5_base}:/{tin_grp_name}/faces"
            h5_points = tin_grp.create_dataset(name="points", data=points, dtype="d")
            h5_points.attrs["projection"] = geometry.crs.to_proj4()
            tin_grp.create_dataset(name="faces", data=faces, dtype="i")
            if face_fields:
                face_grp = tin_grp.create_group("face_fields")
                for field_name, field in face_fields.items():
                    face_grp.create_dataset(name=field_name, data=field)
            if land_cover_repository is not None:
                info_grp = archive.create_group(info_grp_name)
                info_grp.attrs["land_cover_type"] = land_cover_repository.__class__.__name__
                land_covers = [(v.value, v.name.encode("utf-8")) for v in land_cover_repository.land_cover_type]
                land_covers = np.array(land_covers, dtype=[("value", "i4"), ("name", "S100")])
                info_grp.create_dataset(name="land_covers", data=land_covers)
        tree.write(str(xdmf_filename), encoding="utf-8")


    def delete(self, uid: str) -> None:
        if uid in self.content:
            (self.path / f"{uid}.h5").unlink()
            (self.path / f"{uid}.xdmf").unlink()
