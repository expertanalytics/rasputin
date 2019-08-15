from typing import Dict, Any
import numpy as np
from datetime import datetime
from h5py import File
from lxml import etree
from pathlib import Path
from rasputin.mesh import Mesh
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
                mesh = Mesh.from_points_and_faces(points=pts, faces=faces)
                geometries[name] = Geometry(mesh=mesh,
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

    def save(self, *, uid: str, geometries: Dict[str, Geometry]) -> None:
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
            root_grp = archive.create_group("tins")
            for name, geom in geometries.items():
                grp = root_grp.create_group(name)
                points, faces = geom.points, geom.faces
                grid = etree.SubElement(domain, "Grid", Name=name)
                x_geom = etree.SubElement(grid,
                                        "Geometry",
                                        GeometryType="XYZ")
                pts_elm = etree.SubElement(x_geom,
                                           "DataItem",
                                           Format="HDF",
                                           DataType="Float",
                                           Precision="8",
                                           Dimensions=f"{points.shape[0]} {points.shape[1]}")
                pts_elm.text = f"{h5_base}:/tins/{name}/points"
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
                attr_elm.text = f"{h5_base}:/tins/{name}/face_color"
                etree.SubElement(pts_elm,
                                 "Information",
                                 Name="projection",
                                 Value=geom.projection.definition_string())
                f_elm.text = f"{h5_base}:/tins/{name}/faces"
                h5_points = grp.create_dataset(name="points", data=points, dtype="d")
                h5_points.attrs["projection"] = geom.projection.definition_string()
                h5_faces = grp.create_dataset(name="faces", data=faces, dtype="i")
                h5_faces.attrs["color"] = geom.base_color
                color_arr = np.empty((len(faces), 3), dtype="float")
                color_arr[:, 0] = geom.base_color[0]
                color_arr[:, 1] = geom.base_color[1]
                color_arr[:, 2] = geom.base_color[2]
                grp.create_dataset(name="face_color", data=color_arr)
        tree.write(str(xdmf_filename), encoding="utf-8")


    def delete(self, uid: str) -> None:
        if uid in self.content:
            (self.path / f"{uid}.h5").unlink()
            (self.path / f"{uid}.xdmf").unlink()
