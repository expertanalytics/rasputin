from copy import deepcopy
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


def _indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            _indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


class TinRepository:

    def __init__(self, *, path: Path) -> None:
        self.path = path
        self.xdmf_ext = ".xdmf"
        self.h5_ext = ".h5"

    def info(self, *, uid: str) -> Dict[str, Any]:
        info = dict()
        with File(self.path / f"{uid}{self.h5_ext}", "r") as archive:
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
                land_covers = info_grp["land_covers"][:]
                info["info"]["land_covers"] = [(v, n.decode("utf-8")) for (v, n) in land_covers]
            return info

    def read(self, *, uid: str) -> Geometry:
        filename = self.path / f"{uid}.h5"
        if not filename.exists():
            raise FileNotFoundError(f"Mesh with uid '{uid}' not found in tin archive {self.path}.")
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
        files = self.path.glob(f"*{self.h5_ext}")
        meta_info = dict()
        for f in files:
            meta_info[f.stem] = self.info(uid=f.stem)
        return meta_info

    def save(self, *,
             uid: str,
             geometry: Geometry,
             land_cover_repository: LandCoverRepository,
             face_fields: Optional[Dict[str, np.ndarray]] = None) -> None:
        xdmf_filename = self.path / f"{uid}{self.xdmf_ext}"
        h5_base = f"{uid}{self.h5_ext}"
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

            etree.SubElement(pts_elm,
                             "Information",
                             Name="projection",
                             Value=geometry.crs.to_proj4())
            f_elm.text = f"{h5_base}:/{tin_grp_name}/faces"
            h5_points = tin_grp.create_dataset(name="points", data=points, dtype="d")
            h5_points.attrs["projection"] = geometry.crs.to_proj4()
            tin_grp.create_dataset(name="faces", data=faces, dtype="i")
            info_grp = archive.create_group(info_grp_name)
            info_grp.attrs["land_cover_type"] = land_cover_repository.__class__.__name__
            land_covers = [(v.value, v.name.encode("utf-8")) for v in land_cover_repository.land_cover_type]
            if face_fields:
                face_grp = tin_grp.create_group("face_fields")
                for field_name, field in face_fields.items():
                    face_grp.create_dataset(name=field_name, data=field)
                    # Filter land cover fields by the ones found in dataset
                    if field_name == "cover_type":
                        used_fields = set(field)
                        land_covers = [(v, n) for (v, n) in land_covers if v in used_fields]
                    atype = "Vector" if len(field.shape) == 2 and field.shape[1] == 3 else "Scalar"
                    dtype = "Int" if np.issubdtype(field.dtype, np.integer) else "Float"
                    precision = "4" if dtype == "Int" else "8"
                    dims = f"{' '.join([str(d) for d in field.shape])}"
                    attr = etree.SubElement(grid,
                                            "Attribute",
                                            Name=field_name,
                                            Center="Cell",
                                            AttributeType=atype)
                    attr_elm = etree.SubElement(attr,
                                                "DataItem",
                                                Format="HDF",
                                                Precision=precision,
                                                DataType=dtype,
                                                Dimensions=dims)
                    attr_elm.text = f"{h5_base}:/{tin_grp_name}/face_fields/{field_name}"
            info_grp.create_dataset(name="land_covers",
                                    data=np.array(land_covers, dtype=[("value", "i4"), ("name", "S100")]))
        _indent(domain, level=1)
        tree.write(str(xdmf_filename), pretty_print=True, encoding="utf-8")


    def delete(self, uid: str) -> None:
        if uid in self.content:
            (self.path / f"{uid}.h5").unlink()
            (self.path / f"{uid}.xdmf").unlink()


class ShadeRepository:

    def __init__(self, *,
                 tin_repo: TinRepository,
                 tin_uid: str,
                 path: Path,
                 shade_uid: str,
                 overwrite: bool) -> None:
        self.tin_repo = tin_repo
        self.tin_uid = tin_uid
        self.path = path
        self.shade_uid = shade_uid
        self.h5_fh: Optional[File] = None
        self.tree = None
        self.tgrid = None
        self.ncells = ""
        self.overwrite = overwrite

    def __enter__(self):
        assert self.path.is_dir()
        assert self.h5_fh is None
        xdmf_fullpath = (self.path / f"{self.shade_uid}.xdmf")
        hdf5_fullpath = (self.path / f"{self.shade_uid}.h5")
        if not self.overwrite:
            assert not xdmf_fullpath.exists()
            assert not hdf5_fullpath.exists()
        else:
            if xdmf_fullpath.exists():
                xdmf_fullpath.unlink()
            if hdf5_fullpath.exists():
                hdf5_fullpath.unlink()
        self.h5_fh = File(hdf5_fullpath, "w")
        self.h5_fh.create_group(name=self.tin_uid)
        self.tree = etree.parse(str(self.tin_repo.path / f"{self.tin_uid}{self.tin_repo.xdmf_ext}"))
        domain = self.tree.find("Domain")
        grid = domain.find("Grid")
        elm = grid.find("Geometry").find("DataItem")
        elm.text = str(self.tin_repo.path / elm.text)
        elm = grid.find("Topology").find("DataItem")
        elm.text = str(self.tin_repo.path / elm.text)
        self.ncells = [c for c in grid if c.tag == "Topology"][0].get("NumberOfElements")
        template_grid = deepcopy(grid)
        domain.remove(grid)
        collection = etree.SubElement(domain,
                                      "Grid",
                                      Name="shadow_times",
                                      GridType="Collection",
                                      CollectionType="Temporal")
        for child in template_grid.iterchildren():
            if child.tag == "Attribute":
                template_grid.remove(child)

        self.tgrid = template_grid
        return self

    def __exit__(self, type, value, traceback):
        domain = self.tree.find("Domain")
        _indent(domain, level=1)
        xdmf_fullpath = (self.path / f"{self.shade_uid}.xdmf")
        self.tree.write(str(xdmf_fullpath), pretty_print=True, encoding="utf-8")
        self.h5_fh.close()
        self.xdmf_fh = None
        self.h5_fh = None
        self.tree = self.domain = None

    def save(self, timestamp: float, data: np.ndarray) -> None:
        assert self.h5_fh is not None
        grid = deepcopy(self.tgrid)
        h5_base = (self.path / f"{self.shade_uid}.h5")
        etree.SubElement(grid,
                         "Time",
                         Value=f"{timestamp:.3f}")
        attr = etree.SubElement(grid,
                                "Attribute",
                                Name="shade",
                                Center="Cell",
                                AttributeType="Scalar")
        attr_elm = etree.SubElement(attr,
                                    "DataItem",
                                    Format="HDF",
                                    Precision="4",
                                    DataType="Int",
                                    Dimensions=self.ncells)
        shade_name = f"{timestamp:.4f}"
        attr_elm.text = f"{h5_base}:/{self.tin_uid}/{shade_name}/"
        collection = self.tree.find("Domain").find("Grid")
        collection.append(grid)
        self.h5_fh[self.tin_uid].create_dataset(name=f"{timestamp:.4f}", data=data, dtype='i')


