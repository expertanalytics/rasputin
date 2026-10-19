"""tester.md section 3C for a cut run, which reads external input: a domain in
another CRS, a corrupted raster, and schema drift in the index.

`docs/increments/23-basin-scale.md`, "Output" (`MeshIndex`, frozen Pydantic)
and K10 (the core sees numbers in metres; CRS stays in Python).

- CRS: a domain in EPSG:4326 over the projected scene is moved by pyproj as
  today (15b) and then cut; the union passes conformity and both section 3D
  oracles against the domain moved by pyproj's own `Transformer` (DC3's
  area against the run's noded outline, and that outline within the
  noder's snap of the moved domain).
- A truncated DEM is a usage error in a cut run, as uncut (today: exit 2,
  "tifffile could not decode the file", no traceback), and leaves no
  `index.json`: a half-written run is not an index.
- Schema drift: `MeshIndex` refuses an unknown key, a missing `tolerance`,
  and a piece id that is not three integers (`ValidationError`); it is
  frozen, with `extra="forbid"` (confirmed under "Settled after 23c's red step (65e3990)", item 5).
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from pydantic import ValidationError
from shapely.geometry import Polygon

import pieces_fixtures as pf
from cli_helpers import USAGE
from geographic_fixtures import project_ring
from geotiff_fixtures import micro_tiff

WGS84 = "EPSG:4326"
UTM = "EPSG:25833"


@pytest.fixture
def mi() -> ModuleType:
    return importlib.import_module("tin_engine.io.mesh_index")


@pytest.fixture(scope="module")
def index_text(tmp_path_factory: pytest.TempPathFactory) -> str:
    where = tmp_path_factory.mktemp("index")
    dem = pf.write_dem(where / "dem.tif", pf.hills())
    domain = pf.write_domain(where / "sq.geojson", pf.SQUARE)
    done = pf.cut(
        where / "run",
        "--dem",
        str(dem),
        "--domain",
        str(domain),
        "--tolerance",
        "1",
        "--pieces",
        "4",
    )
    assert done.dir is not None
    return (done.dir / "index.json").read_text()


def test_a_domain_in_another_crs_is_moved_then_cut(tmp_path: Path) -> None:
    array = pf.hills()
    dem = pf.write_dem(tmp_path / "dem.tif", array)
    ring = project_ring(UTM, WGS84, pf.SQUARE)
    path = tmp_path / "wgs84.geojson"
    path.write_text(
        json.dumps(
            {
                "type": "Polygon",
                "coordinates": [[*map(list, ring), list(ring[0])]],
                "crs": {"type": "name", "properties": {"name": WGS84}},
            }
        )
    )
    done = pf.cut(
        tmp_path / "run",
        "--dem",
        str(dem),
        "--domain",
        str(path),
        "--tolerance",
        "1",
        "--pieces",
        "4",
    )
    moved = Polygon(project_ring(WGS84, UTM, ring))
    assert len(done.pieces) == 4
    assert pf.conformity_problems(done.pieces) == []
    assert pf.off_lattice_seams(done.pieces) == []
    u = pf.union(done.pieces)
    assert pf.validity_problems(u, moved, done.noded) == []
    assert pf.delaunay_violations(u.points[:, :2], u.triangles, u.constrained) == []
    found = pf.tolerance_check(u, array, moved, 1.0)
    assert (found.over, found.outside_mesh) == (0, 0), found


def test_a_truncated_dem_is_a_usage_error_and_no_index(tmp_path: Path) -> None:
    array = pf.hills()
    data = micro_tiff(array, scale=(pf.H, pf.H, 0.0)).getvalue()
    dem = tmp_path / "cut.tif"
    dem.write_bytes(data[: len(data) - array.nbytes // 2])
    domain = pf.write_domain(tmp_path / "sq.geojson", pf.SQUARE)
    done = pf.run(
        tmp_path / "run",
        "--dem",
        str(dem),
        "--domain",
        str(domain),
        "--tolerance",
        "1",
        "--pieces",
        "4",
    )
    assert done.code == USAGE, done.output
    assert "No such option" not in done.output
    assert "Traceback" not in done.output
    assert "could not decode" in done.output
    assert done.dir is None or not (done.dir / "index.json").exists()


class TestSchemaDrift:
    def test_the_runs_index_validates(self, mi: ModuleType, index_text: str) -> None:
        mi.MeshIndex.model_validate_json(index_text)

    @pytest.mark.parametrize(
        "drift",
        ["an unknown key", "no tolerance", "a two-part piece id", "a piece id of text"],
    )
    def test_refused(self, mi: ModuleType, index_text: str, drift: str) -> None:
        doc: dict[str, Any] = json.loads(index_text)
        if drift == "an unknown key":
            doc["pieces_total"] = 4
        elif drift == "no tolerance":
            del doc["tolerance"]
        elif drift == "a two-part piece id":
            doc["pieces"][0]["id"] = doc["pieces"][0]["id"][:2]
        else:
            doc["pieces"][0]["id"] = ["a", "b", "c"]
        with pytest.raises(ValidationError):
            mi.MeshIndex.model_validate_json(json.dumps(doc))

    def test_frozen(self, mi: ModuleType, index_text: str) -> None:
        model = mi.MeshIndex.model_validate_json(index_text)
        with pytest.raises(ValidationError):
            model.tolerance = 2.0
