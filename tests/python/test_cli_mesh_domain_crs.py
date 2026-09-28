"""`rasputin mesh --dem ... --domain PATH` with the domain in its own CRS (increment 15b).

`docs/increments/15-dem-mosaic.md` R9, R11 and "Tests for @tester" (15b), and
the 15b Acceptance: a catchment polygon in EPSG:4326 over a 25833 mosaic
meshes, with `domain_crs` and `domain_transform` recorded. Pinned by this suite
(see "Pinned by the red suite (15b)"):

- `domain_crs` is `EPSG:n` when pyproj finds an exact EPSG code for the
  domain's CRS, and otherwise text pyproj parses back to that CRS.
- `domain_transform` is pyproj's `description` of the `always_xy` transformer
  from the domain's CRS to the DEM's; for a domain already in the DEM's CRS it
  says `none` (any case), and no transformer is made.
- A domain already in the DEM's CRS meshes bit for bit as increment 16 meshed
  it: the digests below were RECORDED FROM THE TREE BEFORE 15b, at `d34d79d`
  (15a, whose single-file domain path is 16's), before any 15b production
  change, with `_core` rebuilt from that tree. Both are unchanged when the
  same run adds `--bbox` at the domain's bounds, which is what R6 makes the
  domain do, so a window cut to the domain may not move them either. No
  commit may update them to agree with new code.
- `--bbox` with `--domain` is a usage error naming both.
- `dem_seams` is written on the domain path as without a domain: the
  disagreeing pair with a domain in EPSG:4326, `none` when the overlaps agree
  (test amendment after the 15b review).
- A transformed domain refused for its extent (outside the DEM, in no tile, or
  with no image in the DEM's CRS) is a usage error naming the domain's CRS and
  the DEM's EPSG code, and writes nothing.

The axis-order test is able to fail: it patches `pyproj.Transformer.from_crs`
to drop `always_xy`, and the same run that meshes unpatched must then be
refused by the extent check (R9, "Degeneracy policy": axis order).

HOW THIS FILE GOES RED: it imports nothing new. Before 15b a domain in
another CRS is refused as a CRS mismatch, `domain_crs` is not written, and
`--bbox` with `--domain` is accepted, so those tests fail on their assertions.
Six are guards that pass before 15b and must stay green after it: the two
same-CRS digests, the two unknown `--domain-crs` refusals, and the far-away and
latitude-first refusals (16's mismatch message already names both CRSs).
`test_axis_order_is_able_to_fail` is red before 15b through its control run.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import pytest
import shapely
from pyproj import CRS, Transformer
from shapely.geometry import Polygon

import test_cli_mesh_mosaic
from geotiff_fixtures import KARTVERKET, micro_tiff, needs_codecs
from mosaic_fixtures import X0, Y0, blocks, quadrants, whole
from test_cli_mesh_dem import write_tiff
from test_cli_mesh_domain import COLS, HOLE, ROWS, SQUARE, quarter_circle
from test_cli_mesh_domain import geojson as utm33_geojson
from test_cli_mesh_mosaic import USAGE, invoke, terrain, write_tiles
from test_cli_mesh_refine import NUMBER, sentence
from test_cli_mesh_refine import field as sentence_field
from vtkread import VtkFile, read_vtk

Ring = list[tuple[float, float]]
SNAP = 1e-3  # DEFAULT_SNAP_SPACING: the noder's grid, applied in the DEM's CRS (R9)
LCC = "+proj=lcc +lat_1=60 +lat_2=65 +lat_0=62 +lon_0=15 +ellps=GRS80 +units=m +no_defs"

# Recorded from the tree before 15b; see the module docstring.
GOLDEN = {
    "square": "04b8495098bba930b85425aad36c7d4eef1f66f85585ac64c82a764c2d3b0c15",
    "quarter_circle_10m": "696056fc9a7566cb3e9e17029182f19603de66d3f663b20e074b43ea145acf58",
}


def digest(vtk: VtkFile) -> str:
    """SHA-256 over the mesh itself: points, cells and every cell and point
    array, each prefixed with its name, dtype and shape. Field data is left out,
    since 15b adds `domain_crs` and `domain_transform` to it."""
    h = hashlib.sha256()

    def put(name: str, a: Any) -> None:
        arr = np.ascontiguousarray(np.asarray(a))
        h.update(f"{name}{arr.dtype.str}{arr.shape}".encode())
        h.update(arr.tobytes())

    put("points", vtk.points)
    for i, cell in enumerate(vtk.lines):
        put(f"line{i}", cell)
    for i, cell in enumerate(vtk.polygons):
        put(f"polygon{i}", cell)
    for group in (vtk.point_scalars, vtk.scalars):
        for name in sorted(group):
            put(name, group[name].values)
    for block in sorted(vtk.cell_fields):
        for name in sorted(vtk.cell_fields[block]):
            values = vtk.cell_fields[block][name].values
            if isinstance(values, np.ndarray):
                put(f"{block}.{name}", values)
    return h.hexdigest()


def to_crs(src: str, dst: str, ring: Ring) -> Ring:
    """pyproj's own `always_xy` transform of `ring`: the oracle."""
    t = Transformer.from_crs(src, dst, always_xy=True)
    xs, ys = t.transform(np.array([p[0] for p in ring]), np.array([p[1] for p in ring]))
    return [(float(x), float(y)) for x, y in zip(xs, ys, strict=True)]


def description(src: str, dst: str = "EPSG:25833") -> str:
    return str(Transformer.from_crs(src, dst, always_xy=True).description)


def write_geojson(
    path: Path, outer: Ring, holes: tuple[Ring, ...] = (), crs: str | None = None
) -> Path:
    doc: dict[str, Any] = {"type": "Polygon", "coordinates": [[*r, r[0]] for r in (outer, *holes)]}
    if crs is not None:
        doc["crs"] = {"type": "name", "properties": {"name": crs}}
    path.write_text(json.dumps(doc))
    return path


def the_field(vtk: VtkFile, name: str) -> str:
    assert name in vtk.field_data, sorted(vtk.field_data)
    (value,) = vtk.field_data[name].values
    return str(value)


def run(tmp_path: Path, *args: str, out: str = "x.vtk") -> VtkFile:
    target = tmp_path / out
    code, output = invoke(*args, "--out", str(target))
    assert code == 0, output
    return read_vtk(target.read_bytes())


def refused(tmp_path: Path, *args: str, says: tuple[str, ...]) -> str:
    target = tmp_path / "refused.vtk"
    code, output = invoke(*args, "--out", str(target))
    assert code == USAGE, output
    assert "Traceback" not in output
    assert "No such option" not in output, "refused for the wrong reason"
    for word in says:
        assert word in output, f"{word!r} not in {output!r}"
    assert not target.exists()
    return output


def assert_ring_in_output(vtk: VtkFile, ring: Ring) -> None:
    """Every expected vertex has an output point within the noder's snap (U6)."""
    points = vtk.points[:, :2]
    for x, y in ring:
        i = int(np.argmin(np.hypot(points[:, 0] - x, points[:, 1] - y)))
        assert abs(points[i, 0] - x) <= SNAP / 2 + 1e-9, (x, y, points[i])
        assert abs(points[i, 1] - y) <= SNAP / 2 + 1e-9, (x, y, points[i])


def assert_inside(vtk: VtkFile, ring: Ring) -> None:
    """No output vertex outside the domain, in the DEM's CRS, beyond the snap."""
    grown = Polygon(ring).buffer(SNAP)
    inside = shapely.intersects_xy(grown, vtk.points[:, 0], vtk.points[:, 1])
    assert inside.all(), vtk.points[~inside][:5]


# ------------------------------------------------------------------ tiles

# `whole(9, 13)` of `terrain`, in quadrants: nodes x 500 000 .. 500 120 (dx 10),
# y 6 600 000 .. 6 599 960 (dy 5), EPSG:25833, point-registered.
ACROSS: Ring = [
    (X0 + 12.3, Y0 - 36.3),
    (X0 + 107.7, Y0 - 35.9),
    (X0 + 106.1, Y0 - 3.7),
    (X0 + 13.9, Y0 - 4.1),
]


@pytest.fixture
def quad_dir(tmp_path: Path) -> Path:
    tiles = quadrants(whole(9, 13, array=terrain(9, 13)), row_cut=4, col_cut=6, overlap=1)
    write_tiles(tmp_path / "quad", tiles)
    return tmp_path / "quad"


def mesh_args(dem: Path, domain: Path, *extra: str) -> tuple[str, ...]:
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", "1", *extra)


# ------------------------------------------------------------------ tests


class TestTheDomainInItsOwnCrs:
    """ "Tests for @tester", 15b: 4326 GeoJSON, 25832, WKT with --domain-crs
    EPSG:3035, each over a 25833 mosaic; plus a CRS with no EPSG code and
    OGC's CRS84."""

    def test_wgs84_geojson_without_a_crs_member(self, tmp_path: Path, quad_dir: Path) -> None:
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", ACROSS)
        domain = write_geojson(tmp_path / "catchment.geojson", lon_lat)
        vtk = run(tmp_path, *mesh_args(quad_dir, domain))
        expected = to_crs("EPSG:4326", "EPSG:25833", lon_lat)
        assert_ring_in_output(vtk, expected)
        assert_inside(vtk, expected)
        assert the_field(vtk, "domain_crs") == "EPSG:4326"
        assert the_field(vtk, "domain_transform") == description("EPSG:4326")
        assert the_field(vtk, "crs") == "EPSG:25833"
        assert the_field(vtk, "dem_tiles") == "ne.tif; nw.tif; se.tif; sw.tif"
        assert sentence_field(sentence(vtk), rf"achieved max error {NUMBER} m") <= 1.0

    def test_utm32_geojson(self, tmp_path: Path, quad_dir: Path) -> None:
        ring = to_crs("EPSG:25833", "EPSG:25832", ACROSS)
        domain = write_geojson(tmp_path / "d.geojson", ring, crs="urn:ogc:def:crs:EPSG::25832")
        vtk = run(tmp_path, *mesh_args(quad_dir, domain))
        assert_ring_in_output(vtk, to_crs("EPSG:25832", "EPSG:25833", ring))
        assert the_field(vtk, "domain_crs") == "EPSG:25832"
        assert the_field(vtk, "domain_transform") == description("EPSG:25832")

    def test_wkt_in_laea_europe_with_domain_crs(self, tmp_path: Path, quad_dir: Path) -> None:
        ring = to_crs("EPSG:25833", "EPSG:3035", ACROSS)
        domain = tmp_path / "d.wkt"
        domain.write_text(Polygon(ring).wkt)
        vtk = run(tmp_path, *mesh_args(quad_dir, domain, "--domain-crs", "EPSG:3035"))
        assert_ring_in_output(vtk, to_crs("EPSG:3035", "EPSG:25833", ring))
        assert the_field(vtk, "domain_crs") == "EPSG:3035"
        assert the_field(vtk, "domain_transform") == description("EPSG:3035")

    def test_a_crs_with_no_epsg_code(self, tmp_path: Path, quad_dir: Path) -> None:
        """16 refused it ("has no EPSG code"); R9 makes `crs` a string for it."""
        ring = to_crs("EPSG:25833", LCC, ACROSS)
        domain = tmp_path / "d.wkt"
        domain.write_text(Polygon(ring).wkt)
        vtk = run(tmp_path, *mesh_args(quad_dir, domain, "--domain-crs", LCC))
        assert_ring_in_output(vtk, to_crs(LCC, "EPSG:25833", ring))
        recorded = the_field(vtk, "domain_crs")
        assert recorded.isascii()
        assert CRS.from_user_input(recorded) == CRS.from_user_input(LCC)
        assert the_field(vtk, "domain_transform") == description(LCC)

    def test_ogc_crs84_member(self, tmp_path: Path, quad_dir: Path) -> None:
        lon_lat = to_crs("EPSG:25833", "OGC:CRS84", ACROSS)
        domain = write_geojson(tmp_path / "d.geojson", lon_lat, crs="urn:ogc:def:crs:OGC:1.3:CRS84")
        vtk = run(tmp_path, *mesh_args(quad_dir, domain))
        assert_ring_in_output(vtk, to_crs("OGC:CRS84", "EPSG:25833", lon_lat))
        assert CRS.from_user_input(the_field(vtk, "domain_crs")) == CRS.from_user_input("OGC:CRS84")

    def test_a_single_file_records_no_tile_list(self, tmp_path: Path, quad_dir: Path) -> None:
        inside: Ring = [
            (X0 + 5.3, Y0 - 12.9),
            (X0 + 33.7, Y0 - 12.1),
            (X0 + 32.9, Y0 - 3.1),
            (X0 + 6.1, Y0 - 3.7),
        ]
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", inside)
        domain = write_geojson(tmp_path / "d.geojson", lon_lat)
        vtk = run(tmp_path, *mesh_args(quad_dir / "nw.tif", domain))
        assert "dem_tiles" not in vtk.field_data
        assert the_field(vtk, "domain_crs") == "EPSG:4326"

    def test_a_domain_that_avoids_a_missing_tile_meshes(self, tmp_path: Path) -> None:
        """R4 point 5: NaN filler outside the needed region is never read."""
        source = whole(12, 12, dy=10.0, array=terrain(12, 12))
        write_tiles(tmp_path / "ell", blocks(source, 6, 6, skip=[(1, 1)]))
        ell: Ring = [
            (X0 + 5.5, Y0 - 5.5),
            (X0 + 105.5, Y0 - 5.5),
            (X0 + 105.5, Y0 - 35.0),
            (X0 + 35.0, Y0 - 35.0),
            (X0 + 35.0, Y0 - 105.5),
            (X0 + 5.5, Y0 - 105.5),
        ]
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", ell)
        domain = write_geojson(tmp_path / "ell.geojson", lon_lat)
        vtk = run(tmp_path, *mesh_args(tmp_path / "ell", domain))
        assert np.isfinite(vtk.points).all()
        assert_inside(vtk, to_crs("EPSG:4326", "EPSG:25833", lon_lat))


class TestTheSameCrs:
    """ "A domain already in the DEM's CRS is not transformed, and the mesh is
    bit-identical to 16's." """

    @pytest.fixture
    def bumpy(self, tmp_path: Path) -> Path:
        array = np.random.default_rng(16).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
        return write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))

    @pytest.fixture
    def no_transformer(self, monkeypatch: pytest.MonkeyPatch) -> None:
        def refuse(*args: Any, **kwargs: Any) -> Any:
            raise AssertionError(f"Transformer.from_crs{args} was called")

        monkeypatch.setattr(Transformer, "from_crs", staticmethod(refuse))

    def test_the_mesh_is_16s_bit_for_bit(
        self, tmp_path: Path, bumpy: Path, no_transformer: None
    ) -> None:
        square = utm33_geojson(tmp_path / "square.geojson", SQUARE, (HOLE,))
        vtk = run(tmp_path, *mesh_args(bumpy, square))
        assert digest(vtk) == GOLDEN["square"]

    def test_the_fields_say_so(self, tmp_path: Path, bumpy: Path) -> None:
        square = utm33_geojson(tmp_path / "square.geojson", SQUARE, (HOLE,))
        vtk = run(tmp_path, *mesh_args(bumpy, square))
        assert the_field(vtk, "domain_crs") == "EPSG:25833"
        assert "none" in the_field(vtk, "domain_transform").lower()

    @needs_codecs
    def test_the_quarter_circle_is_16s_bit_for_bit(
        self, tmp_path: Path, no_transformer: None
    ) -> None:
        domain = utm33_geojson(tmp_path / "quarter.geojson", quarter_circle())
        vtk = run(tmp_path, "--dem", str(KARTVERKET), "--domain", str(domain), "--tolerance", "10")
        assert digest(vtk) == GOLDEN["quarter_circle_10m"]


class TestUsage:
    def test_bbox_with_domain_is_refused(self, tmp_path: Path, quad_dir: Path) -> None:
        """R11: `--bbox` excludes `--domain` [15b]. The domain is in the DEM's
        own CRS, so nothing but the flag pair can refuse it."""
        domain = utm33_geojson(tmp_path / "d.geojson", ACROSS)
        refused(
            tmp_path,
            *mesh_args(quad_dir, domain),
            "--bbox", "500000", "6599960", "500120", "6600000",
            says=("--bbox", "--domain"),
        )  # fmt: skip

    @pytest.mark.parametrize("bad", ["EPSG:999999", "not-a-crs"])
    def test_an_unknown_domain_crs_is_named(self, tmp_path: Path, quad_dir: Path, bad: str) -> None:
        domain = tmp_path / "d.wkt"
        domain.write_text(Polygon(ACROSS).wkt)
        refused(tmp_path, *mesh_args(quad_dir, domain, "--domain-crs", bad), says=(bad,))


class TestRefusedAfterTheTransform:
    """R9 and "Degeneracy policy": what lands outside the DEM is refused by the
    extent check, naming the domain's CRS and the DEM's."""

    def test_a_domain_far_from_the_dem(self, tmp_path: Path, quad_dir: Path) -> None:
        far = [(x - 200_000.0, y + 50_000.0) for x, y in ACROSS]
        domain = write_geojson(tmp_path / "d.geojson", to_crs("EPSG:25833", "EPSG:4326", far))
        refused(tmp_path, *mesh_args(quad_dir, domain), says=("4326", "25833"))

    def test_a_domain_reaching_past_the_tiles(self, tmp_path: Path, quad_dir: Path) -> None:
        past = [*ACROSS[:2], (X0 + 123.0, Y0 - 3.7), ACROSS[3]]
        domain = write_geojson(tmp_path / "d.geojson", to_crs("EPSG:25833", "EPSG:4326", past))
        refused(tmp_path, *mesh_args(quad_dir, domain), says=("--domain", "outside", "4326"))

    def test_latitude_first_is_a_wrong_file(self, tmp_path: Path, quad_dir: Path) -> None:
        """ "A GeoJSON with latitude first is a wrong file, not a case, and is
        refused by the extent check." """
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", ACROSS)
        domain = write_geojson(tmp_path / "d.geojson", [(y, x) for x, y in lon_lat])
        refused(tmp_path, *mesh_args(quad_dir, domain), says=("4326", "25833"))

    def test_axis_order_is_able_to_fail(
        self, tmp_path: Path, quad_dir: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The same run as `test_wgs84_geojson_without_a_crs_member`, with every
        transformer built latitude first: the vertices land thousands of km
        away, and the extent check must refuse them."""
        real = Transformer.from_crs

        def latitude_first(*args: Any, **kwargs: Any) -> Any:
            kwargs["always_xy"] = False
            return real(*args, **kwargs)

        lon_lat = to_crs("EPSG:25833", "EPSG:4326", ACROSS)
        domain = write_geojson(tmp_path / "catchment.geojson", lon_lat)
        run(tmp_path, *mesh_args(quad_dir, domain))  # the control: unpatched, it meshes
        monkeypatch.setattr(Transformer, "from_crs", staticmethod(latitude_first))
        refused(tmp_path, *mesh_args(quad_dir, domain), says=("4326", "25833"))


class TestSeamsWithADomain:
    """Ola's Q1 revised through `--domain`: `dem_seams` is recorded as without
    one. `test_cli_mesh_mosaic.TestQ1Seams.disagreeing` (reached through its
    module, so pytest does not collect that class twice) plants +4 at global node (2, 6), which is
    (500 060, 6 599 990), inside `ACROSS`."""

    def test_a_disagreeing_pair_is_recorded(self, tmp_path: Path) -> None:
        source = whole(9, 13, array=terrain(9, 13))
        dem = test_cli_mesh_mosaic.TestQ1Seams.disagreeing(tmp_path, source)
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", ACROSS)
        vtk = run(tmp_path, *mesh_args(dem, write_geojson(tmp_path / "c.geojson", lon_lat)))
        assert the_field(vtk, "domain_crs") == "EPSG:4326"
        assert the_field(vtk, "dem_seams") == "ne.tif | nw.tif: nodes 1, max 4, median 4"

    def test_agreeing_tiles_record_none(self, tmp_path: Path, quad_dir: Path) -> None:
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", ACROSS)
        vtk = run(tmp_path, *mesh_args(quad_dir, write_geojson(tmp_path / "c.geojson", lon_lat)))
        assert the_field(vtk, "dem_seams") == "none"


class TestRealSeam:
    """The 15b Acceptance on the committed DTM10 extract: a catchment-like
    polygon in EPSG:4326 over the 6400_4 | 6400_1 seam, across its overlap."""

    UTM33_RING: ClassVar[Ring] = [
        (47_503.3, 6_467_207.7),
        (52_496.1, 6_467_301.9),
        (52_402.7, 6_469_298.3),
        (47_601.9, 6_469_203.1),
    ]

    def test_meshes_with_the_domain_fields(self, tmp_path: Path) -> None:
        seam = Path(__file__).resolve().parents[1] / "fixtures" / "dtm10" / "seam"
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", self.UTM33_RING)
        domain = write_geojson(tmp_path / "catchment.geojson", lon_lat)
        vtk = run(tmp_path, *mesh_args(seam, domain))
        assert the_field(vtk, "dem_tiles") == "6400_1_10m_z33.tif; 6400_4_10m_z33.tif"
        assert the_field(vtk, "domain_crs") == "EPSG:4326"
        assert the_field(vtk, "domain_transform") == description("EPSG:4326")
        assert_inside(vtk, to_crs("EPSG:4326", "EPSG:25833", lon_lat))
        assert np.isfinite(vtk.points).all()
