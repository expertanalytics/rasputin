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
  it, and no transformer is made. Relational, on the same machine: the run
  equals, in every point, cell and array, the same run with `cli.open_dem`
  replaced by 16's data flow (the whole file, and the domain object exactly as
  `read_domain` returned it, never through `DomainPolygon.to_crs`), and the
  domain `_dem_mesh` receives has the read vertices bit for bit. Since
  15f-3 (`docs/increments/15f-edge-strip.md`, S4) the edge strip's vertices
  agree to rounding only, because their lattice coordinates are measured from
  the window's corner and the two flows cut different windows; everything a
  window cannot change stays bit for bit (`assert_16s_but_the_strip`). This
  replaced two digests recorded on macOS arm64 at `d34d79d` (test amendment
  after PR #106's CI): Linux x86 with GCC gives the square a different digest,
  so a recorded digest pins a platform, not a behaviour. The platform-stable
  absolute anchor for a same-CRS domain through the CLI is increment 18's
  `test_refine_golden.py::test_the_cli_with_start_min_angle_0_matches_the_digest`
  (the quarter circle, whose GeoJSON names EPSG:25833); it is referenced, not
  duplicated.
- `--bbox` with `--domain` is a usage error naming both.
- Since increment 25 (`docs/increments/25-plain-output.md`, D2 and D4)
  `domain_crs`, `domain_transform`, `dem_tiles` and `dem_seams` are `--stats`
  rows, read here through `inputs`, and absent from the mesh file.
- `dem_seams` is written on the domain path as without a domain: the
  disagreeing pair with a domain in EPSG:4326, `none: the tiles agree where
  they overlap` when the overlaps agree
  (test amendment after the 15b review). With a domain it counts only the
  nodes inside the needed region (Ola, 2026-09-28): `none` for a disagreement
  wholly outside it, however much of the plan's rectangle it fills.
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
same-CRS runs (then recorded as digests, now relational), the two unknown
`--domain-crs` refusals, and the far-away and latitude-first refusals (16's
mismatch message already names both CRSs).
`test_axis_order_is_able_to_fail` is red before 15b through its control run.
"""

from __future__ import annotations

import dataclasses
import io
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import pytest
import shapely
from pyproj import CRS, Transformer
from shapely.geometry import Polygon

import test_cli_mesh_mosaic
import test_dem_input_domain
import tin_engine.cli as cli
from cli_helpers import HOLE, SQUARE, Ring, geojson, mesh_to_vtk, refused, rough_dem
from geotiff_fixtures import KARTVERKET, needs_codecs
from mosaic_fixtures import X0, Y0, blocks, quadrants, whole
from recordread import sizes_row
from test_cli_mesh_domain import quarter_circle
from test_cli_mesh_mosaic import terrain, write_tiles
from test_cli_mesh_refine import SEAMS_AGREE, file_field, stats_row
from tin_engine.dem_input import DemInput, DemRequest, open_dem
from tin_engine.domain import DomainPolygon, read_domain
from tin_engine.io.geotiff import decode_dem
from vtkread import VtkFile

SNAP = 1e-3  # DEFAULT_SNAP_SPACING: the noder's grid, applied in the DEM's CRS (R9)
LCC = "+proj=lcc +lat_1=60 +lat_2=65 +lat_0=62 +lon_0=15 +ellps=GRS80 +units=m +no_defs"


def to_crs(src: str, dst: str, ring: Ring) -> Ring:
    """pyproj's own `always_xy` transform of `ring`: the oracle."""
    t = Transformer.from_crs(src, dst, always_xy=True)
    xs, ys = t.transform(np.array([p[0] for p in ring]), np.array([p[1] for p in ring]))
    return [(float(x), float(y)) for x, y in zip(xs, ys, strict=True)]


def description(src: str, dst: str = "EPSG:25833") -> str:
    return str(Transformer.from_crs(src, dst, always_xy=True).description)


def the_field(vtk: VtkFile, name: str) -> str:
    return file_field(vtk, name)


#: Fields increment 25 moved from the mesh file to ``--stats`` (D2, D4).
MOVED = ("domain", "domain_crs", "domain_transform", "dem_tiles", "dem_seams")


def inputs(tmp_path: Path, name: str, out: str = "x.vtk") -> str:
    """The ``--stats`` row ``name`` of the last ``mesh_to_vtk`` writing ``out``."""
    return stats_row((tmp_path / out).with_suffix(".md").read_text(encoding="utf-8"), name)


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


#: `TestTheSameCrs`'s DEM, bound here: a class-level binding would get the instance.
bumpy = rough_dem(16)


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
        domain = geojson(tmp_path / "catchment.geojson", lon_lat, crs=None)
        vtk = mesh_to_vtk(tmp_path, *mesh_args(quad_dir, domain))
        expected = to_crs("EPSG:4326", "EPSG:25833", lon_lat)
        assert_ring_in_output(vtk, expected)
        assert_inside(vtk, expected)
        assert inputs(tmp_path, "domain_crs") == "EPSG:4326"
        assert inputs(tmp_path, "domain_transform") == description("EPSG:4326")
        assert the_field(vtk, "crs") == "EPSG:25833"
        assert inputs(tmp_path, "dem_tiles") == "ne.tif; nw.tif; se.tif; sw.tif"
        assert float(the_field(vtk, "max_error_m")) <= 1.0
        for name in MOVED:
            assert name not in vtk.field_data, f"{name} moved to --stats (increment 25, D2)"

    def test_utm32_geojson(self, tmp_path: Path, quad_dir: Path) -> None:
        ring = to_crs("EPSG:25833", "EPSG:25832", ACROSS)
        domain = geojson(tmp_path / "d.geojson", ring, crs="urn:ogc:def:crs:EPSG::25832")
        vtk = mesh_to_vtk(tmp_path, *mesh_args(quad_dir, domain))
        assert_ring_in_output(vtk, to_crs("EPSG:25832", "EPSG:25833", ring))
        assert inputs(tmp_path, "domain_crs") == "EPSG:25832"
        assert inputs(tmp_path, "domain_transform") == description("EPSG:25832")

    def test_wkt_in_laea_europe_with_domain_crs(self, tmp_path: Path, quad_dir: Path) -> None:
        ring = to_crs("EPSG:25833", "EPSG:3035", ACROSS)
        domain = tmp_path / "d.wkt"
        domain.write_text(Polygon(ring).wkt)
        vtk = mesh_to_vtk(tmp_path, *mesh_args(quad_dir, domain, "--domain-crs", "EPSG:3035"))
        assert_ring_in_output(vtk, to_crs("EPSG:3035", "EPSG:25833", ring))
        assert inputs(tmp_path, "domain_crs") == "EPSG:3035"
        assert inputs(tmp_path, "domain_transform") == description("EPSG:3035")

    def test_a_crs_with_no_epsg_code(self, tmp_path: Path, quad_dir: Path) -> None:
        """16 refused it ("has no EPSG code"); R9 makes `crs` a string for it."""
        ring = to_crs("EPSG:25833", LCC, ACROSS)
        domain = tmp_path / "d.wkt"
        domain.write_text(Polygon(ring).wkt)
        vtk = mesh_to_vtk(tmp_path, *mesh_args(quad_dir, domain, "--domain-crs", LCC))
        assert_ring_in_output(vtk, to_crs(LCC, "EPSG:25833", ring))
        recorded = inputs(tmp_path, "domain_crs")
        assert recorded.isascii()
        assert CRS.from_user_input(recorded) == CRS.from_user_input(LCC)
        assert inputs(tmp_path, "domain_transform") == description(LCC)

    def test_ogc_crs84_member(self, tmp_path: Path, quad_dir: Path) -> None:
        lon_lat = to_crs("EPSG:25833", "OGC:CRS84", ACROSS)
        domain = geojson(tmp_path / "d.geojson", lon_lat, crs="urn:ogc:def:crs:OGC:1.3:CRS84")
        vtk = mesh_to_vtk(tmp_path, *mesh_args(quad_dir, domain))
        assert_ring_in_output(vtk, to_crs("OGC:CRS84", "EPSG:25833", lon_lat))
        recorded = CRS.from_user_input(inputs(tmp_path, "domain_crs"))
        assert recorded == CRS.from_user_input("OGC:CRS84")

    def test_a_single_file_records_no_tile_list(self, tmp_path: Path, quad_dir: Path) -> None:
        inside: Ring = [
            (X0 + 5.3, Y0 - 12.9),
            (X0 + 33.7, Y0 - 12.1),
            (X0 + 32.9, Y0 - 3.1),
            (X0 + 6.1, Y0 - 3.7),
        ]
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", inside)
        domain = geojson(tmp_path / "d.geojson", lon_lat, crs=None)
        vtk = mesh_to_vtk(tmp_path, *mesh_args(quad_dir / "nw.tif", domain))
        assert the_field(vtk, "dem_source") == "nw.tif"
        assert inputs(tmp_path, "domain_crs") == "EPSG:4326"

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
        domain = geojson(tmp_path / "ell.geojson", lon_lat, crs=None)
        vtk = mesh_to_vtk(tmp_path, *mesh_args(tmp_path / "ell", domain))
        assert np.isfinite(vtk.points).all()
        assert_inside(vtk, to_crs("EPSG:4326", "EPSG:25833", lon_lat))


class TestTheSameCrs:
    """ "A domain already in the DEM's CRS is not transformed, and the mesh is
    bit-identical to 16's." Relational, on one machine (module docstring): the
    run as it is against the run with 16's data flow put back."""

    @pytest.fixture
    def no_transformer(self, monkeypatch: pytest.MonkeyPatch) -> None:
        def refuse(*args: Any, **kwargs: Any) -> Any:
            raise AssertionError(f"Transformer.from_crs{args} was called")

        monkeypatch.setattr(Transformer, "from_crs", staticmethod(refuse))

    @pytest.fixture
    def handed_on(self, monkeypatch: pytest.MonkeyPatch) -> list[DomainPolygon | None]:
        """Every domain the CLI hands `_dem_mesh`, in call order."""
        seen: list[DomainPolygon | None] = []
        real = cli._dem_mesh

        def spy(*args: Any, **kwargs: Any) -> Any:
            seen.append(kwargs.get("domain", args[7] if len(args) > 7 else None))
            return real(*args, **kwargs)

        monkeypatch.setattr(cli, "_dem_mesh", spy)
        return seen

    @staticmethod
    def as_16_opened_it(request: DemRequest) -> DemInput:
        """16's data flow for one file and a domain: the whole file, opened as
        without a domain, and the domain object as `read_domain` returned it.
        Neither `DomainPolygon.to_crs` nor the window cut to the domain runs."""
        opened = open_dem(request.model_copy(update={"domain": None}))
        return dataclasses.replace(opened, domain=request.domain)

    def assert_as_16(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        handed_on: list[DomainPolygon | None],
        domain: Path,
        *args: str,
    ) -> None:
        now = mesh_to_vtk(tmp_path, *args, out="now.vtk")
        with monkeypatch.context() as m:
            m.setattr(cli, "open_dem", self.as_16_opened_it)
            before = mesh_to_vtk(tmp_path, *args, out="before.vtk")
        assert len(handed_on) == 2, handed_on
        read = read_domain(domain)
        for given in handed_on:
            assert given is not None
            assert given.crs == read.crs
            for got, want in zip(rings(given), rings(read), strict=True):
                assert got.dtype == want.dtype and got.tobytes() == want.tobytes()
        assert len(now.points) > 0
        dem = Path(args[args.index("--dem") + 1])
        report = (tmp_path / "now.md").read_text(encoding="utf-8")
        start = int(sizes_row(report, "start vertices"))
        assert_16s_but_the_strip(now, before, dem, start)

    def test_the_mesh_is_16s_up_to_the_strips_rounding(
        self,
        tmp_path: Path,
        bumpy: Path,
        no_transformer: None,
        handed_on: list[DomainPolygon | None],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        square = geojson(tmp_path / "square.geojson", SQUARE, (HOLE,))
        self.assert_as_16(tmp_path, monkeypatch, handed_on, square, *mesh_args(bumpy, square))

    def test_the_fields_say_so(self, tmp_path: Path, bumpy: Path) -> None:
        square = geojson(tmp_path / "square.geojson", SQUARE, (HOLE,))
        mesh_to_vtk(tmp_path, *mesh_args(bumpy, square))
        assert inputs(tmp_path, "domain_crs") == "EPSG:25833"
        assert "none" in inputs(tmp_path, "domain_transform").lower()

    @needs_codecs
    def test_the_quarter_circle_is_16s_up_to_the_strips_rounding(
        self,
        tmp_path: Path,
        no_transformer: None,
        handed_on: list[DomainPolygon | None],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        domain = geojson(tmp_path / "quarter.geojson", quarter_circle())
        args = ("--dem", str(KARTVERKET), "--domain", str(domain), "--tolerance", "10")
        self.assert_as_16(tmp_path, monkeypatch, handed_on, domain, *args)


def assert_16s_but_the_strip(now: VtkFile, before: VtkFile, dem: Path, start: int) -> None:
    """15f's S4: 16's mesh bit for bit, except the edge strip's vertices.

    A strip point is a fraction computed in lattice coordinates measured
    from the corner of the window the core is given (15f, L16); 16's data
    flow opens the whole file, today's cuts the window to the domain, so
    such a vertex can differ in its last bits. What a window cannot change
    stays bit-identical: the triangles, the constraint edges and every cell
    array (their masks), the start vertices (the first `start`; trim keeps
    order), and every vertex that is a DEM node, with its z. Taking the first
    `start` vertices as the start's rests on refine writing its start
    vertices first and trim keeping order; if either ever changed, a strip
    vertex could fall in that prefix and be held to bit-identity, so the
    test would fail rather than pass wrongly. Every other
    vertex agrees within 1e-9 lattice units in each coordinate and within
    `1e-9 max(1, |z|)` in z. A difference in connectivity means a predicate
    flipped on a rounding-level difference, and fails here."""
    assert [c.tolist() for c in now.polygons] == [c.tolist() for c in before.polygons]
    assert [c.tolist() for c in now.lines] == [c.tolist() for c in before.lines]
    assert sorted(now.cell_fields) == sorted(before.cell_fields)
    for block in now.cell_fields:
        for name, array in now.cell_fields[block].items():
            other = before.cell_fields[block][name].values
            assert np.array_equal(np.asarray(array.values), np.asarray(other)), (block, name)
    for name in now.scalars:
        assert np.array_equal(now.scalars[name].values, before.scalars[name].values), name
    a, b = now.points, before.points
    assert a.shape == b.shape
    m = decode_dem(io.BytesIO(dem.read_bytes())).meta
    col, row = (a[:, 0] - m.x_min) / m.delta_x, (m.y_max - a[:, 1]) / m.delta_y
    node = (col == np.round(col)) & (row == np.round(row))
    node &= m.x_min + np.round(col) * m.delta_x == a[:, 0]
    node &= m.y_max - np.round(row) * m.delta_y == a[:, 1]
    exact = node.copy()
    exact[:start] = True
    assert np.array_equal(a[exact], b[exact]), "a start or DEM-node vertex moved"
    for name in now.point_scalars:
        ours, theirs = now.point_scalars[name].values, before.point_scalars[name].values
        assert np.array_equal(np.asarray(ours)[exact], np.asarray(theirs)[exact]), name
    rest = ~exact
    d_col = np.abs(a[rest, 0] - b[rest, 0]) / m.delta_x
    d_row = np.abs(a[rest, 1] - b[rest, 1]) / m.delta_y
    assert (d_col <= 1e-9).all() and (d_row <= 1e-9).all(), (d_col.max(), d_row.max())
    bound = 1e-9 * np.maximum(1.0, np.abs(a[rest, 2]))
    assert (np.abs(a[rest, 2] - b[rest, 2]) <= bound).all()


def rings(domain: DomainPolygon) -> list[np.ndarray]:
    p = domain.polygon
    return [np.asarray(r.coords) for r in (p.exterior, *p.interiors)]


class TestUsage:
    def test_bbox_with_domain_is_refused(self, tmp_path: Path, quad_dir: Path) -> None:
        """R11: `--bbox` excludes `--domain` [15b]. The domain is in the DEM's
        own CRS, so nothing but the flag pair can refuse it."""
        domain = geojson(tmp_path / "d.geojson", ACROSS)
        refused(
            tmp_path,
            "mesh",
            *mesh_args(quad_dir, domain),
            "--bbox", "500000", "6599960", "500120", "6600000",
            says=("--bbox", "--domain"),
        )  # fmt: skip

    @pytest.mark.parametrize("bad", ["EPSG:999999", "not-a-crs"])
    def test_an_unknown_domain_crs_is_named(self, tmp_path: Path, quad_dir: Path, bad: str) -> None:
        domain = tmp_path / "d.wkt"
        domain.write_text(Polygon(ACROSS).wkt)
        refused(tmp_path, "mesh", *mesh_args(quad_dir, domain, "--domain-crs", bad), says=(bad,))


class TestRefusedAfterTheTransform:
    """R9 and "Degeneracy policy": what lands outside the DEM is refused by the
    extent check, naming the domain's CRS and the DEM's."""

    def test_a_domain_far_from_the_dem(self, tmp_path: Path, quad_dir: Path) -> None:
        far = [(x - 200_000.0, y + 50_000.0) for x, y in ACROSS]
        domain = geojson(tmp_path / "d.geojson", to_crs("EPSG:25833", "EPSG:4326", far), crs=None)
        refused(tmp_path, "mesh", *mesh_args(quad_dir, domain), says=("4326", "25833"))

    def test_a_domain_reaching_past_the_tiles(self, tmp_path: Path, quad_dir: Path) -> None:
        past = [*ACROSS[:2], (X0 + 123.0, Y0 - 3.7), ACROSS[3]]
        domain = geojson(tmp_path / "d.geojson", to_crs("EPSG:25833", "EPSG:4326", past), crs=None)
        refused(
            tmp_path, "mesh", *mesh_args(quad_dir, domain), says=("--domain", "outside", "4326")
        )

    def test_latitude_first_is_a_wrong_file(self, tmp_path: Path, quad_dir: Path) -> None:
        """ "A GeoJSON with latitude first is a wrong file, not a case, and is
        refused by the extent check." """
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", ACROSS)
        domain = geojson(tmp_path / "d.geojson", [(y, x) for x, y in lon_lat], crs=None)
        refused(tmp_path, "mesh", *mesh_args(quad_dir, domain), says=("4326", "25833"))

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
        domain = geojson(tmp_path / "catchment.geojson", lon_lat, crs=None)
        mesh_to_vtk(tmp_path, *mesh_args(quad_dir, domain))  # the control: unpatched, it meshes
        monkeypatch.setattr(Transformer, "from_crs", staticmethod(latitude_first))
        refused(tmp_path, "mesh", *mesh_args(quad_dir, domain), says=("4326", "25833"))


class TestSeamsWithADomain:
    """Ola's Q1 revised through `--domain`: `dem_seams` is recorded as without
    one. `test_cli_mesh_mosaic.TestQ1Seams.disagreeing` (reached through its
    module, so pytest does not collect that class twice) plants +4 at global node (2, 6), which is
    (500 060, 6 599 990), inside `ACROSS`."""

    def test_a_disagreeing_pair_is_recorded(self, tmp_path: Path) -> None:
        source = whole(9, 13, array=terrain(9, 13))
        dem = test_cli_mesh_mosaic.TestQ1Seams.disagreeing(tmp_path, source)
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", ACROSS)
        mesh_to_vtk(tmp_path, *mesh_args(dem, geojson(tmp_path / "c.geojson", lon_lat, crs=None)))
        assert inputs(tmp_path, "domain_crs") == "EPSG:4326"
        assert (
            inputs(tmp_path, "dem_seams")
            == "ne.tif and nw.tif disagree at 1 node, by up to 4 m (median 4 m)"
        )

    def test_agreeing_tiles_record_none(self, tmp_path: Path, quad_dir: Path) -> None:
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", ACROSS)
        mesh_to_vtk(
            tmp_path, *mesh_args(quad_dir, geojson(tmp_path / "c.geojson", lon_lat, crs=None))
        )
        assert inputs(tmp_path, "dem_seams") == SEAMS_AGREE

    @pytest.mark.parametrize(
        ("shape", "recorded"),
        [
            ("STRIP", SEAMS_AGREE),
            ("ELL", "ne.tif and nw.tif disagree at 2 nodes, by up to 2 m (median 1.25 m)"),
        ],
    )
    def test_only_the_needed_region_is_counted(
        self, tmp_path: Path, shape: str, recorded: str
    ) -> None:
        """Ola, 2026-09-28. `test_dem_input_domain.TestSeamsInsideTheNeededRegion`
        (reached through its module, so pytest does not collect it twice): the
        strip's region misses every planted node, the L's takes two of four."""
        planted = test_dem_input_domain.TestSeamsInsideTheNeededRegion
        dem = planted.disagreeing(tmp_path)
        lon_lat = to_crs("EPSG:25833", "EPSG:4326", planted.utm33(getattr(planted, shape)))
        mesh_to_vtk(tmp_path, *mesh_args(dem, geojson(tmp_path / "c.geojson", lon_lat, crs=None)))
        assert inputs(tmp_path, "dem_seams") == recorded


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
        domain = geojson(tmp_path / "catchment.geojson", lon_lat, crs=None)
        vtk = mesh_to_vtk(tmp_path, *mesh_args(seam, domain))
        assert inputs(tmp_path, "dem_tiles") == "6400_1_10m_z33.tif; 6400_4_10m_z33.tif"
        assert inputs(tmp_path, "domain_crs") == "EPSG:4326"
        assert inputs(tmp_path, "domain_transform") == description("EPSG:4326")
        assert_inside(vtk, to_crs("EPSG:4326", "EPSG:25833", lon_lat))
        assert np.isfinite(vtk.points).all()
