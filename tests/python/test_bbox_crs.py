"""Every `--bbox` names its CRS: `CrsBox`, `--bbox-crs`, and the box moved into the frame.

`docs/increments/114-bbox-crs.md`, section 4, tests 1 to 5. Tests 6 and 8 (a
fetch dry run, and hints that run as printed) are in `test_cli_fetch.py`,
beside its range server; test 7 (the MCP schema and instructions) is in
`test_mcp_schema.py` and `test_mcp_datasets.py`.

The figures are `docs/increments/114-probes/box_probe.py`'s, which applies the
design's move (`EDGE_STEPS = 2000` steps a side), with one exception: the
share of the box over the coverage's east edge. The probe measured it against
the coverage's degree bounds as a rectangle, 0.5183; the listing measures it
against the union of the tiles' own lon/lat outlines (`source_report._lonlat`,
8 segments a side), and there the box's outline gives 0.4964, its bounds
0.4961. That figure was measured here, independently of the code under test
(`EDGE_SHARE` below says how).

Fixtures: the quad (`test_cli_mesh_mosaic.mosaic_dir`, 13 x 9 nodes, x
500000..500120 every 10 m, y 6599960..6600000 every 5 m, EPSG:25833) and
`dataset_fixtures`' bowl (EPSG:25833).

HOW THIS FILE GOES RED: `tin_engine.domain` has no `CrsBox` and no CLI command
has `--bbox-crs`. The unit tests reach `CrsBox` through `crs_fixtures.crs_box`
or the `domain` fixture, so each fails on its own; the CLI tests fail on
"No such option" or on their assertions.
"""

from __future__ import annotations

import importlib
import math
import os
import socket
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from pyproj import CRS, Transformer

from cli_helpers import USAGE, invoke, squashed
from crs_fixtures import crs_box, refuse_point_moves
from dataset_fixtures import LOCAL, by_id, local_entry, write_local, write_sources
from mosaic_fixtures import quadrants, whole
from test_cli_mesh_mosaic import run_vtk, same_mesh, terrain, write_tiles
from test_cli_sources import listing
from tin_engine.crs import transform_bounds

UTM33 = "EPSG:25833"
DEGREES = "EPSG:4326"

#: The quad's box, off its node lines: cols 2..9 and rows 1..7 once snapped
#: outward (x 500020..500090, y 6599965..6599995), 8 x 7 = 56 nodes.
QUAD_BOX = (500_022.0, 6_599_969.0, 500_087.0, 6_599_991.0)
#: `QUAD_BOX` moved into EPSG:4326 (`box_probe.txt`), to the 9 decimals printed.
QUAD_DEGREES = (15.000389013, 59.538070744, 15.001538378, 59.53826831)

#: A UTM zone 60 box whose image in degrees spans 179.54 E to 177.84 W.
ACROSS_180 = ("EPSG:32660", (700_000.0, 5_000_000.0, 900_000.0, 5_100_000.0))
#: A polar stereographic box round the north pole.
ROUND_THE_POLE = ("EPSG:3413", (-500_000.0, -500_000.0, 500_000.0, 500_000.0))

#: The bowl's boxes in EPSG:25833 (`box_probe.txt`): wholly inside, wholly
#: outside, and over the coverage's east edge.
INSIDE_UTM = (508_522.7, 6_576_808.3, 511_380.3, 6_582_383.5)
OUTSIDE_UTM = (555_607.8, 6_651_832.7, 561_353.6, 6_663_057.1)
EDGE_UTM = (517_060.4, 6_573_496.4, 522_780.6, 6_579_094.1)
#: The share of `EDGE_UTM` the bowl covers, by area: the union of the four
#: tiles' `_lonlat` outlines intersected with `EDGE_UTM`'s outline in degrees
#: (2000 equal steps a side through `crs.reprojector`), over the outline's
#: area; 0.496418 unrounded. The same share of the outline's bounds is
#: 0.496067, which rounds to 0.4961: 3e-4 away, so `abs=1e-4` tells the two
#: apart.
EDGE_SHARE = 0.4964

NO_BOX_CRS = (
    "--bbox needs --bbox-crs, the CRS its four numbers are in (EPSG:4326 for longitude and"
    " latitude in degrees); a box's CRS is never assumed"
)
ORDER = "XMIN must not exceed XMAX, nor YMIN YMAX"


@pytest.fixture(scope="module")
def domain() -> ModuleType:
    return importlib.import_module("tin_engine.domain")


def no_network(*args: Any, **kwargs: Any) -> Any:
    raise AssertionError("the command opened a socket")


def closed_sockets(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(socket, "socket", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)


def dense_bounds(src: str, dst: str, box: tuple[float, ...], n: int = 200_001) -> np.ndarray:
    """The bounds of `box`'s edges in `dst`, `n` points an edge: the reference
    the edge probe used (`114-probes/edge_probe.py`, `true_bounds`)."""
    move = Transformer.from_crs(src, dst, always_xy=True)
    x0, y0, x1, y1 = box
    s = np.linspace(0.0, 1.0, n)
    xs = np.concatenate([x0 + (x1 - x0) * s, np.full(n, x1), x1 - (x1 - x0) * s, np.full(n, x0)])
    ys = np.concatenate([np.full(n, y0), y0 + (y1 - y0) * s, np.full(n, y1), y1 - (y1 - y0) * s])
    u, v = move.transform(xs, ys)
    return np.array([np.min(u), np.min(v), np.max(u), np.max(v)])


# ------------------------------------------------------------ 1. CrsBox


class TestCrsBox:
    def test_edge_steps_is_the_designs(self, domain: ModuleType) -> None:
        assert domain.EDGE_STEPS == 2000

    @pytest.mark.parametrize(
        ("crs", "dst"),
        [(UTM33, UTM33), (DEGREES, "OGC:CRS84"), (UTM33, CRS.from_epsg(25833))],
        ids=["same-text", "4326-as-crs84", "a-crs-object"],
    )
    def test_the_same_crs_returns_the_numbers_bit_for_bit_and_moves_no_point(
        self, monkeypatch: pytest.MonkeyPatch, crs: str, dst: Any
    ) -> None:
        numbers = QUAD_BOX if crs == UTM33 else QUAD_DEGREES
        box = crs_box(crs, numbers)
        refuse_point_moves(monkeypatch)
        assert box.to_crs(dst) == numbers

    def test_utm_to_degrees(self) -> None:
        # Degrees near 59 N; the largest value checked is 59.54, to 1e-9 deg (0.1 mm).
        got = crs_box(UTM33, QUAD_BOX).to_crs(DEGREES)
        np.testing.assert_allclose(got, QUAD_DEGREES, rtol=0, atol=1e-9)

    def test_a_point_box_gives_a_point(self) -> None:
        """Equal corners are allowed: `sources` takes a point."""
        got = crs_box(UTM33, (500_000.0, 6_600_000.0, 500_000.0, 6_600_000.0)).to_crs(DEGREES)
        np.testing.assert_allclose(got, (15.0, 59.53834913, 15.0, 59.53834913), rtol=0, atol=1e-8)
        assert got[0] == got[2] and got[1] == got[3]

    def test_the_outline_samples_each_side_at_edge_steps_plus_one_points(
        self, domain: ModuleType
    ) -> None:
        box = crs_box(UTM33, QUAD_BOX)
        xy = np.asarray(box.outline(UTM33))
        assert xy.shape == (4 * (domain.EDGE_STEPS + 1), 2)
        x0, y0, x1, y1 = QUAD_BOX
        corners = {(x0, y0), (x1, y0), (x1, y1), (x0, y1)}
        assert corners <= {(float(x), float(y)) for x, y in xy}
        assert xy[:, 0].min() == x0 and xy[:, 0].max() == x1
        assert xy[:, 1].min() == y0 and xy[:, 1].max() == y1

    def test_the_outline_in_another_crs_is_its_points_moved(self) -> None:
        xy = np.asarray(crs_box(UTM33, QUAD_BOX).outline(DEGREES))
        here = np.asarray(crs_box(UTM33, QUAD_BOX).outline(UTM33))
        lon, lat = Transformer.from_crs(UTM33, DEGREES, always_xy=True).transform(
            here[:, 0], here[:, 1]
        )
        np.testing.assert_allclose(xy, np.column_stack([lon, lat]), rtol=1e-15, atol=0)

    @pytest.mark.parametrize(
        ("src", "dst"), [("EPSG:32723", DEGREES), (DEGREES, "EPSG:32723")], ids=["utm", "degrees"]
    )
    def test_a_sao_francisco_box_falls_inside_the_true_image_by_at_most_2_mm(
        self, src: str, dst: str
    ) -> None:
        """A 12 x 14 degree box (about 1300 x 1550 km, the basin's) against
        200 001 points an edge. Measured: 0.0 m one way, 1.6 mm the other
        (`edge_probe.txt`). Bound 2 mm, in metres (degrees times 111 320 m);
        no edge lies outside the reference by more than 1e-6 m, rounding."""
        sf = (-48.0, -21.0, -36.0, -7.0)
        box = sf if src == DEGREES else transform_bounds(DEGREES, src, sf)
        got = np.array(crs_box(src, box).to_crs(dst))
        ref = dense_bounds(src, dst, box)
        inside = np.array([got[0] - ref[0], got[1] - ref[1], ref[2] - got[2], ref[3] - got[3]])
        metres = inside * (111_320.0 if dst == DEGREES else 1.0)
        assert metres.max() <= 2e-3, metres
        assert metres.min() >= -1e-6, metres

    def test_a_pole_box_into_a_projected_crs_is_not_refused(self) -> None:
        """The ±180° check applies only into a geographic CRS: round the pole
        into Arctic polar stereographic (EPSG:3995) the image is one piece."""
        crs, box = ROUND_THE_POLE
        got = crs_box(crs, box).to_crs("EPSG:3995")
        assert all(map(math.isfinite, got)) and got[0] < got[2] and got[1] < got[3]

    def test_latitude_90_is_allowed(self) -> None:
        assert crs_box(DEGREES, (10.0, 80.0, 11.0, 90.0)).box == (10.0, 80.0, 11.0, 90.0)


class TestCrsBoxRefusals:
    @pytest.mark.parametrize(
        "box",
        [
            (math.nan, 0.0, 1.0, 1.0),
            (0.0, -math.inf, 1.0, 1.0),
            (0.0, 0.0, math.inf, 1.0),
            (0.0, 0.0, 1.0, math.nan),
        ],
        ids=["nan-xmin", "inf-ymin", "inf-xmax", "nan-ymax"],
    )
    def test_a_non_finite_value(self, box: tuple[float, ...]) -> None:
        with pytest.raises(ValueError, match="finite"):
            crs_box(UTM33, box)

    @pytest.mark.parametrize(
        "box", [(1.0, 0.0, 0.0, 1.0), (0.0, 1.0, 1.0, 0.0)], ids=["x-inverted", "y-inverted"]
    )
    def test_an_inverted_box(self, box: tuple[float, ...]) -> None:
        with pytest.raises(ValueError, match=ORDER):
            crs_box(UTM33, box)

    def test_an_unreadable_crs(self) -> None:
        with pytest.raises(ValueError, match="cannot read the CRS 'EPSG:0'"):
            crs_box("EPSG:0", QUAD_BOX)

    @pytest.mark.parametrize(
        "box", [(10.0, 59.0, 11.0, 90.5), (10.0, -90.5, 11.0, 59.0)], ids=["north", "south"]
    )
    def test_a_latitude_past_90_in_a_geographic_crs(self, box: tuple[float, ...]) -> None:
        with pytest.raises(ValueError, match="latitude must lie within ±90° in EPSG:4326"):
            crs_box(DEGREES, box)

    @pytest.mark.parametrize(
        ("crs", "box"), [ACROSS_180, ROUND_THE_POLE], ids=["across-180", "round-the-pole"]
    )
    def test_an_image_across_180_or_round_a_pole(self, crs: str, box: tuple[float, ...]) -> None:
        """The outline's neighbours step 360° of longitude (`box_probe.txt`);
        pyproj's own bounds would come back west > east, or as the whole circle."""
        words = rf"^the box \(.+\) in {crs} crosses ±180° or holds a pole in EPSG:4326; split it$"
        with pytest.raises(ValueError, match=words):
            crs_box(crs, box).to_crs(DEGREES)

    def test_no_image(self) -> None:
        """UTM 33 at 1e8 m has no image in LAEA Europe: pyproj gives inf."""
        box = crs_box("EPSG:32633", (500_000.0, 6_600_000.0, 1e8, 1e8))
        words = r"^the box \(.+\) in EPSG:32633 has no image in EPSG:3035$"
        with pytest.raises(ValueError, match=words):
            box.to_crs("EPSG:3035")


# ------------------------------------------------------------ DemRequest


class TestDemRequest:
    """`box` joins `bounds` and `domain`; at most one of the three."""

    def test_a_box_alone_is_accepted(self, tmp_path: Path) -> None:
        from tin_engine.dem_input import DemRequest

        box = crs_box(DEGREES, QUAD_DEGREES)
        assert DemRequest(sources=(tmp_path / "a.tif",), box=box).box == box

    def test_bounds_and_box_are_refused(self, tmp_path: Path) -> None:
        from tin_engine.dem_input import DemRequest
        from tin_engine.io.models import Bounds

        bounds = Bounds(x_min=QUAD_BOX[0], y_min=QUAD_BOX[1], x_max=QUAD_BOX[2], y_max=QUAD_BOX[3])
        with pytest.raises(ValueError):
            DemRequest(sources=(tmp_path / "a.tif",), bounds=bounds, box=crs_box(UTM33, QUAD_BOX))

    def test_a_box_and_a_domain_are_refused(self, tmp_path: Path) -> None:
        import shapely

        from tin_engine.dem_input import DemRequest
        from tin_engine.domain import DomainPolygon

        ring = shapely.box(*QUAD_BOX)
        with pytest.raises(ValueError):
            DemRequest(
                sources=(tmp_path / "a.tif",),
                box=crs_box(UTM33, QUAD_BOX),
                domain=DomainPolygon(polygon=ring, crs=UTM33),
            )


# ------------------------------------------------------------ 2, 3. mesh


@pytest.fixture
def quad(tmp_path: Path) -> Path:
    return write_quad(tmp_path / "tiles")


def write_quad(directory: Path) -> Path:
    """The quad, as `test_cli_mesh_mosaic.mosaic_dir` writes it, in `directory`."""
    write_tiles(directory, quadrants(whole(9, 13, array=terrain(9, 13)), 4, 6, 1))
    return directory


def quad_nodes() -> set[tuple[float, float]]:
    """Cols 2..9, rows 1..7 of the quad: every node of `QUAD_BOX` snapped outward."""
    return {(500_000.0 + 10 * c, 6_600_000.0 - 5 * r) for c in range(2, 10) for r in range(1, 8)}


class TestMesh:
    def test_a_box_in_the_dems_crs_meshes_its_snapped_window(
        self, tmp_path: Path, quad: Path
    ) -> None:
        """Test 2: the numbers are used as given, so the mesh is the one
        without `--bbox-crs` before it existed: all 56 nodes, at stride 1."""
        vtk = run_vtk(
            tmp_path / "utm.vtk", "--dem", str(quad), "--stride", "1",
            "--bbox", *map(repr, QUAD_BOX), "--bbox-crs", UTM33,
        )  # fmt: skip
        assert len(vtk.points) == 56
        assert {(float(x), float(y)) for x, y, _ in vtk.points} == quad_nodes()

    def test_a_box_in_degrees_meshes_the_same_bit_for_bit(self, tmp_path: Path, quad: Path) -> None:
        """Test 3. The margin: the round trip moves the box's edges by at most
        0.9 mm (`box_probe.txt`), and they lie 2 to 9 m from a node line."""
        utm = run_vtk(
            tmp_path / "utm.vtk", "--dem", str(quad), "--stride", "1",
            "--bbox", *map(repr, QUAD_BOX), "--bbox-crs", UTM33,
        )  # fmt: skip
        degrees = run_vtk(
            tmp_path / "deg.vtk", "--dem", str(quad), "--stride", "1",
            "--bbox", *map(repr, QUAD_DEGREES), "--bbox-crs", DEGREES,
        )  # fmt: skip
        same_mesh(degrees, utm)

    def test_with_out_crs_a_box_in_degrees_meshes_as_the_box_in_the_target_crs(
        self, tmp_path: Path, quad: Path
    ) -> None:
        """The reprojected path: the box is moved into `--out-crs` and becomes
        the target grid's domain. EPSG:3006 maps the quad's numbers onto
        themselves (SWEREF 99 TM and UTM 33 on ETRS89 share their projection),
        so `QUAD_BOX` in EPSG:3006 is the same box. Measured before this
        change, with the box then read in `--out-crs`: 162 points on the
        target grid's 5 m lines, x 500010..500095, y 6599960..6600000. The
        box's edges lie at least 1 m from a 5 m line; the move, under 1 mm."""
        target = run_vtk(
            tmp_path / "t.vtk", "--dem", str(quad), "--stride", "1", "--out-crs", "EPSG:3006",
            "--bbox", *map(repr, QUAD_BOX), "--bbox-crs", "EPSG:3006",
        )  # fmt: skip
        xs, ys = target.points[:, 0], target.points[:, 1]
        assert len(target.points) == 162
        assert (xs.min(), xs.max(), ys.min(), ys.max()) == (
            500010.0,
            500095.0,
            6599960.0,
            6600000.0,
        )
        degrees = run_vtk(
            tmp_path / "d.vtk", "--dem", str(quad), "--stride", "1", "--out-crs", "EPSG:3006",
            "--bbox", *map(repr, QUAD_DEGREES), "--bbox-crs", DEGREES,
        )  # fmt: skip
        same_mesh(degrees, target)


# ------------------------------------------------------------ 4. refusals


@pytest.fixture
def config(monkeypatch: pytest.MonkeyPatch) -> Path:
    """The user's sources file, as `test_cli_sources.config` finds it."""
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)
    return Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"


def command(name: str, root: Path) -> tuple[str, ...]:
    """`name` with what it needs besides a box; the quad for `mesh`, a shipped
    source for `fetch` (sockets are closed, so nothing is requested)."""
    if name == "mesh":
        return ("mesh", "--dem", str(write_quad(root / "tiles")), "--out", str(root / "m.vtk"))
    if name == "fetch":
        return ("fetch", "glo30", "--cache", str(root / "cache"))
    return ("sources",)


def nothing_done(root: Path) -> None:
    assert not (root / "m.vtk").exists()
    cache = root / "cache"
    assert not cache.exists() or not any(cache.iterdir()), sorted(cache.rglob("*"))


COMMANDS = ["mesh", "fetch", "sources"]


class TestRefusals:
    @pytest.mark.parametrize("name", COMMANDS)
    def test_bbox_without_bbox_crs(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str
    ) -> None:
        args = command(name, tmp_path)
        closed_sockets(monkeypatch)
        code, output = invoke(*args, "--bbox", *map(repr, QUAD_BOX))
        assert code == USAGE, output
        assert squashed(f"Invalid value for --bbox-crs: {NO_BOX_CRS}") in squashed(output), output
        nothing_done(tmp_path)

    @pytest.mark.parametrize("name", COMMANDS)
    def test_bbox_crs_without_bbox(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str
    ) -> None:
        args = command(name, tmp_path)
        closed_sockets(monkeypatch)
        code, output = invoke(*args, "--bbox-crs", UTM33)
        assert code == USAGE, output
        wanted = "Invalid value for --bbox-crs: applies only with --bbox"
        assert squashed(wanted) in squashed(output), output
        nothing_done(tmp_path)

    @pytest.mark.parametrize("name", COMMANDS)
    def test_an_unreadable_bbox_crs(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str
    ) -> None:
        args = command(name, tmp_path)
        closed_sockets(monkeypatch)
        code, output = invoke(*args, "--bbox", *map(repr, QUAD_BOX), "--bbox-crs", "EPSG:0")
        assert code == USAGE, output
        wanted = "Invalid value for --bbox-crs: cannot read the CRS 'EPSG:0'"
        assert squashed(wanted) in squashed(output), output
        nothing_done(tmp_path)

    @pytest.mark.parametrize("name", COMMANDS)
    @pytest.mark.parametrize(
        ("box", "crs", "words"),
        [
            ((1.0, 0.0, 0.0, 1.0), UTM33, ORDER),
            ((0.0, 0.0, 1.0, math.nan), UTM33, "finite"),
            ((10.0, 59.0, 11.0, 90.5), DEGREES, "latitude must lie within ±90° in EPSG:4326"),
        ],
        ids=["inverted", "nan", "latitude"],
    )
    def test_a_broken_invariant_names_bbox(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        name: str,
        box: tuple[float, ...],
        crs: str,
        words: str,
    ) -> None:
        args = command(name, tmp_path)
        closed_sockets(monkeypatch)
        code, output = invoke(*args, "--bbox", *map(repr, box), "--bbox-crs", crs)
        assert code == USAGE, output
        assert "Invalid value for --bbox:" in output, output
        assert squashed(words) in squashed(output), output
        assert "validation error" not in output, output
        nothing_done(tmp_path)

    def test_sources_refuses_a_box_whose_image_crosses_180(self) -> None:
        crs, box = ACROSS_180
        code, output = invoke("sources", "--bbox", *map(repr, box), "--bbox-crs", crs)
        assert code == USAGE, output
        wanted = "crosses ±180° or holds a pole in EPSG:4326; split it"
        assert squashed("Invalid value for --bbox-crs: the box") in squashed(output), output
        assert squashed(wanted) in squashed(output), output


# ------------------------------------------------------------ 5. sources


class TestSourcesShare:
    """Test 5: the share of a box in EPSG:25833 the bowl covers, by area."""

    @pytest.fixture
    def registered(self, config: Path, tmp_path: Path) -> Path:
        return write_sources(config, {LOCAL: local_entry(write_local(tmp_path / "tiles"))})

    @pytest.mark.parametrize(
        ("box", "share"),
        [(INSIDE_UTM, 1.0), (OUTSIDE_UTM, 0.0), (EDGE_UTM, EDGE_SHARE)],
        ids=["inside", "outside", "over-the-east-edge"],
    )
    def test_the_share(self, registered: Path, box: tuple[float, ...], share: float) -> None:
        found = by_id(listing("--bbox", *map(repr, box), "--bbox-crs", UTM33))
        # The listing rounds to 1e-4.
        assert found[LOCAL]["box_covered"] == pytest.approx(share, abs=1e-4)

    def test_a_point_in_utm_is_covered_or_not(self, registered: Path) -> None:
        x, y = INSIDE_UTM[:2]
        found = by_id(listing("--bbox", *map(repr, (x, y, x, y)), "--bbox-crs", UTM33))
        assert found[LOCAL]["box_covered"] == 1.0

    def test_a_box_in_degrees_is_measured_as_before(self, registered: Path) -> None:
        """`PARTLY` in EPSG:4326, `dataset_fixtures.SHARE_PARTLY`: unchanged."""
        from dataset_fixtures import PARTLY, SHARE_PARTLY

        found = by_id(listing("--bbox", *map(repr, PARTLY), "--bbox-crs", DEGREES))
        assert found[LOCAL]["box_covered"] == pytest.approx(SHARE_PARTLY, abs=1e-4)
