"""The core takes valid polygons only: an invalid land-cover, `--prefer` or
`--reservoirs` polygon is refused by feature id with GEOS's reason.

`docs/increments/119-core-assumes-clean.md`, section 3 (the contract and the
refusal), section 6 (degeneracies) and section 7, tests 1 to 9, 11 and 12.
A polygon whose bounding box meets the read box (closed intervals) is
checked with `shapely.is_valid` after its move into the DEM's CRS; one that
is not valid stops the run with `feature_input.InvalidPolygons`, a
`FeatureError` that lists every invalid polygon, sorted by fid. Nothing is
repaired. `rasputin mesh` turns it into a usage error (exit 2) naming
`--land-cover`, `--prefer` or `--reservoirs`, and when the source is a
catalogue id whose `<id>-clean-<EPSG code>` sibling is in the catalogue,
the message names that sibling first.

THE FIXTURES, with GEOS's reasons as `@architect` measured them (shapely
2.2.0, GEOS 3.14.1) and as this suite re-checks in its premises:

- the bow-tie `Polygon([(0,0),(4,4),(4,0),(0,4)])`: `Self-intersection[2 2]`;
- the hole across its shell, shell (0,0)-(10,10), hole (5,5)-(15,6):
  `Self-intersection[10 5]`;
- AR50 Innlandet fid 191902 (NIBIO, NLOD 1.0; `tests/fixtures/ar50/NOTICE`),
  12 vertices, valid in EPSG:25832 and invalid after the move to
  EPSG:25833; its reason's figures depend on PROJ, so only the prefix
  `Self-intersection[` is matched;
- for the CLI, the bow-tie moved by (506500, 8316100) into the scene of
  `test_cli_land_cover_polygons.py` (EPSG:31983):
  `Self-intersection[506502 8316102]`.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- The refusal's head is compared with its whitespace collapsed, so where
  its lines break is free; each invalid polygon is on a line of its own
  that, stripped, reads `feature <fid>: <reason>`.
- `<name>` in the head is the file's name (`path.name`), as the readers'
  other refusals have it.
- The `(valid in the file's CRS, <own>: ...)` note shows the file's own
  CRS text as read (a GeoJSON `crs` member `EPSG:25832` gives `EPSG:25832`).
- An empty row is not one of the `m` polygons counted in the head.
- `read_polygons` returns the valid polygons far from ``xy`` too, in file
  order, and drops only the invalid far ones.
- The CLI's sibling sentence is matched as `use <sibling> in place of <id>`
  with whitespace removed (Rich wraps), before the first reason.
- The `--reservoirs` refusal comes before `cli._dem_mesh` is called.
"""

from __future__ import annotations

import inspect
import json
import os
import sqlite3
import traceback
from collections.abc import Sequence
from contextlib import closing
from dataclasses import fields
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
import shapely.affinity
from shapely.errors import GEOSException
from shapely.geometry import MultiPolygon, Polygon, box
from shapely.geometry.base import BaseGeometry

import feature_fixtures as ff
import tin_engine.cli as cli
from cli_helpers import USAGE, runner, squashed
from dataset_fixtures import write_sources
from feature_fixtures import UTM33, Feat, write_geojson
from fraction_fixtures import MESH_CRS, write_dem, write_domain
from test_cli_land_cover_polygons import cover_polygons, write_cover_geojson, write_cover_gpkg
from tin_engine.cells import overlay
from tin_engine.cli import app
from tin_engine.polygon_fractions import assign
from vector_fixtures import vector_entry

BOWTIE = Polygon([(0, 0), (4, 4), (4, 0), (0, 4)])
BOWTIE_REASON = "Self-intersection[2 2]"
HOLE_ACROSS = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)], [[(5, 5), (15, 5), (15, 6), (5, 6)]])
HOLE_REASON = "Self-intersection[10 5]"
#: A valid square clear of both, and far: the bow-tie 10 km east.
SQUARE = box(-50, -50, -40, -40)
FAR_BOWTIE = shapely.affinity.translate(BOWTIE, 10_000, 0)
#: The read box the readers are given: the far bow-tie's box misses it.
NEAR = (-100.0, -100.0, 100.0, 100.0)
#: The bow-tie with its box touching `NEAR`'s east side only (x = 100).
TOUCHING_BOWTIE = shapely.affinity.translate(BOWTIE, 100, 0)
#: Two squares overlapping each other: invalid as one MultiPolygon.
OVERLAPPING_PARTS = MultiPolygon([box(0, 0, 4, 4), box(2, 2, 6, 6)])
#: Points whose box widened by `feature_input.MARGIN` (100 m) is about `NEAR`.
XY = np.array([[0.0, 0.0], [4.0, 4.0]])
CODE = "311"

#: AR50 Innlandet fid 191902 as read from NIBIO's file, in EPSG:25832.
AR50_FID = 191902
AR50_25832 = Polygon(
    [
        (678885.524997621, 6826937.85570249),
        (678887.364248655, 6826917.93785453),
        (678898.243112952, 6826908.89825946),
        (678898.702942203, 6826903.91879882),
        (678899.622602121, 6826893.95987771),
        (678902.572129994, 6826889.210109),
        (678907.551523773, 6826889.66949791),
        (678910.730995975, 6826882.429867),
        (678912.110550418, 6826867.49121125),
        (678908.050702045, 6826857.07233377),
        (678893.801773013, 6826848.22372425),
        (678885.524997621, 6826937.85570249),
    ]
)

#: The CLI scene's bow-tie: inside the domain of `test_cli_land_cover_polygons`.
SCENE_BOWTIE = shapely.affinity.translate(BOWTIE, 506_500, 8_316_100)
SCENE_REASON = "Self-intersection[506502 8316102]"


@pytest.fixture(scope="module")
def fi() -> Any:
    return ff.feature_input()


def cover(path: Path, *rows: tuple[Any, BaseGeometry | None], crs: str = UTM33) -> Path:
    """A GeoJSON land-cover file: one feature per ``(fid, geometry)``, class 311."""
    return write_geojson(path, [Feat(fid, g, {"Code_18": CODE}) for fid, g in rows], crs=crs)


def read_near(
    path: Path,
    near: tuple[float, float, float, float] | None = NEAR,
    mesh_crs: str = UTM33,
    threads: int = 1,
) -> Any:
    """`read_class_polygons` under `corine`, a GeoJSON read whole, checked near ``near``."""
    return ff.feature_input().read_class_polygons(
        path, None, ff.classes("corine"), mesh_crs, lambda _own: None,
        threads=threads, near=near,
    )  # fmt: skip


def flat(text: str) -> str:
    return " ".join(text.split())


def feature_lines(message: str) -> list[str]:
    return [s for line in message.splitlines() if (s := line.strip()).startswith("feature ")]


def head(name: str, n: int, m: int, crs: str) -> str:
    """The refusal's words before its feature lines, whitespace collapsed."""
    return (
        f"{name}: {n} of the {m} polygons reaching the read region are not valid in {crs}, "
        f"the DEM's CRS; the mesh takes valid polygons only. Repair them in {crs} "
        "(shapely.make_valid, or QGIS's Fix geometries) and run again:"
    )


# ------------------------------------------------------------ the premises


class TestPremises:
    """The fixtures are what the design measured, and an invalid polygon
    breaks the overlays: why the core takes valid polygons only."""

    @pytest.mark.parametrize(
        ("polygon", "reason"),
        [
            (BOWTIE, BOWTIE_REASON),
            (HOLE_ACROSS, HOLE_REASON),
            (SCENE_BOWTIE, SCENE_REASON),
            (FAR_BOWTIE, "Self-intersection[10002 2]"),
            (TOUCHING_BOWTIE, "Self-intersection[102 2]"),
            (OVERLAPPING_PARTS, "Self-intersection[2 4]"),
        ],
        ids=["bowtie", "hole", "scene", "far", "touching", "overlapping-parts"],
    )
    def test_the_reasons(self, polygon: BaseGeometry, reason: str) -> None:
        assert shapely.is_valid_reason(polygon) == reason

    def test_the_square_is_valid(self) -> None:
        assert SQUARE.is_valid

    def test_the_far_and_touching_boxes(self) -> None:
        _, y0, x1, y1 = NEAR
        fx0, _, _, _ = FAR_BOWTIE.bounds
        assert fx0 > x1
        tx0, ty0, _, ty1 = TOUCHING_BOWTIE.bounds
        assert tx0 == x1 and y0 <= ty0 and ty1 <= y1

    def test_the_ar50_polygon_is_valid_as_read_and_broken_by_the_move(self) -> None:
        from tin_engine.crs import reprojector

        assert AR50_25832.is_valid
        moved = shapely.transform(AR50_25832, reprojector("EPSG:25832", "EPSG:25833"))
        assert shapely.is_valid_reason(moved).startswith("Self-intersection[")

    def test_assign_throws_on_the_bowtie(self) -> None:
        triangle = np.array([[0.0, 0.0], [4.0, 0.0], [0.0, 4.0]])
        with pytest.raises(GEOSException, match="TopologyException: side location conflict at 2 2"):
            assign(triangle, np.array([[0, 1, 2]]), [(311, BOWTIE)], known=[311])

    def test_overlay_throws_on_the_bowtie(self) -> None:
        triangle = np.array([[[0.0, 0.0], [4.0, 0.0], [0.0, 4.0]]])
        with pytest.raises(GEOSException, match="TopologyException: side location conflict at 2 2"):
            overlay(triangle, [BOWTIE])


# ------------------------------------------------------------ tests 1 to 4: read_class_polygons


class TestReadClassPolygons:
    def test_invalid_polygons_is_a_feature_error(self, fi: Any) -> None:
        assert issubclass(fi.InvalidPolygons, fi.FeatureError)

    def test_1_the_bowtie_is_refused_by_fid_and_reason(self, fi: Any, tmp_path: Path) -> None:
        path = cover(tmp_path / "cover.geojson", (3, SQUARE), (7, BOWTIE))
        with pytest.raises(fi.InvalidPolygons) as caught:
            read_near(path)
        message = str(caught.value)
        assert flat(message) == f"{head('cover.geojson', 1, 2, UTM33)} feature 7: {BOWTIE_REASON}"
        assert feature_lines(message) == [f"feature 7: {BOWTIE_REASON}"], message
        assert "feature 3" not in message, message
        assert "given CRS" not in message, message  # the CLI sends that to --land-cover-crs

    @pytest.mark.parametrize("threads", [1, 3])
    def test_2_every_invalid_polygon_is_named_in_fid_order(
        self, fi: Any, tmp_path: Path, threads: int
    ) -> None:
        """No stop at the first: the bow-tie, read first, has the larger fid."""
        path = cover(tmp_path / "cover.geojson", (5, BOWTIE), (1, SQUARE), (2, HOLE_ACROSS))
        with pytest.raises(fi.InvalidPolygons) as caught:
            read_near(path, threads=threads)
        message = str(caught.value)
        assert flat(message).startswith(head("cover.geojson", 2, 3, UTM33)), message
        assert feature_lines(message) == [
            f"feature 2: {HOLE_REASON}",
            f"feature 5: {BOWTIE_REASON}",
        ], message

    def test_3_a_far_bowtie_is_not_checked_and_built_outside_as_read(self, tmp_path: Path) -> None:
        path = cover(tmp_path / "cover.geojson", (1, SQUARE), (2, FAR_BOWTIE))
        read = read_near(path)
        assert [g.wkt for _, g in read.areas] == [SQUARE.wkt, FAR_BOWTIE.wkt]
        far = read.areas[1][1]
        built = ff.land_cover_map().build(
            [read.areas], [], ff.classes("corine"), box(-40, -40, 0, 0)
        )
        assert [code for code, _ in built.areas] == [311]
        ((code, kept),) = built.outside
        assert code == 311 and kept is far

    def test_a_box_that_only_touches_the_read_box_is_near_and_checked(
        self, fi: Any, tmp_path: Path
    ) -> None:
        """Section 6: closed intervals, as `_fixed`'s ``<=`` and ``>=``."""
        path = cover(tmp_path / "cover.geojson", (4, TOUCHING_BOWTIE))
        with pytest.raises(fi.InvalidPolygons) as caught:
            read_near(path)
        assert feature_lines(str(caught.value)) == ["feature 4: Self-intersection[102 2]"]

    def test_a_multipolygon_is_checked_whole(self, fi: Any, tmp_path: Path) -> None:
        """Section 6: parts overlapping each other are invalid together."""
        assert all(p.is_valid for p in OVERLAPPING_PARTS.geoms)  # each part alone is valid
        path = cover(tmp_path / "cover.geojson", (6, OVERLAPPING_PARTS))
        with pytest.raises(fi.InvalidPolygons) as caught:
            read_near(path)
        assert feature_lines(str(caught.value)) == ["feature 6: Self-intersection[2 4]"]

    def test_an_empty_row_is_skipped_counted_and_not_checked(self, fi: Any, tmp_path: Path) -> None:
        path = cover(tmp_path / "cover.geojson", (1, Polygon()), (2, SQUARE), (3, BOWTIE))
        with pytest.raises(fi.InvalidPolygons) as caught:
            read_near(path)
        assert flat(str(caught.value)).startswith(head("cover.geojson", 1, 2, UTM33))
        assert read_near(cover(tmp_path / "ok.geojson", (1, Polygon()), (2, SQUARE))).empty == 1

    def test_near_none_checks_every_polygon_read(self, fi: Any, tmp_path: Path) -> None:
        path = cover(tmp_path / "cover.geojson", (1, SQUARE), (2, FAR_BOWTIE))
        with pytest.raises(fi.InvalidPolygons) as caught:
            read_near(path, near=None)
        assert feature_lines(str(caught.value)) == ["feature 2: Self-intersection[10002 2]"]

    def test_a_vertex_with_no_image_is_refused_before_the_check(
        self, fi: Any, tmp_path: Path
    ) -> None:
        """Section 6: `_moved`'s refusal, as today, not this one."""
        laea = "EPSG:3035"
        bowtie = ff.moved(shapely.affinity.translate(BOWTIE, *ff.at(0, 0)), UTM33, laea)
        lost = list(ff.moved(ff.square(0, 0, 10, 10), UTM33, laea).exterior.coords)[:-1]
        lost[1] = (1e300, lost[1][1])
        path = cover(tmp_path / "cover.geojson", (1, bowtie), ("lost", Polygon(lost)), crs=laea)
        with pytest.raises(fi.FeatureError) as caught:
            read_near(path, near=None)
        assert not isinstance(caught.value, fi.InvalidPolygons)
        lost_text = "cover.geojson: feature lost has a vertex with no image"
        assert str(caught.value).startswith(lost_text), caught.value

    def test_4_a_polygon_broken_by_the_move_says_so(self, fi: Any, tmp_path: Path) -> None:
        path = cover(tmp_path / "ar50.geojson", (AR50_FID, AR50_25832), crs="EPSG:25832")
        with pytest.raises(fi.InvalidPolygons) as caught:
            read_near(path, near=None, mesh_crs="EPSG:25833")
        message = str(caught.value)
        assert flat(message).startswith(head("ar50.geojson", 1, 1, "EPSG:25833")), message
        (line,) = feature_lines(message)
        assert line.startswith(f"feature {AR50_FID}: Self-intersection["), line
        assert line.endswith("(valid in the file's CRS, EPSG:25832: broken by the move)"), line

    def test_4_the_same_polygon_in_its_own_crs_is_accepted(self, tmp_path: Path) -> None:
        path = cover(tmp_path / "ar50.geojson", (AR50_FID, AR50_25832), crs="EPSG:25832")
        read = read_near(path, near=None, mesh_crs="EPSG:25832")
        assert len(read.areas) == 1


# ------------------------------------------------------------ test 5: read_polygons


class TestReadPolygons:
    def read(self, path: Path) -> Any:
        return ff.feature_input().read_polygons(path, None, UTM33, XY)

    def test_5_tell_is_no_parameter(self, fi: Any) -> None:
        assert "tell" not in inspect.signature(fi.read_polygons).parameters

    def test_5_a_near_bowtie_is_refused_naming_file_and_fid(self, fi: Any, tmp_path: Path) -> None:
        path = cover(tmp_path / "lake.geojson", (1, SQUARE), (7, BOWTIE))
        with pytest.raises(fi.InvalidPolygons) as caught:
            self.read(path)
        message = str(caught.value)
        assert flat(message).startswith(head("lake.geojson", 1, 2, UTM33)), message
        assert feature_lines(message) == [f"feature 7: {BOWTIE_REASON}"], message

    def test_5_a_far_bowtie_is_dropped_and_the_rest_returned_as_a_list(
        self, tmp_path: Path
    ) -> None:
        far_square = shapely.affinity.translate(SQUARE, 20_000, 0)
        path = cover(tmp_path / "lake.geojson", (1, SQUARE), (2, FAR_BOWTIE), (3, far_square))
        found = self.read(path)
        assert isinstance(found, list), type(found)
        assert [g.wkt for g in found] == [SQUARE.wkt, far_square.wkt]


# ------------------------------------------------------------ test 9: the repair is gone


class TestTheRepairIsGone:
    def test_9_no_repair_function(self, fi: Any) -> None:
        assert not hasattr(fi, "_repaired")
        assert not hasattr(cli, "_say_repaired")

    def test_9_the_map_counts_no_repair_and_build_cuts(self) -> None:
        lm = ff.land_cover_map()
        assert "repaired" not in {f.name for f in fields(lm.LandCoverMap)}
        parameters = inspect.signature(lm.build).parameters
        assert "fix" not in parameters
        assert parameters["cut"].kind is inspect.Parameter.KEYWORD_ONLY
        assert parameters["cut"].default is True

    def test_9_the_record_describes_no_repair(self) -> None:
        from tin_engine.run_record import WORDING

        assert "land_cover_repaired" not in WORDING
        assert "reservoirs_repaired" not in WORDING


# ------------------------------------------------------------ tests 6 to 9 and 12: rasputin mesh


@pytest.fixture
def dem(tmp_path: Path) -> Path:
    return write_dem(tmp_path / "dem.tif")


@pytest.fixture
def domain(tmp_path: Path) -> Path:
    return write_domain(tmp_path / "domain.geojson")


def run(dem: Path, domain: Path, target: Path, *extra: str) -> tuple[int, str]:
    """`rasputin mesh` on the scene; no exception escapes, none is printed."""
    given = ["--dem", str(dem), "--domain", str(domain), "--tolerance", "3", "--out", str(target)]
    result = runner.invoke(app, ["mesh", *given, *extra])
    if result.exit_code not in (0, USAGE):
        raise AssertionError(
            result.output + "".join(traceback.format_exception(result.exception))
            if result.exception
            else result.output
        )
    assert "Traceback" not in result.output, result.output
    return result.exit_code, result.output


def bad_cover(path: Path) -> Path:
    return write_cover_geojson(path, [*cover_polygons(), (311, SCENE_BOWTIE)])


def good_cover(path: Path) -> Path:
    return write_cover_geojson(path, cover_polygons())


def corine(path: Path) -> tuple[str, ...]:
    return ("--land-cover", str(path), "--land-cover-classes", "corine")


def refused(output: str, option: str, fid: int) -> str:
    """The squashed output, after checking it names ``option``, ``fid`` and the reason."""
    text = squashed(output)
    assert f"Invalidvaluefor{option}" in text, output
    assert squashed(f"feature {fid}: {SCENE_REASON}") in text, output
    return text


class TestTheRun:
    def test_6_land_cover(self, tmp_path: Path, dem: Path, domain: Path) -> None:
        target = tmp_path / "x.vtk"
        code, output = run(dem, domain, target, *corine(bad_cover(tmp_path / "bad.geojson")))
        assert code == USAGE, output
        refused(output, "--land-cover", 4)
        assert not target.exists()

    def test_7_prefer(self, tmp_path: Path, dem: Path, domain: Path) -> None:
        lake = write_geojson(tmp_path / "lake.geojson", [Feat(2, SCENE_BOWTIE, {})], crs=MESH_CRS)
        target = tmp_path / "x.vtk"
        code, output = run(
            dem, domain, target, *corine(good_cover(tmp_path / "cover.geojson")),
            "--prefer", f"lake={lake}",
        )  # fmt: skip
        assert code == USAGE, output
        refused(output, "--prefer", 2)
        assert not target.exists()

    def test_8_reservoirs_are_refused_before_the_mesh(
        self, tmp_path: Path, dem: Path, domain: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        reached: list[bool] = []

        def no_mesh(*_args: Any, **_kwargs: Any) -> Any:
            reached.append(True)
            raise AssertionError("the mesh was built before the reservoirs were checked")

        monkeypatch.setattr(cli, "_dem_mesh", no_mesh)
        res = write_geojson(
            tmp_path / "res.geojson", [Feat(3, SCENE_BOWTIE, {"magasinNr": "1"})], crs=MESH_CRS
        )
        target = tmp_path / "x.gpkg"
        code, output = run(
            dem, domain, target, *corine(good_cover(tmp_path / "cover.geojson")),
            "--reservoirs", str(res),
        )  # fmt: skip
        assert reached == [], output
        assert code == USAGE, output
        refused(output, "--reservoirs", 3)
        assert not target.exists()

    def test_9_a_land_cover_record_has_no_repair_count(
        self, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        record = tmp_path / "x.json"
        code, output = run(
            dem, domain, tmp_path / "x.vtk", *corine(good_cover(tmp_path / "cover.geojson")),
            "--record", str(record),
        )  # fmt: skip
        assert code == 0, output
        assert "land_cover_repaired" not in json.loads(record.read_text(encoding="ascii"))

    def test_9_a_reservoirs_record_has_no_repair_count(
        self, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        res = write_geojson(
            tmp_path / "res.geojson",
            [Feat(1, box(*SCENE_BOWTIE.bounds), {"magasinNr": "1"})],
            crs=MESH_CRS,
        )
        target = tmp_path / "x.gpkg"
        code, output = run(
            dem, domain, target, *corine(good_cover(tmp_path / "cover.geojson")),
            "--reservoirs", str(res),
        )  # fmt: skip
        assert code == 0, output
        with closing(sqlite3.connect(target)) as conn:
            names = {n for (n,) in conn.execute("SELECT name FROM run")}
        assert "reservoirs" in names, names
        assert "reservoirs_repaired" not in names, names


# ------------------------------------------------------------ tests 11 and 12: the clean sibling


def entry(key: str) -> Any:
    from tin_engine.sources import RemoteSource

    return RemoteSource(
        id=key, kind="local", path=Path(f"/data/{key}.gpkg"), values="features",
        credit="Test", licence_note="Test",
    )  # fmt: skip


class TestCleanSibling:
    def sibling(self, key: str | None, crs: str, known: Sequence[str]) -> Any:
        from tin_engine.sources import clean_sibling

        return clean_sibling(key, crs, {k: entry(k) for k in known})

    def test_11_the_sibling_for_the_dems_epsg_code(self) -> None:
        assert self.sibling("bad", "EPSG:31983", ["bad", "bad-clean-31983"]) == "bad-clean-31983"

    def test_11_no_key(self) -> None:
        assert self.sibling(None, "EPSG:31983", ["bad", "bad-clean-31983"]) is None

    def test_11_a_sibling_for_another_crs_only(self) -> None:
        assert self.sibling("bad", "EPSG:31983", ["bad", "bad-clean-25833"]) is None

    def test_11_a_crs_without_an_epsg_code(self) -> None:
        crs = "+proj=tmerc +lon_0=7.3 +ellps=GRS80"
        assert self.sibling("bad", crs, ["bad", "bad-clean-31983"]) is None

    def test_11_an_id_ending_in_clean_is_not_stripped(self) -> None:
        """Section 3: `ar50-innlandet-clean-25832` on an EPSG:25833 DEM finds none."""
        known = ["x-clean-25832", "x-clean-25833"]
        assert self.sibling("x-clean-25832", "EPSG:25833", known) is None


class TestTheSiblingNamed:
    """Test 12: a user catalogue holding `bad`, a local GeoPackage with the
    scene's bow-tie under `classes = "corine"`, and `bad-clean-31983`."""

    @pytest.fixture
    def catalogue(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Any:
        monkeypatch.delenv("RASPUTIN_DATA", raising=False)
        bad = write_cover_gpkg(tmp_path / "bad.gpkg", [*cover_polygons(), (311, SCENE_BOWTIE)])
        good = write_cover_gpkg(tmp_path / "good.gpkg", cover_polygons())
        config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"

        def write(with_sibling: bool) -> Path:
            entries = {"bad": vector_entry(bad, "Bad cover", classes="corine")}
            if with_sibling:
                entries["bad-clean-31983"] = vector_entry(good, "Clean cover", classes="corine")
            write_sources(config, entries)
            return bad

        return write

    def test_12_land_cover_names_the_sibling_first(
        self, catalogue: Any, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        catalogue(True)
        target = tmp_path / "x.vtk"
        code, output = run(dem, domain, target, "--land-cover", "bad")
        assert code == USAGE, output
        text = refused(output, "--land-cover", 4)
        use = text.find("usebad-clean-31983inplaceofbad")
        assert 0 <= use < text.find(squashed(SCENE_REASON)), output
        assert not target.exists()

    def test_12_no_sibling_in_the_catalogue(
        self, catalogue: Any, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        catalogue(False)
        code, output = run(dem, domain, tmp_path / "x.vtk", "--land-cover", "bad")
        assert code == USAGE, output
        refused(output, "--land-cover", 4)
        assert "-clean-" not in squashed(output), output

    def test_12_the_file_by_path_names_no_sibling(
        self, catalogue: Any, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        bad = catalogue(True)
        code, output = run(dem, domain, tmp_path / "x.vtk", *corine(bad))
        assert code == USAGE, output
        refused(output, "--land-cover", 4)
        assert "-clean-" not in squashed(output), output

    def test_12_reservoirs_name_the_sibling(
        self, catalogue: Any, tmp_path: Path, dem: Path, domain: Path
    ) -> None:
        catalogue(True)
        target = tmp_path / "x.gpkg"
        code, output = run(
            dem, domain, target, *corine(good_cover(tmp_path / "cover.geojson")),
            "--reservoirs", "bad",
        )  # fmt: skip
        assert code == USAGE, output
        text = refused(output, "--reservoirs", 4)
        assert "usebad-clean-31983inplaceofbad" in text, output
        assert not target.exists()
