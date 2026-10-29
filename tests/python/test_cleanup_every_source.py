"""The land-cover clean-up reaches every source but ``property``'s.

`docs/increments/109-cleanup-every-source.md`, *5. Prediction and tests*,
tests 1 to 6, and *3. The design (tonight)*, items 1 to 5. Not named
invariant-critical, so no mutants.

Tests 1 to 5 run ``open_features`` on a 1,000 m square domain in EPSG:25833
at plain coordinates (0 to 1,000 m), as the design measured its figures.
Areas are compared at ``rel=1e-9`` against the input polygon's own area, at
coordinates of at most 1,000 m, where a double resolves about 1e-13 m.

The simplifier's vertex counts are not exact across platforms: they move with
the coordinates' last bits. The same fixtures translated by 11 offsets (sub-
metre ones and up to 1e6 m, the UTM offset of ``feature_fixtures.at`` among
them; shapely 2.2.0, GEOS 3.14.1, darwin/arm64) gave 50 to 58 / 26 to 28 / 14
to 18 for the design's 54 / 28 / 15, and 29 to 33 for test 4's 33. So a count
is bounded (``AT_MOST``: the design's figure plus about 20 %, at least 8
times fewer than the 514 in) and falls strictly with the band; the polygon
counts the simplifier is given (2 or 3) and every area stay exact.

Test 6 runs ``rasputin mesh`` on ``cli_helpers``' 10 m x 5 m micro DEM (cell
10 m) inside its ``SQUARE``.

PINNED HERE, where the design leaves it open (listed for ``@architect``):

- An uncoded rules file is read with its own ``[[class]]`` values (``lake``,
  ``land``); a polygon of it whose clean-up ran still has ``code`` and
  ``polygon`` None, and its lines carry the class's edge names.
- Test 1's and 2's counts are ``FeatureSet.cover_vertices`` (before, after),
  which counts line coordinates; with lakes away from the outline these are
  the design's polygon counts. Test 4's count is the simplifier's own output
  polygons (the source's two, less the rest piece), read through a spy on
  ``feature_input.simplify_borders``; ``cover_vertices`` there is lower (30
  on darwin/arm64), since the land polygon's outline edges are cut off its
  lines.
- Test 5's message is a :class:`FeatureError` whose text is exactly
  ``"<file>: border simplification refused: x; --features-band 0 meshes
  this source with its borders as read"``: ``<why>`` is the ``ValueError``'s
  own text.
- Test 6: each ``features`` entry ends "N lines, M vertices" with that
  source's own counts, sliced from the ``FeatureSet`` by ``counts``;
  ``features_joined`` is "N lines, M vertices after joining shared edges and
  adding crossings"; ``features_joined`` and ``features_accuracy_m`` each
  have a ``run_record.WORDING`` entry; with one source lacking an accuracy,
  ``features_accuracy_m`` is item 3's sentence alone.
- ``MeshingRules.cleaned`` is the field name (the design's); an uncoded file
  without ``positional_accuracy_m`` loads with ``accuracy_m`` None.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import shapely
import typer
from shapely.geometry import Point, Polygon, box
from shapely.geometry.base import BaseGeometry

import feature_fixtures as ff
import tin_engine.cli as cli
import tin_engine.run_record as rr
from cli_helpers import SQUARE, geojson, invoke, rough_dem
from feature_fixtures import UTM33, Feat, domain_of, write_geojson
from legend_fixtures import coded_map, write
from test_cli_mesh_landcover import EAST, WEST, coded, rel

SIDE = 1000.0
DOMAIN = box(0, 0, SIDE, SIDE)
RADIUS = 100.0
QUAD_SEGS = 64  # 257 coordinates per circle, closing one included
LAKE_AREA = 31_412.773  # the 256-gon's area, m2, as the design states it
REPAIR = 0.01
#: The design's counts (54 / 28 / 15 by band; 33 for test 4) plus about 20 %.
AT_MOST = {1.0: 65, 5.0: 34, 20.0: 22}
FILTER_AT_MOST = 40

#: An uncoded rules file with no positional accuracy, the shape of Sol's
#: lake file: a lake class and a land class by value.
LAKES = (
    'name = "lakes"\n'
    'attribute = "kind"\n'
    'unlisted = "drop"\n\n'
    '[[class]]\nvalue = "lake"\nedges = ["water"]\n\n'
    '[[class]]\nvalue = "land"\nedges = ["land_cover"]\n'
)
#: The same classes, coded: land 1, lake 2 (water).
LAKES_CODED = coded_map(
    "lakes-coded",
    attribute="kind",
    unlisted="drop",
    classes=(
        '[[class]]\ncode = 1\nlabel = "Land"\ncolour = "#1f8d49"\n\n'
        '[[class]]\ncode = 2\nlabel = "Lake"\ncolour = "#2532e4"\nwater = true\n'
    ),
)
UNCODED = {"lake": "lake", "land": "land"}
CODED = {"lake": 2, "land": 1}


@pytest.fixture(scope="module")
def fi() -> ModuleType:
    return ff.feature_input()


def lake(x: float, y: float = 500.0) -> Polygon:
    return Point(x, y).buffer(RADIUS, quad_segs=QUAD_SEGS)


def rules_file(tmp_path: Path, text: str, name: str = "lakes.toml") -> Path:
    return write(tmp_path / name, text)


def source(
    fi: ModuleType, tmp_path: Path, shapes: list[tuple[str, BaseGeometry]], coded_rules: bool
) -> Any:
    """One GeoJSON source of ``shapes`` (class, polygon), under the uncoded or
    the coded lake rules."""
    values = CODED if coded_rules else UNCODED
    feats = [Feat(i + 1, g, {"kind": values[kind]}) for i, (kind, g) in enumerate(shapes)]
    path = write_geojson(tmp_path / "lakes.geojson", feats)
    text = LAKES_CODED if coded_rules else LAKES
    loaded = ff.rules(str(rules_file(tmp_path, text)))
    return fi.FeatureSource(path=path, rules=loaded)


def request(fi: ModuleType, src: Any, band: float) -> Any:
    return fi.FeatureRequest(
        sources=(src,), repair_m=REPAIR, merge_same_class=True, band_m=band, outline_snap_m=0.0
    )


class Spy:
    """Every ``feature_input.simplify_borders`` call: its polygons and result."""

    def __init__(self, fi: ModuleType, monkeypatch: pytest.MonkeyPatch) -> None:
        self.calls: list[tuple[list[BaseGeometry], Any]] = []
        real = fi.simplify_borders

        def spy(polygons: Any, *args: Any, **kwargs: Any) -> Any:
            result = real(polygons, *args, **kwargs)
            self.calls.append((list(polygons), result))
            return result

        monkeypatch.setattr(fi, "simplify_borders", spy)


def vertices(geometries: Any) -> int:
    return int(sum(shapely.get_num_coordinates(g) for g in geometries))


# ---------------------------------------------------------------- tests 1 and 2


class TestLakes:
    """Tests 1 and 2: two lakes, which no polygon of their source borders,
    are simplified once the domain's uncovered rest borders them, each
    keeping its area; coded or not, the same counts."""

    @pytest.mark.parametrize("coded_rules", [False, True], ids=["uncoded", "coded"])
    def test_two_lakes_are_simplified(
        self,
        fi: ModuleType,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        coded_rules: bool,
    ) -> None:
        lakes = [lake(300.0), lake(700.0)]
        assert vertices(lakes) == 514  # the premise
        src = source(fi, tmp_path, [("lake", g) for g in lakes], coded_rules)
        counts = []
        for band, most in AT_MOST.items():
            with monkeypatch.context() as patch:
                spy = Spy(fi, patch)
                fs = fi.open_features(request(fi, src, band), domain_of(DOMAIN), UTM33)
            before, after = fs.cover_vertices
            assert before == 514 and after <= most, (band, after)
            counts.append(after)
            self.kept(fs, spy, lakes, coded_rules)
        assert counts == sorted(set(counts), reverse=True), counts  # falls strictly

    @staticmethod
    def kept(fs: Any, spy: Spy, lakes: list[Polygon], coded_rules: bool) -> None:
        """Each lake keeps its area, and its code and label only if coded."""
        (call,) = spy.calls
        rings = [g for g in call[1].polygons if g.area < DOMAIN.area / 2]  # not the rest
        assert sum(g.area for g in rings) == pytest.approx(sum(g.area for g in lakes), rel=1e-9)
        assert lakes[0].area == pytest.approx(LAKE_AREA, abs=5e-4)
        for f in fs.features:  # each lake's ring keeps its area
            for line in f.lines:
                assert Polygon(line).area == pytest.approx(lakes[0].area, rel=1e-9)
            if coded_rules:
                assert f.code == 2 and f.polygon is not None
            else:
                assert f.code is None and f.polygon is None

    def test_the_uncoded_lines_carry_the_class_names(self, fi: ModuleType, tmp_path: Path) -> None:
        src = source(fi, tmp_path, [("lake", lake(300.0))], coded_rules=False)
        fs = fi.open_features(request(fi, src, 5.0), domain_of(DOMAIN), UTM33)
        (f,) = fs.features
        assert f.mask == fi.DEFAULT_VOCABULARY.mask("water")


# ---------------------------------------------------------------- test 3


class TestTiling:
    """Test 3: two halves of the square tile it; the rest is empty, so the
    simplifier gets the source's two polygons and nothing more."""

    @pytest.mark.parametrize("coded_rules", [False, True], ids=["uncoded", "coded"])
    def test_no_rest_piece(
        self,
        fi: ModuleType,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        coded_rules: bool,
    ) -> None:
        halves = [("land", box(0, 0, 500, SIDE)), ("lake", box(500, 0, SIDE, SIDE))]
        spy = Spy(fi, monkeypatch)
        src = source(fi, tmp_path, halves, coded_rules)
        fi.open_features(request(fi, src, 5.0), domain_of(DOMAIN), UTM33)
        (call,) = spy.calls
        assert len(call[0]) == 2
        assert sum(g.area for g in call[0]) == pytest.approx(DOMAIN.area, rel=1e-12)


# ---------------------------------------------------------------- test 4


class TestTheCentimetreFilter:
    """Test 4: a rest part narrower than 1 cm is left out of the simplifier's
    input, a wider one goes in; the source simplifies alike either way."""

    @pytest.mark.parametrize(("width", "given"), [(0.008, 2), (0.012, 3)])
    def test_a_strip_below_1_cm_is_left_out(
        self,
        fi: ModuleType,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        width: float,
        given: int,
    ) -> None:
        water = lake(300.0)
        land = box(0, 0, SIDE, SIDE - width).difference(water)
        assert vertices([water, land]) == 519  # the premise
        strip = DOMAIN.difference(shapely.coverage_union_all([water, land]))
        assert strip.area == pytest.approx(SIDE * width, rel=1e-6)
        assert shapely.get_num_coordinates(strip) == 5
        spy = Spy(fi, monkeypatch)
        src = source(fi, tmp_path, [("lake", water), ("land", land)], coded_rules=False)
        fs = fi.open_features(request(fi, src, 5.0), domain_of(DOMAIN), UTM33)
        (call,) = spy.calls
        assert len(call[0]) == given
        out = call[1].polygons[:2]
        assert vertices(out) <= FILTER_AT_MOST
        assert out[0].area == pytest.approx(water.area, rel=1e-9)
        assert out[1].area == pytest.approx(land.area, rel=1e-9)
        assert water.area == pytest.approx(LAKE_AREA, abs=5e-4)
        before, after = fs.cover_vertices
        assert before == 519 and after <= vertices(out)


# ---------------------------------------------------------------- test 5


class TestRefusal:
    """Test 5: the simplifier's refusal of an uncoded source refuses the run,
    naming the file and the way out."""

    def test_the_run_is_refused(
        self, fi: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        def refused(*_: Any, **__: Any) -> Any:
            raise ValueError("x")

        monkeypatch.setattr(fi, "simplify_borders", refused)
        src = source(fi, tmp_path, [("lake", lake(300.0))], coded_rules=False)
        with pytest.raises(fi.FeatureError) as caught:
            fi.open_features(request(fi, src, 5.0), domain_of(DOMAIN), UTM33)
        assert str(caught.value) == (
            "lakes.geojson: border simplification refused: x; "
            "--features-band 0 meshes this source with its borders as read"
        )


# ---------------------------------------------------------------- item 1: property


class TestProperty:
    """Item 1: ``property``'s polygons are never cleaned; every other rules
    file's are. An uncoded file without an accuracy states none."""

    def test_the_cleaned_flag(self, tmp_path: Path) -> None:
        assert ff.rules("property").cleaned is False
        assert ff.rules("corine").cleaned is True
        assert ff.rules("corine-water").cleaned is True
        uncoded = ff.rules(str(rules_file(tmp_path, LAKES)))
        assert uncoded.cleaned is True
        assert uncoded.accuracy_m is None

    def test_property_lakes_are_not_simplified(
        self, fi: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        spy = Spy(fi, monkeypatch)
        feats = [Feat(1, lake(300.0), {"property": "water"})]
        path = write_geojson(tmp_path / "p.geojson", feats)
        src = fi.FeatureSource(path=path, rules=ff.rules("property"))
        fs = fi.open_features(request(fi, src, 5.0), domain_of(DOMAIN), UTM33)
        assert spy.calls == []
        assert fs.cover_vertices is None
        (f,) = fs.features
        assert vertices(f.lines) == 257

    def test_the_help_says_so(self) -> None:
        params = {p.name: p for p in typer.main.get_command(cli.app).commands["mesh"].params}
        help_ = params["features_rules"].help
        assert "property: each feature's own edge names, never cleaned or simplified." in help_


# ---------------------------------------------------------------- test 6

bumpy = rough_dem(16)
#: Item 3's sentence for ``lakes.toml`` on the micro DEM (cell 10 m).
NO_ACCURACY = "lakes.toml states no positional_accuracy_m; the DEM cell, 10 m, is used"
JOINED = re.compile(
    r"^\d+ lines?, \d+ vert(ex|ices) after joining shared edges and adding crossings$"
)


@pytest.fixture
def square(tmp_path: Path) -> Path:
    return geojson(tmp_path / "square.geojson", SQUARE)


@pytest.fixture
def small_lake(tmp_path: Path) -> Path:
    """An uncoded lake inside ``SQUARE``, more than 5 m from its outline."""
    (x, y) = rel(60.3, -40.1)
    feats = [Feat(1, Point(x, y).buffer(15.0, quad_segs=8), {"kind": "lake"})]
    return write_geojson(tmp_path / "lake.geojson", feats)


@pytest.fixture
def corine(tmp_path: Path) -> Path:
    return write_geojson(
        tmp_path / "corine.geojson", [coded(1, WEST, "311"), coded(2, EAST, "512")]
    )


class Opened:
    """Every ``cli.open_features`` result."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self.sets: list[Any] = []
        real = cli.open_features

        def spy(*args: Any, **kwargs: Any) -> Any:
            found = real(*args, **kwargs)
            self.sets.append(found)
            return found

        monkeypatch.setattr(cli, "open_features", spy)


def mesh(tmp_path: Path, *args: str) -> tuple[dict[str, Any], str]:
    """The record and output of one successful run."""
    rec = tmp_path / "x.json"
    code, output = invoke(
        "mesh",
        *("--tolerance", "1", "--out", str(tmp_path / "x.vtk"), "--record", str(rec)),
        *args,
    )
    assert code == 0, output
    return json.loads(rec.read_text(encoding="ascii")), output


def own_counts(found: Any, i: int) -> str:
    start = sum(found.counts[:i])
    mine = found.features[start : start + found.counts[i]]
    lines = sum(len(f.lines) for f in mine)
    coords = sum(int(shapely.get_num_coordinates(line)) for f in mine for line in f.lines)
    return f"{rr.plural(lines, 'line', 'lines')}, {rr.plural(coords, 'vertex', 'vertices')}"


class TestTheRecord:
    """Test 6: per-source counts, one joined total, and the stated default
    for a missing positional accuracy."""

    def test_each_source_counts_its_own(
        self,
        tmp_path: Path,
        bumpy: Path,
        square: Path,
        small_lake: Path,
        corine: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        rules = rules_file(tmp_path, LAKES)
        spy = Opened(monkeypatch)
        record, _ = mesh(
            tmp_path,
            *("--dem", str(bumpy), "--domain", str(square)),
            *("--features", str(small_lake), "--features-rules", str(rules)),
            *("--features", str(corine), "--features-rules", "corine"),
        )
        (found,) = spy.sets
        entries = record["features"].split("; ")
        assert len(entries) == 2
        own = [own_counts(found, i) for i in range(2)]
        assert own[0] != own[1]  # the premise: the two sources differ
        for entry, counts in zip(entries, own, strict=True):
            assert entry.endswith(f": {counts}") or entry.endswith(f", {counts}"), entry
        assert entries[0].startswith("lake.geojson, rules lakes: 1 feature, ")
        assert entries[1].startswith("corine.geojson, rules corine: 2 features, ")
        assert JOINED.match(record["features_joined"]), record["features_joined"]
        assert "features_joined" in rr.WORDING
        assert "features_accuracy_m" in rr.WORDING

    def test_a_missing_accuracy_is_the_dem_cell(
        self, tmp_path: Path, bumpy: Path, square: Path, small_lake: Path
    ) -> None:
        rules = rules_file(tmp_path, LAKES)
        record, output = mesh(
            tmp_path,
            *("--dem", str(bumpy), "--domain", str(square)),
            *("--features", str(small_lake), "--features-rules", str(rules)),
        )
        assert output.count(NO_ACCURACY) == 1
        assert record["features_accuracy_m"] == NO_ACCURACY
        assert record["features_band_m"] == 5.0  # s/2, s = max(10, 10)
        assert record["features_repair_m"] == pytest.approx(0.1, rel=1e-12)  # s/100
