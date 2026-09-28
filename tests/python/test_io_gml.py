"""`tin_engine.io.gml`: the legacy CORINE GML, read with the standard library (16b-1, Q6).

`docs/increments/16b-terrain-polygons.md`, "Ruled by Ola", Q6 (b): the 30 MB
`tests/fixtures/corine/0000_4326_corine2018_4e6064_GML.gml` stays and 16b
reads it "as a second real fixture: a small GML2 reader (standard-library XML,
no GDAL, about 40 lines) and a class map for `clc18_kode`".

**What the file was, and its repair.** OGR wrote it as GML2
(`gml:coordinates`, `x,y` pairs, EPSG:4326 longitude first). As committed in
2020 (`8e30af4`) it was not well-formed XML: `</gml:featureMember>` was
missing after feature `sql_statement.40965` (line 1097) and the document
ended without closing its root. Ola ruled on 2026-09-28: "Repear the broken
file." `tests/fixtures/corine/repair_gml.py` adds the two closing tags, and the
committed file is its output. Every one of its 200 features
(`ogr:sql_statement`) is a `Polygon` (counted below from the file, not typed
here).

Pinned by this suite (see "Pinned by the red suite (16b-1/2)" and
"Consequences for the red suite"):

- `read_gml(stream, attribute) -> GmlDocument`, `stream` a binary stream (the
  reader opens nothing). `GmlDocument` has `crs` (the geometries' `srsName`,
  or `None` when none gives one) and `features`.
- A feature is an element with a `fid` attribute directly inside a
  `gml:featureMember`; its `fid` is that attribute, its `geometry` the shapely
  geometry of its one geometry property (`Polygon` with `outerBoundaryIs` and
  `innerBoundaryIs`, `MultiPolygon`, `LineString`, `MultiLineString`,
  `Point`; a third coordinate is dropped), and its `value` the text of its
  child element named `attribute` (`None` when absent). Features in document
  order.
- **The reader is strict** (R1): a document that is not well-formed XML is
  refused with `GmlError` (a `ValueError`), including one that ends with
  elements still open, even after complete features.
- Geometries naming different `srsName`s are refused.

HOW THIS FILE GOES RED: `tin_engine.io.gml` is imported inside a fixture.
"""

from __future__ import annotations

import ast
import difflib
import importlib
import importlib.util
import io
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import shapely
from pyproj import Transformer
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon, box

from gpkg_fixtures import LEGACY_GML

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src_python" / "tin_engine"
REPAIR_SCRIPT = LEGACY_GML.parent / "repair_gml.py"
#: Where the legacy GML lived when it was committed in 2020 (`8e30af4`).
GML_2020_PATH = ".circleci/rasputin_data/corine/0000_4326_corine2018_4e6064_GML.gml"
#: The benchmark tile's node rectangle, EPSG:25833.
TILE_7908_3 = box(799_750.0, 7_899_750.0, 850_250.0, 7_950_250.0)

HEAD = (
    '<?xml version="1.0" encoding="utf-8" ?>\n'
    "<ogr:FeatureCollection\n"
    '     xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n'
    '     xmlns:ogr="http://ogr.maptools.org/"\n'
    '     xmlns:gml="http://www.opengis.net/gml">\n'
    "  <gml:boundedBy><gml:Box><gml:coord><gml:X>0</gml:X><gml:Y>0</gml:Y></gml:coord>"
    "<gml:coord><gml:X>1</gml:X><gml:Y>1</gml:Y></gml:coord></gml:Box></gml:boundedBy>\n"
)
TAIL = "</ogr:FeatureCollection>\n"


@pytest.fixture(scope="module")
def gml() -> ModuleType:
    return importlib.import_module("tin_engine.io.gml")


def coords(points: Any) -> str:
    return " ".join(",".join(repr(float(v)) for v in p) for p in points)


def ring(tag: str, points: Any) -> str:
    return (
        f"<gml:{tag}><gml:LinearRing><gml:coordinates>{coords(points)}"
        f"</gml:coordinates></gml:LinearRing></gml:{tag}>"
    )


def polygon(shape: Polygon, srs: str | None = "EPSG:4326") -> str:
    attr = f' srsName="{srs}"' if srs else ""
    inner = "".join(ring("innerBoundaryIs", r.coords) for r in shape.interiors)
    return (
        f"<gml:Polygon{attr}>{ring('outerBoundaryIs', shape.exterior.coords)}{inner}</gml:Polygon>"
    )


def geometry_xml(shape: Any, srs: str | None = "EPSG:4326") -> str:
    attr = f' srsName="{srs}"' if srs else ""
    if isinstance(shape, Polygon):
        return polygon(shape, srs)
    if isinstance(shape, MultiPolygon):
        members = "".join(
            f"<gml:polygonMember>{polygon(p, None)}</gml:polygonMember>" for p in shape.geoms
        )
        return f"<gml:MultiPolygon{attr}>{members}</gml:MultiPolygon>"
    if isinstance(shape, LineString):
        body = f"<gml:coordinates>{coords(shape.coords)}</gml:coordinates>"
        return f"<gml:LineString{attr}>{body}</gml:LineString>"
    if isinstance(shape, MultiLineString):
        members = "".join(
            f"<gml:lineStringMember><gml:LineString><gml:coordinates>{coords(g.coords)}"
            "</gml:coordinates></gml:LineString></gml:lineStringMember>"
            for g in shape.geoms
        )
        return f"<gml:MultiLineString{attr}>{members}</gml:MultiLineString>"
    assert isinstance(shape, Point)
    return f"<gml:Point{attr}><gml:coordinates>{coords(shape.coords)}</gml:coordinates></gml:Point>"


def member(fid: str, shape: Any, code: str | None = "311", srs: str | None = "EPSG:4326") -> str:
    value = f"      <ogr:clc18_kode>{code}</ogr:clc18_kode>\n" if code is not None else ""
    return (
        "  <gml:featureMember>\n"
        f'    <ogr:sql_statement fid="{fid}">\n'
        f"      <ogr:geometryProperty>{geometry_xml(shape, srs)}</ogr:geometryProperty>\n"
        f"{value}"
        "      <ogr:sl_sdeid>99</ogr:sl_sdeid>\n"
        "    </ogr:sql_statement>\n"
        "  </gml:featureMember>\n"
    )


def document(*members: str, tail: str = TAIL) -> io.BytesIO:
    return io.BytesIO((HEAD + "".join(members) + tail).encode("utf-8"))


LAND = Polygon(
    [(21.0, 70.0), (21.5, 70.0), (21.5, 70.2), (21.0, 70.2)],
    [[(21.1, 70.05), (21.2, 70.05), (21.2, 70.1), (21.1, 70.1)]],
)
LAKE = Polygon([(22.0, 70.0), (22.1, 70.0), (22.1, 70.1)])


class TestReadGml:
    def test_polygons_with_holes_their_fids_and_values(self, gml: ModuleType) -> None:
        doc = gml.read_gml(
            document(member("f.1", LAND, "311"), member("f.2", LAKE, "512")), "clc18_kode"
        )
        assert doc.crs == "EPSG:4326"
        assert [f.fid for f in doc.features] == ["f.1", "f.2"]
        assert [f.value for f in doc.features] == ["311", "512"]
        assert shapely.equals_exact(doc.features[0].geometry, LAND, tolerance=0.0)
        assert len(doc.features[0].geometry.interiors) == 1

    def test_coordinates_are_x_then_y_bit_for_bit(self, gml: ModuleType) -> None:
        """GML2's `gml:coordinates` as OGR writes it: `x,y` pairs, longitude
        first even for EPSG:4326. The reader keeps that order and the doubles."""
        shape = Polygon([(21.123456789012345, 70.98765432109876), (21.5, 70.0), (21.0, 70.5)])
        (feature,) = gml.read_gml(document(member("a", shape)), "clc18_kode").features
        assert list(feature.geometry.exterior.coords) == list(shape.exterior.coords)

    @pytest.mark.parametrize(
        "shape",
        [
            pytest.param(MultiPolygon([LAKE, box(23, 70, 23.1, 70.1)]), id="multipolygon"),
            pytest.param(LineString([(21, 70), (21.1, 70.1), (21.3, 70.0)]), id="linestring"),
            pytest.param(
                MultiLineString([[(21, 70), (21.1, 70.1)], [(22, 70), (22.2, 70.1)]]),
                id="multilinestring",
            ),
            pytest.param(Point(21.0, 70.0), id="point"),
        ],
    )
    def test_every_accepted_geometry_type(self, gml: ModuleType, shape: Any) -> None:
        """Points are read, not refused here: refusing input that is not
        geometry is `feature_input`'s, which names the feature and the file."""
        (feature,) = gml.read_gml(document(member("a", shape)), "clc18_kode").features
        assert shapely.equals_exact(feature.geometry, shape, tolerance=0.0)

    def test_a_third_coordinate_is_dropped(self, gml: ModuleType) -> None:
        text = (
            '  <gml:featureMember><ogr:sql_statement fid="z"><ogr:geometryProperty>'
            '<gml:LineString srsName="EPSG:4326"><gml:coordinates>21,70,5 21.1,70.1,6'
            "</gml:coordinates></gml:LineString></ogr:geometryProperty>"
            "<ogr:clc18_kode>311</ogr:clc18_kode></ogr:sql_statement></gml:featureMember>\n"
        )
        (feature,) = gml.read_gml(document(text), "clc18_kode").features
        assert not feature.geometry.has_z
        assert list(feature.geometry.coords) == [(21.0, 70.0), (21.1, 70.1)]

    def test_a_missing_attribute_is_none(self, gml: ModuleType) -> None:
        (feature,) = gml.read_gml(document(member("a", LAKE, None)), "clc18_kode").features
        assert feature.value is None

    def test_another_attribute_is_read_by_name(self, gml: ModuleType) -> None:
        (feature,) = gml.read_gml(document(member("a", LAKE)), "sl_sdeid").features
        assert feature.value == "99"

    def test_no_srs_name_gives_no_crs(self, gml: ModuleType) -> None:
        doc = gml.read_gml(document(member("a", LAKE, srs=None)), "clc18_kode")
        assert doc.crs is None

    def test_two_srs_names_are_refused(self, gml: ModuleType) -> None:
        with pytest.raises(gml.GmlError):
            gml.read_gml(
                document(member("a", LAKE), member("b", LAKE, srs="EPSG:4258")), "clc18_kode"
            )

    def test_an_empty_collection_reads_as_nothing(self, gml: ModuleType) -> None:
        doc = gml.read_gml(document(), "clc18_kode")
        assert doc.features == () or list(doc.features) == []

    def test_a_document_has_no_completeness_flag(self, gml: ModuleType) -> None:
        """A strict reader reads a whole document or refuses it, so the first
        design's `complete` is gone (Ola, 2026-09-28)."""
        doc = gml.read_gml(document(member("a", LAKE)), "clc18_kode")
        assert not hasattr(doc, "complete")


class TestMalformed:
    """The reader is strict: what is not well-formed XML is refused."""

    def test_a_document_whose_root_is_never_closed_is_refused(self, gml: ModuleType) -> None:
        """The 2020 file's second defect, after complete features."""
        with pytest.raises(gml.GmlError):
            gml.read_gml(document(member("a", LAKE), member("b", LAND), tail=""), "clc18_kode")

    def test_a_member_left_open_is_refused(self, gml: ModuleType) -> None:
        """The 2020 file's first defect: one `</gml:featureMember>` missing,
        with the root closed. The members after it nest inside it, so the
        root's closing tag does not match and the parser stops there."""
        broken = member("a", LAKE).replace("  </gml:featureMember>\n", "")
        with pytest.raises(gml.GmlError):
            gml.read_gml(document(broken, member("b", LAND)), "clc18_kode")

    def test_both_defects_together_are_refused(self, gml: ModuleType) -> None:
        """The 2020 file exactly, in small."""
        broken = member("a", LAKE).replace("  </gml:featureMember>\n", "")
        with pytest.raises(gml.GmlError):
            gml.read_gml(document(broken, member("b", LAND), tail=""), "clc18_kode")

    def test_the_refusal_names_the_parsers_line(self, gml: ModuleType) -> None:
        """R1: the refusal names the line. Here the head is 6 lines, the
        broken member 6 and the next 7, so the root's closing tag, where the
        parser stops, is on line 20."""
        broken = member("a", LAKE).replace("  </gml:featureMember>\n", "")
        with pytest.raises(gml.GmlError, match=r"\bline 20\b"):
            gml.read_gml(document(broken, member("b", LAND)), "clc18_kode")

    def test_a_document_ending_inside_a_feature_is_refused(self, gml: ModuleType) -> None:
        text = HEAD + member("a", LAKE) + member("b", LAND)[:200]
        with pytest.raises(gml.GmlError):
            gml.read_gml(io.BytesIO(text.encode()), "clc18_kode")

    def test_a_syntax_error_is_refused(self, gml: ModuleType) -> None:
        text = (
            HEAD + member("a", LAKE).replace("</ogr:sql_statement>", "</ogr:sql_statemen>") + TAIL
        )
        with pytest.raises(gml.GmlError):
            gml.read_gml(io.BytesIO(text.encode()), "clc18_kode")

    @pytest.mark.parametrize(
        "text",
        [
            pytest.param(b"", id="empty"),
            pytest.param(b"not xml at all", id="not-xml"),
            pytest.param(b'{"type": "FeatureCollection"}', id="json"),
        ],
    )
    def test_what_is_not_gml_is_refused(self, gml: ModuleType, text: bytes) -> None:
        with pytest.raises(gml.GmlError):
            gml.read_gml(io.BytesIO(text), "clc18_kode")

    def test_unparseable_coordinates_are_refused_naming_the_feature(self, gml: ModuleType) -> None:
        bad = member("f.bad", LAKE).replace("22.1,70.0", "22.1;70.0")
        with pytest.raises(gml.GmlError, match=r"f\.bad"):
            gml.read_gml(document(bad), "clc18_kode")

    def test_the_error_is_a_value_error(self, gml: ModuleType) -> None:
        assert issubclass(gml.GmlError, ValueError)


class TestLegacyFixture:
    """The committed file, repaired, read whole (30 MB; about a second)."""

    def test_the_committed_file_is_well_formed(self) -> None:
        root = ET.parse(LEGACY_GML).getroot()
        assert len(root.findall("{http://www.opengis.net/gml}featureMember")) == 200

    @pytest.fixture(scope="class")
    def legacy(self, gml: ModuleType) -> Any:
        with LEGACY_GML.open("rb") as stream:
            return gml.read_gml(stream, "clc18_kode")

    def test_every_feature_in_the_file_is_read(self, legacy: Any) -> None:
        text = LEGACY_GML.read_text(encoding="utf-8")
        fids = re.findall(r'<ogr:sql_statement fid="([^"]+)">', text)
        assert len(fids) == 200
        assert [f.fid for f in legacy.features] == fids

    def test_the_values_are_the_files_clc_codes(self, legacy: Any) -> None:
        text = LEGACY_GML.read_text(encoding="utf-8")
        codes = re.findall(r"<ogr:clc18_kode>(\d+)</ogr:clc18_kode>", text)
        assert [f.value for f in legacy.features] == codes

    def test_the_crs_is_wgs84_and_every_geometry_a_polygon(self, legacy: Any) -> None:
        assert legacy.crs == "EPSG:4326"
        assert len(legacy.features) == 200
        assert {f.geometry.geom_type for f in legacy.features} == {"Polygon"}

    def test_holes_are_read(self, legacy: Any) -> None:
        text = LEGACY_GML.read_text(encoding="utf-8")
        assert sum(len(f.geometry.interiors) for f in legacy.features) == text.count(
            "<gml:innerBoundaryIs>"
        )

    def test_forty_two_features_meet_the_benchmark_tile(self, legacy: Any) -> None:
        """Longitude first: read the other way round, no feature lands in
        Finnmark at all. The count is this file's, measured when the suite was
        written (`sql_statement.11029`, the sea polygon, among them)."""
        to_utm = Transformer.from_crs("EPSG:4326", "EPSG:25833", always_xy=True)
        hits = [
            f.fid
            for f in legacy.features
            if shapely.transform(f.geometry, lambda xy: _moved(to_utm, xy)).intersects(TILE_7908_3)
        ]
        assert len(hits) == 42
        assert "sql_statement.11029" in hits


class TestRepairScript:
    """`tests/fixtures/corine/repair_gml.py`: the committed file is the 2020
    file repaired, and the script runs on nothing else."""

    @pytest.fixture(scope="class")
    def repair(self) -> ModuleType:
        spec = importlib.util.spec_from_file_location("repair_gml", REPAIR_SCRIPT)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    @pytest.fixture(scope="class")
    def original(self) -> bytes:
        """The file as committed in 2020, from git; skipped outside a clone
        that has the commit (a shallow CI checkout, a source archive)."""
        try:
            done = subprocess.run(
                ["git", "show", f"8e30af4:{GML_2020_PATH}"],
                cwd=REPO,
                capture_output=True,
                check=True,
            )
        except (OSError, subprocess.CalledProcessError):
            pytest.skip("commit 8e30af4 is not in this checkout")
        return done.stdout

    def test_the_committed_file_is_the_2020_file_repaired(
        self, repair: ModuleType, original: bytes
    ) -> None:
        assert repair.repair(original) == LEGACY_GML.read_bytes()

    def test_the_repair_adds_the_two_closing_tags_and_nothing_else(
        self, repair: ModuleType, original: bytes
    ) -> None:
        before = original.decode().splitlines()
        after = repair.repair(original).decode().splitlines()
        added = list(difflib.ndiff(before, after))
        assert [d for d in added if d.startswith("- ")] == []
        assert [d for d in added if d.startswith("+ ")] == [
            "+   </gml:featureMember>",
            "+ </ogr:FeatureCollection>",
        ]
        assert after[1097] == "  </gml:featureMember>"
        assert "sql_statement.40965" in after[1092]

    def test_the_repaired_file_is_refused(self, repair: ModuleType) -> None:
        with pytest.raises(repair.RepairError, match="not the 2020 file"):
            repair.repair(LEGACY_GML.read_bytes())

    def test_a_file_differing_by_one_byte_is_refused(
        self, repair: ModuleType, original: bytes
    ) -> None:
        with pytest.raises(repair.RepairError, match="not the 2020 file"):
            repair.repair(original + b"\n")


def _moved(transformer: Transformer, xy: Any) -> Any:
    x, y = transformer.transform(xy[:, 0], xy[:, 1])
    return np.column_stack([x, y])


class TestBoundaries:
    def test_the_reader_opens_nothing_and_uses_only_the_standard_library_xml(
        self, gml: ModuleType
    ) -> None:
        path = SRC / "io" / "gml.py"
        tree = ast.parse(path.read_text())
        imported: set[str] = set()
        called: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported |= {a.name for a in node.names}
            elif isinstance(node, ast.ImportFrom):
                imported.add(node.module or "")
            elif isinstance(node, ast.Call):
                f = node.func
                called.add(f.id if isinstance(f, ast.Name) else getattr(f, "attr", ""))
        assert "open" not in called
        assert not any(m.split(".")[0] in {"lxml", "osgeo", "fiona", "pathlib"} for m in imported)
        assert any(m.startswith("xml.") for m in imported)
