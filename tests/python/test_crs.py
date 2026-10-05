"""`tin_engine.crs`: `parse_crs` and `reprojector`, the one transform site (increment 15b).

`docs/increments/15-dem-mosaic.md` R2 ([15b] `tin_engine/crs.py`), R9 and I8.
Pinned by this suite (see "Pinned by the red suite (15b)"):

- `parse_crs(text) -> pyproj.CRS` accepts anything `CRS.from_user_input`
  accepts, and refuses anything else with a `ValueError` naming the text.
  pyproj's own `CRSError` is a `RuntimeError`, not a `ValueError`, so letting it
  through would escape every `except ValueError` that turns a refusal into a
  usage error.
- `reprojector(src, dst)` takes two CRSs (text or `pyproj.CRS`) and returns a
  callable mapping an `(N, 2)` array-like of `(x, y)` to an `(N, 2)` float64
  array. `x` is easting or longitude, `y` northing or latitude, whatever the
  CRS's own axis order (`always_xy`), and the result is pyproj's own
  `always_xy` transform bit for bit.
- The only `Transformer.from_crs` call in `src_python/` is in
  `tin_engine/crs.py` (I8). `test_always_xy.py` still guards that the call
  passes a literal `always_xy=True`.

Real coordinates are checked against pyproj directly, never against numbers
typed into this file. The oracle is `Transformer.from_crs(..., always_xy=True)`
built here, in the test, which the grep guard does not scan.

HOW THIS FILE GOES RED: `tin_engine.crs` is imported inside a fixture, so each
test fails on its own with `ModuleNotFoundError` and collection is unaffected.

PR B of the Python audit (`docs/increments/python-audit.md`, section 9) adds
the one CRS rule, pinned by `TestSameCrs`, `TestTransformLabel` and
`TestSingleCrs`:

- `same_crs(a, b)`: the same once both are in x-then-y order (PROJ's
  equivalence), or else one EPSG code in common at PROJ's identify confidence
  70 ("equivalent, names differ"). So a CRS spelt as a PROJ string, a WKT
  without its ID or GDAL's WKT1 is the EPSG code it defines. Its limits: a
  `+towgs84` spelling is not the same (PROJ matches it on the ellipsoid
  alone); `OGC:CRS84` is `EPSG:4326`; `EPSG:3045` is `EPSG:25833`. A pair PROJ
  finds no operation between is not the same, never an exception.
- `transform_label(src, dst)`: "none" when the two are the same, else
  `transform_description(src, dst)`.
- `single_crs(texts, refusal=ValueError)`: the one text, or `refusal` naming
  them all in plain words.

Those three go red with `AttributeError`, each test on its own.
"""

from __future__ import annotations

import ast
import importlib
import warnings
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from pyproj import CRS, Transformer

SRC_PYTHON = Path(__file__).resolve().parents[2] / "src_python"

# Points in (longitude, latitude): the committed DTM10 seam extract near
# Flekkefjord, the benchmark tile's corner in Finnmark, Oslo, and the UTM 33
# central meridian. Longitude and latitude differ in every row, so swapping
# them cannot give the same answer.
LON_LAT = np.array(
    [
        [7.35061409316692, 58.12324881275789],
        [21.8, 71.1],
        [10.75, 59.91],
        [15.0, 59.5],
    ]
)


@pytest.fixture(scope="module")
def crs() -> ModuleType:
    return importlib.import_module("tin_engine.crs")


def oracle(src: str, dst: str, xy: np.ndarray) -> np.ndarray:
    t = Transformer.from_crs(src, dst, always_xy=True)
    x, y = t.transform(xy[:, 0], xy[:, 1])
    return np.column_stack([x, y])


class TestParseCrs:
    @pytest.mark.parametrize(
        ("text", "epsg"),
        [
            ("EPSG:25833", 25833),
            ("epsg:25833", 25833),
            ("urn:ogc:def:crs:EPSG::25833", 25833),
            ("EPSG:4326", 4326),
        ],
    )
    def test_epsg_spellings(self, crs: ModuleType, text: str, epsg: int) -> None:
        parsed = crs.parse_crs(text)
        assert isinstance(parsed, CRS)
        assert parsed == CRS.from_epsg(epsg)

    def test_a_crs_without_an_epsg_code_is_accepted(self, crs: ModuleType) -> None:
        """R9: `crs: str`, since a domain may have no EPSG code."""
        text = "+proj=lcc +lat_1=60 +lat_2=65 +lat_0=62 +lon_0=15 +ellps=GRS80 +units=m +no_defs"
        parsed = crs.parse_crs(text)
        assert parsed == CRS.from_user_input(text)
        assert parsed.to_epsg() is None

    def test_wkt2_is_accepted(self, crs: ModuleType) -> None:
        wkt = CRS.from_epsg(3035).to_wkt()
        assert crs.parse_crs(wkt) == CRS.from_epsg(3035)

    @pytest.mark.parametrize("text", ["EPSG:999999", "not a crs", "urn:ogc:def:crs:EPSG::0"])
    def test_an_unknown_crs_is_a_value_error_naming_it(self, crs: ModuleType, text: str) -> None:
        with pytest.raises(ValueError) as info:
            crs.parse_crs(text)
        assert text in str(info.value), info.value

    def test_empty_text_is_a_value_error(self, crs: ModuleType) -> None:
        with pytest.raises(ValueError):
            crs.parse_crs("")


class TestReprojector:
    def test_geographic_to_utm33_is_pyprojs_always_xy_bit_for_bit(self, crs: ModuleType) -> None:
        out = crs.reprojector("EPSG:4326", "EPSG:25833")(LON_LAT)
        assert isinstance(out, np.ndarray)
        assert out.dtype == np.float64
        assert out.shape == LON_LAT.shape
        assert out.tobytes() == oracle("EPSG:4326", "EPSG:25833", LON_LAT).tobytes()

    def test_axis_order_is_x_then_y_for_a_latitude_first_crs(self, crs: ModuleType) -> None:
        """EPSG:4326's own axis order is latitude first. The input here is
        (longitude, latitude), as in a GeoJSON file and a GeoTIFF's model space,
        and the central meridian point lands on UTM 33's false easting."""
        out = crs.reprojector("EPSG:4326", "EPSG:25833")(LON_LAT[3:])
        assert out[0, 0] == pytest.approx(500_000.0, abs=1e-6)
        assert 6_500_000.0 < out[0, 1] < 6_700_000.0
        # The same numbers read latitude first land thousands of km away: the
        # mutant this test exists to kill.
        swapped = Transformer.from_crs("EPSG:4326", "EPSG:25833").transform(*LON_LAT[3])
        assert abs(swapped[0] - out[0, 0]) > 1_000_000.0

    def test_output_is_x_then_y(self, crs: ModuleType) -> None:
        """Northing is the larger number everywhere in Norway, so a column swap
        on the way out cannot pass."""
        out = crs.reprojector("EPSG:4326", "EPSG:25833")(LON_LAT)
        assert (out[:, 1] > 6_000_000.0).all()
        assert (out[:, 0] < 1_200_000.0).all()

    @pytest.mark.parametrize(
        ("src", "dst"),
        [
            ("EPSG:25832", "EPSG:25833"),
            ("EPSG:3035", "EPSG:25833"),
            ("EPSG:25833", "EPSG:4326"),
            ("OGC:CRS84", "EPSG:25833"),
        ],
    )
    def test_other_pairs_are_pyprojs_always_xy(self, crs: ModuleType, src: str, dst: str) -> None:
        points = oracle("EPSG:4326", src, LON_LAT)
        out = crs.reprojector(src, dst)(points)
        assert out.tobytes() == oracle(src, dst, points).tobytes()

    def test_crs_objects_are_accepted(self, crs: ModuleType) -> None:
        out = crs.reprojector(CRS.from_epsg(4326), crs.parse_crs("EPSG:25833"))(LON_LAT)
        assert out.tobytes() == oracle("EPSG:4326", "EPSG:25833", LON_LAT).tobytes()

    def test_a_list_of_pairs_is_accepted(self, crs: ModuleType) -> None:
        pairs = [(float(x), float(y)) for x, y in LON_LAT]
        out = crs.reprojector("EPSG:4326", "EPSG:25833")(pairs)
        assert out.tobytes() == oracle("EPSG:4326", "EPSG:25833", LON_LAT).tobytes()

    def test_no_vertex_is_added_or_dropped(self, crs: ModuleType) -> None:
        """R9: vertices only, not densified."""
        for n in (1, 2, 7):
            assert crs.reprojector("EPSG:4326", "EPSG:25833")(LON_LAT[:1].repeat(n, 0)).shape == (
                n,
                2,
            )


# ---------------------------------------------------------------- PR B: one CRS rule

#: EPSG:31287 (MGI / Austria Lambert) written by its parameters, as the
#: GeoKeys of the Austrian openDEM file give it.
LAMBERT = (
    "+proj=lcc +lat_1=46 +lat_2=49 +lat_0=47.5 +lon_0=13.33333333333333 "
    "+x_0=400000 +y_0=400000 +ellps=bessel +units=m +no_defs"
)
#: EPSG:31287 with the datum shift PROJ's own EPSG:31287 to WGS 84 uses.
TOWGS84 = "+towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232"
MARS = "+proj=longlat +a=3396190 +b=3376200"


def without_id(epsg: int) -> str:
    """The WKT2 of `epsg` with its own `ID` removed (its base CRS's and its
    parameters' stay), so the CRS itself is named by definition only."""
    doc = CRS.from_epsg(epsg).to_json_dict()
    del doc["id"]
    text = CRS.from_json_dict(doc).to_wkt()
    assert f'ID["EPSG",{epsg}]' not in text, text
    return text


def proj4_of(epsg: int) -> str:
    """pyproj's PROJ string of `epsg`, without its lossy-conversion warning."""
    with warnings.catch_warnings(action="ignore", category=UserWarning):
        return CRS.from_epsg(epsg).to_proj4()


SAME = {
    "Lambert by parameters": (LAMBERT, "EPSG:31287"),
    "WKT without ID": (without_id(31287), "EPSG:31287"),
    "GDAL WKT1, no axis order": (CRS.from_epsg(3035).to_wkt("WKT1_GDAL"), "EPSG:3035"),
    "PROJ string of 25833": (proj4_of(25833), "EPSG:25833"),
    "CRS84 is 4326 in x-then-y order": ("OGC:CRS84", "EPSG:4326"),
    "one definition, two codes": ("EPSG:3045", "EPSG:25833"),
    "control, case only": ("EPSG:31287", "epsg:31287"),
}
NOT_SAME = {
    "Lambert in US feet": (LAMBERT.replace("+units=m", "+units=us-ft"), "EPSG:31287"),
    "Lambert at 13.5 E": (LAMBERT.replace("13.33333333333333", "13.5"), "EPSG:31287"),
    "another UTM zone": ("EPSG:25832", "EPSG:25833"),
    "another datum": ("EPSG:4258", "EPSG:4326"),
    "the +towgs84 limit": (LAMBERT.replace("+no_defs", f"{TOWGS84} +no_defs"), "EPSG:31287"),
    "Mars, no operation": (MARS, "EPSG:4326"),
}


def both_orders(pairs: dict[str, tuple[str, str]]) -> list[Any]:
    return [
        pytest.param(*pair, id=f"{name}, {order}")
        for name, (a, b) in pairs.items()
        for order, pair in (("given", (a, b)), ("swapped", (b, a)))
    ]


class TestSameCrs:
    @pytest.mark.parametrize(("a", "b"), both_orders(SAME))
    def test_the_same_points(self, crs: ModuleType, a: str, b: str) -> None:
        assert crs.same_crs(a, b) is True

    @pytest.mark.parametrize(("a", "b"), both_orders(NOT_SAME))
    def test_not_the_same(self, crs: ModuleType, a: str, b: str) -> None:
        assert crs.same_crs(a, b) is False

    def test_crs_objects_are_accepted(self, crs: ModuleType) -> None:
        assert crs.same_crs(CRS.from_epsg(31287), CRS.from_user_input(LAMBERT)) is True
        assert crs.same_crs(CRS.from_epsg(25832), "EPSG:25833") is False

    @pytest.mark.parametrize(("a", "b"), [("not a crs", "EPSG:4326"), ("EPSG:4326", "not a crs")])
    def test_unreadable_text_is_parse_crss_value_error(
        self, crs: ModuleType, a: str, b: str
    ) -> None:
        with pytest.raises(ValueError, match=r"cannot read the CRS 'not a crs'"):
            crs.same_crs(a, b)


class TestTransformLabel:
    @pytest.mark.parametrize(
        ("src", "dst"), both_orders({"Lambert by parameters": (LAMBERT, "EPSG:31287")})
    )
    def test_none_for_the_same_crs_by_definition(self, crs: ModuleType, src: str, dst: str) -> None:
        assert crs.transform_label(src, dst) == "none"

    def test_proj_s_description_for_a_real_transform(self, crs: ModuleType) -> None:
        label = crs.transform_label("EPSG:4326", "EPSG:25833")
        assert label == crs.transform_description("EPSG:4326", "EPSG:25833")
        assert label != "none"


class TestSingleCrs:
    MIXED = ("EPSG:25833", "EPSG:25832", "EPSG:25833")
    WORDING = (
        "the DEM files are in 2 different CRSs (EPSG:25832, EPSG:25833); all must be in one CRS"
    )

    class RefusalError(ValueError):
        pass

    def test_one_crs_is_returned(self, crs: ModuleType) -> None:
        assert crs.single_crs(["EPSG:25833"] * 3) == "EPSG:25833"

    def test_any_iterable_is_read_once(self, crs: ModuleType) -> None:
        """The sites pass a generator over the footprints."""
        assert crs.single_crs(t for t in ["EPSG:25833"] * 3) == "EPSG:25833"

    def test_two_crss_raise_the_given_type_in_plain_words(self, crs: ModuleType) -> None:
        with pytest.raises(self.RefusalError) as info:
            crs.single_crs(iter(self.MIXED), self.RefusalError)
        assert str(info.value) == self.WORDING

    def test_two_crss_raise_value_error_by_default(self, crs: ModuleType) -> None:
        with pytest.raises(ValueError) as info:
            crs.single_crs(list(self.MIXED))
        assert type(info.value) is ValueError
        assert str(info.value) == self.WORDING


# ---------------------------------------------------------------- I8


def from_crs_sites(root: Path) -> list[str]:
    """Every call whose callee is named `from_crs`, as `relative/path.py:line`."""
    found: list[str] = []
    for path in sorted(root.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                callee = node.func
                name = callee.attr if isinstance(callee, ast.Attribute) else None
                name = callee.id if isinstance(callee, ast.Name) else name
                if name == "from_crs":
                    found.append(f"{path.relative_to(root).as_posix()}:{node.lineno}")
    return found


def test_exactly_one_from_crs_site_and_it_is_in_crs_py() -> None:
    """I8 and R9: `crs.reprojector` holds the only `Transformer.from_crs` in
    `src_python/`. `test_always_xy.py` checks its keywords."""
    sites = from_crs_sites(SRC_PYTHON)
    assert len(sites) == 1, sites
    assert sites[0].startswith("tin_engine/crs.py:"), sites


def test_the_site_finder_finds_planted_sites(tmp_path: Path) -> None:
    (tmp_path / "a.py").write_text(
        "from pyproj import Transformer\n"
        "t = Transformer.from_crs('EPSG:4326', 'EPSG:25833', always_xy=True)\n"
        "u = from_crs(1, 2)\n",
        encoding="utf-8",
    )
    (tmp_path / "b.py").write_text("x = 1\n", encoding="utf-8")
    assert from_crs_sites(tmp_path) == ["a.py:2", "a.py:3"]
