"""`land_cover_map.build(..., preferred, ...)` then `edges_of`: the land
cover's own polygons take the preferred classes before they become lines,
so a lake has one shoreline, the preferred source's.

`docs/increments/110-prefer-lines.md`, *3. Tests*, test 2. The scene is
`test_preferred.py`'s case 1 in a 100 m square domain at plain coordinates in
EPSG:25833: heath (322) over the square less [40,60]², water (512) on
[40,60]², the preferred lake [30,50] x [40,60]. The source is read with
``corine`` (`read_class_polygons`), the map built with the lake preferred,
and every class chosen (`docs/increments/111-one-land-cover-map.md`,
section 8: `Prefer` and `open_features(prefer=)` become `build`'s
``preferred``).

THE ORACLE: the source's lines, unioned, less the domain outline, against the
preferred lake's boundary (80 m) and the old water border off it (60 m: x = 40
inside the lake and the part beyond x = 50). Every corner is on a whole metre
at most 100 m from the origin, and with repair and simplification at 0 the
clean-up moves no vertex, so lengths and distances are exact up to rounding:
bound `EPS` = 1e-9 m. A crossing point of two lines gives `LINE` (1e-6 m) of
overlap through a band of `EPS`.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- The given lake's feature is the one feature of code 512, with the mask
  ``land_cover`` and ``water``.
- The third way, "repair 0, merge off, tolerance 0", runs with outline snap
  5 m (the CLI's default), so the clean-up stage is asked for and differs
  from the first way; the lake lies 30 m from the outline, beyond its reach.
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import shapely
from shapely.geometry import box

import feature_fixtures as ff
from feature_fixtures import UTM33, Feat, domain_of, write_geojson
from tin_engine.features import DEFAULT_VOCABULARY

EPS = 1e-9
LINE = 1e-6
HEATH, LAKE = 322, 512
SQUARE = box(0, 0, 100, 100)
OLD = box(40, 40, 60, 60)
PREFERRED = box(30, 40, 50, 60)
WAYS = {
    "cleanup-off": {},
    "merge-on": {"repair_m": 0.0, "merge_same_class": True, "max_shift_m": 0.0},
    "merge-off": {"repair_m": 0.0, "merge_same_class": False, "max_shift_m": 0.0,
                  "outline_snap_m": 5.0},
}  # fmt: skip


@pytest.fixture(scope="module")
def fi() -> ModuleType:
    return ff.feature_input()


def opened(fi: ModuleType, tmp_path: Path, way: dict[str, Any]) -> Any:
    shapes = [(HEATH, SQUARE.difference(OLD)), (LAKE, OLD)]
    feats = [Feat(i + 1, g, {"Code_18": str(c)}) for i, (c, g) in enumerate(shapes)]
    path = write_geojson(tmp_path / "cover.geojson", feats)
    domain = domain_of(SQUARE)
    corine, region = ff.classes("corine"), ff.region_of(domain)
    read = fi.read_class_polygons(path, None, corine, UTM33, ff.box_of(domain))
    built = ff.land_cover_map().build([read.areas], [(LAKE, [PREFERRED])], corine, region)
    return fi.edges_of(
        built, [], ff.every_class(), ff.cleanup(**way), domain, DEFAULT_VOCABULARY, region=region
    )


@pytest.mark.parametrize("way", WAYS.values(), ids=list(WAYS))
def test_the_lines_are_the_preferred_shoreline(
    fi: ModuleType, tmp_path: Path, way: dict[str, Any]
) -> None:
    fs = opened(fi, tmp_path, way)
    lines = shapely.union_all([line for f in fs.features for line in f.lines])
    off = lines.difference(SQUARE.boundary)
    assert off.length == pytest.approx(80.0, abs=EPS), off.wkt
    xy = shapely.get_coordinates(off)
    assert shapely.distance(PREFERRED.boundary, shapely.points(xy)).max() <= EPS, off.wkt
    old = OLD.boundary.difference(PREFERRED.boundary)
    assert old.length == pytest.approx(60.0, abs=EPS)  # the premise
    left = old.intersection(lines.buffer(EPS))
    assert left.length <= LINE, left.wkt
    water = DEFAULT_VOCABULARY.mask("land_cover", "water")
    (lake,) = [f for f in fs.features if f.code == LAKE]
    assert lake.mask == water, lake
    ring = shapely.get_coordinates(shapely.union_all(list(lake.lines)))
    assert shapely.distance(PREFERRED.boundary, shapely.points(ring)).max() <= EPS
