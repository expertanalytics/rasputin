"""The clean-up's clip to the domain keeps the land cover a coverage.

``data/clip_noding_17.wkt`` holds three polygons cut from catchment 17 (CORINE
classes 6 and 12 after the outline rule, EPSG:25833): lines 1 and 2 are two
neighbours that tile, line 3 a piece of the domain. Clipped one at a time,
the neighbours' shared border stops matching vertex for vertex where the
outline rule's vertex and the clip's crossing lie about 1e-10 m apart, and
``shapely.coverage_union_all`` refuses the pair. The fixture is from
``docs/research/2026-10-29-outline-order.md``, *Why the clean-up crashes in 17, 23, 32 and 39*.

A test through the whole clean-up cannot use it: the repair and the outline
rule move its vertices first, and the clip then agrees. So the clip is
tested on its own.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import shapely

from tin_engine.feature_input import clip_to_domain

FIXTURE = Path(__file__).parent / "data" / "clip_noding_17.wkt"


@pytest.fixture(scope="module")
def three() -> list[shapely.Geometry]:
    lines = [line for line in FIXTURE.read_text().split("\n") if line]
    return [shapely.from_wkt(line) for line in lines]


def test_the_fixture_shows_the_fault(three: list[shapely.Geometry]) -> None:
    """The pair tiles and the plain one-by-one clip does not; if a later GEOS
    clips it cleanly, this fixture no longer tests anything."""
    a, b, dom = three
    assert [int(shapely.get_num_coordinates(g)) for g in three] == [30, 18, 25]
    assert shapely.coverage_is_valid(np.array([a, b], dtype=object))
    plain = np.array([shapely.intersection(g, dom) for g in (a, b)], dtype=object)
    assert not shapely.coverage_is_valid(plain)


def test_the_clipped_pair_tiles_and_unions(three: list[shapely.Geometry]) -> None:
    a, b, dom = three
    out = clip_to_domain([a, b], dom)
    assert len(out) == 2
    assert shapely.coverage_is_valid(np.array(out, dtype=object))
    shapely.coverage_union_all(np.array(out, dtype=object))  # raises if not noded


def test_the_clip_moves_next_to_no_area(three: list[shapely.Geometry]) -> None:
    """Each polygon keeps its plain clip's area within 1e-6 m2 (measured 2.6e-10)."""
    a, b, dom = three
    out = clip_to_domain([a, b], dom)
    for g, c in zip((a, b), out, strict=True):
        assert abs(shapely.intersection(g, dom).area - c.area) < 1e-6
        assert c.difference(dom).area < 1e-6


def test_a_polygon_outside_the_domain_is_kept_empty(three: list[shapely.Geometry]) -> None:
    """Its slot stays, empty, so the clipped list lines up with its input."""
    a, b, dom = three
    far = shapely.box(0.0, 0.0, 1.0, 1.0)
    out = clip_to_domain([a, far, b], dom)
    assert len(out) == 3
    assert out[1].is_empty and not out[0].is_empty and not out[2].is_empty
