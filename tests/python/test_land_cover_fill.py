"""`landcover.fill_no_class`, `land_cover_fractions.merge_zero_share` and
`landcover.label_triangles(..., rest=)`: no triangle of a finished mesh keeps
code 0, and code 0 means "in no polygon of the map".

`docs/increments/122-lakes-overrule.md`, section 3.2 (the fill's rule: groups,
donors, nearest, ties), section 3.3 (the shares merge), section 3.7
(`rest=`), section 6 (the degeneracies) and section 7, tests 1 to 6, 12 and
23. Each grid is unit squares, `c` columns by 2 rows, each square two
triangles, codes per column (CORINE: 311 broad-leaved forest, 322 heath, 333
sparse, 511 river, 512 lake; water is the class file's `water_codes`).

THE ORACLE is the design's own figures, worked by hand on the grid: a
column's triangles are at distance `|i - j| - 1` from column `j`'s, and two
touching columns share an edge of length 2 (or 1, a row each). Distances and
shared lengths are whole numbers in the grid's unit, exact in float64, so
codes compare exactly; the one area compared (`largest`) is a sum of halves,
exact too. The UTM-scale case shifts the grid by (136216, 6830714), the
largest coordinates checked: float64 spacing there is 1e-9, far under the
unit distances compared.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- `Filled.codes` is int32 and a new array; `codes` given is not changed.
- With nothing to fill, `largest` is `(0.0, nan, nan)` and `groups` 0.
- `asking` replaces the default (`codes == 0`) whole: an asking triangle of
  a nonzero code is filled too, and is not a donor.
- `merge_zero_share` takes `into` as an int array of one code per triangle;
  a row with no 0 share is unchanged; the merged rows keep classes
  ascending.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from shapely.geometry import box

from tin_engine import land_cover_classes
from tin_engine.io.models import ClassFractions

HEATH, SPARSE, RIVER, LAKE, ARABLE = 322, 333, 511, 512, 211
#: CORINE broad-leaved forest, test 23's middle class.
BROADLEAF = 311
WATER = land_cover_classes.load("corine").water_codes
#: The step-2 refusal's words (section 3.2).
NO_DONOR = (
    "the land cover gives no triangle of the mesh a class, so code 0 cannot be "
    "filled; it does not cover the domain"
)

Column = int | tuple[int, int]


def grid(
    columns: Sequence[Column], origin: tuple[float, float] = (0.0, 0.0)
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vertices, triangles and codes of `len(columns)` by 2 unit squares, two
    triangles per square; a column is one code, or (below, above)."""
    c = len(columns)
    ix, iy = np.meshgrid(np.arange(c + 1), np.arange(3), indexing="ij")
    vertices = np.column_stack([ix.ravel(), iy.ravel()]).astype(np.float64)
    vertices += np.asarray(origin)

    def v(i: int, j: int) -> int:
        return i * 3 + j

    triangles, codes = [], []
    for i, column in enumerate(columns):
        below, above = column if isinstance(column, tuple) else (column, column)
        for j, code in ((0, below), (1, above)):
            a, b, cc, d = v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)
            triangles += [[a, b, cc], [a, cc, d]]
            codes += [code, code]
    return vertices, np.array(triangles, np.int64), np.array(codes, np.int32)


def column_of(triangles: np.ndarray) -> np.ndarray:
    """Each triangle's column in `grid`'s order (four per column)."""
    return np.arange(len(triangles)) // 4


def fill(columns: Sequence[Column], **kwargs: object) -> object:
    from tin_engine.landcover import fill_no_class

    return fill_no_class(*grid(columns), water=WATER, **kwargs)


# ------------------------------------------------------------ tests 1 to 5


def test_1_one_group_takes_the_nearest_non_water_class() -> None:
    """[333, 512, 0, 512, 512, 322]: the group (column 2) is 1 from 333 and 2
    from 322, across water."""
    v, t, codes = grid([SPARSE, LAKE, 0, LAKE, LAKE, HEATH])
    from tin_engine.landcover import Filled, fill_no_class

    got = fill_no_class(v, t, codes, water=WATER)
    assert isinstance(got, Filled)
    want = codes.copy()
    want[column_of(t) == 2] = SPARSE
    assert got.codes.dtype == np.int32
    assert_array_equal(got.codes, want)
    assert_array_equal(got.filled, column_of(t) == 2)
    assert got.groups == 1
    assert not (codes != grid([SPARSE, LAKE, 0, LAKE, LAKE, HEATH])[2]).any(), "codes changed"


def test_1_largest_is_the_groups_area_at_its_largest_triangles_centroid() -> None:
    """The group's 4 triangles of 1/2 each: area 2; all equal, so the lowest
    index, triangle 8 = (2,0), (3,0), (3,1): centroid (8/3, 1/3)."""
    got = fill([SPARSE, LAKE, 0, LAKE, LAKE, HEATH])
    area, x, y = got.largest
    assert area == 2.0
    assert (x, y) == pytest.approx((8 / 3, 1 / 3), abs=1e-12)


def test_2_a_tie_at_one_distance_takes_the_lower_code() -> None:
    """[333, 512, 0, 512, 322]: 333 and 322 both at 1, neither shares an edge."""
    _, t, _ = grid([SPARSE, LAKE, 0, LAKE, HEATH])
    got = fill([SPARSE, LAKE, 0, LAKE, HEATH])
    assert set(got.codes[column_of(t) == 2]) == {HEATH}


def test_3_the_longest_shared_border_beats_the_lower_code() -> None:
    """[333, 0, (322 below, 512 above)]: both donors touch the group (0);
    333 shares 2, 322 shares 1; 512 is water. 333, over the lower 322."""
    _, t, _ = grid([SPARSE, 0, (HEATH, LAKE)])
    got = fill([SPARSE, 0, (HEATH, LAKE)])
    assert set(got.codes[column_of(t) == 1]) == {SPARSE}


def test_4_water_only_donors_give_water_marked() -> None:
    """[512, 0]: no non-water donor, so every coded triangle is one."""
    _, t, _ = grid([LAKE, 0])
    got = fill([LAKE, 0])
    assert_array_equal(got.codes, [LAKE] * 8)
    assert_array_equal(got.filled, column_of(t) == 1)


def test_5_nothing_to_fill_from_is_refused_with_the_step_2_words() -> None:
    with pytest.raises(ValueError) as caught:
        fill([0, 0])
    assert str(caught.value) == NO_DONOR


# ------------------------------------------------------------ test 6


def test_6_two_groups_each_take_their_own_class() -> None:
    """[0, 333, 512, 0, 322]: the west group touches 333; the east touches
    322 and is 2 from 333."""
    _, t, codes = grid([0, SPARSE, LAKE, 0, HEATH])
    got = fill([0, SPARSE, LAKE, 0, HEATH])
    col = column_of(t)
    assert set(got.codes[col == 0]) == {SPARSE}
    assert set(got.codes[col == 3]) == {HEATH}
    assert got.filled.sum() == 8
    assert_array_equal(got.filled, codes == 0)
    assert got.groups == 2


def test_6_largest_of_two_equal_groups_is_the_one_with_the_lowest_index() -> None:
    """Both groups 2 in area: the west holds triangle 0, (0,0), (1,0), (1,1)."""
    area, x, y = fill([0, SPARSE, LAKE, 0, HEATH]).largest
    assert area == 2.0
    assert (x, y) == pytest.approx((2 / 3, 1 / 3), abs=1e-12)


CASES = [
    [SPARSE, LAKE, 0, LAKE, LAKE, HEATH],
    [SPARSE, LAKE, 0, LAKE, HEATH],
    [SPARSE, 0, (HEATH, LAKE)],
    [LAKE, 0],
    [0, SPARSE, LAKE, 0, HEATH],
]


@pytest.mark.parametrize("columns", CASES, ids=["1", "2", "3", "4", "6"])
def test_6_no_order_matters_triangles_vertices_and_corners_permuted(
    columns: Sequence[Column],
) -> None:
    """Section 3.2: the result depends on no triangle or vertex order (seed 0;
    each triangle's corners rolled too, which keeps its orientation)."""
    from tin_engine.landcover import fill_no_class

    v, t, codes = grid(columns)
    want = fill_no_class(v, t, codes, water=WATER)
    rng = np.random.default_rng(0)
    order, relabel = rng.permutation(len(t)), rng.permutation(len(v))
    moved = np.empty_like(relabel)
    moved[relabel] = np.arange(len(v))  # old vertex k is now vertex relabel[k]
    vertices = v[moved]
    triangles = relabel[t][order]
    triangles = np.array(
        [np.roll(tri, int(k)) for tri, k in zip(triangles, order % 3, strict=True)]
    )
    got = fill_no_class(vertices, triangles, codes[order], water=WATER)
    assert_array_equal(got.codes, want.codes[order])
    assert_array_equal(got.filled, want.filled[order])
    assert got.groups == want.groups


def test_5_nothing_to_fill_from_is_refused_in_any_order() -> None:
    from tin_engine.landcover import fill_no_class

    v, t, codes = grid([0, 0])
    with pytest.raises(ValueError, match="cannot be filled"):
        fill_no_class(v[::-1], (len(v) - 1 - t)[::-1], codes, water=WATER)


# ------------------------------------------------------------ section 6's degeneracies


def test_no_code_0_the_fill_does_nothing_and_marks_nothing() -> None:
    _, _, codes = grid([SPARSE, LAKE, HEATH])
    got = fill([SPARSE, LAKE, HEATH])
    assert_array_equal(got.codes, codes)
    assert not got.filled.any()
    assert got.groups == 0
    area, x, y = got.largest
    assert area == 0.0 and np.isnan(x) and np.isnan(y)


def test_without_a_class_file_water_is_empty_and_every_code_donates() -> None:
    """Section 6: `water` empty, so the nearer 512 donates."""
    from tin_engine.landcover import fill_no_class

    v, t, codes = grid([HEATH, LAKE, 0])
    got = fill_no_class(v, t, codes, water=frozenset())
    assert set(got.codes[column_of(t) == 2]) == {LAKE}


def test_a_donor_touching_only_at_a_vertex_is_at_distance_0() -> None:
    """Two columns, codes per square: the asking square (0,0) touches 322 at
    (1,1) only; the other squares are water. Distance 0, shared length 0."""
    from tin_engine.landcover import fill_no_class

    v, t, codes = grid([(0, LAKE), (LAKE, HEATH)])
    got = fill_no_class(v, t, codes, water=WATER)
    assert set(got.codes[:2]) == {HEATH}


def test_an_edge_beats_a_vertex_at_distance_0_though_its_code_is_higher() -> None:
    """As above, but the square east of the asking one is 333: both touch it
    (0); 333 shares 1, 322 shares 0. 333, over the lower 322."""
    from tin_engine.landcover import fill_no_class

    v, t, codes = grid([(0, LAKE), (SPARSE, HEATH)])
    got = fill_no_class(v, t, codes, water=WATER)
    assert set(got.codes[:2]) == {SPARSE}


def test_a_nearly_equal_distance_is_no_tie() -> None:
    """Test 2's grid with the line x = 4 moved east by 1e-9: 322 is then
    1 + 1e-9 away against 333's 1, so 333, though 322 is lower (exact float
    distance, section 3.2 step 3)."""
    from tin_engine.landcover import fill_no_class

    v, t, codes = grid([SPARSE, LAKE, 0, LAKE, HEATH])
    v = v.copy()
    v[v[:, 0] == 4.0, 0] += 1e-9
    got = fill_no_class(v, t, codes, water=WATER)
    assert set(got.codes[column_of(t) == 2]) == {SPARSE}


@pytest.mark.parametrize("columns", CASES, ids=["1", "2", "3", "4", "6"])
def test_utm_scale_coordinates_give_the_same_codes(columns: Sequence[Column]) -> None:
    """Shifted to (136216, 6830714) in a projected CRS: the same codes."""
    from tin_engine.landcover import fill_no_class

    v, t, codes = grid(columns)
    want = fill_no_class(v, t, codes, water=WATER)
    got = fill_no_class(*grid(columns, origin=(136216.0, 6830714.0)), water=WATER)
    assert_array_equal(got.codes, want.codes)
    assert_array_equal(got.filled, want.filled)


def test_asking_fills_the_triangles_it_marks_whatever_their_code() -> None:
    """[333, 512, 322] asking the water column: it touches 333 and 322 with 2
    each, so the lower, 322; the asking triangles are no donors."""
    from tin_engine.landcover import fill_no_class

    v, t, codes = grid([SPARSE, LAKE, HEATH])
    asking = column_of(t) == 1
    got = fill_no_class(v, t, codes, water=WATER, asking=asking)
    want = codes.copy()
    want[asking] = HEATH
    assert_array_equal(got.codes, want)
    assert_array_equal(got.filled, asking)


# ------------------------------------------------------------ test 12, the shares


def table(rows: Sequence[dict[int, float]]) -> ClassFractions:
    offsets = np.cumsum([0, *(len(r) for r in rows)]).astype(np.uint32)
    classes = np.array([k for r in rows for k in sorted(r)], np.int32)
    fractions = np.array([r[k] for r in rows for k in sorted(r)], np.float32)
    return ClassFractions(offsets=offsets, classes=classes, fractions=fractions)


def test_12_a_zero_share_joins_the_class_it_is_given() -> None:
    from tin_engine.land_cover_fractions import merge_zero_share

    given = table([{0: 0.25, LAKE: 0.75}, {0: 0.4, HEATH: 0.6}])
    got = merge_zero_share(given, np.array([HEATH, HEATH], np.int32))
    assert isinstance(got, ClassFractions)
    assert_array_equal(got.offsets, np.array([0, 2, 3], np.uint32))
    assert_array_equal(got.classes, np.array([HEATH, LAKE, HEATH], np.int32))
    assert got.fractions.dtype == np.float32
    assert_array_equal(got.fractions, np.array([0.25, 0.75, 1.0], np.float32))


def test_12_a_table_with_no_zero_share_comes_back_equal() -> None:
    from tin_engine.land_cover_fractions import merge_zero_share

    given = table([{HEATH: 0.5, LAKE: 0.5}, {SPARSE: 1.0}])
    got = merge_zero_share(given, np.array([HEATH, SPARSE], np.int32))
    for name in ("offsets", "classes", "fractions"):
        assert_array_equal(getattr(got, name), getattr(given, name))
        assert getattr(got, name).dtype == getattr(given, name).dtype


def test_12_a_row_all_zero_becomes_its_class_alone() -> None:
    from tin_engine.land_cover_fractions import merge_zero_share

    got = merge_zero_share(table([{0: 1.0}, {SPARSE: 1.0}]), np.array([ARABLE, SPARSE], np.int32))
    assert_array_equal(got.offsets, np.array([0, 1, 2], np.uint32))
    assert_array_equal(got.classes, np.array([ARABLE, SPARSE], np.int32))
    assert_array_equal(got.fractions, np.array([1.0, 1.0], np.float32))


# ------------------------------------------------------------ test 23, label_triangles(rest=)

#: Test 23's polygons: one region east of x = 3, the map's areas over all four columns.
EAST_LAKE = [(box(3, 0, 4, 2), LAKE)]
MAP_AREAS = [(box(0, 0, 2, 2), HEATH), (box(2, 0, 3, 2), BROADLEAF), (box(3, 0, 4, 2), LAKE)]


def labelled(**rest: object) -> object:
    """`label_triangles` on the 4 by 2 grid, the two unit segments on x = 3
    the constraint edges, `polygons` the lake east of them."""
    from tin_engine.landcover import label_triangles

    vertices, triangles, _ = grid([0, 0, 0, 0])
    edges = np.array([[9, 10], [10, 11]], np.int64)  # (3, 0)-(3, 1)-(3, 2)
    return label_triangles(vertices, triangles, edges, polygons=EAST_LAKE, margin=1e-3, **rest)


def test_23_the_label_points_pin_the_figures() -> None:
    """The figures below are by construction: each half-square's incentre
    lies strictly inside its own unit square, so off every polygon border."""
    vertices, triangles, _ = grid([0, 0, 0, 0])
    a, b, c = (vertices[triangles[:, k]] for k in range(3))
    la, lb, lc = (np.hypot(*(q - p).T) for p, q in ((b, c), (c, a), (a, b)))
    centre = (la[:, None] * a + lb[:, None] * b + lc[:, None] * c) / (la + lb + lc)[:, None]
    column = column_of(triangles)
    assert ((centre[:, 0] > column) & (centre[:, 0] < column + 1)).all()
    assert ((centre[:, 1] % 1 > 0) & (centre[:, 1] % 1 < 1)).all()


def test_23_a_region_the_borders_leave_0_takes_the_map_per_triangle() -> None:
    """Section 3.7: the west region holds heath and broad-leaved forest with no
    line between them; each triangle takes the map's code at its own label
    point, so two classes in one region, and `CoverLabels.rest` counts 12."""
    got = labelled(rest=MAP_AREAS)
    column = column_of(np.zeros((16, 3)))
    assert_array_equal(got.codes[column <= 1], np.full(8, HEATH, np.int32))
    assert_array_equal(got.codes[column == 2], np.full(4, BROADLEAF, np.int32))
    assert_array_equal(got.codes[column == 3], np.full(4, LAKE, np.int32))
    assert got.rest == 12


def test_23_a_triangle_in_no_polygon_of_the_map_stays_0() -> None:
    got = labelled(rest=[a for a in MAP_AREAS if a[1] != BROADLEAF])
    column = column_of(np.zeros((16, 3)))
    assert_array_equal(got.codes[column <= 1], np.full(8, HEATH, np.int32))
    assert_array_equal(got.codes[column == 2], np.zeros(4, np.int32))
    assert_array_equal(got.codes[column == 3], np.full(4, LAKE, np.int32))
    assert got.rest == 8


def test_23_without_rest_the_labels_are_todays() -> None:
    from tin_engine.landcover import label_triangles

    vertices, triangles, _ = grid([0, 0, 0, 0])
    edges = np.array([[9, 10], [10, 11]], np.int64)
    today = label_triangles(vertices, triangles, edges, polygons=EAST_LAKE, margin=1e-3)
    got = labelled(rest=())
    assert_array_equal(got.codes, today.codes)
    column = column_of(np.zeros((16, 3)))
    assert_array_equal(got.codes, np.where(column == 3, LAKE, 0).astype(np.int32))
    assert got.rest == 0
    assert (got.regions, got.outside, got.overlapped, got.thin) == (
        today.regions,
        today.outside,
        today.overlapped,
        today.thin,
    )
