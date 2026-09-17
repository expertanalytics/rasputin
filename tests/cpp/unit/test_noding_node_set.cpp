// Unit tests for terrain/noding/node_set.hpp.
//
// NOT invariant-critical, and the increment file says so explicitly: this is
// sort, unique and lower_bound over a totally ordered integral key with a
// defaulted operator<=>. There is no topology decision inside it -- the
// decision "equal keys collapse" is made by the key, and the key is 5a's
// snapping, which is where the mutation budget goes. The failure mode here is a
// typo or an off-by-one at a boundary, which an exhaustive small-input suite
// catches. Do not spend a mutation round on this file.
//
// It still owes one adversarial test, and it is the last one: two points
// straddling a cell boundary by a single ulp must remain two nodes. That is the
// case where a reader expects a collapse.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/node_set.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

using terrain::GridPoint;
using terrain::Point2;
using terrain::SnapGrid;
using terrain::noding::NodeSet;

namespace {

[[nodiscard]] std::vector<GridPoint> collect(const NodeSet& set) {
    return std::vector<GridPoint>{set.points().begin(), set.points().end()};
}

}  // namespace

TEST_CASE("an empty NodeSet has no nodes", "[noding][node_set]") {
    const NodeSet set{std::vector<GridPoint>{}};

    REQUIRE(set.size() == 0);
    REQUIRE(set.points().empty());
}

TEST_CASE("a one-point NodeSet holds exactly that point", "[noding][node_set]") {
    const NodeSet set{std::vector<GridPoint>{GridPoint{7, -3}}};

    REQUIRE(set.size() == 1);
    REQUIRE(set.points()[0] == GridPoint{7, -3});
    REQUIRE(set.id_of(GridPoint{7, -3}) == 0);
    REQUIRE(set[0] == GridPoint{7, -3});
}

TEST_CASE("NodeSet sorts lexicographically and dedups", "[noding][node_set]") {
    const NodeSet set{std::vector<GridPoint>{
        GridPoint{2, 5}, GridPoint{-1, 9}, GridPoint{2, -4},
        GridPoint{2, 5}, GridPoint{-1, 9}, GridPoint{-1, 9},
        GridPoint{0, 0},
    }};

    const std::vector<GridPoint> expected{
        GridPoint{-1, 9}, GridPoint{0, 0}, GridPoint{2, -4}, GridPoint{2, 5},
    };

    REQUIRE(set.size() == 4);
    REQUIRE(collect(set) == expected);
    REQUIRE(std::ranges::is_sorted(set.points()));
}

// The ordering is a function of the SET of points and of nothing else -- not of
// input order, not of which thread found a crossing first. This is what makes
// the node numbering, and every downstream mesh index, reproducible across runs
// and across thread counts. Mutant 10 is first-appearance order, and the
// shuffle property in prop_noding_snap_invariants.cpp is what kills it over
// generated input; this case states the rule on a fixture a reader can see.
TEST_CASE("node ids are in grid order, not first-appearance order", "[noding][node_set]") {
    const NodeSet set{std::vector<GridPoint>{
        GridPoint{9, 9}, GridPoint{1, 1}, GridPoint{5, 5},
    }};

    REQUIRE(set.id_of(GridPoint{1, 1}) == 0);
    REQUIRE(set.id_of(GridPoint{5, 5}) == 1);
    REQUIRE(set.id_of(GridPoint{9, 9}) == 2);
}

TEST_CASE("id_of and operator[] are inverse on every node", "[noding][node_set]") {
    std::vector<GridPoint> points;
    for (std::int64_t ix = -3; ix <= 3; ++ix) {
        for (std::int64_t iy = -3; iy <= 3; ++iy) {
            points.push_back(GridPoint{ix, iy});
            points.push_back(GridPoint{ix, iy});  // every key duplicated
        }
    }
    const NodeSet set{points};

    REQUIRE(set.size() == 49);
    for (std::uint32_t id = 0; id < set.size(); ++id) {
        INFO("id = " << id);
        REQUIRE(set.id_of(set[id]) == id);
    }
}

// The boundaries of the binary search: the first and last keys, and keys that
// are adjacent in the ordering, are where an off-by-one lives.
TEST_CASE("id_of finds the extreme and adjacent keys", "[noding][node_set]") {
    constexpr std::int64_t big = std::int64_t{1} << 51;
    const NodeSet set{std::vector<GridPoint>{
        GridPoint{-big, -big}, GridPoint{big, big}, GridPoint{0, -1},
        GridPoint{0, 0}, GridPoint{0, 1}, GridPoint{-1, big},
    }};

    REQUIRE(set.size() == 6);
    REQUIRE(set.id_of(GridPoint{-big, -big}) == 0);
    REQUIRE(set.id_of(GridPoint{big, big}) == set.size() - 1);
    REQUIRE(set.id_of(GridPoint{0, -1}) + 1 == set.id_of(GridPoint{0, 0}));
    REQUIRE(set.id_of(GridPoint{0, 0}) + 1 == set.id_of(GridPoint{0, 1}));
}

// Dedup is exact integer equality and nothing else: a key differing by one in
// either coordinate is a different node, however close the two coordinates are
// in metres.
TEST_CASE("dedup collapses only exactly equal keys", "[noding][node_set]") {
    const NodeSet set{std::vector<GridPoint>{
        GridPoint{5000000, 79000000},
        GridPoint{5000000, 79000000},
        GridPoint{5000001, 79000000},
        GridPoint{5000000, 79000001},
    }};

    REQUIRE(set.size() == 3);
}

// THE ADVERSARIAL ONE. Two coordinates one ulp apart, either side of the
// boundary between cells 5000000 and 5000001 at 0.1 m spacing. A reader expects
// points this close to collapse; they must not, because they snap to different
// keys, and the snapping -- not the NodeSet -- is what decides coincidence.
//
// PROVENANCE: the pair was found by walking nextafter across the half-index
// coordinate (5000000 + 1/2) * 0.1 in the UTM33 easting band.
TEST_CASE("points one ulp apart across a cell boundary stay two nodes", "[noding][node_set]") {
    const SnapGrid grid{0.1};
    const Point2 under{500000.04999999993, 7900000.0};
    const Point2 over{500000.05, 7900000.0};

    REQUIRE(std::nextafter(under.x, std::numeric_limits<double>::infinity()) == over.x);

    const NodeSet set{std::vector<GridPoint>{grid.snap(under), grid.snap(over)}};

    REQUIRE(set.size() == 2);
    REQUIRE(set.id_of(grid.snap(under)) != set.id_of(grid.snap(over)));
}

// The other half of the same statement: two points far enough apart to look
// distinct, inside one cell, are one node.
TEST_CASE("distinct points inside one cell are one node", "[noding][node_set]") {
    const SnapGrid grid{0.1};
    const Point2 a{500000.041, 7900000.044};
    const Point2 b{500000.038, 7900000.037};

    REQUIRE(a != b);

    const NodeSet set{std::vector<GridPoint>{grid.snap(a), grid.snap(b)}};

    REQUIRE(set.size() == 1);
}
