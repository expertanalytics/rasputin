// Increment 23b (docs/increments/23-basin-scale.md, "The seam protocol" step 2,
// "Tolerance, the final check and the constraint check points per piece", K3,
// K4, the degeneracy policy, and "Tests @tester can write red", SP1 to SP4):
// the seam pass, refine_seam, through the C++ API, and a piece refined with its
// seam frozen after it. SP1 and SP2 are invariant-critical (mutation runs at
// green); the Python half of SP1 and SP2 (NumPy, fractions) is
// tests/python/test_core_seam.py.
//
// Interface, pinned at the red step and ruled by N8 to N15 ("seam.hpp:
// refine_seam (one-dimensional greedy over constraint_check_points)"; "Output:
// the inserted points in order from a to b, world (x, y) and z, plus z for a
// and b themselves"; "Settled after 23b's red step"):
//
//   #include <terrain/refinement/seam.hpp>
//   struct SeamPoint   { Point2 at; double z; double s; };   world, vertex_z, parameter from a
//   struct SeamOutcome {
//       Point2 a, b;                     // the edge's ends, ordered so a < b by (x, y)
//       std::optional<double> z_a, z_b;  // vertex_z at each end's lattice position; nullopt on NoData
//       std::vector<SeamPoint> points;   // inserted, s strictly increasing in (0, 1)
//       std::size_t check_points;        // constraint_check_points(...).size() for the one edge
//       std::size_t no_data;             // its no_data(): check points skipped for a NoData stencil
//       double max_error;                // largest |z_p - lerp(p)| over check points at the end,
//   };                                   // pieces with an invalid end excluded
//   template <raster::RasterSource R>
//   SeamOutcome refine_seam(const R& dem, Point2 a, Point2 b, double tolerance);
//
//   - the check points are exactly constraint_check_points' for the edge (15f
//     D2): crossings, node crossings snapped to the node, and midpoints with
//     the ends as neighbours. On a grid-line seam that includes midpoints on
//     cell sides, which the greedy nearly never needs (N14: 2c + 1 check
//     points for c nodes, not the nodes alone);
//   - an inserted point is output at (x_min + col dx, y_max - row dy) of its
//     lattice position, with the check point's own z;
//   - a piece with an invalid end (NoData under a or b) has no lerp. The
//     check point on it nearest the invalid end is inserted (14's carving
//     rule, as 15f D4 applies it along an edge), until no such piece holds a
//     check point;
//   - std::invalid_argument, its message naming refine_seam, for a tolerance
//     that is negative or not finite, for a == b, and for an end outside the
//     node rectangle.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/seam.hpp>

#include "frozen_oracle.hpp"
#include "strip_oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using frozen_oracle::kSeam;
using frozen_oracle::Lat;
using frozen_oracle::Polyline;
using frozen_oracle::Side;
using strip_oracle::Mesh;
using strip_oracle::Terrain;
using terrain::Point2;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::refine;
using terrain::refinement::refine_seam;
using terrain::refinement::RefineOptions;
using terrain::refinement::RefineOutcome;
using terrain::refinement::SeamOutcome;

namespace {

bool before(Point2 p, Point2 q) { return p.x < q.x || (p.x == q.x && p.y < q.y); }

// The polyline the pass leaves: a, its points, b, in the lattice, parameter
// by projection onto a -> b (the oracle's own, not the output's s).
Polyline polyline(const RasterGeometry& g, const SeamOutcome& o) {
    const Lat a = strip_oracle::lat(g, o.a), b = strip_oracle::lat(g, o.b);
    Polyline line;
    line.t.push_back(0.0);
    line.z.push_back(o.z_a.value_or(0.0));
    line.valid.push_back(o.z_a.has_value());
    for (const auto& p : o.points) {
        line.t.push_back(strip_oracle::param_dist(a, b, strip_oracle::lat(g, p.at)).first);
        line.z.push_back(p.z);
        line.valid.push_back(true);
    }
    line.t.push_back(1.0);
    line.z.push_back(o.z_b.value_or(0.0));
    line.valid.push_back(o.z_b.has_value());
    return line;
}

// The oracle's check points for the one edge a -> b (15f's rulings, generated
// independently: strip_oracle::ruled_points).
std::vector<strip_oracle::OraclePoint> ruled(const Raster<float>& dem, Point2 a, Point2 b) {
    const std::vector<Point2> v{before(a, b) ? a : b, before(a, b) ? b : a};
    return strip_oracle::ruled_points(dem, v, {{0, 1}});
}

// What every outcome must satisfy, whatever the tolerance: the ends in order
// and unmoved, s strictly increasing inside (0, 1), each point on the edge
// (within the lattice slack), at one of the oracle's check points with its z,
// no more insertions than check points, and every check point within the
// tolerance of the polyline (or on a void piece, only if an end is invalid).
void check(const Raster<float>& dem, Point2 a, Point2 b, const SeamOutcome& o, double tol) {
    const RasterGeometry& g = dem.geometry();
    REQUIRE(o.a == (before(a, b) ? a : b));
    REQUIRE(o.b == (before(a, b) ? b : a));
    const Lat la = strip_oracle::lat(g, o.a), lb = strip_oracle::lat(g, o.b);
    const auto pts = ruled(dem, a, b);
    const double slack = strip_oracle::on_edge_slack(g);
    REQUIRE(o.check_points == pts.size());
    REQUIRE(o.points.size() <= o.check_points);
    double s_last = 0.0;
    for (const auto& p : o.points) {
        CAPTURE(p.at.x, p.at.y, p.s, p.z);
        REQUIRE(p.s > s_last);
        REQUIRE(p.s < 1.0);
        s_last = p.s;
        const Lat lp = strip_oracle::lat(g, p.at);
        REQUIRE(strip_oracle::param_dist(la, lb, lp).second <= slack);
        const auto it = std::find_if(pts.begin(), pts.end(), [&](const auto& q) {
            return std::hypot(q.at.col - lp.col, q.at.row - lp.row) <= slack;
        });
        REQUIRE(it != pts.end());
        REQUIRE(std::abs(it->z - p.z) <= 1e-9 * std::max(1.0, std::abs(p.z)));
    }
    const auto f = frozen_oracle::seam_findings(la, lb, polyline(g, o), pts, tol);
    CAPTURE(f.checked, f.over, f.on_void, f.worst, o.max_error);
    REQUIRE(f.over == 0);
    if (o.z_a && o.z_b) REQUIRE(f.on_void == 0);
    REQUIRE(o.max_error <= tol);
    REQUIRE(std::abs(o.max_error - f.worst) <= 1e-9 * std::max(1.0, f.worst));
}

// SP1's every-point property: on a grid-line seam the bilinear surface is
// linear between nodes, so the tolerance holds at every point of the seam,
// sampled at 1,000 points with the oracle's own bilinear height.
void every_point(const Raster<float>& dem, const SeamOutcome& o, double tol) {
    const RasterGeometry& g = dem.geometry();
    const Lat a = strip_oracle::lat(g, o.a), b = strip_oracle::lat(g, o.b);
    const Polyline line = polyline(g, o);
    std::size_t over = 0;
    for (int i = 0; i <= 1000; ++i) {
        const double s = i / 1000.0;
        const Lat p{a.col + s * (b.col - a.col), a.row + s * (b.row - a.row)};
        const auto z = strip_oracle::height(dem, p);
        const auto v = line.at(s);
        if (!z || !v) continue;
        if (std::abs(*z - *v) > tol + 1e-9 * std::max(1.0, std::abs(*z))) ++over;
    }
    REQUIRE(over == 0);
}

const RasterGeometry kG17 = strip_oracle::exact_geometry(17, 17);

Point2 at(const RasterGeometry& g, double col, double row) { return strip_oracle::world(g, col, row); }

Raster<float> profile(const std::vector<float>& z_along_row_0) {
    // A 2-row lattice whose row 0 carries the profile and row 1 copies it.
    const auto n = z_along_row_0.size();
    std::vector<float> v(z_along_row_0);
    v.insert(v.end(), z_along_row_0.begin(), z_along_row_0.end());
    return Raster<float>{strip_oracle::exact_geometry(n, 2), std::move(v)};
}

RefineOutcome refine_frozen(const Raster<float>& dem, const strip_oracle::Start& s, double tol) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = 1;
    o.frozen_mask = kSeam;
    return refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                  std::span<const std::uint32_t>{s.masks}, o);
}

}  // namespace

// ------------------------------------------------------------------------ SP1

TEST_CASE("SP1: on a grid-line seam the pass meets the tolerance at every point", "[refinement][seam][sp1]") {
    // Vertical (col 8) and horizontal (row 8) seams across the whole lattice,
    // each end on the outline; tester.md §3A, seams along grid lines.
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth, Terrain::Plane);
    const double tol = GENERATE(0.0, 0.5, 2.0, 1e9);
    const bool vertical = GENERATE(true, false);
    const std::uint32_t seed = GENERATE(1u, 2u);
    CAPTURE(static_cast<int>(kind), tol, vertical, seed);
    const auto dem = strip_oracle::terrain_dem(kG17, kind, seed);
    const Point2 a = vertical ? at(kG17, 8, 0) : at(kG17, 0, 8), b = vertical ? at(kG17, 8, 16) : at(kG17, 16, 8);
    const SeamOutcome o = refine_seam(dem, a, b, tol);
    check(dem, a, b, o, tol);
    every_point(dem, o, tol);
    // The 15 nodes strictly inside the seam are check points, with the 16
    // midpoints between them and the ends (15f's generator). An inserted node
    // is the DEM's node exactly, at its own value: "exact heights" on a
    // lattice-line seam.
    REQUIRE(o.check_points == 2 * 15 + 1);
    for (const auto& p : o.points) {
        const Lat l = strip_oracle::lat(kG17, p.at);
        if (!strip_oracle::is_node(l)) continue;
        const terrain::raster::CellIndex c{static_cast<std::size_t>(l.row), static_cast<std::size_t>(l.col)};
        REQUIRE(p.at == kG17.node(c));
        REQUIRE(p.z == static_cast<double>(dem.value_at(c)));
    }
    if (tol == 1e9 || (kind == Terrain::Plane && tol > 0.0)) REQUIRE(o.points.empty());
    if (tol == 0.0 && kind == Terrain::Rough) REQUIRE(!o.points.empty());
}

TEST_CASE("SP1: the ends' heights are vertex_z at them", "[refinement][seam][sp1]") {
    const auto dem = strip_oracle::terrain_dem(kG17, Terrain::Smooth, 1);
    const SeamOutcome o = refine_seam(dem, at(kG17, 8, 16), at(kG17, 8, 0), 0.5);
    // a < b by (x, y): world y = -row, so row 16 comes first.
    REQUIRE(o.a == at(kG17, 8, 16));
    REQUIRE(o.z_a == static_cast<double>(dem.value_at({16, 8})));
    REQUIRE(o.z_b == static_cast<double>(dem.value_at({0, 8})));
}

// ------------------------------------------------------------------------ SP2

TEST_CASE("SP2: on a general seam every check point is within tolerance", "[refinement][seam][sp2]") {
    // Rational (dyadic) ends, off-node, in every direction, short and long.
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth);
    const double tol = GENERATE(0.0, 0.25, 1.0, 5.0);
    const auto ends = GENERATE(std::array<double, 4>{1.25, 0.5, 14.75, 15.5},   // steep, down-right
                               std::array<double, 4>{15.5, 2.125, 0.375, 9.0},  // shallow, down-left
                               std::array<double, 4>{3.5, 3.5, 12.5, 12.5},     // a diagonal through nodes, ends off-node
                               std::array<double, 4>{2.0, 2.0, 14.0, 14.0},     // the diagonal through nodes
                               std::array<double, 4>{4.0, 1.0, 12.0, 15.0},     // through nodes every other row
                               std::array<double, 4>{7.3, 7.1, 8.9, 7.6});      // under two cells long
    const std::uint32_t seed = GENERATE(1u, 2u);
    CAPTURE(static_cast<int>(kind), tol, ends[0], ends[1], ends[2], ends[3], seed);
    const auto dem = strip_oracle::terrain_dem(kG17, kind, seed);
    const Point2 a = at(kG17, ends[0], ends[1]), b = at(kG17, ends[2], ends[3]);
    check(dem, a, b, refine_seam(dem, a, b, tol), tol);
}

TEST_CASE("SP2: ties go to the smallest parameter from a", "[refinement][seam][sp2]") {
    // Along row 0 from col 0 to col 4, z 0 10 0 10 0: nodes 1 and 3 are both
    // 10 off the chord, exactly. At tolerance 7 the first insertion settles the
    // other (6.67 off the new chord), so the tie decides the output: col 1,
    // the smaller parameter from a = (col 0), whichever end is given first.
    const auto dem = profile({0, 10, 0, 10, 0});
    const auto& g = dem.geometry();
    for (const bool reversed : {false, true}) {
        CAPTURE(reversed);
        const Point2 p0 = at(g, 0, 0), p4 = at(g, 4, 0);
        const SeamOutcome o = reversed ? refine_seam(dem, p4, p0, 7.0) : refine_seam(dem, p0, p4, 7.0);
        REQUIRE(o.a == p0);
        REQUIRE(o.points.size() == 1);
        REQUIRE(o.points[0].at == at(g, 1, 0));
        REQUIRE(o.points[0].z == 10.0);
        REQUIRE(o.points[0].s == 0.25);
        check(dem, p0, p4, o, 7.0);
    }
}

TEST_CASE("SP2: the greedy inserts the worst point first, not every point over", "[refinement][seam][sp2]") {
    // z 0 4 9 4 0 along row 0: node 2 is the worst (9); once it is in, nodes 1
    // and 3 are 0.5 off the new chords, within tolerance 1. A pass that
    // inserted every point over the tolerance would add all three.
    const auto dem = profile({0, 4, 9, 4, 0});
    const auto& g = dem.geometry();
    const SeamOutcome o = refine_seam(dem, at(g, 0, 0), at(g, 4, 0), 1.0);
    REQUIRE(o.points.size() == 1);
    REQUIRE(o.points[0].at == at(g, 2, 0));
    REQUIRE(o.max_error == 0.5);
}

// ------------------------------------------------------------------------ SP3

TEST_CASE("SP3: the edge given reversed gives the same output bit for bit", "[refinement][seam][sp3]") {
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth);
    const double tol = GENERATE(0.0, 0.5);
    const auto ends = GENERATE(std::array<double, 4>{8, 0, 8, 16}, std::array<double, 4>{1.25, 0.5, 14.75, 15.5},
                               std::array<double, 4>{15.5, 2.125, 0.375, 9.0});
    CAPTURE(static_cast<int>(kind), tol, ends[0], ends[1], ends[2], ends[3]);
    const auto dem = strip_oracle::terrain_dem(kG17, kind, 3);
    const Point2 a = at(kG17, ends[0], ends[1]), b = at(kG17, ends[2], ends[3]);
    const SeamOutcome f = refine_seam(dem, a, b, tol), r = refine_seam(dem, b, a, tol);
    REQUIRE(f.a == r.a);
    REQUIRE(f.b == r.b);
    REQUIRE(f.z_a == r.z_a);
    REQUIRE(f.z_b == r.z_b);
    REQUIRE(f.points.size() == r.points.size());
    for (std::size_t i = 0; i < f.points.size(); ++i) {
        REQUIRE(f.points[i].at == r.points[i].at);
        REQUIRE(f.points[i].z == r.points[i].z);
        REQUIRE(f.points[i].s == r.points[i].s);
    }
    REQUIRE(f.check_points == r.check_points);
    REQUIRE(f.max_error == r.max_error);
}

TEST_CASE("SP3: a strip window grown by whole cells gives the same output bit for bit", "[refinement][seam][sp3]") {
    // Two windows of one lattice: 17 x 17 from the origin, and one grown by 4
    // columns left, 3 rows up, 2 right and 5 down, every node's value a
    // function of its world position. Two neighbours that cut their strips
    // differently must still agree. The ends are chosen so that every
    // crossing is exact in both frames (a column span of 4 and a row span of
    // 8 cells), so this pins the window's irrelevance, not rounding luck;
    // general seams agree only to rounding (N17).
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth);
    const double tol = GENERATE(0.0, 0.5);
    const auto ends = GENERATE(std::array<double, 4>{8, 0, 8, 16}, std::array<double, 4>{0, 8, 16, 8},
                               std::array<double, 4>{1.5, 2.25, 5.5, 10.25});
    CAPTURE(static_cast<int>(kind), tol, ends[0], ends[1], ends[2], ends[3]);
    const auto small = strip_oracle::terrain_dem(kG17, kind, 5);
    const RasterGeometry big_g{-4 * 2.0, 3 * 1.0, 2.0, 1.0, 17 + 4 + 2, 17 + 3 + 5};
    std::vector<float> v(big_g.rows() * big_g.cols(), 7.0f);
    for (std::size_t r = 0; r < 17; ++r)
        for (std::size_t c = 0; c < 17; ++c) v[(r + 3) * big_g.cols() + (c + 4)] = small.value_at({r, c});
    const Raster<float> big{big_g, std::move(v)};
    REQUIRE(big_g.node({3, 4}) == kG17.node({0, 0}));

    const Point2 a = at(kG17, ends[0], ends[1]), b = at(kG17, ends[2], ends[3]);
    const SeamOutcome s = refine_seam(small, a, b, tol), l = refine_seam(big, a, b, tol);
    REQUIRE(s.z_a == l.z_a);
    REQUIRE(s.z_b == l.z_b);
    REQUIRE(s.points.size() == l.points.size());
    for (std::size_t i = 0; i < s.points.size(); ++i) {
        REQUIRE(s.points[i].at == l.points[i].at);
        REQUIRE(s.points[i].z == l.points[i].z);
    }
}

// ------------------------------------------------------------------------ SP4

TEST_CASE("SP4: NoData stencils are skipped and counted", "[refinement][seam][sp4]") {
    // A NoData node one column off the seam (row 8, col 9): every check point
    // whose cell touches it is dropped (vertex_z's rule), counted in no_data
    // as the oracle counts it, and never inserted.
    const auto ends = GENERATE(std::array<double, 4>{8, 0, 8, 16}, std::array<double, 4>{1.25, 0.5, 14.75, 15.5});
    const double tol = GENERATE(0.0, 0.5);
    CAPTURE(ends[0], ends[1], tol);
    const auto clean_dem = strip_oracle::terrain_dem(kG17, Terrain::Rough, 4);
    std::vector<float> v(kG17.rows() * kG17.cols());
    for (std::size_t r = 0; r < 17; ++r)
        for (std::size_t c = 0; c < 17; ++c) v[r * 17 + c] = clean_dem.value_at({r, c});
    v[8 * 17 + 9] = -9999.0f;
    const Raster<float> dem{kG17, std::move(v), -9999.0f};
    const Point2 a = at(kG17, ends[0], ends[1]), b = at(kG17, ends[2], ends[3]);
    const SeamOutcome o = refine_seam(dem, a, b, tol);
    const std::size_t all = ruled(clean_dem, a, b).size(), kept = ruled(dem, a, b).size();
    REQUIRE(kept < all);
    REQUIRE(o.no_data == all - kept);
    check(dem, a, b, o, tol);
    for (const auto& p : o.points) REQUIRE(strip_oracle::height(dem, strip_oracle::lat(kG17, p.at)).has_value());
}

TEST_CASE("SP4: an end on NoData carves along the seam from that end", "[refinement][seam][sp4][pinned]") {
    // Ruled by N10 (pinned here at the red step): node (row 0, col 8), the top
    // end, is NoData, so that end has no height and the piece touching it is
    // void. The check point nearest it on a void piece is inserted, again and
    // again, until no void piece holds a check point; nothing on a void piece
    // is measured. The midpoint at row 0.5 is dropped (its cell touches the
    // NoData node), the node at row 1 is kept (a node's height is its own
    // value), so at a tolerance nothing else exceeds, the pass inserts exactly
    // that node.
    const auto clean_dem = strip_oracle::terrain_dem(kG17, Terrain::Smooth, 6);
    std::vector<float> v(kG17.rows() * kG17.cols());
    for (std::size_t r = 0; r < 17; ++r)
        for (std::size_t c = 0; c < 17; ++c) v[r * 17 + c] = clean_dem.value_at({r, c});
    v[0 * 17 + 8] = -9999.0f;
    const Raster<float> dem{kG17, std::move(v), -9999.0f};
    const Point2 top = at(kG17, 8, 0), bottom = at(kG17, 8, 16);
    const SeamOutcome o = refine_seam(dem, top, bottom, 1e9);
    REQUIRE(o.a == bottom);  // a < b by (x, y): the bottom (y = -16) first
    REQUIRE(o.b == top);
    REQUIRE(o.z_a.has_value());
    REQUIRE_FALSE(o.z_b.has_value());
    const auto pts = ruled(dem, top, bottom);
    REQUIRE_FALSE(pts.empty());
    const Lat nearest = std::min_element(pts.begin(), pts.end(), [](const auto& x, const auto& y) {
                            return x.at.row < y.at.row;
                        })->at;
    REQUIRE(nearest == Lat{8, 1});
    REQUIRE(o.points.size() == 1);
    REQUIRE(strip_oracle::lat(kG17, o.points.back().at) == nearest);
    const auto f = frozen_oracle::seam_findings(strip_oracle::lat(kG17, o.a), strip_oracle::lat(kG17, o.b),
                                                polyline(kG17, o), pts, 1e9);
    REQUIRE(f.on_void == 0);  // the void piece holds no check point
}

// ------------------------------------------------------- degeneracy and scale

TEST_CASE("SP5: a seam inside one cell has its one midpoint; inside a NoData cell none", "[refinement][seam]") {
    // Degeneracy policy, "A seam edge shorter than a cell may have no check
    // point". 15f's generator, which the pass reuses, gives such an edge one
    // midpoint (the ends count as neighbours); the NoData cell drops it.
    const auto dem = strip_oracle::terrain_dem(kG17, Terrain::Rough, 7);
    const Point2 a = at(kG17, 3.25, 4.125), b = at(kG17, 3.75, 4.875);
    const SeamOutcome o = refine_seam(dem, a, b, 0.0);
    REQUIRE(o.check_points == 1);
    REQUIRE(o.no_data == 0);
    check(dem, a, b, o, 0.0);

    std::vector<float> v(17 * 17);
    for (std::size_t r = 0; r < 17; ++r)
        for (std::size_t c = 0; c < 17; ++c) v[r * 17 + c] = dem.value_at({r, c});
    v[4 * 17 + 4] = -9999.0f;  // a corner of the cell (row 4, col 3)
    const Raster<float> holed{kG17, std::move(v), -9999.0f};
    const SeamOutcome h = refine_seam(holed, a, b, 0.0);
    REQUIRE(h.check_points == 0);
    REQUIRE(h.no_data == 1);
    REQUIRE(h.points.empty());
}

TEST_CASE("SP5: a sub-millimetre seam at UTM scale", "[refinement][seam][scale]") {
    // tester.md §3A, extreme scales: 30 m cells at (5e5, 7e6), a seam 0.4 mm
    // long inside one cell.
    const RasterGeometry g{500000.0, 7000000.0, 30.0, 30.0, 17, 17};
    const auto dem = strip_oracle::terrain_dem(g, Terrain::Rough, 8);
    const Point2 a{500000.0 + 4.5 * 30.0, 7000000.0 - 6.25 * 30.0};
    const Point2 b{a.x + 0.0003, a.y - 0.0002};
    const SeamOutcome o = refine_seam(dem, a, b, 0.0);
    REQUIRE(o.check_points == 1);
    REQUIRE(o.points.size() <= 1);
    check(dem, a, b, o, 0.0);
}

TEST_CASE("SP5: a grid-line seam at the far edge of a 16,385-node lattice", "[refinement][seam][scale]") {
    // A lattice wider than 8,191 nodes at UTM coordinates: the seam runs
    // along row 2 from col 16,000 to the node rectangle's last column, and
    // along the last row (rows - 1, vertex_z's last-cell rule).
    const RasterGeometry g{500000.0, 7000000.0, 10.0, 5.0, 16385, 5};
    std::vector<float> v(g.rows() * g.cols());
    std::uint32_t x = 12345u;
    for (auto& z : v) {
        x = x * 1664525u + 1013904223u;
        z = 300.0f + static_cast<float>(x >> 20) / 64.0f;
    }
    const Raster<float> dem{g, std::move(v)};
    const double tol = GENERATE(0.0, 0.5, 4.0);
    const double row = GENERATE(2.0, 4.0);
    CAPTURE(tol, row);
    const Point2 a = at(g, 16000, row), b = at(g, 16384, row);
    const SeamOutcome o = refine_seam(dem, a, b, tol);
    check(dem, a, b, o, tol);
    every_point(dem, o, tol);
}

TEST_CASE("SP5: refusals name refine_seam", "[refinement][seam]") {
    const auto dem = strip_oracle::terrain_dem(kG17, Terrain::Smooth, 1);
    const Point2 a = at(kG17, 2, 2), b = at(kG17, 9, 5);
    for (const double bad : {-1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
        REQUIRE_THROWS_MATCHES(refine_seam(dem, a, b, bad), std::invalid_argument,
                               Catch::Matchers::MessageMatches(ContainsSubstring("refine_seam")
                                                               && ContainsSubstring("tolerance")));
    REQUIRE_THROWS_MATCHES(refine_seam(dem, a, a, 0.5), std::invalid_argument,
                           Catch::Matchers::MessageMatches(ContainsSubstring("refine_seam")));
    REQUIRE_THROWS_MATCHES(refine_seam(dem, a, at(kG17, 17.5, 3), 0.5), std::invalid_argument,
                           Catch::Matchers::MessageMatches(ContainsSubstring("refine_seam")
                                                           && ContainsSubstring("outside")));
}

// ------------------------------------------------------ the pass, then refine

TEST_CASE("SP6: a domain cut by a seam meets the tolerance at every node, seam included (K3, K4)",
          "[refinement][seam][k3]") {
    // The seam pass, then the start with its points in the seam (fans), then
    // refine with the seam frozen: as one run over the whole domain and as two
    // pieces run apart. Every valid node of the domain is within tolerance
    // (strip_oracle::node_findings, nothing excluded: the seam's nodes are the
    // seam pass's, and it met them); the two pieces write the same seam vertex
    // sequence, (x, y, z) bit for bit; the frozen seam is intact; Delaunay.
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth);
    const double tol = GENERATE(0.0, 0.5, 2.0);
    const auto seam = GENERATE(std::array<double, 2>{8, 8},     // grid line
                               std::array<double, 2>{4, 12},    // through nodes on every other row
                               std::array<double, 2>{5.25, 9.5});  // general
    CAPTURE(static_cast<int>(kind), tol, seam[0], seam[1]);
    const auto dem = strip_oracle::terrain_dem(kG17, kind, 9);
    const Point2 top = at(kG17, seam[0], 0), bottom = at(kG17, seam[1], 16);
    const SeamOutcome pass = refine_seam(dem, top, bottom, tol);
    check(dem, top, bottom, pass, tol);

    // The inner points from top to bottom, in the lattice.
    std::vector<Lat> inner;
    for (const auto& p : pass.points) inner.push_back(strip_oracle::lat(kG17, p.at));
    if (pass.a == bottom) std::reverse(inner.begin(), inner.end());

    std::vector<std::vector<std::pair<Point2, double>>> sequences;
    for (const auto side : {Side::Both, Side::Left, Side::Right}) {
        CAPTURE(static_cast<int>(side));
        const auto start = frozen_oracle::seam_domain(kG17, 16, 16, seam[0], seam[1], inner, side);
        const RefineOutcome out = refine_frozen(dem, start, tol);
        REQUIRE(out.ok());
        const Mesh m = strip_oracle::mesh_of(out);
        const auto nf = strip_oracle::node_findings(dem, m, tol);
        CAPTURE(nf.over, nf.worst);
        REQUIRE(nf.not_ccw == 0);
        REQUIRE(nf.over == 0);
        REQUIRE(strip_oracle::delaunay_violations(kG17, m) == 0);
        REQUIRE(frozen_oracle::clean(
            frozen_oracle::frozen_findings(kG17, start.mesh.vertices(), start.edges, start.masks, kSeam, m)));
        if (side != Side::Both)
            sequences.push_back(frozen_oracle::seam_sequence(kG17, m, Lat{seam[0], 0}, Lat{seam[1], 16}));
    }
    REQUIRE(sequences.size() == 2);
    REQUIRE(sequences[0].size() == pass.points.size() + 2);
    REQUIRE(sequences[0].size() == sequences[1].size());
    for (std::size_t i = 0; i < sequences[0].size(); ++i) {
        REQUIRE(sequences[0][i].first == sequences[1][i].first);
        REQUIRE(sequences[0][i].second == sequences[1][i].second);
    }
}

TEST_CASE("SP2: an error equal to the tolerance is within it", "[refinement][seam][sp2]") {
    // "While some check point p has |z_p - lerp(p)| > tolerance": strictly
    // greater. Along row 0, a plane at tolerance 0 (every error exactly 0), and
    // z 0 2 0 at tolerance 2 (node 1 exactly 2 off the chord, the midpoints 1),
    // both need no point.
    const auto plane = profile({0, 1, 2, 3, 4});
    const auto& gp = plane.geometry();
    const SeamOutcome p = refine_seam(plane, at(gp, 0, 0), at(gp, 4, 0), 0.0);
    REQUIRE(p.check_points == 7);
    REQUIRE(p.points.empty());
    REQUIRE(p.max_error == 0.0);

    const auto tent = profile({0, 2, 0});
    const auto& gt = tent.geometry();
    const SeamOutcome t = refine_seam(tent, at(gt, 0, 0), at(gt, 2, 0), 2.0);
    REQUIRE(t.check_points == 3);
    REQUIRE(t.points.empty());
    REQUIRE(t.max_error == 2.0);
}
