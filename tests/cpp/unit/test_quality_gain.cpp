// Increment 20c, PR 20c-2 (docs/increments/20c-soft-quality.md, R7, R8 and
// "Tests @tester writes red first", 20c-2): the quality start inserts a
// candidate only when it does not lower the worst angle of the triangles it
// replaces (R7), and splits the constraint line a blocked walk stops at (R8).
// The design's invariant-critical suite test_quality_gain: T-P1, T-P2 and LS1.
// T-P3 and determinism are in tests/cpp/property/prop_quality_gain_refine.cpp.
//
// Interface, as R7 and R8 name it:
//
//   QualityOptions::min_gain_deg   P, degrees; negative is off (the hard rule)
//   QualityOutcome::skipped_no_gain  candidates refused by R7
//   QualityOutcome::line_splits      R8's splits
//
// PINNED HERE, where the design leaves it open (listed for @architect):
//   - QualityOptions::min_gain_deg defaults to a negative value (off), so
//     every caller that does not set it, the existing suites included, keeps
//     20c-1's behaviour; set here by member assignment;
//   - R7's read-only cavity is a function of its own, so T-P1 can compare the
//     prediction with what legalise_around writes:
//
//       namespace terrain::mesh::detail {
//       struct QualityCavity {
//           std::vector<std::uint32_t> removed;                // the cavity's slots
//           std::vector<std::array<MeshVertex, 3>> created;    // p joined to each
//       };                                                     // boundary edge, CCW
//       template <pred::GeometryKernel K>
//       QualityCavity quality_cavity(const LatticeMesh& m, std::uint32_t t, unsigned on,
//                                    MeshVertex p, const LatticeFrame& f);
//       }
//
//     `on` 0..2: p goes in on t's edge `on` (split_edge; the triangle across
//     it, if any, seeds the cavity too, constrained or not); 3: strictly
//     inside t (split_inside). improve asks it for every judged candidate;
//     the order of `removed` and `created` is free, each created triangle is
//     counter-clockwise on (col, -row) like a mesh triangle;
//   - the angles are in the frame (col dx, -(row dy)) improve is given, each
//     capped at theta, and "by more than rounding" (T-P2) is R7's slack,
//     1e-6 degrees: R7 accepts new >= old + P - s (ruled by @architect,
//     docs/increments/20c-soft-quality.md, "Pins ruled for 20c-2's red
//     step", pin 3);
//   - R8 runs only when R7 is on and constraint_feet is on (same ruling,
//     pin 4); every LS1 case that expects no split turns the feet on, so
//     that what refuses the split is the rule under test, not the switch;
//   - a split by R8 is not a skip; a blocked walk R8 declines is counted in
//     skipped_blocked, once, as today;
//   - the validity callable is asked at least once between two insertions
//     (it is today, about every node, and R2 and R8 ask it about the foot),
//     which the recorder below relies on.
//
// Oracles are this file's and tests/cpp/support/quality_gain_support.hpp's: the
// triangles an insertion wrote are the difference of the meshes before and
// after it, angles are recomputed from the corners, and topology, Delaunay
// and the constraint lines are checked from the mesh alone.
//
// Planted faults each case targets (the design's list):
//   T-P1  the cavity crossing a constrained edge        "stops at a constrained edge"
//         an inexact incircle (a plain double            "near-cocircular"
//         determinant instead of the kernel's)
//         the foot's second seed missing                 "on a constrained edge, both sides"
//   T-P2  the acceptance test reversed                   "the only candidate ... refused",
//                                                        "every accepted candidate"
//         the slack dropped and the test made strict     "a candidate that keeps the worst
//                                                        angle exactly is accepted",
//                                                        "a foot that keeps its edge's end
//                                                        angle is accepted"
//   LS1   the split at the midpoint                      "splits the line at the node's foot"
//         the end check removed                          "within half a cell of an end"
//         R8 left on at gain -1                          "gain -1 turns the split off"
//         R8 left on with the feet off                   "with the feet off, no split"
// Added by the mutation round, for faults the list above did not reach:
//   T-P1  quad_flips' integer path flipping on           "on the integer frame a candidate
//         Cocircular (shared with Lawson)                on the circle ..."
//         quad_flips' frame-collinear branch dropped     "a side collinear in the frame ..."
//         (shared with Lawson)
//         the slot across a split edge not seeded        "a line split flips beyond ..."
//         by improve's insert routine
// Each case's comment names the fault it is for.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/mesh/quality.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/kernel.hpp>

#include "constraint_foot_oracles.hpp"
#include "quality_fixtures.hpp"
#include "quality_gain_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <numbers>
#include <optional>
#include <set>
#include <utility>
#include <vector>

using terrain::Point2;
using terrain::TriangleIndices;
using terrain::mesh::improve;
using terrain::mesh::kNoNeighbour;
using terrain::mesh::LatticeFrame;
using terrain::mesh::LatticeMesh;
using terrain::mesh::legalise_all;
using terrain::mesh::MeshVertex;
using terrain::mesh::orient_sign;
using terrain::mesh::QualityOptions;
using terrain::mesh::QualityOutcome;
using terrain::pred::DefaultKernel;
using terrain::pred::FastKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;

namespace qg = quality_gain;
namespace cfo = constraint_foot_oracles;

namespace {

constexpr double kTheta = 25.0;
// Degrees; T-P2's "by more than rounding", R7's slack s (docs/increments/20c-soft-quality.md,
// R7, "Why the slack"). The test keeps its own literal, not quality.hpp's constant.
// Scale: angles under theta = 25 degrees on grids up to 33 cells a side.
constexpr double kRounding = 1e-6;
constexpr double kOnLine = 1e-9;    // cells, for grids up to 33 nodes a side

QualityOptions options(const qg::Fixture& fx, double gain, bool feet = false) {
    QualityOptions o{kTheta, fx.rows, fx.cols};
    o.constraint_feet = feet;
    o.min_gain_deg = gain;
    return o;
}

// The fixture's mesh, legalised in f as refine does before the pass. A hand
// fixture is a CDT already (no flip); the Q7 fan is not.
LatticeMesh built(const qg::Fixture& fx, const LatticeFrame& f, bool cdt = true) {
    auto m = qg::build(fx);
    REQUIRE(m.has_value());
    const auto flips = legalise_all<DefaultKernel>(*m, f, [](std::uint32_t) {});
    if (cdt) REQUIRE(flips == 0);
    return std::move(*m);
}

struct Run {
    QualityOutcome out;
    std::vector<qg::Step> steps;
};

// improve with the recorder in its validity callable; `valid` decides as given.
template <class Valid>
Run recorded(LatticeMesh& m, const LatticeFrame& f, const QualityOptions& o, Valid valid) {
    qg::Recorder rec{m};
    Run r;
    r.out = improve<DefaultKernel>(m, f, o, [&](const MeshVertex& v) {
        rec.snap();
        return valid(v);
    });
    rec.snap();
    bool ok = false;
    r.steps = rec.steps(&ok);
    REQUIRE(ok);  // one vertex per insertion, each located
    return r;
}

Run recorded(LatticeMesh& m, const LatticeFrame& f, const QualityOptions& o) {
    return recorded(m, f, o, [](const MeshVertex&) { return true; });
}

// T-P1 for one insertion: R7's prediction on the mesh before it equals what
// the insertion wrote, as sets of corner triples.
void check_prediction(const qg::Step& s, const LatticeFrame& f) {
    const auto c = terrain::mesh::detail::quality_cavity<DefaultKernel>(s.before, s.at.t, s.at.on, s.p, f);
    std::set<qg::Tri> removed, created;
    for (const auto x : c.removed) removed.insert(qg::tri_of(s.before, x));
    for (const auto& t : c.created) created.insert(qg::tri(t[0], t[1], t[2]));
    REQUIRE(removed.size() == c.removed.size());  // no slot twice
    REQUIRE(created.size() == c.created.size());
    const auto d = qg::diff(s.before, s.after);
    CAPTURE(s.p.col, s.p.row, s.at.t, s.at.on, d.removed.size(), d.created.size(), removed.size(), created.size());
    CHECK(removed == d.removed);
    CHECK(created == d.created);
}

// T-P2 for one judged insertion: the worst angle, capped at theta, did not
// fall by more than rounding below old + gain.
void check_gain(const qg::Step& s, const LatticeFrame& f, double gain) {
    const auto d = qg::diff(s.before, s.after);
    const double old_w = qg::worst_capped(d.removed, f, kTheta), new_w = qg::worst_capped(d.created, f, kTheta);
    CAPTURE(s.p.col, s.p.row, old_w, new_w, gain);
    CHECK(new_w >= old_w + gain - kRounding);
}

// ------------------------------------------------------------------ mesh checks

using Pair = std::pair<std::uint32_t, std::uint32_t>;

void check_topology(const LatticeMesh& m) {
    std::map<Pair, std::uint32_t> directed;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        CAPTURE(t);
        REQUIRE(orient_sign(m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)) > 0);
        for (unsigned k = 0; k < 3; ++k)
            REQUIRE(directed.emplace(Pair{m.triangles()[t][k], m.triangles()[t][(k + 1) % 3]}, t).second);
    }
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto& tri = m.triangles()[t];
            const auto it = directed.find({tri[(k + 1) % 3], tri[k]});
            CAPTURE(t, k);
            REQUIRE(m.neighbours(t)[k] == (it == directed.end() ? kNoNeighbour : it->second));
            if (it == directed.end()) continue;
            unsigned j = 0;
            while (m.triangles()[it->second][j] != tri[(k + 1) % 3]) ++j;
            REQUIRE(m.is_constrained(t, k) == m.is_constrained(it->second, j));
            REQUIRE(m.mask(t, k) == m.mask(it->second, j));
        }
}

std::size_t delaunay_violations(const LatticeMesh& m, const LatticeFrame& f) {
    std::size_t bad = 0;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto u = m.neighbours(t)[k];
            if (u == kNoNeighbour || m.is_constrained(t, k)) continue;
            unsigned j = 0;
            while (m.neighbours(u)[j] != t) ++j;
            if (DefaultKernel::incircle(f.at(m.corner(t, 0)), f.at(m.corner(t, 1)), f.at(m.corner(t, 2)),
                                        f.at(m.corner(u, (j + 2) % 3)))
                == Incircle::Inside)
                ++bad;
        }
    return bad;
}

// 20 Q3's check, to kOnLine: every constraint edge on exactly one input
// segment with its mask, and each segment's pieces chaining end to end.
void check_lines(const LatticeMesh& m, const qg::Fixture& fx) {
    const auto v = m.vertices();
    const auto [edges, masks] = m.constraint_edges();
    std::map<Pair, std::vector<std::pair<double, double>>> pieces;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        CAPTURE(i);
        std::size_t owners = 0;
        for (const auto& [seg, mask] : fx.constraints) {
            const MeshVertex a = fx.vertices[seg.first], b = fx.vertices[seg.second];
            const cfo::Frac fa{a.col, a.row}, fb{b.col, b.row};
            const auto [tp, dp] = cfo::param_dist(fa, fb, cfo::Frac{v[edges[i][0]].col, v[edges[i][0]].row});
            const auto [tq, dq] = cfo::param_dist(fa, fb, cfo::Frac{v[edges[i][1]].col, v[edges[i][1]].row});
            if (dp > kOnLine || dq > kOnLine || std::min(tp, tq) < -kOnLine || std::max(tp, tq) > 1 + kOnLine) continue;
            ++owners;
            REQUIRE(masks[i] == mask);  // bits and masks on both halves
            pieces[seg].push_back(std::minmax(tp, tq));
        }
        REQUIRE(owners == 1);
    }
    REQUIRE(pieces.size() == fx.constraints.size());
    for (auto& [seg, ps] : pieces) {
        std::sort(ps.begin(), ps.end());
        REQUIRE(std::abs(ps.front().first) <= kOnLine);
        REQUIRE(std::abs(ps.back().second - 1.0) <= kOnLine);
        for (std::size_t j = 1; j < ps.size(); ++j) REQUIRE(std::abs(ps[j].first - ps[j - 1].second) <= kOnLine);
    }
}

void check_mesh(const LatticeMesh& m, const qg::Fixture& fx, const LatticeFrame& f) {
    check_topology(m);
    REQUIRE(delaunay_violations(m, f) == 0);
    check_lines(m, fx);
}

bool has_vertex(const LatticeMesh& m, MeshVertex p, double eps) {
    return std::any_of(m.vertices().begin(), m.vertices().end(), [&](MeshVertex v) {
        return std::abs(v.col - p.col) <= eps && std::abs(v.row - p.row) <= eps;
    });
}

// The constraint edges of the mesh that carry `mask`.
std::size_t edges_with_mask(const LatticeMesh& m, std::uint32_t mask) {
    const auto masks = m.constraint_edges().second;
    return static_cast<std::size_t>(std::count(masks.begin(), masks.end(), mask));
}

// The circumcentre of (a, b, c) in (col, row), frame (1, 1).
MeshVertex centre(MeshVertex a, MeshVertex b, MeshVertex c) {
    const double bx = b.col - a.col, by = b.row - a.row, cx = c.col - a.col, cy = c.row - a.row;
    const double d = 2.0 * (bx * cy - by * cx), b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    return MeshVertex{a.col + (cy * b2 - by * c2) / d, a.row + (bx * c2 - cx * b2) / d};
}

// The angle at o between the rays to u and v, degrees, in frame f.
double angle_at(MeshVertex o, MeshVertex u, MeshVertex v, const LatticeFrame& f) {
    const Point2 po = f.at(o), pu = f.at(u), pv = f.at(v);
    const double ux = pu.x - po.x, uy = pu.y - po.y, vx = pv.x - po.x, vy = pv.y - po.y;
    return std::atan2(std::abs(ux * vy - uy * vx), ux * vx + uy * vy) * 180.0 / std::numbers::pi;
}

// ------------------------------------------------------------------ fixtures

// 20c-1's own-edge fixture (tests/cpp/unit/test_constraint_foot_quality.cpp): A (0, 0.99),
// C (2, 6.99), B (20, 0.99), one triangle, every edge an outline. Its node
// N = (10, 1) lies 0.01 cells inside A-B, so the node's three triangles hold
// a sliver on A-B: a candidate that lowers the worst angle.
qg::Fixture own_edge() {
    qg::Fixture f;
    f.vertices = {{0, 0.99}, {2, 6.99}, {20, 0.99}};
    f.triangles = {{0, 1, 2}};
    f.constraints = {{{0u, 2u}, 1u}, {{0u, 1u}, 2u}, {{1u, 2u}, 4u}};
    f.rows = 8;
    f.cols = 21;
    return f;
}

// A Delaunay mesh of nine nodes in a 12 x 12 square (found by a seeded
// search, the arrays as it printed them; no constraints, the square's sides
// have no neighbour). The first candidate improve inserts is slot 7's node
// (7, 3); the worst angle of its cavity, atan(1/3) = 18.43 degrees, is also
// the worst angle of its new triangles, in a triangle that is the old one
// moved by a lattice vector, so the two are equal to the bit under any angle
// formula built from edge vectors. Gain 0 accepts it; a strict test refuses.
qg::Fixture equal_gain() {
    qg::Fixture f;
    f.vertices = {{0, 0}, {12, 0}, {12, 12}, {0, 12}, {6, 1}, {7, 1}, {8, 3}, {9, 5}, {3, 5}};
    f.triangles = {{7, 3, 2}, {4, 0, 8}, {2, 1, 7}, {1, 0, 4}, {1, 4, 5}, {8, 6, 4},
                   {1, 5, 6}, {6, 5, 4}, {1, 6, 7}, {8, 7, 6}, {0, 3, 8}, {8, 3, 7}};
    f.rows = f.cols = 13;
    return f;
}

// Two triangles across a constrained edge a-b: Y = (a, d, b) below it, and
// X = (a, b, c) above, flat, its circle reaching far below a-b. p = (4, 5) in
// Y lies inside X's circle, but the constraint is between them.
qg::Fixture across_a_line() {
    qg::Fixture f;
    f.vertices = {{0, 2}, {8, 2}, {4, 1}, {4, 8}};  // a, b, c, d
    f.triangles = {{0, 3, 1}, {0, 1, 2}};          // Y, X
    f.constraints = {{{0u, 1u}, 2u}, {{0u, 3u}, 1u}, {{1u, 3u}, 1u}, {{0u, 2u}, 1u}, {{1u, 2u}, 1u}};
    f.rows = f.cols = 9;
    return f;
}

// The rectangle a (1.5, 0), b (1.9, 0), c (1.9, 0.3), p (1.5, 0.3), each
// coordinate k * 0.1 as a double, so no coordinate is exact: X = (a, c, b)
// and the point p are cocircular on the integers (15, 0), (19, 0), (19, 3),
// (15, 3), and in these doubles the exact kernel says Cocircular while the
// plain double determinant says Inside, in every rotation of X (checked
// below). Y = (a, e, c) with e (1.3, 0.6) holds p. Frame (1, 1), so the
// frame point is the vertex itself.
qg::Fixture near_cocircular() {
    qg::Fixture f;
    f.vertices = {{15 * 0.1, 0 * 0.1}, {19 * 0.1, 0 * 0.1}, {19 * 0.1, 3 * 0.1}, {13 * 0.1, 6 * 0.1}};  // a b c e
    f.triangles = {{0, 2, 1}, {0, 3, 2}};  // X, Y
    f.constraints = {{{0u, 1u}, 1u}, {{1u, 2u}, 1u}, {{2u, 3u}, 1u}, {{0u, 3u}, 1u}};
    f.rows = f.cols = 3;
    return f;
}
const MeshVertex kCocircularP{15 * 0.1, 3 * 0.1};

// LS1: qg::line_beyond(), a long line e from E1 (2, 10) to E2 (30, 10.5)
// whose worst triangle's node (15, 18) lies beyond it (see the support header).
using qg::line_beyond;
constexpr std::uint32_t kLineMask = 2;
const MeshVertex kBeyondNode{15, 18};

// The node's foot on E1-E2, the orthogonal projection, frame (1, 1).
MeshVertex foot_on_e(const qg::Fixture& fx, MeshVertex n) {
    const MeshVertex a = fx.vertices[0], b = fx.vertices[1];
    const double ux = b.col - a.col, uy = b.row - a.row;
    const double s = ((n.col - a.col) * ux + (n.row - a.row) * uy) / (ux * ux + uy * uy);
    return MeshVertex{a.col + s * ux, a.row + s * uy};
}

// The same, with e an outline: the triangle below it, (E2, E1, D), removed.
qg::Fixture line_beyond_outline() {
    qg::Fixture f = line_beyond();
    f.triangles.erase(f.triangles.begin() + 1);
    f.constraints.erase({1u, 3u});
    f.constraints.erase({0u, 3u});
    f.vertices.erase(f.vertices.begin() + 3);  // D: renumber 4, 5, 6 to 3, 4, 5
    for (auto& t : f.triangles)
        for (auto& i : t)
            if (i > 3) --i;
    return f;
}

// LS1's end case: a short line E1 (0.75, 10) to E2 (2.25, 10), mask 2, with A
// (1.5, 9.875) just above it and T (1.5, 8.5), D (1.5, 11.5) closing the
// outline. B = (E1, E2, A), angles 9.5 degrees at E1 and E2, is the only bad
// triangle; its circumcentre (1.5, 12.1875) rounds to the node (2, 12), whose
// foot on e is (2, 10), 0.25 cells from E2: within half a cell of an end.
qg::Fixture line_near_end() {
    qg::Fixture f;
    f.vertices = {{0.75, 10}, {2.25, 10}, {1.5, 8.5}, {1.5, 11.5}, {1.5, 9.875}};
    f.triangles = {{0, 1, 4}, {1, 0, 3}, {1, 2, 4}, {2, 0, 4}};
    f.constraints = {{{0u, 1u}, 2u}, {{0u, 3u}, 1u}, {{1u, 3u}, 1u}, {{1u, 2u}, 1u}, {{0u, 2u}, 1u}};
    f.rows = 25;
    f.cols = 33;
    return f;
}

// A candidate on the circle of the triangle across, every vertex a node, for
// quad_flips' integer path. In (x, y) on the circle x^2 + y^2 = 25: U = (5, 0)
// (3, 4) (-5, 0) has the circle; t = (3, 4) (5, 0) (8, 6) holds p = (4, 3)
// strictly, and E = (8, 6) is outside U's circle (distance 10). Mapped to
// (col, row) = (x + 5, 10 - y), which keeps the orientation.
qg::Fixture cocircular_nodes() {
    qg::Fixture f;
    f.vertices = {{10, 10}, {8, 6}, {0, 10}, {13, 4}};  // P1 P2 P3 E
    f.triangles = {{0, 1, 2}, {1, 0, 3}};                // U, t
    f.rows = 11;
    f.cols = 14;
    return f;
}
const MeshVertex kOnUsCircle{9, 7};

// tests/cpp/unit/test_mesh_lawson.cpp's L4b quad, as a candidate's cavity: A (0, 0),
// B (1, 1), D (0, 1), E (1, 0) in (col, row), t = (A, B, E) and u = (B, A, D).
// C = (0.34, the double below 0.34) lies inside t, one ulp of row off A-B, so
// (A, B, C) is counter-clockwise in the mesh but collinear in the frame
// (3, 1.5); seen from u, C is strictly inside (B, A, D)'s circle.
qg::Fixture frame_collinear() {
    qg::Fixture f;
    f.vertices = {{0, 0}, {1, 1}, {0, 1}, {1, 0}};  // A B D E
    f.triangles = {{0, 1, 3}, {1, 0, 2}};           // t, u
    f.rows = f.cols = 2;
    return f;
}
const MeshVertex kFrameCollinearC{0.34, std::nextafter(0.34, 0.0)};

// LS1's line with a triangle (D, G, E2) beyond the far side's edge D-E2 (now
// interior and free), G = (25, 19). Its circle holds the foot of the node
// (15, 18) on e, about 4.9 cells from the centre (17.25, 9.25) against a
// radius of about 12.8, and not E1: so the split of e must flip D-E2 in the
// slot across e, which only the far side's seed reaches.
qg::Fixture line_beyond_far_flip() {
    qg::Fixture f = line_beyond();
    f.vertices.push_back({25, 19});  // G, index 7
    f.triangles.push_back({3, 7, 1});  // (D, G, E2)
    f.constraints.erase({1u, 3u});
    f.constraints.emplace(std::pair{1u, 7u}, 1u);
    f.constraints.emplace(std::pair{3u, 7u}, 1u);
    return f;
}

// 20's Q7 ring: a 23-gon fanned from one vertex on a 33 x 33 grid, every
// side an outline with mask 1, so R8 can never fire (nothing lies beyond).
qg::Fixture q7_ring() {
    qg::Fixture f;
    const auto ring = quality_fixtures::q7_ring();
    for (const auto& [c, r] : ring) f.vertices.push_back(MeshVertex{c, r});
    const auto n = static_cast<std::uint32_t>(ring.size());
    for (std::uint32_t i = 1; i + 1 < n; ++i) f.triangles.push_back({0, i, i + 1});
    for (std::uint32_t i = 0; i < n; ++i) f.constraints.emplace(std::minmax(i, (i + 1) % n), 1u);
    f.rows = f.cols = 33;
    return f;
}

}  // namespace

// ================================================================== the option

TEST_CASE("R7: QualityOptions leaves the gain test off by default", "[mesh][quality][gain][R7]") {
    REQUIRE(QualityOptions{}.min_gain_deg < 0.0);  // pinned: off unless set
}

// ================================================================== T-P1, the helper on hand fixtures

TEST_CASE("T-P1: the cavity stops at a constrained edge, though the circle beyond holds p",
          "[mesh][quality][gain][T-P1]") {
    // Planted fault: the cavity crossing a constrained edge.
    const LatticeFrame f{1.0, 1.0};
    const auto fx = across_a_line();
    const LatticeMesh m = built(fx, f);
    const MeshVertex p{4, 5};
    {  // the premise: p strictly inside Y, strictly inside X's circle, and a-b constrained
        const auto at = qg::locate(m, p);
        REQUIRE(at.has_value());
        REQUIRE(at->t == 0);
        REQUIRE(at->on == 3);
        REQUIRE(DefaultKernel::incircle(f.at(m.corner(1, 0)), f.at(m.corner(1, 1)), f.at(m.corner(1, 2)), f.at(p))
                == Incircle::Inside);
        REQUIRE(m.is_constrained(0, 2));  // Y = (a, d, b): edge 2 is b-a
    }
    LatticeMesh after = m;
    qg::replay<DefaultKernel>(after, {0, 3}, p, f);
    const auto c = terrain::mesh::detail::quality_cavity<DefaultKernel>(m, 0, 3, p, f);
    REQUIRE(c.removed == std::vector<std::uint32_t>{0});
    REQUIRE(c.created.size() == 3);
    check_prediction(qg::Step{m, after, p, {0, 3}}, f);
}

TEST_CASE("T-P1: on a near-cocircular quad the cavity follows the exact incircle, not the double one",
          "[mesh][quality][gain][T-P1]") {
    // Planted fault: an inexact incircle, a plain double determinant instead
    // of the kernel's. The premise checks that the doubles disagree here.
    const LatticeFrame f{1.0, 1.0};
    const auto fx = near_cocircular();
    const LatticeMesh m = built(fx, f);
    const MeshVertex p = kCocircularP;
    {
        const auto at = qg::locate(m, p);
        REQUIRE(at.has_value());
        REQUIRE(at->t == 1);  // inside Y
        REQUIRE(at->on == 3);
        REQUIRE_FALSE(m.is_constrained(1, 2));  // Y's edge c-a is X's: the cavity may grow there
        const std::array<Point2, 3> x{f.at(m.corner(0, 0)), f.at(m.corner(0, 1)), f.at(m.corner(0, 2))};
        REQUIRE(DefaultKernel::incircle(x[0], x[1], x[2], f.at(p)) == Incircle::Cocircular);
        for (std::size_t k = 0; k < 3; ++k) {
            CAPTURE(k);
            REQUIRE(FastKernel::incircle(x[k], x[(k + 1) % 3], x[(k + 2) % 3], f.at(p)) == Incircle::Inside);
        }
        // And in the form legalise_around asks it, p's new triangle (a, p, c)
        // against X's apex b: exact says no flip, the double says flip.
        const Point2 a = f.at(m.corner(0, 0)), c = f.at(m.corner(0, 1)), b = f.at(m.corner(0, 2));
        const std::array<Point2, 3> y{a, f.at(p), c};
        for (std::size_t k = 0; k < 3; ++k) {
            CAPTURE(k);
            REQUIRE(DefaultKernel::incircle(y[k], y[(k + 1) % 3], y[(k + 2) % 3], b) != Incircle::Inside);
            REQUIRE(FastKernel::incircle(y[k], y[(k + 1) % 3], y[(k + 2) % 3], b) == Incircle::Inside);
        }
    }
    LatticeMesh after = m;
    qg::replay<DefaultKernel>(after, {1, 3}, p, f);
    REQUIRE(qg::diff(m, after).removed.size() == 1);  // Lawson left X alone
    const auto c = terrain::mesh::detail::quality_cavity<DefaultKernel>(m, 1, 3, p, f);
    REQUIRE(c.removed == std::vector<std::uint32_t>{1});
    check_prediction(qg::Step{m, after, p, {1, 3}}, f);
}

TEST_CASE("T-P1: a point on a constrained edge takes the triangles on both sides",
          "[mesh][quality][gain][T-P1]") {
    // Planted fault: the foot's second seed missing (the triangle across the
    // split edge, which the growth cannot reach across the constraint).
    const LatticeFrame f{1.0, 1.0};
    const auto fx = across_a_line();
    const LatticeMesh m = built(fx, f);
    const MeshVertex p{3, 2};  // on a-b, off its midpoint
    REQUIRE(orient_sign(fx.vertices[0], fx.vertices[1], p) == 0);
    // From either side: Y's edge 2 is (b, a), X's edge 0 is (a, b).
    const auto [t, on] = GENERATE(std::pair<std::uint32_t, unsigned>{0, 2}, std::pair<std::uint32_t, unsigned>{1, 0});
    CAPTURE(t, on);
    REQUIRE(m.is_constrained(t, on));
    LatticeMesh after = m;
    qg::replay<DefaultKernel>(after, {t, on}, p, f);
    REQUIRE(qg::diff(m, after).removed.size() == 2);
    const auto c = terrain::mesh::detail::quality_cavity<DefaultKernel>(m, t, on, p, f);
    REQUIRE(c.removed.size() == 2);
    REQUIRE(c.created.size() == 4);
    check_prediction(qg::Step{m, after, p, {t, on}}, f);
}

TEST_CASE("T-P1: a point on an outline edge makes two triangles, not a flat third",
          "[mesh][quality][gain][T-P1]") {
    const LatticeFrame f{1.0, 1.0};
    const auto fx = across_a_line();
    const LatticeMesh m = built(fx, f);
    const MeshVertex p{2, 5};  // the midpoint of a-d, Y's edge 0, an outline
    REQUIRE(orient_sign(fx.vertices[0], fx.vertices[3], p) == 0);
    REQUIRE(m.neighbours(0)[0] == kNoNeighbour);
    LatticeMesh after = m;
    qg::replay<DefaultKernel>(after, {0, 0}, p, f);
    const auto c = terrain::mesh::detail::quality_cavity<DefaultKernel>(m, 0, 0, p, f);
    REQUIRE(c.created.size() == 2);
    check_prediction(qg::Step{m, after, p, {0, 0}}, f);
}

// The next two cases name the cavity they expect, because the cavity and
// legalise_around share detail::quad_flips: a fault there moves both, and
// check_prediction alone cannot see it.

TEST_CASE("T-P1: on the integer frame a candidate on the circle of the triangle across does not grow the cavity",
          "[mesh][quality][gain][T-P1]") {
    // Planted fault: quad_flips' integer path flipping on Cocircular.
    const auto fx = cocircular_nodes();
    const LatticeFrame f = terrain::mesh::lattice_frame(1.0, 1.0, fx.rows, fx.cols);
    REQUIRE(f.integer());
    const LatticeMesh m = built(fx, f);
    const MeshVertex p = kOnUsCircle;
    {  // the premise: p strictly inside t, on U's circle, and the integer path answers it
        const auto at = qg::locate(m, p);
        REQUIRE(at.has_value());
        REQUIRE(at->t == 1);
        REQUIRE(at->on == 3);
        REQUIRE(DefaultKernel::incircle(f.at(m.corner(0, 0)), f.at(m.corner(0, 1)), f.at(m.corner(0, 2)), f.at(p))
                == Incircle::Cocircular);
        REQUIRE(terrain::mesh::lattice_incircle(m.corner(0, 0), m.corner(0, 1), m.corner(0, 2), p, f)
                == std::optional<Incircle>{Incircle::Cocircular});
    }
    LatticeMesh after = m;
    qg::replay<DefaultKernel>(after, {1, 3}, p, f);
    REQUIRE(qg::diff(m, after).removed.size() == 1);  // Lawson leaves U alone
    const auto c = terrain::mesh::detail::quality_cavity<DefaultKernel>(m, 1, 3, p, f);
    REQUIRE(c.removed == std::vector<std::uint32_t>{1});
    REQUIRE(c.created.size() == 3);
    check_prediction(qg::Step{m, after, p, {1, 3}}, f);
}

TEST_CASE("T-P1: a side collinear in the frame is decided from the triangle across, and the cavity grows there",
          "[mesh][quality][gain][T-P1]") {
    // Planted fault: quad_flips' second branch (the side that is not
    // counter-clockwise in the frame) dropped.
    const LatticeFrame f{3.0, 1.5};
    const auto fx = frame_collinear();
    const LatticeMesh m = built(fx, f);
    const MeshVertex a = fx.vertices[0], b = fx.vertices[1], d = fx.vertices[2], c = kFrameCollinearC;
    {  // the premise, as L4b has it
        const auto at = qg::locate(m, c);
        REQUIRE(at.has_value());
        REQUIRE(at->t == 0);
        REQUIRE(at->on == 3);
        REQUIRE(orient_sign(a, b, c) > 0);
        REQUIRE(DefaultKernel::orient2d(f.at(a), f.at(b), f.at(c)) == Orientation::Collinear);
        REQUIRE(DefaultKernel::orient2d(f.at(b), f.at(a), f.at(d)) == Orientation::CounterClockwise);
        REQUIRE(DefaultKernel::incircle(f.at(b), f.at(a), f.at(d), f.at(c)) == Incircle::Inside);
    }
    LatticeMesh after = m;
    qg::replay<DefaultKernel>(after, {0, 3}, c, f);
    REQUIRE(qg::diff(m, after).removed.size() == 2);  // Lawson flips A-B
    const auto cav = terrain::mesh::detail::quality_cavity<DefaultKernel>(m, 0, 3, c, f);
    REQUIRE(cav.removed.size() == 2);
    REQUIRE(cav.created.size() == 4);
    check_prediction(qg::Step{m, after, c, {0, 3}}, f);
    REQUIRE(delaunay_violations(after, f) == 0);
}

TEST_CASE("T-P1: a line split flips beyond the triangle across the line, as predicted",
          "[mesh][quality][gain][T-P1][LS1]") {
    // Planted fault: improve's insert routine leaving the slot across a split
    // edge out of legalise_around's seeds.
    const LatticeFrame f{1.0, 1.0};
    const auto fx = line_beyond_far_flip();
    const MeshVertex foot = foot_on_e(fx, kBeyondNode);
    const MeshVertex d = fx.vertices[3], e2 = fx.vertices[1], g = fx.vertices[7];
    {  // the premise: the foot inside (D, G, E2)'s circle, E1 outside it
        REQUIRE(DefaultKernel::incircle(f.at(d), f.at(g), f.at(e2), f.at(foot)) == Incircle::Inside);
        REQUIRE(DefaultKernel::incircle(f.at(d), f.at(g), f.at(e2), f.at(fx.vertices[0])) == Incircle::Outside);
    }
    LatticeMesh m = built(fx, f);
    const auto r = recorded(m, f, options(fx, 0.0, true));
    REQUIRE(r.out.line_splits >= 1);
    REQUIRE_FALSE(r.steps.empty());
    const auto& first = r.steps.front();  // P Q R's foot on e, as in "splits the line at the node's foot"
    CAPTURE(first.p.col, first.p.row, foot.col, foot.row);
    REQUIRE(std::abs(first.p.col - foot.col) <= 1e-12);  // cells
    REQUIRE(std::abs(first.p.row - foot.row) <= 1e-12);
    REQUIRE(qg::diff(first.before, first.after).removed.count(qg::tri(d, g, e2)) == 1);  // D-E2 flipped
    for (const auto& s : r.steps) check_prediction(s, f);
    check_mesh(m, fx, f);
}

// ================================================================== T-P1 and T-P2 over whole runs

TEST_CASE("T-P1 and T-P2: on whole runs every insertion is the predicted cavity and none lowers the worst angle",
          "[mesh][quality][gain][T-P1][T-P2]") {
    // Planted faults: any T-P1 fault that a run reaches; the acceptance test
    // reversed (an accepted candidate that lowers the worst angle). Two
    // fixtures with no interior constraint, so R8 never fires and every
    // insertion is judged: the Q7 ring (off-node outline) and the nine-node
    // mesh of the equality case. Counted over both, for non-vacuity.
    const bool integer = GENERATE(true, false);
    const double gain = GENERATE(0.0, 2.0);
    const bool feet = GENERATE(false, true);
    CAPTURE(integer, gain, feet);
    std::size_t accepted = 0, refused = 0;
    for (const auto& fx : {q7_ring(), equal_gain()}) {
        const LatticeFrame f = integer ? terrain::mesh::lattice_frame(1.0, 1.0, fx.rows, fx.cols) : LatticeFrame{10.0, 5.0};
        LatticeMesh m = built(fx, f, false);
        const auto r = recorded(m, f, options(fx, gain, feet));
        REQUIRE(r.out.line_splits == 0);  // the premise: no interior line, so every insertion is judged
        REQUIRE(r.steps.size() == r.out.inserted + r.out.feet);
        accepted += r.out.inserted;
        refused += r.out.skipped_no_gain;
        for (const auto& s : r.steps) {
            check_prediction(s, f);
            check_gain(s, f, gain);
        }
        check_mesh(m, fx, f);
    }
    REQUIRE(accepted > 0);  // not vacuous: candidates accepted ...
    REQUIRE(refused > 0);   // ... and refused
}

TEST_CASE("T-P1: across the line fixture every insertion, R8's splits included, is the predicted cavity",
          "[mesh][quality][gain][T-P1][LS1]") {
    // R8 runs only with the feet on (pin 4): splits with them, none without.
    const LatticeFrame f{1.0, 1.0};
    const auto fx = line_beyond();
    const bool feet = GENERATE(false, true);
    CAPTURE(feet);
    LatticeMesh m = built(fx, f);
    const auto r = recorded(m, f, options(fx, 0.0, feet));
    if (feet)
        REQUIRE(r.out.line_splits > 0);
    else
        REQUIRE(r.out.line_splits == 0);
    REQUIRE(r.steps.size() == r.out.inserted + r.out.feet + r.out.line_splits);
    for (const auto& s : r.steps) check_prediction(s, f);
    check_mesh(m, fx, f);
}

// ================================================================== T-P2 on hand fixtures

TEST_CASE("T-P2: the only candidate, which would lower the worst angle, is refused and counted",
          "[mesh][quality][gain][T-P2]") {
    // Planted fault: the acceptance test reversed (it would accept this one).
    const LatticeFrame f{1.0, 1.0};
    const auto fx = own_edge();
    const MeshVertex node{10, 1};
    {  // the premise: the node's insertion lowers the worst angle, 18.43 degrees, to a sliver
        LatticeMesh probe = built(fx, f);
        const auto at = qg::locate(probe, node);
        REQUIRE(at.has_value());
        const LatticeMesh before = probe;
        qg::replay<DefaultKernel>(probe, *at, node, f);
        const auto d = qg::diff(before, probe);
        REQUIRE(qg::worst_capped(d.removed, f, kTheta) > 18.0);
        REQUIRE(qg::worst_capped(d.created, f, kTheta) < 1.0);
    }
    LatticeMesh m = built(fx, f);
    const auto before = qg::triangle_set(m);
    const QualityOutcome q = improve<DefaultKernel>(m, f, options(fx, 0.0));
    REQUIRE(q.skipped_no_gain == 1);
    REQUIRE(q.inserted == 0);
    REQUIRE(q.feet == 0);
    REQUIRE(q.line_splits == 0);
    REQUIRE(m.vertices().size() == fx.vertices.size());
    REQUIRE(qg::triangle_set(m) == before);
}

TEST_CASE("T-P2: with gain -1 the same candidate goes in, as the hard rule has it",
          "[mesh][quality][gain][T-P2][T-P3]") {
    const LatticeFrame f{1.0, 1.0};
    const auto fx = own_edge();
    LatticeMesh m = built(fx, f);
    const QualityOutcome q = improve<DefaultKernel>(m, f, options(fx, -1.0));
    REQUIRE(q.inserted == 1);
    REQUIRE(q.skipped_no_gain == 0);
    REQUIRE(has_vertex(m, MeshVertex{10, 1}, 0.0));
}

TEST_CASE("T-P2: a candidate that keeps the worst angle exactly is accepted at gain 0",
          "[mesh][quality][gain][T-P2]") {
    // Planted fault: the slack dropped and the test made strict (new > old + P).
    const LatticeFrame f{1.0, 1.0};
    const auto fx = equal_gain();
    const MeshVertex node{7, 3};
    {  // the premise: new == old to the bit, by a moved copy of the worst triangle
        LatticeMesh probe = built(fx, f);
        const auto shape = terrain::mesh::detail::quality_shape(probe, 7, f);
        REQUIRE(std::round(shape.centre.x) == node.col);
        REQUIRE(std::round(-shape.centre.y) == node.row);
        const auto at = qg::locate(probe, node);
        REQUIRE(at.has_value());
        const LatticeMesh before = probe;
        qg::replay<DefaultKernel>(probe, *at, node, f);
        const auto d = qg::diff(before, probe);
        const auto worst = [&](const std::set<qg::Tri>& ts) {
            return *std::min_element(ts.begin(), ts.end(), [&](const qg::Tri& a, const qg::Tri& b) {
                return qg::min_angle_deg(a, f) < qg::min_angle_deg(b, f);
            });
        };
        const qg::Tri old_t = worst(d.removed), new_t = worst(d.created);
        REQUIRE(qg::min_angle_deg(old_t, f) < kTheta);
        REQUIRE(qg::translate_of(old_t, new_t));
        REQUIRE(qg::min_angle_deg(old_t, f) == qg::min_angle_deg(new_t, f));
        for (const auto& t : d.removed)  // no other angle within a degree of it
            if (!qg::translate_of(t, old_t)) REQUIRE(qg::min_angle_deg(t, f) > qg::min_angle_deg(old_t, f) + 1.0);
        for (const auto& t : d.created)
            if (!qg::translate_of(t, new_t)) REQUIRE(qg::min_angle_deg(t, f) > qg::min_angle_deg(new_t, f) + 1.0);
    }
    LatticeMesh m = built(fx, f);
    const auto r = recorded(m, f, options(fx, 0.0));
    REQUIRE_FALSE(r.steps.empty());
    CHECK(r.steps.front().p == node);  // improve's first insertion is this candidate
    for (const auto& s : r.steps) {
        check_prediction(s, f);
        check_gain(s, f, 0.0);
    }
}

TEST_CASE("T-P2: a foot that keeps its edge's end angle is accepted at gain 0", "[mesh][quality][gain][T-P2]") {
    // Planted fault: the slack dropped and the test made strict (pin 3). CF2's
    // own edge with the feet on: the node N = (10, 1) lies 0.01 cells inside
    // the outline A-B, so R2 inserts its foot F = (10, 0.99) instead. (F, B, C)
    // keeps the angle at B (the same two rays), and F is 10 cells from both B
    // and C, so the triangle is isosceles with the same angle at C: new equals
    // old, 18.43 degrees, in exact arithmetic, and only the slack keeps the
    // decision off F's last bit.
    const LatticeFrame f{1.0, 1.0};
    const auto fx = own_edge();
    const MeshVertex c = fx.vertices[1], b = fx.vertices[2], foot{10, 0.99};
    const double end_angle = std::atan(1.0 / 3.0) * 180.0 / std::numbers::pi;  // 18.43 degrees
    {  // the premise. 1e-9 degrees: coordinates up to 20 cells, checked on this fixture only
        LatticeMesh probe = built(fx, f);
        const auto at = qg::locate(probe, foot);
        REQUIRE(at.has_value());
        REQUIRE(at->on < 3);  // on A-B
        const LatticeMesh before = probe;
        qg::replay<DefaultKernel>(probe, *at, foot, f);
        const auto d = qg::diff(before, probe);
        REQUIRE(d.removed.size() == 1);
        REQUIRE(d.created.size() == 2);
        REQUIRE(d.created.count(qg::tri(foot, c, b)) == 1);  // (F, C, B), counter-clockwise as A C B
        CAPTURE(angle_at(b, foot, c, f), angle_at(c, b, foot, f), end_angle);
        REQUIRE(std::abs(angle_at(b, foot, c, f) - end_angle) <= 1e-9);
        REQUIRE(std::abs(angle_at(c, b, foot, f) - end_angle) <= 1e-9);
        REQUIRE(std::abs(qg::worst_capped(d.removed, f, kTheta) - end_angle) <= 1e-9);
        REQUIRE(std::abs(qg::worst_capped(d.created, f, kTheta) - end_angle) <= 1e-9);
    }
    LatticeMesh m = built(fx, f);
    const auto r = recorded(m, f, options(fx, 0.0, true));
    REQUIRE(r.out.feet >= 1);
    REQUIRE_FALSE(r.steps.empty());
    const auto& first = r.steps.front();  // improve's first insertion is the foot
    CAPTURE(first.p.col, first.p.row);
    CHECK(std::abs(first.p.col - foot.col) <= 1e-12);  // cells
    CHECK(std::abs(first.p.row - foot.row) <= 1e-12);
    for (const auto& s : r.steps) {
        check_prediction(s, f);
        check_gain(s, f, 0.0);
    }
}

// ================================================================== LS1, R8

TEST_CASE("LS1: a node beyond a long constraint line splits the line at the node's foot",
          "[mesh][quality][gain][LS1]") {
    // Planted fault: the split at the midpoint. The feet on: R8 runs only
    // with them (pin 4).
    const LatticeFrame f{1.0, 1.0};
    const auto fx = line_beyond();
    const MeshVertex foot = foot_on_e(fx, kBeyondNode);
    const MeshVertex mid{(fx.vertices[0].col + fx.vertices[1].col) / 2, (fx.vertices[0].row + fx.vertices[1].row) / 2};
    {  // the premise: P Q R's node is (15, 18), beyond e, its foot well inside e and off its midpoint
        const auto c = centre(fx.vertices[4], fx.vertices[5], fx.vertices[6]);
        REQUIRE(std::round(c.col) == kBeyondNode.col);
        REQUIRE(std::round(c.row) == kBeyondNode.row);
        REQUIRE(orient_sign(fx.vertices[0], fx.vertices[1], kBeyondNode) != orient_sign(fx.vertices[0], fx.vertices[1], fx.vertices[6]));
        REQUIRE(std::hypot(foot.col - mid.col, foot.row - mid.row) > 0.5);
        REQUIRE(std::hypot(foot.col - fx.vertices[0].col, foot.row - fx.vertices[0].row) > 10.0);
    }
    LatticeMesh m = built(fx, f);
    const auto r = recorded(m, f, options(fx, 0.0, true));
    REQUIRE(r.out.line_splits >= 1);
    REQUIRE_FALSE(r.steps.empty());
    const auto& first = r.steps.front();  // P Q R is the worst triangle, so the first candidate
    CAPTURE(first.p.col, first.p.row, foot.col, foot.row);
    CHECK(std::abs(first.p.col - foot.col) <= 1e-12);  // cells; the orthogonal projection
    CHECK(std::abs(first.p.row - foot.row) <= 1e-12);
    REQUIRE(first.at.on < 3);  // a split of e itself
    REQUIRE(std::minmax(first.before.triangles()[first.at.t][first.at.on],
                        first.before.triangles()[first.at.t][(first.at.on + 1) % 3])
            == std::minmax(std::uint32_t{0}, std::uint32_t{1}));
    REQUIRE(m.vertices().size() == fx.vertices.size() + r.out.inserted + r.out.feet + r.out.line_splits);
    check_mesh(m, fx, f);  // the lines are the same set, bits and masks on both halves
}

TEST_CASE("LS1: gain -1 turns the split off", "[mesh][quality][gain][LS1][T-P3]") {
    // Planted fault: R8 left on at gain -1. The feet on, so the gain is what refuses.
    const LatticeFrame f{1.0, 1.0};
    const auto fx = line_beyond();
    LatticeMesh m = built(fx, f);
    const QualityOutcome q = improve<DefaultKernel>(m, f, options(fx, -1.0, true));
    REQUIRE(q.line_splits == 0);
    REQUIRE(q.skipped_blocked >= 1);
    REQUIRE(edges_with_mask(m, kLineMask) == 1);  // e is whole
    REQUIRE_FALSE(has_vertex(m, foot_on_e(fx, kBeyondNode), 1e-9));
}

TEST_CASE("LS1: with the feet off, no split", "[mesh][quality][gain][LS1]") {
    // Planted fault: R8 left on with the feet off (pin 4: --no-constraint-feet
    // keeps every point off the lines).
    const LatticeFrame f{1.0, 1.0};
    const auto fx = line_beyond();
    LatticeMesh m = built(fx, f);
    const QualityOutcome q = improve<DefaultKernel>(m, f, options(fx, 0.0, false));
    REQUIRE(q.line_splits == 0);
    REQUIRE(q.skipped_blocked >= 1);
    REQUIRE(edges_with_mask(m, kLineMask) == 1);  // e is whole
    REQUIRE_FALSE(has_vertex(m, foot_on_e(fx, kBeyondNode), 1e-9));
}

TEST_CASE("LS1: no split within half a cell of an end of the line", "[mesh][quality][gain][LS1]") {
    // Planted fault: the end check removed (it splits e at (2, 10)).
    const LatticeFrame f{1.0, 1.0};
    const auto fx = line_near_end();
    {  // the premise: B's node (2, 12) is beyond e, its foot (2, 10) inside e, 0.25 from E2
        const auto c = centre(fx.vertices[0], fx.vertices[1], fx.vertices[4]);
        REQUIRE(std::round(c.col) == 2.0);
        REQUIRE(std::round(c.row) == 12.0);
        REQUIRE(orient_sign(fx.vertices[0], fx.vertices[1], MeshVertex{2, 12}) < 0);  // beyond e, seen from B
        const double to_end = fx.vertices[1].col - 2.0;
        REQUIRE(to_end > 0.0);
        REQUIRE(to_end < 0.5);  // delta_q = min(dx, dy) / 2
        LatticeMesh probe = built(fx, f);
        REQUIRE(terrain::mesh::detail::foot_fits(probe, 0, 0, MeshVertex{2, 10}));  // only the end check refuses it
    }
    LatticeMesh m = built(fx, f);
    const QualityOutcome q = improve<DefaultKernel>(m, f, options(fx, 0.0, true));
    REQUIRE(q.line_splits == 0);
    REQUIRE(q.skipped_blocked == 1);  // pinned: a declined split stays one blocked skip
    REQUIRE(m.vertices().size() == fx.vertices.size());
}

TEST_CASE("LS1: an outline edge with nothing beyond is not split", "[mesh][quality][gain][LS1]") {
    const LatticeFrame f{1.0, 1.0};
    const auto fx = line_beyond_outline();
    LatticeMesh m = built(fx, f);
    REQUIRE(m.neighbours(0)[1] == kNoNeighbour);  // the premise: (Q, E1, E2)'s edge E1-E2 has no neighbour
    const QualityOutcome q = improve<DefaultKernel>(m, f, options(fx, 0.0, true));
    REQUIRE(q.line_splits == 0);
    REQUIRE(q.skipped_blocked >= 1);
    REQUIRE(edges_with_mask(m, kLineMask) == 1);
    REQUIRE_FALSE(has_vertex(m, foot_on_e(fx, kBeyondNode), 1e-9));
}

TEST_CASE("LS1: a frozen line is not split", "[mesh][quality][gain][LS1]") {
    const LatticeFrame f{1.0, 1.0};
    const auto fx = line_beyond();
    LatticeMesh m = built(fx, f);
    m.set_frozen_mask(kLineMask);
    const QualityOutcome q = improve<DefaultKernel>(m, f, options(fx, 0.0, true));
    REQUIRE(q.line_splits == 0);
    REQUIRE(q.skipped_blocked >= 1);
    REQUIRE(edges_with_mask(m, kLineMask) == 1);
}

TEST_CASE("LS1: a foot the validity callable refuses is not split in", "[mesh][quality][gain][LS1]") {
    // Nodes only: every foot on the tilted e is off-node, so none goes in.
    const LatticeFrame f{1.0, 1.0};
    const auto fx = line_beyond();
    LatticeMesh m = built(fx, f);
    const auto r = recorded(m, f, options(fx, 0.0, true), [](const MeshVertex& v) { return v.is_node(); });
    REQUIRE(r.out.line_splits == 0);
    REQUIRE(edges_with_mask(m, kLineMask) == 1);
    REQUIRE_FALSE(has_vertex(m, foot_on_e(fx, kBeyondNode), 1e-9));
}
