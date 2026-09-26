// Increment 20 (docs/increments/20-start-quality.md, R2 to R6, R10): the
// minimum-angle pass over a LatticeMesh. INVARIANT-CRITICAL: the
// postcondition and the constraints are decided here, so @reviewer
// mutation-tests this suite. Q1 to Q5 of the design, plus Q6 (frame), Q8
// (stride grid) and Q9 (outside the rectangle) at the header's own level.
//
// Provisional choices this suite encodes, pending Ola: C1 (a) no Steiner point
// on an input segment except a DEM node lying exactly on it, C2 (a) 25°.
//
// Interface: include/terrain/mesh/quality.hpp. The design names QualityOptions,
// QualityOutcome and improve<K>; the field names are left open and chosen here:
//
//   struct QualityOptions {
//       double min_angle_deg;      // declared in this order: the suite uses
//       std::size_t rows;          // designated initialisers
//       std::size_t cols;
//   };
//   struct QualityOutcome {
//       std::size_t inserted;         // DEM nodes added
//       std::size_t skipped_floor;    // R4 step 3, in R4's order
//       std::size_t skipped_outside;  // circumcentre outside the node rectangle
//       std::size_t skipped_vertex;   // the snapped node is already a vertex
//       std::size_t skipped_blocked;  // the walk crossed a constraint or left the mesh
//       std::size_t walk_bound_hits;  // R4 step 4's bound, asserted zero (Q5)
//   };
//   template <pred::GeometryKernel K>
//   QualityOutcome improve(LatticeMesh&, const LatticeFrame&, const QualityOptions&);
//
// improve is called on a mesh that legalise_all has already made constrained
// Delaunay, as refine does (R1).
//
// THE ORACLE IS THIS FILE'S OWN. "Bad" is recomputed from the circumradius
// (a b c / (4 area)) and the shortest edge in the world frame
// (col * dx, -(row * dy)), with the circumcentre from the textbook formula,
// not from the pass. Q1's skip reasons are recomputed from the output: the
// circumcentre outside [0, cols-1] x [0, rows-1]; the rounded node already a
// vertex; the rounded node outside the (convex) domain, by exact orientation
// against the input ring. That last one stands in for "the walk crossed a
// constraint or left the mesh" only because Q1's domain is convex and its
// only constraints are its boundary, so a node inside it is on the inner side
// of every constraint's line and a visibility walk cannot meet one. On
// non-convex constraints (Q4) the walk's reason is a property of its path, not
// of the output, and Q4 checks counts instead; see the handback.
//
// Delaunay and orientation use DefaultKernel on the frame points, the
// producer's predicate but not its records (14b's convention).

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/mesh/quality.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include "refinement_fixtures.hpp"

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
using terrain::mesh::QualityOptions;
using terrain::mesh::QualityOutcome;
using terrain::pred::DefaultKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;

namespace {

constexpr double kTheta = 25.0;
constexpr double kRel = 1e-9;  // slack on the double-valued criterion, never on topology

using Pair = std::pair<std::uint32_t, std::uint32_t>;

Pair undirected(std::uint32_t a, std::uint32_t b) { return std::minmax(a, b); }

// ------------------------------------------------------------ the fixtures

struct Fixture {
    std::vector<MeshVertex> vertices;
    std::vector<TriangleIndices> triangles;
    std::map<Pair, std::uint32_t> constraints;  // undirected vertex pair -> mask
    std::vector<MeshVertex> ring;               // Q1's convex domain, CCW in world
    std::size_t rows = 0;
    std::size_t cols = 0;
};

LatticeMesh build(const Fixture& f) {
    std::vector<std::uint8_t> bits(f.triangles.size(), 0);
    std::vector<std::array<std::uint32_t, 3>> masks(f.triangles.size(), {0, 0, 0});
    for (std::size_t t = 0; t < f.triangles.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            if (const auto it = f.constraints.find(undirected(f.triangles[t][k], f.triangles[t][(k + 1) % 3]));
                it != f.constraints.end()) {
                bits[t] |= static_cast<std::uint8_t>(1u << k);
                masks[t][k] = it->second;
            }
    auto m = LatticeMesh::build(f.vertices, f.triangles, std::move(bits), std::move(masks));
    REQUIRE(m.has_value());
    return std::move(*m);
}

// A small analogue of the quarter circle (design, M1): the centre on a node
// near the bottom-right corner of a 64 x 64 grid, a straight edge up along
// column 62 and one left along row 62 (both lattice lines), and a dense arc of
// radius 60 through the upper left, every interior arc vertex off-node. The
// start is the fan from the centre, the pathology increment 20 removes.
// Masks: east edge 1, arc 2, south edge 4.
Fixture quarter_circle(std::size_t arc_segments = 70) {
    Fixture f;
    f.rows = f.cols = 64;
    const double c0 = 62.0, r0 = 62.0, radius = 60.0;
    f.vertices.push_back(MeshVertex{c0, r0});
    for (std::size_t i = 0; i <= arc_segments; ++i) {
        const double t = std::numbers::pi / 2.0
                       + std::numbers::pi / 2.0 * static_cast<double>(i) / static_cast<double>(arc_segments);
        f.vertices.push_back(MeshVertex{c0 + radius * std::cos(t), r0 - radius * std::sin(t)});
    }
    f.vertices[1] = MeshVertex{c0, r0 - radius};           // on the east lattice line, a node
    f.vertices[arc_segments + 1] = MeshVertex{c0 - radius, r0};  // on the south lattice line, a node
    const auto n = static_cast<std::uint32_t>(f.vertices.size());
    for (std::uint32_t i = 1; i + 1 < n; ++i) {
        f.triangles.push_back({0, i, i + 1});
        f.constraints[undirected(i, i + 1)] = 2;
    }
    f.constraints[undirected(0, 1)] = 1;
    f.constraints[undirected(n - 1, 0)] = 4;
    f.ring = f.vertices;
    return f;
}

// Q4: a 64 x 64-cell square (65 x 65 nodes) with two breakline segments
// meeting at P at 10°, 40 cells long, off-node ends. Triangulated by hand;
// see the handback for the orientation arithmetic. Square sides mask 1, the
// breaklines 16 and 32.
Fixture ten_degree_wedge() {
    Fixture f;
    f.rows = f.cols = 65;
    const double deg = std::numbers::pi / 180.0, len = 40.0;
    const MeshVertex p{10.3, 32.2};
    // Frame y = -row, so a positive frame angle is a smaller row.
    const MeshVertex a{p.col + len * std::cos(5.0 * deg), p.row - len * std::sin(5.0 * deg)};
    const MeshVertex b{p.col + len * std::cos(-5.0 * deg), p.row - len * std::sin(-5.0 * deg)};
    // 0 TL, 1 TR, 2 BR, 3 BL, 4 P, 5 A (upper), 6 B (lower)
    f.vertices = {MeshVertex{0, 0}, MeshVertex{64, 0}, MeshVertex{64, 64}, MeshVertex{0, 64}, p, a, b};
    f.triangles = {{0, 3, 4}, {3, 6, 4}, {3, 2, 6}, {2, 5, 6}, {2, 1, 5}, {1, 0, 5}, {0, 4, 5}, {4, 6, 5}};
    f.constraints = {{undirected(0, 1), 1}, {undirected(1, 2), 1}, {undirected(2, 3), 1},
                     {undirected(3, 0), 1}, {undirected(4, 5), 16}, {undirected(4, 6), 32}};
    return f;
}

// One triangle, every side constrained with mask 1.
Fixture single(MeshVertex a, MeshVertex b, MeshVertex c, std::size_t rows, std::size_t cols) {
    Fixture f;
    f.rows = rows;
    f.cols = cols;
    f.vertices = {a, b, c};
    f.triangles = {{0, 1, 2}};
    f.constraints = {{undirected(0, 1), 1}, {undirected(1, 2), 1}, {undirected(2, 0), 1}};
    f.ring = f.vertices;
    return f;
}

QualityOutcome run(LatticeMesh& m, const LatticeFrame& frame, const Fixture& f, double theta = kTheta) {
    legalise_all<DefaultKernel>(m, frame, [](std::uint32_t) {});
    return improve<DefaultKernel>(m, frame, QualityOptions{.min_angle_deg = theta, .rows = f.rows, .cols = f.cols});
}

// ------------------------------------------------------------- the oracle

Point2 lattice_point(MeshVertex v) { return Point2{v.col, -v.row}; }

Point2 world(const LatticeFrame& f, MeshVertex v) { return Point2{v.col * f.dx, -(v.row * f.dy)}; }

struct Shape {
    double circumradius;
    double shortest;
    Point2 centre;  // world frame
};

Shape shape(const LatticeFrame& f, MeshVertex va, MeshVertex vb, MeshVertex vc) {
    const Point2 a = world(f, va), b = world(f, vb), c = world(f, vc);
    const double ab = std::hypot(b.x - a.x, b.y - a.y), bc = std::hypot(c.x - b.x, c.y - b.y),
                 ca = std::hypot(a.x - c.x, a.y - c.y);
    const double cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    const double d = 2.0 * cross;
    const double a2 = a.x * a.x + a.y * a.y, b2 = b.x * b.x + b.y * b.y, c2 = c.x * c.x + c.y * c.y;
    const Point2 centre{(a2 * (b.y - c.y) + b2 * (c.y - a.y) + c2 * (a.y - b.y)) / d,
                        (a2 * (c.x - b.x) + b2 * (a.x - c.x) + c2 * (b.x - a.x)) / d};
    return Shape{ab * bc * ca / (2.0 * std::abs(cross)), std::min({ab, bc, ca}), centre};
}

// Minimum angle below theta, by the design's R2 ratio, with slack so a
// triangle the producer's double expression calls exactly at the bound is not
// flagged either way.
bool clearly_bad(const Shape& s, double theta) {
    return s.circumradius / s.shortest > (1.0 + kRel) / (2.0 * std::sin(theta * std::numbers::pi / 180.0));
}

double floor_of(const LatticeFrame& f) { return std::hypot(f.dx, f.dy); }

bool clearly_above_floor(const Shape& s, const LatticeFrame& f) {
    return s.circumradius >= floor_of(f) * (1.0 + kRel);
}

// The DEM node nearest the circumcentre, R3; empty when the centre lies
// outside the node rectangle (Q9: skipped, never clamped).
std::optional<MeshVertex> snapped(const Shape& s, const LatticeFrame& f, std::size_t rows, std::size_t cols) {
    const double col = s.centre.x / f.dx, row = -s.centre.y / f.dy;
    const double eps = kRel * static_cast<double>(std::max(rows, cols));
    if (col < eps || row < eps || col > static_cast<double>(cols - 1) - eps
        || row > static_cast<double>(rows - 1) - eps)
        return std::nullopt;  // outside, or near enough to the border to count as either
    return MeshVertex{std::round(col), std::round(row)};
}

bool outside_convex_ring(const std::vector<MeshVertex>& ring, MeshVertex p) {
    for (std::size_t i = 0; i < ring.size(); ++i)
        if (DefaultKernel::orient2d(lattice_point(ring[i]), lattice_point(ring[(i + 1) % ring.size()]),
                                    lattice_point(p))
            == Orientation::Clockwise)
            return true;
    return false;
}

// Topology from the triangles alone: positive orientation, each directed edge
// once, neighbour links and constraint bits and masks consistent, no two
// vertices equal, and the output rebuilds.
void check_topology(const LatticeMesh& m) {
    const auto v = m.vertices();
    std::set<std::pair<double, double>> seen;
    for (const auto& p : v) REQUIRE(seen.insert({p.col, p.row}).second);  // no duplicate vertex

    std::map<Pair, std::uint32_t> directed;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        const auto& tri = m.triangles()[t];
        CAPTURE(t);
        REQUIRE(DefaultKernel::orient2d(lattice_point(v[tri[0]]), lattice_point(v[tri[1]]), lattice_point(v[tri[2]]))
                == Orientation::CounterClockwise);
        for (unsigned k = 0; k < 3; ++k) REQUIRE(directed.emplace(Pair{tri[k], tri[(k + 1) % 3]}, t).second);
    }
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        const auto& tri = m.triangles()[t];
        for (unsigned k = 0; k < 3; ++k) {
            CAPTURE(t, k);
            const auto it = directed.find({tri[(k + 1) % 3], tri[k]});
            REQUIRE(m.neighbours(t)[k] == (it == directed.end() ? kNoNeighbour : it->second));
            if (it == directed.end()) continue;
            const auto& o = m.triangles()[it->second];
            unsigned j = 0;
            while (o[j] != tri[(k + 1) % 3]) ++j;
            REQUIRE(m.is_constrained(t, k) == m.is_constrained(it->second, j));
            REQUIRE(m.mask(t, k) == m.mask(it->second, j));
        }
    }
    std::vector<MeshVertex> vs(v.begin(), v.end());
    std::vector<TriangleIndices> ts(m.triangles().begin(), m.triangles().end());
    REQUIRE(LatticeMesh::build(vs, ts, std::vector<std::uint8_t>(ts.size(), 0),
                               std::vector<std::array<std::uint32_t, 3>>(ts.size(), {0, 0, 0}))
                .has_value());
}

// Q2: every unconstrained interior edge is locally Delaunay in the frame.
std::size_t delaunay_violations(const LatticeMesh& m, const LatticeFrame& f) {
    const auto v = m.vertices();
    std::size_t bad = 0;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto u = m.neighbours(t)[k];
            if (u == kNoNeighbour || m.is_constrained(t, k)) continue;
            const auto& tt = m.triangles()[t];
            std::uint32_t apex = 0;
            for (const auto x : m.triangles()[u])
                if (x != tt[k] && x != tt[(k + 1) % 3]) apex = x;
            if (DefaultKernel::incircle(f.at(v[tt[0]]), f.at(v[tt[1]]), f.at(v[tt[2]]), f.at(v[apex]))
                == Incircle::Inside)
                ++bad;
        }
    return bad;
}

// Q2's second half: every vertex the pass added is a node in the rectangle,
// and there are exactly `inserted` of them.
void check_inserted_are_nodes(const LatticeMesh& m, const Fixture& f, const QualityOutcome& q) {
    const auto v = m.vertices();
    REQUIRE(v.size() == f.vertices.size() + q.inserted);
    for (std::size_t i = 0; i < f.vertices.size(); ++i) REQUIRE(v[i] == f.vertices[i]);  // start vertices untouched
    for (std::size_t i = f.vertices.size(); i < v.size(); ++i) {
        CAPTURE(i, v[i].col, v[i].row);
        REQUIRE(v[i].is_node());
        REQUIRE(v[i].col >= 0.0);
        REQUIRE(v[i].row >= 0.0);
        REQUIRE(v[i].col <= static_cast<double>(f.cols - 1));
        REQUIRE(v[i].row <= static_cast<double>(f.rows - 1));
    }
}

// Q3: the constraints as a set of lines are unchanged. Each output constraint
// edge lies on exactly one input constraint segment, exactly (orient2d
// Collinear, within its extent), with that segment's mask; the pieces of each
// input segment chain from one end to the other without gaps or overlaps; a
// piece endpoint that is not an input vertex is a node.
void check_constraints(const LatticeMesh& m, const Fixture& f) {
    const auto v = m.vertices();
    const auto [edges, masks] = m.constraint_edges();
    auto on_closed = [&](MeshVertex a, MeshVertex b, MeshVertex p) {
        return DefaultKernel::orient2d(lattice_point(a), lattice_point(b), lattice_point(p)) == Orientation::Collinear
            && std::min(a.col, b.col) <= p.col && p.col <= std::max(a.col, b.col)
            && std::min(a.row, b.row) <= p.row && p.row <= std::max(a.row, b.row);
    };
    std::map<Pair, std::vector<Pair>> pieces;  // input segment -> output pieces
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const MeshVertex p = v[edges[i][0]], q = v[edges[i][1]];
        CAPTURE(i, p.col, p.row, q.col, q.row);
        std::size_t owners = 0;
        for (const auto& [seg, mask] : f.constraints) {
            const MeshVertex a = f.vertices[seg.first], b = f.vertices[seg.second];
            if (!on_closed(a, b, p) || !on_closed(a, b, q)) continue;
            ++owners;
            REQUIRE(masks[i] == mask);
            pieces[seg].push_back({edges[i][0], edges[i][1]});
        }
        REQUIRE(owners == 1);
        for (const auto x : edges[i])
            if (x >= f.vertices.size()) REQUIRE(v[x].is_node());
    }
    REQUIRE(pieces.size() == f.constraints.size());  // no input segment lost
    for (const auto& [seg, list] : pieces) {
        // Walk from seg.first to seg.second along the pieces.
        std::map<std::uint32_t, std::vector<std::uint32_t>> adj;
        for (const auto& [a, b] : list) {
            adj[a].push_back(b);
            adj[b].push_back(a);
        }
        std::uint32_t at = seg.first, prev = kNoNeighbour;
        std::size_t steps = 0;
        while (at != seg.second) {
            REQUIRE(adj[at].size() == (at == seg.first ? 1u : 2u));
            const auto next = adj[at][0] != prev ? adj[at][0] : adj[at][1];
            prev = at;
            at = next;
            REQUIRE(++steps <= list.size());
        }
        REQUIRE(steps == list.size());
    }
}

// Q1: every triangle meets theta, or has R below the floor, or is explained by
// one of R4's skip reasons, recomputed here. Returns how many needed a reason.
std::size_t check_postcondition(const LatticeMesh& m, const LatticeFrame& f, const Fixture& fx, double theta) {
    const auto v = m.vertices();
    std::set<std::pair<double, double>> vertex_at;
    for (const auto& p : v) vertex_at.insert({p.col, p.row});
    std::size_t explained = 0;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        const auto& tri = m.triangles()[t];
        const Shape s = shape(f, v[tri[0]], v[tri[1]], v[tri[2]]);
        if (!clearly_bad(s, theta) || !clearly_above_floor(s, f)) continue;
        const auto node = snapped(s, f, fx.rows, fx.cols);
        const bool outside = !node;
        const bool vertex = node && vertex_at.count({node->col, node->row}) != 0;
        const bool blocked = node && outside_convex_ring(fx.ring, *node);
        CAPTURE(t, s.circumradius, s.shortest, s.centre.x, s.centre.y);
        INFO("a bad triangle above the floor with no skip reason: the pass left work undone");
        REQUIRE((outside || vertex || blocked));
        ++explained;
    }
    return explained;
}

std::size_t skips(const QualityOutcome& q) {
    return q.skipped_outside + q.skipped_vertex + q.skipped_blocked + q.walk_bound_hits;
}

double worst_angle_deg(const LatticeMesh& m, const LatticeFrame& f) {
    const auto v = m.vertices();
    double worst = 180.0;
    for (const auto& tri : m.triangles()) {
        const Shape s = shape(f, v[tri[0]], v[tri[1]], v[tri[2]]);
        worst = std::min(worst, std::asin(std::min(1.0, s.shortest / (2.0 * s.circumradius))) * 180.0 / std::numbers::pi);
    }
    return worst;
}

}  // namespace

// ------------------------------------------------------------------ Q1 to Q3

TEST_CASE("Q1 Q2 Q3: the quarter-circle analogue meets the postcondition", "[mesh][quality]") {
    const double dy = GENERATE(10.0, 5.0);
    CAPTURE(dy);
    const LatticeFrame frame{10.0, dy};
    const Fixture fx = quarter_circle();
    LatticeMesh m = build(fx);
    const double before = worst_angle_deg(m, frame);
    const QualityOutcome q = run(m, frame, fx);
    CAPTURE(q.inserted, q.skipped_floor, q.skipped_outside, q.skipped_vertex, q.skipped_blocked);

    REQUIRE(before < 5.0);   // the fans are there to be removed
    REQUIRE(q.inserted > 0);
    REQUIRE(q.walk_bound_hits == 0);  // Q5

    check_topology(m);
    REQUIRE(delaunay_violations(m, frame) == 0);  // Q2
    check_inserted_are_nodes(m, fx, q);           // Q2
    check_constraints(m, fx);                      // Q3
    const std::size_t explained = check_postcondition(m, frame, fx, kTheta);  // Q1
    REQUIRE(explained <= skips(q));
}

TEST_CASE("Q3: a node exactly on a constrained edge splits it and keeps bit and mask", "[mesh][quality]") {
    // A right triangle's circumcentre is its hypotenuse's midpoint. Hypotenuse
    // (0, 1)-(20, 1) along row 1, right angle at (2, 7): 18.4° at (20, 1), so
    // bad; R = 10 cells; the circumcentre is node (10, 1), exactly on the
    // constrained edge (Collinear in the exact kernel). R6: split_edge there.
    const LatticeFrame frame{1.0, 1.0};
    const Fixture fx = single(MeshVertex{0, 1}, MeshVertex{2, 7}, MeshVertex{20, 1}, 8, 21);
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, frame, fx);
    REQUIRE(q.inserted > 0);
    REQUIRE(q.walk_bound_hits == 0);
    const auto v = m.vertices();
    REQUIRE(std::find(v.begin(), v.end(), MeshVertex{10, 1}) != v.end());
    const auto [edges, masks] = m.constraint_edges();
    std::size_t on_row_1 = 0;
    for (std::size_t i = 0; i < edges.size(); ++i)
        if (v[edges[i][0]].row == 1.0 && v[edges[i][1]].row == 1.0) {
            ++on_row_1;
            REQUIRE(masks[i] == 1u);
        }
    REQUIRE(on_row_1 >= 2);
    check_topology(m);
    REQUIRE(delaunay_violations(m, frame) == 0);
    check_inserted_are_nodes(m, fx, q);
    check_constraints(m, fx);
    REQUIRE(check_postcondition(m, frame, fx, kTheta) <= skips(q));
}

TEST_CASE("Q1: the arc keeps exactly its input vertices", "[mesh][quality]") {
    // C1 (a): the pass never computes a point on a segment. The arc's chords
    // pass through no node exactly here, so every arc constraint edge in the
    // output is an input chord, unsplit.
    const LatticeFrame frame{10.0, 10.0};
    const Fixture fx = quarter_circle();
    LatticeMesh m = build(fx);
    run(m, frame, fx);
    const auto [edges, masks] = m.constraint_edges();
    std::size_t arc = 0;
    for (std::size_t i = 0; i < edges.size(); ++i)
        if (masks[i] == 2u) {
            ++arc;
            REQUIRE(edges[i][0] < fx.vertices.size());
            REQUIRE(edges[i][1] < fx.vertices.size());
        }
    REQUIRE(arc == fx.vertices.size() - 2);
}

TEST_CASE("Q1: theta 0 is off and leaves the mesh untouched", "[mesh][quality]") {
    const LatticeFrame frame{10.0, 10.0};
    const Fixture fx = quarter_circle();
    LatticeMesh legalised = build(fx);
    legalise_all<DefaultKernel>(legalised, frame, [](std::uint32_t) {});
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, frame, fx, 0.0);
    REQUIRE(q.inserted == 0);
    REQUIRE(skips(q) + q.skipped_floor == 0);
    REQUIRE(std::vector<TriangleIndices>(m.triangles().begin(), m.triangles().end())
            == std::vector<TriangleIndices>(legalised.triangles().begin(), legalised.triangles().end()));
    REQUIRE(m.vertices().size() == legalised.vertices().size());
}

TEST_CASE("Q1: the postcondition holds at 20 and 30 degrees too", "[mesh][quality]") {
    // 20° and 30° as well as 25°: the postcondition is stated for any theta.
    const double theta = GENERATE(20.0, 30.0);
    CAPTURE(theta);
    const LatticeFrame frame{10.0, 10.0};
    const Fixture fx = quarter_circle();
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, frame, fx, theta);
    REQUIRE(q.walk_bound_hits == 0);
    check_topology(m);
    REQUIRE(delaunay_violations(m, frame) == 0);
    check_inserted_are_nodes(m, fx, q);
    check_constraints(m, fx);
    REQUIRE(check_postcondition(m, frame, fx, theta) <= skips(q));
}

// ------------------------------------------------------------------------ Q4

TEST_CASE("Q4: two breaklines at 10 degrees terminate within the bound", "[mesh][quality]") {
    // R10: every insertion is a new DEM node, so the pass ends. Stated bound:
    // one insertion per four lattice nodes, 1056 here; the floor is what keeps
    // the wedge's apex from consuming every node near it.
    const LatticeFrame frame{1.0, 1.0};
    const Fixture fx = ten_degree_wedge();
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, frame, fx);
    CAPTURE(q.inserted, q.skipped_floor, q.skipped_outside, q.skipped_vertex, q.skipped_blocked);
    REQUIRE(q.inserted > 0);
    REQUIRE(q.inserted <= fx.rows * fx.cols / 4);
    REQUIRE(q.walk_bound_hits == 0);  // Q5
    check_topology(m);
    REQUIRE(delaunay_violations(m, frame) == 0);
    check_inserted_are_nodes(m, fx, q);
    check_constraints(m, fx);

    // C1 (a): the 10° corner at P cannot be improved without a vertex on a
    // breakline, so some triangle there keeps an angle of at most 10°.
    REQUIRE(worst_angle_deg(m, frame) <= 10.0 + 1e-6);

    // What remains bad and above the floor is at most what the pass says it
    // skipped (the wedge's triangles are recorded skips).
    const auto v = m.vertices();
    std::size_t remaining = 0;
    for (const auto& tri : m.triangles()) {
        const Shape s = shape(frame, v[tri[0]], v[tri[1]], v[tri[2]]);
        if (clearly_bad(s, kTheta) && clearly_above_floor(s, frame)) ++remaining;
    }
    REQUIRE(remaining <= skips(q));
}

// ------------------------------------------------------------------------ Q6

TEST_CASE("Q6: a triangle good in the world but thin in the lattice is left alone", "[mesh][quality][frame]") {
    // dx = 1, dy = 3. Lattice (0,0), (15,3), (30,0) is 11.3° at its base in
    // (col, -row) but 31° in the world, so at 25° there is nothing to do.
    const LatticeFrame frame{1.0, 3.0};
    const Fixture fx = single(MeshVertex{0, 0}, MeshVertex{15, 3}, MeshVertex{30, 0}, 4, 31);
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, frame, fx);
    REQUIRE(q.inserted == 0);
    REQUIRE(skips(q) + q.skipped_floor == 0);
    REQUIRE(m.triangle_count() == 1);
}

TEST_CASE("Q6: a triangle good in the lattice but thin in the world is improved", "[mesh][quality][frame]") {
    // dx = 1, dy = 3. Lattice (0,0), (15,26), (30,0) is about 60° everywhere in
    // (col, -row) and 21.8° at its apex in the world. World circumcentre at
    // row 12.52, so node (15, 13), inside the triangle; R = 40.4 >= sqrt(10).
    const LatticeFrame frame{1.0, 3.0};
    Fixture fx = single(MeshVertex{0, 0}, MeshVertex{15, 26}, MeshVertex{30, 0}, 27, 31);
    fx.ring = fx.vertices;
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, frame, fx);
    REQUIRE(q.inserted > 0);
    REQUIRE(q.walk_bound_hits == 0);
    const auto v = m.vertices();
    REQUIRE(std::find(v.begin(), v.end(), MeshVertex{15, 13}) != v.end());
    check_topology(m);
    REQUIRE(delaunay_violations(m, frame) == 0);
    check_inserted_are_nodes(m, fx, q);
    check_constraints(m, fx);
    REQUIRE(check_postcondition(m, frame, fx, kTheta) <= skips(q));
}

// ------------------------------------------------------------------------ Q8

TEST_CASE("Q8: a stride grid with square cells gets nothing", "[mesh][quality]") {
    // Right isosceles triangles: R / shortest = sqrt(2) / 2, under 1 / (2 sin 25°).
    const std::size_t stride = GENERATE(4u, 8u);
    CAPTURE(stride);
    const auto g = terrain::raster::RasterGeometry{500000.0, 7000000.0, 10.0, 10.0, 33, 33};
    const auto s = refinement_fixtures::grid_mesh(g, stride);
    Fixture fx;
    fx.rows = fx.cols = 33;
    for (const auto& p : s.lattice)
        fx.vertices.push_back(MeshVertex{static_cast<double>(p.col), static_cast<double>(p.row)});
    fx.triangles.assign(s.mesh.triangles().begin(), s.mesh.triangles().end());
    for (std::size_t i = 0; i < s.edges.size(); ++i) fx.constraints[undirected(s.edges[i][0], s.edges[i][1])] = s.masks[i];
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, LatticeFrame{10.0, 10.0}, fx);
    REQUIRE(q.inserted == 0);
    REQUIRE(skips(q) + q.skipped_floor == 0);
    REQUIRE(std::vector<TriangleIndices>(m.triangles().begin(), m.triangles().end()) == fx.triangles);
}

// ------------------------------------------------------------------------ Q9

TEST_CASE("Q9: a circumcentre outside the node rectangle is skipped not clamped", "[mesh][quality]") {
    // Base along row 0 from col 0 to col 40, apex at (20, 2): about 5.7° at the
    // base, circumcentre at row -99, outside. Clamping would put a node on
    // row 0 (on the constrained base) and split it.
    const LatticeFrame frame{10.0, 10.0};
    const Fixture fx = single(MeshVertex{0, 0}, MeshVertex{20, 2}, MeshVertex{40, 0}, 3, 41);
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, frame, fx);
    REQUIRE(q.inserted == 0);
    REQUIRE(q.skipped_outside == 1);
    REQUIRE(q.skipped_vertex + q.skipped_blocked + q.walk_bound_hits == 0);
    REQUIRE(m.vertices().size() == 3);
    REQUIRE(m.triangle_count() == 1);
}

TEST_CASE("Q1: a bad triangle below the floor is left alone", "[mesh][quality]") {
    // Base 6 m, height 1 m on a 10 m grid: 18.4° at the base, R = 5 m against
    // a floor of 14.1 m, circumcentre inside the grid. Not nodes.
    const LatticeFrame frame{10.0, 10.0};
    const Fixture fx = single(MeshVertex{2.0, 5.0}, MeshVertex{2.3, 5.1}, MeshVertex{2.6, 5.0}, 11, 11);
    LatticeMesh m = build(fx);
    const Shape s = shape(frame, fx.vertices[0], fx.vertices[1], fx.vertices[2]);
    REQUIRE(clearly_bad(s, kTheta));
    REQUIRE(s.circumradius < floor_of(frame));
    const QualityOutcome q = run(m, frame, fx);
    REQUIRE(q.inserted == 0);
    REQUIRE(skips(q) == 0);
    REQUIRE(m.triangle_count() == 1);
}
