// Increment 21b (docs/increments/21-parallel-refine.md, section 3 QW2, section
// 7's 21b row, and "Pinned by the red suite (21b)"): an integer incircle for
// quads whose four corners are DEM nodes, answered at the top of must_flip.
// This is 21b's invariant-critical suite.
//
// 21b is bit-identical. Where lattice_incircle answers, its answer must be the
// sign DetriaExact computes on the frame doubles (col * dx, -(row * dy)),
// because that is the number must_flip decides on today. Where any of QW2's
// conditions fails it must refuse, and must_flip must then run today's path.
//
// Interface pinned here (reachable through include/terrain/mesh/lawson.hpp):
//
//   namespace terrain::mesh {
//   // A frame for this grid that also records, once, whether QW2 may answer
//   // on it: dx == dy, finite and > 0, and col * dx, row * dy exact for
//   // every node with col < cols, row < rows.
//   [[nodiscard]] LatticeFrame lattice_frame(double dx, double dy,
//                                            std::size_t rows, std::size_t cols) noexcept;
//   // Precondition: a, b, c strictly counter-clockwise on (col, -row).
//   [[nodiscard]] std::optional<pred::Incircle> lattice_incircle(
//       MeshVertex a, MeshVertex b, MeshVertex c, MeshVertex d, const LatticeFrame& f) noexcept;
//   }
//
// A LatticeFrame built directly, LatticeFrame{dx, dy} as every caller builds
// it today, never enables the integer path.
//
// The oracle is DetriaExact::incircle_ccw on the frame points, computed here
// from col * dx and -(row * dy), not through LatticeFrame::at. For today's
// must_flip the oracle is a copy of it (lawson.hpp at 93a8066).

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/detria_exact.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <span>
#include <string>
#include <utility>
#include <vector>

using terrain::Point2;
using terrain::TriangleIndices;
using terrain::mesh::kNoNeighbour;
using terrain::mesh::lattice_frame;
using terrain::mesh::lattice_incircle;
using terrain::mesh::LatticeFrame;
using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using terrain::mesh::MeshVertex;
using terrain::pred::DefaultKernel;
using terrain::pred::DetriaExact;
using terrain::pred::Incircle;
using terrain::pred::Orientation;

namespace {

constexpr std::int64_t kBound = std::int64_t{1} << 14;  // QW2's spread bound, in nodes
constexpr std::size_t kBench = 5051;                    // the 1 m benchmark's rows and cols

// A node at (col, row). MeshVertex is (col, row); LatticeVertex is {row, col}.
MeshVertex node(std::int64_t col, std::int64_t row) {
    return MeshVertex{static_cast<double>(col), static_cast<double>(row)};
}

// The frame point, written out rather than taken from LatticeFrame::at.
Point2 fp(double dx, double dy, MeshVertex v) { return Point2{v.col * dx, -(v.row * dy)}; }

// Twice the signed area on (col, -row), exact for nodes of the sizes used here.
std::int64_t iorient(MeshVertex a, MeshVertex b, MeshVertex c) {
    const auto ax = static_cast<std::int64_t>(a.col), ay = -static_cast<std::int64_t>(a.row);
    const auto bx = static_cast<std::int64_t>(b.col), by = -static_cast<std::int64_t>(b.row);
    const auto cx = static_cast<std::int64_t>(c.col), cy = -static_cast<std::int64_t>(c.row);
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

// The oracle: DetriaExact on the frame doubles. Empty if a, b, c is not
// counter-clockwise in the frame (DetriaExact's precondition), which under an
// exact frame and a counter-clockwise lattice triangle cannot happen.
std::optional<Incircle> oracle(double dx, double dy, MeshVertex a, MeshVertex b, MeshVertex c,
                               MeshVertex d) {
    const Point2 pa = fp(dx, dy, a), pb = fp(dx, dy, b), pc = fp(dx, dy, c), pd = fp(dx, dy, d);
    if (DetriaExact::orient2d(pa, pb, pc) != Orientation::CounterClockwise)
        return std::nullopt;
    return DetriaExact::incircle_ccw(pa, pb, pc, pd);
}

std::string name(Incircle s) {
    return s == Incircle::Inside ? "Inside" : s == Incircle::Outside ? "Outside" : "Cocircular";
}
std::string describe(MeshVertex v) {
    return "(col " + std::to_string(v.col) + ", row " + std::to_string(v.row) + ")";
}

// Agreement bookkeeping: every quad the caller hands in must be answered, with
// the oracle's sign. The first failure is kept for the report, and the three
// signs are counted so a run that saw only one of them is visible.
struct Agreement {
    std::size_t quads = 0, refused = 0, disagreed = 0;
    std::array<std::size_t, 3> by_sign{};  // Outside, Cocircular, Inside
    std::string first;

    void check(const LatticeFrame& f, double dx, MeshVertex a, MeshVertex b, MeshVertex c,
               MeshVertex d) {
        ++quads;
        const auto want = oracle(dx, dx, a, b, c, d);
        const auto got = lattice_incircle(a, b, c, d, f);
        if (!want || !got || *got != *want) {
            (got ? disagreed : refused) += 1;
            if (first.empty())
                first = describe(a) + " " + describe(b) + " " + describe(c) + " d " + describe(d)
                      + " dx " + std::to_string(dx) + ": got "
                      + (got ? name(*got) : std::string{"nullopt"}) + ", oracle "
                      + (want ? name(*want) : std::string{"not ccw in frame"});
            return;
        }
        ++by_sign[static_cast<std::size_t>(static_cast<int>(*want) + 1)];
    }
    void require_all_agree() const {
        CAPTURE(quads, refused, disagreed, first);
        REQUIRE(quads > 0);
        REQUIRE(refused == 0);
        REQUIRE(disagreed == 0);
    }
};

// The four points in cyclic order, turned counter-clockwise on (col, -row).
std::array<MeshVertex, 4> ccw_cycle(std::array<MeshVertex, 4> p) {
    if (iorient(p[0], p[1], p[2]) < 0)
        std::reverse(p.begin(), p.end());
    return p;
}

// Every lattice point on the circle of squared radius r2 about (cx, cy), in
// angular order.
std::vector<MeshVertex> circle(std::int64_t cx, std::int64_t cy, std::int64_t r2) {
    std::vector<std::pair<double, MeshVertex>> pts;
    const auto root = [](std::int64_t v) {  // floor(sqrt(v)), exact for the sizes here
        auto s = static_cast<std::int64_t>(std::sqrt(static_cast<double>(v)));
        while (s * s > v)
            --s;
        while ((s + 1) * (s + 1) <= v)
            ++s;
        return s;
    };
    const auto r = root(r2);
    for (std::int64_t x = -r; x <= r; ++x) {
        const std::int64_t y = root(r2 - x * x);
        if (x * x + y * y != r2)
            continue;
        for (const std::int64_t sy : {y, -y}) {
            pts.emplace_back(std::atan2(static_cast<double>(sy), static_cast<double>(x)), node(cx + x, cy + sy));
            if (y == 0)
                break;
        }
    }
    std::sort(pts.begin(), pts.end(), [](const auto& l, const auto& r) { return l.first < r.first; });
    std::vector<MeshVertex> out;
    for (const auto& [angle, v] : pts)
        out.push_back(v);
    return out;
}

// Every choice of four points from a cyclically ordered set, as quads in
// cyclic order.
std::vector<std::array<MeshVertex, 4>> quads_of(const std::vector<MeshVertex>& ring) {
    std::vector<std::array<MeshVertex, 4>> out;
    const std::size_t n = ring.size();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j)
            for (std::size_t k = j + 1; k < n; ++k)
                for (std::size_t l = k + 1; l < n; ++l)
                    out.push_back(ccw_cycle({ring[i], ring[j], ring[k], ring[l]}));
    return out;
}

// A cocircular quad in each of its four rotations, (a, b, c) the triangle and
// d the fourth corner, as must_flip hands them over.
template <class Fn>
void for_each_rotation(const std::array<MeshVertex, 4>& q, Fn&& fn) {
    for (std::size_t k = 0; k < 4; ++k)
        fn(q[k], q[(k + 1) % 4], q[(k + 2) % 4], q[(k + 3) % 4]);
}

}  // namespace

// ---------------------------------------------------------------- agreement

TEST_CASE("lattice_incircle agrees with DetriaExact on every quad of a 5 x 5 neighbourhood",
          "[lattice_incircle][agreement]") {
    // Every ordered (a, b, c, d) of distinct nodes with a, b, c strictly
    // counter-clockwise: all ties, all near misses and all clear cases that
    // fit. At the origin and at the far corner of the benchmark's grid, in
    // four exact frames (dx = 10 is the benchmark's).
    const double dx = GENERATE(1.0, 10.0, 0.5, 0.375);
    const std::int64_t base = GENERATE(std::int64_t{0}, std::int64_t{5046});
    CAPTURE(dx, base);
    const LatticeFrame f = lattice_frame(dx, dx, kBench, kBench);
    std::vector<MeshVertex> nodes;
    for (std::int64_t r = 0; r < 5; ++r)
        for (std::int64_t c = 0; c < 5; ++c)
            nodes.push_back(node(base + c, base + r));
    Agreement agree;
    for (const auto& a : nodes)
        for (const auto& b : nodes)
            for (const auto& c : nodes) {
                if (iorient(a, b, c) <= 0)
                    continue;
                for (const auto& d : nodes)
                    if (d != a && d != b && d != c)
                        agree.check(f, dx, a, b, c, d);
            }
    agree.require_all_agree();
    // All three answers occur, ties included, so a constant answer cannot pass.
    REQUIRE(agree.by_sign[0] > 1000);
    REQUIRE(agree.by_sign[1] > 1000);
    REQUIRE(agree.by_sign[2] > 1000);
}

TEST_CASE("lattice_incircle answers Cocircular on the tie shapes 21b-ties found, in every rotation",
          "[lattice_incircle][agreement][ties]") {
    // docs/benchmarks/2026-09-27/21b-ties/README.md, "The shapes of the
    // exact-path quads": axis-aligned squares and rectangles (1x1, 1x2, 2x2,
    // 1x3, 2x3, and legalise_all's 40x40 and 10x40), rotated squares and
    // rectangles, isosceles trapezoids with parallel sides on rows or columns
    // and in another direction, and other cyclic quads. Each is a tie: the
    // answer must be Cocircular, and that must be the oracle's answer too.
    const double dx = GENERATE(1.0, 10.0);
    CAPTURE(dx);
    const LatticeFrame f = lattice_frame(dx, dx, kBench, kBench);
    const auto rect = [](std::int64_t c, std::int64_t r, std::int64_t w, std::int64_t h) {
        return ccw_cycle({node(c, r), node(c + w, r), node(c + w, r + h), node(c, r + h)});
    };
    const auto quad = [](std::array<std::array<std::int64_t, 2>, 4> p) {
        return ccw_cycle({node(p[0][0], p[0][1]), node(p[1][0], p[1][1]), node(p[2][0], p[2][1]),
                          node(p[3][0], p[3][1])});
    };
    const std::vector<std::pair<const char*, std::array<MeshVertex, 4>>> shapes{
        {"axis square 1x1", rect(7, 9, 1, 1)},
        {"axis square 2x2", rect(7, 9, 2, 2)},
        {"axis square 40x40", rect(40, 80, 40, 40)},
        {"axis rectangle 1x2", rect(7, 9, 1, 2)},
        {"axis rectangle 2x1", rect(7, 9, 2, 1)},
        {"axis rectangle 1x3", rect(7, 9, 1, 3)},
        {"axis rectangle 2x3", rect(7, 9, 2, 3)},
        {"axis rectangle 10x40", rect(40, 80, 10, 40)},
        {"rotated square", quad({{{10, 11}, {11, 10}, {12, 11}, {11, 12}}})},
        {"rotated square, side (3, 2)", quad({{{10, 12}, {13, 10}, {15, 13}, {12, 15}}})},
        {"rotated rectangle", quad({{{10, 11}, {11, 10}, {13, 12}, {12, 13}}})},
        {"isosceles trapezoid, parallel sides on rows", quad({{{10, 10}, {14, 10}, {13, 12}, {11, 12}}})},
        {"isosceles trapezoid, parallel sides on columns", quad({{{10, 10}, {12, 11}, {12, 13}, {10, 14}}})},
        // On the radius-5 circle about (20, 20): chords (5,0)-(0,5) and
        // (4,-3)-(-3,4) are both along (-1, 1).
        {"isosceles trapezoid, parallel sides along a diagonal",
         quad({{{25, 20}, {20, 25}, {17, 24}, {24, 17}}})},
        // Radius 5 again: no two sides parallel.
        {"other cyclic quad", quad({{{25, 20}, {23, 24}, {15, 20}, {16, 17}}})},
    };
    for (const auto& [label, q] : shapes) {
        CAPTURE(std::string{label});
        for_each_rotation(q, [&](MeshVertex a, MeshVertex b, MeshVertex c, MeshVertex d) {
            CAPTURE(describe(a), describe(b), describe(c), describe(d));
            REQUIRE(iorient(a, b, c) > 0);
            REQUIRE(oracle(dx, dx, a, b, c, d) == Incircle::Cocircular);
            REQUIRE(lattice_incircle(a, b, c, d, f) == Incircle::Cocircular);
        });
    }
}

TEST_CASE("lattice_incircle agrees with DetriaExact on every quad of lattice circles, and one node off them",
          "[lattice_incircle][agreement][ties]") {
    // Every four of the 12 lattice points on the radius-5 circle and of the 16
    // on the radius-sqrt(65) circle, in every rotation: all ties. Then d moved
    // one node in each axis direction, which puts it strictly inside or outside.
    // The scaled radius-5 circle (x 1617 = 3 * 7^2 * 11, radius 8085, still 12
    // lattice points) has spreads up to 16170 nodes, just under 2^14: there a
    // one-node move of d changes the determinant by about 2^40 in terms near
    // 2^56, beyond what plain double arithmetic on the differences resolves.
    const double dx = GENERATE(1.0, 10.0);
    const std::int64_t r2 = GENERATE(std::int64_t{25}, std::int64_t{65}, std::int64_t{25} * 1617 * 1617);
    CAPTURE(dx, r2);
    const std::int64_t centre = 9000;  // keeps every point at col, row >= 0
    const auto ring = circle(centre, centre, r2);
    REQUIRE(ring.size() == (r2 == 65 ? 16u : 12u));
    const LatticeFrame f = lattice_frame(dx, dx, 20000, 20000);
    Agreement ties, moved;
    std::size_t cocircular = 0;
    for (const auto& q : quads_of(ring))
        for_each_rotation(q, [&](MeshVertex a, MeshVertex b, MeshVertex c, MeshVertex d) {
            ties.check(f, dx, a, b, c, d);
            cocircular += lattice_incircle(a, b, c, d, f) == Incircle::Cocircular;
            for (const auto& [mc, mr] : {std::pair{1, 0}, std::pair{-1, 0}, std::pair{0, 1}, std::pair{0, -1}}) {
                const MeshVertex e{d.col + mc, d.row + mr};
                if (e != a && e != b && e != c)
                    moved.check(f, dx, a, b, c, e);
            }
        });
    ties.require_all_agree();
    moved.require_all_agree();
    REQUIRE(cocircular == ties.quads);
    REQUIRE(ties.by_sign[1] == ties.quads);
    REQUIRE(moved.by_sign[0] > 0);
    REQUIRE(moved.by_sign[2] > 0);
}

TEST_CASE("lattice_incircle agrees with DetriaExact on random quads up to the spread bound",
          "[lattice_incircle][agreement][random]") {
    // Two families, spreads drawn log-uniformly from 1 to 2^14 nodes:
    //   - uniform: a, b, c anywhere within the spread of d;
    //   - near-circle: four points rounded from one circle, so the
    //     determinant is small against its terms, where rounding would show.
    const std::uint32_t seed = GENERATE(range(1u, 5u));
    const double dx = GENERATE(1.0, 10.0, 0.25, 3.0);
    CAPTURE(seed, dx);
    const LatticeFrame f = lattice_frame(dx, dx, 1u << 16, 1u << 16);
    std::mt19937_64 rng{seed};
    std::uniform_real_distribution<double> unit{0.0, 1.0};
    const std::int64_t mid = std::int64_t{1} << 15;
    Agreement uniform, near;
    const auto within = [](MeshVertex p, MeshVertex d) {
        return std::fabs(p.col - d.col) <= static_cast<double>(kBound)
            && std::fabs(p.row - d.row) <= static_cast<double>(kBound);
    };
    for (int i = 0; i < 4000; ++i) {
        const auto s = static_cast<std::int64_t>(std::exp2(14.0 * unit(rng)));
        const auto off = [&] { return static_cast<std::int64_t>(std::llround((2.0 * unit(rng) - 1.0) * s)); };
        // One draw per statement: argument evaluation order is unspecified.
        std::array<std::int64_t, 8> o{};
        for (auto& x : o)
            x = off();
        const MeshVertex d = node(mid + o[0], mid + o[1]);
        MeshVertex a = node(mid + o[0] + o[2], mid + o[1] + o[3]);
        MeshVertex b = node(mid + o[0] + o[4], mid + o[1] + o[5]);
        const MeshVertex c = node(mid + o[0] + o[6], mid + o[1] + o[7]);
        if (iorient(a, b, c) == 0 || d == a || d == b || d == c)
            continue;
        if (iorient(a, b, c) < 0)
            std::swap(a, b);
        uniform.check(f, dx, a, b, c, d);
    }
    for (int i = 0; i < 4000; ++i) {
        const double radius = std::exp2(1.0 + 12.9 * unit(rng));  // diameter <= 2^14
        const double cx = static_cast<double>(mid) + 1000.0 * unit(rng);
        const double cy = static_cast<double>(mid) + 1000.0 * unit(rng);
        std::array<MeshVertex, 4> p{};
        for (auto& v : p) {
            const double t = 2.0 * std::numbers::pi * unit(rng);
            v = MeshVertex{std::round(cx + radius * std::cos(t)), std::round(cy + radius * std::sin(t))};
        }
        auto [a, b, c, d] = p;
        if (iorient(a, b, c) == 0 || d == a || d == b || d == c || !within(a, d) || !within(b, d)
            || !within(c, d))
            continue;
        if (iorient(a, b, c) < 0)
            std::swap(a, b);
        near.check(f, dx, a, b, c, d);
    }
    uniform.require_all_agree();
    near.require_all_agree();
    REQUIRE(uniform.quads > 3000);
    REQUIRE(near.quads > 3000);
    REQUIRE(near.by_sign[0] > 100);
    REQUIRE(near.by_sign[2] > 100);
}

// ---------------------------------------------------------------- the spread bound

namespace {

// a, b, c with one corner `far` from d along one axis, placed at `pos` in the
// counter-clockwise triple; the other two close to d.
std::array<MeshVertex, 3> far_triple(MeshVertex d, std::int64_t dc, std::int64_t dr, unsigned pos) {
    const auto dcol = static_cast<std::int64_t>(d.col), drow = static_cast<std::int64_t>(d.row);
    std::array<MeshVertex, 3> t{node(dcol + dc, drow + dr), node(dcol + 1, drow + 2), node(dcol - 2, drow + 1)};
    if (iorient(t[0], t[1], t[2]) < 0)
        std::swap(t[1], t[2]);
    std::rotate(t.begin(), t.begin() + (3 - pos) % 3, t.end());
    return t;
}

}  // namespace

TEST_CASE("lattice_incircle answers at a spread of exactly 2^14 nodes and refuses at 2^14 + 1",
          "[lattice_incircle][bound]") {
    // QW2: every coordinate difference from d at most 2^14 nodes, in absolute
    // value. For each corner position, each axis and each direction: at 2^14
    // it answers, with the oracle's sign; one node further it refuses.
    const std::int64_t spread = GENERATE(kBound, kBound + 1);
    const unsigned pos = GENERATE(0u, 1u, 2u);
    const auto [dc, dr] = GENERATE(std::pair{1, 0}, std::pair{-1, 0}, std::pair{0, 1}, std::pair{0, -1});
    CAPTURE(spread, pos, dc, dr);
    const double dx = 10.0;
    const LatticeFrame f = lattice_frame(dx, dx, 1u << 16, 1u << 16);
    const MeshVertex d = node(1 << 15, 1 << 15);
    const auto [a, b, c] = far_triple(d, dc * spread, dr * spread, pos);
    REQUIRE(iorient(a, b, c) > 0);
    const auto got = lattice_incircle(a, b, c, d, f);
    if (spread == kBound) {
        REQUIRE(got.has_value());
        REQUIRE(got == oracle(dx, dx, a, b, c, d));
    } else {
        REQUIRE_FALSE(got.has_value());
    }
}

TEST_CASE("lattice_incircle measures the spread from d, not between the other corners",
          "[lattice_incircle][bound]") {
    // a and b are 2^14 from d on either side, so 2^15 apart; the bound is on
    // differences from d (QW2), and this quad is within it.
    const LatticeFrame f = lattice_frame(1.0, 1.0, 1u << 16, 1u << 16);
    const MeshVertex d = node(1 << 15, 1 << 15);
    MeshVertex a = node((1 << 15) + kBound, (1 << 15) + 3), b = node((1 << 15) - kBound, (1 << 15) + 1),
               c = node(1 << 15, (1 << 15) - kBound);
    if (iorient(a, b, c) < 0)
        std::swap(a, b);
    const auto got = lattice_incircle(a, b, c, d, f);
    REQUIRE(got.has_value());
    REQUIRE(got == oracle(1.0, 1.0, a, b, c, d));
}

// ---------------------------------------------------------------- orientation

TEST_CASE("lattice_incircle takes the determinant on (col, -row), as the frame is oriented",
          "[lattice_incircle][orientation]") {
    // In the frame (col, -row): a (0, -2), b (2, -2), c (1, 0) turn
    // counter-clockwise; their circle has centre (1, -1.25) and radius 1.25,
    // and d (1, -1) is 0.25 from the centre, strictly inside. On (col, row)
    // the four points are the mirror image: a, b, c turn clockwise, so the
    // same determinant formula has the opposite sign and says Outside. That is
    // the answer an implementation taking (col, row) gives here.
    const LatticeFrame f = lattice_frame(1.0, 1.0, 16, 16);
    const MeshVertex a = node(0, 2), b = node(2, 2), c = node(1, 0), d = node(1, 1);
    REQUIRE(iorient(a, b, c) > 0);
    REQUIRE(oracle(1.0, 1.0, a, b, c, d) == Incircle::Inside);
    REQUIRE(lattice_incircle(a, b, c, d, f) == Incircle::Inside);
    // And a point clearly outside, which (col, row) would call Inside.
    const MeshVertex e = node(5, 5);
    REQUIRE(oracle(1.0, 1.0, a, b, c, e) == Incircle::Outside);
    REQUIRE(lattice_incircle(a, b, c, e, f) == Incircle::Outside);
}

// ---------------------------------------------------------------- refusals

TEST_CASE("lattice_incircle refuses a quad with an off-node corner, in each position",
          "[lattice_incircle][refuse]") {
    // QW2 answers only when all four corners are nodes (MeshVertex::is_node).
    // Off-node corners are the domain's arc vertices (16 R2) and 20b's feet.
    const LatticeFrame f = lattice_frame(10.0, 10.0, kBench, kBench);
    std::array<MeshVertex, 4> q{node(10, 12), node(12, 12), node(11, 10), node(11, 11)};
    REQUIRE(iorient(q[0], q[1], q[2]) > 0);
    REQUIRE(lattice_incircle(q[0], q[1], q[2], q[3], f).has_value());  // all nodes: answers
    const std::size_t which = GENERATE(0u, 1u, 2u, 3u);
    const auto [dc, dr] = GENERATE(std::pair{0.5, 0.0}, std::pair{0.0, 0.25},
                                   std::pair{std::ldexp(1.0, -40), 0.0}, std::pair{0.0, std::ldexp(1.0, -40)});
    CAPTURE(which, dc, dr);
    q[which] = MeshVertex{q[which].col + dc, q[which].row + dr};
    REQUIRE_FALSE(q[which].is_node());
    REQUIRE(terrain::mesh::orient_sign(q[0], q[1], q[2]) > 0);  // still counter-clockwise
    REQUIRE_FALSE(lattice_incircle(q[0], q[1], q[2], q[3], f).has_value());
}

TEST_CASE("lattice_incircle refuses every frame QW2 excludes", "[lattice_incircle][refuse]") {
    // The quad is a plain all-node quad the benchmark frame answers; only the
    // frame changes.
    const MeshVertex a = node(3, 7), b = node(6, 7), c = node(4, 3), d = node(4, 5);
    REQUIRE(iorient(a, b, c) > 0);
    REQUIRE(lattice_incircle(a, b, c, d, lattice_frame(10.0, 10.0, 100, 100)).has_value());

    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    SECTION("dx != dy: the frame is not the lattice times one constant") {
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, lattice_frame(1.0, 2.0, 100, 100)).has_value());
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, lattice_frame(10.0, 5.0, 100, 100)).has_value());
        REQUIRE_FALSE(
            lattice_incircle(a, b, c, d, lattice_frame(1.0, std::nextafter(1.0, 2.0), 100, 100)).has_value());
    }
    SECTION("dx not finite and positive") {
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, lattice_frame(0.0, 0.0, 100, 100)).has_value());
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, lattice_frame(-1.0, -1.0, 100, 100)).has_value());
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, lattice_frame(inf, inf, 100, 100)).has_value());
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, lattice_frame(nan, nan, 100, 100)).has_value());
    }
    SECTION("an inexact frame: some col * dx or row * dy of the grid rounds") {
        // 3 * 0.1 != 0.3 exactly; the grid holds col and row 3.
        REQUIRE(std::fma(3.0, 0.1, -(3.0 * 0.1)) != 0.0);
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, lattice_frame(0.1, 0.1, 100, 100)).has_value());
        // 0.1 has 52 significant bits, so a grid two nodes wide is exact on
        // that axis; the other axis is long enough to round. Rows and cols
        // must both count. The quads lie inside those grids.
        const MeshVertex ta = node(0, 5), tb = node(1, 5), tc = node(0, 3), td = node(1, 4);  // cols 0..1
        const MeshVertex wa = node(3, 1), wb = node(6, 1), wc = node(4, 0), wd = node(5, 0);  // rows 0..1
        REQUIRE(iorient(ta, tb, tc) > 0);
        REQUIRE(iorient(wa, wb, wc) > 0);
        REQUIRE_FALSE(lattice_incircle(ta, tb, tc, td, lattice_frame(0.1, 0.1, 100, 2)).has_value());
        REQUIRE_FALSE(lattice_incircle(wa, wb, wc, wd, lattice_frame(0.1, 0.1, 2, 100)).has_value());
        // A 53-bit dx: 3 * dx needs 54 bits.
        const double full = 1.0 + std::ldexp(1.0, -52);
        REQUIRE(std::fma(3.0, full, -(3.0 * full)) != 0.0);
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, lattice_frame(full, full, 100, 100)).has_value());
    }
    SECTION("a LatticeFrame built directly never enables the integer path") {
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, LatticeFrame{10.0, 10.0}).has_value());
        REQUIRE_FALSE(lattice_incircle(a, b, c, d, LatticeFrame{1.0, 1.0}).has_value());
    }
}

TEST_CASE("lattice_frame enables the integer path wherever QW2's sufficient condition holds",
          "[lattice_incircle][frame]") {
    // QW2: the frame is exact when the significant bits of dx plus
    // bit_width(max(rows, cols) - 1) are at most 53. Inside that, it must
    // answer (that is the speed-up, and must_flip's zero-kernel-call test
    // below depends on it).
    const MeshVertex a = node(0, 1), b = node(1, 1), c = node(0, 0), d = node(1, 0);  // unit square
    REQUIRE(iorient(a, b, c) > 0);
    // The benchmark: 3 + 13 = 16.
    REQUIRE(lattice_incircle(a, b, c, d, lattice_frame(10.0, 10.0, kBench, kBench)) == Incircle::Cocircular);
    // Exactly 53: a 52-bit dx on a 2 x 2 grid (bit_width(1) = 1).
    const double dx52 = 1.0 + std::ldexp(1.0, -51);
    REQUIRE(lattice_incircle(a, b, c, d, lattice_frame(dx52, dx52, 2, 2)) == Incircle::Cocircular);
    // A power of two on a 2^20 grid: 1 + 20.
    const double tiny = std::ldexp(1.0, -30);
    REQUIRE(lattice_incircle(a, b, c, d, lattice_frame(tiny, tiny, 1u << 20, 1u << 20)) == Incircle::Cocircular);
    // 0.1 on a 2 x 2 grid (52 + 1 = 53): col and row are 0 or 1, exact.
    REQUIRE(lattice_incircle(a, b, c, d, lattice_frame(0.1, 0.1, 2, 2)) == Incircle::Cocircular);
}

TEST_CASE("an inexact frame changes the answer on some lattice tie, so its refusal is load-bearing",
          "[lattice_incircle][refuse]") {
    // QW2's determinism argument: without the exact-frame condition the
    // integer answer is the true lattice answer, and the rounded frame's
    // answer can differ on a tie. Find one: a cocircular lattice quad on which
    // DetriaExact at dx = dy = 0.1 does not say Cocircular. lattice_incircle
    // must refuse it (it would otherwise say Cocircular and change a flip).
    const double dx = 0.1;
    const LatticeFrame f = lattice_frame(dx, dx, 100, 100);
    std::optional<std::array<MeshVertex, 4>> found;
    for (const std::int64_t r2 : {25, 50, 65, 85, 125})
        for (const auto& q : quads_of(circle(40, 40, r2)))
            for_each_rotation(q, [&](MeshVertex a, MeshVertex b, MeshVertex c, MeshVertex d) {
                const auto o = oracle(dx, dx, a, b, c, d);
                if (!found && o && *o != Incircle::Cocircular)
                    found = std::array<MeshVertex, 4>{a, b, c, d};
            });
    REQUIRE(found.has_value());
    const auto [a, b, c, d] = *found;
    CAPTURE(describe(a), describe(b), describe(c), describe(d));
    REQUIRE_FALSE(lattice_incircle(a, b, c, d, f).has_value());
}

// ---------------------------------------------------------------- must_flip

namespace {

// must_flip as of 93a8066, the oracle for decisions.
bool reference_must_flip(const LatticeMesh& m, std::uint32_t t, unsigned e, const LatticeFrame& f) {
    const auto u = m.neighbours(t)[e];
    if (u == kNoNeighbour || m.is_constrained(t, e))
        return false;
    const auto& tri = m.triangles()[t];
    unsigned j = 0;
    while (m.triangles()[u][j] != tri[(e + 1) % 3])
        ++j;
    const auto v = m.vertices();
    const Point2 a = fp(f.dx, f.dy, v[tri[e]]), b = fp(f.dx, f.dy, v[tri[(e + 1) % 3]]),
                 c = fp(f.dx, f.dy, v[tri[(e + 2) % 3]]), d = fp(f.dx, f.dy, v[m.triangles()[u][(j + 2) % 3]]);
    if (DefaultKernel::orient2d(a, b, c) == Orientation::CounterClockwise)
        return DefaultKernel::incircle(a, b, c, d) == Incircle::Inside;
    return DefaultKernel::orient2d(b, a, d) == Orientation::CounterClockwise
        && DefaultKernel::incircle(b, a, d, c) == Incircle::Inside;
}

// DefaultKernel, counting its calls: must_flip asks it nothing when the
// integer path answers.
struct CountingKernel {
    static inline std::size_t calls = 0;
    static Orientation orient2d(const Point2& a, const Point2& b, const Point2& c) {
        ++calls;
        return DefaultKernel::orient2d(a, b, c);
    }
    static Incircle incircle(const Point2& a, const Point2& b, const Point2& c, const Point2& d) {
        ++calls;
        return DefaultKernel::incircle(a, b, c, d);
    }
};
static_assert(terrain::pred::GeometryKernel<CountingKernel>);

// An (n+1) x (n+1) grid of nodes `step` apart from (base, base), two triangles
// per cell, the ring constrained. Mixed diagonals give ties both ways.
LatticeMesh grid(std::uint32_t n, std::uint32_t step, std::uint32_t base) {
    std::vector<LatticeVertex> v;
    for (std::uint32_t r = 0; r <= n; ++r)
        for (std::uint32_t c = 0; c <= n; ++c)
            v.push_back(LatticeVertex{base + r * step, base + c * step});
    const auto at = [n](std::uint32_t r, std::uint32_t c) { return r * (n + 1) + c; };
    std::vector<TriangleIndices> t;
    for (std::uint32_t r = 0; r < n; ++r)
        for (std::uint32_t c = 0; c < n; ++c) {
            const auto tl = at(r, c), tr = at(r, c + 1), bl = at(r + 1, c), br = at(r + 1, c + 1);
            if ((r + c) % 2 == 0) {
                t.push_back({tl, bl, br});
                t.push_back({tl, br, tr});
            } else {
                t.push_back({tl, bl, tr});
                t.push_back({bl, br, tr});
            }
        }
    std::vector<std::uint8_t> bits;
    std::vector<std::array<std::uint32_t, 3>> masks;
    const std::uint32_t lo = base, hi = base + n * step;
    for (const auto& tri : t) {
        std::uint8_t b = 0;
        for (unsigned k = 0; k < 3; ++k) {
            const auto p = v[tri[k]], q = v[tri[(k + 1) % 3]];
            if ((p.row == q.row && (p.row == lo || p.row == hi)) || (p.col == q.col && (p.col == lo || p.col == hi)))
                b = static_cast<std::uint8_t>(b | (1u << k));
        }
        bits.push_back(b);
        masks.push_back({b & 1u, (b >> 1) & 1u, (b >> 2) & 1u});
    }
    auto m = LatticeMesh::build(std::move(v), std::move(t), std::move(bits), std::move(masks));
    REQUIRE(m.has_value());
    return *m;
}

// Insert node p strictly inside a triangle or on an edge; false when p is a
// vertex. Exact orientation on (col, -row), so off-node vertices are fine.
bool insert_node(LatticeMesh& m, LatticeVertex p) {
    const MeshVertex x{p};
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        const auto& tri = m.triangles()[t];
        const auto v = m.vertices();
        std::array<int, 3> o{};
        for (unsigned k = 0; k < 3; ++k)
            o[k] = terrain::mesh::orient_sign(v[tri[k]], v[tri[(k + 1) % 3]], x);
        if (o[0] < 0 || o[1] < 0 || o[2] < 0)
            continue;
        const int zeros = (o[0] == 0) + (o[1] == 0) + (o[2] == 0);
        if (zeros >= 2)
            return false;
        if (zeros == 0)
            m.split_inside(t, p);
        else
            m.split_edge(t, o[0] == 0 ? 0u : o[1] == 0 ? 1u : 2u, x);
        return true;
    }
    return false;
}

// Split a constrained ring edge at an off-node point strictly inside it (a
// foot, as 20b inserts them). false if no ring edge is long enough.
bool insert_foot(LatticeMesh& m, std::mt19937& rng) {
    const auto n = m.triangle_count();
    const auto start = static_cast<std::uint32_t>(rng() % n);
    for (std::uint32_t i = 0; i < n; ++i) {
        const auto t = (start + i) % n;
        for (unsigned k = 0; k < 3; ++k) {
            if (!m.is_constrained(t, k))
                continue;
            const auto p = m.vertices()[m.triangles()[t][k]], q = m.vertices()[m.triangles()[t][(k + 1) % 3]];
            const double len = std::fabs(p.col - q.col) + std::fabs(p.row - q.row);  // axis-aligned
            if (len < 2.0)
                continue;
            const double s = (std::floor(len / 2.0) + 0.5) / len;  // off-node, strictly inside
            m.split_edge(t, k, MeshVertex{p.col + s * (q.col - p.col), p.row + s * (q.row - p.row)});
            return true;
        }
    }
    return false;
}

struct Decisions {
    std::size_t edges = 0, flips = 0, ties = 0, differ = 0;
    std::string first;
};

// Every (t, e): must_flip under the lattice frame against today's must_flip.
void compare_decisions(const LatticeMesh& m, double dx, double dy, std::size_t grid_n, Decisions& out) {
    const LatticeFrame with = lattice_frame(dx, dy, grid_n, grid_n);
    const LatticeFrame without{dx, dy};
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned e = 0; e < 3; ++e) {
            if (m.neighbours(t)[e] == kNoNeighbour || m.is_constrained(t, e))
                continue;
            ++out.edges;
            const bool want = reference_must_flip(m, t, e, without);
            const bool got = terrain::mesh::detail::must_flip<DefaultKernel>(m, t, e, with);
            out.flips += want;
            const auto& tri = m.triangles()[t];
            const auto u = m.neighbours(t)[e];
            unsigned j = 0;
            while (m.triangles()[u][j] != tri[(e + 1) % 3])
                ++j;
            const auto v = m.vertices();
            out.ties += DefaultKernel::incircle(fp(dx, dy, v[tri[0]]), fp(dx, dy, v[tri[1]]), fp(dx, dy, v[tri[2]]),
                                                fp(dx, dy, v[m.triangles()[u][(j + 2) % 3]]))
                     == Incircle::Cocircular;
            if (got != want && out.differ++ == 0)
                out.first = "t " + std::to_string(t) + " e " + std::to_string(e) + ": got "
                          + std::to_string(got) + ", today " + std::to_string(want);
        }
}

struct FrameCase {
    double dx, dy;
    bool integer;  // whether lattice_frame enables the integer path on it
};
const std::array<FrameCase, 5> kFrameCases{FrameCase{1.0, 1.0, true}, FrameCase{10.0, 10.0, true},
                                           FrameCase{0.5, 0.5, true}, FrameCase{0.1, 0.1, false},
                                           FrameCase{10.0, 5.0, false}};

}  // namespace

TEST_CASE("must_flip decides as today on random lattice meshes, with ties and off-node feet",
          "[lattice_incircle][must_flip]") {
    // Random nodes inserted into a mixed-diagonal grid without legalising, so
    // many edges must flip and many are ties; every few insertions a foot
    // splits a ring edge off the lattice. After every insertion, every
    // interior unconstrained edge is decided by must_flip under
    // lattice_frame(dx, dy, ...) and by a copy of today's must_flip under
    // LatticeFrame{dx, dy}: the decisions must be equal. Near the origin and
    // near the benchmark's far corner.
    const std::uint32_t seed = GENERATE(range(1u, 7u));
    const std::size_t frame = GENERATE(range(std::size_t{0}, kFrameCases.size()));
    const std::uint32_t base = GENERATE(0u, 5000u);
    const auto [dx, dy, integer] = kFrameCases[frame];
    CAPTURE(seed, dx, dy, base);
    auto m = grid(4, 6, base);  // nodes base .. base + 24
    std::mt19937 rng{seed};
    Decisions dec;
    std::size_t inserted = 0, feet = 0;
    for (int attempt = 0; attempt < 120; ++attempt) {
        if (attempt % 15 == 7 && insert_foot(m, rng))
            ++feet;
        const LatticeVertex p{base + static_cast<std::uint32_t>(rng() % 25),
                              base + static_cast<std::uint32_t>(rng() % 25)};
        if (!insert_node(m, p))
            continue;
        ++inserted;
        compare_decisions(m, dx, dy, kBench, dec);
    }
    CAPTURE(dec.edges, dec.flips, dec.ties, dec.first);
    REQUIRE(dec.differ == 0);
    REQUIRE(inserted >= 40);
    REQUIRE(feet >= 3);
    REQUIRE(dec.flips >= 100);
    REQUIRE(dec.ties >= 100);
    REQUIRE(dec.edges - dec.flips - dec.ties >= 100);
}

TEST_CASE("legalise_all under lattice_frame makes the same flips, writes and mesh as under LatticeFrame{dx, dy}",
          "[lattice_incircle][must_flip]") {
    // The mesh-level form of bit-identity: the whole Lawson run, not one
    // decision at a time. LatticeFrame{dx, dy} never enables the integer path,
    // so it is today's run.
    const std::uint32_t seed = GENERATE(range(1u, 7u));
    const std::size_t frame = GENERATE(range(std::size_t{0}, kFrameCases.size()));
    const auto [dx, dy, integer] = kFrameCases[frame];
    CAPTURE(seed, dx, dy);
    auto m = grid(5, 5, 0);
    std::mt19937 rng{seed};
    for (int attempt = 0; attempt < 150; ++attempt) {
        if (attempt % 20 == 3)
            insert_foot(m, rng);
        insert_node(m, LatticeVertex{static_cast<std::uint32_t>(rng() % 26), static_cast<std::uint32_t>(rng() % 26)});
    }
    LatticeMesh today = m;
    std::vector<std::uint32_t> w_today, w_new;
    const auto f_today = terrain::mesh::legalise_all<DefaultKernel>(
        today, LatticeFrame{dx, dy}, [&](std::uint32_t s) { w_today.push_back(s); });
    const auto f_new = terrain::mesh::legalise_all<DefaultKernel>(
        m, lattice_frame(dx, dy, 26, 26), [&](std::uint32_t s) { w_new.push_back(s); });
    REQUIRE(f_today >= 20);
    REQUIRE(f_new == f_today);
    REQUIRE(w_new == w_today);
    REQUIRE(m.triangle_count() == today.triangle_count());
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        CAPTURE(t);
        REQUIRE(m.triangles()[t] == today.triangles()[t]);
        REQUIRE(m.neighbours(t) == today.neighbours(t));
    }
}

TEST_CASE("must_flip asks the kernel nothing when lattice_incircle answers, and today's questions when it refuses",
          "[lattice_incircle][must_flip]") {
    // QW2: lattice_incircle is called at the top of must_flip, and when it
    // answers, must_flip also skips the frame orient2d (the mesh triangle is
    // counter-clockwise in integers, and so in an exact frame). That is the
    // whole speed-up, so it is pinned: zero kernel calls on an all-node quad
    // under an enabling frame; kernel calls under a direct LatticeFrame, an
    // excluded frame, or with an off-node corner.
    auto m = grid(4, 3, 0);
    std::mt19937 rng{7};
    for (int i = 0; i < 40; ++i)
        insert_node(m, LatticeVertex{static_cast<std::uint32_t>(rng() % 13), static_cast<std::uint32_t>(rng() % 13)});
    const auto run = [&m](const LatticeFrame& f, bool node_quads_only, std::size_t min_edges) {
        std::size_t asked = 0, edges = 0;
        for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
            for (unsigned e = 0; e < 3; ++e) {
                const auto u = m.neighbours(t)[e];
                if (u == kNoNeighbour || m.is_constrained(t, e))
                    continue;
                bool nodes = true;
                for (const auto i : m.triangles()[t])
                    nodes = nodes && m.vertices()[i].is_node();
                for (const auto i : m.triangles()[u])
                    nodes = nodes && m.vertices()[i].is_node();
                if (nodes != node_quads_only)
                    continue;
                ++edges;
                CountingKernel::calls = 0;
                const bool got = terrain::mesh::detail::must_flip<CountingKernel>(m, t, e, f);
                REQUIRE(got == reference_must_flip(m, t, e, LatticeFrame{f.dx, f.dy}));
                asked += CountingKernel::calls;
            }
        REQUIRE(edges >= min_edges);
        return asked;
    };
    REQUIRE(run(lattice_frame(10.0, 10.0, kBench, kBench), true, 20) == 0);
    REQUIRE(run(LatticeFrame{10.0, 10.0}, true, 20) > 0);
    REQUIRE(run(lattice_frame(0.1, 0.1, 100, 100), true, 20) > 0);
    REQUIRE(run(lattice_frame(10.0, 5.0, kBench, kBench), true, 20) > 0);
    for (int i = 0; i < 6; ++i)
        REQUIRE(insert_foot(m, rng));
    REQUIRE(run(lattice_frame(10.0, 10.0, kBench, kBench), false, 3) > 0);
}
