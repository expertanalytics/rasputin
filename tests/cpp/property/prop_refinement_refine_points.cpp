// Increment 15c-1 (docs/increments/15c-geographic-dem.md, J2, J3, J8, J9 and
// "Tests for @tester", RP3 and RP4): refine_points against scattered check
// points over a rough function, through the C++ API. RP3's NumPy twin, the
// one the record names, is tests/python/test_core_refine_points.py; this file
// adds the constrained-Delaunay oracle that the QA rules (section D) require
// of every refinement property test, and RP4's thread and add-order sweep.
// Invariant-critical (the record: RP3, and RP2, RP4, RP5 for scan_points).
//
// Both oracles are independent of check_points.hpp and refine_points.hpp:
//
//   tolerance oracle (J2): every given check point that is not a start vertex,
//     in every output triangle whose CLOSED area holds it (exact orientation
//     on (col, -row), recovered from the output's world points), is within
//     tolerance of that triangle's plane, the plane recomputed here from the
//     output's z. It walks the given points, never the store, and reads no
//     scan record. It returns its violations, so it can be shown to fail: the
//     control plants the output's z shifted by twice the tolerance.
//   Delaunay oracle: every interior edge that is not a constraint edge has
//     neither apex strictly inside the other triangle's circumcircle, by the
//     exact incircle in the producer's frame, lattice_frame(h, h, ...), i.e.
//     (col h, -(row h)).
//
// Positions are dyadic (k / 1024 of a cell) on an integral frame and z is a
// float of few bits, so each given point is its stored position bit for bit
// (D4) and the oracle may use the given one.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/refine_points.hpp>

#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <set>
#include <span>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::pred::DefaultKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::raster::RasterGeometry;
using terrain::refinement::CheckPoints;
using terrain::refinement::PointRefineOptions;
using terrain::refinement::PointRefineOutcome;
using terrain::refinement::refine_points;

namespace {

constexpr std::size_t kN = 33;  // nodes a side; 32 x 32 cells
constexpr double kH = 30.0, kX = 1000.0, kY = 2000.0;

RasterGeometry grid() { return RasterGeometry{kX, kY, kH, kH, kN, kN}; }
Point2 world(double col, double row) { return Point2{kX + col * kH, kY - row * kH}; }
Point2 frame(Point2 w) { return Point2{(w.x - kX) / kH, -((kY - w.y) / kH)}; }  // (col, -row)

double surface(double c, double r) { return 50.0 + 20.0 * std::sin(0.3 * c) + 15.0 * std::cos(0.25 * r); }

struct Start {
    IndexedMesh2 mesh;
    std::vector<double> z;
    std::vector<std::uint8_t> valid;
    std::vector<std::array<std::uint32_t, 2>> edges;
    std::vector<std::uint32_t> masks;
};

// A coarse grid mesh standing in for phase 1's output: z from the smooth
// surface at its nodes, the outer ring constrained with masks 1, 2, 4, 8.
Start phase1() {
    auto s = refinement_fixtures::grid_mesh(grid(), 8);
    Start out{std::move(s.mesh), {}, {}, std::move(s.edges), std::move(s.masks)};
    for (const auto& rc : s.lattice) {
        out.z.push_back(surface(static_cast<double>(rc.col), static_cast<double>(rc.row)));
        out.valid.push_back(1);
    }
    return out;
}

struct Points {
    std::vector<Point2> xy;
    std::vector<float> z;
};

// Two points per cell at seeded dyadic offsets in (0, 1), so none is a node; z the surface plus noise of
// up to +-4 m, rounded to 1/64 m so it is a float exactly.
Points scattered(std::uint32_t seed) {
    std::mt19937 gen{seed};
    Points p;
    for (std::size_t r = 0; r + 1 < kN; ++r)
        for (std::size_t c = 0; c + 1 < kN; ++c)
            for (int k = 0; k < 2; ++k) {
                const double col = static_cast<double>(c) + static_cast<double>(1u + gen() % 1023u) / 1024.0;
                const double row = static_cast<double>(r) + static_cast<double>(1u + gen() % 1023u) / 1024.0;
                const double noise = (static_cast<double>(gen() % 513u) - 256.0) / 64.0;
                p.xy.push_back(world(col, row));
                p.z.push_back(static_cast<float>(std::round((surface(col, row) + noise) * 64.0) / 64.0));
            }
    return p;
}

CheckPoints store(const Points& p, bool reversed_in_chunks = false) {
    CheckPoints cp{grid()};
    if (!reversed_in_chunks) {
        cp.add(std::span<const Point2>{p.xy}, std::span<const float>{p.z});
    } else {
        std::vector<Point2> xy(p.xy.rbegin(), p.xy.rend());
        std::vector<float> z(p.z.rbegin(), p.z.rend());
        const std::size_t cut = xy.size() / 3;
        cp.add(std::span<const Point2>{xy}.subspan(cut), std::span<const float>{z}.subspan(cut));
        cp.add(std::span<const Point2>{xy}.first(cut), std::span<const float>{z}.first(cut));
    }
    cp.freeze();
    return cp;
}

PointRefineOutcome run(const CheckPoints& cp, const Start& s, double tol, unsigned threads) {
    PointRefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    return refine_points(cp, s.mesh, std::span<const double>{s.z}, std::span<const std::uint8_t>{s.valid},
                         std::span<const std::array<std::uint32_t, 2>>{s.edges},
                         std::span<const std::uint32_t>{s.masks}, o);
}

double cross(Point2 a, Point2 b, Point2 c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); }

// J2's oracle: the number of (point, triangle) pairs over tolerance, with z
// taken from `z` (the output's, or a planted copy).
std::size_t violations(const Points& pts, const Start& start, const PointRefineOutcome& out,
                       const std::vector<double>& z, double tol) {
    std::set<std::pair<double, double>> start_xy;
    for (const Point2 v : start.mesh.vertices()) start_xy.insert({v.x, v.y});
    std::vector<Point2> fp;
    double zmax = 1.0;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        fp.push_back(frame(out.vertices[i]));
        zmax = std::max(zmax, std::abs(z[i]));
    }
    std::size_t bad = 0;
    for (const auto& t : out.triangles) {
        const Point2 a = fp[t[0]], b = fp[t[1]], c = fp[t[2]];
        REQUIRE(DefaultKernel::orient2d(a, b, c) == Orientation::CounterClockwise);
        if (!(out.valid[t[0]] && out.valid[t[1]] && out.valid[t[2]])) continue;
        const double lo_x = std::min({a.x, b.x, c.x}), hi_x = std::max({a.x, b.x, c.x});
        const double lo_y = std::min({a.y, b.y, c.y}), hi_y = std::max({a.y, b.y, c.y});
        const double two_a = cross(a, b, c);
        for (std::size_t i = 0; i < pts.xy.size(); ++i) {
            if (start_xy.contains({pts.xy[i].x, pts.xy[i].y})) continue;  // J2 excludes them
            const Point2 p = frame(pts.xy[i]);
            if (p.x < lo_x || p.x > hi_x || p.y < lo_y || p.y > hi_y) continue;
            if (DefaultKernel::orient2d(a, b, p) == Orientation::Clockwise
                || DefaultKernel::orient2d(b, c, p) == Orientation::Clockwise
                || DefaultKernel::orient2d(c, a, p) == Orientation::Clockwise)
                continue;
            const double plane = (cross(p, b, c) * z[t[0]] + cross(a, p, c) * z[t[1]] + cross(a, b, p) * z[t[2]]) / two_a;
            if (std::abs(plane - static_cast<double>(pts.z[i])) > tol + 1e-9 * zmax) ++bad;
        }
    }
    return bad;
}

void delaunay_oracle(const PointRefineOutcome& out) {
    auto lf = [&](std::uint32_t i) {
        const Point2 f = frame(out.vertices[i]);
        return Point2{f.x * kH, f.y * kH};
    };
    std::set<std::pair<std::uint32_t, std::uint32_t>> constrained;
    for (const auto& e : out.edges) constrained.insert(std::minmax(e[0], e[1]));
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> sides;
    for (std::size_t t = 0; t < out.triangles.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            sides[std::minmax(out.triangles[t][k], out.triangles[t][(k + 1) % 3])].push_back(t);
    std::size_t bad = 0;
    for (const auto& [e, ts] : sides) {
        REQUIRE(ts.size() <= 2);
        if (ts.size() != 2 || constrained.contains(e)) continue;
        for (unsigned s = 0; s < 2; ++s) {
            const auto& tri = out.triangles[ts[s]];
            std::uint32_t apex = 0;
            for (const auto x : out.triangles[ts[1 - s]])
                if (x != e.first && x != e.second) apex = x;
            if (DefaultKernel::incircle(lf(tri[0]), lf(tri[1]), lf(tri[2]), lf(apex)) == Incircle::Inside) {
                UNSCOPED_INFO("edge " << e.first << "-" << e.second << " apex " << apex);
                ++bad;
            }
        }
    }
    REQUIRE(bad == 0);
}

// J9: every output constraint edge lies on one side of the outer ring and
// carries that side's mask; per side the pieces add up to the side.
void constraints_oracle(const PointRefineOutcome& out) {
    const double last = static_cast<double>(kN - 1);
    std::map<std::uint32_t, double> length;
    for (std::size_t k = 0; k < out.edges.size(); ++k) {
        const Point2 p = frame(out.vertices[out.edges[k][0]]), q = frame(out.vertices[out.edges[k][1]]);
        CAPTURE(k, p.x, p.y, q.x, q.y);
        std::uint32_t side = 0;
        if (p.y == 0.0 && q.y == 0.0) side = 1;
        else if (p.x == last && q.x == last) side = 2;
        else if (p.y == -last && q.y == -last) side = 4;
        else if (p.x == 0.0 && q.x == 0.0) side = 8;
        REQUIRE(side != 0);
        REQUIRE(out.masks[k] == side);
        length[side] += std::hypot(q.x - p.x, q.y - p.y);
    }
    for (const std::uint32_t side : {1u, 2u, 4u, 8u}) REQUIRE(length[side] == last);
}

// J8: every inserted vertex is a given check point, carrying its own z.
void inserted_are_check_points(const Points& pts, const Start& start, const PointRefineOutcome& out) {
    std::map<std::pair<double, double>, float> given;
    for (std::size_t i = 0; i < pts.xy.size(); ++i) given[{pts.xy[i].x, pts.xy[i].y}] = pts.z[i];
    const std::size_t n0 = start.mesh.vertices().size();
    REQUIRE(out.vertices.size() == n0 + out.inserted);
    for (std::size_t i = n0; i < out.vertices.size(); ++i) {
        const auto it = given.find({out.vertices[i].x, out.vertices[i].y});
        REQUIRE(it != given.end());
        REQUIRE(out.z[i] == static_cast<double>(it->second));
        REQUIRE(out.valid[i] == 1);
    }
}

}  // namespace

TEST_CASE("RP3: J2 holds at every check point, and the mesh is constrained Delaunay",
          "[refine_points][RP3][property]") {
    const double tol = GENERATE(0.0, 0.5, 2.0, 8.0);
    const std::uint32_t seed = GENERATE(1u, 2u);
    CAPTURE(tol, seed);
    const Start s = phase1();
    const Points pts = scattered(seed);
    const auto out = run(store(pts), s, tol, 0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted > 0);
    REQUIRE(out.max_error <= tol);
    REQUIRE(out.uncovered == 0);
    REQUIRE(out.coincident == 0);  // no offset is 0, so no point is a node
    REQUIRE(violations(pts, s, out, out.z, tol) == 0);
    delaunay_oracle(out);
    constraints_oracle(out);
    inserted_are_check_points(pts, s, out);
}

TEST_CASE("RP3: the oracle fails on a mesh whose z is shifted by twice the tolerance",
          "[refine_points][RP3][control]") {
    const double tol = GENERATE(0.5, 2.0, 8.0);
    CAPTURE(tol);
    const Start s = phase1();
    const Points pts = scattered(1);
    const auto out = run(store(pts), s, tol, 0);
    REQUIRE(out.ok());
    std::vector<double> planted = out.z;
    for (auto& z : planted) z += 2.0 * tol;
    REQUIRE(violations(pts, s, out, planted, tol) > 0);
}

TEST_CASE("RP5: a split shared edge skips the neighbour's stale result",
          "[refine_points][RP5]") {
    // Each 4 x 4 square is A = (tl, bl, br) and, after it in index order,
    // B = (tl, br, tr), sharing the diagonal tl-br. A's worst point p is the
    // square's centre, on the diagonal; B's worst is q, inside B but in the
    // half (p, br, tr) that split_edge(A, p) appends, not the half it leaves
    // in B's slot. In round 1 A splits the shared edge; B must then be skipped
    // as touched, since its scan names q in a triangle its slot no longer
    // holds. No flip follows p's insertion (every circle through p and an
    // outer side has that side as diameter), so only the skip protects B.
    auto g = refinement_fixtures::grid_mesh(grid(), 4);
    Start s{std::move(g.mesh), std::vector<double>(g.lattice.size(), 0.0),
            std::vector<std::uint8_t>(g.lattice.size(), 1), std::move(g.edges), std::move(g.masks)};
    Points pts;
    for (std::size_t r = 0; r + 4 < kN; r += 4)
        for (std::size_t c = 0; c + 4 < kN; c += 4) {
            const auto rr = static_cast<double>(r), cc = static_cast<double>(c);
            pts.xy.push_back(world(cc + 2.0, rr + 2.0));  // p, |error| 5
            pts.z.push_back(5.0f);
            pts.xy.push_back(world(cc + 3.5, rr + 2.25));  // q, |error| 9
            pts.z.push_back(9.0f);
        }
    const auto out = run(store(pts), s, 1.0, 1);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == pts.xy.size());
    REQUIRE(out.triangles.size() == s.mesh.triangle_count() + 2 * pts.xy.size());
    double area = 0.0;
    for (const auto& t : out.triangles) {
        const Point2 a = frame(out.vertices[t[0]]), b = frame(out.vertices[t[1]]), c = frame(out.vertices[t[2]]);
        REQUIRE(DefaultKernel::orient2d(a, b, c) == Orientation::CounterClockwise);
        area += cross(a, b, c) / 2.0;
    }
    REQUIRE(area == static_cast<double>((kN - 1) * (kN - 1)));
    REQUIRE(violations(pts, s, out, out.z, 1.0) == 0);
    delaunay_oracle(out);
    inserted_are_check_points(pts, s, out);
}

TEST_CASE("RP5: interior edges split with work on both sides in one round stay a triangulation",
          "[refine_points][RP5][property]") {
    // Every 4 x 4 square of the start is split by its diagonal, an interior
    // edge shared by two triangles. Check points sit exactly on every diagonal
    // (and two inside each square), so in the first round both triangles of a
    // square name a point on their shared edge: the first to split must mark
    // its neighbour touched, or the neighbour splits a triangle that no longer
    // exists. threads 1, so the order is the serial one.
    auto g = refinement_fixtures::grid_mesh(grid(), 4);
    Start s{std::move(g.mesh), std::vector<double>(g.lattice.size(), 0.0),
            std::vector<std::uint8_t>(g.lattice.size(), 1), std::move(g.edges), std::move(g.masks)};
    const std::uint32_t seed = GENERATE(0u, 1u, 2u);
    CAPTURE(seed);
    std::mt19937 gen{seed};
    Points pts;
    for (std::size_t r = 0; r + 4 < kN; r += 4)
        for (std::size_t c = 0; c + 4 < kN; c += 4) {
            const auto rr = static_cast<double>(r), cc = static_cast<double>(c);
            std::vector<std::pair<double, double>> at;
            for (const double t : {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.25}) at.emplace_back(cc + t, rr + t);
            at.emplace_back(cc + 3.0, rr + 1.0);
            at.emplace_back(cc + 1.0, rr + 3.0);
            for (const auto& [col, row] : at) {
                pts.xy.push_back(world(col, row));
                pts.z.push_back(static_cast<float>(static_cast<double>(gen() % 1537u) - 768.0) / 64.0f);
            }
        }
    const double tol = 1.0;
    const auto out = run(store(pts), s, tol, 1);
    REQUIRE(out.ok());
    REQUIRE(out.inserted > 0);
    double area = 0.0;
    for (const auto& t : out.triangles) {
        const Point2 a = frame(out.vertices[t[0]]), b = frame(out.vertices[t[1]]), c = frame(out.vertices[t[2]]);
        REQUIRE(DefaultKernel::orient2d(a, b, c) == Orientation::CounterClockwise);
        area += cross(a, b, c) / 2.0;
    }
    const double side = static_cast<double>(kN - 1);
    REQUIRE(area == side * side);  // dyadic corners: every term is exact
    REQUIRE(violations(pts, s, out, out.z, tol) == 0);
    delaunay_oracle(out);
    constraints_oracle(out);
    inserted_are_check_points(pts, s, out);
}

TEST_CASE("RP4: output is bit-identical for 1, 2 and 8 threads and two add orders",
          "[refine_points][RP4][determinism]") {
    const Start s = phase1();
    const Points pts = scattered(3);
    const auto ref = run(store(pts), s, 2.0, 1);
    REQUIRE(ref.ok());
    REQUIRE(ref.inserted > 0);
    for (const bool reversed : {false, true})
        for (const unsigned threads : {1u, 2u, 8u}) {
            CAPTURE(reversed, threads);
            const auto out = run(store(pts, reversed), s, 2.0, threads);
            REQUIRE(out.ok());
            REQUIRE(out.vertices == ref.vertices);
            REQUIRE(out.z == ref.z);
            REQUIRE(out.valid == ref.valid);
            REQUIRE(out.triangles == ref.triangles);
            REQUIRE(out.edges == ref.edges);
            REQUIRE(out.masks == ref.masks);
            REQUIRE(out.rounds == ref.rounds);
            REQUIRE(out.inserted == ref.inserted);
            REQUIRE(out.flips == ref.flips);
            REQUIRE(out.max_error == ref.max_error);
            REQUIRE(out.uncovered == ref.uncovered);
            REQUIRE(out.carved == ref.carved);
            REQUIRE(out.coincident == ref.coincident);
            REQUIRE(out.coincident_max_error == ref.coincident_max_error);
        }
}
