// Increment 18 (docs/increments/18-row-span-scan.md, R1 to R1c, "Tests for
// @tester"): the row-span iterator. INVARIANT-CRITICAL, with T2: @reviewer
// runs the mutation round here first, because T1 compares node SETS and an
// argmax comparison caught a dropped vertex exclusion on 8 of 199 313
// triangles in the prototype.
//
// Interface assumed (R1; names the design leaves open are chosen here):
//   namespace terrain::mesh, in <terrain/mesh/row_spans.hpp>:
//     struct RowSpan { std::uint32_t row, c0, c1; std::int8_t flat_edge; };
//     template <std::invocable<RowSpan> F>
//     void for_each_row_span(const std::array<MeshVertex, 3>& v, F&& f);
//     std::int64_t floor_div(std::int64_t n, std::int64_t d);   // CHOSEN
//     std::int64_t ceil_div(std::int64_t n, std::int64_t d);    // CHOSEN
//   The design names floor_div and ceil_div in row_spans.hpp but not their
//   namespace or signature; they are taken as terrain::mesh free functions on
//   int64.
//
// T1: the union of the spans equals scan_oracle::bbox_node_set, the frozen box
// walk (tests/cpp/support/scan_oracle.hpp), as a sorted list of (row, col);
// rows strictly ascending; no empty span; every span inside the node box and
// the grid; flat_edge = k exactly when edge k is horizontal on that row.
//
// dx != dy: the iterator's signature takes lattice vertices only, so it cannot
// see a frame; that independence is structural. The frame is exercised where
// it can matter, in prop_refinement_scan_equivalence (T2), which runs under
// two geometries with dx != dy.

#include <catch2/catch_test_macros.hpp>

#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/row_spans.hpp>

#include "scan_oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <tuple>
#include <vector>

using terrain::mesh::ceil_div;
using terrain::mesh::floor_div;
using terrain::mesh::for_each_row_span;
using terrain::mesh::LatticeVertex;
using terrain::mesh::MeshVertex;
using terrain::mesh::orient_sign;
using terrain::mesh::RowSpan;

namespace {

constexpr std::uint32_t kN = 64;  // the grid is kN x kN nodes
using Tri = std::array<MeshVertex, 3>;

MeshVertex mv(double col, double row) { return MeshVertex{col, row}; }

std::vector<RowSpan> spans_of(const Tri& v) {
    std::vector<RowSpan> out;
    for_each_row_span(v, [&](RowSpan s) { out.push_back(s); });
    return out;
}

std::string show(const Tri& v) {
    std::string s;
    for (const auto& p : v) s += "(" + std::to_string(p.col) + ", " + std::to_string(p.row) + ") ";
    return s;
}

// Checks every T1 property of one counter-clockwise triangle. Returns false
// on the first failure so a random sweep reports one triangle, not thousands.
bool agrees_with_oracle(const Tri& v) {
    INFO("triangle (col, row): " << show(v));
    const auto spans = spans_of(v);
    const auto [rlo, rhi] = std::minmax({v[0].row, v[1].row, v[2].row});
    const auto [clo, chi] = std::minmax({v[0].col, v[1].col, v[2].col});
    std::vector<LatticeVertex> got;
    bool ok = true;
    for (std::size_t i = 0; i < spans.size(); ++i) {
        const RowSpan s = spans[i];
        INFO("span " << i << ": row " << s.row << " [" << s.c0 << ", " << s.c1 << "] flat "
                     << int{s.flat_edge});
        if (i > 0 && !(spans[i - 1].row < s.row)) { FAIL_CHECK("rows not strictly ascending"); ok = false; }
        if (s.c0 > s.c1) { FAIL_CHECK("empty span reported"); ok = false; }
        if (s.row < std::ceil(rlo) || s.row > std::floor(rhi) || s.c0 < std::ceil(clo)
            || s.c1 > std::floor(chi) || s.c1 >= kN || s.row >= kN) {
            FAIL_CHECK("span outside the node box or the grid");
            ok = false;
        }
        int flat = -1;
        for (unsigned k = 0; k < 3; ++k)
            if (v[k].row == v[(k + 1) % 3].row && v[k].row == static_cast<double>(s.row))
                flat = static_cast<int>(k);
        if (int{s.flat_edge} != flat) { FAIL_CHECK("flat_edge is " << int{s.flat_edge} << ", expected " << flat); ok = false; }
        for (std::uint32_t c = s.c0; c <= s.c1 && s.c0 <= s.c1; ++c) got.push_back(LatticeVertex{s.row, c});
    }
    const auto want = scan_oracle::bbox_node_set(v);
    if (got != want) {
        FAIL_CHECK("node set differs from the oracle: " << got.size() << " nodes against "
                                                        << want.size());
        ok = false;
    }
    return ok;
}

// A counter-clockwise triangle from three points, or nothing if collinear.
std::optional<Tri> ccw(MeshVertex a, MeshVertex b, MeshVertex c) {
    const int s = orient_sign(a, b, c);
    if (s > 0) return Tri{a, b, c};
    if (s < 0) return Tri{a, c, b};
    return std::nullopt;
}

// mt19937's raw output is specified by the standard; distributions are not.
struct Rng {
    std::mt19937 gen;
    std::uint32_t node() { return gen() % kN; }
    double unit() { return static_cast<double>(gen() % (1u << 20)) / double(1u << 20); }
    double coord() { return std::min(static_cast<double>(node()) + unit(), double(kN - 1)); }
    bool coin() { return (gen() & 1u) != 0; }
};

// Nudge a node by +-eps in col and/or row, staying inside the grid.
double nudge(Rng& g, std::uint32_t x, double eps) {
    const double d = g.coin() ? eps : -eps;
    const double out = static_cast<double>(x) + d;
    return std::clamp(out, 0.0, double(kN - 1));
}

using Gen = std::optional<Tri> (*)(Rng&);

std::optional<Tri> node_only(Rng& g) {
    return ccw(mv(g.node(), g.node()), mv(g.node(), g.node()), mv(g.node(), g.node()));
}
std::optional<Tri> mixed(Rng& g) {
    auto p = [&] { return g.coin() ? mv(g.node(), g.node()) : mv(g.coord(), g.coord()); };
    return ccw(p(), p(), p());
}
std::optional<Tri> sliver(Rng& g) {
    const auto a = mv(g.node(), g.node()), b = mv(g.node(), g.node());
    const double eps = g.coin() ? 1e-7 : 1e-12;
    const double mc = (a.col + b.col) / 2, mr = (a.row + b.row) / 2;
    const double len = std::hypot(b.col - a.col, b.row - a.row);
    if (len == 0) return std::nullopt;
    const double s = g.coin() ? eps : -eps;  // off the midpoint, perpendicular
    const double c = std::clamp(mc - s * (b.row - a.row) / len, 0.0, double(kN - 1));
    const double r = std::clamp(mr + s * (b.col - a.col) / len, 0.0, double(kN - 1));
    return ccw(a, b, mv(c, r));
}
std::optional<Tri> collinear_but_one(Rng& g) {
    // a and b nodes on a long lattice line; c the lattice node one step off it.
    const std::uint32_t r0 = g.node() / 4, dr = 1 + g.gen() % 3, dc = 1 + g.gen() % 3;
    const std::uint32_t k = std::min((kN - 1 - r0) / dr, (kN - 2) / dc);
    if (k < 2) return std::nullopt;
    const auto a = mv(0, r0), b = mv(double(k * dc), double(r0 + k * dr));
    const std::uint32_t m = k / 2;
    return ccw(a, b, mv(double(m * dc + (g.coin() ? 1 : 0)), double(r0 + m * dr + (g.coin() ? 0 : 1))));
}
std::optional<Tri> near_node(Rng& g) {
    const double eps = g.coin() ? 1e-12 : 1e-9;
    auto p = [&] {
        const auto c = g.node(), r = g.node();
        switch (g.gen() % 3) {
            case 0: return mv(nudge(g, c, eps), r);
            case 1: return mv(c, nudge(g, r, eps));
            default: return mv(nudge(g, c, eps), nudge(g, r, eps));
        }
    };
    return ccw(p(), p(), g.coin() ? mv(g.node(), g.node()) : p());
}
std::optional<Tri> horizontal(Rng& g) {
    // An edge on a row (integer) or between rows (fractional), node or off-node ends.
    const double r = g.coin() ? double(g.node()) : std::min(g.node() + 0.5, kN - 1.0);
    auto x = [&] { return g.coin() ? double(g.node()) : g.coord(); };
    return ccw(mv(x(), r), mv(x(), r), g.coin() ? mv(g.node(), g.node()) : mv(g.coord(), g.coord()));
}
std::optional<Tri> vertical(Rng& g) {
    const double c = g.coin() ? double(g.node()) : g.coord();
    auto y = [&] { return g.coin() ? double(g.node()) : g.coord(); };
    return ccw(mv(c, y()), mv(c, y()), g.coin() ? mv(g.node(), g.node()) : mv(g.coord(), g.coord()));
}
std::optional<Tri> on_border(Rng& g) {
    auto edge = [&] {
        const double e = g.coin() ? 0.0 : double(kN - 1);
        const double t = g.coin() ? double(g.node()) : g.coord();
        return g.coin() ? mv(e, t) : mv(t, e);
    };
    return ccw(edge(), edge(), g.coin() ? mv(g.node(), g.node()) : mv(g.coord(), g.coord()));
}

void sweep(Gen gen, std::uint32_t seed, int count) {
    Rng g{std::mt19937{seed}};
    int tested = 0, failed = 0;
    for (int i = 0; i < count && failed < 3; ++i)
        if (const auto t = gen(g)) {
            ++tested;
            if (!agrees_with_oracle(*t)) ++failed;
        }
    CHECK(failed == 0);
    CHECK(tested > count / 2);  // the generator is not silently all-degenerate
}

// The quarter-circle shape: one off-node centre, an off-node arc, huge slivers.
std::vector<Tri> fan(double cx, double cy, double radius, int n) {
    std::vector<Tri> out;
    const double pi = std::acos(-1.0);
    auto arc = [&](int i) {
        const double t = pi / 2 * i / n;
        return mv(cx - radius * std::cos(t), cy - radius * std::sin(t));
    };
    for (int i = 0; i < n; ++i)
        if (const auto t = ccw(mv(cx, cy), arc(i), arc(i + 1))) out.push_back(*t);
    return out;
}

}  // namespace

// ---------------------------------------------------------------- T1

TEST_CASE("T1 node-only triangles give the oracle's node set", "[mesh][row_spans][T1]") {
    sweep(node_only, 1801, 20000);
}

TEST_CASE("T1 mixed node and off-node triangles give the oracle's node set", "[mesh][row_spans][T1]") {
    sweep(mixed, 1802, 20000);
}

TEST_CASE("T1 slivers a vertex 1e-7 or 1e-12 off an edge midpoint", "[mesh][row_spans][T1]") {
    sweep(sliver, 1803, 20000);
}

TEST_CASE("T1 collinear-but-one lattice slivers", "[mesh][row_spans][T1]") {
    sweep(collinear_but_one, 1804, 20000);
}

TEST_CASE("T1 off-node vertices 1e-12 and 1e-9 from a node on either side", "[mesh][row_spans][T1]") {
    sweep(near_node, 1805, 20000);
}

TEST_CASE("T1 horizontal edges on a row and between rows", "[mesh][row_spans][T1]") {
    sweep(horizontal, 1806, 20000);
}

TEST_CASE("T1 vertical edges", "[mesh][row_spans][T1]") { sweep(vertical, 1807, 20000); }

TEST_CASE("T1 triangles with a vertex on the grid border", "[mesh][row_spans][T1]") {
    sweep(on_border, 1808, 20000);
}

TEST_CASE("T1 huge fans from an off-node centre to an off-node arc", "[mesh][row_spans][T1]") {
    for (const auto& [cx, cy, n] : {std::tuple{62.7, 62.3, 534}, std::tuple{63.0, 63.0, 40},
                                    std::tuple{61.999999999999, 62.000000001, 200}}) {
        const auto tris = fan(cx, cy, 60.4, n);
        REQUIRE(tris.size() == static_cast<std::size_t>(n));
        int failed = 0;
        for (const auto& t : tris)
            if (failed < 3 && !agrees_with_oracle(t)) ++failed;
        CHECK(failed == 0);
    }
}

// ---------------------------------------------------------------- T1-deg

TEST_CASE("T1-deg a row touching only a node vertex gives no span", "[mesh][row_spans][T1-deg]") {
    // Rows 0..4; row 0 holds only the vertex (0, 0); row 4 is edge 1, flat.
    const Tri v{mv(0, 0), mv(0, 4), mv(4, 4)};
    const std::vector<std::array<int, 4>> want{{1, 0, 1, -1}, {2, 0, 2, -1}, {3, 0, 3, -1}, {4, 1, 3, 1}};
    const auto got = spans_of(v);
    REQUIRE(got.size() == want.size());
    for (std::size_t i = 0; i < want.size(); ++i) {
        INFO("span " << i);
        CHECK(int(got[i].row) == want[i][0]);
        CHECK(int(got[i].c0) == want[i][1]);
        CHECK(int(got[i].c1) == want[i][2]);
        CHECK(int{got[i].flat_edge} == want[i][3]);
    }
    CHECK(agrees_with_oracle(v));
}

TEST_CASE("T1-deg a horizontal node edge on a row drops both endpoints and sets flat_edge",
          "[mesh][row_spans][T1-deg]") {
    // Edge 0 from (1, 5) to (6, 5) is flat on row 5; the apex is above it.
    const Tri v{mv(1, 5), mv(6, 5), mv(3, 2)};
    REQUIRE(orient_sign(v[0], v[1], v[2]) > 0);
    const auto got = spans_of(v);
    REQUIRE_FALSE(got.empty());
    const RowSpan last = got.back();
    CHECK(last.row == 5u);
    CHECK(last.c0 == 2u);
    CHECK(last.c1 == 5u);
    CHECK(int{last.flat_edge} == 0);
    CHECK(agrees_with_oracle(v));
}

TEST_CASE("T1-deg a row touching only an off-node vertex gives no span", "[mesh][row_spans][T1-deg]") {
    // Apex (2.5, 0) is on row 0 but not a node; edge 1 is flat on row 3 with
    // off-node ends 0.5 and 4.5, so row 3 holds the nodes between them.
    const Tri v{mv(2.5, 0), mv(0.5, 3), mv(4.5, 3)};
    const std::vector<std::array<int, 4>> want{{1, 2, 3, -1}, {2, 2, 3, -1}, {3, 1, 4, 1}};
    const auto got = spans_of(v);
    REQUIRE(got.size() == want.size());
    for (std::size_t i = 0; i < want.size(); ++i) {
        INFO("span " << i);
        CHECK(int(got[i].row) == want[i][0]);
        CHECK(int(got[i].c0) == want[i][1]);
        CHECK(int(got[i].c1) == want[i][2]);
        CHECK(int{got[i].flat_edge} == want[i][3]);
    }
    CHECK(agrees_with_oracle(v));
}

TEST_CASE("T1-deg a triangle whose only node is one interior node", "[mesh][row_spans][T1-deg]") {
    const Tri v{mv(0.4, 0.4), mv(0.4, 1.7), mv(1.7, 1.0)};
    REQUIRE(orient_sign(v[0], v[1], v[2]) > 0);
    const auto got = spans_of(v);
    REQUIRE(got.size() == 1);
    CHECK(got[0].row == 1u);
    CHECK(got[0].c0 == 1u);
    CHECK(got[0].c1 == 1u);
    CHECK(int{got[0].flat_edge} == -1);
}

TEST_CASE("T1-deg a triangle with an empty node set reports nothing", "[mesh][row_spans][T1-deg]") {
    CHECK(spans_of(Tri{mv(0.2, 0.2), mv(0.2, 0.8), mv(0.8, 0.5)}).empty());
    // A node-only sliver of area 1/2: by Pick, its only lattice points are its vertices.
    const Tri unimodular{mv(0, 0), mv(3, 2), mv(2, 1)};
    REQUIRE(scan_oracle::bbox_node_set(unimodular).empty());
    CHECK(spans_of(unimodular).empty());
}

TEST_CASE("T1-deg floor_div and ceil_div for every sign and remainder", "[mesh][row_spans][T1-deg]") {
    // (n, d, floor, ceil): negative numerators, negative divisors, remainders of each sign, exact.
    const std::vector<std::array<std::int64_t, 4>> cases{
        {7, 2, 3, 4},    {-7, 2, -4, -3}, {7, -2, -4, -3}, {-7, -2, 3, 4},  {6, 3, 2, 2},
        {-6, 3, -2, -2}, {6, -3, -2, -2}, {0, 5, 0, 0},    {0, -5, 0, 0},   {1, 3, 0, 1},
        {-1, 3, -1, 0},  {1, -3, -1, 0},  {-1, -3, 0, 1},
        {(std::int64_t{1} << 48) + 1, 2, std::int64_t{1} << 47, (std::int64_t{1} << 47) + 1},
        {-((std::int64_t{1} << 48) + 1), 2, -(std::int64_t{1} << 47) - 1, -(std::int64_t{1} << 47)},
    };
    for (const auto& [n, d, fl, ce] : cases) {
        INFO(n << " / " << d);
        CHECK(floor_div(n, d) == fl);
        CHECK(ceil_div(n, d) == ce);
    }
}
