// Increment 32 (docs/increments/32-landcover-simplify.md, sections 6 to 9):
// land-cover borders simplified by area-preserving segment collapse, each
// border within a two-sided band of its source, every face's area kept, the
// coverage's topology kept. THE INVARIANT-CRITICAL SUITE of the increment;
// its kill record covers M1 to M6 (section 9).
//
// Interface, as section 6 fixes it:
//
//   #include <terrain/vector_simplify/border_collapse.hpp>
//   namespace terrain::vector_simplify {
//   enum class BorderStatus : std::uint8_t { Ok, InvalidBand, BadRings };
//   struct BorderCounts { junctions, borders, fixed_borders, collinear,
//                         collapses, rejected_crossing, rejected_side; };
//   struct BorderOutcome { std::vector<Point2> points;
//                          std::vector<std::uint64_t> ring_starts;
//                          BorderStatus status; BorderCounts counts; };
//   template <pred::GeometryKernel K>
//   BorderOutcome simplify_borders(std::span<const Point2> points,
//                                  std::span<const std::uint64_t> ring_starts,
//                                  double band);
//   }
//
// THE ORACLES are computed here from the input and the output alone, never
// from the producer's records; geometric decisions use classify<DefaultKernel>
// and orient2d, the exact predicates the producer uses:
//
// - areas: each ring's signed area kept, relative 1e-12 (test 1) or 1e-9
//   unless a test says otherwise (coordinates under 2e3 m: E's rounding is
//   about 1e-13 m there); same sign, so the same orientation;
// - fixed: every input junction (section 6's edge rule, recomputed here) in
//   the same rings, bit for bit; every edge used by one ring or by three or
//   more, present in the same rings in the same direction;
// - a valid coverage: brute force over every pair of distinct output edges,
//   Touching if they share an end point, Disjoint otherwise; every output
//   edge used by at most two rings; the pairs of rings that share an edge
//   unchanged;
// - the same faces: every input vertex still in the output has the same
//   winding number about every ring it is not on, before and after;
// - the band, per border: each ring cut at its junctions into borders, and
//   each border's source and output compared both ways, sampled at every
//   vertex and at every 1/1000 of the polyline's length (a lower estimate of
//   the Hausdorff distance, so it cannot be red on a correct output). Slack
//   1e-9 max(1, band) m unless a test says otherwise.
//
// PINNED HERE where the design leaves it open (listed in the handback):
// - BorderCounts::borders counts every border, the fixed ones (the outline's
//   chains between junctions, a fixed closed ring) included; fixed_borders
//   counts those; junctions counts distinct junction points (test 3, 4).
// - Zero rings (`ring_starts == {0}`, no points) is Ok with zero rings out;
//   an empty `ring_starts` is BadRings.
// - collinear: at least the number of collinear vertices in the border,
//   whether a shared vertex counts once or per ring.
// - A refused collapse is counted in rejected_side when test (ii) refuses it
//   (test 5) and in rejected_crossing when test (i) does (test 6).
// - A collinear vertex on the outline is kept: its edges are fixed (test 3).
//
// RED at the commit that adds this file: the header does not exist, so the
// target test_border_collapse does not compile ("'terrain/vector_simplify/
// border_collapse.hpp' file not found"); every other target builds.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/vector_simplify/border_collapse.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <span>
#include <thread>
#include <utility>
#include <vector>

using terrain::Point2;
using terrain::Segment2;
using terrain::noding::classify;
using terrain::noding::SegmentRelation;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;
using terrain::vector_simplify::BorderOutcome;
using terrain::vector_simplify::BorderStatus;
using terrain::vector_simplify::simplify_borders;

namespace {

using Ring = std::vector<Point2>;
using Rings = std::vector<Ring>;
constexpr double nan_v = std::numeric_limits<double>::quiet_NaN();
constexpr double inf_v = std::numeric_limits<double>::infinity();

struct Flat {
    std::vector<Point2> points;
    std::vector<std::uint64_t> starts;
};

Flat flatten(const Rings& rings) {
    Flat f{{}, {0}};
    for (const Ring& r : rings) {
        f.points.insert(f.points.end(), r.begin(), r.end());
        f.starts.push_back(f.points.size());
    }
    return f;
}

BorderOutcome run_flat(const std::vector<Point2>& points, const std::vector<std::uint64_t>& starts,
                       double band) {
    return simplify_borders<DefaultKernel>(std::span<const Point2>{points},
                                           std::span<const std::uint64_t>{starts}, band);
}

BorderOutcome run(const Rings& rings, double band) {
    const Flat f = flatten(rings);
    return run_flat(f.points, f.starts, band);
}

Rings rings_of(const BorderOutcome& out) {
    Rings rings;
    REQUIRE(!out.ring_starts.empty());
    REQUIRE(out.ring_starts.back() == out.points.size());
    for (std::size_t k = 0; k + 1 < out.ring_starts.size(); ++k) {
        REQUIRE(out.ring_starts[k] < out.ring_starts[k + 1]);
        rings.emplace_back(out.points.begin() + static_cast<std::ptrdiff_t>(out.ring_starts[k]),
                           out.points.begin() + static_cast<std::ptrdiff_t>(out.ring_starts[k + 1]));
    }
    return rings;
}

// ------------------------------------------------------------- geometry

// Signed, relative to the first vertex: at UTM-sized coordinates the plain
// shoelace's products (3e12 m2) round by 5e-4 m2, more than a small ring's
// area; the differences here are exact (Sterbenz) for nearby vertices.
double area(const Ring& r) {
    double twice = 0.0;
    for (std::size_t i = 1; i + 1 < r.size(); ++i)
        twice += cross(r[i] - r[0], r[i + 1] - r[0]);
    return 0.5 * twice;
}

double distance(const Point2& p, const Point2& a, const Point2& b) {
    const double dx = b.x - a.x, dy = b.y - a.y;
    const double len2 = dx * dx + dy * dy;
    const double t = len2 > 0.0 ? std::clamp(((p.x - a.x) * dx + (p.y - a.y) * dy) / len2, 0.0, 1.0) : 0.0;
    return std::hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
}

double distance(const Point2& p, const Ring& line) {
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i + 1 < line.size(); ++i)
        best = std::min(best, distance(p, line[i], line[i + 1]));
    return best;
}

// Every vertex and every 1/1000 of the polyline's length.
std::vector<Point2> samples(const Ring& line) {
    std::vector<Point2> out(line.begin(), line.end());
    double total = 0.0;
    for (std::size_t i = 0; i + 1 < line.size(); ++i)
        total += std::hypot(line[i + 1].x - line[i].x, line[i + 1].y - line[i].y);
    std::size_t seg = 0;
    double before = 0.0;
    for (int k = 1; k < 1000 && line.size() > 1; ++k) {
        const double s = total * k / 1000.0;
        while (seg + 2 < line.size()
               && before + std::hypot(line[seg + 1].x - line[seg].x, line[seg + 1].y - line[seg].y) < s) {
            before += std::hypot(line[seg + 1].x - line[seg].x, line[seg + 1].y - line[seg].y);
            ++seg;
        }
        const double len = std::hypot(line[seg + 1].x - line[seg].x, line[seg + 1].y - line[seg].y);
        const double t = len > 0.0 ? std::clamp((s - before) / len, 0.0, 1.0) : 0.0;
        out.push_back(line[seg] + t * (line[seg + 1] - line[seg]));
    }
    return out;
}

// The two directed sampled distances, source to output and output to source.
std::pair<double, double> both_ways(const Ring& source, const Ring& output) {
    double fwd = 0.0, back = 0.0;
    for (const Point2& p : samples(source))
        fwd = std::max(fwd, distance(p, output));
    for (const Point2& p : samples(output))
        back = std::max(back, distance(p, source));
    return {fwd, back};
}

bool on_ring(const Point2& p, const Ring& r) {
    for (std::size_t i = 0; i < r.size(); ++i)
        if (classify<DefaultKernel>(Segment2{p, p}, Segment2{r[i], r[(i + 1) % r.size()]})
            != SegmentRelation::Disjoint)
            return true;
    return false;
}

int winding(const Point2& p, const Ring& r) {
    int w = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
        const Point2& a = r[i];
        const Point2& b = r[(i + 1) % r.size()];
        const Orientation o = DefaultKernel::orient2d(a, b, p);
        if (a.y <= p.y && b.y > p.y && o == Orientation::CounterClockwise)
            ++w;
        if (a.y > p.y && b.y <= p.y && o == Orientation::Clockwise)
            --w;
    }
    return w;
}

bool same_bits(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }
bool same_bits(const Point2& a, const Point2& b) { return same_bits(a.x, b.x) && same_bits(a.y, b.y); }
bool same_bits(const std::vector<Point2>& a, const std::vector<Point2>& b) {
    return a.size() == b.size()
           && std::equal(a.begin(), a.end(), b.begin(), [](const Point2& p, const Point2& q) { return same_bits(p, q); });
}

// The same cycle, bit for bit, up to where it starts.
bool same_cycle(const Ring& a, const Ring& b) {
    if (a.size() != b.size())
        return false;
    for (std::size_t s = 0; s < a.size(); ++s) {
        bool all = true;
        for (std::size_t i = 0; i < a.size() && all; ++i)
            all = same_bits(a[(s + i) % a.size()], b[i]);
        if (all)
            return true;
    }
    return a.empty();
}

Ring reversed(Ring r) {
    std::reverse(r.begin(), r.end());
    return r;
}

// ------------------------------------------------- the design's edge rule

// A vertex by its exact coordinates, -0.0 and 0.0 the same (section 8).
using Key = std::pair<double, double>;
Key key(const Point2& p) { return {p.x + 0.0, p.y + 0.0}; }
using EdgeKey = std::pair<Key, Key>;
EdgeKey edge_key(const Point2& a, const Point2& b) {
    const Key u = key(a), v = key(b);
    return u < v ? EdgeKey{u, v} : EdgeKey{v, u};
}

struct Topology {
    std::map<EdgeKey, std::multiset<std::size_t>> rings_of_edge;
    std::set<Key> junctions;
};

// Section 6: an edge used by exactly two rings is shared; by one or by three
// or more, fixed. A junction has other than two distinct incident edges, or
// two whose rings differ.
Topology topology(const Rings& rings) {
    Topology t;
    for (std::size_t k = 0; k < rings.size(); ++k)
        for (std::size_t i = 0; i < rings[k].size(); ++i)
            t.rings_of_edge[edge_key(rings[k][i], rings[k][(i + 1) % rings[k].size()])].insert(k);
    std::map<Key, std::vector<EdgeKey>> incident;
    for (const auto& [e, _] : t.rings_of_edge) {
        incident[e.first].push_back(e);
        incident[e.second].push_back(e);
    }
    for (const auto& [v, edges] : incident)
        if (edges.size() != 2 || t.rings_of_edge.at(edges[0]) != t.rings_of_edge.at(edges[1]))
            t.junctions.insert(v);
    return t;
}

// Ring `r` cut at its junctions into borders, each from one junction to the
// next, both ends included; a ring with no junction is one closed border
// (first vertex repeated at the end). The cut starts at junction `first` if
// given, else at the ring's first junction.
std::vector<Ring> borders_of(const Ring& r, const std::set<Key>& junctions, const Key* first = nullptr) {
    std::vector<std::size_t> at;
    for (std::size_t i = 0; i < r.size(); ++i)
        if (junctions.contains(key(r[i])))
            at.push_back(i);
    if (at.empty()) {
        Ring closed = r;
        closed.push_back(r.front());
        return {closed};
    }
    std::size_t s0 = 0;
    if (first != nullptr)
        while (s0 < at.size() && key(r[at[s0]]) != *first)
            ++s0;
    REQUIRE(s0 < at.size());
    std::vector<Ring> out;
    for (std::size_t j = 0; j < at.size(); ++j) {
        const std::size_t from = at[(s0 + j) % at.size()], to = at[(s0 + j + 1) % at.size()];
        Ring piece{r[from]};
        for (std::size_t i = (from + 1) % r.size();; i = (i + 1) % r.size()) {
            piece.push_back(r[i]);
            if (i == to)
                break;
        }
        out.push_back(piece);
    }
    return out;
}

// Every guarantee of section 7 that the output shows, computed from input and
// output alone. `slack` below 0 means 1e-9 max(1, band).
Rings check(const Rings& in, const BorderOutcome& out, double band, double area_rel = 1e-9,
            double slack = -1.0, bool valid_coverage = true) {
    REQUIRE(out.status == BorderStatus::Ok);
    const Rings got = rings_of(out);
    REQUIRE(got.size() == in.size());
    slack = slack >= 0.0 ? slack : 1e-9 * std::max(1.0, band);
    const Topology before = topology(in);
    const Topology after = topology(got);

    for (std::size_t k = 0; k < in.size(); ++k) {
        INFO("ring " << k);
        REQUIRE(got[k].size() >= 3);
        for (const Point2& p : got[k])
            REQUIRE((std::isfinite(p.x) && std::isfinite(p.y)));
        // 1. Area and orientation.
        const double a0 = area(in[k]), a1 = area(got[k]);
        CHECK(std::abs(a1 - a0) <= area_rel * std::abs(a0));
        CHECK((a0 > 0.0) == (a1 > 0.0));
        // 4. Junctions, bit for bit up to the sign of zero, in the same rings
        // (section 8: -0.0 and 0.0 are one vertex; which sign a vertex shared
        // by two rings comes out with is not pinned).
        for (const Point2& p : in[k])
            if (before.junctions.contains(key(p)))
                CHECK(std::any_of(got[k].begin(), got[k].end(), [&](const Point2& q) { return key(p) == key(q); }));
        // 4. Fixed edges: same rings, same direction, as the junctions.
        for (std::size_t i = 0; i < in[k].size(); ++i) {
            const Point2& a = in[k][i];
            const Point2& b = in[k][(i + 1) % in[k].size()];
            if (before.rings_of_edge.at(edge_key(a, b)).size() == 2)
                continue;
            bool found = false;
            for (std::size_t j = 0; j < got[k].size() && !found; ++j)
                found = key(got[k][j]) == key(a) && key(got[k][(j + 1) % got[k].size()]) == key(b);
            CHECK(found);
        }
    }

    // 3. A valid coverage.
    if (valid_coverage) {
        std::vector<std::pair<Point2, Point2>> edges;
        for (const auto& [e, users] : after.rings_of_edge) {
            CHECK(users.size() <= 2);
            edges.push_back({Point2{e.first.first, e.first.second}, Point2{e.second.first, e.second.second}});
        }
        for (std::size_t i = 0; i < edges.size(); ++i)
            for (std::size_t j = i + 1; j < edges.size(); ++j) {
                const auto& [a, b] = edges[i];
                const auto& [c, d] = edges[j];
                const bool shared = a == c || a == d || b == c || b == d;
                const SegmentRelation rel = classify<DefaultKernel>(Segment2{a, b}, Segment2{c, d});
                INFO("edges " << i << " and " << j);
                CHECK(rel == (shared ? SegmentRelation::Touching : SegmentRelation::Disjoint));
            }
        const auto pairs = [](const Topology& t) {
            std::set<std::pair<std::size_t, std::size_t>> s;
            for (const auto& [_, users] : t.rings_of_edge)
                if (users.size() == 2)
                    s.insert({*users.begin(), *users.rbegin()});
            return s;
        };
        CHECK(pairs(before) == pairs(after));
        // The same faces: an input vertex still present keeps its winding
        // number about every ring it is on neither before nor after.
        for (std::size_t r = 0; r < in.size(); ++r)
            for (const Point2& p : in[r]) {
                if (!std::any_of(got[r].begin(), got[r].end(), [&](const Point2& q) { return same_bits(p, q); }))
                    continue;
                for (std::size_t q = 0; q < in.size(); ++q) {
                    if (q == r || on_ring(p, in[q]) || on_ring(p, got[q]))
                        continue;
                    INFO("vertex of ring " << r << " about ring " << q);
                    CHECK(winding(p, in[q]) == winding(p, got[q]));
                }
            }
    }

    // 2. The band, per border, both ways.
    for (std::size_t k = 0; k < in.size(); ++k) {
        const std::vector<Ring> src = borders_of(in[k], before.junctions);
        const Key* first = nullptr;
        Key start{};
        if (std::any_of(in[k].begin(), in[k].end(), [&](const Point2& p) { return before.junctions.contains(key(p)); })) {
            start = key(src.front().front());
            first = &start;
        }
        const std::vector<Ring> dst = borders_of(got[k], before.junctions, first);
        REQUIRE(src.size() == dst.size());
        for (std::size_t b = 0; b < src.size(); ++b) {
            if (first != nullptr) {
                REQUIRE(key(src[b].front()) == key(dst[b].front()));
                REQUIRE(key(src[b].back()) == key(dst[b].back()));
            }
            const auto [fwd, back] = both_ways(src[b], dst[b]);
            INFO("ring " << k << ", border " << b << ": source to output " << fwd << ", output to source " << back);
            CHECK(fwd <= band + slack);
            CHECK(back <= band + slack);
        }
    }
    return got;
}

std::size_t vertices(const Rings& rings) {
    std::size_t n = 0;
    for (const Ring& r : rings)
        n += r.size();
    return n;
}

bool has(const Ring& r, const Point2& p) {
    return std::any_of(r.begin(), r.end(), [&](const Point2& q) { return same_bits(p, q); });
}

// ----------------------------------------------------------- fixtures

// A zig-zag from `from` to `to` (exclusive of both), `n` steps, offsets of
// +-`amp` across the line, alternating, starting with +.
Ring zigzag(Point2 from, Point2 to, int n, double amp) {
    const Point2 d = to - from;
    const double len = std::hypot(d.x, d.y);
    const Point2 across{-d.y / len, d.x / len};
    Ring out;
    for (int k = 1; k < n; ++k)
        out.push_back(from + (static_cast<double>(k) / n) * d + ((k % 2 == 1) ? amp : -amp) * across);
    return out;
}

Ring cat(std::initializer_list<Ring> parts) {
    Ring out;
    for (const Ring& p : parts)
        out.insert(out.end(), p.begin(), p.end());
    return out;
}

// Test 1: two squares 100 m on a side sharing a zig-zag border at x = 100,
// 49 inner vertices, 3 m either side; junctions (100, 0) and (100, 100).
// `o` shifts everything; `s` scales.
Rings two_squares(Point2 o = {0.0, 0.0}, double s = 1.0, double amp = 3.0, int n = 50) {
    const auto at = [&](double x, double y) { return Point2{o.x + s * x, o.y + s * y}; };
    const Ring up = zigzag(at(100, 0), at(100, 100), n, s * amp);
    const Ring left = cat({{at(0, 0), at(100, 0)}, up, {at(100, 100), at(0, 100)}});
    const Ring right = cat({{at(100, 0), at(200, 0), at(200, 100), at(100, 100)}, reversed(up)});
    return {left, right};
}

// Test 2: the CORINE border of the German case saved as
// docs/benchmarks/2026-10-08/clc-simplify/design-probe/overshoot_fixture.json
// (177 vertices, local origin, rounded to 1 cm), copied here verbatim so the
// suite reads no file.
const Ring overshoot_border{
    {0.0, 0.0}, {31.13, -91.87}, {10.24, -192.43}, {-39.3, -254.85}, {-83.33, -297.07},
    {-91.75, -383.69}, {-180.58, -335.63}, {-176.53, -408.3}, {-146.05, -479.5}, {-147.12, -578.97},
    {-194.81, -674.41}, {-121.4, -683.56}, {-71.13, -634.37}, {-74.82, -568.3}, {-53.19, -480.95},
    {-6.57, -484.97}, {2.24, -524.24}, {17.3, -556.54}, {58.05, -574.15}, {54.73, -514.68},
    {70.5, -440.91}, {114.53, -398.69}, {147.93, -403.45}, {113.09, -491.55}, {173.65, -508.05},
    {166.34, -614.49}, {127.8, -636.53}, {80.82, -625.9}, {98.08, -697.83}, {117.93, -816.01},
    {104.01, -922.82}, {26.18, -834.39}, {26.89, -728.31}, {-1.38, -696.75}, {-59.36, -726.5},
    {-35.12, -804.67}, {22.53, -887.61}, {53.38, -965.41}, {70.65, -1037.34}, {25.87, -1066.35},
    {-10.45, -1128.02}, {-66.22, -1197.4}, {-151.0, -1222.02}, {-215.22, -1258.74},
    {-268.81, -1248.48}, {-248.65, -1134.69}, {-177.46, -1104.2}, {-184.09, -985.29},
    {-239.89, -935.39}, {-305.96, -939.07}, {-325.01, -1072.68}, {-341.49, -1252.53},
    {-409.77, -1216.58}, {-458.96, -1166.31}, {-402.08, -1116.74}, {-374.95, -1009.2},
    {-334.97, -894.3}, {-187.81, -799.93}, {-219.4, -708.92}, {-350.76, -848.79},
    {-416.45, -859.09}, {-341.6, -775.38}, {-375.41, -644.73}, {-389.73, -625.64},
    {-451.73, -702.01}, {-509.72, -731.75}, {-568.81, -741.67}, {-525.85, -798.92},
    {-597.05, -829.4}, {-643.3, -831.98}, {-641.85, -739.13}, {-702.79, -716.01},
    {-752.31, -778.42}, {-821.67, -841.93}, {-850.65, -916.46}, {-743.47, -936.98},
    {-773.56, -991.68}, {-845.11, -1015.55}, {-924.03, -1026.58}, {-990.09, -1030.27},
    {-1101.3, -1056.35}, {-1099.82, -1082.79}, {-1157.43, -1119.13}, {-1251.4, -1097.87},
    {-1263.15, -1125.03}, {-1295.07, -1146.7}, {-1272.67, -1191.84}, {-1262.37, -1257.53},
    {-1156.67, -1251.63}, {-1194.47, -1286.88}, {-1246.21, -1309.65}, {-1270.79, -1344.15},
    {-1209.12, -1380.47}, {-1111.5, -1348.52}, {-1071.87, -1346.3}, {-984.14, -1374.54},
    {-912.2, -1357.27}, {-862.63, -1414.16}, {-804.65, -1384.41}, {-775.54, -1484.69},
    {-703.07, -1490.96}, {-659.03, -1448.74}, {-614.25, -1419.73}, {-551.85, -1469.26},
    {-500.84, -1433.28}, {-465.99, -1345.19}, {-438.09, -1370.13}, {-394.4, -1440.6},
    {-354.02, -1451.6}, {-342.65, -1417.83}, {-276.59, -1414.14}, {-300.8, -1455.26},
    {-221.15, -1457.44}, {-170.14, -1421.46}, {-130.5, -1419.24}, {-77.65, -1416.29},
    {-8.3, -1352.78}, {8.24, -1411.5}, {4.95, -1471.32}, {19.27, -1490.41}, {65.88, -1494.44},
    {117.62, -1471.67}, {160.19, -1403.02}, {204.59, -1367.41}, {250.46, -1358.22},
    {240.91, -1305.73}, {262.54, -1218.38}, {299.97, -1176.52}, {341.82, -1213.95},
    {344.21, -1256.81}, {372.96, -1260.52}, {416.25, -1294.12}, {467.1, -1253.8},
    {511.87, -1224.79}, {561.77, -1168.99}, {592.22, -1120.9}, {623.04, -1079.42},
    {606.52, -1020.69}, {544.84, -984.37}, {478.77, -988.06}, {493.47, -1013.75},
    {543.39, -1077.23}, {450.54, -1075.79}, {391.08, -1079.1}, {330.14, -1056.0}, {274.34, -1006.1},
    {218.2, -1068.87}, {228.49, -1134.57}, {180.77, -1110.72}, {163.87, -1045.4}, {167.16, -985.56},
    {182.19, -898.58}, {206.4, -857.46}, {207.11, -751.39}, {210.02, -684.96}, {259.56, -622.54},
    {328.91, -559.03}, {259.16, -496.65}, {227.94, -412.24}, {228.28, -299.55}, {180.9, -163.03},
    {162.89, -77.88}, {187.48, -43.37}, {217.59, -107.96}, {242.94, -205.96}, {266.08, -264.32},
    {278.95, -376.26}, {301.35, -421.4}, {329.6, -333.67}, {278.89, -137.68}, {338.69, -21.68},
    {357.77, -7.36}, {364.04, -119.67}, {383.86, -118.57}, {377.65, -244.84}, {465.01, -266.46},
    {466.1, -274.26},
};

// The border closed into two polygons by a frame at least 93 m from it away
// from its ends (604 m beyond the two notches): the frame's top at y = 600
// meets the border's start S = (0, 0) through a notch to (+-300, 600), and
// its right side meets the end T through a notch to (1300, T.y) and
// (1300, T.y - 300), so S and T are junctions (three distinct edges) and the
// border is the whole fixture. Both rings counter-clockwise.
Rings overshoot_coverage() {
    const Ring& b = overshoot_border;
    const Point2 T = b.back();
    const Ring a = cat({b, {{1300.0, T.y}, {1300.0, 600.0}, {300.0, 600.0}}});
    const Ring c = cat({reversed(b),
                        {{-300.0, 600.0}, {-1900.0, 600.0}, {-1900.0, -2100.0}, {1300.0, -2100.0}, {1300.0, T.y - 300.0}}});
    return {area(a) > 0.0 ? a : reversed(a), area(c) > 0.0 ? c : reversed(c)};
}

// Test 2b: an out-and-back border, a tongue of Q reaching into P from the
// outline at x = 0, its two arms 10.69 m apart at their closest against a
// band of 20 m. (Drawn from the M6 probe's generator, seed 130, rounded to
// 1 cm; see the handback.)
const Ring tongue_border{
    {0.0, 0.0}, {11.12, 2.87}, {24.26, -3.79}, {44.33, 2.77}, {58.45, -2.3}, {63.64, 4.28},
    {78.78, -3.59}, {98.08, 0.07}, {103.09, -1.12}, {116.41, -1.51}, {136.62, -3.64}, {155.59, 6.52},
    {164.08, -6.53}, {188.21, -6.33}, {202.81, -1.32}, {217.39, 9.91}, {229.2, 0.14}, {244.2, 10.31},
    {226.35, 16.37}, {219.82, 24.19}, {201.28, 20.48}, {188.95, 21.45}, {162.17, 18.0}, {153.33, 25.32},
    {135.46, 25.12}, {119.81, 17.21}, {103.57, 23.87}, {101.31, 17.05}, {80.72, 20.59}, {62.9, 18.48},
    {58.86, 19.13}, {45.69, 22.06}, {25.17, 22.05}, {9.03, 23.06}, {0.0, 20.62},
};

Rings tongue_coverage() {
    const Ring q = tongue_border; // counter-clockwise: out along y ~ 0, back along y ~ 20
    const Ring p = cat({{{0.0, -200.0}, {600.0, -200.0}, {600.0, 200.0}, {0.0, 200.0}}, reversed(tongue_border)});
    return {p, q};
}

// Tests 5 to 7: a border J1-B-C-J2 with a trapezoid bulge into Q, J1 = (0, 0)
// and J2 = (100, 0) on the outline. Its one collapse puts E at (32, 16) (on
// line J1-B) or (68, 16) (on line C-J2): the first sweeps the ground near
// (75, 8) out of P, the second the ground near (25, 8).
const Point2 J1{0.0, 0.0}, B{20.0, 10.0}, C{80.0, 10.0}, J2{100.0, 0.0};

Ring rect(double x0, double y0, double x1, double y1) { return {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}}; }

} // namespace

// ---------------------------------------------------------------- tests

TEST_CASE("1. two squares sharing a zig-zag: areas kept, the band both ways, fewer vertices",
          "[vector_simplify][border_collapse][area][band]") {
    const Rings in = two_squares();
    const BorderOutcome out = run(in, 10.0);
    const Rings got = check(in, out, 10.0, 1e-12);
    CHECK(vertices(got) < vertices(in));
    CHECK(out.counts.collapses > 0);
}

TEST_CASE("1. extreme scale: a sub-millimetre zig-zag at UTM-sized coordinates",
          "[vector_simplify][border_collapse][scale]") {
    // Scale: coordinates near 5e5 and 6.6e6 m, where one ulp is 9.3e-10 m;
    // two 0.1 m squares, the border between them 0.1 m long with a zig-zag
    // 0.5 mm either side in 1 mm steps (99 inner vertices), the band 2 mm.
    // Area bound 1e-5 relative (1e-7 m2 on 0.01 m2): E's rounding (one ulp)
    // times the border's length times the 99 collapses at most is 9.2e-9
    // m2, with a margin of 10. Band slack 1e-8 m, ten ulps. Checked at this
    // fixture only.
    const Rings in = two_squares({500'000.3, 6'600'000.7}, 0.001, 0.5, 100);
    const BorderOutcome out = run(in, 0.002);
    const Rings got = check(in, out, 0.002, 1e-5, 1e-8);
    CHECK(vertices(got) < vertices(in));
}

TEST_CASE("2. the overshoot border: within 50 m both ways under the anchored check (M4)",
          "[vector_simplify][border_collapse][band][overshoot]") {
    // On this border alone today's check (deviation_of) leaves the source
    // 58.7 m from the output, the anchored one 48.3 m (iso.py beside the
    // fixture); with this frame the M4 probe gave 58.3 m by this oracle.
    const Rings in = overshoot_coverage();
    REQUIRE(in[0].size() == overshoot_border.size() + 3);
    const BorderOutcome out = run(in, 50.0);
    const Rings got = check(in, out, 50.0);
    CHECK(vertices(got) < vertices(in));
}

TEST_CASE("2b. an out-and-back border whose arms are closer than the band: the band both ways",
          "[vector_simplify][border_collapse][band][anchor_order]") {
    // Section 9 asks this chain to kill M6 (the anchor-order condition
    // dropped). The probe in the handback found M6 changes the output here
    // (and on 44 of 1 021 such chains) but never past the band: the order
    // can break only on F_A's or F_D's own segment, where the reversed piece
    // holds no source vertex. Kept as a band test of the out-and-back shape.
    const Rings in = tongue_coverage();
    const BorderOutcome out = run(in, 20.0);
    const Rings got = check(in, out, 20.0);
    CHECK(vertices(got) < vertices(in));
}

TEST_CASE("3. junctions and the outline stay bit for bit (M3); the counts of borders",
          "[vector_simplify][border_collapse][junction]") {
    // Three polygons in a 300 m square meet at J = (150, 150); their borders,
    // zig-zags 2 m either side, run from J to (150, 0), (150, 300) and
    // (300, 150) on the outline. The two vertical ones are in line through J,
    // so J is the best collapse there is if it may move. The outline carries
    // collinear vertices at (75, 0), (225, 0), (300, 75), (300, 225),
    // (0, 200) and (0, 100): its edges are fixed, so they stay.
    const Point2 J{150.0, 150.0}, S{150.0, 0.0}, N{150.0, 300.0}, E{300.0, 150.0};
    const Ring down = zigzag(J, S, 15, 2.0), up = zigzag(J, N, 15, 2.0), right = zigzag(J, E, 15, 2.0);
    const Ring west = cat({{{0.0, 0.0}, {75.0, 0.0}, S}, reversed(down), {J}, up, {N, {0.0, 300.0}, {0.0, 200.0}, {0.0, 100.0}}});
    const Ring south_east = cat({{S, {225.0, 0.0}, {300.0, 0.0}, {300.0, 75.0}, E}, reversed(right), {J}, down});
    const Ring north_east = cat({{J}, right, {E, {300.0, 225.0}, {300.0, 300.0}, N}, reversed(up)});
    const Rings in{west, south_east, north_east};
    for (const Ring& r : in)
        REQUIRE(area(r) > 0.0);

    const BorderOutcome out = run(in, 20.0);
    const Rings got = check(in, out, 20.0);
    CHECK(out.counts.collapses > 0);
    for (const Ring& r : got)
        CHECK(has(r, J));
    CHECK(has(got[0], S));
    CHECK(has(got[1], S));
    CHECK(has(got[1], E));
    CHECK(has(got[2], E));
    CHECK(has(got[0], N));
    CHECK(has(got[2], N));
    for (const Point2 p : {Point2{75.0, 0.0}, Point2{0.0, 200.0}, Point2{0.0, 100.0}})
        CHECK(has(got[0], p));
    for (const Point2 p : {Point2{225.0, 0.0}, Point2{300.0, 75.0}})
        CHECK(has(got[1], p));
    CHECK(has(got[2], Point2{300.0, 225.0}));
    // PINNED: J, S, N, E; three inner borders and three outline chains.
    CHECK(out.counts.junctions == 4);
    CHECK(out.counts.borders == 6);
    CHECK(out.counts.fixed_borders == 3);
}

TEST_CASE("3. a loop through one junction keeps it: an island touching its hole's outline at J (M3)",
          "[vector_simplify][border_collapse][junction][loop]") {
    // A 100 m square whose hole is an island Q touching the outline at one
    // point, J = (50, 0), a vertex of both. J has four edges, so it is a
    // junction, and Q's ring (and the hole's) is one border from J round to
    // J (section 8: "a loop's junction never moves"). J is the tip of a
    // 10 m spike: the collapse that removes it puts E near (46.7, 3.3),
    // within the band of 10 m and crossing nothing, so only the rule that
    // a junction is never B or C keeps it. The mutation round's M3 (the loop
    // made a closed border, so J is an inner node) moves it.
    const Point2 J{50.0, 0.0};
    const Ring Q{J, {60.0, 10.0}, {70.0, 10.0}, {70.0, 30.0}, {30.0, 30.0}, {30.0, 10.0}, {40.0, 10.0}};
    const Ring shell{{0.0, 0.0}, J, {100.0, 0.0}, {100.0, 100.0}, {0.0, 100.0}};
    const Rings in{shell, reversed(Q), Q};
    REQUIRE(area(Q) > 0.0);
    const BorderOutcome out = run(in, 10.0);
    const Rings got = check(in, out, 10.0);
    for (const Ring& r : got)
        CHECK(has(r, J));
    CHECK(same_cycle(got[0], shell));
    CHECK(same_cycle(got[1], reversed(got[2])));
}

TEST_CASE("4. an island: its ring and the hole's are the same points in opposite order",
          "[vector_simplify][border_collapse][island]") {
    // A 64-gon on a circle of radius 100 at (500, 500) (cocircular up to
    // rounding), a hole in a 1000 m square.
    Ring island;
    for (int k = 0; k < 64; ++k) {
        const double t = 2.0 * std::numbers::pi * k / 64.0;
        island.push_back({500.0 + 100.0 * std::cos(t), 500.0 + 100.0 * std::sin(t)});
    }
    const Rings in{rect(0.0, 0.0, 1000.0, 1000.0), reversed(island), island};
    const BorderOutcome out = run(in, 20.0);
    const Rings got = check(in, out, 20.0);
    CHECK(out.counts.collapses > 0);
    CHECK(got[2].size() >= 4);
    CHECK(same_cycle(got[1], reversed(got[2])));
    // PINNED: no junction; the outline's closed ring and the island's border.
    CHECK(out.counts.junctions == 0);
    CHECK(out.counts.borders == 2);
    CHECK(out.counts.fixed_borders == 1);
}

TEST_CASE("5. a bulge whose collapse would sweep a whole island leaves it in its face (M1)",
          "[vector_simplify][border_collapse][side]") {
    // Two 2 x 1 m islands in P, at (75, 8) and (25, 8), 1.5 m below the bulge's
    // top: either E sweeps one of them out of P without crossing it.
    const Ring P{{0.0, -50.0}, {100.0, -50.0}, J2, C, B, J1};
    const Ring Q{J1, B, C, J2, {100.0, 60.0}, {0.0, 60.0}};
    const Ring i1 = rect(74.0, 7.5, 76.0, 8.5), i2 = rect(24.0, 7.5, 26.0, 8.5);
    const Rings in{P, reversed(i1), reversed(i2), Q, i1, i2};
    const BorderOutcome out = run(in, 10.0);
    const Rings got = check(in, out, 10.0);
    for (const Ring& island : {got[4], got[5]})
        for (const Point2& p : island)
            CHECK(winding(p, got[0]) == 1);
    CHECK(out.counts.rejected_side >= 1);
}

TEST_CASE("6. borders of two polygon pairs side by side closer than the band do not cross (M2)",
          "[vector_simplify][border_collapse][crossing]") {
    // Y, the P2|P3 border, is straight at y = 13, 3 m above the bulge's top;
    // either E (y = 16) lies beyond it, and Y has no vertex near.
    const Ring P1{{0.0, -50.0}, {100.0, -50.0}, J2, C, B, J1};
    const Ring P2{J1, B, C, J2, {100.0, 13.0}, {0.0, 13.0}};
    const Ring P3{{0.0, 13.0}, {100.0, 13.0}, {100.0, 60.0}, {0.0, 60.0}};
    const Rings in{P1, P2, P3};
    const BorderOutcome out = run(in, 10.0);
    check(in, out, 10.0);
    CHECK(out.counts.rejected_crossing >= 1);
}

TEST_CASE("7. a junction with another border leaving into the region a collapse would sweep: refused",
          "[vector_simplify][border_collapse][junction][side]") {
    // Borders leave J1 to (15, 5.5) and J2 to (85, 5.5), each into the wedge
    // between the old and the new edge at its end, then run down to the
    // outline. Every collapse there is refused; nothing else can collapse.
    const Point2 z{15.0, 5.5}, zp{85.0, 5.5};
    const Ring left{{0.0, -50.0}, {15.0, -50.0}, z, J1};
    const Ring mid{{15.0, -50.0}, {85.0, -50.0}, zp, J2, C, B, J1, z};
    const Ring right{{85.0, -50.0}, {100.0, -50.0}, J2, zp};
    const Ring Q{J1, B, C, J2, {100.0, 60.0}, {0.0, 60.0}};
    const Rings in{left, mid, right, Q};
    const BorderOutcome out = run(in, 10.0);
    const Rings got = check(in, out, 10.0);
    CHECK(out.counts.collapses == 0);
    for (std::size_t k = 0; k < in.size(); ++k)
        CHECK(same_cycle(got[k], in[k]));
}

TEST_CASE("8. area sign: an asymmetric chain whose E on the wrong side would change the area (M5)",
          "[vector_simplify][border_collapse][area]") {
    // A 60 m border along y = 0 with offsets 1.0, -0.3, 0.8, -0.5 m: each
    // chain's area is small and one-sided, so E on the wrong side of the
    // equal-area line is within the band of 2 m and the area changes by
    // twice the chain's (about 1e-3 relative in the probe).
    const Ring inner{{10.0, 1.0}, {20.0, -0.3}, {30.0, 0.8}, {40.0, -0.5}, {50.0, 0.0}};
    const Ring lower = cat({{{0.0, -40.0}, {60.0, -40.0}, {60.0, 0.0}}, reversed(inner), {{0.0, 0.0}}});
    const Ring upper = cat({{{0.0, 0.0}}, inner, {{60.0, 0.0}, {60.0, 40.0}, {0.0, 40.0}}});
    const Rings in{lower, upper};
    const BorderOutcome out = run(in, 2.0);
    check(in, out, 2.0, 1e-12);
    CHECK(out.counts.collapses > 0);
}

TEST_CASE("9. refusals: InvalidBand, BadRings, empty output", "[vector_simplify][border_collapse][status]") {
    const Flat ok = flatten(two_squares());
    const auto refused = [](const BorderOutcome& out, BorderStatus status) {
        CHECK(out.status == status);
        CHECK(out.points.empty());
        CHECK(out.ring_starts.empty());
    };
    for (const double band : {-1.0, nan_v, inf_v, -inf_v})
        refused(run_flat(ok.points, ok.starts, band), BorderStatus::InvalidBand);

    const std::vector<Point2> square{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    // Starts not increasing.
    refused(run_flat(square, {0, 4, 2, 4}, 1.0), BorderStatus::BadRings);
    refused(run_flat(square, {0, 2, 2, 4}, 1.0), BorderStatus::BadRings);
    // A ring under three vertices.
    refused(run_flat(square, {0, 2, 4}, 1.0), BorderStatus::BadRings);
    refused(run_flat(square, {0, 1, 4}, 1.0), BorderStatus::BadRings);
    // The last start is not the point count.
    refused(run_flat(square, {0, 3}, 1.0), BorderStatus::BadRings);
    refused(run_flat(square, {0, 5}, 1.0), BorderStatus::BadRings);
    // PINNED: no starts at all.
    refused(run_flat(square, {}, 1.0), BorderStatus::BadRings);
    // A non-finite coordinate.
    for (const double bad : {nan_v, inf_v, -inf_v}) {
        std::vector<Point2> p = ok.points;
        p[3].x = bad;
        refused(run_flat(p, ok.starts, 1.0), BorderStatus::BadRings);
        p = ok.points;
        p[7].y = bad;
        refused(run_flat(p, ok.starts, 1.0), BorderStatus::BadRings);
    }
    // A good input is Ok; PINNED: zero rings are Ok, with zero rings out.
    CHECK(run_flat(ok.points, ok.starts, 1.0).status == BorderStatus::Ok);
    const BorderOutcome none = run_flat({}, {0}, 1.0);
    CHECK(none.status == BorderStatus::Ok);
    CHECK(none.points.empty());
    CHECK(none.ring_starts == std::vector<std::uint64_t>{0});
}

TEST_CASE("10. band 0 is the input bit for bit; the same input twice gives the same bits",
          "[vector_simplify][border_collapse][determinism]") {
    // A collinear run and -0.0 coordinates, both of which a band above 0
    // would act on.
    Rings in = two_squares({-100.0, 0.0});
    in[0].insert(in[0].begin() + 1, Point2{-50.0, 0.0});
    in[1][0] = Point2{-0.0, 0.0};
    const Flat f = flatten(in);
    const BorderOutcome zero = run_flat(f.points, f.starts, 0.0);
    REQUIRE(zero.status == BorderStatus::Ok);
    CHECK(same_bits(zero.points, f.points));
    CHECK(zero.ring_starts == f.starts);
    CHECK(zero.counts.collapses == 0);

    const Flat o = flatten(overshoot_coverage());
    const BorderOutcome a = run_flat(o.points, o.starts, 50.0), b = run_flat(o.points, o.starts, 50.0);
    REQUIRE(a.status == BorderStatus::Ok);
    CHECK(same_bits(a.points, b.points));
    CHECK(a.ring_starts == b.ring_starts);
    CHECK(a.counts.collapses == b.counts.collapses);
}

TEST_CASE("10. pure: four threads at once give the first run's bits",
          "[vector_simplify][border_collapse][determinism][threads]") {
    const Flat o = flatten(overshoot_coverage());
    const BorderOutcome first = run_flat(o.points, o.starts, 50.0);
    std::vector<BorderOutcome> outs(4);
    {
        std::vector<std::jthread> threads;
        for (std::size_t t = 0; t < outs.size(); ++t)
            threads.emplace_back([&, t] { outs[t] = run_flat(o.points, o.starts, 50.0); });
    }
    for (const BorderOutcome& out : outs) {
        CHECK(same_bits(out.points, first.points));
        CHECK(out.ring_starts == first.ring_starts);
    }
}

TEST_CASE("11. a collinear run inside a border is dropped, never at a junction; nearly collinear too",
          "[vector_simplify][border_collapse][degenerate][collinear]") {
    // The shared border runs straight along x = 100 from (100, 0) through
    // five exactly collinear vertices to the junction K = (100, 50), where a
    // third polygon's border leaves at right angles, then on through K
    // (collinear with its neighbours) as a run of vertices 1e-9 m off the
    // line to (100, 100).
    const Point2 K{100.0, 50.0};
    Ring run_a, run_b;
    for (int k = 1; k <= 5; ++k)
        run_a.push_back({100.0, 8.0 * k});
    for (int k = 1; k <= 5; ++k)
        run_b.push_back({100.0 + ((k % 2 == 1) ? 1e-9 : -1e-9), 50.0 + 8.0 * k});
    const Ring left = cat({{{0.0, 0.0}, {100.0, 0.0}}, run_a, {K}, run_b, {{100.0, 100.0}, {0.0, 100.0}}});
    const Ring low_right = cat({{{100.0, 0.0}, {200.0, 0.0}, {200.0, 50.0}, K}, reversed(run_a)});
    const Ring high_right = cat({{K, {200.0, 50.0}, {200.0, 100.0}, {100.0, 100.0}}, reversed(run_b)});
    const Rings in{left, low_right, high_right};
    const BorderOutcome out = run(in, 1.0);
    const Rings got = check(in, out, 1.0);
    CHECK(out.counts.collinear >= 5);
    for (const Point2& p : run_a)
        CHECK_FALSE(has(got[0], p));
    CHECK(has(got[0], K));
    CHECK(has(got[1], K));
    CHECK(has(got[2], K));
}

TEST_CASE("11. two- and three-vertex borders stay; a lens keeps three vertices a side",
          "[vector_simplify][border_collapse][degenerate][lens]") {
    // A 100 m square cut at y ~ 50 into T and B, with a lens L between the
    // junctions (20, 50) and (80, 50): its upper side borders T, its lower
    // side B, seven vertices each. T|B is (0, 50)-(10, 51)-(20, 50), three
    // vertices, and (80, 50)-(100, 50), two.
    const Point2 a{20.0, 50.0}, b{80.0, 50.0}, mid{10.0, 51.0};
    const Ring top_arc{{30.0, 55.0}, {40.0, 58.0}, {50.0, 59.0}, {60.0, 58.0}, {70.0, 55.0}};
    const Ring low_arc{{30.0, 45.0}, {40.0, 42.0}, {50.0, 41.0}, {60.0, 42.0}, {70.0, 45.0}};
    const Ring T = cat({{{0.0, 50.0}, mid, a}, top_arc, {b, {100.0, 50.0}, {100.0, 100.0}, {0.0, 100.0}}});
    const Ring Bt = cat({{{0.0, 0.0}, {100.0, 0.0}, {100.0, 50.0}, b}, reversed(low_arc), {a, mid, {0.0, 50.0}}});
    const Ring L = cat({{a}, low_arc, {b}, reversed(top_arc)});
    const Rings in{T, Bt, L};
    const BorderOutcome out = run(in, 30.0);
    const Rings got = check(in, out, 30.0);
    CHECK(out.counts.collapses > 0);
    CHECK(got[2].size() >= 4); // a, b and at least one vertex a side
    CHECK(has(got[0], mid));
    CHECK(has(got[1], mid));
    const std::vector<Ring> lens = borders_of(got[2], topology(in).junctions);
    REQUIRE(lens.size() == 2);
    for (const Ring& side : lens)
        CHECK(side.size() >= 3);
}

TEST_CASE("11. an edge used by three rings is fixed", "[vector_simplify][border_collapse][degenerate][broken]") {
    // A broken coverage: a triangle R on one edge u-v of test 1's border,
    // overlapping the left square. u-v is used three times: fixed. R's other
    // edges are used once: fixed. The rest of the border still simplifies.
    Rings in = two_squares();
    const Point2 u = in[0][20], v = in[0][21];
    in.push_back({v, u, Point2{60.0, (u.y + v.y) / 2.0}});
    if (area(in[2]) < 0.0)
        in[2] = reversed(in[2]);
    const BorderOutcome out = run(in, 10.0);
    const Rings got = check(in, out, 10.0, 1e-9, -1.0, false);
    CHECK(same_cycle(got[2], in[2]));
    for (const Ring& r : got) {
        bool found = false;
        for (std::size_t j = 0; j < r.size() && !found; ++j) {
            const Point2& p = r[j];
            const Point2& q = r[(j + 1) % r.size()];
            found = (same_bits(p, u) && same_bits(q, v)) || (same_bits(p, v) && same_bits(q, u));
        }
        CHECK(found);
    }
    CHECK(out.counts.collapses > 0);
}

TEST_CASE("11. two rings 1e-9 m apart on a shared border: its edges fixed, nothing moves",
          "[vector_simplify][border_collapse][degenerate][unmatched]") {
    Rings in = two_squares();
    for (Point2& p : in[1])
        if (p.x != 100.0 && p.x != 200.0) // the zig-zag's inner vertices only
            p.x += 1e-9;
    const BorderOutcome out = run(in, 10.0);
    const Rings got = check(in, out, 10.0, 1e-9, -1.0, false);
    CHECK(out.counts.collapses == 0);
    CHECK(same_cycle(got[0], in[0]));
    CHECK(same_cycle(got[1], in[1]));
}

TEST_CASE("11. -0.0 and 0.0 are the same vertex", "[vector_simplify][border_collapse][degenerate][zero]") {
    // Two squares sharing a border along x = 0 from (0, 0) to (0, 100) in
    // 2 m steps, x = 3, 0, -3, 0, ...: every other inner vertex and both
    // junctions have x exactly 0, which the right ring writes as -0.0.
    Ring border;
    for (int k = 1; k < 50; ++k)
        border.push_back({k % 4 == 1 ? 3.0 : k % 4 == 3 ? -3.0 : 0.0, 2.0 * k});
    const Ring left = cat({{{-100.0, 0.0}, {0.0, 0.0}}, border, {{0.0, 100.0}, {-100.0, 100.0}}});
    Ring right = cat({{{0.0, 0.0}, {100.0, 0.0}, {100.0, 100.0}, {0.0, 100.0}}, reversed(border)});
    for (Point2& p : right)
        if (p.x == 0.0)
            p.x = -0.0;
    const Rings in{left, right};
    const BorderOutcome out = run(in, 10.0);
    check(in, out, 10.0);
    CHECK(out.counts.collapses > 0);
}
