// Property tests for terrain/noding/node.hpp -- the noder's guarantees.
//
// INVARIANT-CRITICAL, mutation round (docs/increments/05b-noder-driver.md,
// "What is worth testing"). The file name is the one testing.md reserves;
// `grep -n prop_noding_no_crossings testing.md` resolves it, and the line it
// sits on says "increment 5b: guarantees 14 and 15".
//
// NO BROAD PHASE ANYWHERE IN THE CHECKS. Guarantee 14 is checked brute force
// over every (edge, edge) and every (node, edge) pair. A broad-phase defect and
// a matching verification blind spot cannot cancel against an oracle that never
// consults the index. O(n^2) is why the generated sets are small, and why this
// is a property rather than the production path.
//
// ---------------------------------------------------------------------------
// GUARANTEE 15'S ORACLE, AND THE TWO WAYS IT HAS ALREADY BEEN GOT WRONG.
//
// It is built FROM THE INPUT and never from the dedup's provenance map. An
// oracle that asks the noder which chains it recorded as contributing, and then
// asserts the union over those, is the dedup restating itself: green on any
// internally consistent merge, including one that silently drops a contributor.
// Borrow the producer's PREDICATE, never its RECORDS
// (.claude/skills/computational-geometry/SKILL.md; the worked history is
// docs/increments/05-noder.md, guarantees 14 and 15).
//
// The predicate borrowed is segment_meets_cell<K> plus arc-order betweenness --
// PROXIMITY, not incidence. on_segment appears nowhere in this file. world(g)
// is not affine at a non-dyadic spacing, so a snapped node is NEAR a segment and
// not ON it; an exact-incidence oracle over snapped output is red on correct
// output, on 207 of 2976 measured cases at 0.1 m.
//
// And the assertion is a SUBSET, in one direction only. The converse is FALSE,
// not merely unproven: a chain passing near both nodes and spanning them
// satisfies the relation without having contributed geometry, so the output may
// legitimately be a strict superset of what the oracle can prove. An equality
// assertion would be red on correct output. A spuriously ADDED property escapes
// this assertion by construction and is owed a mutant instead -- mutant 7, which
// lives in tests/cpp/unit/test_noding_node.cpp because it needs a fixture whose
// two contributors are known by name.
// ---------------------------------------------------------------------------
//
// No rapidcheck; Catch2 GENERATE over a seeded range, as elsewhere in this tree.
//
// FastKernel IS FORBIDDEN IN THIS FILE. The oracle is built from the kernel, so
// an oracle less exact than its subject reports false failures on precisely the
// degenerate inputs the property is about -- and every input here is snapped by
// construction, which is where the filter's fall-throughs live (38 % of
// grid-collinear triples at 0.1 m). DefaultKernel only.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "noding_generators.h"

#include <terrain/core/edge_properties.hpp>
#include <terrain/core/noded_pslg.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/noding/node.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

using terrain::Chain;
using terrain::EdgeProperties;
using terrain::GridPoint;
using terrain::NodedPslg;
using terrain::Point2;
using terrain::Pslg;
using terrain::Segment2;
using terrain::SnapGrid;
using terrain::noding::NodeOptions;
using terrain::noding::NodeOutcome;
using terrain::noding::NodeStatus;
using terrain::noding::SegmentRelation;
using terrain::noding::classify;
using terrain::noding::node;
using terrain::noding::segment_meets_cell;
using terrain::pred::DefaultKernel;
using terrain::testing::PolylineOptions;
using terrain::testing::PolylineSet;
using terrain::testing::random_polylines;
using terrain::testing::shuffled;
using terrain::testing::to_pslg;

namespace {

constexpr int seed_count = 24;
constexpr double kSpacing = 0.1;

// The two arms of the density dial. testing.md:220 calls the generator's
// subject "controllable density of intersections" and both ends of the control
// are worth running: the dense arm exercises the split pass and the cascade, the
// sparse arm exercises the path where the noder must change nothing at all and
// where a driver that splits something anyway is visible.
[[nodiscard]] PolylineOptions dense() { return PolylineOptions{}; }

[[nodiscard]] PolylineOptions sparse() {
    PolylineOptions o;
    o.reach = 5.0;
    return o;
}

// Every output edge, as a flat (c, k) position with its node id pair.
struct OutputEdge {
    std::size_t c{};
    std::size_t k{};
    std::size_t flat{};
    std::uint32_t lo{};
    std::uint32_t hi{};
    Segment2 segment{};
};

[[nodiscard]] std::vector<OutputEdge> output_edges(const NodedPslg& p) {
    std::vector<OutputEdge> edges;
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        const auto idx = p.indices_of(c);
        const std::size_t n = idx.size();
        for (std::size_t k = 0; k < p.edge_count(c); ++k) {
            const std::uint32_t u = idx[k];
            const std::uint32_t v = idx[(k + 1) % n];
            edges.push_back(OutputEdge{c, k, p.edge_base(c) + k, std::min(u, v),
                                       std::max(u, v), p.edge(c, k)});
        }
    }
    return edges;
}

[[nodiscard]] std::vector<GridPoint> node_grid_points(const NodedPslg& p) {
    std::vector<GridPoint> nodes;
    nodes.reserve(p.vertices().size());
    for (const Point2& v : p.vertices()) {
        nodes.push_back(p.grid().snap(v));
    }
    return nodes;
}

// ---------------------------------------------------------------------------
// Guarantee 15's oracle, built from the input Pslg.
// ---------------------------------------------------------------------------

// The producer's sort key for step 4: the dot product of (world(g) - s.a) with
// (s.b - s.a), ties broken by GridPoint's defaulted, therefore lexicographic,
// operator<=>. Borrowed as a PREDICATE. Nothing here reads a record the noder
// kept.
struct ArcKey {
    double t{};
    GridPoint g{};

    [[nodiscard]] friend bool operator<(const ArcKey& a, const ArcKey& b) {
        return a.t != b.t ? a.t < b.t : a.g < b.g;
    }
};

[[nodiscard]] ArcKey arc_key(const Segment2& s, const SnapGrid& grid, const GridPoint& g) {
    const Point2 d{s.b.x - s.a.x, s.b.y - s.a.y};
    const Point2 w = grid.world(g);
    return ArcKey{(w.x - s.a.x) * d.x + (w.y - s.a.y) * d.y, g};
}

// For one input segment: the node ids whose cells it meets, in arc order. This
// is exactly the sequence step 4 says the segment is split into, reconstructed
// from the input. Precomputed per segment because the alternative is
// O(edges * segments * nodes) calls into the exact kernel.
[[nodiscard]] std::vector<std::uint32_t> split_sequence(const Segment2& s, const SnapGrid& grid,
                                                        const std::vector<GridPoint>& nodes) {
    std::vector<std::pair<ArcKey, std::uint32_t>> met;
    for (std::uint32_t i = 0; i < nodes.size(); ++i) {
        if (segment_meets_cell<DefaultKernel>(grid, s, nodes[i])) {
            met.emplace_back(arc_key(s, grid, nodes[i]), i);
        }
    }
    std::ranges::sort(met, [](const auto& a, const auto& b) { return a.first < b.first; });

    std::vector<std::uint32_t> ids;
    ids.reserve(met.size());
    for (const auto& [key, id] : met) {
        ids.push_back(id);
    }
    return ids;
}

// Does this input segment contribute geometry to the output edge (lo, hi)? Yes
// when both nodes are in its split sequence and nothing it meets lies strictly
// between them in arc order -- i.e. when (lo, hi) is one of the pieces the
// segment would be cut into.
[[nodiscard]] bool contributes(const std::vector<std::uint32_t>& sequence, std::uint32_t lo,
                               std::uint32_t hi) {
    for (std::size_t i = 0; i + 1 < sequence.size(); ++i) {
        const std::uint32_t a = sequence[i];
        const std::uint32_t b = sequence[i + 1];
        if (std::min(a, b) == lo && std::max(a, b) == hi) {
            return true;
        }
    }
    return false;
}

// ---------------------------------------------------------------------------
// The guarantees, each as a free function so that a failure names one.
// ---------------------------------------------------------------------------

void check_guarantee_11_12_13(const NodedPslg& p) {
    for (const Point2& v : p.vertices()) {
        const Point2 again = p.grid().snapped(v);
        REQUIRE(again.x == v.x);
        REQUIRE(again.y == v.y);
    }

    std::vector<GridPoint> nodes = node_grid_points(p);
    std::vector<GridPoint> sorted = nodes;
    std::ranges::sort(sorted);
    REQUIRE(std::ranges::adjacent_find(sorted) == sorted.end());

    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        const auto idx = p.indices_of(c);
        const std::size_t n = idx.size();
        for (std::size_t k = 0; k < p.edge_count(c); ++k) {
            REQUIRE(idx[k] != idx[(k + 1) % n]);
        }
    }
}

void check_guarantee_14(const NodedPslg& p) {
    const std::vector<OutputEdge> edges = output_edges(p);
    const std::vector<GridPoint> nodes = node_grid_points(p);

    // 14(a), amended: Disjoint or Touching, OR Overlapping with EQUAL node-id
    // pairs. A partial overlap is still a violation; after the split pass there
    // is no third case, because each collinear overlap's endpoints lie in the
    // other's hot pixels and so each is split at the other's nodes.
    for (std::size_t i = 0; i < edges.size(); ++i) {
        for (std::size_t j = i + 1; j < edges.size(); ++j) {
            const SegmentRelation r =
                classify<DefaultKernel>(edges[i].segment, edges[j].segment);
            if (r == SegmentRelation::Overlapping) {
                REQUIRE(edges[i].lo == edges[j].lo);
                REQUIRE(edges[i].hi == edges[j].hi);
            } else {
                REQUIRE(r != SegmentRelation::Crossing);
            }
        }
    }

    // 14(b), the hot-pixel clause, and it is STRICTLY STRONGER than exact
    // incidence. That is why it catches a T-junction the split pass missed
    // instead of blessing it, and it is the reason on_segment is not used here.
    for (const OutputEdge& e : edges) {
        for (std::uint32_t g = 0; g < nodes.size(); ++g) {
            if (g == e.lo || g == e.hi) {
                continue;
            }
            REQUIRE_FALSE(segment_meets_cell<DefaultKernel>(p.grid(), e.segment, nodes[g]));
        }
    }
}

void check_guarantee_15(const Pslg& in, const NodedPslg& out) {
    const std::vector<GridPoint> nodes = node_grid_points(out);

    // One split sequence per input segment, and the chain it came from.
    std::vector<std::pair<std::size_t, std::vector<std::uint32_t>>> sequences;
    for (std::size_t c = 0; c < in.chains().size(); ++c) {
        for (std::size_t k = 0; k < in.edge_count(c); ++k) {
            sequences.emplace_back(c, split_sequence(in.edge(c, k), out.grid(), nodes));
        }
    }

    for (const OutputEdge& e : output_edges(out)) {
        EdgeProperties proven;
        for (const auto& [c, sequence] : sequences) {
            if (contributes(sequence, e.lo, e.hi)) {
                proven = proven | in.chains()[c].properties;
            }
        }
        // SUBSET, one direction. See the head of this file for why the converse
        // is false rather than unproven.
        REQUIRE(out.edge_properties()[e.flat].contains(proven));
    }
}

void check_guarantee_16(const Pslg& in, const NodedPslg& out) {
    REQUIRE(out.chains().size() == in.chains().size());
    for (std::size_t c = 0; c < in.chains().size(); ++c) {
        REQUIRE(out.chains()[c].role == in.chains()[c].role);
        REQUIRE(out.chains()[c].properties == in.chains()[c].properties);
    }
}

void check_edge_base(const NodedPslg& p) {
    std::size_t running = 0;
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        REQUIRE(p.edge_base(c) == running);
        running += p.edge_count(c);
    }
    REQUIRE(p.edge_base(p.chains().size()) == running);
    REQUIRE(p.edge_properties().size() == running);
}

void check_everything(const Pslg& in, const NodedPslg& out) {
    check_guarantee_11_12_13(out);
    check_guarantee_14(out);
    check_guarantee_15(in, out);
    check_guarantee_16(in, out);
    check_edge_base(out);
}

}  // namespace

// ---------------------------------------------------------------------------
// The probe has to be able to fail, so the first thing asserted is that the
// generator produces work and that the noder mostly succeeds on it. A noder
// that returned NotConverged on everything would satisfy every guarantee below
// vacuously, because they are all conditioned on an engaged pslg.
// ---------------------------------------------------------------------------

TEST_CASE("the generated sets are noded successfully and are not already noded",
          "[noding][no_crossings]") {
    int ok = 0;
    int split_something = 0;
    for (int seed = 0; seed < seed_count; ++seed) {
        const PolylineSet set = random_polylines(static_cast<std::uint64_t>(seed), dense());
        const Pslg in = to_pslg(set);
        const NodeOutcome out = node<DefaultKernel>(in, NodeOptions{kSpacing});
        if (out.status != NodeStatus::Ok) {
            continue;
        }
        ++ok;

        std::size_t in_edges = 0;
        for (std::size_t c = 0; c < in.chains().size(); ++c) {
            in_edges += in.edge_count(c);
        }
        if (out.pslg->edge_properties().size() > in_edges) {
            ++split_something;
        }
    }

    REQUIRE(ok >= seed_count * 3 / 4);
    REQUIRE(split_something >= seed_count / 2);
}

TEST_CASE("a non-Ok outcome carries no pslg and a message", "[noding][no_crossings]") {
    // MUTANT 6's property-suite half: max_rounds treated as a success, returning
    // the last candidate instead of NotConverged. The fixture built to need more
    // rounds than a cap of 1 is in test_noding_node.cpp; this is the blanket
    // statement over generated input, and it is what makes every guarantee below
    // safe to condition on ok().
    const int seed = GENERATE(range(0, seed_count));
    const PolylineSet set = random_polylines(static_cast<std::uint64_t>(seed), dense());
    const Pslg in = to_pslg(set);

    const NodeOutcome out = node<DefaultKernel>(in, NodeOptions{kSpacing, 1});
    REQUIRE((out.status == NodeStatus::Ok) == out.ok());
    REQUIRE(out.ok() == out.message.empty());
}

// ---------------------------------------------------------------------------
// Guarantees 11 to 16, on both arms of the density dial.
// ---------------------------------------------------------------------------

TEST_CASE("noded output satisfies every guarantee, dense input", "[noding][no_crossings]") {
    const int seed = GENERATE(range(0, seed_count));
    const PolylineSet set = random_polylines(static_cast<std::uint64_t>(seed), dense());
    const Pslg in = to_pslg(set);

    const NodeOutcome out = node<DefaultKernel>(in, NodeOptions{kSpacing});
    if (!out.ok()) {
        // A legitimate refusal. Its own statuses are pinned by name in
        // test_noding_node.cpp; here it only must not pretend to have succeeded.
        REQUIRE(out.status != NodeStatus::Ok);
        return;
    }
    check_everything(in, *out.pslg);
}

TEST_CASE("noded output satisfies every guarantee, sparse input", "[noding][no_crossings]") {
    const int seed = GENERATE(range(0, seed_count));
    const PolylineSet set = random_polylines(static_cast<std::uint64_t>(seed), sparse());
    const Pslg in = to_pslg(set);

    const NodeOutcome out = node<DefaultKernel>(in, NodeOptions{kSpacing});
    if (!out.ok()) {
        REQUIRE(out.status != NodeStatus::Ok);
        return;
    }
    check_everything(in, *out.pslg);
}

TEST_CASE("noded output satisfies every guarantee at a dyadic spacing",
          "[noding][no_crossings]") {
    // 0.0625 is 2^-4, where world() IS affine and snapping does preserve
    // collinearity. Run because the decimal arm is where the surprises are and a
    // suite that only ran there could not tell a rounding defect from a
    // topological one.
    const int seed = GENERATE(range(0, seed_count));
    const PolylineSet set = random_polylines(static_cast<std::uint64_t>(seed), dense());
    const Pslg in = to_pslg(set);

    const NodeOutcome out = node<DefaultKernel>(in, NodeOptions{0.0625});
    if (!out.ok()) {
        REQUIRE(out.status != NodeStatus::Ok);
        return;
    }
    check_everything(in, *out.pslg);
}

// ---------------------------------------------------------------------------
// Determinism. node<K> is a pure function of (pslg, options), and node ids are
// a function of the point SET rather than of the order anything was found in.
// ---------------------------------------------------------------------------

TEST_CASE("the same input noded twice gives bit-identical output", "[noding][no_crossings]") {
    const int seed = GENERATE(range(0, seed_count));
    const PolylineSet set = random_polylines(static_cast<std::uint64_t>(seed), dense());
    const Pslg in = to_pslg(set);

    const NodeOutcome a = node<DefaultKernel>(in, NodeOptions{kSpacing});
    const NodeOutcome b = node<DefaultKernel>(in, NodeOptions{kSpacing});
    REQUIRE(a.status == b.status);
    if (!a.ok()) {
        return;
    }

    REQUIRE(a.pslg->vertices().size() == b.pslg->vertices().size());
    for (std::size_t i = 0; i < a.pslg->vertices().size(); ++i) {
        REQUIRE(a.pslg->vertices()[i].x == b.pslg->vertices()[i].x);
        REQUIRE(a.pslg->vertices()[i].y == b.pslg->vertices()[i].y);
    }
    REQUIRE(std::ranges::equal(a.pslg->chain_indices(), b.pslg->chain_indices()));
}

TEST_CASE("shuffling the input vertex order does not move a single node",
          "[noding][no_crossings]") {
    // Node ids are NodeSet's sorted order, which is a function of the SET of
    // grid points -- not of input order, not of which thread found a crossing
    // first, not of the order the broad phase visited buckets in. Permuting the
    // input vertex buffer must therefore leave the output identical, chain
    // slices included, because chain order is preserved by guarantee 16.
    const int seed = GENERATE(range(0, seed_count));
    const PolylineSet set = random_polylines(static_cast<std::uint64_t>(seed), dense());
    const PolylineSet permuted = shuffled(set, static_cast<std::uint64_t>(seed) + 991u);

    const NodeOutcome a = node<DefaultKernel>(to_pslg(set), NodeOptions{kSpacing});
    const NodeOutcome b = node<DefaultKernel>(to_pslg(permuted), NodeOptions{kSpacing});
    REQUIRE(a.status == b.status);
    if (!a.ok()) {
        return;
    }

    REQUIRE(a.pslg->vertices().size() == b.pslg->vertices().size());
    for (std::size_t i = 0; i < a.pslg->vertices().size(); ++i) {
        REQUIRE(a.pslg->vertices()[i].x == b.pslg->vertices()[i].x);
        REQUIRE(a.pslg->vertices()[i].y == b.pslg->vertices()[i].y);
    }
    REQUIRE(std::ranges::equal(a.pslg->chain_indices(), b.pslg->chain_indices()));

    for (std::size_t e = 0; e < a.pslg->edge_properties().size(); ++e) {
        REQUIRE(a.pslg->edge_properties()[e] == b.pslg->edge_properties()[e]);
    }
}

TEST_CASE("an unreferenced far-away vertex changes nothing", "[noding][no_crossings]") {
    // MUTANT 12: an input-dependent anchor -- subtracting a bounding-box origin
    // before snapping. 05-noder.md assigns this mutant to 5b's driver, since
    // 5a's snap() has no argument an anchor could arrive through. An extra
    // unreferenced vertex a kilometre away moves the input's bounding box and
    // nothing else; the node coordinates must be bit-identical.
    //
    // Unreferenced vertices are legal in a Pslg -- the validator checks them for
    // finiteness precisely because the CDT wrapper hands detria the whole point
    // array (core/pslg_builder.hpp, stage 2).
    const int seed = GENERATE(range(0, seed_count));
    const PolylineSet set = random_polylines(static_cast<std::uint64_t>(seed), dense());

    PolylineSet moved = set;
    moved.vertices.push_back(Point2{100000.0, -100000.0});

    const NodeOutcome a = node<DefaultKernel>(to_pslg(set), NodeOptions{kSpacing});
    const NodeOutcome b = node<DefaultKernel>(to_pslg(moved), NodeOptions{kSpacing});
    REQUIRE(a.status == b.status);
    if (!a.ok()) {
        return;
    }

    // Step 3 builds the NodeSet over EVERY snapped input vertex, so the extra
    // one becomes a node of its own: b has exactly one more. Its grid index is
    // larger on x than anything in the domain, so it sorts last and no existing
    // node id moves. That reading of step 3 is pinned here rather than assumed,
    // because it is the difference between "one more vertex" and "every id
    // shifted", and only the first leaves the comparison below meaningful.
    REQUIRE(b.pslg->vertices().size() == a.pslg->vertices().size() + 1);

    for (std::size_t i = 0; i < a.pslg->vertices().size(); ++i) {
        REQUIRE(a.pslg->vertices()[i].x == b.pslg->vertices()[i].x);
        REQUIRE(a.pslg->vertices()[i].y == b.pslg->vertices()[i].y);
    }
    REQUIRE(std::ranges::equal(a.pslg->chain_indices(), b.pslg->chain_indices()));
}
