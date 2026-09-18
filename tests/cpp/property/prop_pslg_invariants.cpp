// Property tests for the guarantees of terrain::Pslg.
//
// INVARIANT-CRITICAL, alongside test_pslg_builder.cpp. The unit suites pin
// named configurations; this file asserts the laws that must hold for EVERY
// valid Pslg, over generated constraint sets -- and, for the exhaustiveness
// law, over generated CORRUPTIONS of them.
//
// The guarantees are numbered as the increment numbers them, and each named
// property below cites the one it enforces. That numbering is the contract a
// downstream module reads instead of re-validating, so a property that drifts
// away from its guarantee silently widens what the noder and the CDT wrapper
// have to check for themselves.
//
// The NEGATIVE guarantees are here too and are equally load-bearing: a valid
// Pslg promises no simplicity, no pairwise disjointness and no nesting. A test
// asserting one of those would be testing the noder and belongs in its catalog,
// so what this file asserts instead is that inputs violating them are ACCEPTED.
//
// No rapidcheck. Catch2's GENERATE over a fixed seed range plus a seeded
// std::mt19937_64 gives reproducible input without a second framework; a
// failure prints its seed as the generator index and rerunning that section
// reproduces it exactly.
//
// Single kernel throughout. Pslg is kernel-free data and the only
// kernel-parameterised code is the winding stage, whose dependence on K is
// pinned by the one TEMPLATE_TEST_CASE in the builder suite.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <point_families.hpp>
#include <pslg_cases.hpp>
#include <ring_cases.hpp>

#include <terrain/core/edge_properties.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <random>
#include <span>
#include <thread>
#include <vector>

using terrain::Chain;
using terrain::ChainRole;
using terrain::EdgeProperties;
using terrain::IndexedRing;
using terrain::Point2;
using terrain::Pslg;
using terrain::PslgBuilder;
using terrain::PslgBuildResult;
using terrain::PslgError;
using terrain::Segment2;
using terrain::is_closed;
using terrain::orientation;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;
using terrain::test::ChainSpec;
using terrain::test::builder_from;
using terrain::test::concatenated;
using terrain::test::flipped;
using terrain::test::points;
using terrain::test::render;
using terrain::test::utm33_offset;
using terrain::test::valid_chain_specs;

namespace {

constexpr int seed_count = 24;

[[nodiscard]] std::mt19937_64 seeded(int seed) {
    return std::mt19937_64{0x9513'0000ULL + static_cast<std::uint64_t>(seed)};
}

[[nodiscard]] std::vector<ChainSpec> generated_specs(std::mt19937_64& rng) {
    std::uniform_int_distribution<std::size_t> holes{0, 4};
    std::uniform_int_distribution<std::size_t> lines{0, 3};
    return valid_chain_specs(rng, holes(rng), lines(rng));
}

// Every property that needs a valid Pslg goes through this, so "the generator
// produced something the validator rejects" fails loudly at the generator
// rather than as a confusing downstream assertion.
[[nodiscard]] Pslg built(const std::vector<ChainSpec>& specs) {
    PslgBuildResult r = builder_from(specs).build<DefaultKernel>();
    INFO(render(r));
    REQUIRE(r.ok());
    return std::move(*r.pslg);
}

// The indices of the closed chains, which several properties iterate.
[[nodiscard]] std::vector<std::size_t> closed_chains(const Pslg& p) {
    std::vector<std::size_t> out;
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        if (is_closed(p.chains()[c].role)) out.push_back(c);
    }
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// Guarantees 1-9, each named
// ---------------------------------------------------------------------------

// Guarantee 1: every coordinate in vertices() is finite, INCLUDING
// unreferenced ones -- the CDT wrapper hands detria's setPoints the entire
// point array.
// Guarantee 2: every value in chain_indices() is < vertices().size(). This is
// the guarantee that makes IndexedRing::vertex(i)'s unchecked read safe, and
// increment 2 deferred it here explicitly.
TEST_CASE("guarantees 1 and 2: finite vertices, in-range indices", "[pslg][property]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    const Pslg p = built(generated_specs(rng));

    for (const Point2& v : p.vertices()) {
        REQUIRE(std::isfinite(v.x));
        REQUIRE(std::isfinite(v.y));
    }
    for (const std::uint32_t i : p.chain_indices()) {
        REQUIRE(i < p.vertices().size());
    }
}

// Guarantee 3: chains() is non-empty and at least one chain has role Outer.
// Guarantee 4: every closed chain has count >= 3; every breakline has
// count >= 2.
TEST_CASE("guarantees 3 and 4: a bounded domain, and no chain below its minimum",
          "[pslg][property]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    const Pslg p = built(generated_specs(rng));

    REQUIRE_FALSE(p.chains().empty());
    std::size_t outers = 0;
    for (const Chain& c : p.chains()) {
        if (c.role == ChainRole::Outer) ++outers;
        REQUIRE(c.count >= (is_closed(c.role) ? 3u : 2u));
    }
    REQUIRE(outers >= 1);
}

// Guarantee 5: no closed chain stores its closure -- its first and last
// vertices are distinct POINTS, not merely distinct indices.
TEST_CASE("guarantee 5: closure is implied, never stored", "[pslg][property]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    const Pslg p = built(generated_specs(rng));

    for (const std::size_t c : closed_chains(p)) {
        const std::span<const std::uint32_t> idx = p.indices_of(c);
        REQUIRE_FALSE(p.vertices()[idx.front()] == p.vertices()[idx.back()]);
    }
}

// Guarantee 6: every Outer ring is counterclockwise and every Hole ring is
// clockwise under the kernel that validated it; neither is collinear. Asserted
// under the same kernel the build used, which is the whole content of "valid
// with respect to a kernel".
TEST_CASE("guarantee 6: declared role and observed winding agree", "[pslg][property]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    const Pslg p = built(generated_specs(rng));

    for (const std::size_t c : closed_chains(p)) {
        const Orientation o = orientation<DefaultKernel>(p.ring(c));
        REQUIRE(o != Orientation::Collinear);
        REQUIRE(o == (p.chains()[c].role == ChainRole::Outer ? Orientation::CounterClockwise
                                                             : Orientation::Clockwise));
    }
}

// Guarantee 7: the sub-spans indices_of(0..n) partition chain_indices() in
// chain order, contiguously, with no gap and no overlap. A consumer may walk
// the flat buffer once.
TEST_CASE("guarantee 7: indices_of partitions the flat buffer exactly", "[pslg][property]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    const Pslg p = built(generated_specs(rng));

    REQUIRE(p.chains().front().begin == 0u);
    std::size_t cursor = 0;
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        const Chain& ch = p.chains()[c];
        INFO(std::format("chain {} begin {} count {}", c, ch.begin, ch.count));
        REQUIRE(ch.begin == cursor);
        const std::span<const std::uint32_t> sub = p.indices_of(c);
        REQUIRE(sub.size() == ch.count);
        REQUIRE(sub.data() == p.chain_indices().data() + cursor);
        cursor += ch.count;
    }
    REQUIRE(cursor == p.chain_indices().size());
}

// Guarantee 8: ring(c) never throws for a closed chain. IndexedRing's
// constructor rechecks size, closure and the two boundary indices; the
// validator established all three, so every one of those checks is provably
// unreachable here. This property is what turns "provably" into "observed".
TEST_CASE("guarantee 8: ring(c) never throws for any closed chain", "[pslg][property]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    const Pslg p = built(generated_specs(rng));

    for (const std::size_t c : closed_chains(p)) {
        REQUIRE_NOTHROW(p.ring(c));
        const IndexedRing r = p.ring(c);
        REQUIRE(r.size() == p.chains()[c].count);
        for (std::size_t i = 0; i < r.size(); ++i) {
            REQUIRE(&r.vertex(i) == &p.vertices()[p.indices_of(c)[i]]);
        }
    }
}

// Guarantee 9: the vertex buffer is element-wise equal to what the builder
// accumulated. No dedup, no reordering, no reversal, no insertion. A caller's
// index k means the same point going out as it did going in, which is what lets
// Python hold a parallel attribute array.
TEST_CASE("guarantee 9: the vertex buffer is the builder's input verbatim",
          "[pslg][property]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    const std::vector<ChainSpec> specs = generated_specs(rng);
    const std::vector<Point2> expected = concatenated(specs);
    const Pslg p = built(specs);

    REQUIRE(p.vertices().size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        REQUIRE(p.vertices()[i] == expected[i]);
    }
    // And the index runs the point-taking overload emitted are consecutive.
    std::uint32_t next = 0;
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        for (const std::uint32_t i : p.indices_of(c)) REQUIRE(i == next++);
    }
}

// ---------------------------------------------------------------------------
// The derived accessors
// ---------------------------------------------------------------------------

// edge(c, edge_count(c) - 1) closes a ring and does not close a breakline.
// The whole reason this accessor exists is so that neither the noder's broad
// phase nor the CDT wrapper's setConstrainedEdge loop writes (k + 1) % count
// and gets the open case wrong.
TEST_CASE("the last edge closes a ring and only a ring", "[pslg][property][edge]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    const Pslg p = built(generated_specs(rng));

    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        const Chain& ch = p.chains()[c];
        const std::span<const std::uint32_t> idx = p.indices_of(c);
        INFO(std::format("chain {} count {} closed {}", c, ch.count, is_closed(ch.role)));

        REQUIRE(p.edge_count(c) == (is_closed(ch.role) ? ch.count : ch.count - 1u));
        for (std::size_t k = 0; k + 1 < idx.size(); ++k) {
            REQUIRE(p.edge(c, k) ==
                    Segment2{p.vertices()[idx[k]], p.vertices()[idx[k + 1]]});
        }

        const Segment2 last = p.edge(c, p.edge_count(c) - 1);
        if (is_closed(ch.role)) {
            REQUIRE(last == Segment2{p.vertices()[idx.back()], p.vertices()[idx.front()]});
        } else {
            REQUIRE(last ==
                    Segment2{p.vertices()[idx[idx.size() - 2]], p.vertices()[idx.back()]});
        }
    }
}

// ---------------------------------------------------------------------------
// Value semantics
// ---------------------------------------------------------------------------

// Pslg stores no span and no iterator into itself, which is what makes the
// defaulted copy correct. A cached IndexedRing member would make a copied Pslg
// point at the ORIGINAL's buffers -- a bug that survives every test that never
// copies, and survives a value-only comparison after one that does.
TEST_CASE("a copy compares equal element-wise and points into its own buffers",
          "[pslg][property][value_semantics]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    const Pslg original = built(generated_specs(rng));
    const Pslg copy = original;

    REQUIRE(copy.vertices().data() != original.vertices().data());
    REQUIRE(copy.chain_indices().data() != original.chain_indices().data());
    REQUIRE(copy.chains().data() != original.chains().data());

    REQUIRE(copy.vertices().size() == original.vertices().size());
    for (std::size_t i = 0; i < copy.vertices().size(); ++i)
        REQUIRE(copy.vertices()[i] == original.vertices()[i]);
    REQUIRE(copy.chains().size() == original.chains().size());
    for (std::size_t i = 0; i < copy.chains().size(); ++i)
        REQUIRE(copy.chains()[i] == original.chains()[i]);

    for (const std::size_t c : closed_chains(copy)) {
        const IndexedRing r = copy.ring(c);
        REQUIRE(&r.vertex(0) == &copy.vertices()[copy.indices_of(c)[0]]);
        REQUIRE(&r.vertex(0) != &original.vertices()[original.indices_of(c)[0]]);
    }
}

// Guarantee 10: a const Pslg is safe for concurrent read from any number of
// threads. There is no mutator, no non-const accessor, no lazy cache and no
// mutable member, so this needs no synchronisation -- which is a requirement of
// parallel_refinement.md, where every refinement thread reads the same
// constraint set, and not an aspiration.
//
// The assertion is agreement, not absence of a crash: every thread computes the
// same reduction over the same const Pslg and all of them must produce the same
// answer. Under tsan this is also where a cache added to an accessor shows up.
TEST_CASE("guarantee 10: concurrent readers of a const Pslg agree",
          "[pslg][property][concurrency]") {
    std::mt19937_64 rng = seeded(0);
    const Pslg p = built(valid_chain_specs(rng, 3, 2));

    const auto reduce = [](const Pslg& q) {
        double acc = 0.0;
        for (std::size_t c = 0; c < q.chains().size(); ++c) {
            for (std::size_t k = 0; k < q.edge_count(c); ++k) {
                const Segment2 s = q.edge(c, k);
                acc += s.a.x - s.b.x + 3.0 * (s.a.y - s.b.y);
            }
            if (is_closed(q.chains()[c].role)) {
                acc += orientation<DefaultKernel>(q.ring(c)) == Orientation::CounterClockwise
                           ? 1.0
                           : -1.0;
            }
        }
        return acc;
    };

    const double expected = reduce(p);
    constexpr std::size_t thread_count = 8;
    std::vector<double> results(thread_count, 0.0);
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    for (std::size_t t = 0; t < thread_count; ++t) {
        threads.emplace_back([&p, &results, t, &reduce] { results[t] = reduce(p); });
    }
    for (std::thread& th : threads) th.join();

    for (std::size_t t = 0; t < thread_count; ++t) REQUIRE(results[t] == expected);
}

// ---------------------------------------------------------------------------
// Exhaustiveness, over generated corruptions
// ---------------------------------------------------------------------------

// The named property behind mutant 8: N INDEPENDENTLY BROKEN CHAINS PRODUCE AT
// LEAST N DIAGNOSTICS, IN STAGES 1-6. The unit suite pins three hand-written
// cases; this one varies which chains are broken, how many, and by which
// corruption, so an early return that happens to survive the hand-written
// arrangement does not survive this.
//
// Stages 1-6, not stage 0: the size-overflow check early-returns by
// specification, because once begin/count have been truncated every later
// diagnostic is noise attributed to chains that may be well-formed. Nothing
// here reaches it -- every generated constraint set is orders of magnitude
// inside uint32 -- so this property neither pins nor forbids that behaviour.
//
// "At least" rather than "exactly": one corruption may legitimately raise more
// than one diagnostic -- a chain with two out-of-range indices raises two -- and
// a property that demanded exactly N would be pinning the diagnosis's
// granularity rather than its completeness.
TEST_CASE("N independently broken chains produce at least N diagnostics",
          "[pslg][property][exhaustive]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);

    std::vector<ChainSpec> specs = valid_chain_specs(rng, 4, 2);
    std::vector<std::size_t> closed;
    for (std::size_t i = 0; i < specs.size(); ++i) {
        if (is_closed(specs[i].role)) closed.push_back(i);
    }
    REQUIRE(closed.size() >= 2);

    // Break every closed chain, each in a way independent of the others.
    std::uniform_int_distribution<int> how{0, 2};
    std::size_t broken = 0;
    for (const std::size_t i : closed) {
        switch (how(rng)) {
            case 0:  // reverse the winding
                specs[i].pts = flipped(points(specs[i].pts));
                break;
            case 1:  // store the closure
                specs[i].pts.push_back(specs[i].pts.front());
                break;
            default:  // drop below the minimum vertex count
                specs[i].pts.resize(2);
                break;
        }
        ++broken;
    }

    const PslgBuildResult r = builder_from(specs).build<DefaultKernel>();
    INFO(render(r));
    REQUIRE_FALSE(r.ok());
    REQUIRE_FALSE(r.pslg.has_value());
    REQUIRE(r.diagnostics.size() >= broken);

    // And every broken chain is named by at least one diagnostic. A validator
    // that reported N diagnostics all about the first chain would satisfy the
    // count and still be useless.
    for (const std::size_t i : closed) {
        const bool named = std::any_of(r.diagnostics.begin(), r.diagnostics.end(),
                                       [i](const terrain::PslgDiagnostic& d) {
                                           return d.chain == static_cast<std::uint32_t>(i);
                                       });
        INFO(std::format("chain {} unnamed", i));
        REQUIRE(named);
    }
}

// A non-finite vertex anywhere in the buffer rejects the whole set, wherever it
// sits and whether or not a chain names it.
TEST_CASE("a non-finite vertex anywhere rejects the build", "[pslg][property][exhaustive]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);

    const std::vector<ChainSpec> specs = valid_chain_specs(rng, 2, 1);
    std::vector<Point2> verts = concatenated(specs);
    const std::size_t n = verts.size();

    // One extra vertex that no chain will ever name.
    verts.push_back(Point2{1.0, 1.0});
    std::uniform_int_distribution<std::size_t> which{0, verts.size() - 1};
    const std::size_t victim = which(rng);
    const bool x_axis = (rng() & 1u) != 0;
    const double bad = (rng() & 1u) != 0 ? std::numeric_limits<double>::quiet_NaN()
                                         : std::numeric_limits<double>::infinity();
    if (x_axis) {
        verts[victim].x = bad;
    } else {
        verts[victim].y = bad;
    }

    PslgBuilder b{verts};
    std::uint32_t cursor = 0;
    for (const ChainSpec& s : specs) {
        std::vector<std::uint32_t> idx;
        for (std::size_t i = 0; i < s.pts.size(); ++i) idx.push_back(cursor++);
        b.add_chain(std::span<const std::uint32_t>{idx}, s.role, s.properties);
    }
    REQUIRE(cursor == n);

    const PslgBuildResult r = std::move(b).build<DefaultKernel>();
    INFO(std::format("victim {} of {}\n{}", victim, verts.size(), render(r)));
    REQUIRE_FALSE(r.ok());
    const auto* d = terrain::test::find_error(r, PslgError::NonFiniteVertex);
    REQUIRE(d != nullptr);
    REQUIRE(d->vertex == static_cast<std::uint32_t>(victim));

    // And nothing else. A chain carrying a non-finite vertex is excluded from
    // the winding stage, so no geometric diagnostic may be raised about it --
    // orient2d on a NaN coordinate returns Collinear, and reporting that as
    // DegenerateRing sends the reader to look for collinear geometry that is
    // not there. The generated rings are as small as three vertices, where
    // orientation<K>'s independent prev/next walk cannot route around the bad
    // vertex and the wrong order is forced into the open.
    REQUIRE(terrain::test::count_errors(r, PslgError::DegenerateRing) == 0);
    REQUIRE(terrain::test::count_errors(r, PslgError::WrongWinding) == 0);
}

// ---------------------------------------------------------------------------
// The negative guarantees
// ---------------------------------------------------------------------------

// A valid Pslg promises NO simplicity, NO pairwise disjointness and NO nesting.
// Stated as properties, that means: adding geometry which violates any of them
// to an otherwise valid set does not change the build outcome. If one of these
// ever starts failing, someone has added a check, and the noder's contract --
// and the "rejects legitimate input" failure mode the design argues against --
// moved without anybody deciding to move it.
TEST_CASE("crossing, touching and un-nested geometry does not change the outcome",
          "[pslg][property][negative]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);

    std::vector<ChainSpec> specs = valid_chain_specs(rng, 2, 1);
    REQUIRE(builder_from(specs).build<DefaultKernel>().ok());

    // A hole far outside every outer ring: inside no ring at all.
    specs.push_back(ChainSpec{
        std::vector<Point2>{Point2{10000.0, 10000.0}, Point2{10000.0, 10010.0},
                            Point2{10010.0, 10010.0}},
        ChainRole::Hole, EdgeProperties{}});
    // A second outer ring overlapping the first, so the hole added next lies
    // inside two of them.
    specs.push_back(ChainSpec{
        std::vector<Point2>{Point2{-50.0, -50.0}, Point2{50.0, -50.0}, Point2{50.0, 50.0},
                            Point2{-50.0, 50.0}},
        ChainRole::Outer, EdgeProperties{}});
    // A breakline that crosses everything, including itself.
    specs.push_back(ChainSpec{
        std::vector<Point2>{Point2{-200.0, -200.0}, Point2{200.0, 200.0},
                            Point2{200.0, -200.0}, Point2{-200.0, 200.0}},
        ChainRole::Breakline, EdgeProperties{}});

    const PslgBuildResult r = builder_from(specs).build<DefaultKernel>();
    INFO(render(r));
    REQUIRE(r.ok());
}

// The same statement at UTM33 magnitudes, where a coordinate carries eight
// significant digits before the decimal point and any predicate that reaches
// for a tolerance stops working. Nothing in the validator constructs a point,
// so translating a whole constraint set by an exactly representable UTM33
// offset cannot change a build outcome.
TEST_CASE("a UTM33 translation does not change a build outcome", "[pslg][property][utm33]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);

    std::vector<ChainSpec> local = valid_chain_specs(rng, 3, 2);
    std::vector<ChainSpec> shifted = local;
    for (ChainSpec& s : shifted) {
        for (Point2& v : s.pts) v = utm33_offset(v);
    }

    const PslgBuildResult a = builder_from(local).build<DefaultKernel>();
    const PslgBuildResult b = builder_from(shifted).build<DefaultKernel>();
    INFO(render(a));
    INFO(render(b));
    REQUIRE(a.ok());
    REQUIRE(b.ok() == a.ok());

    // And the windings are the same rings, not merely the same verdict.
    REQUIRE(a.pslg->chains().size() == b.pslg->chains().size());
    for (std::size_t c = 0; c < a.pslg->chains().size(); ++c) {
        if (!is_closed(a.pslg->chains()[c].role)) continue;
        REQUIRE(orientation<DefaultKernel>(a.pslg->ring(c)) ==
                orientation<DefaultKernel>(b.pslg->ring(c)));
    }
}
