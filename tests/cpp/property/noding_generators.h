#pragma once

// Deterministic generators for the noder's property suites -- increment 5b.
//
// testing.md:220 reserves this path and says what lives here: "random sets of
// polylines with controllable density of intersections". It arrives with 5b and
// not with 5a, because 5a's properties are over points, spacings and segment
// PAIRS, and nothing in 5a consumes a polyline set.
//
// SEEDED AND REPRODUCIBLE. Every function here is a pure function of its seed
// and its options; std::mt19937_64 is specified by the standard down to the bit,
// so a failure at seed N reproduces on every platform CI runs.
//
// THE SETS ARE SMALL AND THE SIZE IS FIXED HERE ON PURPOSE.
// prop_noding_no_crossings.cpp checks guarantee 14 by brute force over every
// (edge, edge) and every (node, edge) pair, which is O(n^2) per noding round,
// and CI runs it on two platforms plus an asan+ubsan Debug build of the same
// suites (05b-noder-driver.md, risk 19). A few dozen segments is milliseconds;
// the failure mode the risk names is "a few dozen" drifting upward during the
// round because a bigger generator finds more. If a bigger set is wanted, it
// belongs in a separate, separately-registered suite with its own budget.
//
// This header deliberately knows nothing about NodedPslg, NodeOptions or
// node<K>: it produces INPUT. An oracle built from a generator that also knew
// the output shape would be one step from being built from the output.

#include <terrain/core/edge_properties.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace terrain::testing {

struct PolylineOptions {
    // One counterclockwise outer square is always emitted first, so chain 0 is
    // an Outer chain and the set is a legal Pslg (increment 3, NoOuterChain).
    double extent{100.0};

    std::size_t breakline_count{8};
    std::size_t vertices_per_breakline{3};

    // THE INTERSECTION-DENSITY DIAL. A breakline's first vertex is uniform in
    // the square and each subsequent vertex is a uniform step of up to `reach`
    // from the last. reach == extent gives chords spanning the domain and a
    // crossing between almost every pair; reach == extent / 20 gives short local
    // polylines that mostly miss each other. Both arms matter: the dense arm
    // exercises the split pass and the cascade, the sparse arm exercises the
    // path where the noder must change nothing at all.
    double reach{100.0};

    // How many distinct property bits the breaklines draw from. Three is enough
    // for the union oracle to distinguish "unioned" from "took the first" and
    // from "took its own"; more only slows the suite down.
    unsigned property_bits{3};
};

// The raw polyline set, before it becomes a Pslg. Exposed because the guarantee
// 15 oracle needs the INPUT chains and their property sets, and reading them off
// a Pslg is fine but reading them off the noder's output is the self-confirming
// shape this whole increment is guarding against.
struct PolylineSet {
    std::vector<Point2> vertices;
    std::vector<std::vector<std::uint32_t>> runs;
    std::vector<ChainRole> roles;
    std::vector<EdgeProperties> properties;
};

[[nodiscard]] inline PolylineSet random_polylines(std::uint64_t seed,
                                                  const PolylineOptions& options) {
    std::mt19937_64 rng{seed};
    std::uniform_real_distribution<double> inside{0.0, options.extent};
    std::uniform_real_distribution<double> step{-options.reach, options.reach};
    std::uniform_int_distribution<unsigned> which_bit{0, options.property_bits - 1};
    std::uniform_int_distribution<int> how_many{1, 2};

    PolylineSet set;

    // Chain 0: the outer square, counterclockwise, with a margin so that a
    // breakline stepping outside the square stays a legal Pslg rather than
    // becoming an accident of the generator.
    const double m = options.extent;
    set.vertices = {Point2{-m, -m}, Point2{2.0 * m, -m}, Point2{2.0 * m, 2.0 * m},
                    Point2{-m, 2.0 * m}};
    set.runs.push_back({0, 1, 2, 3});
    set.roles.push_back(ChainRole::Outer);
    set.properties.push_back(EdgeProperties{});

    for (std::size_t c = 0; c < options.breakline_count; ++c) {
        std::vector<std::uint32_t> run;
        Point2 p{inside(rng), inside(rng)};
        for (std::size_t v = 0; v < options.vertices_per_breakline; ++v) {
            run.push_back(static_cast<std::uint32_t>(set.vertices.size()));
            set.vertices.push_back(p);
            p = Point2{p.x + step(rng), p.y + step(rng)};
        }
        set.runs.push_back(std::move(run));
        set.roles.push_back(ChainRole::Breakline);

        EdgeProperties props;
        const int count = how_many(rng);
        for (int i = 0; i < count; ++i) {
            props = props | EdgeProperties::bit(which_bit(rng));
        }
        set.properties.push_back(props);
    }

    return set;
}

// Throws std::runtime_error rather than returning a result type: a generator
// that emits an invalid Pslg is a defect in the generator, and swallowing it
// into an optional would let a suite quietly test nothing.
[[nodiscard]] inline Pslg to_pslg(const PolylineSet& set) {
    PslgBuilder builder{set.vertices};
    for (std::size_t c = 0; c < set.runs.size(); ++c) {
        builder.add_chain(set.runs[c], set.roles[c], set.properties[c]);
    }
    PslgBuildResult result = std::move(builder).build<pred::DefaultKernel>();
    if (!result.ok()) {
        throw std::runtime_error{"noding_generators: generated an invalid Pslg"};
    }
    return std::move(*result.pslg);
}

// The same geometry with the vertex buffer permuted and every chain index
// remapped. Used to assert that node ids are a function of the point SET: they
// are NodeSet's sorted order, so a permutation of the input must produce
// bit-identical output (05b-noder-driver.md, "Determinism").
[[nodiscard]] inline PolylineSet shuffled(const PolylineSet& set, std::uint64_t seed) {
    std::mt19937_64 rng{seed};

    std::vector<std::uint32_t> order(set.vertices.size());
    std::iota(order.begin(), order.end(), 0u);
    std::shuffle(order.begin(), order.end(), rng);

    // order[new] = old; invert it so chain indices can be rewritten.
    std::vector<std::uint32_t> where(set.vertices.size());
    for (std::uint32_t n = 0; n < order.size(); ++n) {
        where[order[n]] = n;
    }

    PolylineSet out;
    out.roles = set.roles;
    out.properties = set.properties;
    out.vertices.reserve(set.vertices.size());
    for (const std::uint32_t old : order) {
        out.vertices.push_back(set.vertices[old]);
    }
    for (const std::vector<std::uint32_t>& run : set.runs) {
        std::vector<std::uint32_t> remapped;
        remapped.reserve(run.size());
        for (const std::uint32_t i : run) {
            remapped.push_back(where[i]);
        }
        out.runs.push_back(std::move(remapped));
    }
    return out;
}

}  // namespace terrain::testing
