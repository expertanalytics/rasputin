// Property tests for terrain/noding/broad_phase.hpp -- increment 5b's index.
//
// INVARIANT-CRITICAL, mutation round (docs/increments/05b-noder-driver.md,
// "What is worth testing"). This is 05-noder.md risk 12's whole mitigation and
// the one function in the increment with an oracle that does not touch the
// noder: the contract mentions no segment pair, no node and no snap grid, which
// is exactly what makes an independent oracle possible. A check written in terms
// of "the pairs the noder needed" would share the defect it is looking for.
//
// THE ORACLE IS ALL-PAIRS BBOX INTERSECTION COMPUTED WITHOUT THE INDEX, and the
// assertion is one-directional:
//
//     every segment whose closed bbox intersects q is visited.
//
// False positives are free and none is asserted against. That asymmetry is the
// contract, not a weakening of the test: any positive bucket count is correct
// and only sqrt(n) is fast, so no fixture can distinguish sqrt(n) from 1 except
// by timing, and a timing assertion is not a test (05b-noder-driver.md, "Not
// mutants").
//
// No rapidcheck. Catch2's GENERATE over a fixed seed range plus a seeded
// std::mt19937_64, the same arrangement as prop_noding_snap_invariants.cpp: a
// failure prints its seed as the generator index and rerunning that section
// reproduces it exactly.
//
// THE GENERATOR HERE IS LOCAL AND IS NOT noding_generators.h's. That header
// produces polyline sets with a controllable intersection density, which is what
// the noder's own property suite needs; this suite needs long segments and
// adversarial query boxes, and sharing a generator between the index and its
// consumer is how a blind spot in one comes to be a blind spot in both.
//
// FastKernel is not banned here because no kernel appears: the index compares
// bounding boxes and never asks an orientation. That is also why this target
// registers with the plain add_terrain_test helper.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/noding/broad_phase.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

using terrain::Box2;
using terrain::Point2;
using terrain::Segment2;
using terrain::bbox;
using terrain::noding::BroadPhase;

namespace {

constexpr int seed_count = 24;

// UTM33 easting/northing: the magnitudes the engine actually runs at, so that a
// bucket arithmetic defect that only shows up when the domain origin is far from
// zero shows up here.
constexpr double origin_x = 430000.0;
constexpr double origin_y = 6900000.0;

[[nodiscard]] std::vector<std::uint32_t> visit(const BroadPhase& index, const Box2& q) {
    std::vector<std::uint32_t> seen;
    index.for_each_candidate(q, [&seen](std::uint32_t i) { seen.push_back(i); });
    std::ranges::sort(seen);
    const auto dup = std::ranges::unique(seen);
    seen.erase(dup.begin(), dup.end());
    return seen;
}

// The oracle. No index, no buckets, no shared code with the thing under test.
[[nodiscard]] std::vector<std::uint32_t> oracle(std::span<const Segment2> segments,
                                                const Box2& q) {
    std::vector<std::uint32_t> required;
    for (std::size_t i = 0; i < segments.size(); ++i) {
        if (bbox(segments[i]).intersects(q)) {
            required.push_back(static_cast<std::uint32_t>(i));
        }
    }
    return required;
}

void require_conservative(std::span<const Segment2> segments, const Box2& q) {
    const BroadPhase index{segments};
    const std::vector<std::uint32_t> seen = visit(index, q);

    for (const std::uint32_t i : seen) {
        REQUIRE(i < segments.size());
    }
    for (const std::uint32_t i : oracle(segments, q)) {
        // A false negative here is the one defect the function's entire job is
        // not to have: the noder would then never compare a pair that crosses,
        // and the verification pass would agree with it if it shared the index.
        REQUIRE(std::ranges::binary_search(seen, i));
    }
}

// Segments in a 1 km square, with one in eight deliberately long -- up to the
// full domain -- so that the interior of a segment routinely leaves the bucket
// its endpoints are in. That ratio is what makes mutant 1 (bucketing by
// endpoints) reachable at all from the random arm.
[[nodiscard]] std::vector<Segment2> random_segments(std::uint64_t seed, std::size_t n) {
    std::mt19937_64 rng{seed};
    std::uniform_real_distribution<double> coord{0.0, 1000.0};
    std::uniform_real_distribution<double> shortish{-20.0, 20.0};
    std::uniform_int_distribution<int> long_one{0, 7};

    std::vector<Segment2> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const Point2 a{origin_x + coord(rng), origin_y + coord(rng)};
        const Point2 b = long_one(rng) == 0
                             ? Point2{origin_x + coord(rng), origin_y + coord(rng)}
                             : Point2{a.x + shortish(rng), a.y + shortish(rng)};
        out.push_back(Segment2{a, b});
    }
    return out;
}

[[nodiscard]] Box2 random_box(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> coord{0.0, 1000.0};
    std::uniform_real_distribution<double> extent{0.0, 50.0};
    const double x = origin_x + coord(rng);
    const double y = origin_y + coord(rng);
    return Box2{Point2{x, y}, Point2{x + extent(rng), y + extent(rng)}};
}

}  // namespace

// ---------------------------------------------------------------------------
// The shape of the query. Compile-time, and mutant 2's first half.
// ---------------------------------------------------------------------------

namespace {

template <class I>
concept HasPointQuery = requires(const I& i, Point2 p) {
    i.for_each_candidate(p, [](std::uint32_t) {});
};

template <class I>
concept HasConstBoxQuery = requires(const I& i, Box2 q) {
    i.for_each_candidate(q, [](std::uint32_t) {});
};

}  // namespace

TEST_CASE("the index is queried by box only, and only through a const object",
          "[noding][broad_phase]") {
    // MUTANT 2, first half. A Point2 overload is the one live way to break
    // conservativeness, because the caller that wants it is the node query in
    // node.hpp step 4 and the right argument there is the node's cell box, never
    // world(g). The overload must not exist; if it appears in a diff, that is
    // the mutant as a reviewable change, and this assertion is the tripwire.
    STATIC_REQUIRE(!HasPointQuery<BroadPhase>);

    // Immutable after construction and const-queried: what lets the split pass
    // and the verification pass share one instance, and what lets a future
    // parallel-for query it without synchronisation.
    STATIC_REQUIRE(HasConstBoxQuery<BroadPhase>);
}

// ---------------------------------------------------------------------------
// The property.
// ---------------------------------------------------------------------------

TEST_CASE("every segment whose bbox meets the query box is visited", "[noding][broad_phase]") {
    const int seed = GENERATE(range(0, seed_count));
    const std::vector<Segment2> segments = random_segments(static_cast<std::uint64_t>(seed), 60);

    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 0x9e3779b97f4a7c15ULL};
    for (int q = 0; q < 16; ++q) {
        require_conservative(segments, random_box(rng));
    }
}

TEST_CASE("each segment's own bbox retrieves at least itself", "[noding][broad_phase]") {
    // The self-query is the one the split pass actually issues in step 1, and
    // it is the query a by-endpoint bucketing gets right by accident. It is here
    // to stop a regression, not to catch mutant 1.
    const int seed = GENERATE(range(0, seed_count));
    const std::vector<Segment2> segments = random_segments(static_cast<std::uint64_t>(seed), 40);
    const BroadPhase index{segments};

    for (std::size_t i = 0; i < segments.size(); ++i) {
        const std::vector<std::uint32_t> seen = visit(index, bbox(segments[i]));
        REQUIRE(std::ranges::binary_search(seen, static_cast<std::uint32_t>(i)));
    }
}

TEST_CASE("a long segment is found by a box its endpoints miss", "[noding][broad_phase]") {
    // MUTANT 1, named and hand-built rather than left to the random arm.
    //
    // Sixty-four short segments establish a bucket grid of about eight per axis
    // over the domain; the sixty-fifth runs corner to corner. A query box in the
    // middle of the domain touches the diagonal's interior and neither of its
    // endpoints, and shares no bucket with either endpoint. An index that
    // inserted a segment into the buckets of its two endpoints reports nothing
    // and the noder never compares the pair.
    std::vector<Segment2> segments;
    std::mt19937_64 rng{7};
    std::uniform_real_distribution<double> coord{0.0, 1000.0};
    for (int i = 0; i < 64; ++i) {
        const Point2 a{origin_x + coord(rng), origin_y + coord(rng)};
        segments.push_back(Segment2{a, Point2{a.x + 1.0, a.y + 1.0}});
    }
    const auto diagonal = static_cast<std::uint32_t>(segments.size());
    segments.push_back(Segment2{Point2{origin_x, origin_y}, Point2{origin_x + 1000.0, origin_y + 1000.0}});

    const Box2 q{Point2{origin_x + 499.0, origin_y + 499.0},
                 Point2{origin_x + 501.0, origin_y + 501.0}};
    REQUIRE(bbox(segments[diagonal]).intersects(q));
    REQUIRE_FALSE(q.contains(segments[diagonal].a));
    REQUIRE_FALSE(q.contains(segments[diagonal].b));

    const BroadPhase index{segments};
    const std::vector<std::uint32_t> seen = visit(index, q);
    REQUIRE(std::ranges::binary_search(seen, diagonal));
}

TEST_CASE("a query box straddling a bucket boundary finds both sides",
          "[noding][broad_phase]") {
    // MUTANT 2, second half, in the form the caller can actually produce: a node
    // whose cell straddles a bucket boundary. Queried by the cell's box every
    // segment on either side is required; queried by the cell's centre the
    // segments on the far side are lost, and a T-junction with them.
    //
    // Four segments, one per quadrant around the domain centre, and a query box
    // centred on that point.
    const double cx = origin_x + 500.0;
    const double cy = origin_y + 500.0;
    std::vector<Segment2> segments;
    std::mt19937_64 rng{11};
    std::uniform_real_distribution<double> coord{0.0, 1000.0};
    for (int i = 0; i < 100; ++i) {
        const Point2 a{origin_x + coord(rng), origin_y + coord(rng)};
        segments.push_back(Segment2{a, Point2{a.x + 0.5, a.y - 0.5}});
    }
    for (const double dx : {-0.02, 0.02}) {
        for (const double dy : {-0.02, 0.02}) {
            segments.push_back(Segment2{Point2{cx + dx, cy + dy},
                                        Point2{cx + dx * 1.5, cy + dy * 1.5}});
        }
    }

    // Half a decimetre cell about the centre, as cell_min/cell_max would give.
    const Box2 cell{Point2{cx - 0.05, cy - 0.05}, Point2{cx + 0.05, cy + 0.05}};
    require_conservative(segments, cell);

    // And it really is the interesting case: all four quadrant segments qualify.
    REQUIRE(oracle(segments, cell).size() >= 4);
}

// ---------------------------------------------------------------------------
// Degenerate extents. The design's claim is that these need no branch in the
// query, which is a claim about the answers and is checked as one.
// ---------------------------------------------------------------------------

TEST_CASE("an empty segment set answers every query with nothing", "[noding][broad_phase]") {
    const std::vector<Segment2> none;
    const BroadPhase index{none};
    REQUIRE(visit(index, Box2{Point2{0.0, 0.0}, Point2{1.0, 1.0}}).empty());
    REQUIRE(visit(index, Box2{Point2{origin_x, origin_y}, Point2{origin_x, origin_y}}).empty());
}

TEST_CASE("a single segment is retrieved by anything that meets its bbox",
          "[noding][broad_phase]") {
    const std::vector<Segment2> one{Segment2{Point2{origin_x, origin_y},
                                             Point2{origin_x + 3.0, origin_y + 4.0}}};
    require_conservative(one, Box2{Point2{origin_x + 1.0, origin_y + 1.0},
                                   Point2{origin_x + 2.0, origin_y + 2.0}});
    require_conservative(one, bbox(one[0]));

    // NOTHING HERE ASSERTS THAT A FAR-AWAY QUERY BOX RETURNS NOTHING, and an
    // earlier revision of this case did. It was wrong against the contract, not
    // merely strict: false positives are free, a one-bucket index is a correct
    // index, and a suite that forbids them is a suite that fails on a legal
    // implementation. It failed on the first one it met. Recorded rather than
    // silently deleted, because the assertion reads perfectly reasonable.
    require_conservative(one, Box2{Point2{0.0, 0.0}, Point2{1.0, 1.0}});
}

TEST_CASE("a set with zero extent on one axis is one bucket on that axis",
          "[noding][broad_phase]") {
    // Vertical collinear segments: the set's bounding box has zero width, and
    // ceil(sqrt(n)) buckets over a zero range is a division by zero waiting to
    // happen. The design says this is handled without a branch in the query; the
    // observable consequence is that the answers are still right.
    std::vector<Segment2> vertical;
    for (int i = 0; i < 20; ++i) {
        const double y = origin_y + static_cast<double>(i);
        vertical.push_back(Segment2{Point2{origin_x, y}, Point2{origin_x, y + 0.5}});
    }
    require_conservative(vertical, Box2{Point2{origin_x - 1.0, origin_y + 4.0},
                                        Point2{origin_x + 1.0, origin_y + 6.0}});
    require_conservative(vertical, bbox(vertical[19]));
}

TEST_CASE("a set of identical segments has zero extent on both axes",
          "[noding][broad_phase]") {
    const Segment2 s{Point2{origin_x, origin_y}, Point2{origin_x + 1.0, origin_y}};
    const std::vector<Segment2> same(16, s);
    require_conservative(same, bbox(s));

    const BroadPhase index{same};
    REQUIRE(visit(index, bbox(s)).size() == 16);
}

TEST_CASE("zero-length segments are indexed like any other", "[noding][broad_phase]") {
    // A degenerate segment gets no special case anywhere else in this module
    // (core/segment.hpp's ruling) and gets none here either. Its bbox is a
    // point and a query box containing that point must retrieve it.
    std::vector<Segment2> points;
    for (int i = 0; i < 12; ++i) {
        const Point2 p{origin_x + static_cast<double>(i), origin_y};
        points.push_back(Segment2{p, p});
    }
    require_conservative(points, Box2{Point2{origin_x + 2.5, origin_y - 1.0},
                                      Point2{origin_x + 5.5, origin_y + 1.0}});

    const BroadPhase index{points};
    const Point2 p{origin_x + 3.0, origin_y};
    REQUIRE(std::ranges::binary_search(visit(index, Box2{p, p}), 3u));
}

TEST_CASE("a degenerate query box is a legal query", "[noding][broad_phase]") {
    const int seed = GENERATE(range(0, seed_count));
    const std::vector<Segment2> segments = random_segments(static_cast<std::uint64_t>(seed), 40);

    // Box2{p, p} is a point, and it is what a node query degenerates to at a
    // subnormal spacing. Nothing may divide by its extent.
    for (const Segment2& s : segments) {
        require_conservative(segments, Box2{s.a, s.a});
    }
}

TEST_CASE("sub-millimetre features inside a hundred-kilometre domain", "[noding][broad_phase]") {
    // The extreme-scale case the skill file requires of every geometry suite.
    // Bucket edges computed by subtracting two coordinates of magnitude 1e6 and
    // dividing lose every bit of a 1e-4 feature; the conservativeness contract
    // still has to hold, which is why the oracle compares boxes and not indices.
    std::vector<Segment2> segments;
    std::mt19937_64 rng{23};
    std::uniform_real_distribution<double> far{0.0, 100000.0};
    for (int i = 0; i < 80; ++i) {
        const Point2 a{origin_x + far(rng), origin_y + far(rng)};
        segments.push_back(Segment2{a, Point2{a.x + 0.0001, a.y + 0.0001}});
    }
    segments.push_back(Segment2{Point2{origin_x, origin_y},
                                Point2{origin_x + 100000.0, origin_y + 100000.0}});

    for (const Segment2& s : segments) {
        require_conservative(segments, bbox(s));
    }
}

// ---------------------------------------------------------------------------
// Determinism. node<K> is a pure function of (pslg, options) and the index is
// the only stateful thing inside it.
// ---------------------------------------------------------------------------

TEST_CASE("two indexes over the same segments answer identically", "[noding][broad_phase]") {
    const int seed = GENERATE(range(0, seed_count));
    const std::vector<Segment2> segments = random_segments(static_cast<std::uint64_t>(seed), 50);

    const BroadPhase a{segments};
    const BroadPhase b{segments};

    std::mt19937_64 rng{static_cast<std::uint64_t>(seed)};
    for (int q = 0; q < 8; ++q) {
        const Box2 box = random_box(rng);
        REQUIRE(visit(a, box) == visit(b, box));
    }
}
