// Increment 21a, QW3 (docs/increments/21-parallel-refine.md, section 3, and
// "Pinned by the red suite (21a)"): the next round's `active` is rebuilt by a
// merge instead of sort + unique. 21a is bit-identical, so the only contract is
// that the vector is the one today's code builds. End to end that is held by
// 14's T6 and 18's golden digests, untouched; this suite holds it on
// constructed inputs, where the edge cases can be named.
//
// Interface pinned here (include/terrain/refinement/refine.hpp):
//
//   namespace terrain::refinement::detail {
//   void rebuild_active(std::span<const char> touched,
//                       std::span<const std::uint32_t> skipped,
//                       std::vector<std::uint32_t>& active);
//   }
//
// active's previous contents are replaced. Precondition: skipped ascending
// (refine fills it while walking the ascending `active`) and every entry below
// touched.size(). A skipped slot may also be touched; it appears once.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <terrain/refinement/refine.hpp>

#include <algorithm>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

using terrain::refinement::detail::rebuild_active;

namespace {

// Today's code, refine.hpp at 09b005b, lines 360-366, verbatim but for names.
std::vector<std::uint32_t> reference(const std::vector<char>& touched,
                                     const std::vector<std::uint32_t>& skipped) {
    std::vector<std::uint32_t> active;
    for (std::uint32_t t = 0; t < touched.size(); ++t)
        if (touched[t] != 0)
            active.push_back(t);
    active.insert(active.end(), skipped.begin(), skipped.end());
    std::sort(active.begin(), active.end());
    active.erase(std::unique(active.begin(), active.end()), active.end());
    return active;
}

std::vector<std::uint32_t> rebuilt(const std::vector<char>& touched,
                                   const std::vector<std::uint32_t>& skipped,
                                   std::vector<std::uint32_t> active = {}) {
    rebuild_active(std::span<const char>{touched}, std::span<const std::uint32_t>{skipped}, active);
    return active;
}

std::vector<char> touched_at(std::size_t size, std::initializer_list<std::uint32_t> at) {
    std::vector<char> touched(size, 0);
    for (const auto t : at) touched[t] = 1;
    return touched;
}

}  // namespace

TEST_CASE("rebuild_active: both empty gives an empty active", "[refinement][active]") {
    REQUIRE(rebuilt({}, {}).empty());
    REQUIRE(rebuilt(std::vector<char>(7, 0), {}).empty());
}

TEST_CASE("rebuild_active replaces the previous contents", "[refinement][active]") {
    const std::vector<std::uint32_t> stale{0, 1, 2, 3, 4, 5, 6, 99};
    REQUIRE(rebuilt({}, {}, stale).empty());
    const auto touched = touched_at(6, {4});
    REQUIRE(rebuilt(touched, {1}, stale) == std::vector<std::uint32_t>{1, 4});
}

TEST_CASE("rebuild_active: only touched, only skipped", "[refinement][active]") {
    REQUIRE(rebuilt(touched_at(5, {0, 2, 4}), {}) == std::vector<std::uint32_t>{0, 2, 4});
    REQUIRE(rebuilt(std::vector<char>(9, 0), {1, 3, 8}) == std::vector<std::uint32_t>{1, 3, 8});
}

TEST_CASE("rebuild_active: a skipped slot touched later appears once", "[refinement][active]") {
    // Duplicates between the two inputs, at the front, middle and back.
    const auto touched = touched_at(10, {0, 4, 5, 9});
    const std::vector<std::uint32_t> skipped{0, 3, 5, 9};
    REQUIRE(rebuilt(touched, skipped) == std::vector<std::uint32_t>{0, 3, 4, 5, 9});
    REQUIRE(rebuilt(touched, skipped) == reference(touched, skipped));
}

TEST_CASE("rebuild_active: interleaved inputs", "[refinement][active]") {
    const auto touched = touched_at(12, {1, 3, 5, 7, 11});
    const std::vector<std::uint32_t> skipped{0, 2, 4, 6, 8, 10};
    REQUIRE(rebuilt(touched, skipped) == std::vector<std::uint32_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11});
}

TEST_CASE("rebuild_active: every slot touched and every slot skipped", "[refinement][active]") {
    std::vector<std::uint32_t> all(33);
    for (std::uint32_t i = 0; i < all.size(); ++i) all[i] = i;
    REQUIRE(rebuilt(std::vector<char>(33, 1), all) == all);
}

TEST_CASE("rebuild_active: touched holds values other than 1", "[refinement][active]") {
    // refine writes 1, but resize(.., 1) and a char's truth are the contract:
    // any nonzero entry is touched.
    std::vector<char> touched{0, 2, 0, -1, 1};
    REQUIRE(rebuilt(touched, {2}) == std::vector<std::uint32_t>{1, 2, 3, 4});
}

TEST_CASE("rebuild_active equals today's sort + unique on random rounds", "[refinement][active]") {
    // Shaped like a round: `size` slots, a touched density, and `skipped` an
    // ascending subset drawn independently, so it overlaps touched at random.
    const std::uint32_t seed = GENERATE(range(1u, 41u));
    std::mt19937 rng{seed};
    const std::size_t size = rng() % 2000;
    const std::uint32_t touched_per_mille = rng() % 1001, skipped_per_mille = rng() % 1001;
    std::vector<char> touched(size, 0);
    std::vector<std::uint32_t> skipped;
    for (std::uint32_t t = 0; t < size; ++t) {
        touched[t] = rng() % 1000 < touched_per_mille ? 1 : 0;
        if (rng() % 1000 < skipped_per_mille) skipped.push_back(t);
    }
    CAPTURE(seed, size, touched_per_mille, skipped_per_mille);
    REQUIRE(rebuilt(touched, skipped) == reference(touched, skipped));
    // And into a reused, non-empty buffer, as refine reuses `active`.
    REQUIRE(rebuilt(touched, skipped, std::vector<std::uint32_t>(size / 2 + 3, 7)) ==
            reference(touched, skipped));
}
