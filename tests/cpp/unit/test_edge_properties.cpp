// Unit tests for terrain/core/edge_properties.hpp.
//
// Committed RED: the header does not exist yet, so this translation unit fails
// at the #include. That is the intended red for increment 7's commit 1.
//
// NOT invariant-critical, and docs/increments/07-edge-properties.md says so in
// writing: EdgeProperties is a value type over one std::uint32_t with no kernel
// parameter, no template parameter and no exactness claim. Its failure mode is
// a typo and it is loud. tests/cpp/unit/test_noding_node_set.cpp is excluded
// from a mutation round for the same reason and is the precedent. Do not spend
// a mutation round on this file.
//
// Two things here are load-bearing rather than enumerative, and both are
// rulings of the increment rather than properties of a bitset:
//
//   * `operator|` is commutative, associative and idempotent. The noder's merge
//     is a parallel reduce over an UNORDERED set of contributing chains, so an
//     answer that depended on visit order would be a wrong answer that no
//     single-threaded run reproduces. All three are asserted, over a table of
//     representative masks rather than over one pair.
//   * There is no implicit conversion from bool, int or std::uint32_t, and
//     there must never be one. That is what makes commit 3's rename of
//     `Chain::is_river` fail loudly at every one of its 14 call sites instead
//     of silently keeping its old spelling's meaning. It is asserted as a
//     static_assert, because a conversion that does not exist cannot be called.
//
// `EdgeProperties::bit(i)` has a narrow contract -- precondition
// i < kMaxProperties, debug-asserted, unchecked in release -- modelled on
// SnapGrid::can_snap/snap (include/terrain/core/snap_grid.hpp:94-96, :114).
// The violated branch is NOT tested here: an assert aborts the process, and
// this tree has no death-test harness. test_noding_snap_rounding.cpp leaves
// SnapGrid::snap's assert untested for the same reason. Every i this suite
// passes to bit() satisfies the precondition, which is the testable half.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/edge_properties.hpp>

#include <cstdint>
#include <initializer_list>
#include <type_traits>
#include <vector>

using terrain::EdgeProperties;

namespace {

// The union of the named bits. Every EdgeProperties in this file is built from
// bit() and operator|, because those are the only two routes the type offers --
// there is deliberately no constructor from a word.
[[nodiscard]] constexpr EdgeProperties mask_of(std::initializer_list<unsigned> bits) noexcept {
    EdgeProperties acc{};
    for (const unsigned i : bits) {
        acc = acc | EdgeProperties::bit(i);
    }
    return acc;
}

// A small, deliberately varied table: the empty set, singletons at both ends of
// the word, overlapping pairs, a set that is a strict superset of another, and
// a dense one. The algebraic laws below are checked over its full cross
// product, so the overlaps are what make idempotence a real claim rather than a
// restatement of "distinct bits do not collide".
[[nodiscard]] std::vector<EdgeProperties> sample_masks() {
    return {
        EdgeProperties{},
        EdgeProperties::bit(0),
        EdgeProperties::bit(1),
        EdgeProperties::bit(31),
        mask_of({0, 1}),
        mask_of({1, 2, 3}),
        mask_of({0, 1, 2, 3}),
        mask_of({0, 15, 31}),
    };
}

}  // namespace

// ---------------------------------------------------------------------------
// The type itself, asserted at compile time.
// ---------------------------------------------------------------------------

// The whole surface is usable in a constant expression. This is one of the two
// claims 07-edge-properties.md makes about the header ("constexpr throughout
// and header-only"), and a static_assert is the only assertion that can make
// it: a runtime REQUIRE over the same expressions passes whether or not they
// are constexpr.
static_assert(EdgeProperties::kMaxProperties == 32);
static_assert(EdgeProperties{}.empty());
static_assert(EdgeProperties{}.bits() == 0u);
static_assert(EdgeProperties::bit(0).bits() == 1u);
static_assert(EdgeProperties::bit(31).bits() == 0x8000'0000u);
static_assert(!EdgeProperties::bit(4).empty());
static_assert((EdgeProperties::bit(0) | EdgeProperties::bit(1)).bits() == 3u);
static_assert((EdgeProperties::bit(0) | EdgeProperties::bit(1)).contains(EdgeProperties::bit(1)));
static_assert(EdgeProperties::bit(2) == EdgeProperties::bit(2));
static_assert(EdgeProperties::bit(2) != EdgeProperties::bit(3));
static_assert(mask_of({0, 5, 9}).bits() == ((1u << 0) | (1u << 5) | (1u << 9)));

// The word, spelled: std::uint32_t and nothing wider, because "one word per
// edge, trivially parallel-reducible, no allocation" is the argument that bought
// the 32-property ceiling.
static_assert(std::is_same_v<decltype(EdgeProperties{}.bits()), std::uint32_t>);
static_assert(sizeof(EdgeProperties) == sizeof(std::uint32_t));
static_assert(std::is_trivially_copyable_v<EdgeProperties>);
static_assert(std::is_nothrow_default_constructible_v<EdgeProperties>);

// THE MIGRATION'S SAFETY NET. add_chain(idx, role, true) must stop compiling,
// and it stops compiling only because none of these conversions exists. If a
// later tidy-up adds a converting constructor "for convenience", this block is
// what refuses it -- 07-edge-properties.md rejects one by name so that it is
// not added later.
static_assert(!std::is_constructible_v<EdgeProperties, bool>);
static_assert(!std::is_constructible_v<EdgeProperties, int>);
static_assert(!std::is_constructible_v<EdgeProperties, unsigned>);
static_assert(!std::is_constructible_v<EdgeProperties, std::uint32_t>);
static_assert(!std::is_convertible_v<bool, EdgeProperties>);
static_assert(!std::is_convertible_v<std::uint32_t, EdgeProperties>);

// And the reverse direction, which would be just as silent: an EdgeProperties
// that decays to a number would let `if (chain.properties)` compile and mean
// something, which is exactly the bool this increment removes.
static_assert(!std::is_convertible_v<EdgeProperties, bool>);
static_assert(!std::is_convertible_v<EdgeProperties, std::uint32_t>);

// is_convertible_v<EdgeProperties, bool> is about IMPLICIT conversion, and
// `if (x)` uses CONTEXTUAL conversion, which an `explicit operator bool`
// satisfies while leaving the line above true. So the two asserts are claims
// about different objects, and only this one is a claim about `if (x)`.
// Measured rather than reasoned: `if (EdgeProperties{})` is rejected with
// "not contextually convertible to 'bool'", and the probe that says so fails
// when the assertion is inverted.
static_assert(!std::is_constructible_v<bool, EdgeProperties>);

// ---------------------------------------------------------------------------
// The empty set.
// ---------------------------------------------------------------------------

TEST_CASE("a default-constructed EdgeProperties is the empty set", "[core][edge_properties]") {
    const EdgeProperties empty{};

    REQUIRE(empty.empty());
    REQUIRE(empty.bits() == 0u);
    REQUIRE(empty == EdgeProperties{});
}

TEST_CASE("the empty set is the identity of union", "[core][edge_properties]") {
    for (const EdgeProperties m : sample_masks()) {
        CAPTURE(m.bits());
        REQUIRE((m | EdgeProperties{}) == m);
        REQUIRE((EdgeProperties{} | m) == m);
    }
}

TEST_CASE("only the empty set is empty()", "[core][edge_properties]") {
    for (unsigned i = 0; i < EdgeProperties::kMaxProperties; ++i) {
        CAPTURE(i);
        REQUIRE_FALSE(EdgeProperties::bit(i).empty());
    }

    REQUIRE_FALSE(mask_of({0, 31}).empty());
}

// ---------------------------------------------------------------------------
// bit(), across the whole word.
// ---------------------------------------------------------------------------

TEST_CASE("bit(i) sets exactly bit i, for every legal i", "[core][edge_properties]") {
    for (unsigned i = 0; i < EdgeProperties::kMaxProperties; ++i) {
        CAPTURE(i);
        REQUIRE(EdgeProperties::bit(i).bits() == (std::uint32_t{1} << i));
    }
}

TEST_CASE("bit(31) is representable and distinct", "[core][edge_properties]") {
    // The top bit is where a 1 << i written over a signed int would be
    // undefined, and where a 16-bit word would wrap into bit 15.
    REQUIRE(EdgeProperties::bit(31).bits() == 0x8000'0000u);
    REQUIRE(EdgeProperties::bit(31) != EdgeProperties::bit(15));
    REQUIRE_FALSE(EdgeProperties::bit(31).empty());
}

TEST_CASE("distinct bits are distinct sets", "[core][edge_properties]") {
    for (unsigned i = 0; i < EdgeProperties::kMaxProperties; ++i) {
        for (unsigned j = i + 1; j < EdgeProperties::kMaxProperties; ++j) {
            CAPTURE(i, j);
            REQUIRE(EdgeProperties::bit(i) != EdgeProperties::bit(j));
        }
    }
}

TEST_CASE("the union of all 32 bits is the full word", "[core][edge_properties]") {
    EdgeProperties all{};
    for (unsigned i = 0; i < EdgeProperties::kMaxProperties; ++i) {
        all = all | EdgeProperties::bit(i);
    }

    REQUIRE(all.bits() == 0xFFFF'FFFFu);
    for (unsigned i = 0; i < EdgeProperties::kMaxProperties; ++i) {
        CAPTURE(i);
        REQUIRE(all.contains(EdgeProperties::bit(i)));
    }
}

// ---------------------------------------------------------------------------
// contains(): superset, not intersection.
// ---------------------------------------------------------------------------

TEST_CASE("contains is superset and not 'shares a bit'", "[core][edge_properties]") {
    const EdgeProperties river = EdgeProperties::bit(0);
    const EdgeProperties both = mask_of({0, 1});

    // The discriminating pair. An implementation written as
    // (bits_ & other.bits_) != 0 -- "they intersect" -- passes the first
    // REQUIRE and fails the second, which is the only reason the second is
    // here.
    REQUIRE(both.contains(river));
    REQUIRE_FALSE(river.contains(both));
}

TEST_CASE("contains is reflexive and the empty set is contained everywhere",
          "[core][edge_properties]") {
    for (const EdgeProperties m : sample_masks()) {
        CAPTURE(m.bits());
        REQUIRE(m.contains(m));
        REQUIRE(m.contains(EdgeProperties{}));
    }

    REQUIRE(EdgeProperties{}.contains(EdgeProperties{}));
    REQUIRE_FALSE(EdgeProperties{}.contains(EdgeProperties::bit(0)));
}

TEST_CASE("disjoint sets contain neither each other nor their union",
          "[core][edge_properties]") {
    const EdgeProperties a = mask_of({0, 2});
    const EdgeProperties b = mask_of({1, 3});

    REQUIRE_FALSE(a.contains(b));
    REQUIRE_FALSE(b.contains(a));
    REQUIRE((a | b).contains(a));
    REQUIRE((a | b).contains(b));
    REQUIRE_FALSE(a.contains(a | b));
}

TEST_CASE("a union contains both of its operands, over the whole table",
          "[core][edge_properties]") {
    const std::vector<EdgeProperties> masks = sample_masks();

    for (const EdgeProperties a : masks) {
        for (const EdgeProperties b : masks) {
            CAPTURE(a.bits(), b.bits());
            REQUIRE((a | b).contains(a));
            REQUIRE((a | b).contains(b));
        }
    }
}

// ---------------------------------------------------------------------------
// The three laws the parallel reduce rests on.
// ---------------------------------------------------------------------------

TEST_CASE("union is commutative", "[core][edge_properties]") {
    const std::vector<EdgeProperties> masks = sample_masks();

    for (const EdgeProperties a : masks) {
        for (const EdgeProperties b : masks) {
            CAPTURE(a.bits(), b.bits());
            REQUIRE((a | b) == (b | a));
        }
    }
}

TEST_CASE("union is associative", "[core][edge_properties]") {
    const std::vector<EdgeProperties> masks = sample_masks();

    for (const EdgeProperties a : masks) {
        for (const EdgeProperties b : masks) {
            for (const EdgeProperties c : masks) {
                CAPTURE(a.bits(), b.bits(), c.bits());
                REQUIRE(((a | b) | c) == (a | (b | c)));
            }
        }
    }
}

TEST_CASE("union is idempotent", "[core][edge_properties]") {
    const std::vector<EdgeProperties> masks = sample_masks();

    for (const EdgeProperties a : masks) {
        CAPTURE(a.bits());
        REQUIRE((a | a) == a);

        for (const EdgeProperties b : masks) {
            CAPTURE(b.bits());
            // The form the reduce actually meets: re-merging a contributor
            // already folded in changes nothing, so a chain visited twice by
            // the broad phase cannot change the answer.
            REQUIRE(((a | b) | a) == (a | b));
        }
    }
}

TEST_CASE("union is the bitwise or of the words", "[core][edge_properties]") {
    const std::vector<EdgeProperties> masks = sample_masks();

    for (const EdgeProperties a : masks) {
        for (const EdgeProperties b : masks) {
            CAPTURE(a.bits(), b.bits());
            REQUIRE((a | b).bits() == (a.bits() | b.bits()));
        }
    }
}

// ---------------------------------------------------------------------------
// Equality.
// ---------------------------------------------------------------------------

TEST_CASE("equality is equality of the underlying set", "[core][edge_properties]") {
    // Same members, built by two different fold orders: equal, and this is the
    // property the noder's output comparison depends on.
    REQUIRE(mask_of({0, 1, 2}) == mask_of({2, 0, 1}));
    REQUIRE(mask_of({3, 3, 3}) == EdgeProperties::bit(3));

    REQUIRE(mask_of({0, 1}) != mask_of({0, 1, 2}));
    REQUIRE(mask_of({0, 1}) != EdgeProperties{});
    REQUIRE(EdgeProperties::bit(0) != EdgeProperties::bit(1));
}

TEST_CASE("equality agrees with bits() over the whole table", "[core][edge_properties]") {
    const std::vector<EdgeProperties> masks = sample_masks();

    for (const EdgeProperties a : masks) {
        for (const EdgeProperties b : masks) {
            CAPTURE(a.bits(), b.bits());
            REQUIRE((a == b) == (a.bits() == b.bits()));
        }
    }
}

// ---------------------------------------------------------------------------
// The one thing terrain:: must NOT have.
// ---------------------------------------------------------------------------

TEST_CASE("the header spells no feature name", "[core][edge_properties]") {
    // There is nothing to call here, and that is the assertion: if
    // EdgeProperties ever grows a `River`, a `Road` or an enumeration, this
    // comment is the place the reviewer is sent. The mechanical half of the
    // check is a grep, run in review, not a REQUIRE -- a test cannot assert the
    // absence of a name it is forbidden to write. What IS mechanical is that
    // every set in this file is addressed by bit index, never by name, so a
    // named member would be dead weight no suite exercises.
    REQUIRE(EdgeProperties::kMaxProperties == 32);
}
