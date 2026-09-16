// Unit tests for terrain::pred::DefaultKernel -- FilteredKernel<DetriaExact>,
// the kernel the triangulation will actually call.
//
// The filtered-kernel suite already proves the filter's logic against RefExact.
// What is new here, and what this file exists for, is that the backend is now
// real code with a real precondition enforced by a real signal. Two claims
// follow from that:
//
//   * DefaultKernel returns what FilteredKernel<RefExact> returns. Swapping the
//     exact backend must not move a single answer, since a backend is by
//     definition not a place where policy lives.
//
//   * DefaultKernel::incircle survives non-counterclockwise input. This is the
//     blocking condition for the increment, and it is tested twice, in two
//     different failure modes -- see the two sections below.

#include <catch2/catch_test_macros.hpp>

#include <exact_reference.hpp>
#include <point_families.hpp>

#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/detria_exact.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <format>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#define TERRAIN_HAS_FORK 1
#include <sys/wait.h>

#include <csignal>
#include <cstdlib>
#include <unistd.h>
#else
#define TERRAIN_HAS_FORK 0
#endif

using terrain::Point2;
using terrain::pred::DefaultKernel;
using terrain::pred::DetriaExact;
using terrain::pred::FilteredKernel;
using terrain::pred::GeometryKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::pred::reversed;
using terrain::test::CountingExact;
using terrain::test::RefExact;
using terrain::test::utm33_offset;

namespace {

using Reference = FilteredKernel<RefExact>;
using Counted = FilteredKernel<CountingExact<DetriaExact>>;

std::string describe(const Point2& a, const Point2& b, const Point2& c) {
    return std::format("a={} b={} c={}", a, b, c);
}

std::string describe(const Point2& a, const Point2& b, const Point2& c, const Point2& d) {
    return std::format("a={} b={} c={} d={}", a, b, c, d);
}

// Triples spanning the families the kernel must survive: general position,
// coincident vertices, exact collinearity at unit and UTM33 magnitude, slivers,
// and a long narrow corridor.
const std::vector<std::array<Point2, 3>>& orientation_corpus() {
    static const std::vector<std::array<Point2, 3>> corpus = {
        {Point2{0, 0}, Point2{1, 0}, Point2{0, 1}},
        {Point2{0, 0}, Point2{0, 1}, Point2{1, 0}},
        {Point2{-3, 2}, Point2{4, 2}, Point2{1, -5}},
        {Point2{0, 0}, Point2{1, 1}, Point2{2, 2}},
        {Point2{5, 5}, Point2{5, 5}, Point2{9, 1}},
        {Point2{5, 5}, Point2{5, 5}, Point2{5, 5}},
        {Point2{0, 0}, Point2{16777216, 1}, Point2{33554432, 1}},
        {utm33_offset({0, 0}), utm33_offset({1, 0}), utm33_offset({0, 1})},
        {utm33_offset({0, 0}), utm33_offset({4, 4}), utm33_offset({9, 9})},
        {utm33_offset({0, 0}), utm33_offset({1000000, 0}), utm33_offset({500000, 1})},
    };
    return corpus;
}

}  // namespace

// ---------------------------------------------------------------------------
// Identity and concept conformance.
// ---------------------------------------------------------------------------

// `DefaultKernel` is an alias, not a new type. Pinning that is what stops a
// later "temporary" redefinition -- to FastKernel while debugging, say -- from
// surviving review unnoticed: FastKernel also models GeometryKernel, so nothing
// else in the suite would object.
TEST_CASE("DefaultKernel is FilteredKernel over DetriaExact", "[predicates][default_kernel]") {
    STATIC_REQUIRE(std::is_same_v<DefaultKernel, FilteredKernel<DetriaExact>>);
    STATIC_REQUIRE(GeometryKernel<DefaultKernel>);
    STATIC_REQUIRE(std::is_empty_v<DefaultKernel>);
}

// ---------------------------------------------------------------------------
// The backend is interchangeable: same answers, different arithmetic.
// ---------------------------------------------------------------------------

TEST_CASE("DefaultKernel::orient2d agrees with FilteredKernel<RefExact>", "[predicates][default_kernel]") {
    for (const auto& t : orientation_corpus()) {
        INFO(describe(t[0], t[1], t[2]));
        REQUIRE(DefaultKernel::orient2d(t[0], t[1], t[2]) == Reference::orient2d(t[0], t[1], t[2]));
    }
}

TEST_CASE("DefaultKernel::orient2d agrees with the reference near degeneracy", "[predicates][default_kernel][adversarial]") {
    std::mt19937_64 rng{616161};
    for (int i = 0; i < 500; ++i) {
        const auto t = terrain::test::near_degenerate_triple(rng);
        INFO(describe(t.a, t.b, t.c));
        REQUIRE(DefaultKernel::orient2d(t.a, t.b, t.c) == t.expected);
        REQUIRE(DefaultKernel::orient2d(t.a, t.b, t.c) == Reference::orient2d(t.a, t.b, t.c));
    }
}

TEST_CASE("DefaultKernel::orient2d keeps its symmetries on degenerate input", "[predicates][default_kernel]") {
    for (const auto& t : orientation_corpus()) {
        INFO(describe(t[0], t[1], t[2]));
        const Orientation o = DefaultKernel::orient2d(t[0], t[1], t[2]);

        REQUIRE(DefaultKernel::orient2d(t[0], t[2], t[1]) == reversed(o));
        REQUIRE(DefaultKernel::orient2d(t[1], t[0], t[2]) == reversed(o));
        REQUIRE(DefaultKernel::orient2d(t[2], t[1], t[0]) == reversed(o));

        REQUIRE(DefaultKernel::orient2d(t[1], t[2], t[0]) == o);
        REQUIRE(DefaultKernel::orient2d(t[2], t[0], t[1]) == o);
    }
}

TEST_CASE("DefaultKernel::incircle agrees with FilteredKernel<RefExact>", "[predicates][default_kernel][incircle]") {
    const auto circle = terrain::test::integer_circle_radius_5();

    for (std::size_t i = 0; i < circle.size(); ++i) {
        // Deliberately *not* normalized: the triple is fed in whatever order
        // the circle walk produces, because `incircle` is specified as total.
        const Point2 a = circle[i];
        const Point2 b = circle[(i + 3) % circle.size()];
        const Point2 c = circle[(i + 6) % circle.size()];

        for (const Point2& d : {circle[(i + 9) % circle.size()], Point2{0, 0}, Point2{400, 400},
                                Point2{1, 1}, a}) {
            INFO(describe(a, b, c, d));
            REQUIRE(DefaultKernel::incircle(a, b, c, d) == Reference::incircle(a, b, c, d));
        }
    }
}

TEST_CASE("DefaultKernel::incircle is exact on a representable circle", "[predicates][default_kernel][adversarial]") {
    const auto circle = terrain::test::integer_circle_radius_5();
    for (std::size_t i = 0; i < circle.size(); ++i) {
        const Point2 a = utm33_offset(circle[i]);
        const Point2 b = utm33_offset(circle[(i + 3) % circle.size()]);
        const Point2 c = utm33_offset(circle[(i + 6) % circle.size()]);

        INFO(describe(a, b, c, utm33_offset(circle[(i + 9) % circle.size()])));
        REQUIRE(DefaultKernel::incircle(a, b, c, utm33_offset(circle[(i + 9) % circle.size()])) ==
                Incircle::Cocircular);
        REQUIRE(DefaultKernel::incircle(a, b, c, utm33_offset({0, 0})) == Incircle::Inside);
        REQUIRE(DefaultKernel::incircle(a, b, c, utm33_offset({40, 40})) == Incircle::Outside);
    }
}

TEST_CASE("DefaultKernel::incircle is insensitive to the order of the first three points", "[predicates][default_kernel][incircle]") {
    const auto circle = terrain::test::integer_circle_radius_5();
    const Point2 a = circle[0];
    const Point2 b = circle[4];
    const Point2 c = circle[8];

    const std::array<Point2, 3> permutations[] = {
        {a, b, c}, {a, c, b}, {b, a, c}, {b, c, a}, {c, a, b}, {c, b, a},
    };

    for (const Point2& d : {Point2{0, 0}, Point2{400, 400}, circle[2], Point2{1, 1}}) {
        const Incircle expected = Reference::incircle(a, b, c, d);
        for (const auto& p : permutations) {
            INFO(describe(p[0], p[1], p[2], d));
            REQUIRE(DefaultKernel::incircle(p[0], p[1], p[2], d) == expected);
        }
    }
}

// ---------------------------------------------------------------------------
// The backend is actually reached.
// ---------------------------------------------------------------------------

// Without this, most of the file is testing the *filter*, not the backend
// behind it. On well-separated input `FilteredKernel` answers from double
// arithmetic alone and never calls `DetriaExact`, so a backend whose enum
// mapping inverted Inside and Outside would leave every agreement test above
// green. The two things that have to be asserted together are that the exact
// path is entered, and that what comes back out of it is right.
TEST_CASE("DefaultKernel consults DetriaExact on input the filter cannot settle", "[predicates][default_kernel][filter]") {
    CountingExact<DetriaExact>::reset();

    REQUIRE(Counted::orient2d(Point2{0, 0}, Point2{1, 1}, Point2{2, 2}) == Orientation::Collinear);
    REQUIRE(CountingExact<DetriaExact>::orient2d_calls > 0);

    // Reset before the incircle block so its counter is read against a known
    // zero rather than against whatever the orient2d block left behind.
    CountingExact<DetriaExact>::reset();

    REQUIRE(Counted::incircle(Point2{5, 0}, Point2{0, 5}, Point2{-5, 0}, Point2{3, -4}) ==
            Incircle::Cocircular);
    REQUIRE(CountingExact<DetriaExact>::incircle_calls > 0);
}

// Every Inside/Outside answer the backend gives, observed through the kernel.
//
// This is harder to arrange than it looks, and the difficulty is worth stating
// because the obvious constructions do not work. On integer coordinates the
// incircle determinant is either exactly zero or large compared to the filter's
// relative bound of about 1e-15, so integer input reaches the exact path only
// when it is *exactly* cocircular -- where the answer is Cocircular and an
// inverted Inside/Outside mapping in the backend is invisible. Scaling an
// integer circle up and nudging it by a unit does not help: the nudge is
// relatively enormous and the filter settles it in double arithmetic.
//
// So the fourth point is placed one ulp off the circle instead. The circle
// through (1,0), (0,1), (-1,0) is the unit circle exactly, and every coordinate
// here is exactly representable, so `nextafter(-1, -2)` is strictly farther
// from the centre than the radius and `nextafter(-1, 0)` strictly nearer --
// known by construction, with no oracle needed, which is just as well since
// RefExact only referees integers. The separation is about 1.1e-16 relative:
// below the filter's bound, hence decided by DetriaExact.
TEST_CASE("DefaultKernel reports every circle location the backend decides", "[predicates][default_kernel][filter][adversarial]") {
    const Point2 a{1, 0};
    const Point2 b{0, 1};
    const Point2 c{-1, 0};

    const Point2 on{0, -1};
    const Point2 just_outside{0, std::nextafter(-1.0, -2.0)};
    const Point2 just_inside{0, std::nextafter(-1.0, 0.0)};

    CountingExact<DetriaExact>::reset();

    INFO(describe(a, b, c, on));
    REQUIRE(Counted::incircle(a, b, c, on) == Incircle::Cocircular);
    INFO(describe(a, b, c, just_outside));
    REQUIRE(Counted::incircle(a, b, c, just_outside) == Incircle::Outside);
    INFO(describe(a, b, c, just_inside));
    REQUIRE(Counted::incircle(a, b, c, just_inside) == Incircle::Inside);

    // All three went through the exact path. Without this the test would still
    // pass if the filter had widened to swallow them, and it would then be
    // covering the double arithmetic rather than the backend.
    REQUIRE(CountingExact<DetriaExact>::incircle_calls == 3);
}

// ---------------------------------------------------------------------------
// The SIGTRAP hazard, part 1: the diagnosable form.
// ---------------------------------------------------------------------------

// Declared before the unguarded test below on purpose. Catch2 runs the cases in
// a translation unit in declaration order, so if `FilteredKernel::incircle`'s
// normalization regresses, this is the first thing that touches the hazard --
// and it reports it as a test failure naming the signal, instead of the process
// vanishing partway through the binary.
//
// The mechanism: fork, and do the dangerous calls in the child. detria's
// `incircle` checks its counterclockwise precondition under `#ifndef NDEBUG`,
// prints to std::cerr and calls std::raise(SIGTRAP). SIGTRAP's default
// disposition terminates the process, and Catch2 does not install a handler for
// it (it handles SIGILL, SIGSEGV, SIGABRT, SIGFPE, SIGINT and SIGTERM), so in
// the parent the failure is legible through waitpid: WIFSIGNALED with
// WTERMSIG == SIGTRAP.
//
// The child uses _exit, not exit or a Catch2 assertion: it must not run atexit
// handlers, flush the parent's buffered test output a second time, or try to
// report into a Catch2 run it is only a copy of. Findings travel back as an
// exit code and nothing else.
#if TERRAIN_HAS_FORK
namespace {

// Exit codes for the child. Anything else means it died before deciding.
enum ChildResult : int {
    child_ok = 0,
    child_wrong_answer = 1,
};

[[nodiscard]] int run_non_ccw_incircle_probe() {
    const Point2 a{5, 0};
    const Point2 b{0, 5};
    const Point2 c{-5, 0};

    // (a, c, b) is clockwise: the case that reaches the backend and therefore
    // the case that can raise.
    if (DefaultKernel::incircle(a, c, b, Point2{0, 0}) != Incircle::Inside) {
        return child_wrong_answer;
    }
    if (DefaultKernel::incircle(a, c, b, Point2{400, 400}) != Incircle::Outside) {
        return child_wrong_answer;
    }
    if (DefaultKernel::incircle(a, c, b, Point2{3, -4}) != Incircle::Cocircular) {
        return child_wrong_answer;
    }

    // Collinear and coincident triples: these must be answered from the
    // orientation alone and never reach the backend at all, since there is no
    // counterclockwise ordering of them to hand it.
    if (DefaultKernel::incircle(Point2{0, 0}, Point2{1, 1}, Point2{2, 2}, Point2{99, -7}) !=
        Incircle::Cocircular) {
        return child_wrong_answer;
    }
    if (DefaultKernel::incircle(Point2{7, 7}, Point2{7, 7}, Point2{7, 7}, Point2{1, 2}) !=
        Incircle::Cocircular) {
        return child_wrong_answer;
    }

    // A broad sweep over a coarse lattice, where clockwise and collinear
    // triples are common, so the guarded probe covers far more than the three
    // named cases.
    std::mt19937_64 rng{20250302};
    std::uniform_int_distribution<int> coord{-4, 4};
    const auto pick = [&] {
        return Point2{static_cast<double>(coord(rng)), static_cast<double>(coord(rng))};
    };
    for (int i = 0; i < 2000; ++i) {
        const Point2 p = pick();
        const Point2 q = pick();
        const Point2 r = pick();
        const Point2 s = pick();
        if (DefaultKernel::incircle(p, q, r, s) != Reference::incircle(p, q, r, s)) {
            return child_wrong_answer;
        }
    }

    return child_ok;
}

}  // namespace

TEST_CASE("incircle survives non-counterclockwise input, diagnosably", "[predicates][default_kernel][incircle][signal]") {
    // Drain the parent's buffered output first. ctest captures stdout to a
    // file, which makes it block-buffered, so an unflushed buffer is inherited
    // by the child and reappears in the log as a duplicated Catch2 header.
    std::fflush(nullptr);

    const pid_t pid = fork();
    REQUIRE(pid >= 0);

    if (pid == 0) {
        _exit(run_non_ccw_incircle_probe());
    }

    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);

    if (WIFSIGNALED(status)) {
        const int sig = WTERMSIG(status);
        INFO("child terminated by signal " << sig
             << (sig == SIGTRAP
                     ? " (SIGTRAP: FilteredKernel::incircle handed a non-counterclockwise"
                       " triple to detria, whose Debug assertion raises)"
                     : ""));
        FAIL("DefaultKernel::incircle killed the process on a signal");
    }

    REQUIRE(WIFEXITED(status));
    INFO("child exit code " << WEXITSTATUS(status));
    REQUIRE(WEXITSTATUS(status) == child_ok);
}
#endif  // TERRAIN_HAS_FORK

// ---------------------------------------------------------------------------
// The SIGTRAP hazard, part 2: the real thing, in process.
// ---------------------------------------------------------------------------

// READ BEFORE SIMPLIFYING. This test passes by not crashing at least as much as
// by its assertions, and the assertions are supposed to look mild.
//
// It runs unguarded -- no `#ifdef NDEBUG`, no fork, no subprocess -- because
// the Debug build is the only build where the hazard exists, and our asan+ubsan
// CI job is a Debug build. detria's `math::incircle` asserts under
// `#ifndef NDEBUG` that its first three points are counterclockwise, and
// `detail::detriaAssert` writes to std::cerr and calls std::raise(SIGTRAP). It
// is called directly rather than through a macro, so there is no hook, no
// handler to install and no override short of -DNDEBUG.
//
// So if `FilteredKernel::incircle` ever stops normalizing before it reaches the
// backend, this test does not fail -- the process dies on a signal, here, on
// this line. That is the point: this is the canary that proves the normalization
// holds against the backend that actually raises, rather than against
// `CcwCheckingExact`, which is a stand-in that merely records. The forked probe
// above exists to turn that death into a readable report; this one exists to
// make sure the real, unmediated call path is the thing being exercised.
//
// Deleting the clockwise, collinear or coincident case because "the assertion
// is trivial" removes the coverage entirely.
TEST_CASE("incircle never passes a non-counterclockwise triple to the real backend", "[predicates][default_kernel][incircle][signal]") {
    const Point2 a{5, 0};
    const Point2 b{0, 5};
    const Point2 c{-5, 0};
    REQUIRE(RefExact::orient2d(a, b, c) == Orientation::CounterClockwise);

    SECTION("a clockwise triple: the case that reaches the backend") {
        // (a, c, b) is the reverse of a counterclockwise triple, so
        // `FilteredKernel::incircle` must swap before calling `incircle_ccw`.
        REQUIRE(RefExact::orient2d(a, c, b) == Orientation::Clockwise);

        REQUIRE(DefaultKernel::incircle(a, c, b, Point2{0, 0}) == Incircle::Inside);
        REQUIRE(DefaultKernel::incircle(a, c, b, Point2{400, 400}) == Incircle::Outside);
        REQUIRE(DefaultKernel::incircle(a, c, b, Point2{3, -4}) == Incircle::Cocircular);

        // And the answer is the same one the counterclockwise ordering gives:
        // swapping must not negate. Reversing the vertex order of a triangle
        // does not change the circle through its vertices.
        REQUIRE(DefaultKernel::incircle(a, c, b, Point2{0, 0}) ==
                DefaultKernel::incircle(a, b, c, Point2{0, 0}));
        REQUIRE(DefaultKernel::incircle(a, c, b, Point2{400, 400}) ==
                DefaultKernel::incircle(a, b, c, Point2{400, 400}));
    }

    SECTION("a collinear triple: a degenerate circle, answered without the backend") {
        REQUIRE(DefaultKernel::incircle(Point2{0, 0}, Point2{1, 1}, Point2{2, 2},
                                        Point2{99, -7}) == Incircle::Cocircular);
        REQUIRE(DefaultKernel::incircle(Point2{0, 0}, Point2{2, 2}, Point2{1, 1},
                                        Point2{0, 0}) == Incircle::Cocircular);
        // At UTM33 magnitude, where the filter cannot settle the orientation
        // itself and the exact backend decides the collinearity.
        REQUIRE(DefaultKernel::incircle(utm33_offset({0, 0}), utm33_offset({1, 1}),
                                        utm33_offset({2, 2}), utm33_offset({9, 1})) ==
                Incircle::Cocircular);
    }

    SECTION("a coincident triple: collinear, hence cocircular, hence not a legal backend call") {
        REQUIRE(DefaultKernel::incircle(Point2{7, 7}, Point2{7, 7}, Point2{7, 7},
                                        Point2{1, 2}) == Incircle::Cocircular);
        REQUIRE(DefaultKernel::incircle(Point2{7, 7}, Point2{7, 7}, Point2{1, 2},
                                        Point2{7, 7}) == Incircle::Cocircular);
        REQUIRE(DefaultKernel::incircle(a, a, a, a) == Incircle::Cocircular);
    }
}

// The same hazard from the other side: a triple that the *filter* cannot
// classify without the exact backend, and that is clockwise. Here the
// orientation itself comes back from DetriaExact, and the normalization has to
// act on that answer rather than on a filter-resolved one.
TEST_CASE("incircle normalizes correctly when the orientation itself needed the backend", "[predicates][default_kernel][incircle][signal][adversarial]") {
    std::mt19937_64 rng{777};
    int clockwise_seen = 0;

    for (int i = 0; i < 300; ++i) {
        const auto t = terrain::test::near_degenerate_triple(rng);
        if (t.expected != Orientation::Clockwise) {
            continue;
        }
        ++clockwise_seen;

        // RefExact's incircle bound is 2^25, and this family's coordinates
        // exceed it, so the reference kernel cannot referee. What is asserted
        // is the property that does not need an oracle: the call returns at
        // all, and returns the same thing as the explicitly swapped call.
        const Point2 d{1.0, 1.0};
        INFO(describe(t.a, t.b, t.c, d));
        REQUIRE(DefaultKernel::incircle(t.a, t.b, t.c, d) ==
                DefaultKernel::incircle(t.a, t.c, t.b, d));
    }

    // Guards against the loop silently never running the interesting branch.
    REQUIRE(clockwise_seen > 0);
}
