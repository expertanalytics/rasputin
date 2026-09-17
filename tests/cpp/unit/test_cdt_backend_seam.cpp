// The backend seam: terrain/cdt/triangulate.hpp, the CdtBackend concept, and a
// fake backend standing in for detria.
//
// NOT invariant-critical, and no mutation round is spent here. A mutant of this
// file is a build failure, not a wrong answer, and mutation spend on a
// link-level proof buys nothing.
//
// READ THE LINK LINE BEFORE CHANGING ANYTHING HERE. This target registers with
// the PLAIN add_terrain_test helper: terrain_headers, the test support headers
// and Catch2, and NOTHING ELSE -- no terrain_cdt, no terrain_predicates, no
// detria. That is the test. terrain_cdt contributes exactly one symbol,
// DetriaBackend::triangulate; everything else in terrain/cdt is header-only.
// So if the generic entry point ever CALLS the shipped backend -- names it in a
// default, forwards to it on a fallback path, instantiates it to "check" the
// concept -- this file fails to LINK, and the failure is the finding. No
// assertion inside a test body can establish what an absent symbol establishes
// for free.
//
// Precisely, so the claim is not larger than the mechanism: merely INCLUDING
// detria_backend.hpp from triangulate.hpp would not fail this link, because
// that header is a declaration and a static_assert. What catches that is the
// #error guard below plus tools/check_detria_boundary.py. The link line catches
// the call.
//
// The same constraint is why the Pslg fixture is built with FastKernel. The
// filtered kernel's exact backend lives in the compiled terrain_predicates
// target; FastKernel is header-only. On a well-separated unit square the two
// cannot disagree -- which is exactly the condition under which the cheaper
// kernel is legitimate, and also the fact that lets this file stay unlinked.
// cdt_cases.hpp names no kernel for this reason.
//
// If someone "fixes" the CMake registration to use add_terrain_cdt_test, the
// suite still passes and stops meaning anything. The CMakeLists comment says so
// too.

#include <catch2/catch_test_macros.hpp>

#include <cdt_cases.hpp>

#include <terrain/cdt/result.hpp>
#include <terrain/cdt/triangulate.hpp>
#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/predicates/kernel.hpp>

// The vendored detria.hpp is an implementation detail of two translation units,
// src/predicates/detria_exact.cpp and src/cdt/detria_backend.cpp. If a build
// change ever puts lib/detria on a broad include path, and something in the cdt
// headers starts including it, the guard macro from the pinned version is
// defined here and this stops the build with a readable message.
#ifdef DETRIA_HPP_INCLUDED
#error "terrain/cdt headers must not include detria.hpp -- see src/cdt/detria_backend.cpp"
#endif

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::Pslg;
using terrain::TriangleIndices;
using terrain::cdt::CdtBackend;
using terrain::cdt::CdtOptions;
using terrain::cdt::CdtOutcome;
using terrain::cdt::CdtStatus;
using terrain::kEdgeMask01;
using terrain::pred::FastKernel;
using terrain::test::build_fixture;
using terrain::test::square_domain;

namespace {

// Ignores its arguments and returns a canned two-triangle mesh. It can do that
// only because IndexedMesh2's constructor is public -- a fake backend that
// cannot construct a mesh cannot prove the seam is real, which is the decisive
// argument for that constructor being public rather than befriended.
struct FakeCdtBackend {
    static CdtOutcome triangulate(const Pslg&, const CdtOptions&) {
        std::vector<Point2> vertices = {Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0},
                                        Point2{0.0, 1.0}};
        std::vector<TriangleIndices> triangles = {TriangleIndices{0, 1, 3},
                                                  TriangleIndices{1, 2, 3}};
        std::vector<std::uint8_t> masks = {kEdgeMask01, kEdgeMask01};
        return CdtOutcome{CdtStatus::Ok, {},
                          IndexedMesh2{std::move(vertices), std::move(triangles),
                                       std::move(masks)}};
    }
};

// A backend that fails, so the seam is exercised on both arms: a status the
// caller must read, a message, and no mesh.
struct FailingCdtBackend {
    static CdtOutcome triangulate(const Pslg&, const CdtOptions&) {
        return CdtOutcome{CdtStatus::NotNoded, "two constraints cross", {}};
    }
};

// Right name, wrong return type. A concept nothing fails is a concept that
// constrains nothing, so this type exists to fail it.
struct NotACdtBackend {
    static bool triangulate(const Pslg&, const CdtOptions&) { return true; }
};

// Right signature, but only callable on an instance. The concept is spelled
// with a qualified static call on purpose -- matching GeometryKernel -- so that
// statelessness is part of the signature rather than a separate purity rule a
// future model would not be covered by. parallel_refinement.md has every
// refinement thread sharing one constraint set, and a stateful backend object
// is where that stops being true.
struct StatefulCdtBackend {
    CdtOutcome triangulate(const Pslg&, const CdtOptions&) const { return CdtOutcome{}; }
};

}  // namespace

TEST_CASE("the concept accepts a well-formed backend and rejects the others", "[cdt][seam]") {
    STATIC_REQUIRE(CdtBackend<FakeCdtBackend>);
    STATIC_REQUIRE(CdtBackend<FailingCdtBackend>);
    STATIC_REQUIRE_FALSE(CdtBackend<NotACdtBackend>);
    STATIC_REQUIRE_FALSE(CdtBackend<StatefulCdtBackend>);
    STATIC_REQUIRE_FALSE(CdtBackend<int>);
}

// The entry point checks the seam's postconditions in debug, for every backend
// including ones not written yet. Both arms are exercised below, and that is
// the only test either assert gets: a test that deliberately violated one would
// abort the asan+ubsan job rather than fail.
//
// THE POLARITY OF THE SECOND ASSERT. The postconditions are
//
//     out.ok() == !out.mesh.empty()      no Ok without a mesh, no mesh with an error
//     out.ok() ==  out.message.empty()   success says nothing, failure says something
//
// docs/increments/04-cdt.md spells the second one `out.ok() != out.message.empty()`,
// which is inverted: on success both sides are true and on failure both are
// false, so the assert fires on EVERY call. Transcribed literally it aborts the
// first test that reaches it. Measured, not deduced -- the reference
// implementation this suite was validated against did exactly that.
TEST_CASE("the generic entry point returns what its backend returned", "[cdt][seam]") {
    const Pslg pslg = build_fixture<FastKernel>(square_domain());

    const CdtOutcome out = terrain::cdt::triangulate<FakeCdtBackend>(pslg);

    REQUIRE(out.ok());
    CHECK(out.status == CdtStatus::Ok);
    CHECK(out.message.empty());
    CHECK(out.mesh.triangle_count() == 2);
}

TEST_CASE("a failing backend's status and message reach the caller", "[cdt][seam]") {
    const Pslg pslg = build_fixture<FastKernel>(square_domain());

    const CdtOutcome out = terrain::cdt::triangulate<FailingCdtBackend>(pslg, CdtOptions{});

    CHECK_FALSE(out.ok());
    CHECK(out.status == CdtStatus::NotNoded);
    CHECK_FALSE(out.message.empty());
    CHECK(out.mesh.empty());
}

// The options parameter is defaulted, so a caller with nothing to say about
// Delaunay writes nothing. CdtOptions exists as a struct -- rather than a bare
// bool -- so that adding a second option later does not change the CdtBackend
// signature and therefore does not touch every backend.
TEST_CASE("options default to a Delaunay triangulation", "[cdt][seam]") {
    STATIC_REQUIRE(CdtOptions{}.delaunay);

    const Pslg pslg = build_fixture<FastKernel>(square_domain());
    CHECK(terrain::cdt::triangulate<FakeCdtBackend>(pslg).ok());
}
