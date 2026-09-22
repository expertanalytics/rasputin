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
//
// AMENDED AT INCREMENT 5c, and this file is NOT in 05c-noder-wiring.md's list of
// suites the churn touches -- the list names the backend, the property and the
// constraint-edge suites. It has to change all the same: CdtBackend is spelled
// in terms of the entry point's parameter, so retyping triangulate retypes the
// concept, and every fake backend here declares that parameter. A suite of fakes
// whose signature no longer matches the concept fails STATIC_REQUIRE, which is
// the finding rather than a compile accident.
//
// The fixtures are noded under FastKernel, which the link line permits: node<K>
// is header-only, exactly as the validator is, so nothing compiled is dragged in
// by the noder any more than by build_fixture.
//
// IT ALSO HOUSES THE ChainedGraph CASE, and that placement is deliberate. The
// design says in one place that test_cdt_constrained_edges.cpp is UNCHANGED --
// the concrete return on templating ConstraintEdgeSet -- and in another that it
// is amended with a case for the concept. The first is the load-bearing claim:
// not one of that file's twenty-odd hand-built Pslg fixtures needs adapting,
// and that is only demonstrated by not touching it. The concept itself is a
// seam between two graph types and one consumer, which is this file's subject.

#include <catch2/catch_test_macros.hpp>

#include <cdt_cases.hpp>

#include <terrain/cdt/constrained_edges.hpp>
#include <terrain/cdt/result.hpp>
#include <terrain/cdt/triangulate.hpp>
#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/noded_pslg.hpp>
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

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::NodedPslg;
using terrain::Point2;
using terrain::Pslg;
using terrain::TriangleIndices;
using terrain::cdt::ChainedGraph;
using terrain::cdt::CdtBackend;
using terrain::cdt::ConstraintEdgeSet;
using terrain::cdt::CdtOptions;
using terrain::cdt::CdtOutcome;
using terrain::cdt::CdtStatus;
using terrain::kEdgeMask01;
using terrain::pred::FastKernel;
using terrain::test::build_fixture;
using terrain::test::node_fixture;
using terrain::test::square_domain;

namespace {

// Ignores its arguments and returns a canned two-triangle mesh. It can do that
// only because IndexedMesh2's constructor is public -- a fake backend that
// cannot construct a mesh cannot prove the seam is real, which is the decisive
// argument for that constructor being public rather than befriended.
struct FakeCdtBackend {
    static CdtOutcome triangulate(const NodedPslg&, const CdtOptions&) {
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
    static CdtOutcome triangulate(const NodedPslg&, const CdtOptions&) {
        // BackendFailure rather than NotNoded: as of 5c the latter is a
        // self-check no NodedPslg can reach, and a fake handing one back would
        // read as a claim that it still can.
        return CdtOutcome{CdtStatus::BackendFailure, "the backend gave up", {}};
    }
};

// Right name, wrong return type. A concept nothing fails is a concept that
// constrains nothing, so this type exists to fail it.
struct NotACdtBackend {
    static bool triangulate(const NodedPslg&, const CdtOptions&) { return true; }
};

// Right signature, but only callable on an instance. The concept is spelled
// with a qualified static call on purpose -- matching GeometryKernel -- so that
// statelessness is part of the signature rather than a separate purity rule a
// future model would not be covered by. parallel_refinement.md has every
// refinement thread sharing one constraint set, and a stateful backend object
// is where that stops being true.
struct StatefulCdtBackend {
    CdtOutcome triangulate(const NodedPslg&, const CdtOptions&) const { return CdtOutcome{}; }
};

// The pre-5c signature, kept alive as a type that must now FAIL the concept.
// This is the architectural product of 5b and 5c stated where a compiler checks
// it: un-noded input is unrepresentable at the entry point rather than
// diagnosed inside it, so a backend that still offers to take a Pslg is not a
// backend. Without this, the retype could be reverted and every other
// assertion in the file would stay green.
struct UnNodedCdtBackend {
    static CdtOutcome triangulate(const Pslg&, const CdtOptions&) { return CdtOutcome{}; }
};

}  // namespace

TEST_CASE("the concept accepts a well-formed backend and rejects the others", "[cdt][seam]") {
    STATIC_REQUIRE(CdtBackend<FakeCdtBackend>);
    STATIC_REQUIRE(CdtBackend<FailingCdtBackend>);
    STATIC_REQUIRE_FALSE(CdtBackend<NotACdtBackend>);
    STATIC_REQUIRE_FALSE(CdtBackend<StatefulCdtBackend>);
    STATIC_REQUIRE_FALSE(CdtBackend<UnNodedCdtBackend>);
    STATIC_REQUIRE_FALSE(CdtBackend<int>);
}

// ---------------------------------------------------------------------------
// ChainedGraph: the other half of the signature change
// ---------------------------------------------------------------------------
//
// ConstraintEdgeSet reads exactly three members -- chains().size(),
// edge_count(c), indices_of(c) -- and Pslg and NodedPslg both have all three.
// Templating the constructor on that is what keeps the mask computation
// testable on hand-built Pslg input forever, and the two directions are
// deliberate: the HELPER is structural, the ENTRY POINT is not. A
// ChainedGraph-templated triangulate would restore Pslg as a legal argument and
// undo the increment while leaving every test green, which is why the case
// above asserts the opposite for CdtBackend.
TEST_CASE("the chained-graph concept holds for both graph types", "[cdt][seam][concept]") {
    STATIC_REQUIRE(ChainedGraph<Pslg>);
    STATIC_REQUIRE(ChainedGraph<NodedPslg>);
    // A concept nothing fails constrains nothing.
    STATIC_REQUIRE_FALSE(ChainedGraph<IndexedMesh2>);
    STATIC_REQUIRE_FALSE(ChainedGraph<int>);
}

// The compiler checking the structural equality NodedPslg was given on purpose:
// the same fixture, through both types, yields the same NUMBER of constraint
// edges. Not the same keys -- node ids are the sorted node set's order, so the
// two index spaces differ, and an assertion that they agreed would be asserting
// that the noder had done nothing.
TEST_CASE("the constraint edge set is constructible from either graph", "[cdt][seam][concept]") {
    const Pslg pslg = build_fixture<FastKernel>(square_domain());
    const NodedPslg noded = node_fixture<FastKernel>(square_domain());

    const ConstraintEdgeSet from_pslg{pslg};
    const ConstraintEdgeSet from_noded{noded};

    CHECK(from_pslg.size() == 4);
    CHECK(from_noded.size() == from_pslg.size());
    for (std::size_t c = 0; c < noded.chains().size(); ++c) {
        const std::span<const std::uint32_t> idx = noded.indices_of(c);
        for (std::size_t k = 0; k < noded.edge_count(c); ++k) {
            CHECK(from_noded.contains(idx[k], idx[(k + 1) % idx.size()]));
        }
    }
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
    const NodedPslg pslg = node_fixture<FastKernel>(square_domain());

    const CdtOutcome out = terrain::cdt::triangulate<FakeCdtBackend>(pslg);

    REQUIRE(out.ok());
    CHECK(out.status == CdtStatus::Ok);
    CHECK(out.message.empty());
    CHECK(out.mesh.triangle_count() == 2);
}

TEST_CASE("a failing backend's status and message reach the caller", "[cdt][seam]") {
    const NodedPslg pslg = node_fixture<FastKernel>(square_domain());

    const CdtOutcome out = terrain::cdt::triangulate<FailingCdtBackend>(pslg, CdtOptions{});

    CHECK_FALSE(out.ok());
    CHECK(out.status == CdtStatus::BackendFailure);
    CHECK_FALSE(out.message.empty());
    CHECK(out.mesh.empty());
}

// The options parameter is defaulted, so a caller with nothing to say about
// Delaunay writes nothing. CdtOptions exists as a struct -- rather than a bare
// bool -- so that adding a second option later does not change the CdtBackend
// signature and therefore does not touch every backend.
TEST_CASE("options default to a Delaunay triangulation", "[cdt][seam]") {
    STATIC_REQUIRE(CdtOptions{}.delaunay);

    const NodedPslg pslg = node_fixture<FastKernel>(square_domain());
    CHECK(terrain::cdt::triangulate<FakeCdtBackend>(pslg).ok());
}
