#pragma once

// The backend seam: the CdtBackend concept and the one generic entry point.
//
// THIS HEADER NAMES NO CONCRETE BACKEND. It does not include
// detria_backend.hpp, it does not default its template parameter to
// DetriaBackend, and it never calls one on a fallback path -- which is what
// lets tests/cpp/unit/test_cdt_backend_seam.cpp link without terrain_cdt at
// all. A seam that survives only because both halves are always present is not
// a seam; that link line is the proof.
//
// Spelled with a QUALIFIED STATIC CALL, matching GeometryKernel: the concept
// itself requires models to be usable without an instance, so statelessness is
// part of the signature rather than a separate purity rule a future model would
// not be covered by. A backend that needs scratch storage allocates it inside
// the call -- parallel_refinement.md has every refinement thread sharing one
// constraint set, and a stateful backend object is where that stops being true.
//
// SEMANTIC OBLIGATIONS OF THE CONCEPT, stated here because a concept cannot
// express them and a model that violates one produces a wrong mesh rather than
// a compile error:
//
//   1. The returned vertex array BEGINS WITH pslg.vertices(), element-wise and
//      in order. Additional vertices may only be appended.
//   2. Only in-domain triangles are returned: nothing inside a Hole, nothing
//      outside every Outer.
//   3. Triangles are counterclockwise under an exact orientation predicate.
//   4. The mask convention of terrain/cdt/constrained_edges.hpp is the
//      backend's to honour.
//   5. The backend decides orientation and incircle AT LEAST AS EXACTLY as the
//      kernel that validated the Pslg. Increment 3 pinned that a Pslg is valid
//      with respect to a kernel; a backend using naive doubles could disagree
//      with the validator about a sliver's winding. No kernel template
//      parameter is threaded through the CDT: detria does not accept one, and
//      threading a dead parameter to document an obligation is how a parameter
//      later acquires a meaning nobody checked.

#include <terrain/cdt/result.hpp>
#include <terrain/core/pslg.hpp>

#include <cassert>
#include <concepts>

namespace terrain::cdt {

template <typename B>
concept CdtBackend = requires(const Pslg& pslg, const CdtOptions& options) {
    { B::triangulate(pslg, options) } -> std::same_as<CdtOutcome>;
};

// NOT A FORWARDING WRAPPER. It checks the seam's postconditions, in debug, for
// every backend including ones not written yet -- which is the difference
// between a concept that names a signature and a concept that means something.
template <CdtBackend B>
[[nodiscard]] CdtOutcome triangulate(const Pslg& pslg, const CdtOptions& options = {}) {
    CdtOutcome out = B::triangulate(pslg, options);

    // No Ok with an empty mesh, no mesh alongside an error. The first half is
    // not hypothetical: a triangulation with points but no outline succeeds and
    // yields zero interior triangles, so "forgot to add the outline" is a
    // silent success without this.
    assert(out.ok() == !out.mesh.empty());
    // Success says nothing, failure says something. NOTE THE OPERATOR: with !=
    // the assertion is false in both cases -- on success both sides are true,
    // on failure both are false -- and fires on every call.
    assert(out.ok() == out.message.empty());

    return out;
}

}  // namespace terrain::cdt
