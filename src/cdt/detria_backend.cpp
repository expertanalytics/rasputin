// The second and last translation unit in this repository that includes
// detria.hpp. See lib/detria/README.md for why that is a rule and how it is
// enforced.

#include <terrain/cdt/detria_backend.hpp>

#include <terrain/cdt/constrained_edges.hpp>
#include <terrain/cdt/result.hpp>
#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/noded_pslg.hpp>

#include <detria.hpp>

#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace terrain::cdt {

namespace {

using Triangulation = detria::Triangulation<Point2, std::uint32_t>;

// Every detria::TriangulationError maps to exactly one CdtStatus, by an
// exhaustive switch and NEVER BY A CAST -- increment 1 established that rule for
// detria's orientation enums and mutation-tested the failures a cast produces.
//
// The switch has a second job here: `return BackendFailure` sits AFTER it rather
// than in a `default:` arm, so -Wswitch under -Werror turns a re-pin that adds
// an enumerator into a compile error while the runtime path stays total.
//
// The MalformedInput rows are the interesting ones: each is something the PSLG
// validator already excludes, so MalformedInput firing means either the Pslg did
// not come from the builder or this wrapper mis-built a span. It is a self-check
// with a name, not a user-facing status.
[[nodiscard]] CdtStatus map_error(detria::TriangulationError error) noexcept {
    switch (error) {
        case detria::TriangulationError::NoError:
            return CdtStatus::Ok;

        case detria::TriangulationError::TriangulationNotStarted:
        case detria::TriangulationError::LessThanThreePoints:
        case detria::TriangulationError::NonFinitePositionFound:
        case detria::TriangulationError::PolylineIndexOutOfBounds:
        case detria::TriangulationError::PolylineTooShort:
            return CdtStatus::MalformedInput;

        case detria::TriangulationError::DuplicatePointsFound:
        case detria::TriangulationError::PolylineDuplicateConsecutivePoints:
        case detria::TriangulationError::AllPointsAreCollinear:
            return CdtStatus::DegenerateGeometry;

        case detria::TriangulationError::PointOnConstrainedEdge:
        case detria::TriangulationError::ConstrainedEdgeIntersection:
            return CdtStatus::NotNoded;

        case detria::TriangulationError::HoleNotInsideOutline:
        case detria::TriangulationError::StackedPolylines:
        case detria::TriangulationError::EdgeWithDifferentConstrainedTypes:
            return CdtStatus::InvalidTopology;

        case detria::TriangulationError::AssertionFailed:
            return CdtStatus::BackendFailure;
    }
    // Unreachable at the pinned SHA: the enumeration is closed and every
    // enumerator is named above.
    return CdtStatus::BackendFailure;
}

}  // namespace

// LIFETIMES. addOutline, addHole and setPoints store a span and DO NOT COPY, so
// the caller's buffers must outlive triangulate(). This function hands detria
// exactly two kinds of span -- pslg.vertices() and pslg.indices_of(c) -- plus
// individual std::uint32_t values to setConstrainedEdge, which are copied into
// detria's own vector. NO BUFFER CONSTRUCTED INSIDE THIS FUNCTION IS EVER PASSED
// TO A DETRIA CALL, and the Triangulation is a local in the same scope as the
// const NodedPslg& parameter, so it cannot outlive the argument.
//
// indices_of(c) goes STRAIGHT into addOutline/addHole with no closing index:
// detria closes an open polyline itself (createConstrainedEdges seeds
// prevVertexIdx with polyline.back()), so no scratch buffer exists here to
// dangle. The backend suite's characterisation test pins that behaviour; if a
// re-pin breaks it, the fix is one scratch buffer for the whole triangulation
// with sub-spans handed out, and nothing else changes.
CdtOutcome DetriaBackend::triangulate(const NodedPslg& pslg, const CdtOptions& options) {
    Triangulation tri;  // by value, on the stack: move-only, and a member would
                        // give this backend the state the concept forbids
    tri.setPoints(pslg.vertices());

    for (std::size_t c = 0; c < pslg.chains().size(); ++c) {
        const std::span<const std::uint32_t> idx = pslg.indices_of(c);
        switch (pslg.chains()[c].role) {
            case ChainRole::Outer:
                tri.addOutline(idx);
                break;
            case ChainRole::Hole:
                tri.addHole(idx);
                break;
            case ChainRole::Breakline:
                // EDGE BY EDGE, never as a polyline, and never through
                // addPolylineAutoDetectType: we know the chain's role from its
                // ChainRole, and letting the library re-derive it geometrically
                // is a silent-wrong-answer channel -- a closed-loop breakline
                // auto-detected as a hole would carve a void out of the domain.
                // edge_count(c) is the authority on how many edges an open
                // chain has.
                for (std::size_t k = 0; k < pslg.edge_count(c); ++k) {
                    tri.setConstrainedEdge(idx[k], idx[k + 1]);
                }
                break;
        }
    }

    if (!tri.triangulate(options.delaunay)) {
        // getErrorMessage() builds through <sstream>, so it is called ONLY on
        // this path. On success `message` stays empty rather than carrying
        // detria's own success prose.
        return CdtOutcome{map_error(tri.getError()), tri.getErrorMessage(), {}};
    }

    // Built AFTER the call, not before: it is the only buffer this function
    // owns, and building it on the far side means it cannot be confused for
    // something detria is holding. It is also not built at all on the failure
    // path above.
    const ConstraintEdgeSet constraints{pslg};

    std::vector<TriangleIndices> triangles;
    std::vector<std::uint8_t> masks;
    // An upper bound -- it counts interior, hole and convex-hull triangles --
    // so the reserve never under-allocates.
    triangles.reserve(tri.getMaxNumTriangles());
    masks.reserve(tri.getMaxNumTriangles());

    // forEachTriangle, NOT forEachTriangleOfEveryLocation and NOT
    // forEachTriangleOfLocation: it filters to TriangleLocation::Interior, so
    // hole triangles and convex-hull triangles never reach the mesh. That is
    // what makes the suite's Euler counts meaningful.
    //
    // cwTriangles = false, WHICH IS NOT THE DEFAULT. Under the default v1 and
    // v2 arrive swapped, which both reverses the winding and transposes mask
    // bits 0 and 2. One argument, two invariants.
    tri.forEachTriangle(
        [&](detria::Triangle<std::uint32_t> t) {
            const TriangleIndices v{t.x, t.y, t.z};
            triangles.push_back(v);
            masks.push_back(constrained_mask(constraints, v));
        },
        /*cwTriangles=*/false);

    return CdtOutcome{CdtStatus::Ok, {},
                      IndexedMesh2{{pslg.vertices().begin(), pslg.vertices().end()},
                                   std::move(triangles), std::move(masks)}};
}

}  // namespace terrain::cdt
