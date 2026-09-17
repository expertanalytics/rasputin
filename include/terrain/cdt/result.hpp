#pragma once

// The CDT's result channel: a status, a message, a mesh -- and no exception
// crossing the backend seam.
//
// Why a status rather than std::expected or a throw: std::expected is C++23 and
// this project is C++20; a future backend may be a C library, and a status is
// the only channel a C library can speak; and the pybind11 layer must translate
// to a Python exception regardless, so an exception here would be caught and
// rebuilt one frame later.
//
// THIS IS A DIFFERENT CHANNEL FROM PslgBuilder'S DIAGNOSTICS VECTOR, and the
// difference is deliberate. The builder validates input authored outside the
// process, where several chains are commonly wrong at once for one underlying
// reason, so reporting one failure turns one fix into N round trips. This is a
// backend seam: one call to a vendored library, which stops at its first
// problem and has exactly one thing to say about it. Collecting a vector here
// would be collecting a vector of length one.

#include <terrain/core/indexed_mesh.hpp>

#include <string>
#include <string_view>

namespace terrain::cdt {

// Grouped by WHAT THE CALLER SHOULD DO, not by where in the backend the failure
// arose: the backend's own diagnosis goes into `message`, so nothing is lost by
// the grouping and the offending indices survive.
enum class CdtStatus : int {
    Ok = 0,
    NotRun,              // default-constructed outcome; no backend was called
    NotNoded,            // constraints cross, or a vertex lies on a constraint
    DegenerateGeometry,  // coincident points, repeated indices, all-collinear input
    InvalidTopology,     // hole outside every outline, stacked or overlapping rings
    MalformedInput,      // a precondition the Pslg validator already excludes
    BackendFailure,      // the library failed in a way we do not classify
};

[[nodiscard]] constexpr std::string_view describe(CdtStatus status) noexcept {
    switch (status) {
        case CdtStatus::Ok:
            return "the triangulation succeeded";
        case CdtStatus::NotRun:
            return "no triangulation was attempted";
        case CdtStatus::NotNoded:
            return "the constraints are not noded: two constraints cross, or a vertex lies "
                   "exactly on a constraint edge it is not part of -- the noder is the fix";
        case CdtStatus::DegenerateGeometry:
            return "degenerate geometry: coincident vertices, a repeated consecutive index, "
                   "or an all-collinear point set";
        case CdtStatus::InvalidTopology:
            return "invalid ring topology: a hole outside every outline, a hole containing "
                   "one, stacked rings, or a hole sharing an edge with its outline";
        case CdtStatus::MalformedInput:
            // A self-check with a name rather than a user-facing status: every
            // condition behind it is one the PSLG validator already excludes,
            // so it firing means the Pslg did not come from the builder or the
            // wrapper mis-built a span.
            return "malformed input the PSLG validator should already have excluded: this is "
                   "a bug in the wrapper or a Pslg that did not come from the builder";
        case CdtStatus::BackendFailure:
            return "the triangulation backend failed in a way we do not classify";
    }
    // Unreachable: the enumeration is closed. Spelled after the switch rather
    // than as a `default:` arm so that -Wswitch under -Werror turns a new
    // enumerator into a compile error while the runtime path stays total.
    return "unknown triangulation status";
}

// One member, and that is not an oversight. `delaunay` is forwarded to the
// backend and earns its place by being the knob the property suite turns to
// prove the Delaunay property is earned rather than vacuous. The struct exists
// -- rather than a bare bool -- so that adding a second option later does not
// change the CdtBackend signature and therefore does not touch every backend.
// Nothing else goes in it: not a kernel, not a tolerance (this project has
// none), not a "compute the mask" flag.
struct CdtOptions {
    bool delaunay{true};
};

// status's default member initializer is NotRun, A FAILURE VALUE, for the same
// reason Chain::role defaults to Breakline: a default-constructed outcome that
// reaches a reader must not be able to claim success.
//
// On Ok, `message` is EMPTY -- not the backend's own success prose. The backend
// formats a message only on the failure path, and message.empty() on Ok is a
// pinned invariant so that the mutant which formats on the happy path is killed
// by an assertion rather than by a profiler.
struct CdtOutcome {
    CdtStatus status{CdtStatus::NotRun};
    std::string message;
    IndexedMesh2 mesh;

    [[nodiscard]] bool ok() const noexcept { return status == CdtStatus::Ok; }
};

}  // namespace terrain::cdt
