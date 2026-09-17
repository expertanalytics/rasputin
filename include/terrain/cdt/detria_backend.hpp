#pragma once

// The shipped CDT backend: a declaration, and nothing else.
//
// detria.hpp is an implementation detail of src/cdt/detria_backend.cpp, the
// second and last translation unit in the project that includes it. This header
// must never include it -- lib/detria is PRIVATE on the terrain_cdt target, the
// suites carry an #ifdef DETRIA_HPP_INCLUDED guard, and
// tools/check_detria_boundary.py is the third line of defence.
//
// The static_assert sits in the header so the shipped backend is re-checked
// against the concept in every translation unit that includes it: the cheapest
// possible version of "the seam is not fake".

#include <terrain/cdt/result.hpp>
#include <terrain/cdt/triangulate.hpp>
#include <terrain/core/pslg.hpp>

namespace terrain::cdt {

struct DetriaBackend {
    [[nodiscard]] static CdtOutcome triangulate(const Pslg& pslg, const CdtOptions& options);
};

static_assert(CdtBackend<DetriaBackend>);

}  // namespace terrain::cdt
