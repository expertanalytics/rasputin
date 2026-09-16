#pragma once

// The kernel the engine asks its geometric questions of.
//
// Everything downstream -- the CDT, the flip loop, the refinement policies --
// should name `DefaultKernel` rather than `FilteredKernel<DetriaExact>`, so
// that swapping the exact backend is a one-line change here. `GeometryKernel`
// is what makes that swap safe: the triangulation is templated on the concept,
// not on this alias.

#include <terrain/predicates/detria_exact.hpp>
#include <terrain/predicates/kernel.hpp>

namespace terrain::pred {

// An alias, not a new type. A distinct `struct DefaultKernel` forwarding to the
// filter would also model `GeometryKernel`, and would silently acquire a place
// for policy to accumulate; a backend swap is the only thing this name is for.
using DefaultKernel = FilteredKernel<DetriaExact>;

static_assert(GeometryKernel<DefaultKernel>);

}  // namespace terrain::pred
