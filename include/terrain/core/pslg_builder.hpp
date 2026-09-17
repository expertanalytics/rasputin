#pragma once

// The PSLG builder and its validator: the one place a constraint set is
// checked, so that terrain::Pslg can be the type whose existence is the proof.
//
// A SEPARATE BUILDER TYPE, NOT A VALIDATING FACTORY. A factory taking
// (vertices, indices, chains) would require the caller to have already built
// the flat index buffer and computed every begin/count pair, which is exactly
// the arithmetic this increment exists to do once and correctly. The builder
// owns it: add_chain appends to the flat buffer and records
// begin = buffer.size() - indices.size(), so begin/count are never written by a
// caller and an off-by-one in them is unrepresentable rather than validated.
//
// THE FAILURE CHANNEL IS DIAGNOSTICS, NOT A STATUS AND NOT AN EXCEPTION, and
// the split from increment 2's channel is by error CLASS, not by taste.
// PointRing's and IndexedRing's constructors throw std::invalid_argument
// because a malformed view is a PROGRAMMER error on a type built per query,
// where a returned status would be ignored at ten call sites. A malformed
// constraint set is a DATA error arriving from outside the process, where a
// status must be handled and a diagnostic must be legible. The obvious
// "simplification" is to make one of them match the other; do not.
//
// Inputs at this boundary are wrong in bulk -- a whole layer digitised
// clockwise, a whole file with a shifted index base -- and a channel that
// reports the first failure turns one fix into N round trips through a pipeline
// whose cheapest stage is a DEM decode. So the validator is exhaustive:
// EXHAUSTIVENESS IS A PROPERTY OF DIAGNOSIS, NOT OF EXECUTION. Stages 1-6 never
// early-return. Stage 0 does, and that is specified rather than a lapse: once
// the sizes exceed uint32 the begin/count fields have ALREADY been truncated by
// the builder, so every later diagnostic would be noise attributed to chains
// that may be perfectly well-formed. The validator never stops because it has
// found enough errors, only when continuing would fabricate them.
//
// No exceptions of our own cross this boundary. std::vector and std::string may
// still throw bad_alloc; nothing else does.
//
// NO DEDUP ANYWHERE. This header includes no <unordered_map>, <unordered_set>
// or <map> and specifies no key, because there is no legal key: std::hash<Point2>
// is for finding a KNOWN point in a hashed container, and anything dedup-shaped
// in this project keys on snapped integer coordinates, which do not exist until
// increment 5. Dedup is also a mutation of caller-declared topology, and this
// increment's posture is check, do not fix -- the same posture that refuses to
// reverse a wrongly wound ring.

#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace terrain {

enum class PslgError : int {
    NoOuterChain,
    ChainTooShort,
    IndexOutOfRange,
    NonFiniteVertex,
    StoredClosure,
    WrongWinding,
    DegenerateRing,
    VertexCountOverflow,
};

inline constexpr std::uint32_t kNoChain = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::uint32_t kNoVertex = std::numeric_limits<std::uint32_t>::max();

struct PslgDiagnostic {
    PslgError error;
    std::uint32_t chain{kNoChain};  // index into chains(), or kNoChain

    // A VALID index into vertices(), or kNoVertex -- never an offending
    // out-of-range value, and safe to dereference whenever it is set.
    //
    // IndexOutOfRange therefore leaves it kNoVertex and puts the offending
    // value in the message. That is a choice, not an omission: a consumer doing
    // the natural thing with a populated field -- p.vertices()[d.vertex] in a
    // reporting tool, a Python binding mapping it back to a source feature --
    // would perform exactly the out-of-bounds read that stage exists to
    // prevent, on the one diagnostic where the value is guaranteed to be out of
    // bounds. Do not "fix" this later.
    std::uint32_t vertex{kNoVertex};

    // Human-readable, std::formatted at diagnosis time, carrying the offending
    // VALUES that do not fit the two index fields. NOT machine-parsed, and no
    // test asserts its wording: pinning prose produces a suite that fails on
    // improvements to its own error messages.
    std::string message;
};

struct PslgBuildResult {
    std::optional<Pslg> pslg;  // engaged iff diagnostics.empty()
    std::vector<PslgDiagnostic> diagnostics;

    [[nodiscard]] bool ok() const noexcept { return pslg.has_value(); }
};

namespace detail {

[[nodiscard]] constexpr std::string_view pslg_error_name(PslgError e) noexcept {
    switch (e) {
        case PslgError::NoOuterChain: return "NoOuterChain";
        case PslgError::ChainTooShort: return "ChainTooShort";
        case PslgError::IndexOutOfRange: return "IndexOutOfRange";
        case PslgError::NonFiniteVertex: return "NonFiniteVertex";
        case PslgError::StoredClosure: return "StoredClosure";
        case PslgError::WrongWinding: return "WrongWinding";
        case PslgError::DegenerateRing: return "DegenerateRing";
        case PslgError::VertexCountOverflow: return "VertexCountOverflow";
    }
    return "UnknownPslgError";
}

[[nodiscard]] constexpr std::string_view chain_role_name(ChainRole r) noexcept {
    switch (r) {
        case ChainRole::Outer: return "Outer";
        case ChainRole::Hole: return "Hole";
        case ChainRole::Breakline: return "Breakline";
    }
    return "UnknownChainRole";
}

[[nodiscard]] constexpr std::string_view orientation_name(pred::Orientation o) noexcept {
    switch (o) {
        case pred::Orientation::Clockwise: return "clockwise";
        case pred::Orientation::Collinear: return "collinear";
        case pred::Orientation::CounterClockwise: return "counterclockwise";
    }
    return "unknown";
}

[[nodiscard]] constexpr bool is_finite(const Point2& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y);
}

// True iff a vertex buffer of `vertex_count` points and a flat index buffer of
// `index_count` indices are both representable in the uint32_t fields of Chain.
// Stage 0 is exactly !sizes_fit_u32(...); it computes nothing else.
//
// `<=`, not `<`: a buffer of exactly max() elements has largest index
// max() - 1, so no valid index can ever collide with the kNoVertex sentinel,
// and a `<` spelling would reject a buffer that fits.
//
// This is a free function rather than an expression buried in the member
// template because VertexCountOverflow is otherwise the one enumerator with no
// honest test -- firing it through build() needs 64 GiB of Point2 or 16 GiB of
// uint32_t. Being constexpr, allocation-free and taking plain sizes, it lets
// the suite pin all four boundary corners as static_asserts at no runtime cost.
// It is detail:: because it is not part of the type's contract; it is the
// testable half of one stage.
[[nodiscard]] constexpr bool sizes_fit_u32(std::size_t vertex_count,
                                           std::size_t index_count) noexcept {
    constexpr std::size_t limit = static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max());
    return vertex_count <= limit && index_count <= limit;
}

// The validator, in exactly the order the design specifies. THE ORDER IS NOT
// STYLISTIC: each stage establishes the precondition of the next, and the last
// two evaluate geometry that is meaningless if an earlier one failed. A chain
// that fails a stage is EXCLUDED from the later stages but does not stop them
// running on other chains.
//
//   0. overflow (early-returns)   3. finiteness       6. domain
//   1. structure                  4. stored closure
//   2. index range                5. winding
//
// The load-bearing part of stage 3 is not merely that it precedes stage 5 but
// that a chain carrying a non-finite vertex is EXCLUDED from stage 5. Those are
// two requirements and the second is the one that actually bites: nobody moves
// a stage, but everybody forgets to filter. orient2d on a NaN coordinate
// returns Collinear by orientation.hpp's documented total behaviour, and
// orientation<K> advances prev and next INDEPENDENTLY, so on a ring of four or
// more the walk routes around a single NaN vertex and returns the correct
// winding -- the faulty implementation then passes silently. On a triangle it
// reports DegenerateRing or WrongWinding, either of which sends the reader
// looking for geometry that is not there.
template <pred::GeometryKernel K>
[[nodiscard]] std::vector<PslgDiagnostic> validate_pslg(
        std::span<const Point2> vertices, std::span<const std::uint32_t> chain_indices,
        std::span<const Chain> chains) {
    std::vector<PslgDiagnostic> out;

    // 0. Overflow. The one stage that early-returns, because begin/count have
    // already been truncated and there is no more information to collect.
    if (!sizes_fit_u32(vertices.size(), chain_indices.size())) {
        out.push_back(PslgDiagnostic{
            PslgError::VertexCountOverflow, kNoChain, kNoVertex,
            std::format("{} vertices and {} chain indices do not both fit the 32-bit "
                        "fields of Chain",
                        vertices.size(), chain_indices.size())});
        return out;
    }

    const std::size_t n = chains.size();
    // Cleared by stages 1, 2 and 4: the chain can no longer be dereferenced or
    // wound meaningfully. Kept separate from `finite` so that stage 3's
    // exclusion applies to stage 5 alone, as specified.
    std::vector<char> usable(n, 1);
    std::vector<char> finite(n, 1);

    const auto span_of = [chain_indices](const Chain& ch) {
        return chain_indices.subspan(ch.begin, ch.count);
    };

    // 1. Structure, per chain, O(1). A one-vertex breakline is a point
    // constraint and this pipeline has no such thing. begin and count need no
    // bounds check against the flat buffer: the builder wrote them.
    for (std::size_t c = 0; c < n; ++c) {
        const Chain& ch = chains[c];
        const std::uint32_t minimum = is_closed(ch.role) ? 3u : 2u;
        if (ch.count < minimum) {
            usable[c] = 0;
            out.push_back(PslgDiagnostic{
                PslgError::ChainTooShort, static_cast<std::uint32_t>(c), kNoVertex,
                std::format("chain {} is declared {} and needs at least {} vertices, but has {}",
                            c, chain_role_name(ch.role), minimum, ch.count)});
        }
    }

    // 2. Index range, once over the whole flat buffer. THIS IS THE CHECK THAT
    // MAKES IndexedRing::vertex(i)'s UNCHECKED READ SAFE, and increment 2
    // deferred it here explicitly. One diagnostic per occurrence, not per
    // chain; the offending value and its position go in the message, because
    // `vertex` is documented as a valid index and this value is not one.
    for (std::size_t c = 0; c < n; ++c) {
        const std::span<const std::uint32_t> idx = span_of(chains[c]);
        for (std::size_t k = 0; k < idx.size(); ++k) {
            if (idx[k] >= vertices.size()) {
                usable[c] = 0;
                out.push_back(PslgDiagnostic{
                    PslgError::IndexOutOfRange, static_cast<std::uint32_t>(c), kNoVertex,
                    std::format("chain {} names vertex {} at position {}, but the buffer holds "
                                "only {} vertices",
                                c, idx[k], k, vertices.size())});
            }
        }
    }

    // 3. Finiteness, once over the WHOLE VERTEX BUFFER -- not per ring. A
    // vertex shared by k chains would otherwise be scanned k times and,
    // decisively, UNREFERENCED VERTICES MUST BE CHECKED TOO: the CDT wrapper
    // hands detria's setPoints the entire point array, so an unreferenced NaN
    // reaches the backend regardless of whether any chain names it.
    bool any_non_finite = false;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (!is_finite(vertices[i])) {
            any_non_finite = true;
            out.push_back(PslgDiagnostic{
                PslgError::NonFiniteVertex, kNoChain, static_cast<std::uint32_t>(i),
                std::format("vertex {} is not finite: {}", i, vertices[i])});
        }
    }
    if (any_non_finite) {
        for (std::size_t c = 0; c < n; ++c) {
            for (const std::uint32_t i : span_of(chains[c])) {
                if (i < vertices.size() && !is_finite(vertices[i])) {
                    finite[c] = 0;
                    break;
                }
            }
        }
    }

    // 4. Stored closure, per closed chain, O(1). Breaklines are EXEMPT: an open
    // chain with coincident endpoints is a closed polyline -- a contour, a ring
    // road that is not a domain boundary -- not an implied closure spelled
    // twice.
    //
    // The comparison is on POINTS, not indices: two distinct indices onto
    // coincident vertices close a ring just as surely as one index used twice.
    // `vertex` is idx[begin], the chain's FIRST index, ALWAYS and with no
    // branch -- in the one-index form the choice is vacuous, and in the
    // two-index form the first index is the one the caller keeps, since the fix
    // is to drop the trailing index and never the leading one. The message
    // carries both, so nothing is lost.
    for (std::size_t c = 0; c < n; ++c) {
        const Chain& ch = chains[c];
        if (usable[c] == 0 || !is_closed(ch.role)) {
            continue;
        }
        const std::span<const std::uint32_t> idx = span_of(ch);
        if (vertices[idx.front()] == vertices[idx.back()]) {
            usable[c] = 0;
            out.push_back(PslgDiagnostic{
                PslgError::StoredClosure, static_cast<std::uint32_t>(c), idx.front(),
                std::format("chain {} stores its closure: indices {} and {} name the same "
                            "point {}; drop the trailing index",
                            c, idx.front(), idx.back(), vertices[idx.front()])});
        }
    }

    // 5. Winding, per closed chain, one orientation<K> call -- never
    // signed_area, whose sign cancels to noise on a sliver at UTM33 magnitudes
    // and which increment 2 forbade from driving a topology decision. A wrong
    // winding is NEVER SILENTLY REVERSED: if the caller's geometry disagrees
    // with the caller's declared role, one of the two is a bug and guessing
    // which hides it forever. Python normalises with shapely, where reversal is
    // cheap and testable; C++ asserts what it was promised.
    //
    // Every chain reaching here passed stages 1, 2 and 4 and carries only
    // finite coordinates, which is exactly why IndexedRing's throwing
    // constructor is unreachable from this line.
    for (std::size_t c = 0; c < n; ++c) {
        const Chain& ch = chains[c];
        if (usable[c] == 0 || finite[c] == 0 || !is_closed(ch.role)) {
            continue;
        }
        const pred::Orientation observed = orientation<K>(IndexedRing{vertices, span_of(ch)});
        if (observed == pred::Orientation::Collinear) {
            out.push_back(PslgDiagnostic{
                PslgError::DegenerateRing, static_cast<std::uint32_t>(c), kNoVertex,
                std::format("chain {} is declared {} but every vertex is collinear",
                            c, chain_role_name(ch.role))});
            continue;
        }
        const pred::Orientation expected = ch.role == ChainRole::Outer
                                               ? pred::Orientation::CounterClockwise
                                               : pred::Orientation::Clockwise;
        if (observed != expected) {
            out.push_back(PslgDiagnostic{
                PslgError::WrongWinding, static_cast<std::uint32_t>(c), kNoVertex,
                std::format("chain {} is declared {}, which must be {}, but it is {}",
                            c, chain_role_name(ch.role), orientation_name(expected),
                            orientation_name(observed))});
        }
    }

    // 6. Domain, O(chains). A constraint set with no outer boundary has no
    // bounded domain to mesh, and addOutline is not optional for us.
    const bool has_outer = std::any_of(chains.begin(), chains.end(), [](const Chain& ch) {
        return ch.role == ChainRole::Outer;
    });
    if (!has_outer) {
        out.push_back(PslgDiagnostic{
            PslgError::NoOuterChain, kNoChain, kNoVertex,
            std::format("none of the {} chains has role Outer, so there is no bounded domain",
                        n)});
    }
    return out;
}

}  // namespace detail

// Formats the error, the chain, the vertex and the message. Nothing else --
// no coordinates beyond what the message already carries, no ring previews, no
// suggested fixes. Unbounded rendering of a ring was excluded in increment 2
// for the same reason.
[[nodiscard]] inline std::string describe(const PslgDiagnostic& d) {
    std::string out{detail::pslg_error_name(d.error)};
    if (d.chain != kNoChain) {
        out += std::format(" [chain {}]", d.chain);
    }
    if (d.vertex != kNoVertex) {
        out += std::format(" [vertex {}]", d.vertex);
    }
    if (!d.message.empty()) {
        out += ": ";
        out += d.message;
    }
    return out;
}

class PslgBuilder {
public:
    PslgBuilder() = default;

    explicit PslgBuilder(std::vector<Point2> vertices) : vertices_{std::move(vertices)} {}

    // Returns the index of the first new vertex. Callers whose geometry shares
    // vertices between chains -- a hole touching its outer ring at a node, two
    // breaklines meeting at a junction -- use this once and then the
    // index-taking add_chain, which is the encoding this representation exists
    // for.
    std::uint32_t append_vertices(std::span<const Point2> pts) {
        const auto first = static_cast<std::uint32_t>(vertices_.size());
        vertices_.insert(vertices_.end(), pts.begin(), pts.end());
        return first;
    }

    // ADD_CHAIN NEVER FAILS. It is total: an empty span, a one-element span and
    // an out-of-range index are all recorded and rejected at build(). A builder
    // that can fail mid-accumulation needs a second failure channel, and two
    // channels is how half the failures end up unreported.
    //
    // THE STORED INDEX ORDER IS THE ORDER RECEIVED, BYTE FOR BYTE; the builder
    // never rotates a chain. That is what makes the rotation hazard inherited
    // from increment 2 harmless here: {A, A, B, C} is a legal ring and its
    // rotation {A, B, C, A} is an illegal stored closure. Any future normaliser
    // -- one that rotates a chain to start at its lowest vertex, say -- MUST
    // COLLAPSE AN ADJACENT DUPLICATE FIRST.
    //
    // This index-taking overload is the primitive.
    PslgBuilder& add_chain(std::span<const std::uint32_t> idx, ChainRole role,
                           bool is_river = false) {
        const auto begin = static_cast<std::uint32_t>(chain_indices_.size());
        chain_indices_.insert(chain_indices_.end(), idx.begin(), idx.end());
        chains_.push_back(Chain{begin, static_cast<std::uint32_t>(idx.size()), role, is_river});
        return *this;
    }

    // Appends its points to the vertex buffer VERBATIM, WITH NO DEDUP, and
    // emits the consecutive index run. It exists so tests and simple Python
    // ingress do not each hand-roll the same loop.
    PslgBuilder& add_chain(std::span<const Point2> pts, ChainRole role, bool is_river = false) {
        const std::uint32_t first = append_vertices(pts);
        const auto count = static_cast<std::uint32_t>(pts.size());
        const auto begin = static_cast<std::uint32_t>(chain_indices_.size());
        for (std::uint32_t i = 0; i < count; ++i) {
            chain_indices_.push_back(first + i);
        }
        chains_.push_back(Chain{begin, count, role, is_river});
        return *this;
    }

    // RVALUE-REF-QUALIFIED: std::move(b).build<DefaultKernel>(). The builder is
    // consumed and its three buffers are moved into the Pslg, so a successful
    // build costs zero copies of the vertex data -- and "build twice and get
    // two Pslgs sharing nothing" is a compile error rather than a subtle
    // question about what the second build sees.
    //
    // K is the sole template parameter and is named at the call site. No
    // default, for the same reason ring.hpp gives none: defaulting it would
    // make this header include default_kernel.hpp and drag the compiled
    // terrain_predicates target into every consumer of a pure-header core type.
    template <pred::GeometryKernel K>
    [[nodiscard]] PslgBuildResult build() && {
        PslgBuildResult result;
        result.diagnostics = detail::validate_pslg<K>(vertices_, chain_indices_, chains_);
        if (result.diagnostics.empty()) {
            // Constructed as a local and moved in: std::optional::emplace would
            // construct from optional's context, which friendship does not
            // reach.
            Pslg p{std::move(vertices_), std::move(chain_indices_), std::move(chains_)};
            result.pslg = std::move(p);
        }
        return result;
    }

private:
    std::vector<Point2> vertices_;
    std::vector<std::uint32_t> chain_indices_;
    std::vector<Chain> chains_;
};

}  // namespace terrain
