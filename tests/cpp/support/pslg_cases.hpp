#pragma once

// Test-only scaffolding for the PSLG suites: diagnostic queries, the named
// constraint sets the unit suites reason about, and the generators the property
// suite draws from.
//
// The diagnostic helpers exist because every assertion in the builder suite has
// the same shape -- "this build failed, and among the diagnostics there is one
// with this error naming this chain" -- and spelling that out per test with a
// std::find_if buries the assertion in boilerplate. They deliberately answer
// questions about `error`, `chain` and `vertex` only: the `message` field is
// std::formatted prose that carries the offending values, and a suite that
// pinned its wording would fail on improvements to its own error messages.
//
// Nothing here asserts. A helper that assertions live inside reports failures
// at the helper's line, not the test's.
//
// Every generator takes an explicit std::mt19937_64, matching point_families
// and ring_cases: a test seeds once and the whole sequence is reproducible.

#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/core/ring.hpp>

#include <ring_cases.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <string>
#include <vector>

namespace terrain::test {

// ---------------------------------------------------------------------------
// Spans, without the std::span<const T>{v} ceremony at every call site
// ---------------------------------------------------------------------------

[[nodiscard]] inline std::span<const Point2> points(const std::vector<Point2>& v) {
    return std::span<const Point2>{v};
}

[[nodiscard]] inline std::span<const std::uint32_t> indices(const std::vector<std::uint32_t>& v) {
    return std::span<const std::uint32_t>{v};
}

// ---------------------------------------------------------------------------
// Diagnostic queries
// ---------------------------------------------------------------------------

[[nodiscard]] inline std::size_t count_errors(const PslgBuildResult& r, PslgError e) {
    return static_cast<std::size_t>(std::count_if(
        r.diagnostics.begin(), r.diagnostics.end(),
        [e](const PslgDiagnostic& d) { return d.error == e; }));
}

[[nodiscard]] inline const PslgDiagnostic* find_error(const PslgBuildResult& r, PslgError e) {
    const auto it = std::find_if(r.diagnostics.begin(), r.diagnostics.end(),
                                 [e](const PslgDiagnostic& d) { return d.error == e; });
    return it == r.diagnostics.end() ? nullptr : &*it;
}

// The error reported against one particular chain, which is what a test with
// several deliberately broken chains needs to ask.
[[nodiscard]] inline const PslgDiagnostic* find_error(const PslgBuildResult& r, PslgError e,
                                                      std::uint32_t chain) {
    const auto it = std::find_if(
        r.diagnostics.begin(), r.diagnostics.end(),
        [e, chain](const PslgDiagnostic& d) { return d.error == e && d.chain == chain; });
    return it == r.diagnostics.end() ? nullptr : &*it;
}

[[nodiscard]] inline bool has_error(const PslgBuildResult& r, PslgError e) {
    return find_error(r, e) != nullptr;
}

[[nodiscard]] inline bool has_error(const PslgBuildResult& r, PslgError e, std::uint32_t chain) {
    return find_error(r, e, chain) != nullptr;
}

// For failure messages: the whole diagnostic list rendered through the
// production formatter, so a surprising result prints what it actually was.
[[nodiscard]] inline std::string render(const PslgBuildResult& r) {
    std::string out;
    for (const PslgDiagnostic& d : r.diagnostics) {
        out += terrain::describe(d);
        out += '\n';
    }
    return out.empty() ? std::string{"<no diagnostics>"} : out;
}

// ---------------------------------------------------------------------------
// Named constraint sets
// ---------------------------------------------------------------------------

// Counterclockwise, no stored closure, exactly representable, and large enough
// to contain every hole and breakline below. Every kernel agrees on its winding
// -- the coordinates are integers under 2^53, so both the subtractions and the
// products inside orient2d are exact.
[[nodiscard]] inline std::vector<Point2> ccw_square(double s = 100.0) {
    return {Point2{0.0, 0.0}, Point2{s, 0.0}, Point2{s, s}, Point2{0.0, s}};
}

[[nodiscard]] inline std::vector<Point2> cw_square(double s = 100.0) {
    return {Point2{0.0, 0.0}, Point2{0.0, s}, Point2{s, s}, Point2{s, 0.0}};
}

// A clockwise square well inside ccw_square(100): a legal hole.
[[nodiscard]] inline std::vector<Point2> cw_hole() {
    return {Point2{20.0, 20.0}, Point2{20.0, 40.0}, Point2{40.0, 40.0}, Point2{40.0, 20.0}};
}

// The same square wound the other way: a hole the validator must reject.
[[nodiscard]] inline std::vector<Point2> ccw_hole() {
    return flipped(points(cw_hole()));
}

// Four vertices on one line. Not a ring under any winding rule, and the only
// input for which orientation<K> returns Collinear on finite coordinates.
[[nodiscard]] inline std::vector<Point2> collinear_spine() {
    return {Point2{1.0, 1.0}, Point2{2.0, 2.0}, Point2{3.0, 3.0}, Point2{4.0, 4.0}};
}

// An open polyline: no winding contract, no implied closure.
[[nodiscard]] inline std::vector<Point2> open_breakline() {
    return {Point2{5.0, 5.0}, Point2{50.0, 12.0}, Point2{80.0, 70.0}};
}

// A polyline whose first and last points coincide -- a contour, or a ring road
// that is not a domain boundary. Legitimate geometry as a Breakline, and a
// stored closure as anything else.
[[nodiscard]] inline std::vector<Point2> closed_polyline() {
    return {Point2{5.0, 5.0}, Point2{50.0, 12.0}, Point2{80.0, 70.0}, Point2{5.0, 5.0}};
}

// ---------------------------------------------------------------------------
// The UTM33 sliver the two kernels disagree about
// ---------------------------------------------------------------------------
//
// Three vertices at real UTM33 easting/northing magnitudes forming a hole
// 541 km long and under 4 metres tall, whose exact orientation is CLOCKWISE --
// so DefaultKernel accepts it as a Hole -- while FastKernel reports
// CounterClockwise and the same builder rejects it WrongWinding.
//
// The construction, and why it is not the increment-2 pattern:
//
//   * The apex P lies in the binade [2^18, 2^19) and carries full 2^-34
//     precision; the two far vertices lie in [2^19, 2^20) and carry full 2^-33
//     precision. The x differences P->v and P->N are therefore multiples of
//     2^-34 at a magnitude whose representable grid is 2^-33, so BOTH
//     SUBTRACTIONS ROUND, each by exactly half an ulp.
//   * All three northings share the binade [2^22, 2^23), so the y differences
//     are exact by Sterbenz and contribute no error of their own.
//
// Siting the error in the subtraction is the whole trick, and it is the lesson
// increment 2 paid for. If the edge vectors come out exact -- which they do for
// any ring whose coordinates share a binade -- then rounding-to-nearest is
// MONOTONE, so a contraction-free naive determinant can lose the sign (return
// 0.0) but can never invert it: FastKernel would be right, or Collinear, and
// this case would prove nothing on the CI leg that emits mulsd/mulsd/subsd.
// Three tests in increment 2 had exactly that defect.
//
// The constants were searched with exact rational arithmetic and accepted only
// if the perturbed determinant comes out positive under ALL THREE forms a
// compiler may emit for `a.x*b.y - a.y*b.x`:
//
//     fl(fl(a.x*b.y) - fl(a.y*b.x))     -- contraction off
//     fma(a.x, b.y, -fl(a.y*b.x))       -- one fma ordering
//     fma(-a.y, b.x,  fl(a.x*b.y))      -- the other
//
// and confirmed against the real headers under -ffp-contract={off,on,fast} and
// the compiler default. Two of those forms are invisible on any one host, which
// is why they are named here: anyone re-tuning this constant has to clear all
// three.
//
//     exact cross of the true coordinates = -236465675 / 2^64 (clockwise)
//
// Vertex order matters as much as the values. orientation<K> evaluates the
// triple at the LEXICOGRAPHICALLY LOWEST vertex, which is element 1 here, so
// the kernel call is orient2d(v[0], v[1], v[2]) -- the triple the search
// modelled. Rotating this ring would move the call and void the search.
[[nodiscard]] inline std::vector<Point2> utm33_sliver_hole() {
    return {
        Point2{331907.91178657551, 7851020.3979149926},
        Point2{873417.72277593287, 7851019.1743138293},
        Point2{679340.50859126227, 7851019.6128527075},
    };
}

// An outer ring large enough to contain the sliver above, with integer
// coordinates so that no kernel can disagree about it. It exists so the sliver
// case fails for exactly one reason.
[[nodiscard]] inline std::vector<Point2> utm33_enclosing_outer() {
    return {
        Point2{200000.0, 7800000.0},
        Point2{1000000.0, 7800000.0},
        Point2{1000000.0, 7900000.0},
        Point2{200000.0, 7900000.0},
    };
}

// ---------------------------------------------------------------------------
// Generators
// ---------------------------------------------------------------------------

// A chain description in the form the property suite manipulates before handing
// it to a builder: points plus the role they are declared under. Keeping the
// points separate from the builder is what lets a test corrupt exactly one
// chain and rebuild.
struct ChainSpec {
    std::vector<Point2> pts;
    ChainRole role{ChainRole::Breakline};
    bool is_river{false};
};

// A valid constraint set: one counterclockwise outer star, `holes` clockwise
// stars at disjoint centres, and `breaklines` open polylines. Star-shaped rings
// are simple by construction, so the winding a generated ring is declared under
// is the winding it has -- no rejection loop and no chance of a generated
// self-intersection smuggling a wrong expectation into a property.
[[nodiscard]] inline std::vector<ChainSpec> valid_chain_specs(std::mt19937_64& rng,
                                                              std::size_t holes,
                                                              std::size_t breaklines) {
    std::uniform_int_distribution<std::size_t> n_verts{3, 9};
    std::uniform_real_distribution<double> coord{-40.0, 40.0};

    std::vector<ChainSpec> specs;
    specs.push_back(ChainSpec{star_ring(rng, 12, Point2{0.0, 0.0}, 80.0, 100.0),
                              ChainRole::Outer, false});

    for (std::size_t h = 0; h < holes; ++h) {
        const Point2 centre{coord(rng), coord(rng)};
        std::vector<Point2> ccw = star_ring(rng, n_verts(rng), centre, 1.0, 3.0);
        specs.push_back(ChainSpec{flipped(points(ccw)), ChainRole::Hole, false});
    }
    for (std::size_t b = 0; b < breaklines; ++b) {
        std::vector<Point2> line;
        const std::size_t n = n_verts(rng);
        for (std::size_t i = 0; i < n; ++i) line.push_back(Point2{coord(rng), coord(rng)});
        specs.push_back(ChainSpec{std::move(line), ChainRole::Breakline, b % 2 == 0});
    }
    return specs;
}

// Feeds specs to a builder through the point-taking overload, which appends
// verbatim: chain i owns vertices [sum of earlier counts, +count).
[[nodiscard]] inline PslgBuilder builder_from(const std::vector<ChainSpec>& specs) {
    PslgBuilder b;
    for (const ChainSpec& s : specs) b.add_chain(points(s.pts), s.role, s.is_river);
    return b;
}

// The concatenation the builder above must produce, element for element. The
// "no dedup, no reordering, no reversal" guarantee is the assertion this
// enables.
[[nodiscard]] inline std::vector<Point2> concatenated(const std::vector<ChainSpec>& specs) {
    std::vector<Point2> out;
    for (const ChainSpec& s : specs) out.insert(out.end(), s.pts.begin(), s.pts.end());
    return out;
}

}  // namespace terrain::test
