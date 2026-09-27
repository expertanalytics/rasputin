#pragma once

// Lawson legalisation of a LatticeMesh (docs/increments/14b-delaunay-insertion.md,
// R1, R3 to R5).
//
// An edge is flipped only when the neighbour's apex is strictly inside the
// triangle's circle (`Incircle::Inside`), never on a cocircular tie, and never
// when it is constrained or has no neighbour. Each such flip strictly lowers
// the lifted triangulation, so the loops end; that needs the sign to be exact
// for fixed points, which `FilteredKernel<DetriaExact>` gives.
//
// The circle test runs in a LatticeFrame, (col * dx, -(row * dy)) on the
// fractional coordinates (docs/increments/16-domain-polygon.md, R2): world
// coordinates without the translation, so the Delaunay property holds in the
// world when dx != dy and the coordinates stay small for the filter.
//
// Both legalisers use an explicit stack popped last-in first-out, so the
// result depends only on the mesh and the seeds. They return the flip count
// and report every slot a flip writes through on_write (repeats allowed).
//
// Depends on core, predicates and lattice_mesh.hpp; knows no raster.

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/predicates/kernel.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace terrain::mesh {

// A frame built directly, LatticeFrame{dx, dy}, never enables the integer
// incircle; only lattice_frame() below can (docs/increments/21-parallel-refine.md,
// "Pinned by the red suite (21b)").
class LatticeFrame {
public:
    double dx = 1.0;
    double dy = 1.0;

    constexpr LatticeFrame() noexcept = default;
    constexpr LatticeFrame(double x, double y) noexcept : dx{x}, dy{y} {}

    [[nodiscard]] Point2 at(MeshVertex v) const noexcept {
        return Point2{v.col * dx, -(v.row * dy)};
    }
    // Whether lattice_incircle may answer on this frame.
    [[nodiscard]] bool integer() const noexcept { return integer_; }

private:
    bool integer_ = false;
    friend LatticeFrame lattice_frame(double, double, std::size_t, std::size_t) noexcept;
};

// The frame of a rows x cols grid, with the integer path enabled when the frame
// is the lattice times one exact constant (QW2): dx == dy, dx a positive normal
// double, and col * dx, row * dy exact for every node. Exactness is decided by
// QW2's sufficient condition, which is also what it takes: the significant bits
// of dx plus bit_width(max(rows, cols) - 1) are at most 53, so every product
// fits the mantissa; the largest one must also be finite.
[[nodiscard]] inline LatticeFrame lattice_frame(double dx, double dy, std::size_t rows,
                                                std::size_t cols) noexcept {
    LatticeFrame f{dx, dy};
    if (dx != dy || !std::isnormal(dx) || dx < 0.0)
        return f;
    const std::uint64_t mantissa =
        (std::bit_cast<std::uint64_t>(dx) & ((std::uint64_t{1} << 52) - 1)) | (std::uint64_t{1} << 52);
    const auto significant = 53 - std::countr_zero(mantissa);
    const std::size_t top = std::max(rows, cols);
    const auto span = top == 0 ? 0 : std::bit_width(top - 1);
    f.integer_ = significant + span <= 53
              && std::isfinite(static_cast<double>(top == 0 ? 0 : top - 1) * dx);
    return f;
}

// The incircle sign of the quad (a, b, c, d) from an int64 determinant on the
// lattice (col, -row), translated to d. Answers only on an enabling frame, four
// node corners and every difference from d at most 2^14 nodes; the determinant
// is then below 3 * 2^58 and exact. Under those conditions it is the sign
// DetriaExact gives on the frame points, since the frame is the lattice scaled
// by dx > 0 without rounding. Precondition: a, b, c counter-clockwise on
// (col, -row). Otherwise empty, and the caller asks the kernel.
[[nodiscard]] inline std::optional<pred::Incircle> lattice_incircle(MeshVertex a, MeshVertex b,
                                                                    MeshVertex c, MeshVertex d,
                                                                    const LatticeFrame& f) noexcept {
    if (!f.integer() || !a.is_node() || !b.is_node() || !c.is_node() || !d.is_node())
        return std::nullopt;
    constexpr double bound = 1 << 14;
    std::array<std::array<std::int64_t, 2>, 3> p{};  // (col, -row) minus d's
    const std::array<MeshVertex, 3> abc{a, b, c};
    for (std::size_t i = 0; i < 3; ++i) {
        const double x = abc[i].col - d.col, y = d.row - abc[i].row;  // exact: integers < 2^53
        if (std::fabs(x) > bound || std::fabs(y) > bound)
            return std::nullopt;
        p[i] = {static_cast<std::int64_t>(x), static_cast<std::int64_t>(y)};
    }
    const auto lift = [](const std::array<std::int64_t, 2>& q) { return q[0] * q[0] + q[1] * q[1]; };
    const auto cross = [](const std::array<std::int64_t, 2>& u, const std::array<std::int64_t, 2>& v) {
        return u[0] * v[1] - u[1] * v[0];
    };
    const std::int64_t det =
        lift(p[0]) * cross(p[1], p[2]) + lift(p[1]) * cross(p[2], p[0]) + lift(p[2]) * cross(p[0], p[1]);
    return det > 0 ? pred::Incircle::Inside : det < 0 ? pred::Incircle::Outside : pred::Incircle::Cocircular;
}

namespace detail {

// Whether t's edge e must flip: interior, unconstrained, and the quad's
// incircle determinant positive in the frame -- the apex across e inside t's
// circle, or t's apex inside the neighbour's circle.
template <pred::GeometryKernel K>
[[nodiscard]] bool must_flip(const LatticeMesh& m, std::uint32_t t, unsigned e,
                             const LatticeFrame& f) {
    const auto u = m.neighbours(t)[e];
    if (u == kNoNeighbour || m.is_constrained(t, e))
        return false;
    const auto& tri = m.triangles()[t];
    unsigned j = 0;
    while (m.triangles()[u][j] != tri[(e + 1) % 3])
        ++j;
    const auto v = m.vertices();
    const MeshVertex d_vertex = v[m.triangles()[u][(j + 2) % 3]];
    // The integer path first; the mesh triangle is counter-clockwise on
    // (col, -row), and it answers only where the kernel would give the same sign.
    if (const auto s = lattice_incircle(v[tri[e]], v[tri[(e + 1) % 3]], v[tri[(e + 2) % 3]], d_vertex, f))
        return *s == pred::Incircle::Inside;
    const Point2 a = f.at(v[tri[e]]), b = f.at(v[tri[(e + 1) % 3]]), c = f.at(v[tri[(e + 2) % 3]]),
                 d = f.at(d_vertex);
    // Orientation is exact on (col, -row), but the frame rounds col * dx and
    // row * dy, so a triangle counter-clockwise in the mesh can be collinear
    // or clockwise here. The kernel answers Cocircular for a collinear triple
    // and reorders a clockwise one, whose answer then has the wrong sign for
    // termination. So the test is asked only of a side counter-clockwise in
    // the frame. Both branches flip exactly when the same exact determinant,
    // on the frame points of the quad in mesh CCW order, is positive; a flip
    // lowers the sum of signed lifted volumes, so the process terminates.
    if (K::orient2d(a, b, c) == pred::Orientation::CounterClockwise)
        return K::incircle(a, b, c, d) == pred::Incircle::Inside;
    return K::orient2d(b, a, d) == pred::Orientation::CounterClockwise
        && K::incircle(b, a, d, c) == pred::Incircle::Inside;
}

}  // namespace detail

// Legalise around the vertex q after it was inserted. For each seed slot that
// contains q, the edge opposite q is tested; a flip leaves q at index 0 of both
// slots it writes, and both go back on the stack.
//
// The stack is the caller's, so a loop of insertions reuses one buffer
// (docs/increments/21-parallel-refine.md, 21a); its contents on entry are
// discarded, and it is empty on return.
using FlipStack = std::vector<std::uint32_t>;

template <pred::GeometryKernel K, class OnWrite>
std::size_t legalise_around(LatticeMesh& m, std::uint32_t q, std::span<const std::uint32_t> seeds,
                            const LatticeFrame& f, FlipStack& stack, OnWrite&& on_write) {
    stack.assign(seeds.begin(), seeds.end());
    std::size_t flips = 0;
    while (!stack.empty()) {
        const auto t = stack.back();
        stack.pop_back();
        const auto& tri = m.triangles()[t];
        unsigned i = 0;
        while (i < 3 && tri[i] != q)
            ++i;
        if (i == 3 || !detail::must_flip<K>(m, t, (i + 1) % 3, f))
            continue;
        const auto u = m.neighbours(t)[(i + 1) % 3];
        m.flip(t, (i + 1) % 3);
        ++flips;
        on_write(t);
        on_write(u);
        stack.push_back(t);
        stack.push_back(u);
    }
    return flips;
}

// The same, with a stack of its own.
template <pred::GeometryKernel K, class OnWrite>
std::size_t legalise_around(LatticeMesh& m, std::uint32_t q, std::span<const std::uint32_t> seeds,
                            const LatticeFrame& f, OnWrite&& on_write) {
    FlipStack stack;
    return legalise_around<K>(m, q, seeds, f, stack, std::forward<OnWrite>(on_write));
}

// Legalise the whole mesh: every interior edge once on the stack, and each
// flip pushes the four outer sides of its quad. A stale entry (its slot was
// rewritten since) names some current edge and is simply tested.
template <pred::GeometryKernel K, class OnWrite>
std::size_t legalise_all(LatticeMesh& m, const LatticeFrame& f, OnWrite&& on_write) {
    std::vector<std::pair<std::uint32_t, unsigned>> stack;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            if (const auto u = m.neighbours(t)[k]; u != kNoNeighbour && u > t)
                stack.emplace_back(t, k);
    std::size_t flips = 0;
    while (!stack.empty()) {
        const auto [t, e] = stack.back();
        stack.pop_back();
        if (!detail::must_flip<K>(m, t, e, f))
            continue;
        const auto u = m.neighbours(t)[e];
        m.flip(t, e);  // t = (c, a, d), u = (c, d, b)
        ++flips;
        on_write(t);
        on_write(u);
        stack.insert(stack.end(), {{t, 0u}, {t, 1u}, {u, 1u}, {u, 2u}});
    }
    return flips;
}

}  // namespace terrain::mesh
