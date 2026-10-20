#pragma once

// The exact overlap of triangles with the cells of a class grid.
//
// Everything is in the grid's index space: x is the column and y the row, so
// cell (r, c) is the unit square [c, c + 1] x [r, r + 1] and the window is
// [0, cols] x [0, rows], closed. Per triangle, per row it spans, the triangle
// is clipped to the strip r <= y <= r + 1 (Sutherland and Hodgman, two
// horizontal lines); in the strip, every column wholly between the strip's
// left and right chains is a full cell (area 1, no clip) and the columns its
// edges cross are clipped by vertical lines. Each piece's shoelace area is
// taken in coordinates local to its own first vertex and multiplied by the
// cell area of its row; pieces of zero area are skipped.
//
// Every crossing of a triangle edge with a grid line is computed from the
// edge's two endpoints taken in one canonical order (the lexicographically
// smaller first), never from an already clipped piece, so the two triangles on
// a shared edge, and the two cells on either side of a grid line, see the same
// crossing bit for bit. The triangle's corners are first rotated to start at
// the smallest vertex, counter-clockwise, so the vertex order of the input
// changes nothing. Triangles are independent, so blocks run in parallel and
// the rows are concatenated in triangle order: the result does not depend on
// the thread count, nor on how the triangles are cut into calls.
//
// Pure: no I/O, no Python. The view's nodata is not consulted (class 0 is a
// class like any other, 1.1).

#include <terrain/core/point.hpp>
#include <terrain/parallel_util/chunks.hpp>
#include <terrain/raster/view.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <span>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace terrain::land_cover {

// A sparse table, one row per triangle: triangle t's entries are
// [offsets[t], offsets[t + 1]), classes ascending within a row, areas in the
// unit of row_area (m^2) and positive.
struct ClassAreas {
    std::vector<std::uint32_t> offsets;
    std::vector<std::uint8_t> classes;
    std::vector<double> areas;
};

// A triangle with a corner past the closed window [0, cols] x [0, rows]. The
// window is built to hold every vertex plus one cell, so this is a defect of
// the caller, refused for the whole batch and never clamped.
struct OutsideWindow : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

namespace detail {

// A polygon vertex and the kind of the edge leaving it: 0..2 the index of the
// triangle edge it lies on, kHorizontal a cut along y = const, kVertical a cut
// along x = const.
constexpr int kHorizontal = -1;
constexpr int kVertical = -2;

struct Vertex {
    Point2 p;
    int kind;
};

// At most 7 vertices: a triangle gains one per half-plane clip, and there are
// four (two horizontal, two vertical).
struct Polygon {
    std::array<Vertex, 8> v{};
    std::size_t n = 0;
};

[[nodiscard]] constexpr bool lex_less(const Point2& a, const Point2& b) noexcept {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}

// The crossing of triangle edge `e` (corners `tri`) with the line y = Y
// (horizontal) or x = X (vertical), from the edge's canonical endpoints. An
// endpoint on the line is the crossing itself, so a corner on a grid line is
// never moved by rounding.
[[nodiscard]] inline Point2 canonical_crossing(const std::array<Point2, 3>& tri, int e,
                                               bool horizontal, double line) noexcept {
    const Point2& p = tri[static_cast<std::size_t>(e)];
    const Point2& q = tri[static_cast<std::size_t>((e + 1) % 3)];
    const Point2& a = lex_less(p, q) ? p : q;
    const Point2& b = lex_less(p, q) ? q : p;
    if (horizontal) {
        if (a.y == line) return a;
        if (b.y == line) return b;
        const double t = std::clamp((line - a.y) / (b.y - a.y), 0.0, 1.0);
        return Point2{a.x + t * (b.x - a.x), line};
    }
    if (a.x == line) return a;
    if (b.x == line) return b;
    const double t = std::clamp((line - a.x) / (b.x - a.x), 0.0, 1.0);
    return Point2{line, a.y + t * (b.y - a.y)};
}

// Sutherland and Hodgman against one half-plane: `horizontal` picks y against
// `line`, else x; `keep_below` keeps coordinates <= line, else >= line. A
// point on the line is inside.
[[nodiscard]] inline Polygon clip(const Polygon& in, const std::array<Point2, 3>& tri,
                                  bool horizontal, double line, bool keep_below) noexcept {
    const auto coordinate = [horizontal](const Point2& p) { return horizontal ? p.y : p.x; };
    const auto inside = [&](const Point2& p) {
        return keep_below ? coordinate(p) <= line : coordinate(p) >= line;
    };
    const auto crossing = [&](const Vertex& s) {
        if (s.kind >= 0)
            return canonical_crossing(tri, s.kind, horizontal, line);
        // A cut edge is axis-parallel: the crossing keeps its constant coordinate.
        return horizontal ? Point2{s.p.x, line} : Point2{line, s.p.y};
    };
    Polygon out;
    for (std::size_t i = 0; i < in.n; ++i) {
        const Vertex& s = in.v[i];
        const Vertex& e = in.v[(i + 1) % in.n];
        const bool s_in = inside(s.p), e_in = inside(e.p);
        if (s_in && e_in) {
            out.v[out.n++] = e;
        } else if (s_in) {
            // Leaving: the new vertex starts the cut edge along the clip line.
            out.v[out.n++] = Vertex{crossing(s), horizontal ? kHorizontal : kVertical};
        } else if (e_in) {
            // Entering: the new vertex continues the edge s lies on.
            out.v[out.n++] = Vertex{crossing(s), s.kind};
            out.v[out.n++] = e;
        }
    }
    return out;
}

// The shoelace area in coordinates local to the piece's first vertex, which
// lies in the cell: the differences are exact to about an ulp of 1 and a piece
// a nanometre across keeps its area. The sign is dropped so a clockwise
// triangle gives the same areas as its mirror image.
[[nodiscard]] inline double area_local(const Polygon& poly) noexcept {
    double twice = 0.0;
    const Point2 origin = poly.v[0].p;
    for (std::size_t i = 1; i + 1 < poly.n; ++i)
        twice += cross(poly.v[i].p - origin, poly.v[i + 1].p - origin);
    return 0.5 * std::fabs(twice);
}

// The per-block accumulator: one slot per class, and the classes touched by
// the current triangle, so a triangle costs its classes and not 256.
struct Accumulator {
    std::array<double, 256> area{};
    std::vector<std::uint8_t> touched;

    void add(std::uint8_t k, double a) {
        if (area[k] == 0.0)
            touched.push_back(k);
        area[k] += a;
    }
};

// One triangle's pieces into `acc`: rows, then strips, then columns.
inline void overlap_one(const std::array<Point2, 3>& tri,
                        const raster::RasterView<std::uint8_t>& classes,
                        std::span<const double> row_area, Accumulator& acc) {
    const double y_lo = std::min({tri[0].y, tri[1].y, tri[2].y});
    const double y_hi = std::max({tri[0].y, tri[1].y, tri[2].y});
    // Clamped to the grid: a crossing's rounding may put x_hi an ulp past cols.
    const std::size_t rows = classes.geometry().rows(), cols = classes.geometry().cols();
    const auto r_begin = static_cast<std::size_t>(std::max(0.0, std::floor(y_lo)));
    const auto r_end = std::min(rows, static_cast<std::size_t>(std::max(0.0, std::ceil(y_hi))));
    Polygon whole;
    for (std::size_t i = 0; i < 3; ++i)
        whole.v[whole.n++] = Vertex{tri[i], static_cast<int>(i)};
    for (std::size_t r = r_begin; r < r_end; ++r) {
        const auto fr = static_cast<double>(r);
        const Polygon strip = clip(clip(whole, tri, true, fr, false), tri, true, fr + 1.0, true);
        if (strip.n < 3)
            continue;
        double x_lo = std::numeric_limits<double>::infinity(), x_hi = -x_lo;
        double y_min = x_lo, y_max = x_hi;
        for (std::size_t i = 0; i < strip.n; ++i) {
            x_lo = std::min(x_lo, strip.v[i].p.x);
            x_hi = std::max(x_hi, strip.v[i].p.x);
            y_min = std::min(y_min, strip.v[i].p.y);
            y_max = std::max(y_max, strip.v[i].p.y);
        }
        // Full columns exist only when the strip spans the row's whole height.
        // The left chain's x is convex in y and the right chain's concave, so
        // their extremes are at the strip's bottom and top.
        std::size_t full_begin = 0, full_end = 0;
        if (y_min == fr && y_max == fr + 1.0) {
            double left_bottom = x_hi, left_top = x_hi, right_bottom = x_lo, right_top = x_lo;
            for (std::size_t i = 0; i < strip.n; ++i) {
                const Point2& p = strip.v[i].p;
                if (p.y == fr) {
                    left_bottom = std::min(left_bottom, p.x);
                    right_bottom = std::max(right_bottom, p.x);
                }
                if (p.y == fr + 1.0) {
                    left_top = std::min(left_top, p.x);
                    right_top = std::max(right_top, p.x);
                }
            }
            const double left = std::max(left_bottom, left_top);
            const double right = std::min(right_bottom, right_top);
            if (left < right) {
                full_begin = static_cast<std::size_t>(std::ceil(left));
                full_end = static_cast<std::size_t>(std::floor(right));
            }
        }
        const auto c_begin = static_cast<std::size_t>(std::max(0.0, std::floor(x_lo)));
        const auto c_end = std::min(cols, static_cast<std::size_t>(std::max(0.0, std::ceil(x_hi))));
        const std::span<const std::uint8_t> row = classes.row(r);
        const double weight = row_area[r];
        for (std::size_t c = c_begin; c < c_end; ++c) {
            double piece = 1.0;
            if (c < full_begin || c >= full_end) {
                const auto fc = static_cast<double>(c);
                const Polygon cell = clip(clip(strip, tri, false, fc, false), tri, false, fc + 1.0, true);
                piece = cell.n < 3 ? 0.0 : area_local(cell);
            }
            if (piece > 0.0)
                acc.add(row[c], piece * weight);
        }
    }
}

}  // namespace detail

// The class areas of every triangle (three corners each) over `classes`, with
// `row_area` the area of one cell of each row. Refusals, all before any work
// and all for the whole batch: a corners span whose length is not a multiple
// of three, a row_area span whose length is not rows(), a non-finite corner,
// a zero-area triangle (std::invalid_argument), and a corner past the closed
// window (OutsideWindow). threads 0 means every core; the result does not
// depend on it.
[[nodiscard]] inline ClassAreas class_areas(std::span<const Point2> corners,
                                            const raster::RasterView<std::uint8_t>& classes,
                                            std::span<const double> row_area,
                                            unsigned threads = 0) {
    const std::size_t rows = classes.geometry().rows(), cols = classes.geometry().cols();
    if (corners.size() % 3 != 0)
        throw std::invalid_argument("class_areas: corners must hold three points per triangle");
    if (row_area.size() != rows)
        throw std::invalid_argument("class_areas: row_area must hold one cell area per row");
    const std::size_t n = corners.size() / 3;
    for (std::size_t t = 0; t < n; ++t) {
        const Point2 a = corners[3 * t], b = corners[3 * t + 1], c = corners[3 * t + 2];
        for (const Point2& p : {a, b, c}) {
            if (!std::isfinite(p.x) || !std::isfinite(p.y))
                throw std::invalid_argument("class_areas: a triangle corner is not finite");
            if (p.x < 0.0 || p.y < 0.0 || p.x > static_cast<double>(cols) || p.y > static_cast<double>(rows))
                throw OutsideWindow("class_areas: a triangle corner lies outside the class window");
        }
        if (cross(b - a, c - a) == 0.0)
            throw std::invalid_argument("class_areas: a triangle has zero area");
    }

    std::vector<std::uint32_t> counts(n);
    struct Entry {
        std::uint8_t k;
        double area;
    };
    const parallel_util::BlockSchedule schedule{};
    const unsigned workers = threads == 0 ? std::max(1u, std::thread::hardware_concurrency()) : threads;
    const std::size_t block = parallel_util::default_block(n, workers);
    const std::size_t blocks = n / block + (n % block != 0 ? 1 : 0);
    std::vector<std::vector<Entry>> per_block(blocks);
    parallel_util::for_each_block(n, threads, {block, schedule.inline_below}, [&](std::size_t begin, std::size_t end) {
        detail::Accumulator acc;
        std::vector<Entry>& out = per_block[begin / block];
        for (std::size_t t = begin; t < end; ++t) {
            // The canonical rotation: start at the smallest corner, counter-clockwise.
            std::array<Point2, 3> tri{corners[3 * t], corners[3 * t + 1], corners[3 * t + 2]};
            if (cross(tri[1] - tri[0], tri[2] - tri[0]) < 0.0)
                std::swap(tri[1], tri[2]);
            while (detail::lex_less(tri[1], tri[0]) || detail::lex_less(tri[2], tri[0]))
                std::rotate(tri.begin(), tri.begin() + 1, tri.end());
            detail::overlap_one(tri, classes, row_area, acc);
            std::sort(acc.touched.begin(), acc.touched.end());
            counts[t] = static_cast<std::uint32_t>(acc.touched.size());
            for (const std::uint8_t k : acc.touched) {
                out.push_back(Entry{k, acc.area[k]});
                acc.area[k] = 0.0;
            }
            acc.touched.clear();
        }
    });

    ClassAreas result;
    result.offsets.resize(n + 1);
    std::uint64_t total = 0;
    for (std::size_t t = 0; t < n; ++t) {
        result.offsets[t] = static_cast<std::uint32_t>(total);
        total += counts[t];
        if (total > std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("class_areas: more than 2^32 - 1 entries");
    }
    result.offsets[n] = static_cast<std::uint32_t>(total);
    result.classes.reserve(total);
    result.areas.reserve(total);
    for (const std::vector<Entry>& entries : per_block)
        for (const Entry& e : entries) {
            result.classes.push_back(e.k);
            result.areas.push_back(e.area);
        }
    return result;
}

}  // namespace terrain::land_cover
