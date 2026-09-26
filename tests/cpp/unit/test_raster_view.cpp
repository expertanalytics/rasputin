// Increment 12 (docs/increments/12-dem-to-mesh.md, R2 and R4): the non-owning
// RasterView, row access on the RasterSource concept, and the batch sampler.
//
// Tests 1-5 of the design's list. Tests 2 and 4 are invariant-critical: the
// view must agree with the owning Raster over the same buffer, and the batch
// sampler's edge semantics must be bilinear's, point for point.
//
// The batch function's name and signature are NOT in the design (R2 says only
// "a batch function in sample.hpp ... writes z and a valid flag per point").
// This suite assumes
//     bilinear_batch(const R&, std::span<const Point2> points,
//                    std::span<double> z, std::span<bool> valid)
// and the handback names that as a gap.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/raster/view.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <vector>

using Catch::Matchers::WithinRel;

using terrain::Point2;
using terrain::raster::bilinear;
using terrain::raster::bilinear_batch;
using terrain::raster::CellIndex;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::raster::RasterSource;
using terrain::raster::RasterView;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kSentinel = -32767.0;

// Non-square with unequal spacing, so a transposed index or swapped delta
// cannot pass. UTM-scale origin, like the real fixture.
RasterGeometry grid_4x5() {
    return RasterGeometry{/*x_min=*/500000.0, /*y_max=*/6600000.0, /*delta_x=*/10.0,
                          /*delta_y=*/5.0, /*cols=*/5, /*rows=*/4};
}

// Cell (r, c) holds 10 r + c: every value distinct.
template <typename T>
std::vector<T> distinct_values(const RasterGeometry& g) {
    std::vector<T> v(g.size());
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c)
            v[g.linear_index(CellIndex{r, c})] = static_cast<T>(10 * r + c);
    return v;
}

// A buffer with one sentinel cell and one NaN cell, apart from each other.
template <typename T>
std::vector<T> with_voids(const RasterGeometry& g) {
    auto v = distinct_values<T>(g);
    v[g.linear_index(CellIndex{0, 4})] = static_cast<T>(kSentinel);
    v[g.linear_index(CellIndex{3, 0})] = std::numeric_limits<T>::quiet_NaN();
    return v;
}

// A RasterSource model WITHOUT row(): after the concept gains row access it
// must no longer satisfy the concept.
class RowlessRaster {
public:
    using value_type = double;
    explicit RowlessRaster(RasterGeometry g) : geometry_{g} {}
    [[nodiscard]] const RasterGeometry& geometry() const noexcept { return geometry_; }
    [[nodiscard]] double value_at(const CellIndex&) const noexcept { return 0.0; }
    [[nodiscard]] bool is_nodata(const CellIndex&) const noexcept { return false; }
    // Present so the missing row() is the only reason it fails (increment 18, C2 (b)).
    [[nodiscard]] const std::optional<double>& nodata() const noexcept { return nodata_; }

private:
    RasterGeometry geometry_;
    std::optional<double> nodata_;
};

// The independent test double from test_raster.cpp, with row() added as the
// design requires (R4).
class ConstantRaster {
public:
    using value_type = double;
    ConstantRaster(RasterGeometry g, double v) : geometry_{g}, row_(g.cols(), v) {}
    [[nodiscard]] const RasterGeometry& geometry() const noexcept { return geometry_; }
    [[nodiscard]] double value_at(const CellIndex&) const noexcept { return row_.front(); }
    [[nodiscard]] bool is_nodata(const CellIndex&) const noexcept { return false; }
    [[nodiscard]] const std::optional<double>& nodata() const noexcept { return nodata_; }
    [[nodiscard]] std::span<const double> row(std::size_t) const noexcept { return row_; }

private:
    RasterGeometry geometry_;
    std::vector<double> row_;
    std::optional<double> nodata_;  // none; RasterSource requires nodata() (increment 18, C2 (b))
};

// Points that exercise every branch of bilinear: nodes, cell interiors, the
// last row and column, next to the sentinel and the NaN, outside, non-finite.
std::vector<Point2> probe_points(const RasterGeometry& g) {
    std::vector<Point2> pts;
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c) {
            pts.push_back(g.node(CellIndex{r, c}));
            pts.push_back(Point2{g.node(CellIndex{r, c}).x + 3.0,
                                 g.node(CellIndex{r, c}).y - 1.25});
        }
    pts.push_back(Point2{g.x_min() - 1.0, g.y_max()});
    pts.push_back(Point2{g.x_max() + 1.0, g.y_min()});
    pts.push_back(Point2{g.x_min(), g.y_max() + 1.0});
    pts.push_back(Point2{kNaN, g.y_max()});
    pts.push_back(Point2{g.x_min(), kInf});
    return pts;
}

} // namespace

// ---------------------------------------------------------------------------
// Test 1: concept satisfaction, and that row() is now required.
// ---------------------------------------------------------------------------

TEST_CASE("RasterView<float> and RasterView<double> satisfy RasterSource",
          "[raster][view][concept]") {
    STATIC_REQUIRE(RasterSource<RasterView<float>>);
    STATIC_REQUIRE(RasterSource<RasterView<double>>);
    STATIC_REQUIRE(RasterSource<Raster<float>>);
    STATIC_REQUIRE(RasterSource<Raster<double>>);
    STATIC_REQUIRE(RasterSource<ConstantRaster>);
}

TEST_CASE("RasterSource requires row access", "[raster][concept]") {
    STATIC_REQUIRE_FALSE(RasterSource<RowlessRaster>);
}

TEST_CASE("RasterView is a view: its value_type is the cell type",
          "[raster][view][concept]") {
    STATIC_REQUIRE(std::same_as<RasterView<float>::value_type, float>);
    STATIC_REQUIRE(std::same_as<RasterView<double>::value_type, double>);
}

// ---------------------------------------------------------------------------
// Test 2 (invariant-critical): the view agrees with the owning raster.
// ---------------------------------------------------------------------------

template <typename T>
void require_view_matches_raster() {
    const auto g = grid_4x5();
    const auto data = with_voids<T>(g);
    const std::optional<T> nodata{static_cast<T>(kSentinel)};
    const Raster<T> owning{g, data, nodata};
    const RasterView<T> view{g, data.data(), nodata};

    REQUIRE(view.geometry().cols() == g.cols());
    REQUIRE(view.geometry().rows() == g.rows());
    REQUIRE(view.geometry().x_min() == g.x_min());
    REQUIRE(view.geometry().y_max() == g.y_max());

    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c) {
            const CellIndex i{r, c};
            CAPTURE(r, c);
            REQUIRE(view.is_nodata(i) == owning.is_nodata(i));
            if (!owning.is_nodata(i))
                REQUIRE(view.value_at(i) == owning.value_at(i));
        }
    REQUIRE(view.is_nodata(CellIndex{0, 4}));  // the sentinel
    REQUIRE(view.is_nodata(CellIndex{3, 0}));  // NaN, with a sentinel also set
    REQUIRE_FALSE(view.is_nodata(CellIndex{1, 1}));

    for (const Point2& p : probe_points(g)) {
        CAPTURE(p.x, p.y);
        const auto a = bilinear(owning, p);
        const auto b = bilinear(view, p);
        REQUIRE(a.has_value() == b.has_value());
        if (a)
            REQUIRE(*a == *b);  // same arithmetic over the same bytes: bit-equal
    }
}

TEST_CASE("RasterView<float> agrees with Raster<float> over the same buffer",
          "[raster][view][invariant]") {
    require_view_matches_raster<float>();
}

TEST_CASE("RasterView<double> agrees with Raster<double> over the same buffer",
          "[raster][view][invariant]") {
    require_view_matches_raster<double>();
}

TEST_CASE("RasterView without a sentinel treats only NaN as NoData",
          "[raster][view][edge]") {
    const auto g = grid_4x5();
    const auto data = with_voids<double>(g);
    const RasterView<double> view{g, data.data(), std::nullopt};
    REQUIRE_FALSE(view.is_nodata(CellIndex{0, 4}));
    REQUIRE(view.value_at(CellIndex{0, 4}) == kSentinel);
    REQUIRE(view.is_nodata(CellIndex{3, 0}));
}

// ---------------------------------------------------------------------------
// Test 3: row access, on the view and on the owning raster.
// ---------------------------------------------------------------------------

template <typename Source>
void require_rows_match_cells(const Source& s) {
    const auto& g = s.geometry();
    for (const std::size_t r : {std::size_t{0}, g.rows() - 1}) {
        CAPTURE(r);
        const auto row = s.row(r);
        REQUIRE(row.size() == g.cols());
        REQUIRE(row.front() == s.value_at(CellIndex{r, 0}));
        REQUIRE(row.back() == s.value_at(CellIndex{r, g.cols() - 1}));
    }
}

TEST_CASE("row(i) has cols elements and starts at value_at({i, 0})", "[raster][row]") {
    const auto g = grid_4x5();
    const auto data = distinct_values<float>(g);
    const RasterView<float> view{g, data.data(), std::nullopt};
    const Raster<float> owning{g, data};

    STATIC_REQUIRE(std::same_as<decltype(view.row(0)), std::span<const float>>);
    STATIC_REQUIRE(std::same_as<decltype(owning.row(0)), std::span<const float>>);

    require_rows_match_cells(view);
    require_rows_match_cells(owning);
    require_rows_match_cells(ConstantRaster{g, 7.0});

    // The view's row points into the caller's buffer, not a copy.
    REQUIRE(view.row(2).data() == data.data() + 2 * g.cols());
}

// ---------------------------------------------------------------------------
// Test 4 (invariant-critical): the batch sampler's edge semantics.
// ---------------------------------------------------------------------------

namespace {

struct Batch {
    std::vector<double> z;
    std::vector<bool> valid;
};

template <typename R>
Batch run_batch(const R& raster, const std::vector<Point2>& pts) {
    std::vector<double> z(pts.size(), -1.0);
    // bool[] rather than vector<bool>: the latter has no contiguous storage.
    auto flags = std::make_unique<bool[]>(pts.size());
    bilinear_batch(raster, std::span<const Point2>{pts}, std::span<double>{z},
                   std::span<bool>{flags.get(), pts.size()});
    return Batch{std::move(z), std::vector<bool>(flags.get(), flags.get() + pts.size())};
}

} // namespace

TEST_CASE("bilinear_batch agrees with bilinear point for point", "[raster][batch][invariant]") {
    const auto g = grid_4x5();
    const auto data = with_voids<float>(g);
    const RasterView<float> view{g, data.data(), static_cast<float>(kSentinel)};
    const auto pts = probe_points(g);
    const Batch out = run_batch(view, pts);

    for (std::size_t i = 0; i < pts.size(); ++i) {
        CAPTURE(i, pts[i].x, pts[i].y);
        const auto expected = bilinear(view, pts[i]);
        REQUIRE(out.valid[i] == expected.has_value());
        if (expected)
            REQUIRE(out.z[i] == *expected);
    }
}

TEST_CASE("bilinear_batch: outside the grid and non-finite points are not valid",
          "[raster][batch][edge]") {
    const auto g = grid_4x5();
    const auto data = distinct_values<double>(g);
    const RasterView<double> view{g, data.data(), std::nullopt};
    const std::vector<Point2> pts{
        {g.x_min() - 1e-6, g.y_max()},  // just west
        {g.x_max() + 1e-6, g.y_min()},  // just east
        {g.x_min(), g.y_max() + 1e-6},  // just north
        {g.x_min(), g.y_min() - 1e-6},  // just south
        {kNaN, g.y_max()},
        {g.x_min(), kNaN},
        {kInf, g.y_max()},
        {g.x_min(), -kInf},
    };
    const Batch out = run_batch(view, pts);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        CAPTURE(i);
        REQUIRE_FALSE(out.valid[i]);
    }
}

TEST_CASE("bilinear_batch: a zero-weight NoData corner still makes the point invalid",
          "[raster][batch][edge]") {
    const auto g = grid_4x5();
    auto data = distinct_values<double>(g);
    data[g.linear_index(CellIndex{1, 2})] = kSentinel;
    const RasterView<double> view{g, data.data(), kSentinel};

    // Node (1, 1) is valid, but its bilinear cell (1, 1) has (1, 2) as a corner
    // with weight zero. Node (2, 3) is far from the void and stays valid.
    const std::vector<Point2> pts{g.node(CellIndex{1, 1}), g.node(CellIndex{2, 3})};
    const Batch out = run_batch(view, pts);
    REQUIRE_FALSE(out.valid[0]);
    REQUIRE(out.valid[1]);
    REQUIRE(out.z[1] == 23.0);
}

TEST_CASE("bilinear_batch: NaN corners are NoData even with no sentinel",
          "[raster][batch][edge]") {
    const auto g = grid_4x5();
    auto data = distinct_values<double>(g);
    // (0, 2) is a zero-weight corner of node (0, 1)'s bilinear cell (0, 1).
    data[g.linear_index(CellIndex{0, 2})] = kNaN;
    const RasterView<double> view{g, data.data(), std::nullopt};
    const std::vector<Point2> pts{g.node(CellIndex{0, 1}), g.node(CellIndex{2, 2})};
    const Batch out = run_batch(view, pts);
    REQUIRE_FALSE(out.valid[0]);
    REQUIRE(out.z[0] == 0.0);  // an invalid point's z is 0.0, never NaN (R2)
    REQUIRE(out.valid[1]);
    REQUIRE_FALSE(std::isnan(out.z[1]));
}

TEST_CASE("bilinear_batch: z is the node value, exactly inside, to 1e-9 on the far edges",
          "[raster][batch][invariant]") {
    const auto g = grid_4x5();
    const auto data = distinct_values<float>(g);
    const RasterView<float> view{g, data.data(), std::nullopt};

    std::vector<Point2> pts;
    std::vector<CellIndex> cells;
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c) {
            cells.push_back(CellIndex{r, c});
            pts.push_back(g.node(cells.back()));
        }
    const Batch out = run_batch(view, pts);

    for (std::size_t i = 0; i < pts.size(); ++i) {
        const CellIndex c = cells[i];
        CAPTURE(c.row, c.col);
        REQUIRE(out.valid[i]);
        const double node = static_cast<double>(view.value_at(c));
        const bool far_edge = c.row == g.rows() - 1 || c.col == g.cols() - 1;
        if (far_edge)
            REQUIRE_THAT(out.z[i], WithinRel(node, 1e-9));
        else
            REQUIRE(out.z[i] == node);
    }
}

TEST_CASE("bilinear_batch over an empty span writes nothing", "[raster][batch][edge]") {
    const auto g = grid_4x5();
    const auto data = distinct_values<double>(g);
    const RasterView<double> view{g, data.data(), std::nullopt};
    const Batch out = run_batch(view, {});
    REQUIRE(out.z.empty());
    REQUIRE(out.valid.empty());
}

TEST_CASE("bilinear_batch accepts any RasterSource, not only the view", "[raster][batch]") {
    const auto g = grid_4x5();
    const Batch out = run_batch(ConstantRaster{g, 7.0}, {g.node(CellIndex{1, 1})});
    REQUIRE(out.valid[0]);
    REQUIRE(out.z[0] == 7.0);
}

// ---------------------------------------------------------------------------
// Test 5: the view does not copy.
// ---------------------------------------------------------------------------

TEST_CASE("RasterView does not copy: a changed buffer value shows through",
          "[raster][view]") {
    const auto g = grid_4x5();
    auto data = distinct_values<double>(g);
    const RasterView<double> view{g, data.data(), kSentinel};

    data[g.linear_index(CellIndex{2, 3})] = 99.5;
    REQUIRE(view.value_at(CellIndex{2, 3}) == 99.5);

    data[g.linear_index(CellIndex{1, 1})] = kSentinel;
    REQUIRE(view.is_nodata(CellIndex{1, 1}));
    REQUIRE_FALSE(bilinear(view, g.node(CellIndex{1, 1})).has_value());
}
