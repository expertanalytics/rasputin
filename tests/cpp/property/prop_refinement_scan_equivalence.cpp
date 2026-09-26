// Increment 18 (docs/increments/18-row-span-scan.md, R2, R3, R6, C2 (b);
// "Tests for @tester" T2 and T2-seg). INVARIANT-CRITICAL, with
// test_mesh_row_spans.
//
// T2: refinement::scan and the frozen box walk scan_oracle::scan_bbox
// (tests/cpp/support/scan_oracle.hpp) give identical ScanResults on every
// triangle: max_error by memcmp, node, where, is_void, uncovered. Exact
// equality, per Ola's C1 (a). The DEMs carry NaN NoData, sentinel NoData,
// NoData at vertices (void triangles, R3's carve point) and repeated values
// (ties, the smallest-(row, col) rule), under two frames with dx != dy.
//
// T2-seg: the same through SegmentedRaster, a test-only RasterSource over the
// same Raster whose row_segments cuts every requested range at fixed columns,
// at every column, exactly at c0 and c1, and at seeded random columns. All
// segments share the source's one nodata() (Ola's C2, simplified: no
// per-segment sentinel). Its row() is POISONED -- every cell 1e30 -- so a scan
// that bypasses for_each_row_segment differs, and it counts the segments it
// hands out, so a scan that never asks for one fails.
//
// Interface assumed (R6 and C2 (b); names the design leaves open are chosen here):
//   namespace terrain::raster, in <terrain/raster/row_segments.hpp>:
//     template <typename T> struct RowSegment {
//         std::span<const T> values; std::uint32_t first_col; };
//     template <RasterSource R, std::invocable<RowSegment<typename R::value_type>> F>
//     void for_each_row_segment(const R&, std::size_t row, std::uint32_t c0,
//                               std::uint32_t c1, F&& f);
//   A source's optional hook is
//     template <typename F> void row_segments(std::size_t, std::uint32_t,
//                                             std::uint32_t, F&&) const;
//   RasterSource requires
//     { r.nodata() } -> std::same_as<const std::optional<value_type>&>.
//   The design's single-grid branch passes a third initialiser, dem.nodata(),
//   to a RowSegment that has no nodata field; this suite builds RowSegment
//   with designated initialisers for the two fields only.

#include <catch2/catch_test_macros.hpp>

#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/row_segments.hpp>
#include <terrain/raster/view.hpp>
#include <terrain/refinement/scan.hpp>

#include "refinement_fixtures.hpp"
#include "scan_oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <span>
#include <string>
#include <vector>

using terrain::TriangleIndices;
using terrain::mesh::LatticeMesh;
using terrain::mesh::MeshVertex;
using terrain::mesh::orient_sign;
using terrain::raster::CellIndex;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::raster::RasterSource;
using terrain::raster::RasterView;
using terrain::raster::RowSegment;
using terrain::refinement::ScanResult;

namespace {

constexpr std::size_t kN = 24;
constexpr float kSentinel = -32767.0f;
const float kNaN = std::numeric_limits<float>::quiet_NaN();

// ---------------------------------------------------------------- the adapter

enum class Cut { Whole, EveryColumn, Every3, Every7, AtEnds, Random };

class SegmentedRaster {
public:
    using value_type = float;

    SegmentedRaster(const Raster<float>& base, Cut cut, std::uint32_t seed)
        : base_{&base}, cut_{cut}, seed_{seed},
          poison_(base.geometry().cols(), 1e30f) {}

    [[nodiscard]] const RasterGeometry& geometry() const noexcept { return base_->geometry(); }
    [[nodiscard]] float value_at(const CellIndex& c) const noexcept { return base_->value_at(c); }
    [[nodiscard]] bool is_nodata(const CellIndex& c) const noexcept { return base_->is_nodata(c); }
    [[nodiscard]] const std::optional<float>& nodata() const noexcept { return base_->nodata(); }
    // Poisoned: the scan must reach cells through row_segments only.
    [[nodiscard]] std::span<const float> row(std::size_t) const noexcept { return poison_; }

    template <typename F>
    void row_segments(std::size_t i, std::uint32_t c0, std::uint32_t c1, F&& f) const {
        ++*segments_;
        const auto whole = base_->row(i);
        std::vector<std::uint32_t> starts{c0};  // each segment's first column
        auto cut_at = [&](std::uint32_t c) { if (c > c0 && c <= c1) starts.push_back(c); };
        switch (cut_) {
            case Cut::Whole: break;
            case Cut::EveryColumn: for (std::uint32_t c = c0 + 1; c <= c1; ++c) cut_at(c); break;
            case Cut::Every3: for (std::uint32_t c = 0; c <= c1; c += 3) cut_at(c); break;
            case Cut::Every7: for (std::uint32_t c = 0; c <= c1; c += 7) cut_at(c); break;
            case Cut::AtEnds: cut_at(c0 + 1); cut_at(c1); break;  // [c0], [c0+1, c1-1], [c1]
            case Cut::Random: {
                std::mt19937 g{seed_ ^ static_cast<std::uint32_t>(i * 7919u + c0 * 104729u + c1)};
                for (std::uint32_t c = c0 + 1; c <= c1; ++c)
                    if (g() % 3 == 0) cut_at(c);
                break;
            }
        }
        std::sort(starts.begin(), starts.end());
        starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
        for (std::size_t k = 0; k < starts.size(); ++k) {
            const std::uint32_t a = starts[k], b = k + 1 < starts.size() ? starts[k + 1] - 1 : c1;
            f(RowSegment<float>{.values = whole.subspan(a, b - a + 1), .first_col = a});
        }
    }

    [[nodiscard]] std::size_t segments_handed_out() const noexcept { return *segments_; }

private:
    const Raster<float>* base_;
    Cut cut_;
    std::uint32_t seed_;
    std::vector<float> poison_;
    std::shared_ptr<std::size_t> segments_ = std::make_shared<std::size_t>(0);
};

// ---------------------------------------------------------------- DEMs

struct Dem {
    std::string name;
    Raster<float> raster;
};

std::vector<float> with_holes(std::vector<float> z, std::uint32_t seed, float hole, int every) {
    std::mt19937 g{seed};
    for (auto& x : z)
        if (g() % static_cast<std::uint32_t>(every) == 0) x = hole;
    return z;
}

std::vector<float> two_levels(std::uint32_t seed) {  // many exact ties
    std::mt19937 g{seed};
    std::vector<float> z(kN * kN);
    for (auto& x : z) x = static_cast<float>(g() % 2);
    return z;
}

std::vector<Dem> dems(const RasterGeometry& geo) {
    using namespace refinement_fixtures;
    std::vector<Dem> out;
    out.push_back({"rough", Raster<float>{geo, rough_dem(kN, kN, 7)}});
    out.push_back({"smooth", Raster<float>{geo, smooth_dem(kN, kN, 8)}});
    out.push_back({"NaN NoData", Raster<float>{geo, with_holes(rough_dem(kN, kN, 9), 10, kNaN, 11)}});
    out.push_back({"sentinel NoData",
                   Raster<float>{geo, with_holes(smooth_dem(kN, kN, 11), 12, kSentinel, 9), kSentinel}});
    out.push_back({"ties", Raster<float>{geo, two_levels(13)}});
    out.push_back({"ties with sentinel", Raster<float>{geo, with_holes(two_levels(14), 15, kSentinel, 13), kSentinel}});
    return out;
}

// Two frames, dx != dy both ways; the scan works in the lattice and must not care.
std::vector<RasterGeometry> frames() {
    return {refinement_fixtures::geometry(kN, kN),
            RasterGeometry{500000.0, 7000000.0, 5.0, 10.0, kN, kN}};
}

// ---------------------------------------------------------------- triangles

MeshVertex mv(double col, double row) { return MeshVertex{col, row}; }

struct Rng {
    std::mt19937 gen;
    double node() { return static_cast<double>(gen() % kN); }
    double coord() {
        return std::min(node() + static_cast<double>(gen() % 1024u) / 1024.0, double(kN - 1));
    }
    bool coin() { return (gen() & 1u) != 0; }
};

// Seeded triangles, node-only, mixed, near-node, slivers and horizontal, each as
// a one-triangle LatticeMesh so every result is one triangle's.
std::vector<LatticeMesh> meshes(std::uint32_t seed, int count) {
    Rng g{std::mt19937{seed}};
    std::vector<LatticeMesh> out;
    for (int i = 0; i < count; ++i) {
        std::array<MeshVertex, 3> v{};
        switch (i % 5) {
            case 0: for (auto& p : v) p = mv(g.node(), g.node()); break;
            case 1: for (auto& p : v) p = g.coin() ? mv(g.node(), g.node()) : mv(g.coord(), g.coord()); break;
            case 2:
                for (auto& p : v) {
                    const double e = g.coin() ? 1e-9 : -1e-12;
                    p = mv(std::clamp(g.node() + e, 0.0, double(kN - 1)), g.node());
                }
                break;
            case 3: {
                v[0] = mv(g.node(), g.node());
                v[1] = mv(g.node(), g.node());
                v[2] = mv(std::clamp((v[0].col + v[1].col) / 2 + 1e-7, 0.0, double(kN - 1)),
                          (v[0].row + v[1].row) / 2);
                break;
            }
            default: {
                const double r = g.node();
                v = {mv(g.coord(), r), mv(g.coord(), r), mv(g.coord(), g.coord())};
            }
        }
        const int s = orient_sign(v[0], v[1], v[2]);
        if (s == 0) continue;
        if (s < 0) std::swap(v[1], v[2]);
        auto m = LatticeMesh::build({v[0], v[1], v[2]}, {TriangleIndices{0, 1, 2}}, {0},
                                    {std::array<std::uint32_t, 3>{}});
        if (m) out.push_back(std::move(*m));
    }
    return out;
}

// ---------------------------------------------------------------- comparison

bool same(const ScanResult& a, const ScanResult& b) {
    return std::memcmp(&a.max_error, &b.max_error, sizeof(double)) == 0 && a.node == b.node
        && a.where == b.where && a.is_void == b.is_void && a.uncovered == b.uncovered;
}

struct Tally {
    int triangles = 0, voids = 0, argmax = 0, mismatches = 0;
};

template <RasterSource R>
Tally compare(const R& subject, const Raster<float>& oracle_dem, const std::vector<LatticeMesh>& ms) {
    Tally t;
    for (const auto& m : ms) {
        const ScanResult want = scan_oracle::scan_bbox(oracle_dem, m, 0);
        const ScanResult got = terrain::refinement::scan(subject, m, 0);
        ++t.triangles;
        t.voids += want.is_void;
        t.argmax += want.node.has_value() && !want.is_void;
        if (!same(got, want) && ++t.mismatches <= 3) {
            const auto c = [&](unsigned k) { return m.corner(0, k); };
            FAIL_CHECK("triangle (col, row) (" << c(0).col << ", " << c(0).row << ") ("
                       << c(1).col << ", " << c(1).row << ") (" << c(2).col << ", " << c(2).row
                       << "): max_error " << got.max_error << " vs " << want.max_error
                       << ", void " << got.is_void << " vs " << want.is_void << ", uncovered "
                       << got.uncovered << " vs " << want.uncovered);
        }
    }
    return t;
}

template <typename R>
concept HasSourceNodata = requires(const R& r) {
    { r.nodata() } -> std::same_as<const std::optional<typename R::value_type>&>;
};

// Everything RasterSource required before C2 (b), and no nodata().
class NodatalessRaster {
public:
    using value_type = double;
    explicit NodatalessRaster(RasterGeometry g) : geometry_{g}, row_(g.cols(), 0.0) {}
    [[nodiscard]] const RasterGeometry& geometry() const noexcept { return geometry_; }
    [[nodiscard]] double value_at(const CellIndex&) const noexcept { return 0.0; }
    [[nodiscard]] bool is_nodata(const CellIndex&) const noexcept { return false; }
    [[nodiscard]] std::span<const double> row(std::size_t) const noexcept { return row_; }

private:
    RasterGeometry geometry_;
    std::vector<double> row_;
};

}  // namespace

// ---------------------------------------------------------------- C2 (b)

TEST_CASE("C2 RasterSource requires one source-wide nodata()", "[raster][concept][C2]") {
    CHECK(HasSourceNodata<Raster<float>>);
    CHECK(HasSourceNodata<RasterView<float>>);
    CHECK(HasSourceNodata<RasterView<double>>);
    CHECK(RasterSource<SegmentedRaster>);
    CHECK_FALSE(RasterSource<NodatalessRaster>);
}

TEST_CASE("C2 RasterView::nodata returns the sentinel it was built with", "[raster][view][C2]") {
    const auto g = refinement_fixtures::geometry(2, 2);
    const std::vector<float> z{1.0f, kSentinel, 3.0f, 4.0f};
    const RasterView<float> with{g, z.data(), kSentinel};
    const RasterView<float> without{g, z.data(), std::nullopt};
    REQUIRE(with.nodata().has_value());
    CHECK(*with.nodata() == kSentinel);
    CHECK_FALSE(without.nodata().has_value());
}

// ---------------------------------------------------------------- R6 directly

TEST_CASE("for_each_row_segment on one grid is one segment of row(r)", "[raster][row_segments]") {
    const Raster<float> dem{refinement_fixtures::geometry(kN, kN), refinement_fixtures::rough_dem(kN, kN, 3)};
    int calls = 0;
    terrain::raster::for_each_row_segment(dem, 5, 3, 17, [&](RowSegment<float> s) {
        ++calls;
        CHECK(s.first_col == 3u);
        CHECK(s.values.data() == dem.row(5).data() + 3);
        CHECK(s.values.size() == 15u);
    });
    CHECK(calls == 1);
}

TEST_CASE("for_each_row_segment forwards to row_segments in ascending order", "[raster][row_segments]") {
    const Raster<float> dem{refinement_fixtures::geometry(kN, kN), refinement_fixtures::rough_dem(kN, kN, 4)};
    for (const Cut cut : {Cut::Whole, Cut::EveryColumn, Cut::Every3, Cut::AtEnds, Cut::Random}) {
        const SegmentedRaster seg{dem, cut, 99};
        std::uint32_t next = 2;
        terrain::raster::for_each_row_segment(seg, 7, 2, 20, [&](RowSegment<float> s) {
            CHECK(s.first_col == next);
            CHECK_FALSE(s.values.empty());
            CHECK(s.values.data() == dem.row(7).data() + s.first_col);
            next = s.first_col + static_cast<std::uint32_t>(s.values.size());
        });
        CHECK(next == 21u);
        CHECK(seg.segments_handed_out() == 1u);
    }
}

// ---------------------------------------------------------------- T2

TEST_CASE("T2 scan equals the frozen box walk bit for bit", "[refinement][scan][T2]") {
    const auto ms = meshes(1811, 3000);
    REQUIRE(ms.size() > 2000);
    for (const auto& geo : frames())
        for (const auto& d : dems(geo)) {
            INFO("DEM " << d.name << ", dx " << geo.delta_x() << ", dy " << geo.delta_y());
            const Tally t = compare(d.raster, d.raster, ms);
            CHECK(t.mismatches == 0);
            CHECK(t.argmax > 0);
            if (d.name.find("NoData") != std::string::npos || d.name.find("sentinel") != std::string::npos)
                CHECK(t.voids > 0);  // the carve point is exercised, not just the argmax
        }
}

// ---------------------------------------------------------------- T2-seg

TEST_CASE("T2-seg scan through segmented rows equals the box walk bit for bit",
          "[refinement][scan][T2-seg]") {
    const auto ms = meshes(1812, 1500);
    for (const auto& d : dems(frames().front()))
        for (const Cut cut : {Cut::Whole, Cut::EveryColumn, Cut::Every3, Cut::Every7, Cut::AtEnds, Cut::Random}) {
            INFO("DEM " << d.name << ", cut " << static_cast<int>(cut));
            const SegmentedRaster seg{d.raster, cut, 1813};
            const Tally t = compare(seg, d.raster, ms);
            CHECK(t.mismatches == 0);
            CHECK(seg.segments_handed_out() > 0);  // the scan took the segmented route
        }
}
