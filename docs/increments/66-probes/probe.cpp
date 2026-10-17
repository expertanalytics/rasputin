// Increment 66 probe (docs/increments/66-coincident-written-vertices.md,
// section 5): fixture W1 through refine, constraint_check_points and
// refine_strip, one thread, feet off then on. Prints each outline with a
// coincident vertex (one whose written point equals an earlier vertex's), a
// flat triangle (flat_as_written) or a refused run, and the totals. With a
// fifth argument it also prints each flat triangle's vertices.
//
// Build from the repository root (the single-rule figures of section 5 come
// from building it against the mutants "R1 removed" and "R2 removed"):
//   c++ -std=c++20 -O2 -Iinclude -Ilib/detria -c src/predicates/detria_exact.cpp -o /tmp/d.o
//   c++ -std=c++20 -O2 -Iinclude docs/increments/66-probes/probe.cpp /tmp/d.o -o /tmp/probe
//   /tmp/probe <tolerance> <k_first> <k_last> [detail]
#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/refine_points.hpp>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
using namespace terrain;
using namespace terrain::refinement;
constexpr double kX0 = 500000.0, kY0 = 6600000.0;
raster::RasterGeometry geo() { return raster::RasterGeometry{kX0, kY0, 10.0, 10.0, 12, 12}; }
raster::Raster<float> dem_() {
    auto g = geo(); std::vector<float> z(144);
    for (std::size_t r = 0; r < 12; ++r) for (std::size_t c = 0; c < 12; ++c)
        z[r * 12 + c] = float(int((r * 7 + c * 13) % 11) - 5);
    return raster::Raster<float>{g, std::move(z)};
}
double h(int k, int i, int lo, int hi) { return double(lo + (std::int64_t(k) * (2 * i + 7) * 7919 + i * 104729) % (hi - lo)) / 10.0; }
int main(int argc, char** argv) {
    const double tol = std::atof(argv[1]);
    const int k0 = std::atoi(argv[2]), k1 = std::atoi(argv[3]);
    const bool detail_ = argc > 4;
    auto dem = dem_();
    std::size_t tot_pairs = 0, tot_flat = 0, tot_ref = 0;
    for (int feet = 0; feet < 2; ++feet) {
      std::vector<int> bad;
      for (int k = k0; k <= k1; ++k) {
        std::vector<Point2> xy{{kX0 + h(k,0,50,300), kY0 - h(k,1,700,1050)}, {kX0 + h(k,2,800,1050), kY0 - h(k,3,750,1050)},
                               {kX0 + h(k,4,700,1050), kY0 - h(k,5,50,300)}, {kX0 + h(k,6,50,400), kY0 - h(k,7,50,300)}};
        IndexedMesh2 start{xy, std::vector<TriangleIndices>{{0,1,2},{0,2,3}}, std::vector<std::uint8_t>(2,0)};
        std::vector<std::array<std::uint32_t,2>> edges{{0,1},{1,2},{2,3},{3,0}}; std::vector<std::uint32_t> masks(4,1);
        RefineOptions o; o.tolerance = tol; o.threads = 1; o.constraint_feet = feet;
        auto r = refine(dem, start, edges, masks, o);
        if (!r.ok()) { std::printf("phase1 fail k=%d\n", k); continue; }
        const std::size_t nt = r.triangles.size();
        IndexedMesh2 mesh{r.vertices, r.triangles, std::vector<std::uint8_t>(nt, 0)};
        auto strip = constraint_check_points(dem, mesh.vertices(), std::span<const std::array<std::uint32_t,2>>{r.edges});
        PointRefineOptions po; po.tolerance = tol; po.threads = 1; po.constraint_feet = feet;
        auto out = refine_strip(dem, strip, mesh, std::span<const double>{r.z}, std::span<const std::uint8_t>{r.valid},
                                std::span<const std::array<std::uint32_t,2>>{r.edges}, std::span<const std::uint32_t>{r.masks}, po);
        if (!out.ok()) { ++tot_ref; bad.push_back(k); std::printf("feet=%d k=%d REFUSED %s\n", feet, k, out.message.c_str()); continue; }
        std::map<std::pair<double,double>, int> seen; std::size_t pairs = 0;
        for (auto& v : out.vertices) pairs += seen[{v.x, v.y}]++ ? 1 : 0;
        tot_pairs += pairs; tot_flat += out.flat_as_written;
        if (pairs || out.flat_as_written) { bad.push_back(k); std::printf("feet=%d k=%d coincident=%zu flat=%zu\n", feet, k, pairs, out.flat_as_written); }
        if (detail_ && out.flat_as_written) {
          std::set<std::pair<std::uint32_t,std::uint32_t>> con; for (auto& e : out.edges) con.insert(std::minmax(e[0], e[1]));
          auto g = geo();
          for (auto& t : out.triangles) {
            auto& v = out.vertices;
            if (pred::DefaultKernel::orient2d(v[t[0]], v[t[1]], v[t[2]]) == pred::Orientation::CounterClockwise) continue;
            std::printf("  flat tri (%u %u %u)\n", t[0], t[1], t[2]);
            for (int i = 0; i < 3; ++i) { auto lp = detail::lattice_position(g, v[t[i]]);
              std::printf("    v%u world (%.17g, %.17g) lattice (%.17g, %.17g) z %.17g\n", t[i], v[t[i]].x, v[t[i]].y, lp.col, lp.row, out.z[t[i]]); }
            for (int i = 0; i < 3; ++i) std::printf("    edge %u-%u constrained=%d\n", t[i], t[(i+1)%3], (int)con.count(std::minmax(t[i], t[(i+1)%3])));
          }
        }
      }
    }
    std::printf("TOTAL coincident=%zu flat=%zu refused=%zu\n", tot_pairs, tot_flat, tot_ref);
}
