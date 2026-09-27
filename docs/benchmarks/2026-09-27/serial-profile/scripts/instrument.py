import sys
root = sys.argv[1]
def patch(rel, pairs):
    p = f"{root}/{rel}"
    s = open(p).read()
    for old, new in pairs:
        assert s.count(old) == 1, (rel, old[:60])
        s = s.replace(old, new)
    open(p, "w").write(s)

patch("include/terrain/predicates/kernel.hpp", [
("namespace terrain::pred {\n", """namespace terrain::pred {

// PROFILE-ONLY (serial-profile, 2026-09-27): predicate counters, not committed.
namespace prof {
inline unsigned long long orient_calls = 0, orient_exact = 0, orient_exact_zero = 0;
inline unsigned long long incircle_calls = 0, incircle_exact = 0, incircle_exact_zero = 0;
}  // namespace prof
"""),
("""        if (std::fabs(det) > orient2d_bound_a * permanent) {
            return orientation_of_sign(det);
        }
        return E::orient2d(a, b, c);
    }

    // Total and precondition-free""", """        ++prof::orient_calls;
        if (std::fabs(det) > orient2d_bound_a * permanent) {
            return orientation_of_sign(det);
        }
        ++prof::orient_exact;
        const auto o = E::orient2d(a, b, c);
        prof::orient_exact_zero += o == Orientation::Collinear ? 1 : 0;
        return o;
    }

    // Total and precondition-free"""),
("""        if (std::fabs(det) > incircle_bound_a * permanent) {
            return incircle_of_sign(det);
        }
        return E::incircle_ccw(a, b, c, d);""", """        ++prof::incircle_calls;
        if (std::fabs(det) > incircle_bound_a * permanent) {
            return incircle_of_sign(det);
        }
        ++prof::incircle_exact;
        const auto r = E::incircle_ccw(a, b, c, d);
        prof::incircle_exact_zero += r == Incircle::Cocircular ? 1 : 0;
        return r;"""),
])

patch("include/terrain/mesh/lawson.hpp", [
("namespace terrain::mesh {\n", """namespace terrain::mesh {

// PROFILE-ONLY: when set, legalise_around appends every slot it reads (the
// popped triangle and its neighbour across the tested edge). Serial use only.
namespace prof {
inline std::vector<std::uint32_t>* reads = nullptr;
}  // namespace prof
"""),
("""        if (i == 3 || !detail::must_flip<K>(m, t, (i + 1) % 3, f))
            continue;""", """        if (prof::reads) {
            prof::reads->push_back(t);
            if (i < 3 && m.neighbours(t)[(i + 1) % 3] != kNoNeighbour)
                prof::reads->push_back(m.neighbours(t)[(i + 1) % 3]);
        }
        if (i == 3 || !detail::must_flip<K>(m, t, (i + 1) % 3, f))
            continue;"""),
])

patch("include/terrain/refinement/scan.hpp", [
("namespace terrain::refinement {\n", """namespace terrain::refinement {

// PROFILE-ONLY: nodes visited by scan on this thread.
namespace prof {
inline thread_local unsigned long long scan_nodes = 0;
}  // namespace prof
"""),
("""            const std::uint32_t row = span.row, c = s.first_col;""", """            const std::uint32_t row = span.row, c = s.first_col;
            prof::scan_nodes += s.values.size();"""),
])

patch("include/terrain/refinement/refine.hpp", [
("#include <algorithm>\n", "#include <algorithm>\n#include <cstdio>\n#include <cstdlib>\n#include <mutex>\n"),
("""        t0 = clock::now();
        parallel_util::for_each_chunk(active.size(), options.threads,
                                      [&](std::size_t begin, std::size_t end) {
                                          for (std::size_t i = begin; i < end; ++i)
                                              results[active[i]] = scan(dem, m, active[i]);
                                      });
        out.scan_seconds += since(t0);""", """        t0 = clock::now();
        struct ChunkRec { std::size_t begin, end; double start, stop; unsigned long long nodes; };
        std::vector<ChunkRec> chunk_recs;
        std::mutex chunk_mu;
        parallel_util::for_each_chunk(active.size(), options.threads,
                                      [&](std::size_t begin, std::size_t end) {
                                          const double c0 = since(t0);
                                          const auto n0 = prof::scan_nodes;
                                          for (std::size_t i = begin; i < end; ++i)
                                              results[active[i]] = scan(dem, m, active[i]);
                                          const double c1 = since(t0);
                                          const std::lock_guard lock{chunk_mu};
                                          chunk_recs.push_back({begin, end, c0, c1, prof::scan_nodes - n0});
                                      });
        const double scan_round = since(t0);
        out.scan_seconds += scan_round;
        if (prof_out) {
            std::fprintf(prof_out, "R %zu %zu %.9f %zu\\n", out.rounds, active.size(), scan_round, chunk_recs.size());
            std::sort(chunk_recs.begin(), chunk_recs.end(), [](auto& a, auto& b) { return a.begin < b.begin; });
            for (const auto& c : chunk_recs)
                std::fprintf(prof_out, "C %zu %zu %.9f %.9f %llu\\n", c.begin, c.end, c.start, c.stop, c.nodes);
        }
        const auto pc0 = pred::prof::incircle_calls, pe0 = pred::prof::incircle_exact,
                   pz0 = pred::prof::incircle_exact_zero, oc0 = pred::prof::orient_calls,
                   oe0 = pred::prof::orient_exact;
        std::size_t marked = 0, deferred_touched = 0, deferred_edge = 0;
        std::vector<std::uint32_t> reads, flip_writes;
        if (prof_out)
            mesh::prof::reads = &reads;"""),
("""    std::vector<ScanResult> results;
    std::set<std::pair""", """    // PROFILE-ONLY: RASPUTIN_PROF_OUT names a file for per-round records.
    std::FILE* prof_out = nullptr;
    if (const char* path = std::getenv("RASPUTIN_PROF_OUT"))
        prof_out = std::fopen(path, "w");
    std::vector<ScanResult> results;
    std::set<std::pair"""),
("""            any = true;
            if (touched[t] != 0)
                continue;""", """            any = true;
            ++marked;
            if (touched[t] != 0) {
                ++deferred_touched;
                continue;
            }
            const std::size_t flips_before = out.flips;
            reads.clear();
            flip_writes.clear();"""),
("""                if (u != mesh::kNoNeighbour && touched[u] != 0) {
                    skipped.push_back(t);
                    continue;
                }""", """                if (u != mesh::kNoNeighbour && touched[u] != 0) {
                    skipped.push_back(t);
                    ++deferred_edge;
                    continue;
                }"""),
("""                [&](std::uint32_t s) { touched[s] = 1; });""", """                [&](std::uint32_t s) { touched[s] = 1; if (prof_out) flip_writes.push_back(s); });"""),
("""            ++out.inserted;
            out.carved += r.is_void ? 1 : 0;""", """            ++out.inserted;
            if (prof_out) {
                // I col row on_edge flips | written slots | read slots | flip-written slots.
                // Written: t, the appended slots [before, count), u; flipped
                // slots are in the read list (a flip writes the pair it tested).
                std::fprintf(prof_out, "I %.3f %.3f %d %zu | %u", p.col, p.row, edge ? 1 : 0,
                             out.flips - flips_before, t);
                for (auto s2 = before; s2 < m.triangle_count(); ++s2)
                    std::fprintf(prof_out, " %u", s2);
                if (n_seeds == 4)
                    std::fprintf(prof_out, " %u", seeds[3]);
                std::fprintf(prof_out, " |");
                for (auto x : reads)
                    std::fprintf(prof_out, " %u", x);
                std::fprintf(prof_out, " |");
                for (auto x : flip_writes)
                    std::fprintf(prof_out, " %u", x);
                std::fprintf(prof_out, "\\n");
            }
            out.carved += r.is_void ? 1 : 0;"""),
("""        out.split_seconds += since(t0);
        if (!any)""", """        const double split_round = since(t0);
        out.split_seconds += split_round;
        mesh::prof::reads = nullptr;
        if (prof_out)
            std::fprintf(prof_out, "S %zu %zu %zu %zu %.9f %llu %llu %llu %llu %llu\\n", marked, deferred_touched,
                         deferred_edge, m.triangle_count(), split_round,
                         pred::prof::incircle_calls - pc0, pred::prof::incircle_exact - pe0,
                         pred::prof::incircle_exact_zero - pz0, pred::prof::orient_calls - oc0,
                         pred::prof::orient_exact - oe0);
        if (!any)"""),
("""    // By the stopping rule a void triangle""", """    if (prof_out)
        std::fclose(prof_out);
    // By the stopping rule a void triangle"""),
])
print("patched")
