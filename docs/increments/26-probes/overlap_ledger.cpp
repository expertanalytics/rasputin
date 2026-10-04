// Throwaway probe for docs/increments/26-land-cover.md, not production code and
// not built by CMake. Compile: c++ -std=c++20 -O2 -o overlap_ledger overlap_ledger.cpp
//
// overlap_ledger overlap DIR   reads DIR/tri.f64 (T x 6: col,row of three vertices,
//                              in the class raster's index space), DIR/win.u8 (H x W),
//                              DIR/win.shape (H W), DIR/roww.f64 (H: m^2 per cell, by row);
//                              writes DIR/csr.u32 (T+1 offsets), DIR/cls.u8, DIR/area.f64
//                              (exact triangle-cell overlap, summed per class, in m^2).
// overlap_ledger ledger DIR C M RULE   reads the CSR and DIR/order.u32 (the visiting order);
//                              an entry is too small when below C (relative to the triangle)
//                              and below M (m^2) both; RULE kept|present|any|renorm;
//                              writes DIR/out_<RULE>_<C>_<M>.{csr.u32,cls.u8,area.f64}.
//                              RASPUTIN_PROBE_CAP=N in the environment keeps at most N
//                              entries per triangle (the N largest), the ledger taking the rest.
// overlap_ledger label DIR CLASSES   8-connected components of the cells whose class is in
//                              the comma list CLASSES (water_probe.py's first step).
//                              and prints the ledger's largest |L_k| and the end remainder.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template <typename T> std::vector<T> slurp(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    const auto n = static_cast<std::size_t>(f.tellg());
    std::vector<T> v(n / sizeof(T));
    f.seekg(0);
    f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n));
    return v;
}
template <typename T> void spit(const std::string& path, const std::vector<T>& v) {
    std::ofstream f(path, std::ios::binary);
    f.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)));
}

struct P { double x, y; };
using Poly = std::vector<P>;

// Sutherland-Hodgman against one axis-parallel half-plane: keep side * (coord - v) >= 0.
static Poly clip(const Poly& in, int axis, double v, double side) {
    Poly out;
    const std::size_t n = in.size();
    for (std::size_t i = 0; i < n; ++i) {
        const P& a = in[i];
        const P& b = in[(i + 1) % n];
        const double da = side * ((axis ? a.y : a.x) - v), db = side * ((axis ? b.y : b.x) - v);
        if (da >= 0) out.push_back(a);
        if ((da >= 0) != (db >= 0)) {
            const double t = da / (da - db);
            out.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)});
        }
    }
    return out;
}
static double area(const Poly& p) {
    double s = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const P& a = p[i];
        const P& b = p[(i + 1) % p.size()];
        s += a.x * b.y - a.y * b.x;
    }
    return 0.5 * std::fabs(s);
}

static int overlap(const std::string& dir) {
    const auto tri = slurp<double>(dir + "/tri.f64");
    const auto win = slurp<std::uint8_t>(dir + "/win.u8");
    const auto shape = slurp<std::int64_t>(dir + "/win.shape");
    const auto roww = slurp<double>(dir + "/roww.f64");
    const std::int64_t H = shape[0], W = shape[1];
    const std::size_t T = tri.size() / 6;
    std::vector<std::uint32_t> off{0};
    std::vector<std::uint8_t> cls;
    std::vector<double> ar;
    std::array<double, 256> acc{};
    std::vector<int> touched;
    for (std::size_t t = 0; t < T; ++t) {
        Poly p{{tri[6 * t], tri[6 * t + 1]}, {tri[6 * t + 2], tri[6 * t + 3]}, {tri[6 * t + 4], tri[6 * t + 5]}};
        double y0 = 1e300, y1 = -1e300;
        for (auto& q : p) y0 = std::min(y0, q.y), y1 = std::max(y1, q.y);
        for (std::int64_t r = static_cast<std::int64_t>(std::floor(y0)); r < static_cast<std::int64_t>(std::ceil(y1)); ++r) {
            const Poly strip = clip(clip(p, 1, double(r), 1), 1, double(r + 1), -1);
            if (strip.size() < 3) continue;
            double x0 = 1e300, x1 = -1e300;
            for (auto& q : strip) x0 = std::min(x0, q.x), x1 = std::max(x1, q.x);
            for (std::int64_t c = static_cast<std::int64_t>(std::floor(x0)); c < static_cast<std::int64_t>(std::ceil(x1)); ++c) {
                const Poly cell = clip(clip(strip, 0, double(c), 1), 0, double(c + 1), -1);
                if (cell.size() < 3) continue;
                const double a = area(cell);
                if (!(a > 0)) continue;
                if (r < 0 || r >= H || c < 0 || c >= W) { std::cerr << "outside window\n"; return 2; }
                const int k = win[static_cast<std::size_t>(r * W + c)];
                if (acc[k] == 0.0) touched.push_back(k);
                acc[k] += a * roww[static_cast<std::size_t>(r)];
            }
        }
        std::sort(touched.begin(), touched.end());
        for (int k : touched) { cls.push_back(static_cast<std::uint8_t>(k)); ar.push_back(acc[k]); acc[k] = 0; }
        touched.clear();
        off.push_back(static_cast<std::uint32_t>(cls.size()));
    }
    spit(dir + "/csr.u32", off);
    spit(dir + "/cls.u8", cls);
    spit(dir + "/area.f64", ar);
    return 0;
}

static int ledger(const std::string& dir, double c, double m, const std::string& rule) {
    const auto off = slurp<std::uint32_t>(dir + "/csr.u32");
    const auto cls = slurp<std::uint8_t>(dir + "/cls.u8");
    const auto ar = slurp<double>(dir + "/area.f64");
    const auto order = slurp<std::uint32_t>(dir + "/order.u32");
    const std::size_t T = off.size() - 1;
    std::vector<std::uint32_t> ooff(T + 1, 0);
    std::vector<std::vector<std::pair<std::uint8_t, double>>> out(T);
    std::array<double, 256> L{}, a{}, v{};
    std::array<double, 256> worst{};
    const bool renorm = rule == "renorm";
    // At most this many entries per triangle (0: no cap), from the environment.
    const char* cap_env = std::getenv("RASPUTIN_PROBE_CAP");
    const std::size_t cap = cap_env ? static_cast<std::size_t>(std::stoul(cap_env)) : 0;
    double sum_area = 0;
    for (const std::uint32_t t : order) {
        double A = 0;
        a.fill(0);
        for (std::uint32_t j = off[t]; j < off[t + 1]; ++j) { a[cls[j]] = ar[j]; A += ar[j]; }
        sum_area += A;
        std::vector<int> S;
        for (int k = 0; k < 256; ++k) {
            v[k] = a[k] + (renorm ? 0.0 : L[k]);
            const auto big = [&](double x) { return x >= c * A || x >= m; };
            bool eligible = rule == "any" || (rule == "present" && a[k] > 0) ||
                            ((rule == "kept" || renorm) && big(a[k]));
            if (eligible && big(v[k]) && v[k] > 0) S.push_back(k);
        }
        if (S.empty()) {  // fallback: the largest v among classes present in the triangle
            int best = -1;
            for (int k = 0; k < 256; ++k)
                if (a[k] > 0 && (best < 0 || v[k] > v[best])) best = k;
            S.push_back(best);
        }
        if (cap > 0 && S.size() > cap) {  // keep the `cap` largest (ties to the smaller code)
            std::stable_sort(S.begin(), S.end(), [&](int x, int y) { return v[x] > v[y]; });
            S.resize(cap);
        }
        std::array<double, 256> o{};
        for (;;) {  // proportional fill of A over S; drop entries pushed below the cutoff
            double s = 0;
            for (int k : S) s += std::max(v[k], 0.0);
            if (!(s > 0)) { o.fill(0); o[S.front()] = A; break; }
            std::vector<int> keep;
            for (int k : S) { o[k] = A * std::max(v[k], 0.0) / s; if (o[k] >= c * A || o[k] >= m) keep.push_back(k); }
            if (keep.size() == S.size() || keep.empty()) break;
            for (int k : S) o[k] = 0;
            S = keep;
        }
        for (int k = 0; k < 256; ++k) {
            if (o[k] > 0) out[t].push_back({static_cast<std::uint8_t>(k), o[k]});
            if (!renorm) { L[k] = v[k] - o[k]; worst[k] = std::max(worst[k], std::fabs(L[k])); }
        }
    }
    std::vector<std::uint8_t> ocls;
    std::vector<double> oar;
    for (std::size_t t = 0; t < T; ++t) {
        for (auto& [k, x] : out[t]) { ocls.push_back(k); oar.push_back(x); }
        ooff[t + 1] = static_cast<std::uint32_t>(ocls.size());
    }
    char tag[64];
    std::snprintf(tag, sizeof tag, "%s_%g_%g", rule.c_str(), c, m);
    spit(dir + "/out_" + tag + ".csr.u32", ooff);
    spit(dir + "/out_" + tag + ".cls.u8", ocls);
    spit(dir + "/out_" + tag + ".area.f64", oar);
    std::printf("{\"rule\":\"%s\",\"cutoff\":%g,\"mean_area\":%.6g,\"worst\":{", rule.c_str(), c, sum_area / double(T));
    bool first = true;
    for (int k = 0; k < 256; ++k)
        if (worst[k] > 0 || L[k] != 0) {
            std::printf("%s\"%d\":[%.6g,%.6g]", first ? "" : ",", k, worst[k], L[k]);
            first = false;
        }
    std::printf("}}\n");
    return 0;
}

// 8-connected components of the cells whose class is in `classes` (a comma list).
// Writes DIR/labels.i32 (H x W, 0 = not water, else 1-based component id in
// raster order of first cell) and DIR/sizes.i64 (cells per component, index = id).
static int label(const std::string& dir, const std::string& classes) {
    const auto win = slurp<std::uint8_t>(dir + "/win.u8");
    const auto shape = slurp<std::int64_t>(dir + "/win.shape");
    const std::int64_t H = shape[0], W = shape[1];
    std::array<bool, 256> water{};
    for (std::size_t p = 0; p < classes.size();) {
        const std::size_t q = classes.find(',', p);
        water[static_cast<std::size_t>(std::stoi(classes.substr(p, q - p)))] = true;
        p = q == std::string::npos ? classes.size() : q + 1;
    }
    std::vector<std::int32_t> lab(static_cast<std::size_t>(H * W), 0);
    std::vector<std::int32_t> parent{0};
    const auto find = [&](std::int32_t x) {
        while (parent[static_cast<std::size_t>(x)] != x) {
            parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    };
    for (std::int64_t r = 0; r < H; ++r)
        for (std::int64_t c = 0; c < W; ++c) {
            if (!water[win[static_cast<std::size_t>(r * W + c)]]) continue;
            std::int32_t best = 0;
            const std::int64_t nb[4][2] = {{r, c - 1}, {r - 1, c - 1}, {r - 1, c}, {r - 1, c + 1}};
            for (auto& n : nb) {
                if (n[0] < 0 || n[1] < 0 || n[1] >= W) continue;
                const std::int32_t l = lab[static_cast<std::size_t>(n[0] * W + n[1])];
                if (!l) continue;
                const std::int32_t root = find(l);
                if (!best) best = root;
                else if (root != best) {
                    const std::int32_t lo = std::min(root, best), hi = std::max(root, best);
                    parent[static_cast<std::size_t>(hi)] = lo;
                    best = lo;
                }
            }
            if (!best) { best = static_cast<std::int32_t>(parent.size()); parent.push_back(best); }
            lab[static_cast<std::size_t>(r * W + c)] = best;
        }
    std::vector<std::int32_t> id(parent.size(), 0);
    std::vector<std::int64_t> sizes{0};
    for (std::size_t i = 0; i < lab.size(); ++i) {
        if (!lab[i]) continue;
        const auto root = static_cast<std::size_t>(find(lab[i]));
        if (!id[root]) { id[root] = static_cast<std::int32_t>(sizes.size()); sizes.push_back(0); }
        lab[i] = id[root];
        ++sizes[static_cast<std::size_t>(lab[i])];
    }
    spit(dir + "/labels.i32", lab);
    spit(dir + "/sizes.i64", sizes);
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 3 && std::string(argv[1]) == "overlap") return overlap(argv[2]);
    if (argc >= 4 && std::string(argv[1]) == "label") return label(argv[2], argv[3]);
    if (argc >= 6 && std::string(argv[1]) == "ledger")
        return ledger(argv[2], std::stod(argv[3]), std::stod(argv[4]), argv[5]);
    std::cerr << "usage: overlap DIR | ledger DIR CUTOFF MIN_M2 RULE\n";
    return 1;
}
