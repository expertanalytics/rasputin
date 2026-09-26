#pragma once

// A 64-bit FNV-1a digest of everything `refine` returns that increment 18
// already returned: the six arrays, bit for bit, and the counters. Used by
// increment 20's Q7 (docs/increments/20-start-quality.md): with the quality
// pass off, the output must be bit-identical to increment 18's, and the
// digests Q7 compares against were recorded from increment 18's refine.hpp
// (master a9a93bc, unchanged at edde638) with this function, before any
// production change on the branch.
//
// `topology_digest` is the part Q7 pins across machines: indices, masks,
// validity and counters, all integers, so neither the compiler's FMA choice
// nor the platform's libm can move it while refine's decisions stand. `digest`
// adds every double bit for bit and is for comparisons inside one binary.
//
// The timings are left out: they are not part of the determinism guarantee.
// So are the fields increment 20 adds, so the digest is defined the same way
// before and after the change.

#include <cstdint>
#include <cstring>
#include <type_traits>

namespace refine_digest {

class Fnv {
public:
    template <typename T>
        requires std::is_trivially_copyable_v<T>
    void add(const T& value) noexcept {
        unsigned char bytes[sizeof(T)];
        std::memcpy(bytes, &value, sizeof(T));
        for (const unsigned char b : bytes) {
            h_ ^= b;
            h_ *= 1099511628211ull;
        }
    }
    [[nodiscard]] std::uint64_t value() const noexcept { return h_; }

private:
    std::uint64_t h_ = 14695981039346656037ull;
};

template <typename Outcome>
void add_topology(Fnv& f, const Outcome& out) {
    f.add(out.vertices.size());
    for (const auto v : out.valid) f.add(v);
    f.add(out.triangles.size());
    for (const auto& t : out.triangles)
        for (const auto i : t) f.add(i);
    f.add(out.edges.size());
    for (const auto& e : out.edges) {
        f.add(e[0]);
        f.add(e[1]);
    }
    for (const auto m : out.masks) f.add(m);
    f.add(out.rounds);
    f.add(out.inserted);
    f.add(out.flips);
    f.add(out.uncovered);
    f.add(out.carved);
}

template <typename Outcome>
[[nodiscard]] std::uint64_t topology_digest(const Outcome& out) {
    Fnv f;
    add_topology(f, out);
    return f.value();
}

template <typename Outcome>
[[nodiscard]] std::uint64_t digest(const Outcome& out) {
    Fnv f;
    add_topology(f, out);
    for (const auto& p : out.vertices) {
        f.add(p.x);
        f.add(p.y);
    }
    for (const double z : out.z) f.add(z);
    f.add(out.max_error);
    return f.value();
}

}  // namespace refine_digest
