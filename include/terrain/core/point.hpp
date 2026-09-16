#pragma once

#include <cstddef>
#include <format>
#include <functional>
#include <string_view>

namespace terrain {

struct Point2 {
    double x{};
    double y{};

    friend constexpr bool operator==(const Point2&, const Point2&) = default;

    constexpr Point2& operator+=(const Point2& rhs) noexcept { x += rhs.x; y += rhs.y; return *this; }
    constexpr Point2& operator-=(const Point2& rhs) noexcept { x -= rhs.x; y -= rhs.y; return *this; }
    constexpr Point2& operator*=(double s) noexcept { x *= s; y *= s; return *this; }
    constexpr Point2& operator/=(double s) noexcept { x /= s; y /= s; return *this; }

    [[nodiscard]] constexpr Point2 operator-() const noexcept { return {-x, -y}; }
};

[[nodiscard]] constexpr Point2 operator+(Point2 a, const Point2& b) noexcept { a += b; return a; }
[[nodiscard]] constexpr Point2 operator-(Point2 a, const Point2& b) noexcept { a -= b; return a; }
[[nodiscard]] constexpr Point2 operator*(Point2 p, double s) noexcept { p *= s; return p; }
[[nodiscard]] constexpr Point2 operator*(double s, Point2 p) noexcept { p *= s; return p; }
[[nodiscard]] constexpr Point2 operator/(Point2 p, double s) noexcept { p /= s; return p; }

[[nodiscard]] constexpr double dot(const Point2& a, const Point2& b) noexcept {
    return a.x * b.x + a.y * b.y;
}

[[nodiscard]] constexpr double cross(const Point2& a, const Point2& b) noexcept {
    return a.x * b.y - a.y * b.x;
}

struct Point3 {
    double x{};
    double y{};
    double z{};

    friend constexpr bool operator==(const Point3&, const Point3&) = default;

    constexpr Point3& operator+=(const Point3& rhs) noexcept { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
    constexpr Point3& operator-=(const Point3& rhs) noexcept { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
    constexpr Point3& operator*=(double s) noexcept { x *= s; y *= s; z *= s; return *this; }
    constexpr Point3& operator/=(double s) noexcept { x /= s; y /= s; z /= s; return *this; }

    [[nodiscard]] constexpr Point3 operator-() const noexcept { return {-x, -y, -z}; }
};

[[nodiscard]] constexpr Point3 operator+(Point3 a, const Point3& b) noexcept { a += b; return a; }
[[nodiscard]] constexpr Point3 operator-(Point3 a, const Point3& b) noexcept { a -= b; return a; }
[[nodiscard]] constexpr Point3 operator*(Point3 p, double s) noexcept { p *= s; return p; }
[[nodiscard]] constexpr Point3 operator*(double s, Point3 p) noexcept { p *= s; return p; }
[[nodiscard]] constexpr Point3 operator/(Point3 p, double s) noexcept { p /= s; return p; }

[[nodiscard]] constexpr double dot(const Point3& a, const Point3& b) noexcept {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

[[nodiscard]] constexpr Point3 cross(const Point3& a, const Point3& b) noexcept {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
}

} // namespace terrain

namespace std {

// Consistent with the defaulted operator==: equal points hash equally.
// Both libc++ and libstdc++ normalize -0.0 in std::hash<double>, so +0.0 and
// -0.0 hash alike -- which is required here, since they also compare equal.
// NaN coordinates never compare equal, so a point containing one is never
// findable in a hashed container. Coordinate dedup that must be robust to
// floating drift should snap-round upstream, not lean on this.
template <>
struct hash<terrain::Point2> {
    [[nodiscard]] size_t operator()(const terrain::Point2& p) const noexcept {
        const size_t hx = hash<double>{}(p.x);
        const size_t hy = hash<double>{}(p.y);
        return hx ^ (hy + 0x9e3779b97f4a7c15ULL + (hx << 6) + (hx >> 2));
    }
};

template <>
struct hash<terrain::Point3> {
    [[nodiscard]] size_t operator()(const terrain::Point3& p) const noexcept {
        const size_t hx = hash<double>{}(p.x);
        const size_t hy = hash<double>{}(p.y);
        const size_t hz = hash<double>{}(p.z);
        size_t h = hx;
        h ^= hy + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= hz + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

template <>
struct formatter<terrain::Point2> {
    // Inheriting formatter<string_view> also inherited its parse(), so a width
    // or align spec parsed successfully and was then discarded by a format()
    // that ignores it -- "{:>24}" silently produced unpadded output. Reject
    // what we do not honour rather than accept it and lie.
    constexpr auto parse(format_parse_context& ctx) {
        auto it = ctx.begin();
        if (it != ctx.end() && *it != '}')
            throw format_error("terrain::Point2 does not accept a format spec");
        return it;
    }

    template <typename FormatContext>
    auto format(const terrain::Point2& p, FormatContext& ctx) const {
        return std::format_to(ctx.out(), "Point2({}, {})", p.x, p.y);
    }
};

template <>
struct formatter<terrain::Point3> {
    // Inheriting formatter<string_view> also inherited its parse(), so a width
    // or align spec parsed successfully and was then discarded by a format()
    // that ignores it -- "{:>24}" silently produced unpadded output. Reject
    // what we do not honour rather than accept it and lie.
    constexpr auto parse(format_parse_context& ctx) {
        auto it = ctx.begin();
        if (it != ctx.end() && *it != '}')
            throw format_error("terrain::Point3 does not accept a format spec");
        return it;
    }

    template <typename FormatContext>
    auto format(const terrain::Point3& p, FormatContext& ctx) const {
        return std::format_to(ctx.out(), "Point3({}, {}, {})", p.x, p.y, p.z);
    }
};

} // namespace std
