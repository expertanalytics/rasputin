#pragma once

#include <cstddef>
#include <format>
#include <functional>

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

} // namespace terrain

namespace std {

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

} // namespace std
