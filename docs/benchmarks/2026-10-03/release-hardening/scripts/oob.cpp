// Probe: does -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST trap on an
// out-of-range std::vector / std::span index at -O3 -DNDEBUG?
#include <cstdio>
#include <span>
#include <vector>
int main(int argc, char**) {
    std::vector<int> v(4, 7);
    std::span<int> s(v);
    volatile int i = 3 + argc;  // 4 with no arguments: one past the end
    std::printf("%d\n", v[static_cast<std::size_t>(i)] + s[static_cast<std::size_t>(i)]);
}
