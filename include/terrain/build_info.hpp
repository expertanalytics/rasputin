#pragma once

// What this build is, as the compiled code sees it: today only whether the
// standard library checks its own preconditions (increment 24,
// docs/increments/24-release-hardening.md section 3).
//
// The answer is read from the macros the standard library itself defines, never
// from the CMake option that asked for them, so a toolchain that ignores the
// request (libc++ before 18 does not define the mode constants) reports "none".
// The Python attribute `_core.hardening` and the guard test
// tests/cpp/unit/test_build_hardening.cpp both read this one function.
// Free of pybind11, so the guard can include it.

#include <string_view>
#include <version>  // _LIBCPP_VERSION / __GLIBCXX__ and the libc++ mode constants

namespace terrain {

[[nodiscard]] constexpr std::string_view stdlib_hardening() noexcept {
#if defined(_LIBCPP_HARDENING_MODE_FAST) && _LIBCPP_HARDENING_MODE == _LIBCPP_HARDENING_MODE_FAST
    return "libc++ fast";
#elif defined(_LIBCPP_HARDENING_MODE_EXTENSIVE) && \
    _LIBCPP_HARDENING_MODE == _LIBCPP_HARDENING_MODE_EXTENSIVE
    return "libc++ extensive";
#elif defined(_LIBCPP_HARDENING_MODE_DEBUG) && _LIBCPP_HARDENING_MODE == _LIBCPP_HARDENING_MODE_DEBUG
    return "libc++ debug";
#elif defined(__GLIBCXX__) && defined(_GLIBCXX_ASSERTIONS)
    return "libstdc++ assertions";
#else
    return "none";
#endif
}

}  // namespace terrain
