// The release-hardening guard: increment 24, docs/increments/24-release-hardening.md
// sections 5 and 12 (T1 to T3).
//
// RASPUTIN_HARDENING (CMake, default ON) puts the standard library's own
// precondition checks into every build: _LIBCPP_HARDENING_MODE_FAST for libc++,
// _GLIBCXX_ASSERTIONS for libstdc++. A toolchain that does not know the macro
// ignores it silently (libc++ before 18), so this suite is what makes CI fail
// when "on" stops meaning on.
//
// tests/cpp/CMakeLists.txt gives this target RASPUTIN_EXPECT_HARDENING, 1 when
// the option is ON and 0 when it is OFF. Three guards read it:
//
//   T1  compile time: what the library reports must agree with what CMake asked
//       for, both ways (ON with an old libc++ fails, OFF with a stray definition
//       fails). Without the definition the suite does not compile at all.
//   T2  behaviour, ON only: a forked child reads one past the end of a vector
//       and must be killed by the check. Not WILL_FAIL: CTest does not invert a
//       signal death (section 2 of the design).
//   T3  the reported string, which the Python attribute `_core.hardening`, the
//       `rasputin version` line and the --stats line all repeat.
//
// Not invariant-critical, so no mutation round. The planting the design asks
// for (T1 red with the definitions removed, T2 red with an in-range index) is
// recorded in the red commit's message.

#include "terrain/build_info.hpp"

#include <catch2/catch_test_macros.hpp>

#include <csignal>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>
#include <version>

#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>

#ifndef RASPUTIN_EXPECT_HARDENING
#error "RASPUTIN_EXPECT_HARDENING is not defined: tests/cpp/CMakeLists.txt must pass it to test_build_hardening"
#endif

// T1: the compile-time guard, as section 5.1 states it.
static_assert((terrain::stdlib_hardening() != "none") == (RASPUTIN_EXPECT_HARDENING != 0),
              "the standard library's hardening mode disagrees with RASPUTIN_HARDENING");

namespace {

// The one-past-the-end read the child performs. The index is volatile so the
// compiler cannot prove it out of range and fold the access away, and the value
// read goes to a volatile sink so the load is not dead.
[[noreturn]] void read_past_the_end_and_exit() {
    // Catch2 installs handlers for SIGABRT (and SIGILL, SIGSEGV, ...) in the
    // parent, and the child inherits them. libstdc++'s check aborts, so on
    // Linux the child would print Catch2's "fatal error" report for a case the
    // parent then passes. Restore the defaults: the child just dies by the
    // signal, which is what the parent inspects.
    std::signal(SIGABRT, SIG_DFL);
    std::signal(SIGILL, SIG_DFL);
    std::signal(SIGTRAP, SIG_DFL);
    std::vector<int> v(4, 7);
    volatile std::size_t index = 4;
    volatile int sink = v[index];
    (void)sink;
    _exit(0);  // reached only if the library did not check
}

bool killed_by_a_library_check(int status) {
    if (!WIFSIGNALED(status)) {
        return false;
    }
    const int sig = WTERMSIG(status);
    // SIGTRAP: libc++'s __builtin_trap on arm64; SIGILL: the same on x86-64;
    // SIGABRT: libstdc++'s __glibcxx_assert_fail.
    return sig == SIGTRAP || sig == SIGILL || sig == SIGABRT;
}

}  // namespace

TEST_CASE("T1: the reported mode is a compile-time constant that agrees with the option",
          "[hardening]") {
    // The static_assert above is the guard; this case makes it visible in the
    // ctest listing and states the expectation at run time as well.
    constexpr std::string_view mode = terrain::stdlib_hardening();
    STATIC_REQUIRE(!mode.empty());
    CHECK((mode != "none") == (RASPUTIN_EXPECT_HARDENING != 0));
}

TEST_CASE("T2: an out-of-range vector read stops the process when hardening is on",
          "[hardening]") {
    if (RASPUTIN_EXPECT_HARDENING == 0) {
        SKIP("RASPUTIN_HARDENING is OFF: the read would be undefined behaviour, not a check");
    }
    const pid_t pid = fork();
    REQUIRE(pid >= 0);
    if (pid == 0) {
        read_past_the_end_and_exit();
    }
    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);
    INFO("child " << (WIFSIGNALED(status)
                          ? "killed by signal " + std::to_string(WTERMSIG(status))
                          : "exited with status " + std::to_string(WEXITSTATUS(status))));
    CHECK(killed_by_a_library_check(status));
}

TEST_CASE("T3: stdlib_hardening names the library and mode this build uses", "[hardening]") {
    constexpr std::string_view mode = terrain::stdlib_hardening();
#if RASPUTIN_EXPECT_HARDENING
#if defined(_LIBCPP_VERSION)
    CHECK(mode == "libc++ fast");
#elif defined(__GLIBCXX__)
    CHECK(mode == "libstdc++ assertions");
#else
    FAIL("a standard library this suite does not know; extend T3");
#endif
#else
    CHECK(mode == "none");
#endif
}

TEST_CASE("T3: stdlib_hardening is noexcept", "[hardening]") {
    STATIC_REQUIRE(noexcept(terrain::stdlib_hardening()));
}
