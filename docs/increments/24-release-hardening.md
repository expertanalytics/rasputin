# Increment 24: release hardening, on by default, off for heavy runs

Status: **designed and ruled**, `@architect`, 2026-10-03. Ola chose option (a)
and answered Q1-Q4 the same day (§11). No code yet. One PR, about 55
production lines. The measurement it rests on is
`docs/benchmarks/2026-10-03/release-hardening/README.md` (`@perf`, #149).

**The ask** (Ola, 2026-10-03): "I think we can live with hardening, but for
certain runs, we might want to turn it off for speed. So normal release, turn
it on. For real heavy lifting, option to turn it off."

**What "hardening" means here.** Two predefined macros, one per standard
library, that make the library check its own preconditions in Release builds:
`_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST` for libc++ (macOS) and
`_GLIBCXX_ASSERTIONS` for libstdc++ (Linux). An out-of-range
`std::vector::operator[]` or `std::span` index then stops the process at the
faulting call (a trap or an abort) instead of reading or writing adjacent
memory and carrying on. It is still a crash of the process, not a Python
exception. Measured cost: refine about 5 % (4.4 to 8.7 % across thread counts
and both libraries), whole CLI process 2 to 3 %, meshes byte-identical.

**The design constraint.** The macros are compile-time. "Off" is a different
build of the C++ extension, never a run-time switch inside one build. The
ruled design makes that second build something the user makes on purpose
(§3); the alternative of shipping both builds in every wheel is kept in §4 as
the rejected option, with its costs.

## 1. Prior art: legacy and literature

No algorithm and no novelty claim. The practice this design follows:

- **Hardened by default in production is now mainstream for C++.** Google
  switched hardened libc++ on across its server fleet and reports 0.30 %
  average overhead, over 1,000 bugs found and 30 % fewer segmentation faults
  ([Google Security Blog, 2024-11, "Retrofitting spatial safety to hundreds
  of millions of lines of C++"](https://security.googleblog.com/2024/11/retrofitting-spatial-safety-to-hundreds.html)).
  Our 5 % on refine is far above their fleet average; refine is a tight
  indexed loop over contiguous arrays, which is where bounds checks cost most.
- **C++26 standardises the idea** as "hardened implementations" (P3471,
  adopted February 2025, [wg21.link/P3471](https://wg21.link/P3471)). How a
  hardened implementation is selected is implementation-defined, in practice
  a build flag, as here.
- **Linux distributions build with `_GLIBCXX_ASSERTIONS` by default.**
  Fedora's `redhat-rpm-config` adds `-Wp,-D_GLIBCXX_ASSERTIONS` to every
  package build ([buildflags.md](https://src.fedoraproject.org/rpms/redhat-rpm-config/blob/8d6c6d0761089dfb7e1eacf35112b3d7a3f0d1af/f/buildflags.md)),
  documented there as "no impact on class layout", suitable for release.
- **Numeric libraries: checks on, opt out at compile time.** Armadillo bounds-
  checks by default and is made a "non-development" build by defining
  `ARMA_NO_DEBUG` ([RcppArmadillo release note](https://www.r-bloggers.com/2012/04/rcpparmadillo-0-3-0-2-released-and-on-cran/)).
  That is the ruled design's shape exactly. Eigen goes the other way (checks
  tied to `NDEBUG`, so off in Release). None of the libraries looked at ships
  a checked and an unchecked extension side by side in one Python wheel;
  Python projects that need build variants publish separate distributions
  (`opencv-python` / `opencv-python-headless`), and a standard for variant
  wheels is still under discussion
  ([PEP 817](https://discuss.python.org/t/pep-817-wheel-variants-beyond-platform-tags/105860),
  split into [PEP 825](https://discuss.python.org/t/pep-825-wheel-variants-package-format-split-from-pep-817/106196)),
  not usable today.
- **ABI.** libc++: "Setting a hardening mode does not affect the ABI";
  functions are ABI-tagged per mode so translation units built in different
  modes link without ODR violations ([libc++ Hardening Modes](https://libcxx.llvm.org/Hardening.html)).
  libstdc++'s `_GLIBCXX_ASSERTIONS` does not change layout (unlike
  `_GLIBCXX_DEBUG`). So Catch2's own translation units, which are built
  outside `terrain_headers` and stay unchecked, link safely with checked test
  code — the open question in `@perf`'s README, answered by the documentation.

*Legacy.* Nothing; the legacy tree never set a hardening or assertion macro:

```
$ git grep -l -iE 'hardening|GLIBCXX_ASSERTIONS|_LIBCPP_|GLIBCXX_DEBUG|bounds.?check' legacy-archive -- legacy
(no output, exit 1)
```

## 2. Facts checked for this design

Each was run on 2026-10-03, macOS arm64, AppleClang 21, CMake 4.4.3.

- **`WILL_FAIL` does not catch a crash.** `@perf`'s README suggests a guard
  test "registered `WILL_FAIL`". CTest does not invert a test killed by a
  signal: a one-past-the-end `v[4]` on a `std::vector<int>(4)`, built with the
  FAST define and registered `WILL_FAIL TRUE`, is reported `SIGTRAP***Exception`
  and counted failed; the same program without the define reads past the end,
  exits 0 and is counted failed by the inversion. So under `WILL_FAIL` the
  hardened build fails and the guard is unusable in that form. Probe: a
  two-target CMake project (one `oob.cpp` as described, built with and without
  the define, both `set_tests_properties(... WILL_FAIL TRUE)`); `ctest` output
  "`1 - oob (SIGTRAP)`, `2 - oob_plain (Failed)`". The guard in §5 forks
  instead.
- **This machine's libc++ knows the macro**: `_LIBCPP_VERSION 210106`, and
  `<vector>` defines `_LIBCPP_HARDENING_MODE_FAST (1 << 2)` and an unset mode
  as `_LIBCPP_HARDENING_MODE_DEFAULT` (`clang++ -dM -E` on `#include <vector>`).
  A libc++ older than 18 does not define `_LIBCPP_HARDENING_MODE_FAST`, which
  is what the compile-time guard tests for.
- **The extension today is 720,984 bytes** (`build-pyext/_core.cpython-314-darwin.so`
  in the main checkout). A second copy was option (b)'s wheel cost.
- **Six modules import `tin_engine._core` at module import time**
  (`__init__.py`, `cli.py`, `catchment.py`, `final_check.py`, `raster.py`, and
  `tools/bench.py`'s child;
  `grep -ln "from tin_engine._core import\|import tin_engine._core" src_python/tin_engine/*.py src_python/tin_engine/*/*.py tools/bench.py`
  returns exactly those six files). `cli.py` does so on line 60, before Typer
  parses any option. That is why option (b) could not offer a CLI flag (§4).

## 3. The design (option (a), ruled): one CMake option, default ON

**What it is.** `RASPUTIN_HARDENING`, a CMake option, default ON. ON adds the
two macros to `terrain_headers`' INTERFACE, so they reach `terrain_predicates`,
`terrain_cdt`, `_core` and every ctest suite. OFF builds exactly today's
binary. A normal `pip install` is hardened; a heavy run uses an unchecked
build that the user makes on purpose.

```cmake
option(RASPUTIN_HARDENING "Bounds-checked standard library in every build type" ON)
if(RASPUTIN_HARDENING)
    target_compile_definitions(terrain_headers INTERFACE
        _LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST _GLIBCXX_ASSERTIONS)
endif()
```

A `-D_LIBCPP_HARDENING_MODE=...` passed in `CMAKE_CXX_FLAGS` on top of the
option ON is a macro redefinition and fails the `-Werror` build. That is
intended: to use another mode (extensive, debug), set the option OFF and pass
the mode in the flags. INSTALL says so.

**How to turn it off.**

```sh
# a wheel or a non-editable install
uv pip install ".[codecs]" -C cmake.define.RASPUTIN_HARDENING=OFF
python -m pip install ".[codecs]" -C cmake.define.RASPUTIN_HARDENING=OFF
# plain CMake (C++ tests, or the build-pyext route in REQUIRED-READING)
cmake -S . -B build-unchecked -DCMAKE_BUILD_TYPE=Release -DRASPUTIN_HARDENING=OFF
```

The documented heavy-run recipe is a **second virtual environment**
(`.venv-unchecked`) with a non-editable OFF install, so the everyday
environment never silently changes mode. Swapping one `.so` in and out of a
single venv is the error the reporting below exists to catch, not a workflow
to recommend.

**A stale-cache trap, and its fix.** scikit-build-core reuses
`build/{wheel_tag}` (`pyproject.toml`, `build-dir`). If an OFF install leaves
`RASPUTIN_HARDENING=OFF` in that directory's `CMakeCache.txt`, the next plain
`pip install .` would reconfigure with the cached OFF, because nothing passes
the option. So the default is passed explicitly on every build, from
`pyproject.toml`:

```toml
[tool.scikit-build.cmake.define]
RASPUTIN_HARDENING = "ON"
```

and `-C cmake.define.RASPUTIN_HARDENING=OFF` overrides that key. That the
command-line setting overrides the `pyproject.toml` one is scikit-build-core's
documented precedence, not checked here; T6 in §12 checks it. `bench.py`
likewise passes the value explicitly on every configure (§7).

**How a user knows which build they have** (Q2 and Q4, ruled). From the
module that is actually loaded, never from a cache file or an option value:

- `tin_engine._core.hardening`, a string set at compile time from what the
  standard library itself defines: `"libc++ fast"`, `"libc++ extensive"`,
  `"libc++ debug"`, `"libstdc++ assertions"`, or `"none"`. It is computed by a
  `constexpr` function in a new, pybind11-free header,
  `include/terrain/build_info.hpp`, which the binding and the guard test both
  include (so the attribute and the guard read one object):

  ```cpp
  #include <string_view>
  #include <version>  // defines _LIBCPP_VERSION / __GLIBCXX__ and the libc++ mode constants
  namespace terrain {
  [[nodiscard]] constexpr std::string_view stdlib_hardening() noexcept {
  #if defined(_LIBCPP_HARDENING_MODE_FAST) && _LIBCPP_HARDENING_MODE == _LIBCPP_HARDENING_MODE_FAST
      return "libc++ fast";
  #elif defined(_LIBCPP_HARDENING_MODE_EXTENSIVE) && _LIBCPP_HARDENING_MODE == _LIBCPP_HARDENING_MODE_EXTENSIVE
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
  ```

  The `defined(..._FAST)` half is load-bearing: a libc++ before 18 does not
  define the mode constants, so the branch is not taken and the function
  reports `"none"`, which the guard test then refuses (§5).
- `rasputin version` prints two lines: the bare version, unchanged, then
  `bounds checks: on (<mode>)`, e.g. `bounds checks: on (libc++ fast)`, or
  `bounds checks: off` when the mode is `"none"`. Scripts that read line one
  are unaffected.
- `rasputin mesh --stats` adds the same `bounds checks: ...` text to the
  report, as one line before its first `##` section, so every stats file from
  a heavy run says it was unchecked. `stats.Report` gets a field
  `bounds_checks: str` (the text after `bounds checks: `), filled by
  `cli.py` from `_core.hardening`; `stats.py` itself stays free of `_core`.

**Costs.** About 55 production lines (§8). Build time and wheel size
unchanged; an OFF build is one more build, on demand. One more C++ CI leg (§6).
What a heavy run costs Ola: a one-off second venv, rebuilt after C++ changes.
The residual risk, a heavy run done in the wrong venv, is caught by
`rasputin version` and the `--stats` line, not prevented.

## 4. Rejected: option (b), two extension modules in one wheel

Kept as the alternative Ola did not choose (Q1), with what it would have
taken, so a later proposal starts from here rather than from scratch. It can
be added on top of the ruled design without undoing any of it.

Every build would compile the C++ twice, `_core_checked` and
`_core_unchecked`, both shipped in the wheel, a Python shim
`tin_engine/_core.py` re-exporting one of them, selected by
`RASPUTIN_UNCHECKED=1` at first import. What that required:

- hardening off `terrain_headers`' INTERFACE, the compiled libraries built
  twice by a CMake function (about 35 lines in place of 6);
- `PYBIND11_MODULE(RASPUTIN_CORE_MODULE, m)` with the module name as a
  per-target definition, relying on pybind11's macro indirection;
- every class's `__module__` naming the variant, so pickles do not cross
  variants;
- never both modules in one process (pybind11's shared type registry refuses
  the second), so variant tests in subprocesses;
- no `--unchecked` CLI flag: `cli.py` imports `_core` before Typer parses
  (§2), so a flag means re-exec'ing the process or making six imports lazy;
- an unchecked variant shipped in every wheel with thinner tests than the
  default.

Costs: about 125 production lines; the extension compiled twice in every
install, editable rebuild and CI Python leg; 0.7 MB more per wheel. What it
bought: a run-time switch worth 2 to 3 % of a whole CLI run.

## 5. The guard test

`@perf` measured that the macro works on today's toolchains; the guard makes
CI fail when a toolchain stops honouring it (older libc++ ignores the unknown
macro silently). One Catch2 suite, `tests/cpp/unit/test_build_hardening.cpp`,
target `test_build_hardening`, registered with `add_terrain_test` and given
`RASPUTIN_EXPECT_HARDENING=$<BOOL:${RASPUTIN_HARDENING}>` (expanding to 1 or
0) by `tests/cpp/CMakeLists.txt`.

1. **Compile time.**
   `static_assert((terrain::stdlib_hardening() != "none") == (RASPUTIN_EXPECT_HARDENING != 0));`
   at namespace scope. ON with an old libc++ fails to compile; OFF with a
   stray definition fails too. This is the cheap guard and catches the
   silent-ignore case. A `#ifndef RASPUTIN_EXPECT_HARDENING #error` above it
   stops the suite compiling if the CMake line is lost.
2. **Behaviour, ON only.** The test `fork()`s; the child indexes one past the
   end of a `std::vector<int>(4)` through a `volatile std::size_t` index, uses
   the value (so it is not optimised away) and `_exit(0)`s if it survives; the
   parent `waitpid`s and requires `WIFSIGNALED` with `SIGTRAP`, `SIGILL`
   (libc++'s trap on x86-64) or `SIGABRT` (libstdc++). This proves the check
   fires, not just that the macro is set. With OFF the case `SKIP`s rather
   than commit undefined behaviour. `fork` rather than `WILL_FAIL`, because of
   §2. POSIX only, which both CI platforms are.

Planting (`.claude/REQUIRED-READING.md`, "Claims"): `@tester` shows that
(1) fails when the `target_compile_definitions` is removed with the option ON,
and that (2) fails when the child's index is made in-range. Whether the forked
child's crash puts a crash-reporter dialog on Ola's Mac (as a Python crash
does, per `@perf`'s README) is checked once by running the suite locally; a
crash report file under `~/Library/Logs/DiagnosticReports` is expected and
harmless.

## 6. CI (Q3, ruled)

- **`cpp`, Ubuntu and macOS legs**: unchanged commands; default ON, so they
  build hardened and run the guard. The macOS leg is where an old Xcode libc++
  would be caught.
- **New leg: `cpp` on Ubuntu with `-DRASPUTIN_HARDENING=OFF`**, a matrix
  `include` entry on the `cpp` job (about 12 lines of YAML), building and
  running every ctest suite. The heavy-run build is one nobody else compiles;
  this keeps it compiling and runs the guard's OFF half. Its name says
  "unchecked" so a red leg is readable.
- **`sanitizers` (ASan+UBSan, Debug)**: unchanged; it builds with the default,
  so hardening is on under the sanitizers. The libraries document no conflict.
  An out-of-range container index now stops at the library check (an abort
  naming the assertion) before ASan reports it; the job still goes red, and
  ASan still covers raw pointers and `new[]` arrays the library checks never
  see. Turning hardening off there would lose the in-capacity overruns ASan
  cannot see without container annotations, so it stays on.
- **`tsan`**: unchanged; it builds a named target list, which does not include
  the guard.
- **`python` legs**: `pip install -e` builds with the default, so the whole
  Python suite runs hardened, which is what ships.

## 7. `tools/bench.py` and the 5 % line

- **`bench.py run --hardening on|off`, default `on`.** `build()` passes
  `-DRASPUTIN_HARDENING=ON|OFF` on every configure, so a cached value never
  decides. (On a tree from before this increment the define is an unused
  variable, CMake warns, and the build is unchecked, which is correct.) With
  `--no-build` the flag is refused when given explicitly, since nothing is
  built; the mode is still read from the loaded module (next point).
- **The record says what was measured, read from the loaded module**: the
  child adds `"hardening": getattr(tin_engine._core, "hardening", "none")` to
  its `BENCH` JSON line; the parent requires every sample of a run to report
  the same value (else a `BenchError`, exit 3), and stores it in a new
  `RunRecord.hardening: str = "none"`. A stored `run.json` without the field
  loads as `"none"`: unchecked, which every run before this increment was as
  `bench.py` builds it. (It would be wrong only for a tree hand-built with the
  macro in a cached `CMAKE_CXX_FLAGS`, as `@perf`'s study did deliberately.)
  When `--hardening` was given, the reported mode must agree with it (`on` ⇔
  not `"none"`), else `BenchError`: the flag and the measured object cannot
  silently differ.
- **`comparable()` refuses across modes**: one more pair, `("hardening",
  a.hardening, b.hardening)`, so the refusal reads `hardening: libc++ fast vs
  none`, as power does. A hardened run is never judged against an unchecked
  baseline, so the 5 % threshold stays 5 % and never absorbs the hardening
  cost.
- **Hardened is the baseline from this increment on.** Its acceptance run is
  the one exception to "compare with the previous increment": `bench.py run
  --hardening off` against the previous merge commit must show no regression
  (the PR changed nothing else), and the `--hardening on` run of the same
  tree, measured back to back, is stored as the first hardened baseline. The
  ON/OFF difference of that pair is recorded as the hardening cost, and should
  reproduce `@perf`'s 4 to 9 %. Later increments measure ON against ON as
  usual. `--hardening off` runs are judged only against OFF runs.
- **The generated run README** prints the mode on its build line.
- `.claude/agents/perf.md` and `docs/benchmarks/bench-py.md` get one sentence
  each saying the above. The first is a rule file: its edit goes through the
  governance guard in the same PR.

## 8. LOC

Production lines under `CLAUDE.md` §2. Tests, YAML and prose not counted.

| item | lines |
|---|---:|
| `CMakeLists.txt` (option, definitions) | 6 |
| `pyproject.toml` (`cmake.define` default) | 2 |
| `include/terrain/build_info.hpp` | 18 |
| `bindings/core.cpp` (`m.attr("hardening")`, include) | 2 |
| `_core.pyi` | 1 |
| `cli.py` (`version` line; `--stats` field) | 6 |
| `stats.py` (one report field and line) | 4 |
| `tools/bench.py` (flag, configure argv, child field, record, agreement checks, `comparable` pair, README line) | 18 |
| **total** | **about 57** |
| tests (not counted) | guard suite ~50, Python ~60 |
| CI YAML (not counted) | ~12 |

One PR, well under the ceiling. The PR touches no refine or mesh source, but
it changes every refine binary, so `@perf`'s acceptance is required all the
same, in the paired form of §7.

## 9. Also rejected

- **One extension with run-time dispatch.** The checks live in inline template
  instantiations throughout the core; switching at run time means compiling
  every hot template twice under different names and dispatching at every
  entry point. That is option (b) inside one `.so`, with more code.
- **Two distributions (`rasputin` and `rasputin-unchecked`).** The precedent
  for Python build variants, and what variant wheels would standardise, but it
  needs a publishing pipeline the project does not have; rasputin is installed
  from source today.
- **`@perf`'s "middle option": checked in CI, unchecked shipped.** Splits what
  is tested from what is shipped; superseded by Ola's ruling (on by default).
- **GCC 14's `-fhardened`.** Bundles `_GLIBCXX_ASSERTIONS` with stack
  protector, FORTIFY and PIE flags whose cost was not measured, and has no
  AppleClang counterpart; the two macros are what `@perf` measured.

## 10. Why (a)

Option (a) is the shape the practice in §1 converges on (checks on, opt out
at build time), it is half the code of (b), it leaves build time, wheel size
and the import path alone, and every build that exists is one CI compiles.
What it costs Ola is a one-off second environment for heavy runs and a
rebuild of it when the C++ changes.

## 11. Rulings (Ola, 2026-10-03)

- **Q1, ruled (a)**: the CMake option `RASPUTIN_HARDENING`, default ON, off
  by rebuilding. Option (b) is rejected (§4).
- **Q2, ruled yes to the names**: `RASPUTIN_HARDENING` (CMake),
  `_core.hardening`, and `rasputin version`'s second line,
  `bounds checks: on (libc++ fast)`.
- **Q3, ruled yes**: the Ubuntu CI leg that builds OFF (§6).
- **Q4, ruled yes**: a line naming the mode in the `--stats` report (§3).

## 12. For `@tester`: the red suite

Every item below is a test that fails before the implementation and passes
after it. No mutation round: nothing here is invariant-critical.

**C++, `tests/cpp/unit/test_build_hardening.cpp`** (§5):

- **T1**, compile-time guard: the `static_assert` and the `#error` exactly as
  §5.1. Red because `include/terrain/build_info.hpp` does not exist (the
  suite does not compile), which is the expected red for a new header.
- **T2**, behavioural guard: the forked out-of-range read of §5.2, with the
  three accepted signals, `SKIP` when `RASPUTIN_EXPECT_HARDENING` is 0.
- **T3**, the reporting function's values: on the build in hand,
  `stdlib_hardening()` is `"libc++ fast"` when `_LIBCPP_VERSION` is defined,
  `"libstdc++ assertions"` when `__GLIBCXX__` is, under ON; `"none"` under OFF.
- Registration in `tests/cpp/CMakeLists.txt` with the
  `RASPUTIN_EXPECT_HARDENING` definition. Planting, recorded in the commit
  message: T1 fails with the `target_compile_definitions` removed and the
  option ON; T2 fails with the child's index in range.

**Python**:

- **T4**, `tin_engine._core.hardening` on the installed (default) build is
  `"libc++ fast"` on macOS and `"libstdc++ assertions"` on Linux (keyed on
  `sys.platform`).
- **T5**, `rasputin version` (`CliRunner`): exit 0, two lines; line one is
  `installed_version()` exactly; line two is `bounds checks: on (<mode>)` with
  `<mode>` equal to `_core.hardening`. With `_core.hardening` monkeypatched to
  `"none"` (patch the name `cli.py` reads), line two is `bounds checks: off`.
- **T6**, the stale build-dir trap (§3). A slow test, marked and run in one
  CI leg only, or a manual check recorded in the review if `@tester` finds it
  too slow for CI (it takes two real non-editable installs into throwaway
  venvs): install with `-C cmake.define.RASPUTIN_HARDENING=OFF` into venv A,
  then into venv B with no setting, from the same source tree and so the same
  `build/{wheel_tag}`; assert A's `_core.hardening == "none"` and B's is the
  platform's hardened mode. It also checks the override precedence §3 relies
  on. `@tester` decides between CI and manual and says which in the red
  commit message.
- **T7**, `--stats` (on an existing small fixture run of `rasputin mesh
  --stats -`): the report contains exactly one `bounds checks: ...` line,
  before the first `##` heading, equal to the `version` line two. A
  `stats.render` unit test with `bounds_checks="off"` and with
  `"on (libc++ fast)"` pins the text without `_core`.
- **T8**, `tools/bench.py`, in `tests/python/test_bench.py`, through its injected `Runner`:
  - `build()`'s configure argv contains `-DRASPUTIN_HARDENING=ON` for
    `--hardening on` and `=OFF` for `off`, and exactly one such define;
  - a `run.json` without `hardening` loads with `hardening == "none"`;
  - `comparable()` returns a reason starting `hardening:` for two records that
    differ only in mode, and `None` when they agree;
  - `find_baseline()` skips a stored run of the other mode;
  - children whose `BENCH` lines report different modes make `run` exit 3;
  - a `BENCH` line without `hardening` (an old tree's `_core` has no
    attribute, so the child reports `"none"`): the parsed child's mode is
    `"none"`;
  - `--hardening on` with a child reporting `"none"` exits 3, and
    `--hardening off` with a child reporting a hardened mode exits 3.

**CI** (§6): the Ubuntu OFF leg is `@developer`'s YAML; `@reviewer` checks it
ran the guard (T1 OFF half compiled, T2 skipped) from the job log.

## 13. Documentation in the same PR

- `INSTALL.md`: a section "Bounds checks (on by default)": what they do, the
  measured cost with a link to `@perf`'s README, the two `-C cmake.define`
  commands and the plain-CMake line, the second-venv recipe for heavy runs,
  `rasputin version` to tell which build is installed, and that a mode other
  than FAST is chosen by setting the option OFF and passing the mode in
  `CMAKE_CXX_FLAGS`. Its "C++ changes have no effect" note also applies to the
  unchecked venv.
- `README.md`: one sentence where installation is mentioned, pointing at that
  section.
- `ROADMAP.md`: row 24 updated to shipped by the merge; the measurement row
  already points here.
