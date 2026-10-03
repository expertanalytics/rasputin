# Increment 24: release hardening, on by default, off for heavy runs

Status: **design, two options for Ola to rule**, `@architect`, 2026-10-03.
No code. One PR whichever option is chosen. The measurement it rests on is
`docs/benchmarks/2026-10-03/release-hardening/README.md` (`@perf`, PR #149;
on branch `worktree-agent-a0ff1ca8678bdc0c3` until that merges).

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
build of the C++ extension, never a run-time switch inside one build. The two
options below differ only in where that second build lives: built on demand
by the user (a), or shipped beside the default in every wheel (b).

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
  That is option (a)'s shape exactly. Eigen goes the other way (checks tied to
  `NDEBUG`, so off in Release). None of the libraries looked at ships a checked
  and an unchecked extension side by side in one Python wheel; Python projects
  that need build variants publish separate distributions
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
  two-target CMake project in the session scratchpad (`oob.cpp` above,
  `set_tests_properties(oob PROPERTIES WILL_FAIL TRUE)`), `ctest` output
  "`1 - oob (SIGTRAP)`, `2 - oob_plain (Failed)`". The guard in §5 forks
  instead.
- **This machine's libc++ knows the macro**: `_LIBCPP_VERSION 210106`, and
  `<vector>` defines `_LIBCPP_HARDENING_MODE_FAST (1 << 2)` and an unset mode
  as `_LIBCPP_HARDENING_MODE_DEFAULT` (`clang++ -dM -E` on `#include <vector>`).
  A libc++ older than 18 does not define `_LIBCPP_HARDENING_MODE_FAST`, which
  is what the compile-time guard tests for.
- **The extension today is 720,984 bytes** (`build-pyext/_core.cpython-314-darwin.so`
  in the main checkout). A second copy is option (b)'s wheel cost.
- **Six modules import `tin_engine._core` at module import time**
  (`__init__.py`, `cli.py`, `catchment.py`, `final_check.py`, `raster.py`, and
  `tools/bench.py`'s child;
  `grep -ln "from tin_engine._core import\|import tin_engine._core" src_python/tin_engine/*.py src_python/tin_engine/*/*.py tools/bench.py`
  returns exactly those six files).
  `cli.py` does so on line 60, before Typer parses any option. That decides
  how option (b) can select a variant (§4).

## 3. Option (a): one CMake option, default ON

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

and `-C cmake.define.RASPUTIN_HARDENING=OFF` overrides that key. `@tester`
plants the trap: install OFF, then install with no setting, and assert
`_core.hardening` reports the hardened mode. The same holds for `bench.py`,
which passes the value explicitly on every configure (§7).

**How a user knows which build they have.** From the module that is actually
loaded, never from a cache file or an option value:

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
  reports `"none"` — which the guard test then refuses (§5).
- `rasputin version` prints a second line: `bounds checks: on (libc++ fast)`
  or `bounds checks: off`. The first line stays the bare version (the
  existing test, `tests/python/test_cli.py::test_version_reports_the_installed_distribution`,
  only asserts non-empty output; scripts that read line one are unaffected).
- `rasputin mesh --stats` adds one line to the report's header, the same
  string, so every stats file from a heavy run says it was unchecked.
  Recommended, not essential (about 5 lines); see Q4.

**Costs.**

| | |
|---|---|
| production LOC | about 55 (table in §8) |
| build time | unchanged for the default; an OFF build is one more build, on demand |
| wheel size | unchanged |
| CI | one more C++ leg (OFF, Ubuntu), a few minutes (§6) |
| what a heavy run costs Ola | a one-off second venv, then `source .venv-unchecked/bin/activate`; a rebuild of that venv after C++ changes |
| risk | a heavy run done in the wrong venv; caught by `rasputin version` and the `--stats` line, not prevented |

## 4. Option (b): two extension modules in one wheel

**What it is.** Every build compiles the C++ twice: `_core_checked` (the
macros on) and `_core_unchecked` (off). Both ship in the wheel. A small Python
module, `tin_engine/_core.py`, imports one of them and re-exports its names, so
the six importers above are unchanged. A run selects the variant with an
environment variable, `RASPUTIN_UNCHECKED=1 rasputin mesh ...`.

**What it takes.**

- **CMake.** Hardening can no longer ride on `terrain_headers`' INTERFACE,
  since one tree must produce both modes. `terrain_headers` stays mode-
  neutral; a function builds the compiled set twice
  (`terrain_predicates_<v>`, `terrain_cdt_<v>`, then the module), each copy
  with its own definitions, and installs both modules. The ctest suites take
  the checked definitions. About 35 lines in place of 6.
- **pybind11 module naming.** `PYBIND11_MODULE(name, m)` needs the literal
  module name, which must match the file name. One `bindings/core.cpp`,
  compiled twice with `RASPUTIN_CORE_MODULE=_core_checked` /
  `_core_unchecked`, and `PYBIND11_MODULE(RASPUTIN_CORE_MODULE, m)`; this
  relies on pybind11's internal macro indirection expanding the argument
  before token pasting, which `@developer` must confirm on both compilers.
- **The shim.** `_core.py` beside the existing `_core.pyi` (mypy reads the
  stub, so typing is unchanged). It reads the variable once, at first import,
  then `globals().update(...)` from the chosen module. Every class's
  `__module__` becomes `tin_engine._core_checked` or `..._unchecked`, so a
  pickle written under one variant does not load under the other.
- **Never both in one process.** The two modules register the same C++ types
  with pybind11's shared registry; importing the second fails
  ("type already registered"). The shim guarantees one per process, but the
  tests that compare variants must run each in a subprocess.
- **No `--unchecked` CLI flag without a re-exec.** `cli.py` imports `_core`
  before Typer sees any option, so a flag can only take effect by setting the
  variable and `os.execv`-ing itself, which breaks under `CliRunner` and is
  invisible to the tests. Making every `_core` import lazy instead touches six
  modules. The design therefore offers the environment variable only.
- **Tests.** The Python suite runs against the checked module as today; the
  unchecked one needs at least a subprocess smoke test of a mesh (and arguably
  the golden refine tests twice), otherwise it ships untested.

**Costs.**

| | |
|---|---|
| production LOC | about 125 (§8) |
| build time | the compiled C++ (three translation units, `bindings/core.cpp` the heavy one) builds twice in every `pip install`, every editable rebuild and every CI Python leg: roughly double the extension's compile time (not measured) |
| wheel size | one more extension, 0.7 MB uncompressed on macOS today |
| CI | longer Python legs (the double build), plus the unchecked smoke tests |
| what a heavy run costs Ola | nothing: one environment variable |
| risk | a variant that ships in every wheel with thinner tests than the default; the shim is a new place for import-order bugs |

## 5. The guard test (both options)

`@perf` measured that the macro works on today's toolchains; the guard makes
CI fail when a toolchain stops honouring it (older libc++ ignores the unknown
macro silently). One Catch2 suite, `tests/cpp/unit/test_build_hardening.cpp`,
target `test_build_hardening`, given
`RASPUTIN_EXPECT_HARDENING=$<BOOL:${RASPUTIN_HARDENING}>` by
`tests/cpp/CMakeLists.txt` (in option (b), always 1: the suites are checked).

1. **Compile time.** `static_assert((terrain::stdlib_hardening() != "none") == RASPUTIN_EXPECT_HARDENING)`.
   ON with an old libc++ fails to compile; OFF with a stray definition fails
   too. This is the cheap guard and catches the silent-ignore case.
2. **Behaviour, ON only.** The test `fork()`s; the child indexes one past the
   end of a `std::vector<int>(4)` through a `volatile` index and `_exit(0)`s
   if it survives; the parent `waitpid`s and requires `WIFSIGNALED` with
   `SIGTRAP`, `SIGILL` (libc++'s trap on x86-64) or `SIGABRT` (libstdc++).
   This proves the check fires, not just that the macro is set. With OFF the
   case `SKIP`s rather than commit undefined behaviour. `fork` rather than
   `WILL_FAIL`, because of §2. POSIX only, which both CI platforms are.

Planting (`.claude/REQUIRED-READING.md`, "Claims"): `@tester` shows that
(1) fails when the `target_compile_definitions` is removed with the option ON,
and that (2) fails when the child's index is made in-range. Whether the
forked child's crash puts a crash-reporter dialog on Ola's Mac (as a Python
crash does, per `@perf`'s README) is checked once by running the suite
locally; a crash report file under `~/Library/Logs/DiagnosticReports` is
expected and harmless.

## 6. CI

- **`cpp`, Ubuntu and macOS legs**: unchanged commands; default ON, so they
  build hardened and run the guard. The macOS leg is where an old Xcode libc++
  would be caught.
- **New leg, `cpp` on Ubuntu with `-DRASPUTIN_HARDENING=OFF`** (option (a);
  a matrix `include` entry, about 12 lines of YAML). The heavy-run build is a
  build nobody else compiles; this keeps it compiling and runs the guard's OFF
  half. In option (b) the wheel build compiles both, so this leg is not needed.
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
  variable, CMake warns, and the build is unchecked, which is correct.)
- **The record says what was measured, read from the loaded module**: the
  child adds `getattr(tin_engine._core, "hardening", "none")` to its `BENCH`
  line, the parent requires every sample to agree and stores it in
  `RunRecord` (a new field, default `"none"` so every stored `run.json`
  without it loads as what it was: unchecked). The absent-attribute default is
  right for every tree before this increment as `bench.py` builds it; it
  would be wrong only for a tree hand-built with the macro in a cached
  `CMAKE_CXX_FLAGS`, as `@perf`'s study did deliberately.
- **`comparable()` refuses across modes** (one more pair, `hardening`), the
  same way it refuses AC against battery. A hardened run is never judged
  against an unchecked baseline, so the 5 % threshold stays 5 % and never
  absorbs the hardening cost.
- **Hardened is the baseline from this increment on.** Its acceptance run is
  the one exception to "compare with the previous increment": `bench.py run`
  with `--hardening off` against the previous merge commit must show no
  regression (the PR changed nothing else), and the `--hardening on` run of
  the same tree, measured back to back, is stored as the first hardened
  baseline. The ON/OFF difference of that pair is recorded as the hardening
  cost, and should reproduce `@perf`'s 4 to 9 %. Later increments measure ON
  against ON as usual. `--hardening off` runs are judged only against OFF
  runs, for when a heavy-run speed question is asked.
- In option (b) the same flag selects the variant through the environment
  variable instead of a rebuild, and `build()` must expect two `.so` files
  where it now refuses anything but one.
- `.claude/agents/perf.md` and `docs/benchmarks/bench-py.md` get one sentence
  each saying the above. The first is a rule file: its edit goes through the
  governance guard in the same PR.

## 8. LOC, per option

Production lines under `CLAUDE.md` §2. Tests, YAML and prose not counted.

| item | (a) | (b) |
|---|---:|---|
| `CMakeLists.txt` | 6 | 35 |
| `pyproject.toml` (`cmake.define` default) | 2 | 2 |
| `include/terrain/build_info.hpp` | 18 | 18 |
| `bindings/core.cpp` (attribute; module-name macro in (b)) | 2 | 4 |
| `_core.pyi` | 1 | 1 |
| `src_python/tin_engine/_core.py` (the shim) | — | 20 |
| `cli.py` (`version` line; `--stats` header) | 6 | 6 |
| `stats.py` (one report field and line) | 4 | 4 |
| `tools/bench.py` | 16 | 35 |
| **total** | **about 55** | **about 125** |
| tests (not counted) | guard suite ~50, Python ~40 | guard suite ~50, Python ~90 (subprocess variants) |
| CI YAML | ~12 (OFF leg) | ~10 (unchecked smoke step) |

Either fits one PR well under the ceiling.

## 9. Rejected

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

## 10. Recommendation: option (a)

Option (a) is the shape the practice in §1 converges on (checks on, opt out
at build time), it is half the code, it leaves build time, wheel size and the
import path alone, and every build that exists is a build CI compiles. What
it costs Ola is a one-off second environment for heavy runs and a rebuild of
it when the C++ changes.

Option (b) buys a free switch at run time. Measured against what the switch
saves — 2 to 3 % of a whole CLI run, about 5 % of refine — it doubles every
build of the extension for every user and every CI leg, adds an import shim
and a variant that ships with thinner tests. If heavy runs become routine
(the basin work, run nightly), (b) can be added later on top of (a) without
undoing anything: the option, the attribute, the guard and the `bench.py`
field are all reused.

## 11. Questions for Ola

- **Q1.** Option (a) or option (b)? Recommended: (a).
- **Q2.** The name. `RASPUTIN_HARDENING=ON/OFF` (CMake), `_core.hardening`,
  and `rasputin version`'s second line, `bounds checks: on (libc++ fast)`.
  Recommended as written.
- **Q3.** (a) only: the extra Ubuntu CI leg that builds OFF. Recommended: yes.
- **Q4.** The `--stats` line naming the mode, so every heavy-run report says it
  was unchecked. Recommended: yes (about 5 lines).

## 12. For `@tester` (after the ruling)

The red suite covers: the guard (§5, both halves, planted); `_core.hardening`
equals `"libc++ fast"` or `"libstdc++ assertions"` on the default build (a
Python test, platform-keyed); `rasputin version` prints the bounds-checks line
and still prints the bare version first; the `--stats` line (if Q4 is yes);
`bench.py`: the configure argv carries `-DRASPUTIN_HARDENING=<value>` for both
flag values, a `run.json` without the field loads as `"none"`,
`comparable()` names `hardening` when the modes differ, and a child whose
samples disagree is refused. The stale-cache plant of §3 is a manual check
recorded in the review, not a CI test (it needs two real installs).
No mutation round: nothing here is invariant-critical.

## 13. Documentation in the same PR

- `INSTALL.md`: a section "Bounds checks (on by default)": what they do, the
  measured cost with a link to `@perf`'s README, the two `-C cmake.define`
  commands and the plain-CMake line, the second-venv recipe for heavy runs,
  and `rasputin version` to tell which build is installed. Its "C++ changes
  have no effect" note also applies to the unchecked venv.
- `README.md`: one sentence where installation is mentioned, pointing at that
  section.
- `ROADMAP.md`: the row for 24, updated to shipped by the merge.
