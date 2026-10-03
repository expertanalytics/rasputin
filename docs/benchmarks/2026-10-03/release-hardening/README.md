# Release hardening: what bounds checks cost on the 1 m benchmark (@perf, 2026-10-03)

The ROADMAP row "Release hardening": measure what libc++'s
`_LIBCPP_HARDENING_MODE_FAST` costs in Release (bounds checks on
`std::vector::operator[]`, `std::span`, and the like), and libstdc++'s
equivalent, `_GLIBCXX_ASSERTIONS`, for CI's GCC legs. Measurement only:
nothing in the build files changed. The define went in through
`CMAKE_CXX_FLAGS` on the configure command line.

## Result

The meshes are identical with and without the checks: same sha256, same
quality and the same refine counters in every sample (see Quality). The cost
is time only.

| | libc++ FAST (AppleClang 21) | libstdc++ ASSERTIONS (GCC 16.2) |
|---|---|---|
| refine, 1 thread, tile / quarter | **+4.9 % / +4.5 %** | **+2.6 % / +2.5 %** |
| refine, CLI default threads, tile / quarter | +5.7 % / +6.3 % | +6.9 % / +7.2 % |
| refine, 2 to 20 threads (range of pooled medians) | +4.4..+8.0 % | +4.1..+8.7 % |
| `rasputin mesh` in process (`app_s`: read, mesh, write), default threads | +4.0 % / +4.1 % | +6.1 % / +5.2 % |
| whole child process (`proc_s`, with interpreter start-up), default threads | +1.9 % / +2.6 % | +3.2 % / +2.6 % |
| refine ceiling (1 thread over best), plain → hardened, tile | 2.49x → 2.46x | 2.44x → 2.34x |

Pooled medians over three back-to-back pairs, 5 samples per pair per thread
count, so 15 samples per side. Per-pair spread at 1 thread: clang tile
+4.3..+5.6 %, quarter +3.9..+5.6 %; GCC tile +2.4..+3.7 %, quarter
+2.0..+2.9 %. Every pair, every domain and every thread count came out
positive for refine and for `app_s`. The full per-thread tables, with each
pair's figure, are in `raw/analysis-clang.md` and `raw/analysis-gcc.md`.

**What this means for acceptance.** About 5 to 6 % of refine is at or just
over `bench.py`'s 5 % regression threshold. A PR that turns the checks on
will get `REGRESSION` verdicts against a baseline built before the switch.
In that PR the switch is the measured, deliberate cause, so the acceptance
run for it re-baselines. The runs here are the evidence for that
(h1/h2: `REGRESSION` at 19 and 20 threads, +5.6 and +5.9 %).

**Where the cost goes is not profiled.** One observation: with GCC the
overhead grows with the thread count, from +2.5 % at 1 thread to about +7 %
at 8 to 20 threads. With clang it stays roughly flat (+4.5 to +6.5 %). That
suggests the checks cost more in the parallel part than in the serial phase
under GCC, but nothing here measures it. Treat it as an observation, not a
cause.

## Recommendation (the decision is Ola's)

Switch both on, for every build type, and accept about 5 % on refine and about
4 % on an in-process `rasputin mesh` of the benchmark tile (5 to 6 % with
GCC). On a whole CLI process the cost is about 2 to 3 %.

The gain is an exact stop where there is now undefined behaviour. An
out-of-range `operator[]` that the tests miss today reads or writes adjacent
memory and carries on. With FAST mode it is a `brk` trap at the faulting call
(SIGTRAP), and with `_GLIBCXX_ASSERTIONS` it is an abort that names the
assertion (SIGABRT, `raw/probes.txt`). It is still a crash of the Python
process, not a Python exception: hardening turns silent corruption into a
deterministic stop. It does not make the error recoverable. That was not
measured inside Python. A deliberate crash in a Python process puts a crash
dialog on Ola's screen, so the probe is a C++ binary.

If 5 % is too much for the basin-scale work, a middle option is open but
unmeasured: hardening on in CI's Release legs and the Python legs, off in a
`RASPUTIN_UNCHECKED` build for production runs. That splits what is tested
from what is shipped, which is why it is not the recommendation.

## What a PR that switches it on touches

- `CMakeLists.txt`, one statement on the header target, so that it reaches
  `terrain_predicates`, `terrain_cdt`, `_core` and every ctest suite through
  `terrain_headers`' INTERFACE:
  `target_compile_definitions(terrain_headers INTERFACE
  _LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST _GLIBCXX_ASSERTIONS)`.
  Each standard library ignores the other's macro, so one line covers macOS
  (libc++) and the Ubuntu legs (libstdc++), including the wheel built by
  scikit-build-core, which runs the same `CMakeLists.txt`
  (`pyproject.toml` `[tool.scikit-build]`). The sanitizer and TSan jobs get it
  too; nothing here ran them that way. Catch2 is built by FetchContent outside
  `terrain_headers`, so its own TUs stay unchecked; whether hardened and
  unhardened TUs in one test binary need care was not checked here.
- `.github/workflows/main.yaml`: nothing needs to change. The `cpp` job's
  `ubuntu-latest` leg (system GCC, libstdc++) and `macos-latest` leg
  (AppleClang, libc++) pick the definition up from CMake. `macos-latest`
  needs an Xcode whose libc++ knows `_LIBCPP_HARDENING_MODE` (LLVM 18 and
  later; older libc++ ignores an unknown macro silently). A guard test
  catches that.
- A test, `@tester`'s: a ctest binary that indexes one past the end of a
  `std::vector` and is registered `WILL_FAIL`, or that `static_assert`s the
  mode macro. Without it, a toolchain that ignores the macro passes CI
  unhardened. A crash in a ctest binary is fine; in the Python process it is
  not (see above).
- The acceptance run of that PR re-baselines (see "What this means for
  acceptance").
- About 3 lines of production code and one test file.

## Method

- Tree `390b516` (master), clean. `tools/bench.py` blob
  `6c5de4c0dee4500fbc1934a607fca22c3d21aa92`, unchanged. Default inputs:
  DEM `tests/fixtures/dem_archive/7908_3_10m_z33.tif`, tolerance 1 m,
  domains `tile` and `quarter` (`docs/benchmarks/2026-09-26/quarter.geojson`),
  threads 0 (the CLI default) and 1 to 20, 5 repeats interleaved over the
  thread counts.
- Apple M1 Max, 8 P + 2 E cores, 32 GiB, macOS 27.0, Python 3.14.7,
  numpy 2.5.3 (the main checkout's `.venv`). **Power: AC, 100 %, charged**,
  before and after every run (`pmset -g batt`, in each `run.json`). No other
  agent ran during the measurements. A `caffeinate` assertion was held.
- Four Release builds of `_core` from the one tree, each configured at
  `<tree>/build-bench` with `bench.py`'s own configure line plus, for the
  hardened ones, `-DCMAKE_CXX_FLAGS=<define>`. CMake keeps a cached
  `CMAKE_CXX_FLAGS` when `bench.py` reconfigures without it. Between runs,
  `scripts/pairs.sh` swaps the build directories in and out of `build-bench`
  with `mv`. `bench.py run` then rebuilt (a no-op) and assembled its `pkg/`
  as usual. Every `run.json` records the `_core` sha256, and each side's hash
  is the same across its three runs (`raw/analysis-*.md`, "runs").
  - clang plain `41a22a61…`, clang FAST `bcfda7ba…`: AppleClang 21.0.0,
    `-O3 -DNDEBUG -flto`.
  - GCC plain `32aa3ad3…`, GCC ASSERTIONS `f058d5ef…`: Homebrew GCC 16.2.0
    (`-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16`), `-O3 -DNDEBUG
    -flto=auto`, linked against Homebrew's libstdc++.
  - `run.json`'s `build.cxx_flags_release` shows only `-O3 -DNDEBUG`, because
    `bench.py` records `CMAKE_CXX_FLAGS_RELEASE`, not `CMAKE_CXX_FLAGS`. The
    flags as compiled are in `raw/probes.txt`, from each build's
    `flags.make`.
- **The probes that the hardened objects are the ones measured and do what
  they say** (`scripts/probes.sh` → `raw/probes.txt`):
  - FAST mode adds 149 `brk` traps to `_core` (168 → 317).
  - The GCC hardened `_core` imports `__glibcxx_assert_fail` and the plain
    one does not.
  - `scripts/oob.cpp` indexes one past the end of a `std::vector` and a
    `std::span` at `-O3 -DNDEBUG`. It prints a value read past the end and exits 0 without the
    define. With it, it exits 133 (SIGTRAP, clang) or 134 (SIGABRT, GCC,
    naming `__n < this->size()`). So the probe can fail, and with the define
    it does stop the program.
- Order: clang `plain fast fast plain plain fast` (12:49 to 13:02), then GCC
  in the same pattern (13:03 to 13:17). The order is balanced so that a drift
  over the batch does not land on one side. Pairs: (h1, h2), (h4, h3),
  (h5, h6), and the same for g. `scripts/analyse.py` takes each pair's
  median ratio and the pooled median over all 15 samples per side.
- Not sanitized first: the perf role's sanitize-before-numbers rule binds a
  patch that changes C++, and no source changed here, only a predefined
  macro.
- The verdict lines inside each `raw/<run>/README.md` are `bench.py`'s own
  pick of a baseline among the runs in the scratch out-root. They are not this
  study's comparison. g1-gplain's verdict, for one, is against a clang run.
  The comparison is `raw/analysis-*.md`.
- End to end through the extension: `app_s` is `cli.app()` in process
  (`rasputin mesh --dem … --tolerance 1 --binary`: read the GeoTIFF, mesh,
  write). `proc_s` is the child process from the parent's clock. No separate
  Bygdin run was made: the benchmark tile already goes through the CLI.

## Quality (identical on all twelve runs)

| domain | worst angle | max degree | within tolerance | Delaunay violations | mesh sha256 |
|---|---:|---:|---|---:|---|
| tile | 0.6296° | 74 | yes | 0 of 692,056 | `11741a81adfa17b34a1ea56875ec9e791414d17625eb1c58305e7c0cf10e005c` |
| quarter | 0.3955° | 18 | yes | 0 of 641,791 | `ccebf96a86c6c5e244e4a0281919de4e866fcfe789b66024a290ac2af33771a1` |

These are the same hashes as the 2026-10-02 runs (`2026-10-02/15c-2-r2`).
GCC 16 produces the same mesh as AppleClang. The (domain, threads, rounds,
inserted, flips, max_error) tuples come to 42 distinct values over all
samples of each batch: one per (domain, thread count). So the checks changed
no decision anywhere in refine.

## Files

- `raw/<run>/`: `bench.py`'s `run.json`, `raw.tsv` and generated `README.md`
  for the twelve runs. `h*` are clang, `g*` are GCC.
- `raw/batch1-clang.out`, `raw/batch2-gcc.out`: the driver's log, with
  `pmset` per run.
- `raw/analysis-clang.md`, `raw/analysis-gcc.md`: per-pair tables, from
  `python scripts/analyse.py . h1-plain:h2-fast h4-plain:h3-fast
  h5-plain:h6-fast`, and the same with `g1-gplain:g2-gassert
  g4-gplain:g3-gassert g5-gplain:g6-gassert`.
- `raw/probes.txt`, from `bash scripts/probes.sh`.
- `scripts/pairs.sh`: the driver. Batch 1 ran with an earlier, two-variant
  version of the same swap logic. The copy here is the generalised one that
  ran batch 2.
- Meshes (the `--ascii` quality runs) were left out of the repository, in
  `../rasputin_scratch/hardening-2026-10-03/meshes/<run>/`. To regenerate
  one, rerun the `bench.py _child` line in any `raw/<run>/README.md` with
  `--ascii --out PATH`.

## Review

### Round 1 (@reviewer, 2026-10-03): CHANGES REQUESTED

Reviewed `390b516..305748c`. Production LOC: 0 (docs and evidence only; the
`CLAUDE.md` §2 ceiling does not apply). CI: none; the branch is not pushed, so
it is not merge-ready until `gh pr checks` is green after a push.

Checked and holding: `scripts/analyse.py` re-run on the committed `run.json`s
reproduces `raw/analysis-clang.md` and `raw/analysis-gcc.md` byte for byte;
every figure in the Result table, the per-pair spreads at 1 thread, the
ceilings, the `_core` hashes, the 149 extra `brk` traps, the mesh hashes and
the 42 counter tuples match them; the only negative per-pair figure is in
`proc_s` (clang quarter, 1 thread, pair 1), so "positive for refine and for
`app_s`" holds; the crash-not-exception point is labelled as not measured in
Python; the CMake and CI claims match `CMakeLists.txt`,
`tests/cpp/CMakeLists.txt`, `pyproject.toml` and `.github/workflows/main.yaml`;
`tools/check_citations.py` passes.

Blocking:

1. "What this means for acceptance": "h1/h2: `REGRESSION` at 19 and 20
   threads, +5.6 and +5.9 %" is false. Those are the last two of 39
   `REGRESSION` lines in `raw/h2-fast/README.md` (`pairs.sh` logs only
   `tail -2`): h2 regresses against h1 at every thread count but t=17 on the
   tile and t=2 and t=14 on the quarter.
2. "Where the cost goes": clang "stays roughly flat (+4.5 to +6.5 %)"
   contradicts the Result table (+4.4..+8.0 over 2 to 20 threads) and
   `raw/analysis-clang.md` (quarter t=12 +8.0, t=13 +7.1; tile t=14 +6.8).
   State the range the tables give; the contrast with GCC is in how much the
   overhead grows from 1 thread (clang about 1 to 3 points, GCC about 4 to 5).
3. Reproducibility: the four builds are not scripted. The configure and build
   commands, the `.variant` markers `pairs.sh` and `probes.sh` read, and the
   renames from `build-park-*` to `build-clang-*` that `probes.sh` expects
   are in no script and only paraphrased in Method. Check the commands in
   (a `build.sh`, or the exact lines in Method), and take the worktree path
   `W` the scripts hard-code from an argument or say it must be edited.

Non-blocking: "A guard test catches that" (the `.github` bullet) describes a
test that does not exist yet; "the guard test below would catch that".
