# Harness h14: the asan+ubsan tests in parallel

Status: ruled; build next. Ola ruled on the principle in h13
(`docs/increments/h13-macos-ci-build.md`, §8 ruling 2: "yes", as its own
increment after h13), and §9 records his go-ahead on §8's three questions.
Mechanics, not a design question: one test-registration change (`@tester`) and
one workflow line (`@developer`), no production code, no new test (§6 says how
the change is checked instead).

Why: since h13 fixed the macOS build, the `C++ sanitizers (asan+ubsan)` job is
the slowest job the required `CI result` check waits for: 7.8 to 13.7 minutes
end to end in the four runs that include h13's change (§2a), against 5.9 to
8.5 for the slowest Python job. Most of that time is its Test step, where `ctest` runs
every test one at a time. Every merge runs it twice (once on the pull request,
once in the merge queue).

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature.* Two sources.

- The CTest manual, `ctest(1)` (CMake source `Help/manual/ctest.1.rst`, read
  2026-10-04): "`-j [<level>], --parallel [<level>]`: Run tests in parallel,
  optionally limited to a given level of parallelism." Since 3.29 the level
  "may be omitted, or `0`, in which case: [...] Otherwise, if the value is
  omitted, parallelism is limited by the number of processors, or 2, whichever
  is larger. Otherwise, if the value is `0`, parallelism is unbounded." The
  `COST` test property (`Help/prop_test/COST.rst`, read 2026-10-04): "When
  parallel testing is enabled, tests in the test set will be run in descending
  order of cost. Projects can explicitly define the cost of a test by setting
  this property to a floating point value. When the cost of a test is not
  defined by the project, ctest will initially use a default cost of `0`. It
  computes a weighted average of the cost each time a test is run and uses
  that as an improved estimate of the cost for the next run." A CI job starts
  from an empty build directory, so it never has that measured cost
  (`Testing/Temporary/CTestCostData.txt`), and every test costs 0.
- Ordering the longest jobs first is Graham's LPT (longest processing time)
  rule: R. L. Graham, "Bounds on multiprocessing timing anomalies", SIAM J.
  Appl. Math. 17(2):416-429, 1969 (doi:10.1137/0117039). Greedy list
  scheduling in an arbitrary order can finish up to 2 - 1/m times later than
  the best schedule on m processors; longest first, at most 4/3 - 1/(3m)
  times later. That bound is all LPT promises in general: it can miss the
  best schedule even when one job is longer than a quarter of the total (four
  workers, jobs of 10, 5, 5, 4, 4, 3, 3 and 3, total 37: longest first gives
  11, the best schedule 10, as 10 | 5+4 | 5+4 | 3+3+3; checked by trying all
  4^8 assignments). On this suite's data
  it is optimal: §2c's simulation finishes at 97 / 171 s, which is the lower
  bound (the run cannot finish before its longest test does). The design
  departs from LPT in one way: it gives a cost to one target, not to every
  test (§4, and why).

The runner: GitHub's runner reference (h13 §1, read 2026-10-04) gives
ubuntu-latest 4 CPUs and 16 GB for public repositories. The image's own
README (`actions/runner-images`, `images/ubuntu/Ubuntu2404-Readme.md`, read
2026-10-04) lists CMake 3.31.6 as the default, so the 3.29 behaviour above
applies.

*Legacy.* Nothing; the legacy tree has no CI configuration and no CTest
registration at all.

```
$ git ls-tree -r --name-only legacy-archive -- legacy | grep -iE '\.ya?ml$|travis|appveyor|\.github|CTest'
(no output, exit 1)
$ git grep -n -iE 'ctest|add_test|catch_discover' legacy-archive -- legacy
(no output)
```

## 2. Diagnosis

### 2a. Step timings

Read with `gh run view <id> --json jobs` (Build and Test step, and the whole
job, in minutes). The first four runs include h13's change (macOS Build back
to 1.4-2.4 minutes): 37232172673 is h13's own pull-request run, and merge-queue
run 37233115840 merged it. The rest are before it.

| run | event | asan Build | asan Test | asan job | tsan Test | tsan job | slowest Python job |
|---|---|---|---|---|---|---|---|
| 37235448420 | merge_group | 2.0 | 5.6 | 7.8 | 9.4 | 11.7 | 5.9 |
| 37234445967 | pull_request | 3.1 | 10.4 | 13.7 | 8.4 | 10.7 | 7.1 |
| 37233115840 | merge_group | 2.5 | 7.9 | 10.6 | 8.4 | 10.7 | 7.1 |
| 37232172673 | pull_request | 3.1 | 10.4 | 13.7 | 8.5 | 10.9 | 8.5 |
| 37230369996 | merge_group | 3.4 | 8.2 | | 9.5 | | |
| 37229043133 | pull_request | 2.9 | 6.5 | | 7.3 | | |
| 37226227284 | push | 3.5 | 10.2 | | 8.5 | | |
| 37224723259 | merge_group | 2.7 | 5.5 | | 4.6 | | |
| 37223176932 | pull_request | 3.7 | 7.7 | | 8.4 | | |

The tsan job is not one the `CI result` check waits for
(`.github/workflows/main.yaml`, the `result` job's `needs:` list), so it does
not hold a merge; §3 F says what to do with it.

### 2b. Where the Test step's time goes

Every Catch2 test case is its own ctest test, run as its own process:
`add_terrain_test` calls `catch_discover_tests(${name})`
(`tests/cpp/CMakeLists.txt`, in the function), and Catch2 v3.6.0's
`extras/CatchAddTests.cmake` (lines 63-64 of the copy CMake fetches into
`build/_deps/catch2-src`) runs each binary with `--list-tests --verbosity
quiet` and emits one `add_test` per name it lists. So 994 processes, each an ASan+UBSan Debug binary.

Per-test times parsed from the Test step's log (`gh run view <run> --log --job
<job>`, the `Test #N: <name> ... Passed <s> sec` lines), for the fastest and
the slowest asan Test step among the four runs that include h13's change:

| | 37235448420 (job 111533535723) | 37234445967 (job 111530695102) |
|---|---|---|
| ctest's `Total Test time (real)` | 334.09 s | 622.06 s |
| sum of the 994 per-test times | 332.9 s | 622.0 s |
| longest test | 97.0 s | 171.4 s |
| share of the longest 1 / 5 / 20 tests | 29 % / 52 % / 81 % | 28 % / 51 % / 82 % |
| tests under 0.1 s | 907 | 869 |
| median test | 0.02 s | 0.03 s |

The two runs differ by a factor of 1.9 with the same tests in the same
proportions, test by test; a factor of about 2 shows in the tsan job too
(§2a, 4.6 against 9.4). That is the runner hardware, which GitHub does not let a job
choose, not this suite.

The longest test, in both runs, is one Catch2 case: "ES9: refine_points with
the strip keeps J2 at the source points and E1 at the strip's", in
`tests/cpp/property/prop_refinement_edge_strip.cpp` (the case starting
`TEST_CASE("ES9: ...`). It is 16 `GENERATE` combinations (two terrains, two
seeds, four tolerances) inside one case, so ctest cannot split it, and it is
single-threaded (it calls `after_refine` and `run_points` with their default
`threads = 1`, same file). It is test #923 of 994 in declaration order. The
next longest: RP3 in `prop_refinement_refine_points` (27.9 / 50.2 s), two
`legalise_around` cases in `test_mesh_lawson_stack` (20.7 / 43.9 s), the T1
cases of `test_mesh_row_spans` (each 6-26 s).

### 2c. What parallel runs would give, from these numbers

Greedy list scheduling of the logged per-test times on N workers (each test
starts on the first free worker; a Python simulation over the parsed log,
assuming a test takes as long in parallel as it did alone):

| order | N = 2 | N = 3 | N = 4 |
|---|---|---|---|
| declaration order (ctest's order with no cost data) | 200 / 368 s | 163 / 300 s | 145 / 265 s |
| `COST` on `prop_refinement_edge_strip` only | 167 / 313 s | 112 / 214 s | 97 / 172 s |
| longest first, every test (LPT) | 166 / 311 s | 111 / 207 s | 97 / 171 s |
| lower bound: max(longest test, sum / N) | 166 / 311 s | 111 / 207 s | 97 / 171 s |

(fast run / slow run.) In declaration order, ES9 has 212 s (fast) or 406 s
(slow) of per-test time ahead of it to share out (the sum of the earlier
tests' times; in a parallel run without the cost it starts at 48 / 94 s in
the simulation), and then
runs alone for its 97-171 s. A cost on its one target puts it first and
reaches the lower bound at N = 4. Costs on more targets do not help and can
hurt, because the cost is per target, not per case: with the five targets
holding the longest cases all costed, N = 4 gives 129 / 234 s, worse than one
target, since their cases then run in declaration order among themselves and
ES9 waits behind them.

**What this has not shown.** The simulation assumes four workers as fast as
one. Whether the runner's 4 CPUs are four cores or two cores with two
hardware threads each is not in the logs, and memory bandwidth is shared.
If each test runs 1.5 times slower when four run at once, the N = 4 costed
estimate becomes about 2.4 / 4.3 minutes. §6 measures it.

### 2d. The precondition: do any two tests share state?

Since each test is its own process (§2b), memory, globals and statics in one
test cannot reach another. What could: files (a fixed name or a shared
temporary path), environment variables, other processes, IPC, the working
directory, ctest fixtures, and timing. All searched in this worktree at
`529613a`, over the C++ tests, the core headers and the core sources:

```
$ grep -rnE '#include <(fstream|cstdio|stdio\.h|filesystem|unistd\.h|fcntl\.h|sys/[a-z]+\.h|cstdlib)>' tests/cpp include src
tests/cpp/unit/test_refinement_scan_frozen.cpp:30:#include <cstdlib>
tests/cpp/unit/test_build_hardening.cpp:38:#include <sys/wait.h>
tests/cpp/unit/test_build_hardening.cpp:39:#include <unistd.h>
tests/cpp/unit/test_predicates_default_kernel.cpp:31:#include <cstdio>
tests/cpp/unit/test_predicates_default_kernel.cpp:40:#include <sys/wait.h>
tests/cpp/unit/test_predicates_default_kernel.cpp:43:#include <cstdlib>
tests/cpp/unit/test_predicates_default_kernel.cpp:44:#include <unistd.h>

$ grep -rnwE 'fopen|freopen|ofstream|ifstream|fstream|tmpnam|tmpfile|mkstemp|mkdtemp|temp_directory_path|unlink|rename|creat' tests/cpp include src
(three hits, all the English word "rename" in comments:
 test_edge_properties.cpp:24, test_pslg_builder.cpp:961, prop_cdt_invariants.cpp:398)

$ grep -rnE '"[^"]*(/tmp|\.txt|\.csv|\.tif|\.bin|\.json|\.dat|\.out|\.log)"' tests/cpp include src
(no output, exit 1)

$ grep -rnE 'getenv|setenv|putenv|unsetenv|environ\b' tests/cpp include src
(no output, exit 1)

$ grep -rnE '\bfork\(|exec[lv]p?\(|system\(|popen|socket\(|bind\(|shm_open|mmap\(|sem_open|mkfifo|flock|chdir' tests/cpp include src
tests/cpp/unit/test_build_hardening.cpp:96:    const pid_t pid = fork();
tests/cpp/unit/test_predicates_default_kernel.cpp:352:    const pid_t pid = fork();

$ grep -rnE 'FIXTURES_|RESOURCE_LOCK|RUN_SERIAL|WORKING_DIRECTORY|ENVIRONMENT|OUTPUT_DIR|REPORTER|PROCESSORS|COST' CMakeLists.txt tests/cpp/CMakeLists.txt
(no output, exit 1)
```

What the hits are:

- `<cstdio>` in the default-kernel suite is for `std::fflush(nullptr)` before
  its fork. `_exit`, which both fork tests call in the child, comes from
  `<unistd.h>`; `test_build_hardening.cpp` has no `<cstdlib>`. In
  `test_refinement_scan_frozen.cpp:30` `<cstdlib>` declares `std::size_t`, the
  only name the file uses that it can come from (it calls no `<cstdlib>`
  function). The `<cstdlib>` at `test_predicates_default_kernel.cpp:43`, in
  its fork block, declares nothing the file uses (`exit` and `atexit` appear
  only in the comment at lines 281-284, and `_exit` is `<unistd.h>`'s). No
  file is
  opened anywhere in the C++ tests or the core (the core never opens a file,
  `CLAUDE.md` §2, I/O boundary).
- The two `fork()` calls: each forks a child, waits for that child's own pid
  (`waitpid(pid, ...)`), and reads only its exit status. Nothing is shared
  with another test's process.
- No ctest property couples tests: no fixtures, no resource locks, no
  `RUN_SERIAL`, no per-test environment, no Catch2 reporter writing a file.
  Every test runs in its target's build directory (Catch2's default
  `WORKING_DIRECTORY`), and none writes there.
- ASan and UBSan report to stderr by default; the workflow sets no
  `ASAN_OPTIONS` or `UBSAN_OPTIONS` (`log_path` would write files).

Timing, the remaining way one test can affect another (by loading the
machine):

```
$ grep -rnE 'steady_clock|high_resolution_clock|system_clock|sleep_for|sleep_until|this_thread::yield|BENCHMARK|wait_for|wait_until' tests/cpp
```

returns `test_refinement_chunks_dynamic.cpp` (lines 168-300) and
`test_refinement_refine.cpp` (lines 12-49). Read:

- `test_refinement_refine.cpp` checks only that the phase timings refine
  reports are non-negative and sum to no more than a `steady_clock` interval
  around the same call. Load stretches both sides alike; it cannot fail
  because the machine is busy.
- `test_refinement_chunks_dynamic.cpp` forces block arrival orders with waits
  that time out after 10 s, plus 20-30 ms sleeps. The sleeps only make a
  *wrong* scheduler more likely to fail; a correct one rethrows the lowest
  block whatever the timing, so load cannot fail a correct one. The 10 s
  timeouts wait for another thread to start or reach a counter; under four
  sanitized processes on four CPUs that is milliseconds, not seconds. This is
  the one residual risk, and it is small: a failure would name the case and
  `timed_out`, not hide.

Threads: the refine suites start up to 8 threads, or as many as the machine
has (`threads == 0` means `hardware_concurrency`, `parallel_util/chunks.hpp`).
With four tests at once that is more threads than CPUs for a while. It costs
time, not correctness, and is already the case on Ola's machine whenever he
runs `ctest -j`.

Memory: not measured. Four ASan processes at once must fit in 16 GB; no test's
peak was recorded in CI, and a local measurement would be AppleClang on a
different allocator. A test killed for memory fails ctest, so §6 item 1
catches it.

**Precondition met**: no two tests write the same file or share any other
state that the searches above can see.

## 3. Options

| option | expected asan Test step | cost and what it gives up |
|---|---|---|
| **A. `ctest --parallel "$(getconf _NPROCESSORS_ONLN)"`** on the asan Test step | 2.4-4.4 min ideal (145 / 265 s), up to about 3.6-6.6 min if each test slows by 1.5 | One line. ES9 still starts late (§2c). |
| **B. A + a `COST` on `prop_refinement_edge_strip`** | 1.6-2.9 min ideal (97 / 172 s), about 2.4-4.3 with the 1.5 slowdown | About 8 lines of `tests/cpp/CMakeLists.txt` (§4, LOC). A hand-placed hint that can go stale: if another test becomes the longest, the run gets slower, never wrong. |
| C. Shard the tests over a matrix of jobs (`ctest -I` strides) | no better than B | Each shard builds the whole tree again (2.0-3.7 min each), or the job ships the sanitized binaries between jobs as an artifact. ES9 is still a floor of 97-171 s in whichever shard has it. |
| D. Split the job (build once, test in several jobs) | no better than B | Same as C's artifact route, more YAML, same floor. |
| E. Split ES9 into several Catch2 cases | B's floor drops to the sum / 4 bound: 83 / 155 s | A test edit for 15 s; ES9's `GENERATE`s are the test's design. Not worth it now. |
| F. The tsan job too (run its 20 binaries in parallel) | tsan Test 4.6-9.5 min, maybe halved | Not a check `CI result` waits for, so no merge gets faster. Its Test step is a shell loop, not ctest; parallel needs `xargs -P` with interleaved output, and the suites start many threads under TSan already. Separate question (§8). |
| G. Compile the sanitizer build at `-O1` instead of Debug (`-O0`) | perhaps a half or a third of today's, unmeasured | Changes what the job builds, not how it runs; outside Ola's ruling. Mentioned only. |
| H. `--parallel` on the C++ core legs' `ctest` too | saves 10-20 s on jobs of 2-3 min | Same precondition, already met. Not the floor; out of the ruling's scope (§8). |

After B, the asan job's expected total is about 4.5-8.5 minutes (set-up and
Configure under a minute, Build 2.0-3.7, Test 2.4-4.3 with the slowdown), so on
a fast runner, or at the ideal, the slowest Python job (5.9-8.5 in §2a)
becomes the floor `CI result` waits on. On a slow runner with the slowdown the
asan job can stay the floor: in run 37234445967 it would be about 3.1 + 4.3 +
0.2 = 7.6 minutes (Build, Test, the rest), above that run's Python job at
7.1. The saving on a run is the asan job's time above whichever job is then
the floor: about 2 minutes in run 37235448420 (7.8 - 5.9), about 6 to 7 in
run 37234445967 (13.7 - 7.6 with the slowdown, 13.7 - 7.1 at the ideal).

## 4. Recommendation

Option B.

**Workflow** (`.github/workflows/main.yaml`, `sanitizers` job, Test step):

```yaml
      - name: Test
        # One test per CPU. Every Catch2 case is its own process and none
        # shares a file or other state (h14 §2d).
        run: ctest --test-dir build-san --output-on-failure --parallel "$(getconf _NPROCESSORS_ONLN)"
```

- The job count is the same expression the build lines use (h13), so the
  workflow has one idiom for it.
- If `getconf` ever printed nothing, the line becomes `--parallel ""`. Probed
  locally (CMake 4.4.3, a test project with no compiler): ctest then runs in
  parallel, limited to the CPU count, and exits 0. That is the 3.29 rule of
  §1, and the runner's CMake is 3.31.6. Unlike the build lines (h13 §4), an
  empty value here is safe, not unlimited.
- `--output-on-failure` still prints a failing test's whole output in one
  block under its result line, not interleaved with tests running beside it
  (probed locally: a failing test that prints across a second, run beside a
  noisy one, `--parallel 4`).
- No per-test timeout is set today, and none is added: `grep -rn
  'TIMEOUT\|include(CTest)' CMakeLists.txt tests/cpp/CMakeLists.txt` prints
  nothing (exit 1). The job's `timeout-minutes: 30` still catches a
  hang.

**Test registration** (`tests/cpp/CMakeLists.txt`): `add_terrain_test`
accepts an optional `COST <value>` keyword and passes it to
`catch_discover_tests(... PROPERTIES COST <value>)`; the helpers that forward
`${ARGN}` (`add_terrain_backend_test`, `add_terrain_cdt_test`) pass it through
unchanged. `prop_refinement_edge_strip` is registered with `COST 100` and a
one-line comment: its ES9 case runs 97-171 s under the sanitizers, and
without a cost ctest starts it 923rd (h14 §2c).

- The interface: one keyword on the existing helper, not a second
  `catch_discover_tests` call (calling it twice on one target registers every
  case twice) and not a per-target variable read inside the helper. A
  suite that is not passed `COST` is registered exactly as today.
- `PROPERTIES` applies to every case of the target, not to ES9 alone (§1,
  Catch2's own documentation of the option). That is why one target and not
  five (§2c): the other 16 cases of `prop_refinement_edge_strip` take 5-10 s
  together and finish while ES9 runs.
- An explicitly set cost of 0 behaves as unset: probed locally, ctest still
  replaces it with the measured cost on a second run. So the unset branch can
  be written either way.
- Locally nothing gets worse: a second `ctest -j` run in the same build
  directory already orders by measured cost (probed), and the costed target
  goes first on the first run too.
- The CMake test file is test infrastructure, not a test: no case, assertion
  or `GENERATE` changes, and the test count stays the same (§6 item 1).

**Constant: `COST 100`.** Only its order relative to the other costs matters:
every other test's cost is 0 in CI. It assumes no other target is given a cost
above 100; checked at this tree, where no test has a `COST` (the last grep
of §2d).

**Constant: one ctest job per CPU.** It assumes four sanitized test processes
fit in the runner's 16 GB, about 4 GB each. Not measured (§2d, Memory); §6
item 1 fails if one is killed.

**LOC:** 0 production lines. About 3 lines of workflow (one changed, two
comment) and about 8 of `tests/cpp/CMakeLists.txt` (the keyword parsing, the
`catch_discover_tests` call, the edge-strip registration and its comment).

## 5. Exclusions

- The tsan job, the C++ core legs' `ctest`, and the sanitizer optimisation
  level stay as they are (§3 F, H, G; §8 questions 1 and 2).
- No test, case or `GENERATE` changes (§3 E).
- No change to which jobs are required, to `timeout-minutes`, or to the build
  lines.

## 6. How the change is checked

No unit test can see a workflow's job count or ctest's start order, so there
is no red step (as in h10 and h13). The check is the pull request's own CI
run, and every item can fail. In the `C++ sanitizers (asan+ubsan)` job's log:

1. ctest reports `100% tests passed` out of the same number of tests as the
   `C++ core (ubuntu-latest)` leg of the same run (994 at this tree; compare
   the two, do not hard-code it). A test killed for memory or failing under
   load breaks this.
2. **Parallelism took effect:** ctest's `Total Test time (real)` is at most
   **0.5 times** the sum of the per-test `Passed ... sec` times in the same
   log. Today the two are equal (334.09 against 332.9 s; 622.06 against
   622.0 s). This is independent of how fast the runner is, which a bare time
   limit is not (§2b's factor of 1.9). Four workers give about 0.25 to 0.3;
   a run with the job count lost gives 1.0.
3. **The cost took effect:** the log line `Start <n>: ES9: refine_points ...`
   is time-stamped less than **30 s** after the first `Start` line of the
   Test step. Today, in serial, ES9 starts 213 s (fast runner) to 406 s
   (slow) after the first test, as the command below prints. Parallel without
   the cost would start it at about 48 / 94 s (§2c's simulation), so the
   30 s limit separates the two cases with room on both runners. With the
   cost, only the other 16 edge-strip cases (5-10 s of work, shared over four
   workers) can start before it.
4. **The time:** the Test step under **5.0 minutes** (today 5.6-10.4 in the
   four runs that include h13's change; the estimate is 2.4-4.3). If items 2 and 3 hold and this fails, the
   runner is slower under load than §2c assumed: the increment goes back to
   `@architect` with the measured numbers, not on to another option.

A command that reads items 2 and 3 from a log saved with `gh run view <run>
--log --job <job> > asan.log`. Run against today's two logs, it prints sums of
332.94 and 622.04 s against totals of 334.09 and 622.06 s, and ES9 starting
213 and 406 s after test 1:

```bash
awk -F'\t' '$2=="Test"' asan.log | grep -oE 'Passed +[0-9.]+ sec' \
  | awk '{s+=$2} END {print "sum of per-test times", s}'
grep -m1 'Total Test time' asan.log
awk -F'\t' '$2=="Test"' asan.log | grep -E ' Start +[0-9]+: ' \
  | sed -En '1p;/Start +[0-9]+: ES9/p' | cut -f3 | cut -c1-90   # first Start, and ES9's
```

`@reviewer` is read-only. Its handback after the PR's CI run carries the
measured numbers for items 1 to 4, and the spawner records them, verbatim, as
a review round in `## Review` below. The numbers exist only after the push.

## 7. Blueprint

One step and one registration change; data flows as before:

```
tests/cpp/CMakeLists.txt
  add_terrain_test(name sources... [COST c])     # new optional keyword
     -> catch_discover_tests(name [PROPERTIES COST c])
     -> CTestTestfile: one add_test per Catch2 case, COST c on each of them
.github/workflows/main.yaml, sanitizers job
  Configure, Build      unchanged
  Test                  ctest ... --parallel "$(getconf _NPROCESSORS_ONLN)"
     -> ctest orders by COST (edge strip first), then declaration order
     -> 4 processes at a time; output per test, whole, on failure
```

Who does what, in this order, one commit each:

1. `@tester`: the `COST` keyword in `add_terrain_test` and `COST 100` on
   `prop_refinement_edge_strip`, both in `tests/cpp/CMakeLists.txt` (under
   `tests/`, `@tester`'s write limit; `@developer`'s `CMakeLists.txt` entry is
   the root file, `docs/increments/h6-role-limits.md` §3). No case, assertion
   or `GENERATE` changes. On its own this commit changes nothing in CI, which
   still runs the tests one at a time.
2. `@developer`: the sanitizer job's `ctest` line in
   `.github/workflows/main.yaml` (`.github/`, `@developer`'s write limit).

There is no red step: no test assertion changes, and no unit test can see the
job count or the start order, so the PR's own CI run is the check (§6).
`@reviewer` audits before the first push; after the PR's CI run, `@reviewer`
again with §6's numbers. `@perf` is not involved: no refine or mesh code is
touched.

## 8. Questions for Ola

1. Run the thread-sanitizer job's 20 suites in parallel too? Default: **no**.
   That job does not hold a merge (the `CI result` check does not wait for
   it), so it would save runner time but no waiting, and its test step is a
   shell loop that would need rewriting.
2. Add the same `--parallel` to the ordinary C++ test jobs (Linux and macOS)?
   Default: **no**. They take 2-3 minutes in all, the tests there take 10-30
   seconds, and your ruling named the sanitizer job.
3. Include the cost hint for the one long test (the edit to
   `tests/cpp/CMakeLists.txt`, §4), or only the workflow line? Default:
   **include it**. Without it the slowest test starts near the end and the
   step takes roughly 1.5 times as long (§2c: 145 against 97 s on a fast
   runner, 265 against 172 s on a slow one).

## 9. Ola's rulings

Ola's words, the whole message: "WHen h14 and h12 C2 should get going."
(2026-10-05T04:57:55Z, the main session's transcript). He did not answer §8's
questions one by one. The main session read the message as a go-ahead with
§8's three defaults, and told Ola so in the same turn; each one below is taken
on that reading, and Ola can overturn any of them.

1. **The thread-sanitizer job stays as it is** (§8 question 1, default no).
2. **No `--parallel` on the ordinary C++ test jobs**, Linux and macOS (§8
   question 2, default no).
3. **The cost hint is included**: `COST 100` on `prop_refinement_edge_strip`
   (§8 question 3, default include it; §4).

## Review

### Round 1: `@reviewer`, design, 529613a..0b2ab5d

Verdict: **CHANGES REQUESTED**, on wording only. Every number and the
recommendation reproduced: the step timings, the per-test parse, ES9's 97.03
and 171.4 s, the simulation table cell for cell, the shared-state searches,
the CMake 3.31.6 source, the quoted documentation, and §6's command; all four
§6 checks can fail.

Findings, all fixed in the commit after 0b2ab5d:

1. §2b named `--list-test-names-only` (Catch2 v2); v3.6.0's
   `CatchAddTests.cmake` runs `--list-tests --verbosity quiet` (its lines
   63-64).
2. Intro and §2a said "since/after h13 merged"; run 37232172673 is h13's own
   pull-request run, before queue run 37233115840 merged it. Now "the four
   runs that include h13's change".
3. Intro and §3 gave the slowest Python job as 5.5-8.5 min; the table gives
   5.9-8.5 (5.5 is the fastest Python job).
4. §6 item 4 gave 5.5-10.4 after h13; those rows give 5.6-10.4 (5.5 is the
   pre-h13 run 37224723259).
5. §6 item 3 gave ES9's start as 212 s; §6's own command prints 213 s (the
   timestamp gap, 213.3 s). 212 s is the sum of the earlier tests' times,
   which §2c now names as such.
6. §3 B said two to five lines of `tests/cpp/CMakeLists.txt`, §4 about 8. Now
   about 8 in both.
7. §1 said longest first is optimal *because* one test exceeds a quarter of
   the run; it is not (the counterexample first given here, jobs 10, 7, 7, 7,
   7, was wrong; round 2, finding 1, replaced it). Now: optimal on this data, as §2c's simulation reaches the lower
   bound (97 / 171 s).

Suggestions, taken: §2d's `<cstdlib>`/`_exit` remark corrected (`_exit` is
`<unistd.h>`; `test_build_hardening.cpp` has no `<cstdlib>`; the one in
`test_refinement_scan_frozen.cpp` is for `std::size_t`); §4's no-timeout
claim now shows its grep; §3's saving "about 6" became "about 7" (13.7 - 7.1), itself
corrected in round 2, finding 2, and
Build 2.0-3.7; §6 item 3 compares against parallel without the cost (48 / 94 s
simulated).

### Round 2: `@reviewer`, design, 0b2ab5d..05c6eaf

Verdict: **CHANGES REQUESTED**. Both errors came from round 1's own
suggestions.

Findings, fixed in the commit after 05c6eaf:

1. §1's counterexample (four workers, jobs 10, 7, 7, 7, 7) is not one: the
   best schedule is also 14 (10 | 7 | 7 | 7+7). Replaced by jobs 10, 5, 5, 4,
   4, 3, 3, 3 on four workers: longest first 11, best 10 (10 | 5+4 | 5+4 |
   3+3+3), with the longest job, 10, above a quarter of the total 37 (9.25).
   Brute-forced by `@reviewer`, the main session and `@architect` (all 4^8
   assignments). Round 1's finding 7 now says its example was wrong.
2. §3's slow-run saving ignored §3's own 1.5 slowdown: the asan job after
   option B would be about 3.1 + 4.3 + 0.2 = 7.6 minutes in run 37234445967,
   above that run's Python job (7.1), so asan stays the floor there. Now
   "about 6 to 7 (13.7 - 7.6 with the slowdown, 13.7 - 7.1 at the ideal)",
   and the Python job is the floor only on a fast runner or at the ideal.
   Round 1's "Suggestions, taken" line is corrected to match.

Suggestion, taken: §2d says what the `<cstdlib>` at
`test_predicates_default_kernel.cpp:43` is for: nothing the file uses.

**Round 3, 2026-10-05 (summarised from `@reviewer`'s handback).** Range `05c6eaf..de233d6`. Verdict: APPROVED. Both round-2 findings closed: the counterexample brute-forced over all 4^8 assignments (best 10, longest-first 11); the slow-run saving recomputed (6.1 with the slowdown, 6.6 without). The §2d claim checked with `git log -S` (2b3ae32) and grep. `check_citations.py` clean. Non-blocking: §2d's "declares nothing the file uses" would read better as "nothing the file does not already get from `<cstddef>`".
