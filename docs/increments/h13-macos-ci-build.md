# Harness h13: the macOS C++ job in CI

Status: edit done (a77879f); approved by `@reviewer`, round 2 (see Review);
awaiting the push, then the PR's CI check (§6), whose numbers are recorded
as a further review round.

Why: the `C++ core (macos-latest)` job in `.github/workflows/main.yaml` takes
20-26 minutes on every run, against 2-4 minutes for the same job on
ubuntu-latest. It is a required check, and every merge runs it twice (once on
the pull request, once in the merge queue), so it sets the floor for how long
a merge takes. It is also within 4-10 minutes of the job's own
`timeout-minutes: 30`, so a slightly slower runner would fail it outright.

Independent of h11 (PR #176: the `CI result` gate, push trigger dropped) and
h12 (prose-only fast lane, maybe heavy jobs only in the queue): this change is
one argument on the build lines, and is correct whichever of those land.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature.* The CMake manual, `cmake --build` (cmake.org, `cmake(1)`, read
2026-10-04): "`-j [<jobs>]`, `--parallel [<jobs>]`: The maximum number of
concurrent processes to use when building. If `<jobs>` is omitted the native
build tool's default number is used. The `CMAKE_BUILD_PARALLEL_LEVEL`
environment variable, if set, specifies a default parallel level when this
option is not given." (§2c shows the second sentence's consequence: a bare
`--parallel` is "given", so the variable is ignored.) For the Makefile
generator, the native tool's default under a bare `-j` is GNU make's: "If
there is nothing looking like an integer after the '-j' option, there is no
limit on the number of job slots." (GNU make manual, *Parallel Execution*,
gnu.org, read 2026-10-04). The `make(1)` man page says the same. GitHub's runner reference (docs.github.com, "GitHub-hosted
runners", read 2026-10-04), public repositories: ubuntu-latest has 4 CPUs and
16 GB of memory; macos-latest (arm64) has 3 CPUs (M1) and 7 GB;
macos-15-intel has 4 CPUs and 14 GB. The repository is public
(`gh repo view --json visibility` returns `PUBLIC`).

*Legacy.* Nothing; the legacy tree has no CI configuration at all.

```
$ git ls-tree -r --name-only legacy-archive -- legacy | grep -iE '\.ya?ml$|travis|appveyor|\.github'
(no output)
$ git grep -n -iE 'parallel|ccache|macos' legacy-archive -- 'legacy/.github' 'legacy/**/*.yml' 'legacy/**/*.yaml'
(no output, exit 1)
```

## 2. Diagnosis

### 2a. The time is all in the Build step

Step timings for the C++ jobs of seven recent runs (`gh run view <id> --json
jobs`, Configure / Build / Test in minutes):

| run | event | macOS Build | macOS Test | ubuntu Build | asan Build+Test | tsan Build+Test |
|---|---|---|---|---|---|---|
| 37229043133 | pull_request | 19.7 | 0.4 | 2.2 | 2.9+6.5 | 1.8+7.3 |
| 37226227284 | push | 25.1 | 0.4 | 2.4 | 3.5+10.2 | 2.2+8.5 |
| 37224723259 | merge_group | 22.2 | 0.4 | 2.5 | 2.7+5.5 | 1.6+4.6 |
| 37223176932 | pull_request | 21.5 | 0.4 | 1.9 | 3.7+7.7 | 2.2+8.4 |
| 37220905429 | merge_group | 22.2 | 0.5 | 3.0 | 3.7+10.3 | 2.2+8.4 |
| 37219150291 | pull_request | 21.3 | 0.5 | 3.8 | 3.7+10.3 | 1.9+7.5 |
| 37210876977 | push | 22.6 | 0.4 | 1.9 | 3.8+10.3 | 1.5+4.5 |

Configure is under 10 seconds everywhere; macOS Test is under 30 seconds.

### 2b. The two legs do the same work

Read from the logs (`gh run view 37229043133 --log --job 111514642326` for
macOS, `--job 111514642329` for ubuntu):

- Same commands: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  -DRASPUTIN_HARDENING=ON`, then `cmake --build build --parallel`, then
  `ctest`. Same targets: everything (`grep -c 'Building CXX object'` gives
  178 on both, `grep -c 'Linking CXX'` 73 on both). The job runs every
  registered test, so it needs every target.
- Compilers: AppleClang 21.0.0 on `macos-26-arm64` image; GCC 13.3.0 on
  ubuntu. Generator: Unix Makefiles on both (no `-G`).
- No compiler cache on either leg; nothing is restored between runs.
- Catch2 v3.6.0 is fetched by `FetchContent` (`tests/cpp/CMakeLists.txt:3`)
  and compiled from source on every run: 107 of the 178 translation units.
  pybind11 is not involved (the `_core` module is not built in this job).

### 2c. `--parallel` with no number starts every job at once

The build line passes no job count. A probe with a target that prints the
environment shows what make receives (CMake 4.4.3, GNU make, locally):

```
cmake --build b --parallel        ->  MAKEFLAGS=sj                        (unlimited)
cmake --build b --parallel 3      ->  MAKEFLAGS=s --jobserver-fds=3,4 -j  (limited)
CMAKE_BUILD_PARALLEL_LEVEL=3 cmake --build b --parallel  ->  MAKEFLAGS=sj (env ignored)
CMAKE_BUILD_PARALLEL_LEVEL=3 cmake --build b             ->  limited
```

So the build runs as `make -j` with no limit, and the logs show it: both legs
launch 108 compiler processes within four seconds of the Build step starting
(the 1st and 108th `Building CXX object` lines are at 19:38:54.3 and 19:38:57.9
on macOS, 19:38:42.7 and 19:38:47.0 on ubuntu).

### 2d. 108 compilers at once do not fit in 7 GB

Every translation unit was compiled locally with the same compiler (AppleClang
21.0.0) from this tree's `compile_commands.json`, four at a time, under
`/usr/bin/time -l` (peak memory and wall time per unit):

| group | units | sum of compile times | peak memory, median | peak memory, max | sum of peaks |
|---|---|---|---|---|---|
| Catch2 | 107 | 60 s | 132 MB | 205 MB | 13.3 GB |
| project | 71 | 142 s | 214 MB | 436 MB (`prop_refinement_edge_strip.cpp`) | 16.5 GB |

The first wave (the 107 Catch2 units plus the two library units) wants about
13 GB at once; the macOS runner has 7 GB, the ubuntu runner 16 GB. The macOS
runner pages; the ubuntu one does not.

The timings say the same. The first unit to finish in either leg is
`src/predicates/detria_exact.cpp` (0.7 s locally), whose library links at
19:39:04.7 on ubuntu (22 s after the step starts) and at 19:45:16.1 on macOS
(6 min 22 s; 6 min 41 s in run 37226227284). The whole first wave is about
60 s of compile time locally, so about 30 s on three cores even allowing the
M1 to be half as fast; overcommitting the CPUs cannot stretch 30 s into six
minutes, because overcommit does not add work. Paging can.

Not the cause: the compiler (the total compile time locally with AppleClang is
about 200 s, which three cores would get through in two to three minutes), the
build type (Release on both), or the tests (under 30 s).

**What this diagnosis has not shown**: the runner's memory was not observed
directly (the log has no memory trace, and a GitHub-hosted runner cannot be
inspected from here). §6's acceptance check is what confirms or refutes it.

## 3. Options

Expected times are estimates from §2d, to be replaced by the measured ones in
§6.

| option | expected macOS Build | cost and what it gives up |
|---|---|---|
| **A. A job count equal to the runner's CPUs** on every `cmake --build` line | 2-4 min (from 20-26) | None found. A few lines of YAML. |
| B. Compiler cache (ccache or sccache, kept with `actions/cache`) | after A: maybe 1-2 min on a hit, A's time on a miss | Cache keys and eviction to maintain; the repository's 10 GB cache limit shared with every job. GitHub scopes caches by branch: a run restores caches from its own branch and the default branch, and a pull request also from its base branch ("Dependency caching" reference, read 2026-10-04). Merge-queue branches are new each time, and h11 drops the push trigger, so nothing would write a cache on `master`; most PR and queue runs would miss. Not worth it at A's time. |
| C. Prebuilt Catch2 (Homebrew `catch2` and `find_package`, or Catch2's single-file amalgamation) | after A: saves maybe 30-60 s | A second way to get Catch2 into the build, version drift between Homebrew and the pinned v3.6.0, or one serial 20-30 s unit. Not worth it at A's time. |
| D. Build only the targets the tests need | none | The job runs every test, so it needs every target. |
| E. `macos-15-intel` (4 CPUs, 14 GB) | probably 3-5 min unchanged | Loses arm64, which is what Ola develops on; macos-15-intel is the last Intel image GitHub publishes, supported until August 2027 (as reported in the search results; not read from a GitHub page). Fixes the symptom, not the cause. |
| F. Larger paid runner (`macos-latest-xlarge`) | a few minutes | Money, for a cause A removes for free. |
| G. Take macOS out of the per-PR gate (nightly, or merge queue only) | 0 on the PR | AppleClang `-Werror` differences then surface after review or after merge, not on the PR. Belongs to h12 if wanted, and is not needed once A works. |

Ninja would also limit jobs by default (CPU count plus two), but it means
installing it on both runners and changes nothing A does not.

## 4. Recommendation

Option A, on all three `cmake --build` lines (the `cpp` matrix, which covers
its three legs, `sanitizers` and `tsan`), not only macOS: the ubuntu legs have the same unlimited `make -j` and are safe today only
because 16 GB happens to hold the current 178 units. The sanitizer builds in
particular use more memory per unit than Release.

```yaml
run: cmake --build build --parallel "$(getconf _NPROCESSORS_ONLN)"
```

- `getconf _NPROCESSORS_ONLN` is POSIX and prints the online CPU count on both
  Linux and macOS (locally it prints 10). It is shell, so it goes on the `run:`
  line; a workflow-level `env: CMAKE_BUILD_PARALLEL_LEVEL` cannot evaluate it,
  and §2c shows a bare `--parallel` ignores that variable anyway.
- If `getconf` ever printed nothing, the line would become a bare `--parallel`
  and fall back, silently, to an unlimited `make -j`. Nothing in the workflow
  fails on that; only §6 item 2 on this PR's run catches it, and a later
  runner image change would go unnoticed until the macOS job slowed again.
- A comment above the first build line says why the number is there (an
  unlimited `make -j` pages the 7 GB macOS runner), so nobody "simplifies" it
  back.
- The `tsan` step's `--target` list keeps its order; only the job count is
  added.
- The Python job's `pip install -e` builds one unit (`_core`) and is
  untouched.

**Constant: one job per CPU.** It assumes the largest unit's peak is well
under the runner's memory divided by its CPUs: 7 GB / 3 is about 2.3 GB per job
on macOS, 4 GB per job on ubuntu. Checked at this tree: largest unit 436 MB
(Release, AppleClang), so three at once peak at about 1.3 GB. It needs revisiting
only if a unit grows past about 2 GB.

**LOC:** about 5 lines of YAML (three build lines changed, a two-line comment),
no production code.

**After A, the next floor.** The longest C++ job becomes the asan+ubsan job,
8-14 minutes, most of it in its Test step (5.5-10.3 minutes): `ctest` runs the 994 tests one
at a time. `ctest --parallel` there is the next lever, but it needs a check
that no two tests share a temporary file, so it is not part of this increment
(question 2).

## 5. Exclusions

- No change to which jobs are required, or to when they run (h11, h12).
- No change to `timeout-minutes`; it guards a hung refinement loop, not the
  build.
- No compiler cache, no prebuilt dependencies, no runner change (§3 B, C, E, F).

## 6. How the change is checked

No unit test can see a workflow's job count, so there is no red step (as in
h10; the effect is observable only in CI). The check is the pull request's own
CI run, and it can fail. That run must show:

1. The macOS `Build` step under **6.0 minutes** (measured: 19.7-25.1 today).
   If it is not, the paging diagnosis of §2d is wrong, and the increment stops
   and goes back to `@architect` rather than trying option B.
2. In that job's log, the 108th `Building CXX object` line more than **10 s**
   after the 1st (today: 3.6 s on macOS, 4.3 s on ubuntu). Make prints the
   line when a compile starts. With three at a time the first 108 units, about
   60 s of compile time locally and more on three slower cores, cannot all
   start within 10 s, so this shows the limit took effect, and it is the only
   check that would catch an empty `getconf` result (§4).
3. The ubuntu legs' Build steps at most the top of §2a's range: `C++ core`
   (ubuntu) at most 3.8 minutes, asan+ubsan at most 3.8, tsan at most 2.2.
4. Every leg green; `ctest` reports 100 % of 994 tests passed on both C++ core
   legs; the tsan loop exits 0.

`@reviewer` is read-only. Its handback after the PR's CI run carries the
measured numbers for items 1 to 4, and the spawner records them, verbatim, as
a review round in `## Review` below. The numbers exist only after the push, so
this is a round after the PR's CI run, not before it.

## 7. Questions for Ola

1. Apply the job count to all three C++ build lines, so the ubuntu and
   sanitizer builds are bounded too (default), or only to the macOS leg?
2. Add `ctest --parallel` to the asan+ubsan Test step, the next-longest job
   after this fix? Default: no, a separate small increment after this one,
   because it first needs a check that tests do not share temporary files.

## 8. Ola's rulings

Ola answered both with one reply: "1: all three. 2: yes." (2026-10-04, this
session's transcript). The main session had put the questions to him as:
(1) put the job limit on all three C++ build lines, ubuntu and the sanitizers
too, or only macOS, default all three; (2) should the sanitizer job's tests
run in parallel, default yes, as a separate small increment afterwards,
because it first needs a check that no two tests write the same temporary
file.

1. **Ruling: all three** build lines (`cpp` matrix, `sanitizers`, `tsan`) get
   the explicit job count.
2. **Ruling: yes**, the asan+ubsan tests run in parallel, as its own
   increment after this one; it is not part of h13.

## Review

**h13, review, round 1, 2026-10-04.** Range `d926644..a77879f` (the design, Ola's rulings, the workflow edit). Verdict: CHANGES REQUESTED, on prose only; the workflow edit `a77879f` stands unchanged. LOC: 0 production lines; workflow 3 changed and 2 comment lines (estimate about 5). Checked by running: a bare `--parallel` gives make `-j` with no limit and ignores `CMAKE_BUILD_PARALLEL_LEVEL`, the new form passes `-j10` and the tsan `--target` list still builds (CMake 4.4.3, GNU make 3.81 probe); all seven §2a timings match `gh run view --json jobs` to 0.1 min; run 37229043133 logs: 178 compiles and 73 links on both legs, the 1st and 108th compile 3.7 s apart on macOS and 4.3 s on ubuntu, the 109th at 19:45:35 on macOS, `libterrain_predicates.a` links at 22.0 s (ubuntu) against 6 min 22 s (macOS), 6 min 41 s in run 37226227284; runner specs as cited (ubuntu 4 CPU / 16 GB, macos-latest 3 CPU M1 / 7 GB), repository public; memory spot-checked: `prop_refinement_edge_strip.cpp` 432 MB, a Catch2 unit 140 MB, `detria_exact.cpp` 143 MB in 0.63 s; the edit matches §4 and ruling 1, the Test steps are byte-identical to master; merges cleanly on h11 (`79a442a`); `check_citations.py`'s at-risk `h10-merge-queue.md:183` read as a quoted record. No red step: acceptable (h10 precedent; the effect is observable only in CI). Blocking: (1) §1's quoted sentence "If the -j option is given without an argument, make will not limit…" is from the `make(1)` man page, not the GNU make manual's *Parallel Execution* page, which says "If there is nothing looking like an integer after the '-j' option, there is no limit on the number of job slots."; (2) §6 says `@reviewer` records numbers in `## Review`, but `@reviewer` cannot write. What the PR's CI run must show: macOS Build step under 6.0 min (else back to `@architect`); 108th compile more than N s after the 1st; ubuntu legs' Build at most 3.8 min, asan+ubsan at most 3.8, tsan at most 2.2; every leg green, ctest 100 % of 994 on both C++ core legs, tsan loop exit 0.

(Recorded as the spawner relayed it, with code formatting added to names and
"vs" written out; N was left open by the review and is set to 10 s in §6.)

**h13, review, round 2, 2026-10-04.** Range `a77879f..14d6adc` (the increment file only); whole branch `d926644..14d6adc`, 0 production lines. Verdict: APPROVED, before the first push. Round 1's blockers closed: §1 quotes gnu.org's *Parallel Execution* page word for word, and the older sentence is credited to `make(1)`; §6 has the spawner record the CI numbers after the PR's run. The `cmake(1)` quotation, both sentences and the option order, matches cmake.org. §4's empty-`getconf` point probed: `cmake --build b --parallel ""` runs `make -j` with no limit and exits 0. Run 37229043133: 1st to 108th compile 3.66 s (macOS) and 4.31 s (ubuntu); ctest 100 % of 994 on both C++ core legs; the 10 s threshold can fail and leaves margin. `check_citations.py` exits 0; `h10-merge-queue.md:183` is a quoted record. Suggestions, not taken before the push: write 3.7 s throughout (§6 item 2 says 3.6 s); call §6 item 2 the only check that tells an empty `getconf` apart from a wrong diagnosis.

