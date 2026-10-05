# Harness h15: a faster CI that proves the same things

Status: design, not yet reviewed; Ola's rulings on §9 open. Design only: no workflow change in this
branch. Mechanics, not a design question for the core: the plan (§7) is
three small workflow PRs and one test edit, each ordered after h14 (#181,
the asan+ubsan tests in parallel) and h12's PR A (option C2, the queue
skips a rebuild of a tree its PR run passed), because all of them edit
`.github/workflows/main.yaml`.

Ola, 2026-10-05: "Can we impove speed (important) without sacrifising
presision?" The one test every option here is judged by (§4): **can this
ever let a change reach master that a full run would have failed?** An
option that can is named as such and is not recommended.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature* (all read 2026-10-05):

- The "not rocket science rule" and GitHub's merge queue: h10 §1. The queue
  is what makes every master commit one that passed; nothing here changes
  what it waits for, only how long the checks take.
- GitHub, "Dependency caching reference": "Workflow runs can restore caches
  created in either the current branch or the default branch (usually
  `main`)." "If a workflow run is triggered for a pull request, it can also
  restore caches created in the base branch". "Workflow runs cannot restore
  caches created for child branches or sibling branches." "When a cache is
  created by a workflow run triggered on a pull request, the cache is
  created for the merge ref (`refs/pull/.../merge`). Because of this, the
  cache will have a limited scope and can only be restored by re-runs of
  the pull request." "By default, the limit is 10 GB per repository."
  "GitHub will remove any cache entries that have not been accessed in over
  7 days." "You cannot change the contents of an existing cache." "If the
  job completes successfully, the action automatically creates a new cache".
- GitHub, "Actions limits": Team plan, "60" concurrent jobs on standard
  runners, "5" macOS. The organisation is on the Team plan
  (`gh api orgs/expertanalytics --jq .plan.name` prints `team`), and the 60
  are shared by every repository in it.
- GitHub, "GitHub-hosted runners", public repositories: ubuntu-latest and
  ubuntu-24.04-arm "4" CPU, "16 GB"; macos-latest "3 (M1)", "7 GB"; "Use of
  the standard GitHub-hosted runners is free and unlimited on public
  repositories." "Larger runners": "when larger runners are in use, they
  will always be billed at the per-minute rate."
- ccache manual (ccache.dev, latest): "Ccache is designed to produce exactly
  the same compiler output as a normal compilation." Its one documented
  correctness gap is in *direct mode*: "header files that were used by the
  compiler are recorded, but header files that were **not** used, but would
  have been used if they existed, are not". *Preprocessor mode* hashes "the
  preprocessor output from running the compiler with `-E`" plus the command
  line, so it has no such gap. The default `compiler_check` is the
  compiler's "mtime and size".
- `ctest(1)`: `-I` runs "tests starting at number Start, ending at number
  End, and incrementing by Stride"; `--no-tests=error`: "Consider running no
  tests to be an error", which is not the command-line default
  ("`ignore` ... the default when running tests via the ctest command
  line").
- pytest-xdist, "Known limitations": "It is not possible to have tests that
  differ in order or their amount across workers" (xdist then stops with a
  collection error, so this fails loudly). The page says nothing about tests
  that depend on state another test left behind; §5 P2 deals with that.
- Distributing a test suite over workers in longest-first order is Graham's
  LPT rule (h14 §1).

*Legacy.* Nothing; the legacy tree has no CI configuration, cache or test
distribution.

```
$ git grep -n -iE 'ccache|sccache|xdist|actions/cache|upload-artifact' legacy-archive -- legacy
(no output, exit 1)
$ git ls-tree -r --name-only legacy-archive -- legacy | grep -iE '\.ya?ml$|\.github|travis'
(no output, exit 1)
```

## 2. What each job does today

Read from `.github/workflows/main.yaml` at `7d13f3f` and the job logs. Every
job starts from an empty runner; nothing is cached or shared between jobs or
runs (`gh api repos/expertanalytics/rasputin/actions/cache/usage` prints
`"active_caches_count":0`; the workflow has no `actions/cache`,
`upload-artifact` or `download-artifact`).

| job | runner | builds | runs | gates `CI result`? |
|---|---|---|---|---|
| Changed files | ubuntu | nothing | `tools/ci_changes.py` (h11; h12 C2 adds `tested=`) | yes |
| Governance gates | ubuntu | nothing | three stdlib scripts | yes |
| C++ core (ubuntu-latest) | ubuntu | all 178 C++ units, Release, hardening ON; 107 of them are Catch2, fetched by `FetchContent` and compiled every run (h13 §2b) | `ctest`, 994 tests, serial | yes |
| C++ core (macos-latest) | macOS | the same, AppleClang | the same | yes |
| C++ core (ubuntu-latest, unchecked) | ubuntu | the same, hardening OFF | the same | no (`continue-on-error`) |
| C++ sanitizers (asan+ubsan) | ubuntu | the same 178 units, Debug, ASan+UBSan | `ctest`, serial today; 4 at once after h14 | yes |
| C++ thread sanitizer | ubuntu | 20 refinement suites, Debug, TSan | the 20 binaries in a loop | no (h11 §9) |
| Python 3.12 / 3.13 / 3.14 | ubuntu | the `_core` extension, three times per leg (one editable install per step: `.[dev]`, `.[dev,codecs]`, `.[dev,codecs,viewer]`), each a full ~26-28 s rebuild because each install uses a fresh isolated build environment; the 3.12 leg builds it twice more in the install trap | `pytest` (4709 passed in run 37265663651; one process, with coverage), the codecs subset (436), the viewer subset (153); on 3.12 only: the install trap, mypy, ruff, ruff format | yes |
| CI result | ubuntu | nothing | reads the others' results (h11) | is the required check |

The three Python legs cannot share one extension build: the wheel is tagged
for one interpreter (`rasputin-0.2.0.dev0-cp312-cp312-linux_x86_64.whl` in
the 3.12 log of run 37265663651), so 3.13 and 3.14 each need their own.
The duplication is *within* a leg: three builds where one would do.

## 3. Measurements

### 3a. Per job, the seven code runs since h13 merged

`gh run view <id> --json jobs`, job start to end in minutes. The time the
queue and Ola wait for is when `CI result` finishes, not when the run ends:
the thread sanitizer job runs on after it and gates nothing.

| run | what | `CI result` done | asan+ubsan | Python 3.12 | 3.13 / 3.14 | C++ ubuntu / macOS | tsan |
|---|---|---|---|---|---|---|---|
| 37232172673 | PR #177 | 13.9 | **13.7** | 8.5 | 5.5 / 6.4 | 3.1 / 3.2 | 10.9 |
| 37233115840 | queue #177 | 10.8 | **10.6** | 7.1 | 7.0 / 6.3 | 3.0 / 2.0 | 10.7 |
| 37234445967 | PR #178 | 13.9 | **13.7** | 7.0 | 7.1 / 6.8 | 3.3 / 3.1 | 10.7 |
| 37235448420 | queue #178 | 8.0 | **7.8** | 5.9 | 5.5 / 5.9 | 2.4 / 3.1 | 11.7 |
| 37265663651 | PR #179 | 14.3 | **14.0** | 9.0 | 6.3 / 4.7 | 3.0 / 3.7 | 10.8 |
| 37267264466 | PR #179 | 9.6 | 7.5 | **9.3** | 6.8 / 6.2 | 3.4 / 2.0 | (running) |
| 37267283373 | PR #181, h14's change | 9.3 | 8.6 | **9.1** | 4.5 / 5.2 | 2.5 / 2.4 | (running) |

The job that sets `CI result`'s time is in bold. In these seven runs every
job starts within 18 s of the run's first job, so runners are not queueing today. Prose-only
runs since h11: PR #180 0.8 min, its queue run 0.4 min (37265667235,
37265931630).

### 3b. Inside the two jobs that set the pace

**asan+ubsan.** Build 113-196 s; Test 315-631 s. h14's run (37267283373)
ran the tests four at a time: Test 315 s, of which the one case ES9 in
`prop_refinement_edge_strip` took 314 s (it took 97-171 s alone in h14 §2b),
and the 994 per-test times summed to 1176 s against 333-622 s in the serial
runs of h14 §2b (parsed from the job log's `Test #N: ... sec` lines). So four at once
on this runner does not make each test as fast as one alone; the step is
bounded by total work, and ES9 sits at that bound. h14's own §6 check is
where that run is judged; here it only sets the next floor.

**Python 3.12** (the same seven runs, seconds): Install 33-50 (of which
~28 the extension build), Tests 164-279, install trap 38-59, codecs step
68-111 (of which ~26 a rebuild), viewer step 36-47 (of which ~26 a
rebuild), mypy 5-8, ruff and format under 5 each. Everything after Tests is
serial in the one job and unique to 3.12, which is why 3.12 usually finishes last
of the three: up to 3.9 min after the slower of 3.13 and 3.14 in §3a (once
0.1 min before it).

### 3c. Push to merge, now and after h14 and C2

With "enqueue right after push", the machine wait from push to merge is
the PR run's `CI result` plus the queue run's.

| state | PR run | queue run | push to merge | source |
|---|---|---|---|---|
| today (h13 merged) | 13.9-14.3 (asan) | 8.0-10.8 (asan) | **22-25** | #177: 13.9 + 10.8; #178: 13.9 + 8.0 |
| after h14 | 9.3 (Python 3.12) | ~6-9.3 | ~16-19 | run 37267283373, h14's one run so far; the queue figure is an estimate |
| after h14 and C2, queue tree equals PR tree | 9.3 | ~0.4-0.8 | **~10** | the queue does only `changes`, governance and `CI result`, like a prose run (§3a) |
| after h14 and C2, master moved first | 9.3 | ~6-9.3 | ~16-19 | C2 then runs everything (h12 §11) |

So after h14 and C2 the floor is the Python 3.12 job (5.9-9.3 min), with
asan+ubsan right behind it (8.6 min in h14's run; 7.5 min in the serial
run 37267264466 on a fast runner). **Shortening one of the
two saves at most the gap to the other**, about 1 minute; h15's saving
comes from shortening both.

## 4. The test, and what passes it

**Can this option ever let a change reach master that a full run would
have failed?** A "full run" is today's workflow, run cold on the queue's
commit. An option passes when every check a full run makes is still made,
on the same inputs, with the same pass/fail rule, and its savings come from
running those checks at the same time or from not redoing work whose output
is provably identical. An option fails when its saving rests on a check not
being made (a skipped test, a test run in fewer environments, a different
build) unless something else makes it.

Two kinds of saving pass by construction:

- **Running the same commands in more jobs at once**, or one test list
  divided among jobs whose union is provably the whole list. Job boundaries
  do not change what a command checks, provided every new job gates `CI result`.
  Its risk is a new job that gates nothing; §7 pins that with a test.
- **Reusing an output whose inputs are hashed in full.** ccache in
  preprocessor mode keys each object on the compiler's whole preprocessed
  input and command line (§1); a hit is the object that compile would give.

Savings that rest on "this already passed once" (C2, the content-keyed
skip of §5 O1) pass only if the earlier run's inputs equal this run's in
full and the tests are deterministic; h12 §6 states that cost for C2 and
Ola accepted it (h12 §10). h15 adds no saving of that kind.

## 5. Options

Minutes are per merge, against the baseline after h14 and C2 (§3c: PR run
9.3, queue run ~0.5 when C2 skips). Because the two slowest jobs are within
a minute of each other, options are worth little alone; the "with" column
is the saving in the package of §7.

| # | option | saves alone | saves in §7's package | passes §4? |
|---|---|---|---|---|
| P | Split each Python leg into a main-suite job and an extras job | ~0.7 (down to asan) | ~3.5-4 together with A | **yes** |
| A | asan+ubsan tests in 3 jobs by ctest stride, with ES9 split into 4 cases | ~0 (Python is the floor) | (counted with P) | **yes** |
| B | ccache on the C++ jobs, warmed from master, preprocessor mode | ~0 now | ~0.5-1 after P and A | **yes**, with the settings below |
| X | pytest-xdist (several pytest workers) | ~0 alone | ~1 more, only if Python is the floor after P | **no**, unless one serial leg stays |
| O1 | Content-keyed whole build, skip `ctest` on a hit | ~0 | ~0 | only with a complete key; dominated |
| O2 | Build the extension once, share it across Python legs | impossible across versions | — | — |
| O3 | sccache instead of ccache | as B | as B | as B, less documented |
| O4 | Other runners | unmeasured / money | — | depends (below) |
| O5 | Concurrency and cancel rules | 0 | 0 | yes |
| O6 | `-O1` sanitizer build | maybe 2-3 min | — | **no** |

### P. Split each Python leg in two (recommended)

Each matrix leg today runs, one after another, the main suite (164-279 s on
3.12) and everything after it (the codecs and viewer subsets, and on 3.12
the install trap and the static gates: 150-230 s more, §3b). Split it into
two matrix jobs per version, with the **same commands in the same order and
environment** as today:

- `Python <v>`: checkout, `setup-python`, `pip install -e ".[dev]"`,
  `pytest` (the coverage gate). Expected 3.5-5.6 min on 3.12.
- `Python <v>, extras`: checkout, `setup-python`, `pip install -e
  ".[dev,codecs]"` and the codecs subset, `pip install -e
  ".[dev,codecs,viewer]"` and the viewer subset; on 3.12 also the install
  trap, mypy, ruff, ruff format, in that environment, as today. Expected
  4.4-5.3 min on 3.12, less elsewhere.

What it gives up: nothing a run checks. The codecs step today runs right
after the main suite in the same runner; afterwards it runs after a fresh
install of `.[dev,codecs]`, which is what its own first line already does.
Each step keeps its exact `run:` text. Cost: three more jobs per run (§5 O5,
the job pool). Its own check: §7 PR 1.

A cheaper half was considered and dropped: installing the extras' packages
without the editable reinstall saves the two ~26 s rebuilds, but it stops
running `pip install -e ".[dev,codecs]"`, the command INSTALL.md gives
users, on 3.13 and 3.14. Not worth 50 s in a job that is not the floor.

### A. The asan+ubsan tests in three jobs (recommended, after h14)

After h14 the job is Build (113-196 s) plus a Test step bounded by total
work: 1176 s of per-test time on 4 CPUs, with the single case ES9 at 314 s
(§3b). Splitting the tests over three jobs divides the work; each job
builds the sanitized tree itself (runner time is free) and runs every third
test:

```yaml
name: C++ sanitizers (asan+ubsan, ${{ matrix.shard }} of 3)
strategy: { fail-fast: false, matrix: { shard: [1, 2, 3] } }
...
run: >
  ctest --test-dir build-san --output-on-failure --no-tests=error
  --parallel "$(getconf _NPROCESSORS_ONLN)" -I ${{ matrix.shard }},,3
```

- `-I k,,3` runs tests k, k+3, k+6, ... (ctest numbers them in a fixed
  order from the build's test list), so the three jobs together run every
  test exactly once: the union of the three strides is every number from 1
  to the count. `--no-tests=error` makes a job that matched nothing fail
  rather than pass (§1: the command-line default is to pass).
- **ES9 split into four Catch2 cases**, one per terrain and seed, each
  keeping the four tolerances and every assertion (`@tester`'s edit, in
  `tests/cpp/property/prop_refinement_edge_strip.cpp`). Without it the shard
  holding ES9 is bounded by ES9 alone, 97-314 s. The 16 combinations and the
  assertions stay the same; only how ctest can spread them changes. h14's
  `COST 100` on the target stays.
- Expected: per shard Build 113-196 s plus Test ~100-180 s (a third of the
  work; a quarter of ES9 is 24-43 s alone and ~80 s at the load of h14's
  run): about 4-6.5 min against
  8.6.

What it gives up: nothing a run checks; each test runs once, in the same
build, with the same flags. Costs: two more jobs per run; and it **touches
h12's C2**. C2's condition 4b looks for a check run named exactly
`C++ sanitizers (asan+ubsan)` (h12 §11), and matrix jobs get expanded names.
After A, C2 would always print `tested=false` (safe, but the saving is
gone), and h12's T7 test fails because the name is no longer in the
workflow. So A's PR changes `tools/ci_changes.py` to require every shard
named `C++ sanitizers (asan+ubsan, <k> of 3)` (k = 1, 2, 3), all `success`,
in the same check suite as `CI result`; `@tester` extends T7 and T8 first.

**Constant: 3 shards.** Assumes 994 tests and ~1176 s of per-test time under
load (h14's run 37267283373), ES9 split in four. Checked only at this
tree's test count. Revisit when a shard's Test step passes 4 minutes.

### B. ccache for the C++ jobs (optional; ask after P and A)

Catch2's 107 units never change, and most project units do not change from
one PR to the next. ccache stores compiled objects keyed by their inputs.

Settings that make it pass §4:

- **Preprocessor mode** (`direct_mode = false`): the key is the compiler's
  full preprocessed output plus the command line, so the documented
  direct-mode gap (a header that *would* have been found had it existed)
  does not apply (§1).
- **`compiler_check = content`**: the key includes the compiler binary's
  contents, not its timestamp (the default), so an image update that
  changes the compiler misses. The cache name also carries the runner's
  `ImageOS` and `ImageVersion` and the job name, so jobs with different
  flags never share a store.
- Every link and every test still runs; a hit replaces a compile, nothing
  else. A failed compile is not cached, so `-Werror` still fails.

Who warms it, given there is no `push` trigger (h11): a run can read caches
of its own ref and of master (§1). Queue refs are new each time, and PR
caches are scoped to the PR's merge ref, so without a writer on master
only a PR's own later runs ever hit (the docs promise re-runs; a later push
uses the same merge ref, which should also hit, to be seen). The writer is
a **separate workflow file, `ccache-warm.yaml`**, triggered on `push` to
master, which builds the five C++ configurations and saves their stores; it
runs no tests and gates nothing, so `main.yaml` keeps h11's no-push design.
The same file carries a **weekly cold run** (`schedule`, ccache off) of the
C++ jobs, which is the check that ccache has not been trusted wrongly: a
master that fails cold fails there.

Cache poisoning: master's scope is written only by runs on master, that
is, by code the queue already passed. A PR's run can write only its own
scope; a queue run reads master's scope and its own (empty) one, so **the
run that lands a change never reads anything a PR wrote**.

Expected: Build 100-196 s down to roughly 30-70 s on a hit, the same plus
~20 s for restore and save on a miss; unmeasured, since no build may run in
this design round. Store size per job is unmeasured; the repository limit
is 10 GB, entries unused for 7 days go (§1), and five stores per master
merge churn it, so each store gets a `max_size` and the first PR measures
it. Costs: a second workflow, an install step (`apt-get` / `brew`), and
cache upkeep. After A it saves maybe 0.5-1 min on the slowest shard;
**ask Ola after A's numbers** (§9, question 3).

### X. pytest-xdist (not now)

Several worker processes share the suite. Each test still runs once, but in
a different process and order than the serial run, so a test that fails
only because an earlier test in the same process left state behind (a
module-level cache, an environment change) can pass under xdist and fail
serially. That is a failure a full run can show and xdist can hide, so X
alone fails §4. It passes with one gating job that keeps the serial run;
after P the serial main suite is 3.5-5.6 min, near A's 4-6.5, so X buys
about a minute at most. Also a new dev dependency and a state audit of
114 test files (as h14 §2d did for C++). Revisit only if Python is the
floor after P and A.

### O1. The content-keyed build cache, with `ctest` skipped on a hit

The main session's proposal: key = hash of `include/`, `src/`, `bindings/`,
`CMakeLists.txt`, `tests/cpp/`, toolchain and runner image; restore the
extension and test binaries; on a hit skip or re-run `ctest`.

- **The key as proposed is incomplete**, and an incomplete key is exactly
  the failure §4 forbids: it misses `lib/` (the vendored detria), the
  workflow file (the configure flags live there: Release or Debug,
  hardening, sanitizer flags), and Catch2, which `FetchContent` takes by the
  movable tag `v3.6.0` from the network (`tests/cpp/CMakeLists.txt`). A safe
  key would be deny-by-default, as h11's classifier is: every tracked file
  except prose and the Python-only paths, plus `ImageOS`, `ImageVersion`
  and the compiler's version. The Catch2 tag cannot be hashed without
  fetching it.
- **Who warms it, and scope**: the same as B; without a master writer only
  a PR's own re-runs hit.
- **What it saves here: almost nothing.** A hit needs unchanged C++ inputs.
  In a PR, that is a later push that changed only Python or tools, and
  there the Python jobs are the floor, so skipping C++ saves no wall time.
  In the queue, C2 already skips an unchanged tree, with stronger evidence
  (GitHub's record of the PR run rather than a hand-made key).
- **Restore binaries and re-run `ctest`** saves only the build, which B
  does per unit with ccache's own key, and B re-links from source-matched
  objects.
- **Poisoning**: as B, master's scope is written only by master; PR runs
  cannot reach the queue's reads.

Dominated by B and C2; not recommended.

### O2. One extension build shared by the Python legs

Not possible across legs: each wheel is for one interpreter (§2). Within a
leg, P keeps today's three builds, by choice (P, last paragraph). An
artifact of the sanitized build shared by A's shards is possible, but the
binaries' size is unmeasured (Debug, ASan; 73 linked targets, h13 §2b), and a
build-then-test chain puts the upload and download on the critical path;
A builds in each shard instead.

### O3. sccache

Equivalent in aim to B; it stores each object as its own entry in the
Actions cache. ccache is preferred because its manual documents the mode
(preprocessor) whose key has no known gap, and one store per job is easier
to size against the 10 GB limit.

### O4. Runners

- `ubuntu-24.04-arm` is free and has the same nominal 4 CPUs and 16 GB
  (§1). Its speed for these jobs is unmeasured; moving the sanitizer or
  Python jobs there would drop x86-64 from them (macOS already tests
  arm64), which fails §4. Adding it gains no speed.
- Larger runners (more cores) would shorten the two floor jobs directly,
  but they are billed per minute even here (§1). Ola's call (§9, question 4).
- `macos-15-intel`: rejected in h13 §3 E; macOS is 2.0-3.7 min, not a floor.

### O5. Concurrency and cancel rules

Already right: a push cancels its PR's run in flight; queue entries have
their own refs and never cancel each other (h10 §2). Two things to watch,
neither a change: the old runs of queue entries rebuilt after a failure run
on (runner time, h10 §2), and **the organisation's 60 concurrent jobs are
shared** (§1). P and A take a code run from 11 jobs to 16; a PR run, a queue
run and a second PR run together would take 48. Every job started within
18 s of its run in §3a; §7's checks repeat that measurement.

### O6. The sanitizer build at `-O1`

The sanitizer job compiles at `-O0` (Debug). At `-O1` it would run faster,
but the optimiser can delete an invalid access that is never used before
the sanitizer's checks see it, so a defect the `-O0` build reports could
pass. Fails §4; not recommended unless an `-O0` run stays (then nothing is
saved).

### Smaller items, not worth a PR of their own

- `setup-python`'s pip download cache: ~5-10 s per Python job; no risk.
- Starting the code jobs without waiting for `Changed files`: ~10-15 s;
  more conditions in every job for it.
- `COVERAGE_CORE=sysmon` on 3.12 and 3.13: unmeasured.
- C1 (the full suite only in the queue): rejected in h12 §6, for moving
  failures past Ola's merge yes; not a speed option under h15's test.

## 6. Order with h14 and h12 C2

All three touch `.github/workflows/main.yaml`; h12's PR A also changes
`tools/ci_changes.py` and its tests, which A here changes again. So:

1. h14 (#181) merges first; A's shards keep its `--parallel` line and
   `COST 100`, and A's numbers are measured against h14's.
2. h12 PR A (C2) merges second; A must then update C2's job-name condition
   in the same PR that renames the job (§5 A).
3. Then h15's PRs, in §7's order, each branched from master after the one
   before it merged.

## 7. Plan: small PRs, in order

Production code (the `CLAUDE.md` §2 count): none in PR 1 and PR 3; about
15 lines of `tools/ci_changes.py` in PR 2. Workflow lines are counted
separately, as h13 did.

**PR 1, Python split (P).** `main.yaml`: the `python` matrix job keeps
install and the main `pytest`; a new matrix job `python-extras` (name
`Python <v>, extras`) takes the codecs and viewer steps, and on 3.12 the
install trap and the three static gates, each step's `run:` text unchanged;
`CI result` needs it and checks it like `python`. About 40 workflow lines.
`@tester` first: in `tests/python/test_ci_changes.py`, which already parses
the workflow (h11 T4/T7), a test that every job except `Changed files`,
`CI result` and the thread sanitizer is in `CI result`'s `needs:` and in its
shell's list, so a job that gates nothing fails the suite (green at master,
and it keeps every later split honest); and a test that a job
`python-extras` exists, holds the codecs and viewer steps, and is in that
list (red at master). That each moved step's `run:` text is unchanged is
`@reviewer`'s check, by `git diff` of the workflow, not a test that would
freeze the commands. **Acceptance, on the PR's
run**: each leg reports the same counts as before (main suite 4709 passed /
119 skipped, codecs 436 / 4, viewer 153 / 2, trap 1 passed, at the counts
of the base commit, which the PR's own diff does not change); the slowest
Python job at most 6.0 min; every job starts within 60 s of the run's first.

**PR 2, the asan+ubsan shards (A).** `@tester` first, two commits: ES9
split into four cases (the 16 combinations and their assertions unchanged;
`@reviewer` compares old and new case by case), and the C2 tests (h12 T7,
T8) for three shard names in one check suite. Then `@developer`: the
`sanitizers` matrix of §5 A, `CI result` unchanged in shape (a matrix job
reports one result), `tools/ci_changes.py`'s condition 4b for the three
names. About 10 workflow lines, ~15 tool lines. **Acceptance, on the PR's
run**: the three shards' `out of N` lines sum to the count that `ctest -N`
prints for the same build (994 at master today, more if the split adds
cases: the split replaces one ctest test by four, so 997); each shard's
Test step at most 4 min; `CI result` at most 6.5 min; and on the queue run,
C2's `tested=true` when the trees match (h12 §11 step 3), so the rename
did not end the skip.

**PR 3, ccache (B), only on Ola's yes after PR 2's numbers.** The settings
of §5 B on the five C++ jobs, `ccache -s` and `ccache -p` printed at the
end of each Build step; `ccache-warm.yaml` (push to master: build and save;
weekly: the C++ jobs cold). About 60 workflow lines. **Acceptance**: on a
PR run after a warm run on master, `ccache -s` shows hits on the Catch2
units and the log shows `direct_mode = false` and `compiler_check =
content`; every test count unchanged; the first weekly cold run green.

**Expected end state** (estimates until each PR's run replaces them):

| | PR run | queue run, C2 skips | queue run, master moved | push to merge |
|---|---|---|---|---|
| today | 13.9-14.3 | — | 8.0-10.8 | 22-25 |
| after h14 and C2 | 9.3 | ~0.5 | ~6-9.3 | ~10 / ~16-19 |
| after PR 1 and PR 2 | ~5-6.5 | ~0.5 | ~5-6.5 | **~6-7 / ~10-13** |
| after PR 3 | ~4.5-6 | ~0.5 | ~4.5-6 | ~5-6.5 / ~9-12 |

## 8. Exclusions

- No change to what any job builds or runs: no flags, no tests removed, no
  environment dropped (O4, O6, X).
- No change to which checks are required, to C2's conditions other than
  the job name, or to the thread sanitizer and unchecked jobs (they gate
  nothing and do not hold a merge).
- No paid runners without Ola's ruling.
- No content-keyed skip (O1).

## 9. Questions for Ola

1. **Keep every Python test suite running in one process, in the same
   order as today, and gain the time by splitting each Python job in two?**
   The alternative, several test processes at once (pytest-xdist), is about
   a minute faster but can hide a test that only fails after another test
   left something behind. Default: **yes, split; no xdist**.
2. **Run the sanitizer tests in three jobs instead of one, and split the
   one very long test case (ES9) into four smaller cases with the same
   checks?** It saves about 3 minutes per run and checks the same things;
   it needs a small change to the queue-skip tool (h12's C2), since the
   job's name changes. Default: **yes**.
3. **A compiler cache (ccache), kept warm by a small extra workflow on each
   merge to master and checked by a weekly run without the cache?** It
   saves perhaps half a minute to a minute once PR 1 and PR 2 have landed,
   at the cost of a second workflow file and cache upkeep. Default: **decide
   after PR 2's measured times**.
4. **Paid larger GitHub runners (more cores) for the two slowest jobs?**
   Billed per minute even on a public repository. Default: **no**.
5. **Compile the sanitizer build with light optimisation (-O1) to run it
   faster?** It can hide a memory error the current build reports.
   Default: **no**.
