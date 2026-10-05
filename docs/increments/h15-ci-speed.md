# Harness h15: a faster CI that proves the same things

Status: design approved by `@reviewer` in round 2; Ola ruled on §9 on 2026-10-05 (all five defaults, below). PR 1 (§5 P, the Python job split) merged in #183: red `8a8af9b`, green `83d0601`; code review round 4 approved it. PR 2 is now designed in `docs/increments/h17-ci-test-time.md` §5, without the ES9 split, and no longer bound to land after h12's PR A if Ola says yes to h17 §8 question 2; PR 3 is not started. Mechanics, not a design question for the core: the plan (§7) is
two small workflow PRs and an optional third, each touching
`.github/workflows/main.yaml`. The order rule is §6: PR 1 needs only
master; PR 2 renames the sanitizer job, so it lands after h12's PR A
(option C2, the queue skips a rebuild of a tree its PR run passed) and
updates C2's job-name condition in the same PR. *(Superseded by
`docs/increments/h17-ci-test-time.md` §5 if Ola says yes to its §8
question 2, since Ola's ruling 2 below put PR 2 after h12's change: PR 2
then no longer waits for h12's PR A; whichever merges second carries the
name change.)*

Ola, 2026-10-05: "Can we impove speed (important) without sacrifising
presision?" The one test every option here is judged by (§4): **can this
ever let a change reach master that a full run would have failed?** An
option that can is named as such and is not recommended. The main
session's brief for this round adds that today's CI breaks Ola's flow, so
the design leads with the change that most shortens the time from a push
to a green `CI result`, and that harness work pauses after this increment
and h16, so the first PR is kept small.

**In one paragraph.** A push now waits 9.1-9.3 min for a green `CI result`
(the two PR runs since h14 merged; queue runs 7.6-9.2). Two jobs set that
time: the sanitizer job (asan+ubsan) and the slowest Python job. The
sanitizer job's time depends on the runner it gets: 8.5-8.9 min on the
slower kind (three of the five runs since h14), 4.6-5.1 on the faster
(two); the slowest Python job took 7.3-9.1. A Python job finished last in
three runs, the sanitizer job in two (§3). So shortening one of them alone
saves time only on some runs. PR 1 splits each Python job in two (§5 P):
the smallest change, and the one that saves most alone (0-3.7 min, 1.3 on
average over the five runs). PR 2 splits the sanitizer tests over three
jobs (§5 A). Together they bring the wait to an estimated 5.2-6.1 min, 2.4-3.7
min less in every one of the five runs. Every test still runs once, in the
same build and environment.

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
  On writers: "Only these workflow triggers can create or overwrite caches
  in the default branch's scope: `push`, `workflow_dispatch`,
  `repository_dispatch`, `delete`, `registry_package`, `page_build`,
  `schedule`". On fallbacks: "When a key doesn't match directly, the action
  searches for keys prefixed with the restore key." "If there are multiple
  partial matches for a restore key, the action returns the most recently
  created cache." The search order is the key, then the restore keys, in
  the run's own branch, then the same two in the default branch.
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
  line, so it has no such gap. On `compiler_check`: the default is the
  compiler's "mtime and size"; `content` is "Hash the content of the
  compiler binary"; `string:value` is "Hash **value**", for instance "a
  compiler revision number or another string that the build system
  generates to identify the compiler". And on wrappers: "when the compiler
  (as seen by ccache) actually isn't the real compiler but another compiler
  wrapper — in that case, the default **mtime** method will hash the mtime
  and size of the other compiler wrapper, which means that ccache won't be
  able to detect a compiler upgrade." The sentence names `mtime`; `content`
  hashes the same file, so it has the same blind spot (§5 B).
- `ctest(1)`: `-I` runs "tests starting at number Start, ending at number
  End, and incrementing by Stride"; `--no-tests=error`: "Consider running no
  tests to be an error", which is not the command-line default
  ("`ignore` ... the default when running tests via the ctest command
  line").
- pytest-xdist, "Known limitations": "It is not possible to have tests that
  differ in order or their amount across workers" (xdist then stops with a
  collection error, so this fails loudly). The page says nothing about tests
  that depend on state another test left behind; §5 X deals with that.
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

Read from `.github/workflows/main.yaml` at `bc01cd8` and the job logs. Every
job starts from an empty runner; nothing is cached or shared between jobs or
runs (`gh api repos/expertanalytics/rasputin/actions/cache/usage` prints
`"active_caches_count":0`; the workflow has no `actions/cache`,
`upload-artifact` or `download-artifact`). Its triggers are `pull_request`,
`merge_group` and `workflow_dispatch`; it has no `push` and no `schedule`.

| job | runner | builds | runs | gates `CI result`? |
|---|---|---|---|---|
| Changed files | ubuntu | nothing | `tools/ci_changes.py` (h11; h12 C2 adds `tested=`) | yes |
| Governance gates | ubuntu | nothing | three stdlib scripts | yes |
| C++ core (ubuntu-latest) | ubuntu | all 178 C++ units, Release, hardening ON; 107 of them are Catch2, fetched by `FetchContent` and compiled every run (h13 §2b) | `ctest`, 994 tests, serial | yes |
| C++ core (macos-latest) | macOS | the same, AppleClang | the same | yes |
| C++ core (ubuntu-latest, unchecked) | ubuntu | the same, hardening OFF | the same | no (`continue-on-error`) |
| C++ sanitizers (asan+ubsan) | ubuntu | the same 178 units, Debug, ASan+UBSan | `ctest`, 4 tests at once (h14) | yes |
| C++ thread sanitizer | ubuntu | 20 refinement suites, Debug, TSan | the 20 binaries in a loop | no (h11 §9) |
| Python 3.12 / 3.13 / 3.14 | ubuntu | the `_core` extension, three times per leg (one editable install per step: `.[dev]`, `.[dev,codecs]`, `.[dev,codecs,viewer]`), each a full ~26-28 s rebuild because each install uses a fresh isolated build environment; the 3.12 leg builds it twice more in the install trap | `pytest` (4709 passed in run 37265663651; one process, with coverage), the codecs subset (436), the viewer subset (153); on 3.12 only: the install trap, mypy, ruff, ruff format | yes |
| CI result | ubuntu | nothing | reads the others' results (h11) | is the required check |

The Python job's steps, in order: checkout; `setup-python`; Install
(`python -m pip install --upgrade pip && python -m pip install -e
".[dev]"`); Tests (`pytest`); on 3.12 the hardening install trap; the
codecs step (its own `pip install -e ".[dev,codecs]"`, then a list of
suites); the viewer step (its own `pip install -e ".[dev,codecs,viewer]"`,
then two suites); on 3.12 mypy, ruff, ruff format.

The three Python legs cannot share one extension build: the wheel is tagged
for one interpreter (`rasputin-0.2.0.dev0-cp312-cp312-linux_x86_64.whl` in
the 3.12 log of run 37265663651, `cp313-cp313` in the 3.13 log), and
pybind11 cannot build one stable-ABI (abi3) wheel for all three, so 3.13
and 3.14 each need their own. The duplication is *within* a leg: three
builds where one would do.

## 3. Measurements

Every figure here is from `gh run view <id> --json jobs` (job and step
start and end times) and, for per-test times, the job log's `Test #N: ...
sec` lines (`gh run view <id> --log --job <job id>`). Minutes are from the
run's first job starting.

### 3a. Per job, the green code runs from h13's PR run on

From h13's PR run (37232172673, created before its queue run merged it) to
today, green runs only: the red run 37232282246 (PR #178, a failing test
fixed in its next push) is left out, since a red run's time says nothing
about how long a green one takes. The time the queue and Ola wait for is
when `CI result` finishes, not when the run ends: the thread sanitizer job
runs on after it and gates nothing.

| run | what | h14 in it? | `CI result` done | asan+ubsan | Python 3.12 | 3.13 / 3.14 | C++ ubuntu / macOS | tsan |
|---|---|---|---|---|---|---|---|---|
| 37232172673 | PR #177 | no | 13.9 | **13.7** | 8.5 | 5.5 / 6.4 | 3.1 / 3.2 | 10.9 |
| 37233115840 | queue #177 | no | 10.8 | **10.6** | 7.1 | 7.0 / 6.3 | 3.0 / 2.0 | 10.7 |
| 37234445967 | PR #178 | no | 13.9 | **13.7** | 7.0 | 7.1 / 6.8 | 3.3 / 3.1 | 10.7 |
| 37235448420 | queue #178 | no | 8.0 | **7.8** | 5.9 | 5.5 / 5.9 | 2.4 / 3.1 | 11.7 |
| 37265663651 | PR #179 | no | 14.3 | **14.0** | 9.0 | 6.3 / 4.7 | 3.0 / 3.7 | 10.8 |
| 37267264466 | PR #179 | no | 9.6 | 7.5 | **9.3** | 6.8 / 6.2 | 3.4 / 2.0 | 10.8 |
| 37267283373 | PR #181, h14's change | yes | 9.3 | 8.6 | **9.1** | 4.5 / 5.2 | 2.5 / 2.4 | 10.7 |
| 37268029719 | queue #181 | yes | 9.2 | 4.6 | **8.9** | 5.0 / 6.5 | 2.1 / 4.3 | 10.8 |
| 37268031183 | queue #179 | yes | 7.6 | 5.1 | **7.3** | 7.2 / 6.2 | 3.2 / 2.3 | 10.8 |
| 37276269474 | PR #182 | yes | 9.1 | **8.9** | 6.2 | 5.0 / 7.3 | 3.0 / 2.9 | 10.8 |
| 37277177951 | queue #182 | yes | 8.8 | **8.5** | 7.2 | 8.4 / 7.5 | 3.2 / 2.8 | 11.0 |

Durations are job start to end; the bold job is the gating job that
finishes last, counting its start offset (in 37277177951 the sanitizer job
ends at 8.7 and Python 3.13 at 8.6). In all eleven
runs every job starts within 20 s of the run's first job, so runners are
not queueing today. Prose-only runs since h11: PR #180 0.8 min, its queue
run 0.4 min (37265667235, 37265931630).

### 3b. Inside the jobs that set the pace, since h14

The five runs with h14's change (the last five rows above), seconds.

**asan+ubsan.** Two clusters, one per runner speed:

| run | Build | Test step | per-test sum | ES9 | sum ÷ 4 |
|---|---|---|---|---|---|
| 37267283373 | 189 | 315 | 1176 | 314.3 | 294 |
| 37276269474 | 204 | 314 | 1183 | 313.0 | 296 |
| 37277177951 | 188 | 311 | 1163 | 310.7 | 291 |
| 37268029719 | 115 | 151 | 558 | 150.8 | 140 |
| 37268031183 | 134 | 158 | 592 | 157.6 | 148 |

All 994 tests pass in each. Build takes about 1.5 times and Test about
twice as long in the first three as in the last two, so the runner, not
the change, sets the cluster; three of five were slow. In every run the
one case ES9 (in `prop_refinement_edge_strip`) runs for the whole Test
step, and the step is also within 10-21 s of its work bound, the per-test
sum divided by the four test slots. The step is bounded by both at once.

**Python** (seconds, ranges over the five runs):

| step | 3.12 | 3.13 | 3.14 |
|---|---|---|---|
| Install (of which ~26-28 the extension build) | 28-45 | 29-49 | 39-53 |
| Tests (`pytest`) | 194-271 | 142-293 | 146-242 |
| install trap | 38-59 | — | — |
| codecs step (of which ~26 a rebuild) | 65-109 | 62-111 | 82-108 |
| viewer step (of which ~26 a rebuild) | 34-46 | 27-42 | 34-42 |
| mypy; ruff and format | 4-8; under 5 each | — | — |

Everything after Tests runs serially in the one job; on 3.12 it adds
141-222 s to the job.

### 3c. Push to green, push to merge

With "enqueue right after push", the machine wait from push to merge is
the PR run's `CI result` plus the queue run's. Push to green, what breaks
Ola's flow, is the PR run alone.

| state | PR run | queue run | push to merge | source |
|---|---|---|---|---|
| before h14 | 9.6-14.3 | 8.0-10.8 | **22-25** | pairs #177 (13.9 + 10.8) and #178 (13.9 + 8.0); PR range from §3a's five PR runs without h14 |
| after h14, measured | 9.1-9.3 | 7.6-9.2 | **17.9-18.5** | pairs #181 (9.3 + 9.2) and #182 (9.1 + 8.8); #179's pair spans h14 (9.6 + 7.6) |
| after h14 and C2, queue tree equals PR tree | 9.1-9.3 | ~0.4-0.8 | ~9.5-10 | the queue does only `changes`, governance and `CI result`, like a prose run (§3a) |
| after h14 and C2, master moved first | 9.1-9.3 | 7.6-9.2 | ~17-18.5 | C2 then runs everything (h12 §11) |

C2 is h12's, not h15's; it is the larger saving for push to merge, and it
does nothing for push to green. Push to green is set by the later of two
jobs (§3b): the sanitizer job (8.5-8.9 min on the three slow runners,
4.6-5.1 on the two fast ones) and the slowest Python job (7.3-9.1). Python
3.12 finished last in three runs (9.1, 8.9, 7.3), the sanitizer job in two
(8.9, 8.5), and in those two the slowest Python job was 7.3 and 8.4.
**Shortening one of the two saves time only on the runs where the other
is not close behind**; h15's saving needs both.

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
  do not change what a command checks, provided every new job gates `CI
  result` and runs its commands on the same state as before; §5 P says what
  "the same state" means there. The risk is a new job that gates nothing;
  §7 pins that with a test.
- **Reusing an output whose inputs are hashed in full.** ccache in
  preprocessor mode keys each object on the compiler's whole preprocessed
  input and command line (§1); a hit is the object that compile would give,
  provided the compiler and the tools after it are also in the key (§5 B).

Savings that rest on "this already passed once" (C2, the content-keyed
skip of §5 O1) pass only if the earlier run's inputs equal this run's in
full and the tests are deterministic; h12 §6 states that cost for C2 and
Ola accepted it (h12 §10). h15 adds no saving of that kind.

## 5. Options

### How the savings are estimated

For each of the five runs since h14 (§3b), every job's finish is its start
offset plus its duration, less the steps it would no longer run; a new
job's finish is its start offset plus the steps it would run. `CI result`
finishes at the last gating job plus the run's measured tail (0.1-0.2
min). A sanitizer shard runs the same Build and a third of the Test step,
plus 10% for an uneven split. These are estimates from measured step
times; each PR's own runs replace them (§7).

| run | now | P alone | A alone | ES9 split, no shards | P + A |
|---|---|---|---|---|---|
| 37267283373 (slow asan) | 9.3 | 8.8 | 9.3 | at least 9.3 | 5.6 |
| 37276269474 (slow asan) | 9.1 | 9.1 | 7.5 | at least 8.8 | 5.8 |
| 37277177951 (slow asan) | 8.8 | 8.8 | 8.7 | at least 8.7 | 6.1 |
| 37268029719 (fast asan) | 9.2 | 5.5 | 9.2 | at least 9.2 | 5.5 |
| 37268031183 (fast asan) | 7.6 | 5.3 | 7.6 | at least 7.6 | 5.2 |
| saving | — | 0-3.7, mean 1.3 | 0-1.6, mean 0.3 | at most 0.3 | **2.4-3.7, mean 3.1** |

| # | option | saves alone (§3b runs) | passes §4? |
|---|---|---|---|
| P | Split each Python leg into a main-suite job and an extras job | 0-3.7 min, mean 1.3 | **yes** |
| A | asan+ubsan tests in 3 jobs by ctest stride, ES9 split into 4 cases | 0-1.6, mean 0.3; 2.4-3.7 together with P | **yes** |
| A0 | ES9 split into 4 cases, one sanitizer job | at most 0.3 | yes |
| B | ccache on the C++ jobs, warmed from master, preprocessor mode | 0-0.8 after P and A, one run in five | **yes**, with the settings below |
| X | pytest-xdist (several pytest workers) | ~1, only once Python is the floor | **no**, unless one serial leg stays |
| O1 | Content-keyed whole build, skip `ctest` on a hit | ~0 | only with a complete key; dominated |
| O2 | Build the extension once, share it across Python legs | impossible across versions | — |
| O3 | sccache instead of ccache | as B | as B, less documented |
| O4 | Other runners | unmeasured / money | depends (below) |
| O5 | Concurrency and cancel rules | 0 | yes |
| O6 | `-O1` sanitizer build | maybe 2-3 min | **no** |

### P. Split each Python leg in two (recommended; PR 1)

Each matrix leg today runs the main suite and then, in the same job, the
codecs and viewer subsets, and on 3.12 the install trap and the static
gates (§2). Split it into two matrix jobs per version. Each job runs a
part of today's step list, **in today's order, each step's `run:` text
unchanged**:

- `Python <v>`: checkout; `setup-python`; Install (`python -m pip install
  --upgrade pip && python -m pip install -e ".[dev]"`); Tests (`pytest`,
  the coverage gate). That is today's job up to and including Tests.
  Estimated 3.2-6.0 min from the run's start (§3b's Install plus Tests,
  plus ~15 s for the job's start and setup); the top of that range is
  3.13's slowest main suite, 293 s.
- `Python <v>, extras`: checkout; `setup-python`; Install (the same step,
  pip upgrade included); on 3.12 the hardening install trap; the codecs
  step; the viewer step; on 3.12 mypy, ruff, ruff format. That is today's
  job with the Tests step taken out. Estimated 3.5-4.7 min from the run's
  start; it is not the floor in any of the five runs.

**What each moved step sees, against today.** The codecs and viewer steps
and the static gates run after the same installs as today, in the same
order, so the environment they see is today's. The one thing missing
before them is the main `pytest` process. What it leaves behind is on
disk in the checkout (coverage data, `.pytest_cache`, `__pycache__`) and
in `/tmp` (pytest's temporary directories, in a fresh runner directory
per job either way); nothing the later steps read.

**Where the install trap runs.** In the extras job, right after Install:
the same place relative to every install as today, without the main suite
before it. It does not depend on the main suite: it builds its own venv
from `[build-system].requires` (`build_environment` in
`tests/python/test_hardening.py`) and installs from its own copy of the
source, made from `git ls-files -co --exclude-standard`
(`copy_source_tree`, same file): tracked files plus untracked files that
`.gitignore` does not exclude. The main suite's leftovers listed above
are all ignored (`.gitignore` names `.coverage`, `__pycache__` and
`build`; `.pytest_cache` ignores itself), so the copy is the same with or
without it. A file the main suite wrote into the checkout outside those
patterns would reach the copy today and not after the split; `@reviewer`
checks there is none on PR 1 (`git status --porcelain
--untracked-files=all` after `pytest` in a clean clone prints nothing).
Keeping the trap in the main 3.12 job instead was weighed: it adds 38-59
s to that job and up to 1.0 min to `CI result` after P and A (estimated
6.5 instead of 5.5 in run 37268029719), for no change in what the trap
sees.

What P gives up: nothing a run checks. Cost: three more jobs per run, and
one more `.[dev]` install per leg (runner time only; the job pool is §5
O5). Its own checks: §7 PR 1.

A cheaper half was considered and dropped: installing the extras'
packages without the editable reinstall saves the two ~26 s rebuilds, but
it stops running `pip install -e ".[dev,codecs]"`, the command
INSTALL.md gives users, on 3.13 and 3.14. Not worth 50 s in a job that is
not the floor.

### A. The asan+ubsan tests in three jobs (recommended; PR 2)

The job is Build (115-204 s) plus a Test step at its work bound, with ES9
running the whole step (§3b). Splitting the tests over three jobs divides
the work; each job builds the sanitized tree itself (runner time is free)
and runs every third test:

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
  to the count. That holds when every shard numbers the tests the same way:
  the same commit, configuration and Catch2 listing, which in declaration
  order is the case here. `--no-tests=error` makes a job that matched
  nothing fail rather than pass (§1: the command-line default is to pass).
- **ES9 split into four Catch2 cases**, one per terrain and seed, each
  keeping the four tolerances and every assertion (`@tester`'s edit, in
  `tests/cpp/property/prop_refinement_edge_strip.cpp`). ES9's `GENERATE`
  lines give 2 terrains × 2 seeds × 4 tolerances, and Catch2 runs the whole
  body once per combination, so the four cases run the same 16 bodies.
  Without the split the shard holding ES9 takes at least ES9's time alone,
  97-171 s in h14 §2b's serial runs, so that shard would end at about
  Build plus 171 s, up to 6.3 min, above P's Python floor. h14's `COST 100`
  on the target stays.
- Estimated per shard: Build 115-204 s plus Test 55-115 s (a third of the
  step, plus 10%; a quarter of ES9 is 38-79 s at the same load), so 3.2-5.7
  min from the run's start, against 4.8-9.0 today.

**A0, ES9 split without the shards**, was compared. It renames no job and
leaves C2 alone. But the Test step is already within 10-21 s of its work
bound (§3b: per-test sum ÷ 4 is 140-296 s against a step of 151-315 s), and
splitting ES9 does not reduce the work. So A0 shortens the step by at most
those 10-21 s, about 0.3 min, and only when the sanitizer job is last. The
four slots are the limit, not ES9; only more CPUs (more jobs) move it.

What A gives up: nothing a run checks; each test runs once, in the same
build, with the same flags. Costs: two more jobs per run; and it **touches
h12's C2**. C2's condition 4b looks for a check run named exactly
`C++ sanitizers (asan+ubsan)` (`SANITIZERS_CHECK` in h12's
`tools/ci_changes.py`), and matrix jobs get expanded names. After A, C2
would always print `tested=false` (safe, but the saving is gone), and
h12's T7 test fails because the name is no longer in the workflow. So A's
PR changes `tools/ci_changes.py` to require every shard named `C++
sanitizers (asan+ubsan, <k> of 3)` (k = 1, 2, 3), all `success`, in the
same check suite as `CI result`; `@tester` extends T7 and T8 first. If
the shard count changes later without the tool, C2 prints `tested=false`:
safe.

**Constant: 3 shards.** Assumes 994 tests, a per-test sum of 558-1183 s at
four tests at once, and ES9 split in four (§3b, the five runs since h14).
Checked only at this tree's test count. Two shards would put the slowest
shard at about 6.6 min, above P's Python floor; four would save ~0.5 min
on a slow runner, where Python's main suite (up to 6.0) is then usually
the floor anyway.
Revisit when a shard's Test step passes 2.5 min.

### B. ccache for the C++ jobs (optional; ask after P and A)

Catch2's 107 units never change, and most project units do not change from
one PR to the next. ccache stores compiled objects keyed by their inputs.
After P and A, the sanitizer shard's Build (115-204 s) is most of its
time, but Python's main suite is then the floor in four of the five runs
(estimate), so B saves under a minute, on about one run in five, and more
only if the Python floor drops.

Settings that make it pass §4:

- **Preprocessor mode** (`direct_mode = false`): the key is the compiler's
  full preprocessed output plus the command line, so the documented
  direct-mode gap (a header that *would* have been found had it existed)
  does not apply (§1).
- **The compiler's identity, from the real compiler.** `compiler_check =
  content` is not enough: it hashes the file ccache runs, which here is a
  driver or a stub, not the compiler proper. On ubuntu that is `g++`, which
  runs `cc1plus`; on macOS `/usr/bin/c++` is an `xcrun` stub. §1's wrapper
  sentence says so for `mtime`; `content` hashes the same file. So one
  step before Build computes an identity string, and ccache gets
  `compiler_check = string:<it>`: the SHA-256 of the real compiler binary
  (`$(c++ -print-prog-name=cc1plus)` on ubuntu, `$(xcrun -f clang++)` on
  macOS) and of `c++ -v`'s output.
- **What that still does not see**: the assembler gcc runs after
  `cc1plus` (binutils), and any other tool on the image between the
  compiler and the object. They are guarded by the image: the cache key
  is `ccache-<ImageOS>-<ImageVersion>-<job>-<sha>`, the two image values
  read from the runner's environment in a step and passed as that step's
  output, and the one restore key is `ccache-<ImageOS>-<ImageVersion>-<job>-`.
  A prefix match therefore never crosses `ImageOS`, `ImageVersion` or the
  job (whose name carries its flags); there is no shorter restore key. A
  new image starts cold.
- Every link and every test still runs; a hit replaces a compile, nothing
  else. A failed compile is not cached, so `-Werror` still fails.

**Who writes, and who reads.** A run can read caches of its own ref and of
master (§1). Queue refs are new each time, and PR caches are scoped to the
PR's merge ref, so without a writer on master only a PR's own re-runs
would hit. The design:

- `main.yaml` only reads: it uses `actions/cache/restore`, never
  `actions/cache` (which saves at the end of a green job) or
  `actions/cache/save`. So none of its runs writes a store, whatever the
  trigger: not a PR run, not a queue run, and not a hand-started
  `workflow_dispatch` run on master, which §1's list allows to write
  master's scope.
- A **separate workflow, `ccache-warm.yaml`**, is the one writer. Its
  `push` job (on push to master) builds the five C++ configurations and
  saves their stores with `actions/cache/save`, under `if:
  github.event_name == 'push'`; it runs no tests and gates nothing, so
  `main.yaml` keeps h11's no-push design. Its **weekly cold run**
  (`schedule`) builds and tests the C++ jobs with ccache not installed and
  has no cache step, so it cannot save; it is the check that ccache has not
  been trusted wrongly: a master that fails cold fails there. The file
  has no `workflow_dispatch` and none of §1's other writing triggers.

Cache poisoning: master's scope is then written only by the `push` job,
which builds master's code, that is, code the queue already passed. PR and
queue runs write nothing. **The run that lands a change reads only stores
built from master.**

Expected: Build 115-204 s down to roughly 30-70 s on a hit, the same plus
~20 s for restore on a miss; unmeasured, since no build may run in this
design round. Store size per job is unmeasured; the repository limit is
10 GB, entries unused for 7 days go (§1), and five stores per master merge
churn it, so each store gets a `max_size` and the PR measures it. Costs: a
second workflow, an install step (`apt-get` / `brew`), and cache upkeep.
**Ask Ola after PR 2's numbers** (§9, question 3).

### X. pytest-xdist (not now)

Several worker processes share the suite. Each test still runs once, but in
a different process and order than the serial run, so a test that fails
only because an earlier test in the same process left state behind (a
module-level cache, an environment change) can pass under xdist and fail
serially. That is a failure a full run can show and xdist can hide, so X
alone fails §4. It passes with one gating job that keeps the serial run;
after P and A the serial main suite is the floor in four of five runs
(up to 6.0 min, estimate), so X with one serial leg buys little: that leg is still
there. Also a new dev dependency and a state audit of the test files (as
h14 §2d did for C++). Revisit only if Python's main suite stays the floor
after P and A.

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
  and the compiler's identity. The Catch2 tag cannot be hashed without
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
- **Poisoning**: as B.

Dominated by B and C2; not recommended.

### O2. One extension build shared by the Python legs

Not possible across legs: each wheel is for one interpreter, and pybind11
cannot build one abi3 wheel (§2). Within a leg, P keeps today's three
builds and adds one, by choice (P). An artifact of the sanitized build
shared by A's shards is possible, but the binaries' size is unmeasured
(Debug, ASan; 73 linked targets, h13 §2b), and a build-then-test chain
puts the upload and download on the critical path; A builds in each shard
instead.

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
  and would end the slow-runner lottery of §3b, but they are billed per
  minute even here (§1). Ola's call (§9, question 4).
- `macos-15-intel`: rejected in h13 §3 E; macOS is 2.0-4.3 min, not a floor.

### O5. Concurrency and cancel rules

Already right: a push cancels its PR's run in flight; queue entries have
their own refs and never cancel each other (h10 §2). Two things to watch,
neither a change: the old runs of queue entries rebuilt after a failure run
on (runner time, h10 §2), and **the organisation's 60 concurrent jobs are
shared** (§1). P and A take a code run from 11 jobs to 16; a PR run, a queue
run and a second PR run together would take 48. Every job started within
20 s of its run in §3a; §7's checks repeat that measurement.

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

All of them touch `.github/workflows/main.yaml`; h12's PR A also changes
`tools/ci_changes.py` and its tests, which A here changes again. The rule:

1. h14 (#181) first: A keeps its `--parallel` line in the sanitizer Test
   step and `COST 100` on `prop_refinement_edge_strip`, and A's numbers are
   measured against h14's. (Check: `git grep -n 'COST 100' origin/master
   -- tests/cpp/CMakeLists.txt` prints a line.)
2. PR 1 (P) needs nothing else: it does not touch the sanitizer job or
   C2's names. If h12's PR A is in flight, whichever of the two merges
   second merges master first and resolves the workflow text.
3. PR 2 (A) after h12's PR A, and updates C2's job-name condition in the
   same PR that renames the job (§5 A). (Check: `git grep -n
   SANITIZERS_CHECK origin/master -- tools/ci_changes.py` prints a line
   once C2 has merged.)
   *Superseded by `docs/increments/h17-ci-test-time.md` §5 if Ola says yes
   to its §8 question 2: the order then no longer holds; whichever of
   PR 2 and h12's PR A merges second carries the job-name change, and
   h12's own test fails until it does.*
4. Each h15 PR branches from master after the h15 PR before it merged.

## 7. Plan: small PRs, in order

Production code (the `CLAUDE.md` §2 count): none in PR 1 and PR 3; about
15 lines of `tools/ci_changes.py` in PR 2. Workflow lines are counted
separately, as h13 did.

**PR 1, Python split (P).** The first and smallest PR, and the one that
saves most alone. `main.yaml`: the `python` matrix job keeps checkout,
`setup-python`, Install and Tests; a new matrix job `python-extras` (name
`Python <v>, extras`) has checkout, `setup-python`, Install, the install
trap (3.12), the codecs step, the viewer step, mypy, ruff and ruff format
(3.12), in that order, each step's `run:` text and `if:` unchanged; `CI
result` needs it and checks it like `python`. About 40 workflow lines.
`@tester` first: in `tests/python/test_ci_changes.py`, which already parses
the workflow (h11 T4/T7), a test that every job except `Changed files`,
`CI result` and the thread sanitizer is in `CI
result`'s `needs:` and in its shell's list, so a job that gates nothing
fails the suite (green at master, and it keeps every later split honest);
and a test that a job `python-extras` exists, holds the codecs and viewer
steps after an Install step, and is in that list (red at master). That
each moved step's `run:` text is unchanged, and that each job's steps are
in today's order, is `@reviewer`'s check, by `git diff` of the workflow,
not a test that would freeze the commands; so is the `git status` check
of §5 P. **Acceptance, on the PR's run**: each step reports the same
passed and skipped counts as the same step in the base commit's run (main
suite, codecs, viewer, trap); the slowest `Python <v>` job at most 6.5
min and the slowest extras job at most 5.0, from the run's start; every
job starts within 60 s of the run's first.

**PR 2, the asan+ubsan shards (A).** `@tester` first, two commits: ES9
split into four cases (the 16 combinations and their assertions unchanged;
`@reviewer` compares old and new case by case), and the C2 tests (h12 T7,
T8) for three shard names in one check suite. Then `@developer`: the
`sanitizers` matrix of §5 A, `CI result` unchanged in shape (a matrix job
reports one result), `tools/ci_changes.py`'s condition 4b for the three
names. About 10 workflow lines, ~15 tool lines. **Acceptance, on the PR's
run**: the three shards' `out of N` lines sum to the count that `ctest -N`
prints for the same build (the base count plus 3, since the split
replaces one ctest test by four); each shard's Test step at most 2.5 min;
`CI result` at most 6.5 min; and on the queue run, C2's `tested=true` when
the trees match (h12 §11 step 3), so the rename did not end the skip.

**PR 3, ccache (B), only on Ola's yes after PR 2's numbers.** The settings
of §5 B on the five C++ jobs, `ccache -s` and `ccache -p` printed at the
end of each Build step; `ccache-warm.yaml` (push to master: build and save;
weekly: the C++ jobs cold, no cache step). About 70 workflow lines.
`@tester` first: a test in `tests/python/test_ci_changes.py` that
`main.yaml` uses no `actions/cache@` and no `actions/cache/save`, that every
`restore-keys` line in either workflow carries both image values, and that
`ccache-warm.yaml`'s only save step is conditioned on `push`.
**Acceptance**: on a PR run after a warm run on master, `ccache -s` shows
hits on the Catch2 units and `ccache -p` shows `direct_mode = false` and
`compiler_check = string:...`; every test count unchanged; the first
weekly cold run green.

**Expected end state** (minutes; estimates until each PR's runs replace
them; "push to merge" is with C2 skipping / with master moved):

| | PR run (push to green) | queue run, C2 skips | queue run, master moved | push to merge |
|---|---|---|---|---|
| before h14, measured | 9.6-14.3 | — | 8.0-10.8 | 22-25 |
| after h14, measured | 9.1-9.3 | — | 7.6-9.2 | 17.9-18.5 |
| after h14 and C2 | 9.1-9.3 | ~0.5 | 7.6-9.2 | ~9.5-10 / ~17-18.5 |
| after PR 1 | 5.3-9.1 | ~0.5 | 5.3-9.1 | ~6-9.5 / ~11-18 |
| after PR 1 and PR 2 | 5.2-6.1 | ~0.5 | 5.2-6.1 | **~6-6.5 / ~10.5-12** |
| after PR 3 | 5.0-6.1 | ~0.5 | 5.0-6.1 | ~5.5-6.5 / ~10-12 |

## 8. Exclusions

- No change to what any job builds or runs: no flags, no tests removed, no
  environment dropped (O4, O6, X).
- No change to which checks are required, to C2's conditions other than
  the job name, or to the thread sanitizer and unchecked jobs (they gate
  nothing and do not hold a merge).
- No paid runners without Ola's ruling.
- No content-keyed skip (O1).

## 9. Questions for Ola

1. **Split each Python job in two, as the first PR?** One job runs the
   main test suite, the other the extra test sets and the checks that run
   after it today, in the same order, so every test still runs once, in one
   process, in today's order. Estimated from the five runs since h14: it
   alone cuts the wait for a green check by 0-3.7 min (1.3 on average); it
   saves nothing on runs where the sanitizer job happens to get a slow
   machine. The alternative, several test processes at once
   (pytest-xdist), can hide a test that only fails after another test
   left something behind. Default: **yes, split, as PR 1; no xdist**.
2. **Then run the sanitizer tests in three jobs instead of one, and split
   the one very long test case (ES9) into four smaller cases with the same
   checks?** On top of PR 1 it saves 2.7-3.3 min on the three of five runs
   where the sanitizer job got a slow machine, and nothing on the other
   two; with PR 1, the wait for a green check goes from 7.6-9.3 to an
   estimated 5.2-6.1 min. Splitting ES9 without the extra jobs saves at
   most 0.3 min. It needs a small change to the queue-skip tool (h12's
   C2), since the job's name changes, so it comes after that. Default:
   **yes, as PR 2**.
3. **A compiler cache (ccache), kept warm by a small extra workflow on each
   merge to master and checked by a weekly run without the cache?** After
   PR 1 and PR 2 it saves under a minute, on about one run in five,
   because the Python main suite is then usually the slowest job, at the cost of a
   second workflow file and cache upkeep. Default: **decide after PR 2's
   measured times**.
4. **Paid larger GitHub runners (more cores) for the slowest jobs?** Billed
   per minute even on a public repository. Default: **no**.
5. **Compile the sanitizer build with light optimisation (-O1) to run it
   faster?** It can hide a memory error the current build reports.
   Default: **no**.

## Ola's rulings

2026-10-05, Ola, verbatim: "Accept all five defaults." So: 1 yes, split the Python jobs as PR 1, no xdist; 2 yes, sanitizer shards and ES9 split as PR 2, after h12's queue-skip change; 3 decide on ccache after PR 2's measured times; 4 no paid larger runners; 5 no -O1 sanitizer build.

## Review

### Round 1: `@reviewer`, design, `fc7346a` on `7d13f3f`

`@reviewer`'s handback, word for word (its Result and Pinned sections):

> ## Result
>
> **Verdict: CHANGES REQUESTED**
>
> **What I reviewed:** commit `fc7346a` on top of `7d13f3f`. It adds one file, `/Users/skavhaug/projects/rasputin/.claude/worktrees/h15-ci/docs/increments/h15-ci-speed.md` (520 lines). It is a design only: 0 lines of production code and no workflow change, so the 700-line ceiling does not apply. CI: the branch has not been pushed, so there are no checks to read. `python3 tools/check_citations.py` passes ("All citations resolve, and none point into a file this branch edits"), with no at-risk lines to re-read.
>
> **What checks out** (each one measured from `gh run view --json jobs`, the job logs and the GitHub and ccache docs):
> - **Job times, §3a:** all seven rows match to 0.1 min. No job started more than 18 s after its run's first job. The prose runs took 0.8 and 0.4 min. "No cache in use" is true: `active_caches_count` is 0.
> - **Python 3.12 step times:** every range in §3b matches the logs.
> - **Sanitizer job (asan+ubsan) in h14's PR run 37267283373:** the 994 per-test times sum to 1176 s, ES9 took 314.31 s, and the Test step took 315 s.
> - **Today's push-to-merge time of 22-25 min** comes from the pairs #177 (13.9 + 10.8) and #178 (13.9 + 8.0).
> - **The wheel cannot be shared across Python versions:** the 3.12 log has the `cp312-cp312` wheel and the 3.13 log has `cp313-cp313`.
> - **Test counts in run 37265663651:** 4709 passed / 119 skipped, codecs 436/4, viewer 153/2, install trap 1 passed.
> - **Splitting ES9:** its `GENERATE` has 2 terrains × 2 seeds × 4 tolerances = 16 combinations. Catch2 re-runs the whole body for each combination, so four cases of one terrain and seed, each keeping the four tolerances, run exactly the same 16 bodies. `catch_discover_tests` makes one ctest test per case, so the count becomes 997.
> - **Sharding with `ctest -I k,,3`:** the three strides together cover 1..N exactly once for any N. If N is below 3, `--no-tests=error` makes the empty shard fail loudly; it never passes wrongly. This assumes every shard numbers the tests the same way, which holds for the same commit, the same configuration and discovery in declaration order.
> - **The queue-skip tool (option C2) after sharding:** h12's check 4b reads `SANITIZERS_CHECK = "C++ sanitizers (asan+ubsan)"` in `tools/ci_changes.py` on h12's branch. Requiring all three shard names, all `success`, in the same check suite as `CI result` keeps 4b's purpose. If the shard count is ever changed without updating the tool, the result is `tested=false`, which is safe. A matrix job still reports one combined result to `CI result`.
> - **GitHub's cache scoping:** the quotes match the docs. A queue run can read only its own new ref and master. A PR's cache is scoped to `refs/pull/N/merge`, which is a sibling of the queue ref, so the run that lands a change cannot read anything a PR wrote.
> - **The rejections:** O1's key really is missing `lib/` (the vendored detria), the workflow file (where the build flags live) and Catch2 (fetched by the movable tag `v3.6.0`, `GIT_SHALLOW`). The reasoning for rejecting `-O1`, moving jobs to arm64 and xdist alone is sound under §4's test.
>
> **Blocking issues**
>
> 1. **Since the design was written, the sanitizer job after h14 has two more measurements, and they overturn §3c's premise and the saving promised to Ola in question 2.** h14's two queue runs are 37268029719 (#181) and 37268031183 (#179):
>    - Sanitizer job: 4.6 and 5.1 min. Test step: 151 and 158 s. Per-test sum: 558 and 592 s. ES9: 150.8 and 157.6 s, which is the entire Test step.
>    - Python 3.12 job: 8.9 and 7.3 min.
>
>    So h14's PR run (8.6 min, 1176 s) was a slow runner, not the typical case. Two of three runs show the sanitizer job 2.2-4.3 min *behind* Python 3.12, not "within a minute". What has to change:
>    - Rewrite §3c and the "Shortening one of the two saves at most the gap" paragraph from all three runs.
>    - Put the measured queue times (7.6 and 9.2 min) in place of the "~6-9.3, estimate".
>    - Re-derive option P's "saves alone" (more than 0.7 min).
>    - **Add the variant "split ES9 into four cases, no sharding."** In all three runs ES9 *is* the critical path. With ES9 split and 4 tests at a time, the Test step falls toward per-test sum ÷ 4 ≈ 140-150 s at the queue runs' load, with no job rename and no change to the queue-skip tool. Sharding then buys only the remainder, mainly on slow runners. Compare the two variants, and make question 2's "about 3 minutes" a measured figure for whichever variant is recommended.
> 2. **"Same commands in the same order and environment as today" (§5 P) is not true as written.**
>    - Today the install trap runs right after the main `pytest`, *before* the codecs install. P and §7 PR 1 move it into the extras job after the codecs/viewer installs, and the order there is not stated.
>    - The extras job's step list leaves out `python -m pip install --upgrade pip`, which today runs before every later install.
>
>    Fix: state the extras job's exact step order, and keep the pip upgrade. Either keep the trap in its old position (after the main suite) or say why moving it changes nothing. The evidence for that exists: `copy_source_tree` and `build_environment` in `tests/python/test_hardening.py` build their own venv from tracked and non-ignored files.
> 3. **Option B's "passes §4 with the settings below" overstates what `compiler_check = content` gives.** The ccache manual says that when the compiler ccache sees "actually isn't the real compiler but another compiler wrapper", hashing it "won't be able to detect a compiler upgrade". `content` hashes only the driver: `g++`, not `cc1plus`, and on macOS `/usr/bin/c++` is an `xcrun` stub. What actually guards a toolchain change is `ImageVersion` in the cache name, and only if no `restore-keys` prefix falls back across image versions. The design is silent on `restore-keys`. Fix: say the restore never crosses `ImageOS`/`ImageVersion`, and/or use the command form (`%compiler% -v`, or `string:` of `xcrun -f clang++` and its version). Also add the docs' write rule, which h15 leaves out: "only these workflow triggers can create or overwrite caches in the default branch's scope: push, workflow_dispatch, ... schedule". `main.yaml` has `workflow_dispatch`, and the weekly `schedule` cold run must not save a store. Neither is a route for poisoning (both run master code), but §5 B's sentence "master's scope is written only by runs on master" should name them.
> 4. **Small factual fixes:**
>    - §3a "the seven code runs since h13 merged" includes h13's own PR run (#177, 37232172673, created before its queue run) and silently leaves out the red run 37232282246. Say "from h13's PR run on, green runs only".
>    - §3c's "today, PR run 13.9-14.3" leaves out 37267264466, a PR run at today's workflow that finished in 9.6 min (§3a lists it). Give the range as 9.6-14.3, or say why it is excluded.
>    - Line 5 and §6 describe h14 (#181) as yet to merge. It merged (`0e946ed` on origin/master). Write the order rule, not the state.
>
> ## Pinned or assumed beyond the design
>
> - I did not run a build (the brief forbids it).
> - Tests numbering the same in every shard is reasoned from Catch2's declaration-order listing, not run.
> - The extras-job time (4.4-5.3 min) looks high against the step sums (my rough sum is 2.9-4.5 min). This is non-blocking, because PR 1's acceptance run will measure it.
> - O2's "impossible" is true for pybind11, which cannot build one stable-ABI (abi3) wheel. The design could say that is the reason.

What changed in response (`@architect`, on `worktree-h15-ci` after merging
master at `bc01cd8`):

1. §3 is re-measured from all eleven green runs, including the five with
   h14's change; two of them (37276269474, 37277177951, PR #182) came after
   the review, and both had a slow sanitizer runner, so three of five runs
   were slow and two fast (§3b). §3c now gives measured queue times
   (7.6-9.2) and push-to-merge (17.9-18.5) after h14, and the "at most the
   gap" paragraph is replaced. P's saving alone is re-derived per run
   (§5's first table: 0-3.7, mean 1.3). The variant without shards is §5
   A0. Its figure differs from the review's: per-test sum ÷ 4 at the queue
   runs' load is 140-148 s against a Test step of 151-158 s there, so the
   split brings the step down 10-21 s at most, because the step is already
   at its work bound; A0 is therefore compared and not recommended.
   Question 2's figure is now derived per run from measured step times
   (2.7-3.3 min on top of PR 1, on the slow-runner runs).
2. §5 P lists both jobs' steps in order, each a part of today's list; the
   extras job keeps the Install step with its pip upgrade, and the trap
   runs right after it. Why the trap's result cannot change, and the
   `git status` check that closes the one gap, are in §5 P; keeping the
   trap in the main job was weighed and costs up to 1.0 min.
3. §5 B now uses `compiler_check = string:` of the real compiler's hash
   and `c++ -v`, names what that still misses (the assembler and other
   image tools), keys the cache on the image with a single restore key that
   includes it, quotes the docs' writer list (§1), and makes `main.yaml`
   restore-only, so `workflow_dispatch` cannot write, while the weekly
   `schedule` run has no cache step. The manual's wrapper sentence is
   about `mtime`; §1 quotes it as such and says why it applies to
   `content` too.
4. §3a's heading and note, §3c's PR range (9.6-14.3) and the order rule
   (line 5, §6) are fixed. O2 now names abi3.

The extras-job time is re-derived as 3.5-4.7 min (non-blocking note); it
had counted the moved steps on the slowest leg twice.

### Round 2: `@reviewer`, design, `fc7346a..6997557`

`@reviewer`'s verdict, word for word:

> **Verdict: APPROVED**
>
> **What I reviewed:** `fc7346a..6997557` on `worktree-h15-ci`. `f25eaa5` merges master at `bc01cd8`. It has no conflicts, and its tree is the same as `git merge-tree --write-tree fc7346a bc01cd8`. `6997557` is the rework: `docs/increments/h15-ci-speed.md` +533 / −247. Against master the branch changes only that one file (806 lines). It is a design only, with 0 lines of production code and no workflow change, so the 700-line ceiling does not apply. CI: the branch has not been pushed, so there are no checks to read. It changes prose only, so a PR run would run just the governance gates. `python3 tools/check_citations.py` passes ("All citations resolve, and none point into a file this branch edits"), with no at-risk lines.
>
> **Round 1's blocking issues: all four fixed.**
> 1. §3 is re-measured from all eleven green runs. I rebuilt every row of §3a from `gh run view --json jobs`: all match to 0.1 min, and every job starts within 20 s of its run's first job. §3c's ranges are measured: PR 9.1-9.3, queue 7.6-9.2, push to merge 17.9-18.5. The variant with no shards is §5 A0, with a figure. Question 2's figure is now derived run by run.
> 2. §5 P gives both jobs' steps in today's order. The extras job keeps `Install` with its pip upgrade, and the trap runs straight after it. The argument for the trap holds. `copy_source_tree` copies `git ls-files -co --exclude-standard`, and `build_environment` builds its own venv from `[build-system].requires` (`tests/python/test_hardening.py`). What the main suite leaves behind is all ignored: `.coverage`, `__pycache__` and `build` are in `.gitignore`, and `.pytest_cache` ignores itself. No main-suite test runs `pip install` into the job's environment. The `git status` check on PR 1 covers the one gap left.
> 3. §5 B now uses `compiler_check = string:`. It hashes the real compiler (`cc1plus` or `xcrun -f clang++`) and `c++ -v`. No job sets `CXX`, so `c++` is the compiler CMake uses. The cache key and its single restore key carry `ImageOS` and `ImageVersion`. `main.yaml` only restores, and the weekly cold run has no cache step. Every quoted ccache passage matches ccache.dev/manual/latest word for word: `content`, `mtime`, `string:value`, the wrapper sentence, preprocessor mode and the direct-mode gap. The manual's Caveats section names only the direct-mode gap. The GitHub quotes, including the list of triggers that can write, match "Dependency caching reference".
> 4. The factual fixes are in: §3a's heading, the PR range 9.6-14.3, the order rule in its place of the state (`git grep 'COST 100' origin/master` prints line 381; `SANITIZERS_CHECK` is not yet on master), and abi3 named in O2.
>
> **The rework's new claims, checked by running them:**
> - **Sanitizer job since h14 (§3b):** I parsed all 994 `Test #N` lines in each of the five logs. Build, Test step, per-test sum and ES9 match every row: 189/315/1176/314.3, 204/314/1183/312.9, 188/311/1163/310.7, 115/151/558/150.8 and 134/158/592/157.6. The Test step minus (per-test sum ÷ 4) is 21, 18, 20, 11 and 10 s. ES9 takes the whole Test step in every run.
> - **A0 saves at most about 0.3 min:** true. The four test slots (`--parallel "$(getconf _NPROCESSORS_ONLN)"`, 4 CPUs) are already about 91-96% busy, and splitting ES9 removes no work. Round 1's suggestion that the step would fall to sum ÷ 4 at no cost was wrong; the design's correction is right.
> - **Savings per run (§5's first table):** I rebuilt them with the stated method (finish = start offset + kept steps; shard = Build + Test ÷ 3 × 1.1; plus each run's measured tail). P alone: 0-3.7 min, mean 1.3. A alone: 0-1.6, mean 0.3. P + A: 2.4-3.7. Nineteen of the twenty cells match to 0.1 min (the exception is under Suggestions). §5 A's per-shard figures (55-115 s, 3.2-5.7 min), two shards at about 6.6 min and four shards saving about 0.5 min all reproduce. So do §5 P's ranges for the main job (3.2-6.0 min) and the extras job (3.5-4.7 min) and §7's end-state table.
> - **Python step ranges (§3b):** every cell matches, as does the 141-222 s that the steps after Tests add on 3.12.
>
> **Blocking issues:** none.
>
> **Suggestions (non-blocking):**
> - In run 37276269474, P + A comes out at 5.72 min by my count, not 5.8. The mean saving is then 3.2, not 3.1, and question 2's range is 2.7-3.4. Question 2's "nothing on the other two" is 0 and 0.1-0.2 min (37268031183: 5.35 → 5.15). Worth one line to fix, since it is a figure put to Ola.
> - §1, the search order: per GitHub's docs, the action tries the exact key, then keys that start with it, then the restore keys, then the same steps on the default branch. The design leaves out the keys-that-start-with-it step. It does no harm with these keys.
> - §4 and §5 B, "the whole preprocessed input and command line": preprocessor mode leaves include and define options (`-I`, `-include`, `-D`) out of the hash, on the grounds that their effect shows up in the preprocessed output. A point of precision only.
> - PR 3: GitHub now has a `cache-mode` workflow key. `cache-mode: read` on `main.yaml` would make GitHub enforce "only restores" instead of relying on how the workflow is written. `@tester`'s PR 3 test could check for it.
> - PR 3: the cache key names `<job>`. After A, the three shards build the same tree, so key on the build configuration, not the shard's job name. Otherwise `ccache-warm.yaml` must write one store per shard name.
> - §5 A, ES9 left whole in one shard: the figure uses h14's serial 97-171 s. On CI under load ES9 took 151-314 s. The conclusion (above Python's floor) only gets stronger.

### Round 3: `@reviewer`, code, PR 1, `eac65cc..83d0601`

`@reviewer`'s verdict, word for word:

> **Verdict: CHANGES REQUESTED**
>
> **What I reviewed:** `eac65cc..83d0601` on `worktree-h15-ci`. Red `8a8af9b` (`@tester`) changes only `tests/python/test_ci_changes.py` (+226/−17). Green `83d0601` (`@developer`) changes only `.github/workflows/main.yaml` (+30/−7). Against master `bc01cd8` the branch changes three files: the workflow, the test file and `docs/increments/h15-ci-speed.md`. CI: the branch has not been pushed and has no PR, so there are no checks to read.
>
> **Size:** 0 lines of production code under `CLAUDE.md` §2 (the change is tests and workflow only). The design counts workflow lines separately: 19 added and 3 removed, not counting blank or comment lines, which nets to 16. The raw diff is +30/−7. §7 estimated "about 40 workflow lines". The change came in under that, so no split seam applies.
>
> **Moved steps, checked by `git diff` and a text comparison of the job blocks, master against `83d0601`:**
> - The `python` job is byte-identical to master's up to and including `Tests` (`pytest`). Only the steps after `Tests` are gone.
> - `python-extras` (name `Python ${{ matrix.python-version }}, extras`) has the same `runs-on`, `needs: changes`, `if:` and matrix (3.12, 3.13, 3.14, `fail-fast: false`) as `python`. Its steps are checkout, `setup-python`, Install and then master's steps after `Tests`:
>   - Install has the same `run:` text, pip upgrade included.
>   - The trap, codecs, viewer, mypy, ruff and ruff format steps follow, in master's order.
>   - Every `run:` and `if:` is byte-identical to master's.
> - The only text change inside a moved step is the codecs comment. I accept the new wording. "The main suite (the python job) covers the absent half" is true: that job installs `.[dev]` only, so it runs the half where the codecs package is absent. The old "the step above" was wrong after the split, since the step above is the trap on 3.12 and Install elsewhere, and it was already loose before.
> - No other job changed.
>
> **`CI result` wiring:**
> - `needs:` includes `python-extras`.
> - `OTHERS` passes `needs.python-extras.result` into a fourth slot, labelled "Python extras" in the loop.
> - `gh api .../required_status_checks` lists only `CI result`, so the new job's checks need no change to branch protection, and gate merges through `CI result`.
>
> **Tests, run in throwaway clones with no C++:**
> - At `8a8af9b`, `tests/python/test_ci_changes.py` gives 6 failed and 68 passed. All six failures are the missing `python-extras` job.
> - At `83d0601`, 74 passed.
> - I planted four faults, and each one failed the suite:
>   - `python-extras` dropped from `needs:`: 2 failed.
>   - "Python extras" dropped from the loop: 3 failed, including "a failure in OTHERS slot 3 passed".
>   - Its result dropped from `OTHERS`: 5 failed.
>   - The job renamed: 8 failed.
> - No red-step scaffolding is left in the test file or the workflow.
> - `ruff check` and `ruff format --check` pass on the test file. `check_prohibited_deps.py` is clean.
> - `check_citations.py` resolves every citation. It lists two lines as at risk: `h10-merge-queue.md:183` and `h11-ci-path-filter.md:451`. Re-read as quotations, both are quoted review records of their own time, and they stay as they are.
>
> **Stash:** `git stash list` is empty and there is no `refs/stash` log, so the developer's stash left the shared stack clean.
>
> **Blocking issues:**
> 1. **A prose claim the change made false.** The status line of `docs/increments/h15-ci-speed.md` says "Design only: no workflow change in this branch" and "next PR 1's red step (`@tester`)". The branch now holds PR 1's red and green commits and a workflow change. Before the push, the status line should say:
>    - PR 1 (§5 P) is implemented on this branch: red `8a8af9b`, green `83d0601`.
>    - The result of this review round, and that the push waits for Ola.
>    - PR 2 and PR 3 are not started.
>
>    Record this round under `## Review`. There is no ROADMAP row for h15; no harness increment (h10 to h15) has one, so nothing there needs changing.
> 2. **The clean-clone `git status` check of §5 P has not been run.** §5 P requires it of `@reviewer` on PR 1: in a clean clone, after `pytest`, `git status --porcelain --untracked-files=all` must print nothing. I could not run it without building C++:
>    - this worktree's venv has no `tin_engine`;
>    - the main checkout's `_core` dates from 2026-09-29, and the bindings have changed since.
>
>    The PR run cannot stand in for it, because no CI step runs `git status`. The route: an agent allowed to build C++ (no `@perf` run and no other C++ build at the same time) runs this at `83d0601`: `git clone`, a Python 3.12 venv, `pip install -e ".[dev]"`, `pytest`, then the `git status` command above. The `git status` must print nothing. If it prints any path, the trap's input differs after the split, and the design must answer for that file before PR 1 goes ahead. The output then goes to `@reviewer` for one more short round.
>
> **After the push (§7 PR 1 acceptance):**
> - Each step reports the same passed and skipped counts as on the base commit's run: the main suite, the codecs step, the viewer step and the trap.
> - The slowest `Python <v>` job takes at most 6.5 min from the run's start, and the slowest extras job at most 5.0 min.
> - Every job starts within 60 s of the run's first job.
> - `CI result` is green.
>
> Red CI turns this into CHANGES REQUESTED.
>
> **Suggestions (non-blocking):**
> - The last line of the reworded codecs comment is 104 characters, against about 80 for its neighbours. Re-wrap it.
> - A step in the `python` job that fails when `git status --porcelain --untracked-files=all` prints anything after `pytest` would make blocking issue 2 a check on every run instead of a one-off. It adds a step, so it is a design change, for `@architect` to weigh.

### Round 4: `@reviewer`, code, PR 1, `eac65cc..c08263b`

`@reviewer`'s verdict, word for word:

> **Verdict: APPROVED**
>
> **What I reviewed:** `83d0601..c08263b` on `worktree-h15-ci`, plus the clean-clone check of §5 P at `83d0601`. `c08263b` changes only `docs/increments/h15-ci-speed.md` (+64/−2): the status line, and round 3 recorded under `## Review`. CI: the branch has not been pushed and has no PR (`gh pr list --head worktree-h15-ci --state all` and `git ls-remote origin worktree-h15-ci` both print nothing), so there are no checks to read.
>
> **Size:** unchanged from round 3. There are 0 lines of production code under `CLAUDE.md` §2. `c08263b` is prose only.
>
> **Blocking issue 1 (status line): fixed.** Each claim on the new status line holds against the branch:
> - PR 1 (§5 P) is implemented: red `8a8af9b`, green `83d0601`.
> - Round 3 asked for this line and for the clean-clone check.
> - The push waits for that check, a short review round and Ola.
> - PR 2 and PR 3 are not started: the only workflow change on the branch is `83d0601`'s Python split.
>
> No "design only" or "no workflow change" wording is left anywhere outside the quoted review records. There is no h15 row in `ROADMAP.md` or `docs/increments/README.md` that needs changing.
>
> **Blocking issue 2 (clean-clone `git status` check of §5 P): run, and it prints nothing.** Commands and output (the scratch directory is my session scratchpad):
> ```
> $ git clone --no-local /Users/skavhaug/projects/rasputin/.claude/worktrees/h15-ci h15-clone
> $ cd h15-clone && git checkout --detach 83d0601
> HEAD is now at 83d0601 h15 PR 1 green step: each Python CI job split into a main-suite job and an extras job (@developer)
> $ git rev-parse HEAD
> 83d0601bc4f0d23fbafa5be37b8df5b3503fe916
> $ python3.12 -m venv ../venv312          # uv's CPython 3.12.12
> $ ../venv312/bin/python -m pip install --upgrade pip
> $ ../venv312/bin/python -m pip install -e ".[dev]"
> Successfully installed ... rasputin-0.2.0.dev0 ...      (exit 0; builds tin_engine/_core.cpython-312-darwin.so)
> $ git status --porcelain --untracked-files=all          # after install, before pytest
> (no output, exit 0)
> $ pytest
> Required test coverage of 85% reached. Total coverage: 98.62%
> ========== 4896 passed, 118 skipped, 7 warnings in 288.92s (0:04:48) ===========
> (exit 0)
> $ git status --porcelain --untracked-files=all
> (no output, exit 0; 0 lines)
> ```
> To make sure this check could have failed, I also listed the ignored files (`git status --porcelain --ignored --untracked-files=all`). After filtering out `build/`, `__pycache__/`, `.pytest_cache/` and `.coverage`, nothing is left (grep exit 1). Those four are exactly the leftovers §5 P names, and `copy_source_tree` excludes them. So the install trap sees the same source copy with or without the main suite before it, and the argument of §5 P holds. The clone and the venv have been removed.
>
> **Red-step scaffolding:** none. `c08263b` touches no test or workflow file.
>
> **After the push (§7 PR 1 acceptance), unchanged from round 3:**
> - Each step reports the same passed and skipped counts as on the base commit's run: the main suite, the codecs step, the viewer step and the trap.
> - The slowest `Python <v>` job takes at most 6.5 min from the run's start, and the slowest extras job at most 5.0 min.
> - Every job starts within 60 s of the run's first job.
> - `CI result` is green.
>
> Red CI turns this into CHANGES REQUESTED.
>
> **Suggestions (non-blocking), carried from round 3 and still open:**
> - Re-wrap the 104-character last line of the codecs comment in `.github/workflows/main.yaml`.
> - A `git status --porcelain --untracked-files=all` step after `pytest` in the `python` job would turn this one-off check into a check on every run. That is a design change, for `@architect` to weigh.
