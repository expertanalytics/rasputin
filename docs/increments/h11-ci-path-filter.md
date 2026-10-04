# Harness h11: skip the code jobs on prose-only pull requests

Status: approved by `@reviewer` (round 2, see Review), on condition that CI
is green after the push. Red step 21d675e, 8db5bfc and
8986083 (the install trap's copy, §6.1); green 8b5c4ee (`@developer`); the
branch has no CI run yet. Rules text for §7 written (`@architect`, 2026-10-04).
After the merge, Ola's settings step (§5, step 2) follows directly; the
standing yes asked for in §10 is answered.

Why: `.github/workflows/main.yaml` has no path filter, so a pull request that
changes only prose (for example #174: `ROADMAP.md` and one increment status
line) runs the whole matrix: C++ on ubuntu, macOS and the unchecked leg,
asan+ubsan, TSan, Python 3.12 to 3.14, and the governance gates. With the merge
queue (h10), every PR runs CI three times: on `pull_request`, on `merge_group`,
and on `push` to `master` after it lands.

## 0. Verdict on the proposal

The main session proposed (a) a first job that lists changed files with plain
`git diff` and gates the C++, sanitizer and Python jobs with a job-level `if:`,
governance always running; (b) dropping the `push: master` trigger.

**(a) is right in aim and wrong in one mechanism; amended below.** A job
skipped by `if:` does report success, but **5 of the 7 required checks are
matrix jobs**, and a skipped matrix job is never expanded: it reports under
its unexpanded name (for example ``Python ${{ matrix.python-version }}``), so
the required checks `Python 3.12`, `Python 3.13`, `Python 3.14`,
`C++ core (ubuntu-latest)` and `C++ core (macos-latest)` would stay "Expected,
waiting for status to be reported" and **every prose PR would be blocked**,
on `pull_request` and in the queue alike. The required set on `master` was
read with
`gh api repos/expertanalytics/rasputin/branches/master/protection/required_status_checks`
on 2026-10-04: those five, plus `Governance gates` and
`C++ sanitizers (asan+ubsan)`. The amendment: one final job, `CI result`,
that needs every gating job and becomes the only required check (§2.3). A
second, smaller defect: the proposed list of "code" paths misses prose files
the Python suite reads (§2.1).

**(b) is sound, on one condition that GitHub's documentation does not settle**:
that the queue lands the very commit it tested. No queued merge has happened
on this repository yet (`gh run list --event merge_group` returned `[]` on
2026-10-04), so this is checked on the first one before (b) takes effect
(§5, step 0). Nothing depends on push runs: no cache, artifact or badge in the
tree (`grep -rnE 'actions/cache|upload-artifact|download-artifact' .github`
and `grep -rn badge README.md docs` return nothing), and
`tools/session_state.py` reads git, not CI (`origin/master`, `git log --merges`).
Branch protection has `enforce_admins: true` and no force pushes, so `master`
moves only through merges that passed the required checks.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature* (all read 2026-10-04):

- GitHub, "Using conditions to control job execution": "A job that is skipped
  will report its status as 'Success'. It will not prevent a pull request from
  merging, even if it is a required check."
- GitHub, "Troubleshooting required status checks": with a workflow-level path
  filter, "associated checks stay in a 'Pending' state and block merging"; the
  advice is to "avoid requiring workflows that can be skipped", and to "use
  `always()` with `needs` for required checks that depend on other jobs". This
  is why the filter lives in jobs, not in `on: paths-ignore`.
- Skipped matrix jobs keep unexpanded names: not stated in GitHub's docs;
  reported in GitHub community discussions #26822 and #72708 and in
  wgtechlabs/build-flow-action issue #51 (a skipped job's check shown as
  `Container flow${{ ... matrix.image.image-name ... }}`). Treated as true;
  the h11 PR is where it is seen on this repository (§6, real-PR checks).
- The `merge_group` payload (octokit/webhooks schema,
  `merge_group/checks_requested`): `base_sha` is "The SHA of the merge
  group's parent commit", `head_sha` "The SHA of the merge group". GitHub,
  "Events that trigger workflows": under `merge_group`, `GITHUB_SHA` is the
  "SHA of the merge group"; under `pull_request`, `GITHUB_SHA` is "the last
  merge commit of the pull request merge branch".
- GitHub, "Managing a merge queue": GitHub "will merge all these changes into
  the `base_branch` once the checks required ... pass". It does **not** say
  that the landed commit is the tested one. h10 §2 states "the queue
  fast-forwards `master` to it" without a source; secondary sources say the
  same, which is why §5 checks it instead of citing it.
- The pattern itself (a path classifier job plus one aggregate required job)
  is common; third-party actions do it (`dorny/paths-filter`,
  `re-actors/alls-green`). h11 writes both in plain shell and stdlib Python
  instead: no third-party action gets read access to the run.

*Legacy.* Nothing; the legacy tree has no CI of this kind.

```
$ git grep -n -iE 'paths-ignore|paths-filter|alls-green|merge_group|workflow' legacy-archive -- legacy
(no output, exit 1)
```

## 2. The change

### 2.1 The classifier: `tools/ci_changes.py` (new, stdlib only)

Deny by default: **a path is prose only if it is on a short allowlist;
everything else is code.** A wrong "code" costs one full run; a wrong "prose"
lets a broken change land green.

```python
NOT_PROSE: frozenset[str]  # Markdown files the build or the suites read
def is_prose(path: str) -> bool
def needs_full_ci(paths: Sequence[str]) -> bool   # empty -> True
def changed_paths(base: str, head: str, repo: Path) -> list[str]  # raises on any git failure
def main(argv: list[str]) -> int   # always prints code=true|false, exits 0
def prose_reads(paths: Iterable[str], root: Path) -> list[str]   # §6, T5
```

`prose_reads` keeps the paths that lie inside `root`, turns each into its
`/`-separated path relative to `root`, and returns those for which `is_prose`
is true; paths outside `root` are dropped. It is pure (no file system access
beyond the path arithmetic) and is what the T5 hook calls (§6.1).

`is_prose(path)` is true when `path` ends in `.md` (case-sensitive), is not
in `NOT_PROSE`, and either has no `/` (a root file) or starts with `docs/`.
So `ROADMAP.md`, `INSTALL.md`, `testing.md`, `docs/increments/*.md`,
`docs/retrospectives/**` are prose; everything under `src/`, `include/`,
`bindings/`, `lib/`, `src_python/`, `tests/`, `tools/`, `.claude/`,
`.github/`, and every non-Markdown file anywhere (`LICENSE`, `.gitignore`,
`CMakeLists.txt`, `pyproject.toml`, `docs/increments/h4-probes/h4_probe.py`,
`docs/benchmarks/**/*.geojson`) is code.

`NOT_PROSE` holds the Markdown the code jobs read, found by grep on this tree:

| File | Read by |
|---|---|
| `README.md` | `pyproject.toml` (`readme = "README.md"`), so `pip install -e` |
| `CLAUDE.md` | `tests/python/test_spawn_rule_text.py`, `test_session_state.py` (rule-size budget) |
| `docs/increments/README.md` | `tests/python/test_session_state.py` (`_copy_real_rule_files`) |
| `docs/PRINCIPLES.md` | the same |

`.claude/**` is code wholesale: the suite reads `REQUIRED-READING.md`, the
agent files, the skills, `briefs/common.md` and `settings.json`, and copies
the hooks. The main session's list would have missed the two `docs/` rows:
a prose PR that pushed `docs/increments/README.md` past the session-state
budget would go green and break the next, unrelated code PR.

`changed_paths` runs `git diff --name-only -z --no-renames <base> <head>`.
`--no-renames` is load-bearing: with rename detection, moving
`src/a.cpp` to `docs/a.md` lists only `docs/a.md`, which is prose, and the
deletion of a source file would skip the build. `-z` keeps unusual file
names unquoted.

`main` takes `<base> <head>` and prints exactly one line, `code=true` or
`code=false`, for `>> "$GITHUB_OUTPUT"`. **Safe default: `code=true`**, with
one plain sentence on stderr saying why, whenever the base is empty (a
`workflow_dispatch` run), git fails (the base is not in the clone, a bad
SHA), or the diff is empty. It exits 0 in all of these, so a classifier
problem shows as a full run, not as a red job. Only a crash of the script
itself fails the `changes` job, and then `CI result` fails (§2.3): also safe.

### 2.2 The workflow: `.github/workflows/main.yaml`

```yaml
on:
  pull_request:
    branches: [ master ]
  merge_group:
    types: [checks_requested]
  workflow_dispatch:        # replaces push: master; a hand-run is a full run
```

A new first job, `changes` (name `Changed files`), on ubuntu: checkout with
`fetch-depth: 0` (the base commit must be in the clone; governance already
fetches the full history), then

```
python3 tools/ci_changes.py \
  "${{ github.event.pull_request.base.sha || github.event.merge_group.base_sha }}" \
  "${{ github.sha }}" >> "$GITHUB_OUTPUT"
```

with `outputs: code: ${{ steps.<id>.outputs.code }}`. Both bases are
ancestors of `github.sha`, and both can be older than the head's actual
parent (a `pull_request` base moved on since the event; a queue entry behind
another entry, if `base_sha` is `master` rather than the entry ahead), so the
diff can only be wider than the PR, never narrower. Wider is safe.

`cpp`, `sanitizers`, `tsan` and `python` gain

```yaml
    needs: changes
    if: needs.changes.outputs.code != 'false'
```

`!= 'false'`, not `== 'true'`: an empty or unexpected output runs the jobs.
`governance` is unchanged: no `needs`, no `if`, it runs on every event. It
uses no other job's artifacts (no job uploads any), so skipping the others
costs it nothing, and `check_citations.py` matters most on prose.

The `unchecked` C++ leg gains `continue-on-error: ${{ matrix.hardening == 'OFF' }}`,
and `tsan` stays out of `CI result`'s `needs`, so the merge gate stays exactly
what it is today: neither is required now (§8 asks whether they should be).

### 2.3 `CI result`, the one required check

```yaml
  result:
    name: CI result
    if: always()
    needs: [changes, cpp, sanitizers, governance, python]
    runs-on: ubuntu-latest
    steps:
      - name: Every gating job passed, or was skipped on a prose-only change
        env:
          CODE: ${{ needs.changes.outputs.code }}
          CHANGES: ${{ needs.changes.result }}
          GOVERNANCE: ${{ needs.governance.result }}
          OTHERS: ${{ needs.cpp.result }} ${{ needs.sanitizers.result }} ${{ needs.python.result }}
        run: |   # plain shell; @developer writes it
```

The rule the step enforces, in words: `changes` and `governance` must be
`success`. If `CODE` is exactly `false`, each of the others must be `success`
or `skipped`; otherwise each must be `success`. Anything else (`failure`,
`cancelled`, a `skipped` job on a code change) exits 1 and says which job
ended how. `if: always()` makes the job run, and so report, even when a
needed job failed; without it a failure upstream would skip `CI result`,
and a skipped job reports success.

## 3. Failure modes and their defaults

| Case | What happens |
|---|---|
| Base SHA missing, empty, or not in the clone | `code=true`, full run |
| Diff empty (base equals head) | `code=true`, full run |
| `workflow_dispatch` (no base) | `code=true`, full run |
| Rename from code to a prose path | both paths listed (`--no-renames`): code |
| `changes` job crashes | code jobs skipped, `CI result` fails on `CHANGES` |
| A code job fails or is cancelled | `CI result` fails |
| A new suite starts reading a prose file | caught by the suite's own read check (§6, T5) |
| `tools/ci_changes.py` missing or broken | every `pytest` session over `tests/python` ends red, saying so (§6.1) |
| A new job added to the workflow | T4 fails until it is in `CI result`'s `needs` or on the named exemption list |
| PR with a merge conflict | no `pull_request` run at all (GitHub's rule), as today |
| Concurrency | group unchanged; `changes` and `result` cancel with the rest |

## 4. LOC

`tools/ci_changes.py`: about 45 counted lines (actual: 52). `main.yaml`:
about 35 added, 2 removed (the `push` trigger) (actual: 52 added, 2 removed);
h10 counted workflow lines apart from production, as here. Tests are excluded.
Total about 104 (actuals from review round 1). Well under 700.

## 5. Sequence, and Ola's settings step

0. **Before (b) takes effect:** on the first queued merge (likely #174),
   check that the queue landed the tested commit:
   `gh run list --event merge_group --limit 5 --json headSha` must contain
   `git rev-parse origin/master` taken right after that merge. If it does
   not, keep `push: master`; (a) stands without (b).
1. The h11 PR changes `.github/` and `tools/`, so its own run is a full run:
   all seven current required checks report under their names, plus
   `CI result`. It merges like any other PR.
2. **Ola, straight after the merge**: Settings, Branches, the `master` rule,
   required status checks: remove the seven, add `CI result` (offered because
   the h11 run reported it), with "GitHub Actions" chosen as its expected
   source, so that a check of the same name posted by another app cannot
   satisfy it. Until this is done, a prose PR waits on the five
   matrix checks that never report; nothing breaks, it just waits.
   `Require branches to be up to date` and the merge queue stay as h10 left
   them.
3. The first prose PR afterwards is the proof (§6, real-PR checks).

## 6. Tests (`@tester`)

`tests/python/test_ci_changes.py`, loading the tool the way
`harness_fixtures.Tool` loads the other `tools/` scripts.

- **T1, `is_prose` table.** Prose: `ROADMAP.md`, `INSTALL.md`, `testing.md`,
  `docs/increments/23c-x.md`, `docs/retrospectives/2026-10-04-x.md`,
  `docs/GLOSSARY.md`. Code: every `NOT_PROSE` entry;
  `docs/increments/h4-probes/h4_probe.py`;
  `docs/benchmarks/2026-09-26/quarter.geojson`; `.claude/agents/tester.md`;
  `.claude/briefs/common.md`; `.github/pull_request_template.md`;
  `.github/workflows/main.yaml`; `lib/detria/README.md`;
  `tests/python/notes.md`; `src_python/tin_engine/cli.py`; `LICENSE`;
  `pyproject.toml`; `CMakeLists.txt`; `.gitignore`; `tools/ci_changes.py`;
  `docs/x.MD`; `Docs/x.md`.
- **T2, `needs_full_ci`.** Empty list: true. All prose: false. One code path
  among prose: true.
- **T3, end to end on a fixture repository** (a temporary `git init`, as in
  `harness_fixtures.make_repo`), running `main` and reading stdout:
  a prose-only commit gives `code=false`; a code commit `code=true`; a rename
  `src/a.cpp` to `docs/a.md` gives `code=true` (red without `--no-renames`);
  deleting `src/a.cpp` alone gives `code=true`; an unknown base SHA, an empty
  base, and base equal to head each give `code=true`, a non-empty stderr, and
  exit 0.
- **T4, the workflow's shape**, read as text (PyYAML is not installed, and
  this is not a reason to add it): every top-level job other than `changes`,
  `governance` and `result` carries `needs: changes` and
  `if: needs.changes.outputs.code != 'false'`; `result` has `if: always()`;
  its `needs` equals the set of jobs minus `result` and minus a named
  exemption tuple (`tsan`), so a job added later fails the test until someone
  decides; `on:` has no `push`. The test fails loudly if the file's layout
  stops matching the shapes it parses, rather than passing on nothing.
- **T5, the suite's own reads.** A new `tests/python/conftest.py` installs a
  `sys.addaudithook` for `open` events, keeps every path inside the
  repository root, and at session end fails if `is_prose` is true for any of
  them, naming the file and saying to add it to `NOT_PROSE`. The pure part
  (`prose_reads(paths, root) -> list[str]`) is unit-tested; the wiring is
  proved by a planted run (a scratch test that reads `ROADMAP.md`, run in a
  subprocess, must fail the session). Limit, stated in the conftest: a tool
  the suite runs as a subprocess against the real tree is not seen; no suite
  does that today (they run tools against fixture repositories).

The red step: T1 to T4 fail on the missing tool and the unchanged workflow;
T5's planted run fails because the hook does not exist.

### 6.1 Two points the red step raised

**1. Who writes the T5 hook: `@tester`.** `tests/python/conftest.py` is a test
file, so it is `@tester`'s, and `@developer` does not edit it; `@tester` wrote
it in 8db5bfc. `@developer`'s share of T5 is `prose_reads` in the tool
(§2.1). The hook **fails closed**: it imports the tool at session end, and if
the tool is missing or has no `prose_reads`, the session ends red with a
sentence naming `tools/ci_changes.py`. So between the red step and the green
one every `pytest` run over `tests/python` is red, locally and on every CI
leg, whatever it tested; that is the red step working, not a regression.

**2. The install trap's copy of the tree: leave prose out; `@tester` makes
the change.** `copy_source_tree` in `tests/python/test_hardening.py` (T6 of
increment 24, run only with `RASPUTIN_INSTALL_TRAP=1`, which CI's "Hardening
install trap" step sets on the Python 3.12 leg) lists the checkout with
`git ls-files -co --exclude-standard` and copies each file with
`shutil.copy2`, which opens it in the test process. The hook would see
`ROADMAP.md`, `INSTALL.md`, `testing.md` and every `docs/**/*.md`, and that CI
step would turn red as soon as the tool exists. The copy does not depend on
what those files say.

The fix: `copy_source_tree` skips every listed path for which the tool's
`is_prose` is true (loading the tool the way `test_ci_changes.py` does, through
`harness_fixtures.Tool`). The `NOT_PROSE` files, `README.md` among them, are
still copied, so the install still finds the readme `pyproject.toml` names.
The test's docstring says why prose is left out, and the install's failure
message adds one sentence: if the build now needs a file that is not copied,
that file belongs in `NOT_PROSE`.

Weighed against the other two routes:

- *Exempt the copy from the hook* (a marker or an allowlist in the conftest).
  Rejected: it adds an exemption mechanism to the one check T5 rests on, and
  the next test that wants a pass would reach for it.
- *Copy with a subprocess* (`git checkout-index --prefix`, `cp`), whose reads
  the hook does not see. Rejected: it passes only by routing around the hook.
- *The cost of skipping*, the point `@tester` raised: a build that silently
  needs a prose file would break only in the trap. That is a gain, not a
  loss. The installs run in a `pip` subprocess, which the hook cannot see
  (its stated limit), so today nothing anywhere would notice the build
  starting to need, say, `docs/x.md`, and a prose-only PR could then change
  that file with no build run at all. With prose left out of the copy, the
  trap is the one place that build dependency shows up, as a red install on
  the 3.12 leg. The message is a build error rather than T5's sentence, which
  is why the failure message above names `NOT_PROSE`.

What the rest of the suite reads, checked for this ruling: the suite was run
on this branch with a throwaway pytest plugin that records every file opened
in the process and applies §2.1's rule, using the main checkout's virtual
environment (no C++ build allowed in that run, and its extension is older than
the tree). It saw no prose read in the 3,458 tests that passed. The stale
extension kept 32 modules from being collected and made 32 tests fail and 75
error, so those were not seen; in the 32 modules a `.md` name appears only in
docstrings, in files under `tmp_path`, and in the trap copy above
(`git grep -n '\.md'` over them). The probe was shown able to fail by a
planted test reading `ROADMAP.md`, which it named. The conftest itself is the
authority: the green run on CI covers every module.

**Only a real PR can show** (the main session, with `gh pr checks`):

- On the h11 PR: all nine jobs plus `Changed files` and `CI result` green,
  `code=true` in the `Changed files` log.
- On the first prose PR after §5 step 2: the C++, sanitizer and Python checks
  show as skipped (their names unexpanded, which confirms the reason for
  `CI result`), `Governance gates` and `CI result` green, on both the
  `pull_request` and the `merge_group` run, and the PR lands.
- After any queued merge: no `push` run on `master`
  (`gh run list --event push --branch master --limit 1` shows nothing newer
  than the merge).

## 7. Rules text

- `CLAUDE.md` §4, CI: what `CI result` checks; that it is meant to be the one
  required check, with the `gh api` command that shows whether it is yet (so
  the text is not false between the merge and Ola's settings step, §5 step 2),
  and that a prose-only PR waits on the never-reported matrix checks until
  then; on a PR that changes only prose (as `tools/ci_changes.py` decides),
  the C++, sanitizer and Python jobs are skipped, so that PR's green says
  nothing about code. Eight lines.
- `testing.md`, the "Live today" block and the sanitizer bullet ("run on
  every PR"): "on every PR that changes code". (Aside, not h11's: the tsan
  bullet there still says planned; the job is live.)
- `.claude/agents/reviewer.md` §5, after "Local green is not green": on a
  prose-only PR the C++, sanitizer and Python checks show as skipped; that is
  the design, not a red, and `CI result` decides. Two lines added, one
  rewrapped (b2e0a7b). It already said to read `.github/` when a change
  touches the build; the new sentence keeps a skipped check from being read
  as red CI.
- h10 is a record and stays as written; its "the queue fast-forwards
  `master`" is checked by §5 step 0, not relied on.

## 8. Question for Ola (answered: §9, ruling 1)

Should the TSan job and the unchecked C++ leg block a merge? Today they do
not (they are not required checks). Default: keep it so; adding them later is
one line each (`tsan` into `CI result`'s `needs`, and dropping the
`continue-on-error`).

## 9. Ola's rulings

Ola answered all three below with one reply: "yes to all, go ahead"
(2026-10-04T17:38:45.628Z, UTC). The main session put the questions to him;
each ruling is the default he said yes to.

1. **Question (§8):** should the TSan job and the unchecked C++ build (the
   leg built with hardening off) block a merge? **Ruling: no.** Both stay out
   of the merge gate: `CI result` does not list `tsan` in its `needs`, and the
   unchecked leg keeps `continue-on-error` (§2.2, §2.3). Making either one
   gating later is one line each.
2. **Question (§5, step 2):** will Ola change the required checks on `master`
   straight after h11 merges? **Ruling: yes.** Right after the merge, Ola
   removes the seven required checks and adds the single `CI result` check.
3. **Question:** should the prose-only fast lane go further, as its own
   increment? **Ruling: yes**, as h12, designed after h11 merges (§10).

## 10. Follow-on: h12, a fast lane for prose-only changes

h12 is designed after h11 merges and reuses h11's file classifier
(`tools/ci_changes.py`). Nothing below is designed or ruled; it is the main
session's sketch, kept for the record only.

The sketch: a tool sorts a change into one of three tiers from its changed
files. *Bookkeeping* (ROADMAP rows, status lines, review records copied in,
citation fixes): checked by a tool, no `@reviewer` round, governance CI only.
*Design and docs*: `@reviewer` as now, governance CI only. *Rule files*:
unchanged, full review. Open for Ola when h12 is designed: a standing yes for
the bookkeeping tier.

**That question is answered.** Ola gave a standing yes for bookkeeping changes
in the main session of 2026-10-04 (transcript 86af816a), in his words:
"Standing yes to bookkeepings changes." h12's design takes it from there:
what counts as bookkeeping, and which acts the yes covers, are h12's to
define and Ola's to confirm. One rule h12 must settle with it: today a grant
does not outlive the session it was given in (`.claude/REQUIRED-READING.md`,
the approval boundary), so until h12 lands, the standing yes is not yet a
rule any session can act on.

## Review

**h11, review, round 1, 2026-10-04.** Range `master...2faa3c4` (merge base `d926644`). Verdict: CHANGES REQUESTED. LOC: `tools/ci_changes.py` 52 production lines (estimate about 45); workflow 52 added / 2 removed (estimate about 35), counted apart as h10 did; total about 104. Gates clean; `check_citations.py` lists `h10-merge-queue.md:183` as at risk, read as a quoted record, stays. Full suite 4603 passed, 15 skipped, coverage 98.47 % on a rebuilt `_core`. Red: 65 failed / 1 passed at each red commit; 66 passed at `8b5c4ee`. Seven planted workflow and classifier faults, each caught. The `CI result` step, run under `bash -eo pipefail` over 14 result combinations, passes only on all green, or on `code=false` with the three code jobs skipped. §5 step 0 confirmed: the `merge_group` run on `d926644` equals `origin/master`. Blocking: §7 said `.claude/agents/reviewer.md` §5 needs no change, but `b2e0a7b` changed it; §7 said the `CLAUDE.md` change is two lines, it is five; say what was written. Suggestions: word `CLAUDE.md` §4's "the one required check" so it is not false before Ola's settings step (§5 step 2); in §5 step 2, add `CI result` with "GitHub Actions" as its expected source; in §4, record the actual line counts next to the estimates. After the push, CI must show all nine code jobs green plus `Changed files` and `CI result`, `code=true` in the `Changed files` log, and the 3.12 leg's "Hardening install trap" step green.

Fixes for round 1 (`@architect`): §7 now says what was written in `CLAUDE.md` and `.claude/agents/reviewer.md`; `CLAUDE.md` §4 says `CI result` is meant to be the one required check and names the `gh api` command that shows whether it is yet (eight lines now, not five); §5 step 2 names the expected source; §4 carries the actual counts.

**h11, review, round 2, 2026-10-04.** Range `2faa3c4..95b68da` (one docs commit, `@architect`; `CLAUDE.md` +5/-2, the increment file +29/-9; 0 production lines). Verdict: APPROVED for the first push, on condition that CI is green after it. Round 1's blocking finding is fixed: §7's `reviewer.md` bullet matches `git show b2e0a7b`, and its "eight lines" matches `git diff d926644..95b68da -- CLAUDE.md`. Claims checked: `gh api repos/expertanalytics/rasputin/branches/master/protection/required_status_checks` lists the seven per-job checks and not `CI result`, all with `app_id` 15368 (GitHub Actions); a skipped matrix job posts no check under its per-leg names, so a prose-only PR waits on the five matrix checks §5 step 2 names; `.github/workflows/main.yaml:306-307` names the job `CI result`; §4's counts recomputed. `check_citations.py` exits 0; its one at-risk line, `h10-merge-queue.md:183`, is a quoted record and stays. Blocking: none. After the push: all nine code jobs green plus `Changed files` and `CI result`, `code=true` in the `Changed files` log, the 3.12 leg's "Hardening install trap" step green; red CI turns this into CHANGES REQUESTED.
