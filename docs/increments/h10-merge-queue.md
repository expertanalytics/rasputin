# Harness h10: merge through GitHub's merge queue on `master`

Status: approved by `@reviewer` (round 2, 2026-10-04), pending green CI on the
PR; then Ola ticks Require merge queue, then the rules commit (§4). Mechanics, not a design question. One workflow edit (2 lines,
`@developer`, done at `d6b0b1c`), one guard change (`@tester` then
`@developer`, §4a), one settings change (Ola's), one rules commit (wording in
§4).

Why: master's branch protection requires a branch to be up to date before it
merges. Merging several PRs is therefore a serial loop: update branch, wait for
full CI (the macOS job alone took over 20 minutes on #163), merge, and the next
PR is behind again. A merge queue tests each PR on top of the PRs queued ahead
of it, on a temporary branch, and merges them in order without the loop.

## 1. Prior art

Tooling; no novelty claimed.

*Literature.* The idea is the "not rocket science rule" (Graydon Hoare, 2014:
keep a branch whose every commit passed its tests, by testing the merge result
before it lands), as run by bors for Rust. GitHub's form of it
(docs.github.com, "Managing a merge queue", read 2026-10-04): a workflow "must
use the `merge_group` event to trigger your GitHub Actions workflow when a pull
request is added to a merge queue"; without it, "the merge will fail as the
required status check will not be reported". The queue's temporary branches
carry "the special prefix `gh-readonly-queue/{base_branch}`", and the queue
"provides the same benefits as the **Require branches to be up to date before
merging** branch protection, but does not require a pull request author to
update their pull request branch".

*Legacy.* Nothing; the legacy tree has no CI of this kind.

```
$ git grep -n -iE 'merge_group|merge queue|update-branch' legacy-archive -- legacy
(no output, exit 1)
```

## 2. The workflow change (`@developer`)

In `.github/workflows/main.yaml`, `on:` gains:

```yaml
  merge_group:
    types: [checks_requested]
```

Nothing else in the file needs to change. Read in full at 45acf22:

- **Concurrency group** `${{ github.workflow }}-${{ github.ref }}`. Under
  `merge_group`, `github.ref` is the queue entry's own temporary branch
  (`refs/heads/gh-readonly-queue/master/pr-<n>-<sha>`), so queue runs never
  cancel each other or a PR's run. Unchanged: two queue merges landing close
  together still cancel the first `push` run on `master`, as two quick merges
  do today.
- **No** `github.event.pull_request.*`, `github.base_ref` or `github.head_ref`
  anywhere; **no** path filters; the only `if:` conditions are step-level and
  test `matrix.python-version`. Nothing is skipped or empty under the new event.
- `actions/checkout@v4` checks out `github.sha`, which under `merge_group` is
  the queue's test commit; `fetch-depth: 0` (Governance gates) behaves as on a
  PR.
- The non-required jobs (thread sanitizer, the `unchecked` C++ leg) also run in
  the queue. The queue waits only for the required checks, so they cost runner
  time, not queue time. Left as is: a red one on a queue run is still worth
  seeing.
- When a queue entry fails, the queue rebuilds the entries behind it on new
  temporary branches. Their old runs are not cancelled, because each new
  branch is a different `github.ref` and so a different concurrency group.
  They run to the end for nothing: runner time again, not queue time, since
  the queue no longer waits on them.
- The `push` trigger on `master` stays. With the queue it re-tests a commit the
  queue has already tested (the queue fast-forwards `master` to it), so it is
  redundant but cheap insurance, and minutes are free on a public repository.

**What the PR's own CI proves:** the YAML parses and the existing triggers
still run every required check. **What only the first queued merge proves:**
that `merge_group` fires, that all seven required checks report on the queue
branch under their exact names, and that the queue merges with a merge commit.
The h10 PR itself merges the old way (update and merge), because the queue is
switched on only after the trigger is on `master`.

## 3. Ola's step, after the workflow change is on `master`

Settings → Branches → the `master` rule → tick **Require merge queue**. The
order matters: switched on before the trigger is on `master`, every queued PR
waits for checks that never report and fails at the timeout.

| Setting | Value | Why |
|---|---|---|
| Merge method | **Merge commit** | `docs/increments/README.md`: increment PRs never squash; the red commit must stay ahead of the green one in history. |
| Build concurrency | **2** | Rarely more than two or three PRs are ready together, and each queue entry starts a full CI including a macOS job; 2 keeps macOS runners free for PR runs. |
| Minimum group size | **1** | A lone PR merges at once instead of waiting for company. |
| Maximum group size | **1** | `master` moves forward one PR at a time, so each merge commit on it is one PR. This does not change how entries are tested: merge limits do not combine `merge_group` builds (GitHub, "Managing a merge queue", Merge limits, read 2026-10-04), so every entry gets its own CI run whatever the group size. |
| Wait time to meet minimum | **0** (or the lowest offered) | Moot at minimum 1. |
| Only merge non-failing pull requests | **On** | Every merge commit on `master` is then one that passed by itself; with it off, a later passing entry can carry an earlier failing PR in. |
| Status check timeout | **60 minutes** | The C++ jobs cap at 30 minutes and the slowest (macOS) has taken over 20; 60 leaves room for runner queueing without letting a hung run hold the queue for hours. |

"Require branches to be up to date" may stay ticked; the queue supersedes it.

"Allow auto-merge" (repository settings, General) is **on**: Ola ticked it on
2026-10-04, before the queue, and the main session read back
`allow_auto_merge: true` from the repository API. It matters because gh
2.101.0 always sends an auto-merge request when the base branch has a queue
(`pkg/cmd/pr/merge/merge.go:298-304`, found by `@reviewer`), whether or not
the PR's checks have passed; with the setting off, that request could be
refused. With it on, that risk is closed. Check:
`gh api repos/expertanalytics/rasputin --jq .allow_auto_merge` prints `true`.

## 4. Rules text that changes (one rules commit, after Ola's step)

`git grep -n -E 'update-branch|gh pr merge'` over `tools/`, `.claude/` and the
rules files finds **no script or hook that runs `gh pr update-branch` or
`gh pr merge`**; the update-then-merge loop lived only in the main session's
practice. The places that describe merging:

- **`CLAUDE.md` §4, CI.** After the `gh pr checks <pr>` block add: "`master`
  merges through a merge queue: `gh pr merge <pr>` enqueues the PR (the queue's
  method, merge commit, applies), and the queue tests it on top of the PRs
  ahead and merges it. No update-branch loop. A method flag such as `--merge`
  only prints a warning (the queue's method wins); `-d`/`--delete-branch` is
  refused. Never pass `--admin`: it bypasses the queue."
- **`.claude/REQUIRED-READING.md`, approval boundary** (the "A fresh yes"
  paragraph) and **`docs/PRINCIPLES.md` E1**: after `gh pr merge` add "(which
  enqueues; one yes covers one enqueue)". `gh pr merge` on a PR whose checks
  have not passed turns on auto-merge into the queue (`gh pr merge --help`,
  gh 2.101.0); that is the same act and needs the same yes.
- **`.claude/REQUIRED-READING.md`, "The harness"**: the guard's list reads
  `gh pr create/merge/ready/edit/update-branch` once §4a has landed (not
  before: until then the sentence would be false).

The guard already asks before every `gh pr merge`, with or without `--auto`
(the verb is the third word either way;
`.claude/hooks/guard_push.py@d6b0b1c:148`), so no act that
enqueues a PR passes silently.

## 4a. The guard asks before `gh pr update-branch` (`@tester`, `@developer`)

Ola's ruling, 2026-10-04: "yes, add update-branch to the guard". The command
writes a merge commit (or, with `--rebase`, a force-push) to the PR's branch on
the remote, yet the guard passes it: probed at 45acf22, `guard_push.py` exits 0
with no ask on `gh pr update-branch 170` and asks on `gh pr merge 170`. With
the queue the command has no routine use, so asking before it costs nothing.

*Failing test first (`@tester`)*, in `tests/python/test_guard_push.py`, whose
`ASKED` table drives `test_an_asked_command_asks_while_attended` and
`test_an_asked_command_is_denied_and_queued_while_unattended`. New rows, each
with the existing reason `"gh pr changes a pull request"`:

- `"gh pr update-branch 1"`
- `"gh pr update-branch --rebase 1"`
- `"gh pr merge --auto 1"` (passes today; pins the claim above)

and in `PASSED`, so the change cannot widen to every `gh pr` verb:
`"gh pr checks 1"`, `"gh pr view 1"`. The two `update-branch` rows are red at
`d6b0b1c`; the rest are green.

*The change (`@developer`)*, in `.claude/hooks/guard_push.py`: add
`update-branch` to the verb set in both places that list it, the text pattern
(`.claude/hooks/guard_push.py@d6b0b1c:39`, `create|merge|ready|edit`, used for
a line `shell_scan` cannot read) and the parsed check
(`.claude/hooks/guard_push.py@d6b0b1c:148`, `("create", "merge", "ready",
"edit")`). The reason
text is unchanged. A hook edit; Ola's ruling above is the approval. The module
docstring's spec reference gains `h10 §4a`.

Not changing: `tools/session_state.py` finds the last landing with
`git log --merges`; queue merges are merge commits, so it still finds them.

## 5. Tests

The guard change has its tests (§4a). For the rest none are written: a
workflow trigger and a repository setting have no unit to test, and a suite
cannot create a merge queue. What stands in: the h10 PR's own
CI (§2, the YAML and existing triggers), then the **first queued merge**,
watched by the main session: `gh pr checks <pr>` shows the seven required
checks on the `gh-readonly-queue/master/...` run, and `git log -1 --merges
origin/master` shows a merge commit for that PR. If the queue run never
starts, the trigger is missing from `master`; the fix is to untick Require
merge queue, not to bypass with `--admin`.

## Review

**h10, review, round 1, 2026-10-04.** Range `45acf22..d6b0b1c` (the note and 2 YAML lines). Verdict: CHANGES REQUESTED. LOC: 0 production (2 workflow lines, estimate about 2). The trigger form, the seven required checks and their names, the concurrency group, the absence of pull-request-only expressions and path filters, checkout at `github.sha`, gh 2.101.0's queue behaviour and the guard finding (`gh pr update-branch` exits 0 without asking; `.claude/hooks/guard_push.py@45acf22:39`, `:148`) are checked and true. Blocking: (1) `docs/increments/06-cdt-viewer.md:966` cites `.github/workflows/main.yaml:59-83` for the `sanitizers` job, which the two inserted lines moved to 61-85; cite the job by name. (2) The §3 reason for "Maximum group size = 1" is false: merge limits do not combine `merge_group` builds (GitHub, "Managing a merge queue", Merge limits); the value stays, the reason is that master moves one PR at a time. Not pushed; no CI. The PR's CI must show all nine jobs green; only the first queued merge can show that a `merge_group` run starts, that the seven required checks report under the same names, that the queue lands a merge commit, that "Require branches to be up to date" left ticked does not block it, and whether `gh pr merge` works with "Allow auto-merge" off.

(Recorded verbatim except one edit: the guard's line 39, cited by bare file
name and line, is pinned to `45acf22`, because `tools/check_citations.py`
fails on an unpinned line citation into a rule file.)

**h10, review, round 2, 2026-10-04.** Range `d6b0b1c..35673c8` (round-1 revision, guard red `fad1843`, guard green `391578c`, merge of master `72d4608`). Verdict: APPROVED, pending green CI on the PR. LOC: 0 production (2 workflow lines, estimate about 2; about 3 hook lines). Round 1's two blockers are fixed: `06-cdt-viewer.md` names the `sanitizers` job instead of a line range, and §3 gives the true reason for maximum group size 1. `allow_auto_merge` reads `true`. Guard tests: 6 failed / 130 passed at `fad1843`, all "update-branch passed silently"; 136 passed at the head. Mutant check: removing `update-branch` from the parsed verb tuple (`.claude/hooks/guard_push.py@35673c8:151`) is killed; removing it from the text pattern (`.claude/hooks/guard_push.py@35673c8:40`) survives, because no row reaches the text fallback (non-blocking). The merge has no conflicts and differs from master only in the 5 h10 files. `check_citations.py` exits 0. Not pushed; no CI. The PR's CI covers the full suite, since the branch changes no C++ or package Python. `.claude/REQUIRED-READING.md@35673c8:137` omits `update-branch` until the §4 rules commit.
