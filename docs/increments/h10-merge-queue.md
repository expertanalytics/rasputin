# Harness h10: merge through GitHub's merge queue on `master`

Status: **designed** (`@architect`, 2026-10-04). Mechanics, not a design
question. One workflow edit (about 2 lines, `@developer`), one settings change
(Ola's), one rules commit (wording in §4).

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
  anywhere; **no** path filters; the only job `if:` conditions test
  `matrix.python-version`. Nothing is skipped or empty under the new event.
- `actions/checkout@v4` checks out `github.sha`, which under `merge_group` is
  the queue's test commit; `fetch-depth: 0` (Governance gates) behaves as on a
  PR.
- The non-required jobs (thread sanitizer, the `unchecked` C++ leg) also run in
  the queue. The queue waits only for the required checks, so they cost runner
  time, not queue time. Left as is: a red one on a queue run is still worth
  seeing.
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
| Maximum group size | **1** | Each PR is tested alone (on top of those ahead), so a red queue run names its PR; volume is too low for batching to save anything. |
| Wait time to meet minimum | **0** (or the lowest offered) | Moot at minimum 1. |
| Only merge non-failing pull requests | **On** | Every merge commit on `master` is then one that passed by itself; with it off, a later passing entry can carry an earlier failing PR in. |
| Status check timeout | **60 minutes** | The C++ jobs cap at 30 minutes and the slowest (macOS) has taken over 20; 60 leaves room for runner queueing without letting a hung run hold the queue for hours. |

"Require branches to be up to date" may stay ticked; the queue supersedes it.
"Allow auto-merge" in the repository's general settings is not needed for
enqueueing a PR whose checks are green; see §4 for a PR whose checks are not.

## 4. Rules text that changes (one rules commit, after Ola's step)

`git grep -n -E 'update-branch|gh pr merge'` over `tools/`, `.claude/` and the
rules files finds **no script or hook that runs `gh pr update-branch` or
`gh pr merge`**; the update-then-merge loop lived only in the main session's
practice. The places that describe merging:

- **`CLAUDE.md` §4, CI.** After the `gh pr checks <pr>` block add: "`master`
  merges through a merge queue: `gh pr merge <pr>` enqueues the PR (the queue's
  method, merge commit, applies), and the queue tests it on top of the PRs
  ahead and merges it. No update-branch loop. Never pass `--admin`: it
  bypasses the queue."
- **`.claude/REQUIRED-READING.md`, approval boundary** (the "A fresh yes"
  paragraph) and **`docs/PRINCIPLES.md` E1**: after `gh pr merge` add "(which
  enqueues; one yes covers one enqueue)". `gh pr merge` on a PR whose checks
  have not passed turns on auto-merge into the queue (`gh pr merge --help`,
  gh 2.101.0); that is the same act and needs the same yes.
- **`.claude/hooks/guard_push.py`** (for `@developer`, optional, needs Ola's
  approval as a hook change): `gh pr update-branch` writes a merge commit to
  the PR's remote branch but is not in the guard's list; probed at 45acf22, the
  guard exits 0 with no ask on `gh pr update-branch 170` and asks on
  `gh pr merge 170`. With the queue the command has no routine use, which
  makes asking before it cheap.

Not changing: `tools/session_state.py` finds the last landing with
`git log --merges`; queue merges are merge commits, so it still finds them.

## 5. Tests

None are written: a workflow trigger and a repository setting have no unit to
test, and a suite cannot create a merge queue. What stands in: the h10 PR's own
CI (§2, the YAML and existing triggers), then the **first queued merge**,
watched by the main session: `gh pr checks <pr>` shows the seven required
checks on the `gh-readonly-queue/master/...` run, and `git log -1 --merges
origin/master` shows a merge commit for that PR. If the queue run never
starts, the trigger is missing from `master`; the fix is to untick Require
merge queue, not to bypass with `--admin`.
