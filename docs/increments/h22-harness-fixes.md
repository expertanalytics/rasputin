# Harness h22: the private remote under a push grant, two false reds after a commit, and scratch copies of the private package

Status: **design**, 2026-10-17, on `worktree-h22-harness-fixes` off master
`db12c3e1`. One PR, by day: every production file it changes is governed.
Implements three fixes Ola approved on 2026-10-17 as items 4 to 6 of the
night's queue, each from a failure seen that night.

Acceptance: the harness suites (`pytest -m harness`) and the live checks of section 8, recorded in this file; no `@perf` run (no refine or mesh code).

Words used here. **Public / private repository**: `expertanalytics/rasputin`
(remote `origin`, merge queue) and `expertanalytics/rasputin-private` (remote
`private`, no merge queue). **Private worktree**: a worktree named `lsm-*`
under `.claude/worktrees/`, holding a branch of the private repository.
**Private-only commit**: a commit reachable from `refs/remotes/private/master`
and not from `refs/remotes/origin/master`. **Grant**: a branch Ola named with
`away.py --grant` (`docs/increments/h20-away-push-grants.md`). **Governed**: a
file `guard_governance.py` asks about; changing it needs Ola's prompt.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature.*

- git, `git-ls-files(1)`, `--unmerged` / `-u`: "Show information about
  unmerged files in the output", and `git-merge(1)`, *How conflicts are
  presented*: on a conflict the index records "up to three versions" of the
  path and `MERGE_HEAD` is kept. Fix 2 reads unmerged index entries, not
  `MERGE_HEAD`, because a stopped `cherry-pick`, `revert`, `rebase` or
  `pull` leaves the same unmerged entries without `MERGE_HEAD`.
- GitHub GraphQL, `PullRequest.statusCheckRollup` (as `gh pr view --json
  statusCheckRollup` returns it): a list of `CheckRun` items (`status`,
  `conclusion`) and `StatusContext` items (`state`). Fix 1's private merge
  reads it because the private master is not protected (section 4.1), so
  `gh pr merge` merges on red checks.
- scikit-build-core, *Editable installs* (scikit-build-core.readthedocs.io,
  configuration/editable, read 2026-10-17): the default redirect mode "uses a
  custom redirecting finder to combine the static CMake install dir with the
  original source code"; "Older interpreters emit only `.pth` files" (3.15
  adds `.start` files). Fix 3 rests on both halves: the compiled `_lsm` sits
  in site-packages, and the `.pth` file adds the worktree's `src` to
  `sys.path` (section 6.1).
- Python, `site` module: a `.pth` line that is a path is appended to
  `sys.path`; `sitecustomize` is imported after the `.pth` files are
  processed. The h16 design already relies on the second fact.

*Legacy.* Nothing: the legacy tree has no harness.

```
$ git grep -l -i -e 'push_grant' -e 'scratch_copy' -e 'gates_after_commit' legacy-archive -- legacy
(no output, exit 1)
```

*Earlier increments.* `h21-harness-fixes-1014.md` item 4 (section 5) designs
a pre hook that records each committing command's HEAD, which covers fix
2(a) and two more cases (a merge with nothing to merge; a fast-forward of
master). Its PR B is not built: `git log --oneline origin/master --
.claude/hooks/gates_after_commit.py` ends at h19 (`f9e4648f`). Fix 2(a) here
is the narrower index check, independent of that hook; h21 PR B, when built,
keeps it or replaces it and says which.

## 2. The ideas file

`../rasputin_scratch/ideas/ideas.md`, read 2026-10-17. No entry is taken up.
Overlapping, left out (same files, other changes): "Make the after-commit hook
formatter-clean"; "After-commit hook: keep the cause when two faults meet";
"Cap on runner and shell words per line in the push guard".

## 3. Scope

| Fix | What | Files |
|---|---|---|
| 1 | The three granted forms for the private repository, and the rule that a branch goes only to its own repository | `tools/push_grant.py`, `tools/away.py` |
| 2 | The after-commit hook: no gates on a tree with unmerged paths; a gate whose script is absent is skipped with one line | `.claude/hooks/gates_after_commit.py` |
| 3 | `tools/scratch_copy.py` for the private layout: `src/`, the `_lsm` extension, the `.pth` path entry | `tools/scratch_copy.py` |

Not in scope: `guard_push.py` (its reasons for the new forms are already
`PUSH` and `PR`, which `judge` accepts; checked by reading `publishes()` in
`.claude/hooks/guard_push.py`, lines 267-285: `git push` gives `PUSH`
whatever the remote, and `gh pr create` and `gh pr merge` give `PR` whatever
`-R` names), `record_review.py` (section 4.1), and the private repository's
own copy of the harness (section 9).

## 4. Fix 1: the private forms in `tools/push_grant.py`

### 4.1 Facts the design rests on (each checked 2026-10-17)

- The private worktrees are worktrees of the main checkout's repository:
  `git worktree list` lists `.claude/worktrees/lsm-*`, and the main
  checkout's remotes are `origin` and `private`. So `judge`'s check that the
  command runs in the harness state's repository passes for them unchanged.
- The verdict for a private branch is read where it is for a public one:
  `record_review.py` appends every `Review: <verdict> <sha>` line to
  `<git-common-dir>/harness/reviews.jsonl` of `harness_mode.ROOT` (the main
  checkout, since the settings run `$CLAUDE_PROJECT_DIR/.claude/hooks/...`),
  keyed by sha only. The common dir is shared by every worktree. What is
  wrong for a private branch is the walk, not the record: `review_of` walks
  `refs/remotes/origin/master..<head>` and `changes` accepts only merges of
  `origin/master`. For `lsm-melt-out`, `origin/master..` holds 94 commits and
  `private/master..` holds 1.
- `gh api repos/expertanalytics/rasputin-private/branches/master/protection`
  answers `Branch not protected` (404), and the repository has
  `allow_auto_merge: false`. So `gh pr merge --merge` merges at once, red
  checks or not. The guard checks the checks itself (4.4).
- `gh pr view lsm-lake-inflow -R expertanalytics/rasputin-private --json
  state,baseRefName,headRefOid,statusCheckRollup` took 0.74 s and returned
  one `CheckRun` with conclusion `SUCCESS`, state `OPEN`, base `master`.
- History is shared up to the split: `git merge-base origin/master
  private/master` is `a058bb6b` (2026-10-14). `private/master` holds 93
  private-only commits. `lsm-melt-out` and `lsm-ddd-routing` hold all 93;
  `worktree-h22-harness-fixes` and `design-terrain-radiation` hold none.

### 4.2 The forms

Two constants: `PRIVATE_REMOTE = "private"`, `PRIVATE_REPO =
"expertanalytics/rasputin-private"`. `Act` gains `remote: str`, `"origin"` or
`"private"`. Every rule of h20 §3.2 still holds for every form: one simple
command, nothing the shell could expand (`unexpandable`), the title alone in
single quotes, a 40-hex sha, each option once, the body only by
`--body-file`/`-F`, a branch name by `branch_name`.

| Form | Exactly |
|---|---|
| Private push | `git [-C <absolute dir>] push private <40-hex>:refs/heads/<branch>` |
| Private pull request | `gh pr create -R expertanalytics/rasputin-private --head <branch> [--base master] [--title '<t>'] [--body-file <f>] [--draft]`; `-R` also spelled `--repo`, `--repo=` |
| Private merge | `gh pr merge <branch> -R expertanalytics/rasputin-private --merge --match-head-commit <40-hex>`, options in any order |

The public forms are unchanged, and still take no `-R`.

Refused as not a granted act (`granted_act` returns None, so the guard asks
as before):

1. A push to any remote but `origin` and `private`; to `private` with
   `--force`, `-f`, `+<sha>`, a short sha, a bare `:<branch>` ref, two refs,
   or `master`.
2. `gh pr create` or `gh pr merge` with `-R`/`--repo` naming anything but
   `expertanalytics/rasputin-private`, `-R` twice, or `-R` with no value.
3. The private merge without `--merge`, or with `--squash`, `--rebase`,
   `--auto`, `--admin`, `-d`/`--delete-branch`, or a second positional; with
   a PR number for the branch; without `--match-head-commit`.
4. The public merge with `--merge` keeps h20's meaning (allowed, the queue's
   method wins); the private form without it is refused because, with no
   queue, gh would prompt for a method.

The act kinds are `push`, `pull request`, `enqueue` (public merge) and
`merge` (private merge). `merge` is once per head per window, as `enqueue`.
The `grants-used.jsonl` line gains `"remote"`.

### 4.3 A branch's own repository

`home(repo, head) -> tuple[str, str]`: `(remote, "")` with remote
`"origin"` or `"private"`, or `("", reason)` when the branch's repository
cannot be told. The check fails closed: whenever the private remote is set
up, a head counts as public only when the private side was read and the
head holds none of it.

Two counts: `git rev-list --count refs/remotes/private/master
^refs/remotes/origin/master` (N) and the same with `^<head>` added (M).
Whether the private remote is set up is read from `git config --get
remote.private.url`, run with a 10 s timeout, by its exit status, three ways
(git 2.55.0, measured 2026-10-17: a key that is there exits 0, a key that is
not exits 1, a `GIT_CONFIG_GLOBAL` file holding `[broken` exits 128 with
`fatal: bad config line 1`):

- exit 0: set;
- exit 1: not set;
- any other exit, or the timeout: `home` refuses, with `("", "git config
  could not be read (it exited <code>), so it cannot be told whether the
  private remote is set up; fix the file git names")`, or for the timeout
  `("", "git config did not answer within 10 s, so it cannot be told whether
  the private remote is set up")`. The `rev-list` counts are not run.

A read error is therefore never taken for "not set": a malformed
`~/.gitconfig` is not a governed file, so without this it would turn every
head public. That anchor cannot be removed under a grant
without a prompt: `git remote remove`/`rm` and `git config` writes ask
(`git_why` in `guard_push.py`, `REMOTE_WRITES` and the `config` branch), and
`.git/config` is governed (`GOVERNED_SUFFIXES` in `guard_governance.py`). The
ref is not such an anchor: `git branch -dr private/master` deletes it with no
prompt. So the ref's absence is never read as "no private history":

| `remote.private.url` | `refs/remotes/private/master` | Result |
|---|---|---|
| unreadable (exit other than 0 or 1, or timeout) | not read | the config reason above |
| set | present, N > 0 | `"private"` when M < N, else `"origin"` |
| set | missing | `("", "the private remote is set up but refs/remotes/private/master is missing, so it cannot be told whether this branch holds private commits; run git fetch private")` |
| set | present, N = 0 | `("", "refs/remotes/private/master holds no commit that origin/master lacks, so it cannot be told whether this branch holds private commits; run git fetch private")` |
| not set | missing | `"origin"` (a clone without the private remote, CI) |
| not set | present, N > 0 | as the first row (a ref left behind; it is still read) |
| not set | present, N = 0 | `"origin"` |

`rev-list` exiting non-zero for any other reason (a damaged ref, a timeout of
10 s) is `("", "the private and public histories could not be compared")`.

In `judge`, after the head is resolved and before the review: a reason from
`home` refuses every form, origin and private alike, with that reason.
Otherwise, when the remote differs from `act.remote`, refuse with

- `head <h7> holds commits of private/master, so it goes only to the private repository` (an origin form on a private branch: the leak refusal), or
- `head <h7> holds no commit of private/master, so it goes only to origin` (a private form on a public branch).

A private form in a clone without the private remote (the fourth row) is
refused by the second sentence.

Known limit, stated, not fixed: a branch forked before the split
(`a058bb6b`) and never merged with `private/master` counts as `"origin"`
(`design-terrain-radiation` today), so its own commits could go to origin
under a grant. Ola grants by name, and the guard cannot tell intent.

### 4.4 The review walk and the checks

`review_of(repo, head, records, remote="origin")` and `changes(repo, commit,
remote="origin")` walk and test against `refs/remotes/<remote>/master`;
`judge` passes `act.remote`. A merge of `private/master` into a private
branch is then the same as a merge of `origin/master` into a public one
(h20 §3.3 item 4). `reviewed` gains the same keyword; `away.py`'s
`grant_lines` passes `home(root, head)` for it, and labels a `merge` act
`merged <h7> at <time>`.

For the private merge only, after the review passes, `checks(view, head) ->
str | None` (pure; the reason, or None when green) on the JSON of
`gh pr view <branch> -R expertanalytics/rasputin-private --json
state,baseRefName,headRefOid,statusCheckRollup`, run with a 20 s timeout
(the call took 0.74 s; one PR, one page). Refused, each with its own words:

1. gh exits non-zero, times out, or prints JSON without those keys: `the pull request's checks could not be read`.
2. Three conditions, each with its own reason: `state` is not `OPEN` (`the pull request is <state>, not open`); `baseRefName` is not `master` (`the pull request is into <base>, not master`); `headRefOid` is not the head (`the pull request's head is <h7>, not <h7>`).
3. `statusCheckRollup` is empty: `no checks reported`.
4. Any `CheckRun` whose `status` is not `COMPLETED`, or whose `conclusion` is not `SUCCESS`, `SKIPPED` or `NEUTRAL`; any `StatusContext` whose `state` is not `SUCCESS`. The reason names the first such check.
5. No item is a success (all skipped or neutral).

`--match-head-commit` still pins the head at the forge, as in h20 §1.

The hook's time. `.claude/settings.json` gives `guard_push.py` no `timeout`
(read 2026-10-17), so Claude Code's default for a command hook applies: 600 s
(code.claude.com/docs/en/hooks, *Common fields*, read 2026-10-17). The 20 s
`gh` call, `home`'s `git config --get` and two `rev-list` calls (10 s each,
new here) and the guard's own two `git --list-cmds` calls (10 s each, `guard_push.py`
line 138) add to at most 70 s, well inside it; `push_grant.py`'s other git calls have no timeout today
and are local reads. The same page says a timed-out `PreToolUse`
hook does not block the call, which then goes through the normal permission
flow; so the 20 s cap is what keeps the guard from ever reaching that limit,
and no `timeout` is added to the settings.

## 5. Fix 2: `.claude/hooks/gates_after_commit.py`

### 5.1 A tree with unmerged paths

In `main`, for each resolved tree (a `Path`), before `run_gates`:
`git ls-files -u` non-empty means a merge, cherry-pick, revert, rebase or
pull stopped on conflicts there. Then no gate runs in that tree, and its block
is

> A git operation stopped on conflicts in `<tree>`: `<n>` paths are unmerged (`<first three, comma-separated>`). No gates ran there. Any commit made before it stopped (an earlier pick of a rebase or cherry-pick, or an earlier command on the same line) is not checked either; the gates run on the next commit there.

The wording does not say no commit was made: a rebase or a range
cherry-pick that stops on its second commit has committed the first.

The opening "The commit is in, and below are ..." is printed only when at
least one block is a gate block or a `NOT CHECKED` block; a run whose only
blocks are unmerged trees prints them alone. Exit 2 either way, so the agent
reads it.

### 5.2 A gate whose script is absent

`gates(root)` is unchanged; `run_gates` skips a gate whose second word starts
with `tools/` and is not a file under `root`, and returns the line
`Skipped in <root>: <script> is not in this tree.` beside the failures. Ruff
gates are never skipped. When a tree's only lines are skips, the hook prints
them on stdout and exits 0: a skip is not red. When there are failures too,
the skip lines go in that tree's block. Today `lsm-melt-out/tools` lacks
`check_detria_boundary.py` and has the other two scripts.

A gate script the commit itself deleted is skipped the same way, exit 0.
That is accepted, for two reasons. Every gate script is under
`tools/check_`, a governed prefix (`GOVERNED_PREFIXES` in
`guard_governance.py`), so deleting one already needed Ola's prompt. And CI
runs each gate by name from `.github/workflows/main.yaml` (`python
tools/check_detria_boundary.py` and the others), so on a public branch the
missing script fails CI. The skip line says which script was missing, so the
deletion is visible at commit time too.

## 6. Fix 3: `tools/scratch_copy.py`

### 6.1 What breaks (reproduced 2026-10-17 on `lsm-lake-inflow` at `cea32bad`)

The private copy of the tool (`tools/scratch_copy.py` in that worktree, the
same code as public but for one docstring line) made a copy, and from it:

- `rasputin_lsm.__file__` was `.../worktrees/lsm-lake-inflow/src/rasputin_lsm/__init__.py`:
  the worktree's, not the copy's. The sitecustomize drops the editable
  finder, but the skbc `.pth` file's second line already put the worktree's
  `src` on `sys.path`, and the printed `PYTHONPATH` names `src_python`,
  which the private layout lacks. A mutant run there tests the worktree.
- `import rasputin_lsm._lsm` failed: the extension is in
  `site-packages/rasputin_lsm/` and was not copied.
- The tool copied `site-packages/tin_engine/_core*.so` into a new
  `src_python/tin_engine/` in the copy (the core is a plain install there,
  pinned by `pyproject.toml`); harmless, since the installed package wins,
  but wrong.
- `pytest -x` stopped at `test_cli_inflow.py` after 421 passes.

### 6.2 The rule

1. **Source roots**: those of `src_python` and `src` that exist in the copy,
   in that order. The printed `PYTHONPATH` is the sitecustomize directory,
   then each source root of the copy.
2. **Extensions**: for each package directory `P` (holding `__init__.py`)
   directly under a source root of the copy, every `*.so` under
   `ROOT/.venv/lib/python3.*/site-packages/P/` is copied to the same relative
   place under the copy's `P`. A site-packages package with no source in the
   copy is left alone. The warning when none is found names the packages
   looked for.
3. **The sitecustomize** drops the editable finders as today, and also
   every `sys.path` entry whose real path is `ROOT/src_python` or `ROOT/src`
   (paths written into the generated file). It drops nothing else: a
   prototype that dropped every entry under `ROOT` removed the venv's
   site-packages and lost `pytest`.
4. **The C++-differs warning** compares `git diff --quiet <rev> HEAD -- '*.cpp'
   '*.hpp' '*.h' '*CMakeLists.txt' ':(exclude)tests/'` instead of the
   directory list, since the private `src/` holds Python only.
5. `copy` reads the module's `ROOT` at call time (it does today), so the
   live check can run the public tool on a private worktree.

The prototype of rules 1 to 3, by hand on `lsm-lake-inflow` at `cea32bad`
into the session scratchpad: `rasputin_lsm` and `_lsm` both from the copy,
`tin_engine` from the venv; the full suite **3374 passed, 27 skipped, 0
failed in 434 s**.

## 7. Tests `@tester` writes first

All `harness`-marked, real git in `tmp_path`, no network.

**Fixture A (push grant)**: a repository with `refs/remotes/origin/master` at
M1 and `refs/remotes/private/master` at P1, both children of M0 (N = 1
private-only commit); branch `pub` = M1 + one code commit; branch `prv` =
P1 + one code commit; `remote.private.url` set to a path in `tmp_path` that
is never fetched (so N = 1 is read from the ref, with no network); a flag
granting both; a `reviews.jsonl` approving each head. `checks` is fed dicts,
never gh.

`tests/python/test_push_grant.py`:

1. Each private form parses to its `Act` with `remote="private"`, including
   `-C /abs`, `--repo`, `--repo=`, and the merge's options in two orders.
2. Each refusal of 4.2 items 1 to 3, one parametrised case each (about 18):
   `granted_act` returns None.
3. The public merge with `-R expertanalytics/rasputin` returns None, and so
   does `gh pr create -R expertanalytics/rasputin --head pub`.
4. `home`, one case per row of the 4.3 table that the fixture can reach:
   (a) `prv` is `"private"` (M = 0 < N = 1) and `pub` is `"origin"` (M = 1);
   (b) URL set, ref deleted (`git branch -dr private/master`): both get the
   missing-ref reason; (c) URL set, ref moved to M0 (`git update-ref`, so
   N = 0): both get the N = 0 reason; (d) URL unset and ref deleted: both are
   `"origin"`; (e) URL set, ref present, and `GIT_CONFIG_GLOBAL` (set with
   `monkeypatch.setenv`, which `home`'s git inherits) naming a `tmp_path` file
   holding `[broken`: both get the config reason with `exited 128`, the
   measured exit (4.3), and no `rev-list` is reached.
5. `judge`: (a) private push of `prv` passes and logs `"remote": "private"`;
   (b) origin push of `prv` is refused with the leak sentence; (c) private
   push of `pub` is refused with the second sentence; (d) the leak case of
   design review round 1: origin push of `prv` with the ref deleted, URL set,
   is refused with the missing-ref reason; (e) the same with the ref at M0
   (N = 0) is refused with the N = 0 reason; (f) a private push of `prv` in a
   clone with no URL and no ref is refused; (g) origin push of `prv` with the
   malformed `GIT_CONFIG_GLOBAL` of 4(e) is refused with the config reason.
6. The walk: on `prv`, a code commit after the review is refused; a clean
   merge of `private/master` after the review passes; a clean merge of
   `origin/master` into `prv` is refused (not its master).
7. `checks`, one parametrised case per condition, each asserting its own
   reason: green passes; `state` `MERGED`; `baseRefName` `main`;
   `headRefOid` another sha (the three of 4.4 item 2, each alone on an
   otherwise green view); an empty rollup (item 3); a `CheckRun` `IN_PROGRESS`,
   a `CheckRun` `FAILURE`, a `StatusContext` `PENDING` (item 4); all
   `SKIPPED` (item 5); and item 1 three ways: the stubbed runner exits 1,
   times out, prints `[]`.
8. A second private merge of the same head in one window is refused.

`tests/python/test_away.py`: `--back` on a granted `prv` reports the review
against `private/master` (`reviewed: yes`), and a logged `merge` act reads
`merged <h7>`.

`tests/python/test_gates_after_commit.py` (its existing `Repos` fixture,
whose gate stubs append `<gate> <cwd>` to a log when run):

9. A worktree without `tools/check_detria_boundary.py`, other gates green: a
   commit exits 0, stdout is one `Skipped in <wt>: tools/check_detria_boundary.py
   is not in this tree.` line, stderr empty.
10. The same with a red ruff stub: exit 2, the ruff block and the skip line,
    no detria entry.
11. A conflicting `git merge` in a worktree: exit 2, `stopped on
    conflicts`, nothing appended to the stubs' log, no `The commit is in`.
12. A cherry-pick of a two-commit range whose second commit conflicts (the
    first is committed): as 11.
13. `git commit -qm a && git merge <conflicting>`: as 11, with the sentence
    on an earlier commit.
14. A clean merge still runs the gates (regression).

`tests/python/test_scratch_copy.py` (**Fixture B**: a repository with
`src/pkg/__init__.py`, a fake venv whose site-packages holds
`pkg/_ext.cpython-3x-darwin.so` (stand-in bytes), `other/_x.so` with no
source, and a `.pth` naming `ROOT/src`):

15. `_ext` lands at `copy/src/pkg/`; nothing is made for `other`; no
    `src_python` appears.
16. The printed `PYTHONPATH` holds `copy/src`; a child Python run with it
    and the fake `.pth`'s entry on `sys.path` imports `pkg` from the copy.
17. The sitecustomize leaves `site-packages` on `sys.path` (the prototype's
    failure, pinned).
18. A change to a `.py` under `src/` does not warn that the extension is
    from HEAD; a change to a `.cpp` does.
19. The existing public-layout tests pass unchanged, except the ones that
    pin the old warning text, which `@tester` updates in the red commit.

No suite here is invariant-critical; no mutation round.

## 8. Live checks (after the merge, on master; hooks run from the settings a session starts with)

1. `python3 -c` importing the merged `tools/scratch_copy.py` with `ROOT` set
   to an `lsm-*` worktree with a `.venv`, copy `HEAD`, run the printed
   command: `rasputin_lsm` and `_lsm` from the copy; the suite's count beside
   section 6.2's 3374 / 27.
2. A commit in an `lsm-*` worktree: the skip line, no red.
3. A conflicting `git merge` in a scratchpad repository: the unmerged block,
   no gates.
4. The private forms are judged in the next away window Ola grants a private
   branch in; the `--back` report shows the acts.

## 9. Delivery to the private repository

The hooks run from the main checkout for every worktree under it, so fixes 1
and 2 act on `lsm-*` worktrees on merge. Fix 3 does not: an agent in a
private worktree runs that tree's own `tools/scratch_copy.py`. It reaches the
private repository by the harness sync of `51-open-core.md` §5.3 (a private
PR that checks the harness paths out of a named public commit), which
`@orchestrator` owns.

## 10. Rule text (lands with Ola's prompt)

`.claude/REQUIRED-READING.md`, *The harness*, the sentence on grants, from
"may be pushed (as `git push origin <sha>:refs/heads/<branch>`), opened as a
pull request and enqueued (once per head)" to:

> may be pushed (as `git push origin <sha>:refs/heads/<branch>`, or `git push
> private ...` for a branch holding commits of `private/master`), opened as a
> pull request and enqueued, or on the private repository merged once its
> checks are green (once per head),

## 11. Size

Net production lines (`CLAUDE.md` §2), estimates:

| File | Fix | Estimate |
|---|---|---|
| `tools/push_grant.py` (forms, `home` with its fail-closed rows, `checks` with a reason per condition, the remote through the walk) | 1 | about 100 |
| `tools/away.py` | 1 | about 5 |
| `.claude/hooks/gates_after_commit.py` | 2 | about 25 |
| `tools/scratch_copy.py` | 3 | about 20 |
| Total | | about 150 |

One PR, full path (over 50 lines): design review, red, green, code review,
Ola's prompts for the governed files and the rule file.

## Review

- **Design review, round 1, 2026-10-17.** Commit `949a6b2b`, 0 counted lines, @reviewer: CHANGES REQUESTED. The leak check fails open: with `refs/remotes/private/master` deleted (no prompt needed) or counting no private-only commits, every head counts as public, so a granted `git push origin` of an `lsm-*` branch would publish private history; refuse the origin forms then, with a test each. Suggestions: refuse `gh pr create -R expertanalytics/rasputin` in test 3; test 7 per condition; a rebase or cherry-pick stopped partway did make commits; a gate deleted by the commit; the guard hook timeout against the 20 s `gh` call.
- **Design review, round 2, 2026-10-17.** Commit `892c4de8`, 0 counted lines, @reviewer: CHANGES REQUESTED. The six rows of 4.3 fail closed and the `remote.private.url` anchor needs a prompt to remove; one hole: "set up" is read as `git config --get` exit 0, so a read error (exit 128 from a malformed `~/.gitconfig`) counts as not set and fails open. Exit 1 means not set; any other exit or a timeout refuses; that call gets a timeout inside the hook budget.
- **Design review, round 3, 2026-10-17.** Commit `b32c527b`, 0 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session: the hook-time sum is 70 s (two `git --list-cmds` calls, `guard_push.py` line 138). Idea taken for the red step: a monkeypatched `TimeoutExpired` on the config call tests the timeout reason.
