# Harness h27: an amend of an unpushed commit is not a history rewrite

Status: **designed**, 2026-10-20, on `worktree-h27-unpushed-amend` off master
`d3e55b0a`. Short path (`docs/increments/README.md`): about 25 counted lines,
no C++ file, no refine or mesh code. The green step edits
`.claude/hooks/guard_push.py`, a governed file, so it runs with Ola at the
prompt.

Source: Ola said yes on 2026-10-19 to "an amend of a commit that was never
pushed is not a history rewrite, so the guard should not refuse it". That
night `git -C <worktree> commit --amend --no-edit -a`, on a commit never
pushed, was refused in away mode as "rewrites history".

**The picture.** The push guard asks (and, while Ola is away, refuses) every
`git commit --amend`, because amending a published commit forces a later
force-push. An amend only replaces the commit HEAD points at. If no
remote-tracking branch and no tag holds that commit, nobody else can have it,
and the amend is ordinary local work. The guard checks that before asking.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature* (git 2.55.0, `git --version`; read 2026-10-20).

- `git help push`: a push to a branch that the remote's fetch refspec maps
  (`+refs/heads/*:refs/remotes/origin/*` here, and the same for `private`,
  `git config --get-all remote.<name>.fetch`) also updates the matching
  remote-tracking ref locally. So a commit pushed from this clone, in any
  worktree (they share refs), is held by some `refs/remotes/*` ref. The
  test in section 3 pins this rather than trusting the page.
- `git for-each-ref --contains <commit> refs/remotes refs/tags` lists the
  refs whose tip descends from the commit. Run in this worktree on HEAD
  (`d3e55b0a`, on master): it printed `refs/remotes/origin/HEAD` and
  `refs/remotes/origin/master`, exit 0. On an unknown commit it exits 129
  ("no such commit").
- git accepts any unambiguous prefix of a long option (`gitcli`). Probed
  read-only with `--dry-run` in this worktree: `git commit --am --dry-run
  --allow-empty` and `--amen` exit 0 as amends; `--amendx` is rejected.
  Today's guard looks only for the full word `--amend`, so
  `git commit --amen` passes it silently. Section 2 closes that too, since
  the exemption must be keyed to the same detection.

*Legacy.* Nothing: the legacy tree is the old meshing product, with no
harness.

```
$ git grep -l -i -e 'amend' -e 'remote-tracking' legacy-archive -- legacy
(no output, exit 1)
```

*The ideas file* (`../rasputin_scratch/ideas/ideas.md`, read 2026-10-20):
no entry on amends or history rewrites; nothing taken up or overlapped.

## 2. The change

All in `.claude/hooks/guard_push.py`.

**(a) Amend detection.** In `publishes`, `"--amend" in args` becomes "some
argument is one of `--am`, `--ame`, `--amen`, `--amend`" (a module constant
`AMENDS`). `--a` is ambiguous to git and fails on its own.

**(b) The exemption.** A new function, called from `main` right after the
reasons are collected and before the grant check:

```python
def unpushed_amend(command: str, cwd: str | None) -> bool:
    """True when the line is one bare amend and its HEAD is held by no remote-tracking ref or tag."""
```

`main` returns 0 (silent, in every mode) when `reasons == [HISTORY]` and
`unpushed_amend(command, event.get("cwd"))` is true. Any other reason on the
line (`--no-verify`, a push, a second rewrite) keeps the ask, exemption or
not.

`unpushed_amend` returns True only if every row holds, checked in order;
anything else, an exception included, returns False (the ask stands, as
today):

| Check | Why |
|---|---|
| `shell_scan` imported, `shell_scan.parse(command)` gives exactly one simple command, and `shlex.split(command)` equals its argv | One command, no `cd`, no `&&`, no assignment in front (`GIT_DIR=... git ...`), no substitution: HEAD cannot move between the check and the amend |
| `EXPANDS_NOTHING.fullmatch(command)` | No word the shell expands; messages go in single quotes or the line asks |
| argv is `git`, then either nothing or exactly `-C <absolute path>`, then `commit`, and an `AMENDS` word among the rest | No `-c`, `--git-dir`, `--work-tree` or second `-C` that could point git elsewhere; a relative `-C` asks |
| the directory is the `-C` path, else `cwd`; it is a non-empty string | No cwd in the event: unknown, so refuse |
| `git -C <dir> rev-parse --verify -q HEAD^{commit}` exits 0 with a sha | No HEAD (a new repository): refuse |
| `git -C <dir> for-each-ref --contains <sha> --count=1 --format=%(refname) refs/remotes refs/tags` exits 0 with empty output | Held by any remote-tracking branch of any remote, or by any tag (tags may have been pushed): refuse |

Each `git` call has `timeout=10`, `capture_output=True`, `check=False`; a
timeout or non-zero exit is False.

**What counts as pushed.** Any ref under `refs/remotes/` or `refs/tags/`
whose tip is the commit or descends from it. Fetch state can be stale; every
way it is stale errs toward refusing except one, listed under
*Left out*.

- A push from this clone updates the tracking ref (section 1), so the
  ordinary case is seen.
- The remote moved on (a merge queue, `gh pr update-branch`): the old tracking
  ref still holds HEAD, so the guard asks.
- A worktree path in `-C`: worktrees share refs, so the check sees pushes
  from every worktree. A `-C` into another repository judges that
  repository's own refs, which is what the amend there touches.

**Text-rule path.** When `shell_scan` cannot parse the line, the text rules
run as before and no exemption applies.

The module docstring gains one sentence: an amend whose HEAD is held by no
remote-tracking ref or tag, on a bare one-command line, passes; anything
unknown asks.

**Estimate**: about 25 counted lines (the function about 20, the constant 1,
the `publishes` line 1, the `main` lines 2), measured with
`python3 tools/count_loc.py`.

## 3. Tests

In `tests/python/test_guard_push.py`. Fixtures are built by the test with
`git` in `tmp_path`, as `make_repo` and `add_worktree` already do. Each
fixture is `make_repo` (one commit), plus a bare repository
(`git init --bare`) added as remote `origin`, with `git push origin master`
run by the test itself so the tracking ref comes from git and not from
`update-ref`.

Pass, silent with modes `off`, `on`, `expired`, and nothing queued:

1. One new commit after the push, then `git commit --amend --no-edit` with no
   `-C` (cwd = repo).
2. The same with `git -C <abs repo> commit --amend --no-edit -a` (tonight's
   form).
3. A worktree from `add_worktree` with one commit of its own, amended with
   `-C <abs worktree>`.
4. `--amen` and `--am` forms of 1.

Ask by day, refuse and queue while away, reason contains the HISTORY text:

5. The pushed HEAD itself amended (no new commit), with and without `-C`.
6. A new commit, then a tag on it (`git tag t`), then the amend.
7. A new commit pushed under another branch name (`git push origin
   HEAD:refs/heads/x`): held by `refs/remotes/origin/x`.
8. A new repository with no commit (`git init` only), amended: no HEAD.
9. Unpushed HEAD, but the line is not bare: `cd <repo> && git commit
   --amend`, `GIT_DIR=<repo>/.git git commit --amend`, `git -c a.b=c commit
   --amend`, `git -C rel commit --amend` (relative), `git commit --amend -m
   "x"` (double quotes), `git commit --amend --no-verify` (two reasons).
10. Unpushed HEAD, event with no `cwd` key, no `-C`.
11. `git commit --amen` and `--am` on a pushed HEAD: today these pass; after
    (a) they ask. This is the one row that fails today for a different reason
    than the exemption.

Existing rows: `ASKED["git commit --amend"]` runs in the `repo` fixture,
which has no remote, so after the change it would pass. `@tester` moves it to
a fixture whose HEAD is pushed (row 5 covers it) and says so in the red
commit.

Not mutation-tested (short path; no invariant-critical suite).

## 4. Left out

- **Rebase and `reset --hard`.** A rebase rewrites every commit between its
  base and HEAD, so "HEAD is unpushed" is not enough: an older commit in the
  range may be pushed. The check needs the rebase's base parsed from its
  arguments (`--onto`, `--root`, the upstream); not the same check, so it
  stays asked. `reset --hard` moves no published commit but throws away work
  in the tree, which is a different question.
- **The one stale case that passes.** A commit pushed from another clone, or
  only to a pull-request ref, and whose tracking ref was then pruned
  (`fetch --prune` after the branch was deleted on GitHub without a merge),
  is held by no local ref. Amending it changes nothing anyone tracks; the
  closed pull request keeps the old commit. Accepted.
- **The hook's own environment.** The guard's `git` calls run with the hook's
  environment, not the Bash tool's shell. A `GIT_DIR` exported in Ola's shell
  profile would make the two disagree. Assumed absent; not checked.
