# Harness h20: away-mode push grants, and wording fixes without a second review

Status: **design**, 2026-10-13, by @architect, on `worktree-night-push` off
master `97c10c4d`. One PR, by day: every production file it changes is
governed. Implements two proposals Ola approved on 2026-10-13 from the morning
retrospective `docs/retrospectives/2026-10-13-night-40-41-42.md` (branch
`worktree-retro-1013`, #244), section 7: P2, the push grant (first approved as
P2 of `docs/retrospectives/2026-10-09-night-32-33.md`, lines 186-192 and
274-280), and P3, a reviewer's named wording fix needs no second review.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature.* Two code-review systems settle when an approval survives new
commits, and this design borrows from both.

- Gerrit's label `copyCondition` carries a vote to a new patch set by change
  kind (gerrit-review.googlesource.com/Documentation/config-labels.html, read
  2026-10-13): `NO_CODE_CHANGE` ("the same parent tree ... and the same code
  diff"), `TRIVIAL_REBASE` ("the same diff (including context lines)") and
  `MERGE_FIRST_PARENT_UPDATE` ("a new merge commit which only differs from
  the merge commit in the previous patch set in its first parent").
- GitHub's branch protection can "dismiss stale pull request approvals when
  commits are pushed that affect the diff", and can require that "the most
  recent reviewable push must be approved" (docs.github.com, *About protected
  branches*, read 2026-10-13).

What differs, and why. This design keeps an approval across (a) commits that
change only prose files and (b) merges of master whose result differs from
git's automatic merge only in prose files. (b) is close to Gerrit's
`MERGE_FIRST_PARENT_UPDATE` but also accepts a hand-resolved conflict in a
prose file, because the conflicts that idle the night are in `ROADMAP.md`
(#239, retrospective 2026-10-13 §1 item 2). (a) is weaker than Gerrit's
`NO_CODE_CHANGE`: a post-review prose commit is not checked against what the
reviewer named. That is the price of a mechanical check; §4.4 names it.

The forge enforces the last step: `gh pr merge --match-head-commit <sha>` sets
the GraphQL `expectedHeadOid`, so GitHub refuses the enqueue if the pull
request's head is not the commit the guard judged. gh 2.101's source passes it
on every merge path, the merge-queue path included (cli/cli
`pkg/cmd/pr/merge/http.go`, `mergePullRequest`, read 2026-10-13).

*Legacy.* Nothing; the legacy tree has no harness.

```
$ git grep -l -i -e 'unattended' -e 'guard_push' -e 'away.py' legacy-archive -- legacy
(no output, exit 1)
```

*Ideas file* (`../rasputin_scratch/ideas/ideas.md`, read 2026-10-13). No entry
is taken up. Overlapping, left out: "Cap on runner and shell words per line in
the push guard" (same file, other rule); "Governed-file asks answered without
a prompt in auto mode" (the grant does not touch asks while Ola is present).
The same retrospective's P5 (a GraphQL read passes the push guard) and P7
(`git merge-tree` before a second enqueue) touch the same files and are not
part of this increment.

## 2. Scope

| # | What | Where |
|---|---|---|
| 1 | `away.py <duration> --grant <branch>` (repeatable) writes the grants into the flag | `tools/away.py` |
| 2 | The push guard passes a granted act on a reviewed head while away | `tools/push_grant.py` (new), `.claude/hooks/guard_push.py` |
| 3 | The review verdict travels as a `Review:` commit trailer | `docs/increments/README.md` step 4, `.claude/agents/reviewer.md` |
| 4 | A reviewer's named wording fix needs no second round (P3) | `docs/increments/README.md` step 4, `.claude/REQUIRED-READING.md` |
| 5 | `away.py --back` reports each grant | `tools/away.py` |
| 6 | The two new imports of a hook are governed | `.claude/hooks/guard_governance.py` |

Out: grants while Ola is present (the guard asks, as today); any act other
than the three forms of §3.2; adding a grant during a window without
re-entering it.

## 3. Design

### 3.1 The grant (away.py)

```
python3 tools/away.py 8h --grant worktree-threads --grant worktree-cite-gate
```

- `--grant` is repeatable and only valid with a duration (with `--back`, or
  alone: usage error, exit 2). A grant is a local branch name, checked with
  `git check-ref-format --branch <name>`; `master`, `main` and `HEAD` are
  refused (exit 2, before the terminal opens). A branch that does not exist
  yet is accepted, and the confirmation says "does not exist yet" beside it.
- The confirmation on `/dev/tty` lists the grants: "Granted while away: push,
  pull request and one enqueue per reviewed head of: a, b." The same `y`
  confirms the window and the grants.
- The flag gains one key, `"grants": ["a", "b"]`, written `[]` without
  `--grant`; `windows.jsonl`'s enter entry gains the same key.
- **Expiry is the window's**: the grants live in the flag, so they end when
  the window ends (`until`) or on `--back`, which deletes the flag.
  Re-entering replaces the flag, so the new call's `--grant` list is the whole
  list. A grant survives a context loss: it is Ola's own write, read by the
  guard, not a yes remembered by a session.
- `harness_mode.read_mode` is unchanged: an unknown key does not break the
  flag. A `grants` value that is not a list of strings counts as no grants
  (fail closed), and does not make the flag broken.

### 3.2 The granted acts (push_grant.py)

The guard passes, while the mode is `on`, exactly these three forms, for a
branch `<b>` in the flag's grants. The command must be **one simple command
and nothing else**: `shell_scan.parse(command)` returns one `Simple`, and
`shlex.split(command) == simple.argv` (no wrapper, assignment, redirection,
operator or substitution). Any other word or option makes it not a granted
act, and it is refused as today.

| Act | Form (options in any order) |
|---|---|
| push | `git [-C <dir>] push [-u \| --set-upstream] origin <b>` |
| pull request | `gh pr create --head <b> [--base master] [--title <t>] [--body-file <f>] [--draft]` (`-H`, `-B`, `-t`, `-F`, `-d` and `--opt=value` spellings too) |
| enqueue | `gh pr merge <b> --match-head-commit <40-hex sha> [--merge]` |

- A granted PR body comes from `--body-file`, because `--body "$(cat <<EOF …)"`
  is a substitution and fails the one-command test.
- `--admin`, `-d`/`--delete-branch`, `--auto`, `--squash`, `--rebase`, `-R`,
  `--force`, `--no-verify`, a refspec, a remote other than `origin`: none is
  in the forms, so each is refused.
- The repository: `git rev-parse --git-common-dir` in `<dir>` (push with
  `-C`) or the event's `cwd` must be the guard's own (the parent of
  `harness_mode.state_dir(ROOT)`). Otherwise not granted.
- The head `H` is `git rev-parse refs/heads/<b>` in that repository. For the
  enqueue, the `--match-head-commit` sha must equal `H`.
- Each act also needs `reviewed(H)` (§3.3) to hold.
- **One enqueue per reviewed head**: an enqueue is refused if the grants log
  (§3.5) already records an enqueue of `<b>` at `H` in this window. A new
  head (a master merge after the queue ejected the PR) may be enqueued once.
  Pushes and PR creates are not counted.

The guard checks the push's head before git runs. If another writer commits
on `<b>` between the check and the push, an unreviewed commit reaches the
remote branch; the enqueue still cannot take it, since its
`--match-head-commit` must be a reviewed `H` and GitHub refuses a mismatched
head. Accepted as a limit.

### 3.3 "Reviewed": the rule the guard runs

The verdict travels as a git trailer on the main session's recording commit
(the commit that copies a review round into the increment file):

```
Review: APPROVED <40-hex sha>
Review: APPROVED WITH WORDINGS <40-hex sha>
Review: CHANGES REQUESTED <40-hex sha>
```

`<sha>` is the head of the range the reviewer reviewed. Every recording
commit carries exactly one `Review:` trailer. A trailer is chosen over parsing
the `## Review` prose: the records' wording varies (compare the `## Review`
sections of `39-se2-catchments.md` and `h19-commit-gate-tree.md`), and a
trailer also works on a branch with no increment file, such as a
retrospective's.

`reviewed(H)` holds when all of these hold (git calls only; no network):

1. Walk `git log --first-parent refs/remotes/origin/master..H`, newest first;
   `V` is the first commit with a `Review:` trailer
   (`%(trailers:key=Review,valueonly)`). None: "no review recorded". The
   first-parent walk keeps master's own recording commits, brought in by a
   master merge, out of the search.
2. `V`'s verdict is `APPROVED` or `APPROVED WITH WORDINGS`. Otherwise: "last
   review round: CHANGES REQUESTED".
3. `R`, the trailer's sha, is a full 40-hex sha on the first-parent chain from
   `H` (`R` is `H` or a first-parent ancestor of it, and reachable in the
   walk of step 1).
4. Every commit `C` on the first-parent chain after `R` up to `H` (so `V`
   too) is one of:
   - **a prose commit**: one parent, and every path of
     `git diff --name-only --no-renames C^ C` is prose by
     `tools/ci_changes.py`'s `is_prose` (Markdown at the root or under
     `docs/`, less `NOT_PROSE`);
   - **a master merge**: two parents, the second an ancestor of
     `refs/remotes/origin/master`, and every path where `C`'s tree differs
     from `git merge-tree --write-tree C^1 C^2`'s tree (the first line of its
     output; exit 0 clean or 1 conflicted, anything else fails) is prose.

   Anything else (a code commit, a merge with a code conflict resolved by
   hand, an octopus merge): "commit <short> changes <path> after the review of
   <R short>".

Checked by hand on 2026-10-13 in a scratch repository (git 2.55): a
conflicted `merge-tree --write-tree` prints the tree first and exits 1, and
the diff between that tree and the hand-resolved merge names only the
conflicted file (`R.md`), not the file master changed cleanly (`c.py`).

### 3.4 What the rule trusts

The trailer is written by the main session, not by `@reviewer` (who is
read-only), and nothing authenticates it. The guard is a tripwire against
pushing an unreviewed or since-changed head, not proof of review; Ola's trust
in a granted night rests on the grant naming the branch, the trailer
discipline, and the `--back` report (§3.6), which names every granted act and
the `V` commit it relied on. A post-review prose commit is accepted whether or
not the reviewer named it (§1). Question 1 (§7) offers a stronger route.

### 3.5 Guard wiring (guard_push.py, push_grant.py)

`tools/push_grant.py`, new, standard library only, imports `harness_mode`,
`shell_scan` and `ci_changes`. One public entry point:

```python
def judge(event: dict[str, Any], command: str, reasons: list[str]) -> Grant | None
```

`None`: not a grant-shaped act for a granted branch while away; the guard goes
on as today. `Grant(passed: bool, branch: str, act: str, head: str | None,
why: str)`: a granted branch's act, passed or not, and the reason in plain
words. Internals, each testable alone: `granted_act(argv) -> Act | None`
(pure, §3.2's forms), `grants(flag) -> tuple[str, ...]` (pure),
`reviewed(repo, head) -> str | None` (None or the reason, §3.3), `log_use`.

`judge` returns `None` unless the mode is `on` and `reasons` is `[PUSH]` or
`[PR]` (guard_push's reason strings: no force, no `--no-verify`, no history
rewrite, no unknown command).

In `guard_push.main`, after `reasons` is computed and before `emit`:

- `Grant(passed=True)`: append one line to `<state>/grants-used.jsonl`
  (`{"at", "window": flag since, "branch", "act", "head", "review": V}`, one
  write, as `harness_mode.append_queue`), then print nothing and return: the
  act passes. If the line cannot be written, the act is refused (an enqueue
  that cannot be counted is not granted).
- `Grant(passed=False)`: `why` gains "; granted branch <b>, but <reason>", so
  the refusal and its queue line say why the grant did not apply.
- Any exception in `push_grant` falls to the existing `except`: deny.

### 3.6 `away.py --back`

After the queue and before the `ASK OLA` lines:

```
Push grants (2):
  worktree-threads: pushed 1a2b3c4 at 23:41 UTC; PR opened; enqueued 1a2b3c4 at 23:44 UTC (review in 5d6e7f8)
  worktree-cite-gate: not used; head 9a8b7c6 reviewed: no (commit 9a8b7c6 changes tools/x.py after the review of 69a4eb8c)
```

Uses this window's lines of `grants-used.jsonl` (`window` equal to the flag's
`since`); for an unused grant, `reviewed()` at the branch's current head. No
grants: `Push grants: (none)`. The log is not archived: its `window` key
separates windows.

### 3.7 Rule text (lands in the same PR, governed, with Ola's prompt)

`.claude/REQUIRED-READING.md@97c10c4d:179`, in the paragraph that opens
"**Unattended mode.**", after "is never redone by any route.":

> While away, a branch Ola named with `away.py --grant <branch>` may be
> pushed, opened as a pull request and enqueued (once per head) without a
> fresh yes, provided its last review round, recorded as a `Review:` trailer
> (`docs/increments/README.md`, step 4), approved its head; the guard checks
> this, and the grant ends with the window.

*A fresh yes, every time* (line 121), after "one yes covers one enqueue":
", except a branch Ola granted for an away window (*Unattended mode*)".

*`@reviewer` runs … again before a push that follows a round of findings*
(line 134), after "round of findings": " (a verdict other than `APPROVED`,
`docs/increments/README.md` step 4)".

`docs/increments/README.md`, step 4, after "Check: `grep -n '^## Review'
docs/increments/NN-name.md` on the branch.":

> The recording commit carries one trailer, `Review: <verdict> <sha>`: the
> verdict as `@reviewer` gave it and the full sha of the head it reviewed.
> When `@reviewer`'s only findings are wordings in prose files (Markdown at
> the root or under `docs/`) that it names line by line with the new text,
> its verdict is `APPROVED WITH WORDINGS`. The main session applies them
> exactly, in the recording commit, and a branch whose only change after the
> review is that commit counts as reviewed: no second round. A wording in a
> code file, or any other finding, is `CHANGES REQUESTED`.

`.claude/agents/reviewer.md`, check list item 1 (line 55), replacing
"`APPROVED` or `CHANGES REQUESTED` (with explicit blocking issues)." by:

> `APPROVED`, `APPROVED WITH WORDINGS` (only prose wordings, each named line
> by line with its new text; `docs/increments/README.md` step 4), or
> `CHANGES REQUESTED` (with explicit blocking issues), and the full sha of
> the head you reviewed.

The frontmatter `description` keeps "returning APPROVED or CHANGES
REQUESTED"; unchanged, it is still true of the two outcomes.

## 4. Tests @tester writes first

All in the harness suites (`-m harness`), on `harness_fixtures.make_repo`
copies; `COPIED` gains `tools/push_grant.py` and `tools/ci_changes.py`. The
fixture repository for §4.2-4.4 (measured by construction, not sampled):
`master` with 1 commit; `refs/remotes/origin/master` set to it with `git
update-ref` (no remote needed); branch `feat` with a code commit
`tools/x.py` (1 line), and a recording commit changing
`docs/increments/h99.md` (1 line) with `Review: APPROVED <code commit sha>`.
Flag on, `grants: ["feat"]`, unless a test says otherwise.

### 4.1 `away.py` (test_away.py)

1. `8h --grant feat --grant other` writes `"grants": ["feat", "other"]` to
   the flag and to the `windows.jsonl` enter entry; the tty text lists both.
2. `8h` alone writes `"grants": []`; a flag with no `grants` key (written
   before h20) reads as no grants.
3. `--grant master`, `--grant main`, `--grant HEAD`, `--grant 'a..b'`: exit 2,
   tty never opened, no flag.
4. `--back --grant feat` and `--grant feat` alone: exit 2.
5. A grant for a branch that does not exist: accepted, tty text says "does
   not exist yet".
6. Re-entering with `--grant b` only: the flag's grants are `["b"]`.
7. `--back` with a 3-line `grants-used.jsonl` (2 lines of this window, 1 of
   an earlier one): prints this window's 2 acts under `Push grants`, not the
   third; an unused grant prints its head and `reviewed: yes` or the reason;
   no grants prints `Push grants: (none)`.

### 4.2 The forms (`granted_act`, pure, parametrised)

8. Each form of §3.2 parses, including `-C <dir>`, `-u`, `--set-upstream`,
   `-H`, `--head=feat`, `-B master`, `--body-file f`, `--merge`.
9. Each of these is `None`: `--force`, `-f`, `--force-with-lease`,
   `--no-verify`, `origin feat:main`, `origin HEAD`, `upstream feat`, a second
   branch, `--admin`, `-d`, `--auto`, `--squash`, `--rebase`, `-R x/y`,
   `--body text`, `--base other`, a short sha (39 hex), a missing
   `--match-head-commit`, `gh pr merge 241` (a number, not the branch).

### 4.3 The guard end to end (test_guard_push.py, through the hook copy)

10. Flag on, `feat` granted and reviewed: `git push -u origin feat`,
    `gh pr create --head feat --title t --body-file f` and `gh pr merge feat
    --match-head-commit <H> --merge` each print no decision; three lines in
    `grants-used.jsonl`; nothing in `queue.jsonl`.
11. The same three with the mode off: ask, as today; nothing logged.
12. With the flag broken: deny, as today.
13. Not granted (`grants: ["other"]`): deny and queue, `why` without the
    grant clause.
14. Each wrapper is refused: `env A=1 git push origin feat`,
    `git push origin feat 2>&1`, `git push origin feat && echo`,
    `sh -c 'git push origin feat'`, `timeout 5 git push origin feat`.
15. A second `gh pr merge feat --match-head-commit <H>` in the window: deny,
    `why` says the head was enqueued already. After a clean master merge
    (new `H2`): one enqueue at `H2` passes.
16. `--match-head-commit` of a sha other than `H`: deny.
17. `git -C <a directory of another repository> push origin feat`: deny.
18. `grants-used.jsonl` unwritable (a directory at its path): deny.

### 4.4 `reviewed` (the invariant-critical suite of this increment)

Each case builds on the fixture; the expected result is pass or the named
reason.

19. As built: pass.
20. No trailer on the branch: "no review recorded".
21. Last trailer `CHANGES REQUESTED`, an earlier one `APPROVED`: refused.
22. `APPROVED WITH WORDINGS` and a prose-only recording commit: pass.
23. A code commit after the recording commit: refused, naming the commit
    and `tools/x.py`.
24. A recording commit that also changes `tools/x.py`: refused.
25. A post-review commit changing `CLAUDE.md`, `README.md` or
    `docs/increments/README.md` (`NOT_PROSE`): refused.
26. A clean master merge (master changed `src/y.py`): pass.
27. A master merge with a `ROADMAP.md` conflict resolved by hand: pass.
28. A master merge with a `tools/x.py` conflict resolved by hand: refused.
29. A clean merge whose second parent is not on `origin/master` (another
    branch): refused.
30. A trailer naming a sha not on the first-parent chain, or a short sha:
    refused.
31. A master recording commit with `Review: APPROVED <sha>` brought in by a
    master merge, and no trailer on the branch's own commits: "no review
    recorded" (the first-parent walk).

**Mutation targets** (the kill record must cover each): the first-parent
flag dropped (kills 31); `is_prose` replaced by "any `.md`" (25); the
second-parent ancestry check dropped (29); the merge-tree diff replaced by
"the merge is clean" (27 or 28, one of them); the verdict test accepting any
verdict (21); the one-command test dropped (14); the per-head enqueue count
dropped (15).

## 5. Size

Against the 700-line limit (`CLAUDE.md` §2), net production lines:

| File | Estimate |
|---|---|
| `tools/push_grant.py` (new) | about 140 |
| `tools/away.py` (`--grant`, flag key, tty text, `--back` section) | about 40 |
| `.claude/hooks/guard_push.py` (call `judge`, log, extend `why`) | about 15 |
| `.claude/hooks/guard_governance.py` (`GOVERNED` gains `tools/push_grant.py` and `tools/ci_changes.py`: the push guard imports them, so an edit to either is an edit to the guard) | 2 |
| Total | about 200 |

`pyproject.toml`'s mypy `files` gains `tools/push_grant.py` (not counted).
Rule text, §3.7: about +150 words over four files.

## 6. Order

Design review; red (tests of §4, all failing); green; code review; the rule
text of §3.7 applied by the main session at Ola's prompts, in the same PR; the
recording commit with its own `Review:` trailer. The guard's edit and the new
`GOVERNED` entries are governed writes, so the round runs by day.

## 7. Questions for Ola

1. **Who vouches for the review verdict the guard reads?** The design has
   the main session write the reviewer's verdict into the recording commit
   as a line `Review: APPROVED <commit>`, and the guard trusts that line. A
   session that wrote it falsely would get an unreviewed branch pushed and
   enqueued while you are away. The stronger route is a new hook that copies
   the verdict straight from `@reviewer`'s handback into the harness state,
   where no agent can write: about 40 more lines and a change to
   `.claude/settings.json`. Default: the commit line, as designed; the
   morning report lists every granted act and the commit it relied on.
2. **How many enqueues does a grant cover?** Today one yes covers one
   enqueue. The design allows one enqueue per reviewed head commit, so if the
   merge queue throws the pull request out (as it did #239 on the night of
   2026-10-12), a merge of master, a push and one more enqueue can happen
   before morning; the same head is never enqueued twice. Default: one per
   head, as designed. The alternative is one enqueue per grant per night.
