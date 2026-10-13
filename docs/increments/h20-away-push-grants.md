# Harness h20: away-mode push grants, and wording fixes without a second review

Status: **green, code review round 1 recorded**, 2026-10-13, on `worktree-night-push` off
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
reviewer named. That is the price of a mechanical check; §3.4 names it.

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
| 3 | A hook copies the verdict from `@reviewer`'s handback into the harness state | `.claude/hooks/record_review.py` (new), `.claude/settings.json` (§3.8), `docs/increments/README.md` step 4, `.claude/agents/reviewer.md` |
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
  alone: usage error, exit 2). A grant is a local branch name that matches
  `[A-Za-z0-9._/-]+` (so no shell expansion can change it: no `~`, `=`, `$`,
  glob or brace character) and passes `git check-ref-format --branch <name>`;
  `master`, `main` and `HEAD` are
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
| push | `git [-C <dir>] push origin <H>:refs/heads/<b>`, `<H>` the 40-hex head the guard judged |
| pull request | `gh pr create --head <b> [--base master] [--title <t>] [--body-file <f>] [--draft]` (`-H`, `-B`, `-t`, `-F`, `-d` and `--opt=value` spellings too) |
| enqueue | `gh pr merge <b> --match-head-commit <40-hex sha> [--merge]` |

- A granted PR body comes from `--body-file`, because `--body "$(cat <<EOF …)"`
  is a substitution and fails the one-command test.
- `--admin`, `-d`/`--delete-branch`, `--auto`, `--squash`, `--rebase`, `-R`,
  `--force`, `-u`, `--no-verify`, a `+` or any other refspec, a remote other
  than `origin`: none is in the forms, so each is refused.
- **No word the shell could expand.** Every word of the argv matches
  `[A-Za-z0-9._/:@-]+` (with `=` allowed only inside a `--opt=value` word)
  and none starts with `~` or `=`; in a `--opt=value` word the value does
  not start with `~` or `=` either, because zsh's `magic_equal_subst`
  expands after the `=` (`--body-file==ls` becomes `--body-file=/bin/ls`,
  @reviewer, design review round 2); the one exception is the `--title` value,
  which must stand in single quotes in the command text and contain no single
  quote. The one-command test alone lets these through, because
  `shlex.split` does not expand: `git -C ~/other push …` and `git -C =other
  push …` pass it and zsh expands both to another directory (probed by
  @reviewer, design review round 1); so does zsh's `$H:r` modifier, which
  turns `$H:refs/heads/b` into another word.
- The repository: `<dir>` of `-C` must be an absolute path (it starts with
  `/`, so it cannot start with `~` or `=`); `git rev-parse --git-common-dir`
  in `<dir>`, or in the event's `cwd` without `-C`, must be the guard's own
  (the parent of `harness_mode.state_dir(ROOT)`). Otherwise not granted.
- The head `H` is `git rev-parse refs/heads/<b>` in that repository. The
  push's source sha and the enqueue's `--match-head-commit` sha must equal
  `H`.
- Each act also needs `reviewed(H)` (§3.3) to hold.
- **One enqueue per reviewed head** (Ola's ruling 2): an enqueue is refused
  if the grants log (§3.5) already records an enqueue of `<b>` at `H` in this
  window. A new head (a master merge after the queue ejected the PR) may be
  enqueued once. Pushes and PR creates are not counted.

**Why the push names the sha.** On a repository with a merge queue, `gh pr
merge` turns on auto-merge when checks are still pending, and GitHub keeps
auto-merge on after a push by anyone with write access (@reviewer, design
review round 1). So any commit that reaches the remote branch after an
enqueue can merge. With `git push origin <b>`, a commit made on `<b>` between
the guard's check and git's read of the ref would be pushed unreviewed. The
form `git push origin <H>:refs/heads/<b>` pushes exactly the commit the guard
judged, whatever the local branch holds by then (checked 2026-10-13, git 2.55:
with the local branch one commit ahead, the remote branch was created at the
named sha). Every commit a granted push publishes is therefore a reviewed
head. The cost: no `-u`, so the local branch gets no upstream; no tool in
`tools/` or `.claude/hooks/` reads one (`grep -n -e upstream -e '@{u}'`
over both finds nothing, 2026-10-13). What stays open is a push the guard never sees
(§3.4).

### 3.3 The verdict: recorded by a hook from `@reviewer`'s handback

Ola's ruling 1: the guard reads a verdict a hook copies from `@reviewer`'s
handback into the harness state, not a line a session writes in a commit.

**The handback line.** `@reviewer` ends every review handback with exactly one
line of the form

```
Review: APPROVED <40-hex sha>
Review: APPROVED WITH WORDINGS <40-hex sha>
Review: CHANGES REQUESTED <40-hex sha>
```

matched by `^Review: (APPROVED|APPROVED WITH WORDINGS|CHANGES REQUESTED) ([0-9a-f]{40})$`
on a line of its own. `<sha>` is the head of the range it reviewed.

**The hook, `.claude/hooks/record_review.py`** (new, standard library, imports
`harness_mode`; under `.claude/hooks/`, so already governed). Two events call
it (the entries are in §3.8):

- `PreToolUse` with matcher `SubagentHandback`: the handback is the tool's
  `tool_input.message`. This is how a persona hands back in this harness (the
  subagent transcripts show the `SubagentHandback` tool call, followed by a
  short closing text that is not the report).
- `SubagentStop` with matcher `reviewer`: the handback is
  `last_assistant_message`, for a run whose report is its last text (a plain
  `claude -p` subagent).

The field that proves the agent is `@reviewer` is the event's `agent_type`,
which Claude Code sets only for an event fired inside a subagent, to the
subagent's type, and which no agent can set. Checked 2026-10-13 on the two
versions in use (the main session's transcript records 2.1.289; `claude
--version` prints 2.1.296): the harness queue logs carry
`agent_type: "reviewer"` on `PreToolUse` events of `@reviewer`'s own Bash and
Edit calls (`<git-common-dir>/harness/queue-*.jsonl`, 5 `guard_push` and 9
`guard_governance` lines, the latest 2026-10-12); a probe on 2.1.296 in a
scratch directory (`claude -p
--settings <scratch file> --agents …`, a subagent of type `reviewer`) logged
`SubagentStop` with `agent_type: "reviewer"`, `agent_id` and
`last_assistant_message` holding the subagent's reply, and `PreToolUse` of
the main session's `Agent` call with no `agent_type`. **Not checked:** that
`PreToolUse` fires for the `SubagentHandback` tool, since no hook here
matches it yet. The first step of the green phase is that check, once Ola
has added the §3.8 entries: one review handback, then `reviews.jsonl` has
its line. If the event does not fire, nothing is recorded and every grant
is refused (fail closed); the design then falls back to `SubagentStop` with
the reviewer told to repeat the line as its last text.

The hook writes one line to `<state>/reviews.jsonl` (`{"at", "sha",
"verdict", "agent_id", "session_id", "event"}`, one write with `O_APPEND`,
as `harness_mode.append_queue`) only when all hold, and otherwise writes
nothing:

1. `agent_type` is present and equals `"reviewer"` (absent, or another
   persona: nothing; fail closed);
2. the handback text has exactly one line matching the pattern (none, two or
   more, or a malformed one such as a short sha or an unknown verdict:
   nothing).

It writes in every mode (off, on, broken), so a review by day counts at
night. It never blocks: it prints nothing and exits 0, also on an error
(then nothing is recorded). The records are never pruned; a sha not on the
branch is ignored by the rule below.

**`reviewed(H)`** holds when all of these hold (git calls only; no network):

1. `L` is `git rev-list --first-parent refs/remotes/origin/master..H`: `H`'s
   first-parent chain back to where it leaves master.
2. `V` is the newest record of `reviews.jsonl` (last in the file) whose `sha`
   is in `L`. None: "no review recorded for this branch". Records for other
   branches and for master's own commits are not in `L`, so they never count.
   A `Review:` line or trailer in a commit message is not read.
3. `V`'s verdict is `APPROVED` or `APPROVED WITH WORDINGS`. Otherwise: "last
   review round: CHANGES REQUESTED".
4. `R` is `V`'s sha. Every commit `C` on the first-parent chain after `R` up
   to `H` is one of:
   - **a prose commit**: one parent, and every path of
     `git diff --name-only --no-renames C^ C` is prose by
     `tools/ci_changes.py`'s `is_prose` (Markdown at the root or under
     `docs/`, less `NOT_PROSE`); the main session's recording commit is one;
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

### 3.4 What the guard trusts

Each item below, if wrong or changed, lets the guard pass an act it should
refuse. None is authenticated; the guard is a tripwire against pushing an
unreviewed or since-changed head, not proof.

| Trusted | Who can change it | What it would let through |
|---|---|---|
| The grants in the flag `<state>/unattended.json` | Ola, through `away.py`; `guard_governance` refuses an agent's write to `<state>/` it can read, but not a write from inside an agent-written script it cannot | a branch Ola did not name |
| `<state>/reviews.jsonl` | `record_review.py`; the same gap as the flag | a verdict `@reviewer` never gave |
| `<state>/grants-used.jsonl` | `guard_push.py`; the same gap | a second enqueue of one head |
| `agent_type` in the hook event | Claude Code only | a verdict from another persona |
| What `@reviewer` writes on its `Review:` line | `@reviewer`, from its brief | a wrong sha, or a verdict it did not reach |
| `refs/remotes/origin/master` | `git fetch`; `git update-ref` (`guard_push` asks before it); a direct write to `.git/refs/` or `.git/packed-refs`, which nothing governs | moved to a code branch's tip, a merge of that branch passes as "a master merge"; moved back, the walk of `reviewed` grows (records on master's own commits then count) |
| `refs/heads/<b>` | any writer in the worktree | nothing beyond the reviewed rule: the guard judges the sha it reads, and the push sends that sha |
| `tools/ci_changes.py`'s `is_prose`, `tools/push_grant.py`, the hooks, `.claude/settings.json` | governed writes (asked by day, refused while away) | anything |
| `git` and `gh` on `PATH` | the machine's owner | anything |
| The forge's `expectedHeadOid` check | GitHub | an enqueue of a head other than the one judged |
| Only the `Bash` tool pushes | — | a push from a program the guard does not see (a script, a Python `subprocess`) is not judged at all; while away, nothing but Ola's grant makes the guard pass a push, but a push it never sees needs no pass |

A post-review prose commit is accepted whether or not the reviewer named it
(§1). The `--back` report (§3.6) names every granted act and the review
record it relied on, so Ola checks the night against it.

### 3.5 Guard wiring (guard_push.py, push_grant.py)

`tools/push_grant.py`, new, standard library only, imports `harness_mode`,
`shell_scan` and `ci_changes`. One public entry point:

```python
def judge(event: dict[str, Any], command: str, reasons: list[str]) -> Grant | None
```

`None`: not a grant-shaped act for a granted branch while away; the guard goes
on as today. `Grant(passed: bool, branch: str, act: str, head: str | None,
why: str)`: a granted branch's act, passed or not, and the reason in plain
words. Internals, each testable alone: `granted_act(command) -> Act | None`
(pure, §3.2's forms over the command text, with the expansion test), `grants(flag) ->
tuple[str, ...]` (pure), `parse_review(text) -> tuple[str, str] | None`
(pure, shared with `record_review.py`: the verdict and sha of the one
`Review:` line, or None), `reviewed(repo, head, records) -> str | None`
(None or the reason, §3.3), `jsonl(path)`; `review_of` also returns the record relied on.

`judge` returns `None` unless the mode is `on` and `reasons` is `[PUSH]` or
`[PR]` (guard_push's reason strings: no force, no `--no-verify`, no history
rewrite, no unknown command).

In `guard_push.main`, after `reasons` is computed and before `emit`:

- `judge`, before returning `Grant(passed=True)`, appends one line to `<state>/grants-used.jsonl`
  (`{"at", "window": flag since, "branch", "act", "head", "review": [V's sha,
  V's at]}`, one
  write, as `harness_mode.append_queue`), then print nothing and return: the
  act passes. If the line cannot be written, `judge` returns `Grant(passed=False)` instead (an enqueue
  that cannot be counted is not granted).
- `Grant(passed=False)`: `why` gains "; granted branch <b>, but <reason>", so
  the refusal and its queue line say why the grant did not apply.
- Any exception in `push_grant` falls to the existing `except`: deny.

### 3.6 `away.py --back`

After the queue and before the `ASK OLA` lines:

```
Push grants (2):
  worktree-threads: pushed 1a2b3c4 at 23:41 UTC; PR opened; enqueued 1a2b3c4 at 23:44 UTC (review of 5d6e7f8, recorded 23:30 UTC)
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
> pushed (as `git push origin <sha>:refs/heads/<branch>`), opened as a pull
> request and enqueued (once per head) without a fresh yes, provided
> `@reviewer`'s last verdict on the branch, recorded by a hook from its
> handback, approved its head; the guard checks this, and the grant ends
> with the window.

*A fresh yes, every time* (line 121), after "one yes covers one enqueue":
", except a branch Ola granted for an away window (*Unattended mode*)".

*`@reviewer` runs … again before a push that follows a round of findings*
(line 134), after "round of findings": " (a verdict other than `APPROVED`,
`docs/increments/README.md` step 4)".

`docs/increments/README.md`, step 4, after "Check: `grep -n '^## Review'
docs/increments/NN-name.md` on the branch.":

> `@reviewer`'s handback ends with one line, `Review: <verdict> <sha>`: its
> verdict and the full sha of the head it reviewed; a hook records it, and
> the away-mode push grant reads only that record. When `@reviewer`'s only
> findings are wordings in prose files (Markdown at the root or under
> `docs/`) that it names line by line with the new text, its verdict is
> `APPROVED WITH WORDINGS`. The main session applies them exactly, in the
> recording commit, and a branch whose only change after the review is that
> commit counts as reviewed: no second round. A wording in a code file, or
> any other finding, is `CHANGES REQUESTED`.

`.claude/agents/reviewer.md`, item 1 of its "Feedback format" list (line 55), replacing
"`APPROVED` or `CHANGES REQUESTED` (with explicit blocking issues)." by:

> `APPROVED`, `APPROVED WITH WORDINGS` (only prose wordings, each named line
> by line with its new text; `docs/increments/README.md` step 4), or
> `CHANGES REQUESTED` (with explicit blocking issues). End the handback with
> exactly one line `Review: <verdict> <full sha of the head you reviewed>`;
> a hook records it, and no other line may start with `Review:`.

The frontmatter `description` keeps "returning APPROVED or CHANGES
REQUESTED"; unchanged, it is still true of the two outcomes.

### 3.8 The settings entries (Ola's approval, by day)

Two entries in `.claude/settings.json`, added by Ola or at his prompt; the
rest of the file is unchanged. In the `"PreToolUse"` array, a new element:

```json
{
  "matcher": "SubagentHandback",
  "hooks": [
    { "type": "command", "command": "$CLAUDE_PROJECT_DIR/.claude/hooks/record_review.py" }
  ]
}
```

and a new top-level key beside `"PostToolUse"`:

```json
"SubagentStop": [
  {
    "matcher": "reviewer",
    "hooks": [
      { "type": "command", "command": "$CLAUDE_PROJECT_DIR/.claude/hooks/record_review.py" }
    ]
  }
]
```

The hook checks `agent_type` itself; the `reviewer` matcher only saves a
start for other personas, and the `SubagentHandback` matcher has none (the
tool's event carries the persona only in `agent_type`).

## 4. Tests @tester writes first

All in the harness suites (`-m harness`), on `harness_fixtures.make_repo`
copies; `COPIED` gains `tools/push_grant.py`, `tools/ci_changes.py` and
`.claude/hooks/record_review.py`. The fixture repository for §4.2-4.4
(measured by construction, not sampled): `master` with 1 commit;
`refs/remotes/origin/master` set to it with `git update-ref` (no remote
needed); branch `feat` with a code commit `X` changing `tools/x.py` (1 line),
then a recording commit changing `docs/increments/h99.md` (1 line, no
trailer); `<state>/reviews.jsonl` with 1 record, `APPROVED` at `X`. Flag on,
`grants: ["feat"]`, unless a test says otherwise. `H` is `feat`'s head.

### 4.1 `away.py` (test_away.py)

1. `8h --grant feat --grant other` writes `"grants": ["feat", "other"]` to
   the flag and to the `windows.jsonl` enter entry; the tty text lists both.
2. `8h` alone writes `"grants": []`; a flag with no `grants` key (written
   before h20) reads as no grants.
3. `--grant master`, `--grant main`, `--grant HEAD`, `--grant 'a..b'`,
   `--grant '~x'`, `--grant '=x'`, `--grant 'a*'`: exit 2, tty never opened,
   no flag.
4. `--back --grant feat` and `--grant feat` alone: exit 2.
5. A grant for a branch that does not exist: accepted, tty text says "does
   not exist yet".
6. Re-entering with `--grant b` only: the flag's grants are `["b"]`.
7. `--back` with a 3-line `grants-used.jsonl` (2 lines of this window, 1 of
   an earlier one): prints this window's 2 acts under `Push grants`, not the
   third; an unused grant prints its head and `reviewed: yes` or the reason;
   no grants prints `Push grants: (none)`.

### 4.2 The forms (`granted_act(command)`, pure, parametrised)

8. Each form of §3.2 parses, including `-C /abs/dir`, `<H>:refs/heads/feat`,
   `-H`, `--head=feat`, `-B master`, `--title 'a title'`, `--body-file f`,
   `--merge`.
9. Each of these is `None`: `--force`, `-f`, `--force-with-lease`,
   `--no-verify`, `-u`, `origin feat` (no sha), `+<H>:refs/heads/feat`,
   `<H>:feat`, `<H>:refs/heads/main`, `origin HEAD`, `upstream
   <H>:refs/heads/feat`, a second refspec, `--admin`, `-d`, `--auto`,
   `--squash`, `--rebase`, `-R x/y`, `--body text`, `--base other`, a short
   sha (39 hex), a missing `--match-head-commit`, `gh pr merge 241` (a
   number, not the branch); and the expansions: `-C ~/other`, `-C =other`,
   `-C rel/dir`, `-C '$HOME/x'`, `--head 'fea*'`, `--body-file ~/f`,
   `--body-file==ls`, `--body-file=~/f`,
   `--title t` (not in single quotes), `--title "t"`.

### 4.3 The guard end to end (test_guard_push.py, through the hook copy)

10. Flag on, `feat` granted and reviewed: `git push origin
    <H>:refs/heads/feat`, `gh pr create --head feat --title 't' --body-file f`
    and `gh pr merge feat --match-head-commit <H> --merge` each print no
    decision; three lines in `grants-used.jsonl`; nothing in `queue.jsonl`.
11. The same three with the mode off: ask, as today; nothing logged.
12. With the flag broken: deny, as today.
13. Not granted (`grants: ["other"]`): deny and queue, `why` without the
    grant clause.
14. Each wrapper is refused: `env A=1 git push origin <H>:refs/heads/feat`,
    the same with `2>&1`, with `&& echo`, inside `sh -c '…'`, and after
    `timeout 5`.
15. A second `gh pr merge feat --match-head-commit <H>` in the window: deny,
    `why` says the head was enqueued already. After a clean master merge
    (new `H2`): one enqueue at `H2` passes.
16. A sha other than `H` as the push's source or as `--match-head-commit`:
    deny.
17. `git -C <absolute path of another repository> push origin
    <H>:refs/heads/feat`: deny; `git -C ~/other …` and `git -C =other …`,
    with `~/other` and `=other` set up to reach the guard's own repository
    (a `HOME` and a `PATH` entry pointing there): deny.
18. `grants-used.jsonl` unwritable (a directory at its path): deny.

### 4.4 `reviewed` (the invariant-critical suite of this increment)

Each case builds on the fixture; the expected result is pass or the named
reason.

19. As built: pass.
20. No record in `reviews.jsonl`, or only a record of a sha on another
    branch: "no review recorded for this branch".
21. Two records on `feat`'s chain, `APPROVED` then (newer, last in the file)
    `CHANGES REQUESTED`: refused; the same two in the other order: pass.
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
30. A branch `other` off master has a code commit `S` (`tools/s.py`, 1
    line); `feat` merges `other` before `X`, so `S` is reachable from `H`
    but not on its first-parent chain. `reviews.jsonl`: `APPROVED` at `X`,
    then a newer `APPROVED` at `S`. Expected: **pass** (the first-parent
    walk does not see `S`, picks `X`'s record, and only the prose recording
    commit follows `X`). A full `origin/master..H` walk picks `S` and
    refuses on `X` ("commit X changes tools/x.py after the review of S").
31. Master gains commit `M` (code), recorded `APPROVED` at `M`, and `feat`
    merges master; no record on `feat`'s own chain: "no review recorded for
    this branch" (the first-parent walk).
32. A code commit `Y` after the recording commit, whose message ends with
    the line `Review: APPROVED <Y>`: refused, naming `Y` and `tools/x.py`
    (a line in a commit message is not read).

### 4.5 The recording hook (test_record_review.py, through the hook copy)

Each event is a `PreToolUse` event for `SubagentHandback` with
`tool_input.message` set, unless a test says otherwise; the expected result
is the lines of `reviews.jsonl`. Every case prints nothing and exits 0.

33. A handback ending `Review: APPROVED <40 hex>` with `agent_type`
    `"developer"`, and the same with no `agent_type`: no line. With
    `"reviewer"`: one line, with mode off, on and broken alike; the same
    from a `SubagentStop` event's `last_assistant_message`: one line.
34. Malformed lines from `"reviewer"`: a 39-hex sha, upper-case hex,
    `Review: APPROVED`, `Review: LGTM <sha>`, `Review: APPROVED <sha> extra`,
    the line indented: no line.
35. Two `Review:` lines from `"reviewer"`, the same or different: no line.
36. `reviews.jsonl`'s directory unwritable: no line, exit 0, nothing
    printed.

**Mutation targets** (the kill record must cover each): the first-parent
walk replaced by a full `origin/master..H` walk (kills 30); "newest
record" replaced by "first record" (21); `is_prose` replaced by "any `.md`"
(25); the second-parent ancestry check dropped (29); the merge-tree diff
replaced by "the merge is clean" (27 or 28, one of them); the verdict test
accepting any verdict (21); the `agent_type` test dropped, or passing when
absent (33); "exactly one line" replaced by "the first line" (35); the
expansion test dropped (9, `--body-file ~/f`; test 17's cases fall to the absolute-path rule); the push's source-sha test dropped (16); the
one-command test dropped (14); the per-head enqueue count dropped (15).

## 5. Size

Against the 700-line limit (`CLAUDE.md` §2), net production lines:

| File | Estimate |
|---|---|
| `tools/push_grant.py` (new; forms with the expansion test, `parse_review`, `reviewed` over records) | about 155 |
| `.claude/hooks/record_review.py` (new) | about 40 |
| `tools/away.py` (`--grant` with the name check, flag key, tty text, `--back` section) | about 40 |
| `.claude/hooks/guard_push.py` (call `judge`, log, extend `why`) | about 15 |
| `.claude/hooks/guard_governance.py` (`GOVERNED` gains `tools/push_grant.py` and `tools/ci_changes.py`: the push guard and the recording hook import them, so an edit to either is an edit to a hook) | 2 |
| Total | about 250 |

`pyproject.toml`'s mypy `files` gains `tools/push_grant.py` and the settings
entries of §3.8 are JSON (neither counted). Rule text, §3.7: about +170 words
over four files.

## 6. Order

Design review; red (tests of §4, all failing); green, with
`record_review.py`; Ola adds the §3.8 entries (a settings change: his yes,
by day); the live check of §3.3 (one `@reviewer` handback, then a line in
`reviews.jsonl`), recorded in this file; code review, whose own handback the
hook records; the rule text of §3.7 applied by the main session at Ola's
prompts, in the same PR. The guard's edit, the new hook and the new
`GOVERNED` entries are governed writes, so the round runs by day.

## 7. Questions for Ola

Both answered (see *Ola's rulings*): question 1, who vouches for the review
verdict, is the hook of §3.3; question 2, one enqueue per reviewed head, is
§3.2.

## Ola's rulings

Ola, 2026-10-13, on section 7's two questions: "We use the safer option on 1, yes to 2." Question 1: the guard reads a verdict a hook copies from `@reviewer`'s handback into harness state, not a commit line. Question 2: one enqueue per reviewed head.

## Review

**Design review, round 1, 2026-10-13.** Range `97c10c4d..fa6fdaa6`, 0 counted lines, @reviewer: CHANGES REQUESTED. Swap the Review commit line for a hook-recorded verdict (Ola's ruling 1: hook event and persona field named and checked on this Claude Code version, fail closed; one `Review: <verdict> <40-hex sha>` line in the handback; `reviews.jsonl` in every mode; `reviewed(H)` the newest record on H's first-parent chain within origin/master..H; tests 20, 21, 31 replaced and four added; size about 250), correct the auto-merge race claim in section 3.2 (gh turns on auto-merge when checks are pending and GitHub keeps it after a push by anyone with write access), refuse shell-expanded words (`-C ~/x`, `-C =x`; grant names `[A-Za-z0-9._/-]+`), list everything the guard trusts in section 3.4 (including `refs/remotes/origin/master` and the harness state), and fix line 36's pointer to section 3.4.

**Design review, round 2, 2026-10-13.** Range `97c10c4d..452b1b28`, 0 counted lines, @reviewer: CHANGES REQUESTED. Round 1's five blockers fixed; the hook design accepted (fails closed; SubagentStop alone is not enough, as the report is in the handback tool's input). The expansion-test kill moves from test 17 to test 9 (`--body-file ~/f`); test 30 rebuilt so the full `origin/master..H` walk is caught (feat merges a branch `other` whose commit S has a newer APPROVED record; expected pass; the mutant picks S and refuses); `--opt==x` refused (zsh's `magic_equal_subst`), added to test 9. Suggestion: section 3.7's "check list item 1" is item 1 of reviewer.md's "Feedback format".

**Design review, round 3, 2026-10-13.** Range `fb2d71e6..8b403eff`, 0 counted lines, @reviewer: APPROVED. Round 2's three items and its suggestion are applied: test 9's `--body-file ~/f` is the expansion mutant's kill, test 30 is rebuilt so the full-walk mutant is caught, `--opt==x` and `--opt=~x` are refused, and section 3.7 points to reviewer.md's "Feedback format" item 1. Suggestion: test 30's quoted reason should say "refuses" without naming X, since the full-walk mutant may refuse the merge of `other` first.

**Ruling, 2026-10-13 (Ola, "1, always between the two parts").** In the report `away.py --back` prints, the push-grants section always sits between the queued acts and the ASK OLA lines, "Push grants: (none)" included. The green step `4ad42e45` put "(none)" after the ASK OLA lines to keep `test_back_prints_the_queue_by_branch_and_archives_it` unchanged; @tester changes that test's expected lines to include "Push grants: (none)" before the ASK OLA heading, and @developer moves the line.

**Code review, round 1, 2026-10-13.** Range `97c10c4d..a5df3a0b`, 320 counted lines, @reviewer: APPROVED WITH WORDINGS. Three design-text wordings (status line, `log_use`, where grants-used is written), applied by the main session; kill record complete (12 targets, the one-command mutant equivalent through `unexpandable()`); the shell-written governed files read in full and match the design. Merge waits on the §3.8 settings entries, the live check and the §3.7 rule text, which need Ola. Suggestions: unpack away.py's `done = [...]` and push_grant `changes()`'s merge test; wrap `grant_lines` in `--back` so a git failure cannot stop the queue being archived.
