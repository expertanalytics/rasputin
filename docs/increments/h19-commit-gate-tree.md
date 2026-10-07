# Harness h19: the after-commit gate checks the tree that was committed

Status: design review round 1 answered on `worktree-commit-gate` off
`5c41238`; next, design review round 2. Not yet red. One PR, about 45 net
production lines (§6), in one governed file,
`.claude/hooks/gates_after_commit.py`, plus one new test file. Ola's ruling
in §7.

What this is. Row T1 of the 2026-10-06 day retrospective
(`docs/retrospectives/2026-10-06-day-bottlenecks-and-merges.md`, table
"Proposals"). Ola approved it twice on 2026-10-07:

- 05:42 UTC, "1-3 default. 4 must wait", answering the main session's list
  of 05:36 UTC, whose item 1 was 20c question 7, item 2 was 20c question 6,
  and item 3 was the 2026-10-06 retrospective: push it, open its PR, and
  its proposals T1-T6, R1, R2 and S1, each with the default yes. So the one
  answer is both the 20c ruling that `20c-soft-quality.md` records (on its
  questions 7 and 6) and the approval of T1 (checked against the session
  transcript).
- 08:39 UTC, "Both binary by default. T1 and guard PR B now/", which asked
  for this work now.

Labels used here:

| Label | What it is |
|---|---|
| the hook | `.claude/hooks/gates_after_commit.py`, run by Claude Code after every Bash call (`PostToolUse`, `.claude/settings.json`) |
| the gates | the five checks the hook runs: `check_prohibited_deps.py`, `check_detria_boundary.py`, `check_citations.py`, `ruff check .`, `ruff format --check .` |
| D3 | finding D3 of that retrospective: a master merge that dropped an import was gated in the main checkout, not in the worktree it was made in, and passed |
| PR B | the guard fixes of `docs/increments/h16-harness-fixes.md`, on `worktree-h16b` |
| h18 | `tools/new_worktree.py` and `tools/merge_master.py`, on `worktree-harness-tools` |
| subagent-shaped event | a hook input as a subagent's Bash call produces it: `agent_id` and `agent_type` set (a session started with `--agent` sets `agent_type` too; only `agent_id` is subagent-only), `cwd` the session's directory (the main checkout), the worktree named only inside the command |

## 1. The two faults, checked

Both checked at `5c41238` by loading the hook from its path and feeding it
a subagent-shaped event (`cwd` = `/Users/skavhaug/projects/rasputin`):

- **Wrong tree.** For `cd <this worktree> && git commit -m x`,
  `committed_tree(event)` returns `/Users/skavhaug/projects/rasputin`, the
  main checkout. The hook takes the tree from `cwd` alone, and a subagent's
  shell returns to the session directory between calls, so it works in a
  worktree only through `cd <wt> && ...` inside the command. D3 counts about
  7,900 subagent hook events in one session, every one with that `cwd`.
- **Missed trigger.** The hook fires only when the command text contains
  `git commit` or `git merge`. `git -C <wt> commit -m x` and
  `git -c a=b commit -m x` both exit 0 in 0.04 s (no gate ran; a run with
  the gates takes about 11 s, §3.3). h18 §4.7 tells a merging persona to
  type `git -C <wt> commit -F <file>` by hand, so this route is the one its
  own instructions use.

The text test also fires where nothing was committed: `echo "git commit"`,
`grep -n "git merge" x`, `git log --grep="git commit"`, each costing a full
gate run.

## 2. Prior art: legacy and literature

*Literature.* Git's own hooks solve "which tree": "Before Git invokes a
hook, it changes its working directory to either $GIT_DIR in a bare
repository or the root of the working tree in a non-bare repository"
(git-scm.com/docs/githooks, read 2026-10-07). A `post-commit` or
`post-merge` git hook would therefore get the tree for free. It is not used
here, for three reasons: it lives in `.git/hooks` or `core.hooksPath`,
which is per clone and untracked (and a `git config` write is guarded); it
"is meant primarily for notification, and cannot affect the outcome", and
its output lands inside the Bash result, where it reads as the command's
own output, not as the hook's exit-2 feedback the agent is told to act on;
and `post-merge` "is not executed, if the merge failed due to conflicts",
nor does `post-commit` cover a merge. So the design keeps the Claude Code
hook and recovers the directory the way the shell would: by following the
command's own `cd` and `git -C`. Claude Code's hook input carries `cwd`,
"Current working directory when the hook is invoked"; `agent_id`, only
when the hook fires inside a subagent; and `agent_type`, inside a subagent
or in a session started with `--agent` (code.claude.com/docs/en/hooks, read
2026-10-07). No novelty is claimed.

*Legacy.* Nothing; the legacy tree has no harness or git hooks.

```
$ git grep -l -i -e "post-commit" -e "after.commit" -e "show-toplevel" -e "PostToolUse" legacy-archive -- legacy
(no files; exit 1)
```

## 3. Design

All in the hook. It already imports nothing of the repository; it gains
`shell_scan` by `sys.path.append(str(<repo>/tools))` and the import in
`try`, `shell_scan = None` if it fails. **Append, never
`sys.path.insert(0, ...)`**: put first, a `tools/` file named after a
standard-library module replaces it for every later import, the hole PR B's
fix G4 closes in the four guards and three tools (the `insert(0, ...)` that
`guard_push.py` still has at `5c41238`). `tools/shell_scan.py` is not
changed (§4).

### 3.1 Which commands trigger it

A Bash call triggers the gates when `shell_scan.parse(command)` finds a
simple command that is git (`base(argv[0]) == "git"`, after `shell_scan`'s
wrapper stripping, so `env`, `timeout` and `/usr/bin/git` count) whose
subcommand is in

```python
COMMITTING = frozenset({"commit", "merge", "pull", "cherry-pick", "revert", "am", "rebase"})
```

the subcommands that make commits on the checked-out branch. The
subcommand is the first word after git's own options; an option in
`GIT_TAKES_ARG = {"-C", "-c", "--git-dir", "--work-tree", "--namespace",
"--attr-source", "--config-env"}` (the set `guard_push.py` uses on PR B)
consumes the next word; the glued form (`--work-tree=/a`) is one word and
consumes none. The set is copied, not imported from
`guard_push.py`, so that one hook failing to load cannot take the other
with it.

Commands inside `$(...)`, backticks and `sh -c` count: `shell_scan` parses
them. Words inside quotes and heredocs do not, so the three false triggers
of §1 stop.

When `parse` returns `None` (a line it cannot read), the hook falls back to
a text test, widened to options:
`re.search(r"\bgit(\s+-\S+(\s+[^-\s]\S*)?)*\s+(commit|merge|pull|cherry-pick|revert|am|rebase)(?![\w-])", command)`.
The end is `(?![\w-])`, not `\b`: with `\b` it matches `git commit-tree`
and `git merge-base` (checked: the `\b` form matches `git commit-tree T`,
this form matches neither, and matches `git -C /a commit` and `git -c a=b
commit -m "x`).

`commit-tree` is not in the set: it writes an object, not a branch. D3's
hand-built merge ended in `git merge --ff-only`, which is.

### 3.2 Which tree

Walk the parsed simple commands in order, keeping a current directory that
starts at the event's `cwd` (or the hook's own checkout if `cwd` is
missing, as today):

| Command | Current directory becomes |
|---|---|
| `cd D` or `pushd D` (options such as `-P` skipped) | `D` if absolute, else current / `D`; `~` expanded |
| `cd`, `cd -`, `popd`, or a `D` containing `$` or a backtick | unknown, holding the word as written |
| a relative `cd D` while the current directory is unknown | stays unknown, holding the first unknown word |
| a committing git | not changed; its own directory is the current one with each `-C X` applied in order, by the rows above (so a relative `-C X` on an unknown directory stays unknown) |
| a committing git with `--git-dir` or `--work-tree`, as two words (`--work-tree /a`) or glued (`--work-tree=/a`) | that command's directory is unknown, holding the option as written |

Each committing git command yields a directory or an unknown word. A
directory maps to `git rev-parse --show-toplevel` run there (the existing
`git()` helper). A linked worktree's toplevel is the worktree. The trees
are deduplicated in first-seen order.

- **Unknown word, or a directory that does not exist or is in no work
  tree:** that tree is NOT CHECKED. The gates are not run in some other
  tree in its place; running them in the wrong tree is the fault being
  fixed.
- **Unreadable line, text test matched:** the gates run in `cwd`'s
  toplevel (today's behaviour), and a NOT CHECKED line says the line could
  not be read, so that tree may not be the one committed.

Accepted limits, each rare in this repository's commands and each failing
towards a gate run, not towards silence where it can:

- A `cd` inside a subshell, `(cd x) && git commit`, is read as if it
  persisted; `shell_scan` flattens subshells. The hook would then gate `x`.
- A variable assigned on the line (`WT=/a; cd $WT`) is not substituted;
  `shell_scan` substitutes such values only into write targets, and
  extending that touches PR B's file (§4). It reports NOT CHECKED.
- A `cd` inside `$(...)`, backticks or `sh -c` leaks the same way, and
  `shell_scan` emits a substitution's inner commands before the command
  that holds it (checked at `5c41238`: `x=$(cd /b && pwd) && git commit`
  parses as `cd /b`, `pwd`, `git commit`). That line gates `/b`, while the
  commit ran in `cwd`. Wrong tree, but a gate run and the tree named in the
  output, so the reader can see it.

### 3.3 What runs there

The same five gates as today, each with `cwd` = the tree, unchanged:

- `python3 tools/check_*.py`: the tree's own copies, under the `python3`
  on `PATH`. They are stdlib only, so no venv is needed.
- `ruff`: `find_ruff(tree)` as today: the tree's `.venv/bin/ruff`, else the
  main checkout's `.venv/bin/ruff`, else `PATH`, else DID NOT RUN. Ruff
  imports none of the project's code, so a worktree without a venv (65 of
  the 109 directories under `.claude/worktrees/` have no
  `.venv/bin/python`, counted 2026-10-07) loses nothing by using the main
  checkout's; the version is pinned (`ruff~=0.16.8`, `pyproject.toml`), and
  h18's `new_worktree.py` installs the same pin.

Not added: `mypy` and `pytest`. They import `tin_engine`, so they are right
only in a venv that runs this tree's code, which most worktrees lack until
h18's `--existing` gives them one; in any other venv they check another
checkout (the "which Python" lessons of the day retrospective). `pytest`
also needs `_core` built, and takes minutes. h18's `merge_master.py` runs
both in the worktree's own venv, and CI runs both. The hook stays "fast,
no build, no network" (its `gates()` docstring).

Time, measured in this worktree at `5c41238` with the main checkout's
ruff: `check_citations.py` 9.9 s, `check_prohibited_deps.py` 0.7 s, the
other three under 0.1 s each; about 11 s per tree. A line that commits in
two trees takes two runs. The per-gate timeout stays 180 s, unchanged; at
this repository's size the slowest gate uses about a twentieth of it.

### 3.4 What it says

Silent and exit 0 when every tree's gates pass, as today. Otherwise exit 2
(stderr goes back to the agent), with each failure or NOT CHECKED line
headed by the tree it concerns, for example:

```
Gates in /Users/skavhaug/projects/rasputin/.claude/worktrees/audit-geojson:
$ .../ruff check .  (exit 1)
...
NOT CHECKED: git ran in `$WT`, which this hook cannot resolve; run the gates there yourself.
```

The closing paragraph (not the full set; `gh pr checks` decides) stays.

### 3.5 Interface for the tests

Two functions, loadable by path, pure apart from the `git rev-parse` call:

```python
def git_dirs(command: str, cwd: Path) -> list[Path | str] | None:
    """The directory each committing git command on the line runs in, in order:
    a Path where §3.2 resolves it, the word as written where it does not.
    [] when no command commits; None when shell_scan cannot read the line
    or failed to import (the hook then uses §3.1's text test)."""

def committed_trees(event: dict) -> list[Path | str]:
    """The distinct work-tree tops to gate, in first-seen order; a str is a
    NOT CHECKED entry (§3.2). [] means the hook does nothing."""
```

`committed_tree(event)` is replaced by `committed_trees`; nothing else
imports it (`grep -rn committed_tree` finds only the hook, and the
retrospective's prose).

## 4. Other branches and merge order

**PR B** (`worktree-h16b`) changes `tools/shell_scan.py` (the shell list,
`mv` as a copier), both guards, `.claude/REQUIRED-READING.md` and
`tests/python/harness_fixtures.py`. This PR changes none of them:

- It only calls `shell_scan.parse` and `base`, and reads `Simple.argv`; PR B
  keeps all three. PR B's wider shell list means `dash -c "cd x && git
  commit"` is parsed too once both are in; either order works.
- The test copies the hook and `shell_scan.py` into its temporary
  repository itself and does **not** add to `harness_fixtures.COPIED`,
  which PR B extends at the same place.
- `.claude/REQUIRED-READING.md` ("puts the gates' own output in the
  transcript after a commit or merge") stays true and is not edited; PR B
  rewrites that paragraph.

No file is shared, so the order is free. Default: **PR B first** (it is
further along: past its review round 9), this PR merged after it with a master
merge before its review, so that its tests run against PR B's
`shell_scan.py`.

**h18** (`worktree-harness-tools`) appends its `ROADMAP.md` row after the
same last row as this PR, so whichever merges second gets a one-hunk
conflict there: keep both rows, h18 first. h18 does not touch the hook.
Its `merge_master.py` commits inside the tool; the Bash line is `python3
tools/merge_master.py ...`, which still does not trigger the hook (no git
command on the line), as h18 §4.4 intends: the tool runs a superset of the
gates in the right tree.

## 5. Red tests for `@tester`

New file `tests/python/test_gates_after_commit.py`, `pytestmark =
pytest.mark.harness` (h17 §4b), no network, temporary repositories only.
Fixture: `harness_fixtures.make_repo` for a main checkout, plus
`add_worktree` for a linked worktree; copy the hook and
`tools/shell_scan.py` into the main repository; in the main repository
plant `.venv/bin/ruff` as a two-line shell script that exits 0, and stub
`tools/check_*.py` (three files) that exit 0; commit them, so the worktree
has them too. In the worktree, replace `tools/check_detria_boundary.py` with
one that prints `RED IN <its cwd>` and exits 1. Events are
subagent-shaped: `hook_event_name` `PostToolUse`, `tool_name` `Bash`,
`cwd` = the main repository, `**SUBAGENT`.

1. **The worktree is gated** (D3). `cd <wt> && git commit -q --allow-empty
   -m x` (really run in `<wt>` first): the hook exits 2 and stderr contains
   `RED IN <wt>` and names `<wt>`. Today: exit 0. The same with `git -C
   <wt> commit ...` and with `cd <wt>; git -c user.name=x merge --no-ff y`
   (a branch `y` made for it).
2. **The main checkout stays gated** for `git commit -m x` with no `cd`
   (main session shape, no `agent_id`): exit 0 here; and with the main
   repository's own detria stub made red, exit 2 naming the main
   repository.
3. **Triggers**, as `git_dirs(command, cwd)` against expected lists
   (parametrised, no gate run): `git -C /a commit`, `git -c a=b commit`,
   `git --no-pager merge x`, `/usr/bin/git commit`, `env X=1 git commit`,
   `timeout 9 git commit`, `bash -c "cd /a && git commit"`, each of
   `pull`, `cherry-pick`, `revert`, `am`, `rebase`; and `[]` for `echo "git
   commit"`, `grep -n "git merge" f`, `git log --grep="git commit"`, `git
   merge-base a b`, `git commit-tree T`, `git status`.
4. **Directories**, same form: `cd /a && git commit` gives `/a`; `cd /a &&
   cd b && git commit` gives `/a/b`; `git -C /a -C b commit` gives `/a/b`;
   `cd /a && git -C b commit` gives `/a/b`; a relative `cd b` from `cwd`;
   `cd ~/x` expands; `cd /a && git commit && cd /b && git commit` gives
   both in order; `cd "$WT" && git commit`, `cd - && git commit` and `git
   --work-tree=/a commit` give the word as written (a `str`); so do
   `cd - && cd b && git commit` and `cd - && git -C b commit` (a relative
   step from an unknown directory stays unknown), while `cd - && cd /a &&
   git commit` gives `/a`.
5. **NOT CHECKED, end to end**: `cd "$WT" && git commit -m x` exits 2,
   stderr has `NOT CHECKED` and `$WT`, and no gate ran (the stub
   gates write a marker file when run; none exists). The same for `cd
   <tmp>/nonexistent && git commit`.
6. **Two trees**: `git -C <main> commit ... && git -C <wt> commit ...`
   gates both; stderr names `<wt>`; the main checkout's stubs ran once
   (markers), so dedup holds when the same tree appears twice.
7. **Unreadable line**: `git commit -m "unclosed` makes `git_dirs` return
   `None`, and the hook gates `cwd`'s tree and prints `NOT CHECKED` with
   the words "could not be read".
8. **No trigger is fast**: `git status` and a non-Bash event exit 0 with no
   marker written.
9. **Import form**: the hook's source has `sys.path.append(` and no
   `sys.path.insert(0,` (the check PR B's `APPENDERS` test makes for seven
   other files, written here so neither PR edits the other's test file).
   And with `shell_scan` made unimportable (its copy in the temporary
   repository deleted), `git_dirs` returns `None` and `git commit -m x`
   still gates `cwd`'s tree through the text test.

Mutation targets, if `@tester` chooses to spend a round (not required: this
is not an invariant-critical suite): drop the `cd` row, drop `-C`
composition, drop dedup, fall back to `cwd` on an unknown word.

## 6. LOC estimate

`.claude/hooks/gates_after_commit.py`: about 45 net lines (the import, two
constants, `git_dirs` about 25, `committed_trees` about 10, the loop over
trees and the NOT CHECKED lines in `main` about 8, less `committed_tree`'s
4). The retrospective's "about 10 lines" covered the `cd` case only; the
trigger, `-C` and the unknown-directory rule are the rest. Tests about 200
lines, not counted (`CLAUDE.md` §2).

## 7. Ola's ruling

Question 1, closed. Asked: should a commit whose directory the hook cannot
work out (for example `cd "$WT" && git commit`) be reported back to the
agent as "not checked" (exit 2, the same channel as a red gate), or pass
silently as today? Default: report it. Ola, 2026-10-07 09:38 UTC: "Still
open for you: h19's question, telling the agent "not checked" when the hook
can't tell the folder: yes." So §3.2 and §3.4 stand as written: NOT
CHECKED, exit 2, no gate run in another tree in its place.

## Review

h19 design review round 1 (@reviewer, 5c412383..abeff234 on worktree-commit-gate, docs only, 0 counted lines): CHANGES REQUESTED, with two blocking items. B1: the design says to import `shell_scan` "the way guard_push.py imports it" (/Users/skavhaug/projects/rasputin/.claude/worktrees/commit-gate/docs/increments/h19-commit-gate-tree.md@abeff234:73-76). At the design's own base that way is `sys.path.insert(0, ...)`. PR B's fix G4 removes exactly that form as a live hole (/Users/skavhaug/projects/rasputin/.claude/worktrees/h16b/.claude/hooks/guard_push.py@c83dfe0e:32-35). The design must say `sys.path.append`. B2: the quote "1-3 default" (docs/increments/h19-commit-gate-tree.md@abeff234:10) is recorded in the tree only as Ola's answer on 20c questions 7 and 6 (/Users/skavhaug/projects/rasputin/docs/increments/20c-soft-quality.md@5c412383:1912-1913). Either check it against the transcript as an answer on the retrospective's questions, or drop it and keep the verbatim "T1 and guard PR B now". Everything else I checked holds. 6 suggestions: S1 Ola's "no" answer to question 1 spelled out, S2 unknown directory stays unknown, S3 `=` forms, S4 nested `cd` leaks, S5 `git_dirs` without shell_scan, S6 two wording fixes.

`@architect` answer to round 1: B1 taken, §3 names `sys.path.append` and why not `insert(0, ...)`, with test 9; B2 checked against the session transcript (the 05:36 UTC list's item 3 was the retrospective with T1-T6, R1, R2, S1), both quotes now in the opening; S1 moot, Ola ruled yes at 09:38 UTC (§7); S2 taken (§3.2 row, test 4 cases); S3 taken (§3.1, §3.2); S4 taken as a limit, checked with `shell_scan.parse` at `5c41238`; S5 taken (docstring, test 9's second half; the `APPENDERS` idea done as test 9 so neither PR edits the other's file); S6 taken (`agent_type` under `--agent`, PR B past round 9). The round 1 line above is word for word but one change: its short citation of this file gains the `docs/increments/` prefix, so that `check_citations.py` resolves it.
