# Harness h8: the unattended window, the recap and the size table

Status: **design, ruled, round 1 fixes**, @architect, 2026-10-03; Ola
answered §7's four questions the same day. One PR, by day: every production
file it changes is governed. Implements Ola's rulings of 2026-10-03 on items
1 to 6 of `docs/retrospectives/next.md`, section "The window of 2026-10-02, the restart
of 2026-10-03, h7" (recorded in #155), and his request of the same day for a
size table in the recap. Evidence for items 1 to 6:
`docs/retrospectives/2026-10-03-night-restart-h7.md` §1a to §1d and §2-§3.

## 1. Prior art

Tooling; no novelty claimed.

*Literature.* The size table answers Ola's concern that rule text only grows.
Anthropic's guidance for Claude Code says the same of its own instruction
files: "target under 200 lines per CLAUDE.md file. Longer files consume more
context and reduce adherence", and "The more specific and concise your
instructions, the more consistently Claude follows them"
(code.claude.com/docs/en/memory, read 2026-10-03). Nothing there measures
growth over time; the table does that, in words, against the last
retrospective. The idle measure is the method §1a of the evidence file used
by hand (last commit on any ref inside the window), made a print.

*Legacy.* Nothing; the legacy tree has no harness.

```
$ git grep -l -iE 'session_state|away\.py|ask ola|current-task|shell-snapshots' legacy-archive -- legacy
(no output, exit 1)
```

## 2. Scope

| # | Ruling or request | Where |
|---|---|---|
| 1 | Fallback rule before Ola leaves; `--back` prints the longest stretch without a commit | `REQUIRED-READING.md`, `tools/away.py` |
| 2 | `--back` and the recap find the main checkout from the common git dir and list `ASK OLA:` lines from every worktree | `tools/session_state.py`, `tools/away.py` |
| 3 | A line counts only if it starts, after an optional bullet, with `ASK OLA:`; an empty one warns | `tools/session_state.py` |
| 4 | `session.md`: one `NOW:`, one `QUEUE:`, one `ASK OLA:` per decision, nothing else; the recap warns on any other line | `REQUIRED-READING.md`, `tools/session_state.py` |
| 5 | List running background jobs before starting any; the recap prints them | `REQUIRED-READING.md`, `tools/session_state.py` |
| 6 | After merging a change to `CLAUDE.md` or `.claude/agents/`, restart before spawning the changed persona; until then the brief says to read the persona file from disk | `REQUIRED-READING.md` |
| 7 | The recap prints words per agent, skill and rule file, with the change since the last retrospective | new `tools/rule_sizes.py`, `tools/session_state.py` |

Out: item 17 of "The merges of 2026-10-03" (a fallback names the files it
writes; the recap flags governed ones) is not ruled. Where it would go is in
§3.4 and §6. Items 7 and 8 of the first section and 9 to 16 of the second are
not ruled either and are not touched; 9(a) (`/compact` re-reads the project
`CLAUDE.md`) would add a cheaper route to item 6's rule, but not for persona
files.

## 3. Behaviour

The three tools keep their present boundaries: `away.py` is Ola's, prints, and
imports from `session_state.py` only inside a `try` (the recap is optional to
`--back`, h3 §3.3); `session_state.py` prints the recap and never fails as a
whole, each new section costs one line on a fault, as the harness section does
now (h3 §3.8); the new `rule_sizes.py` is pure counting plus one git read, and
`session_state.py` imports it inside a `try` too. Every function below that
reads no process state is pure and is tested on strings.

### 3.1 The main checkout and every worktree (item 2)

In `session_state.py`:

- `main_checkout(repo: Path) -> Path`: from
  `git -C repo rev-parse --path-format=absolute --git-common-dir`; when that
  directory's name is `.git`, its parent; otherwise (bare repository, git
  failure) `repo` itself.
- `checkouts(repo: Path) -> list[Path]`: the `worktree <path>` lines of
  `git -C repo worktree list --porcelain`, main checkout first, keeping only
  paths that exist as directories (a pruned worktree is skipped, not an
  error). On git failure, `[repo]`.

What reads which checkout:

| Recap part | Read from |
|---|---|
| Last landed, In flight | the checkout the script runs in (unchanged) |
| `session.md` and subagent files (`print_current_task`) | the main checkout |
| Waiting on Ola | every checkout in `checkouts()` (§3.2) |
| `session.md` format warnings | the main checkout's `session.md` only (§3.3; item 4 ruled that file) |
| Predecessor turns | unchanged: the transcript folder of the checkout the script runs in. Not ruled; merging the main checkout's folder in would print the running main session's own turns as pending |

The `.claude/current-task` path stops being an import-time constant derived
from `__file__`; it is computed in `main()` from `main_checkout()`, so tests
can point it at a fixture repository.

`away.py --back` takes its decisions from `all_decisions(root)` (below), so
run from a worktree it lists the main checkout's lines, which is the fault
§1b of the evidence file reproduced. Its heading becomes
`ASK OLA lines, every worktree's .claude/current-task/:`.

### 3.2 `ASK OLA:` matching (item 3)

```python
ASK = re.compile(r"^[ \t]*(?:[-*+][ \t]+)?ASK OLA:(.*)$")   # case-sensitive
```

- `pending_decisions(tasks: Path) -> list[str]` keeps its signature. A line
  counts only on a full `ASK.match`. Not counted: `ask ola: x` (case), `Note:
  ASK OLA: x` (not at the start), `**ASK OLA:** x`, `1. ASK OLA: x` (a
  number is not a bullet), `ASK OLA x` (no colon).
- A counted line whose remainder is blank (`ASK OLA:`, `- ASK OLA:   `)
  is listed as `<file>:<lineno>: WARNING: empty ASK OLA line` in place of the
  line. Others are listed as now: `<file>: <line stripped>`, in full.
  The recap's display of the list is capped (§3.6: 5 lines, each cut to 160
  characters); `away.py --back` prints the list uncapped, as it is not
  hook output.
- `all_decisions(repo: Path) -> list[str]`: `pending_decisions` over
  `<c>/.claude/current-task` for each `c` in `checkouts(repo)`. Lines from the
  main checkout carry no prefix; lines from another checkout are prefixed with
  its path relative to the main checkout when it is inside it, else absolute,
  and `: ` (`.claude/worktrees/15c-1: session.md: ASK OLA: ...`).

### 3.3 The `session.md` format (item 4)

`session_format(text: str, name: str = "session.md") -> list[str]`, pure:

- Kind of a line: `^[ \t]*(?:[-*+][ \t]+)?(NOW|QUEUE|ASK OLA):(.*)$`,
  case-sensitive. Blank lines are ignored.
- Warnings, in this order and with these texts:
  - `<name>: <n> NOW lines; exactly one expected` when n is not 1; the same
    for `QUEUE`;
  - `<name>:<lineno>: empty NOW line` (or `QUEUE`); an empty `ASK OLA:` is
    already warned by §3.2 and is not repeated here;
  - `<name>:<lineno>: not a NOW, QUEUE or ASK OLA line: <first 60 characters>`;
  - `<name>:<lineno>: <n> characters; at most 300` for any non-blank line,
    of any kind, longer than 300 characters (`len()` of the line without its
    newline and trailing whitespace; `MAX_LINE = 300`). 300 passes, 301 warns.
- An absent `session.md` gives no warning (it is deleted when a round lands).

The recap runs this on the main checkout's `session.md` only, prints the
warnings under `session.md format:` directly after "Waiting on Ola", at most
4 then `... <n> more` (§3.6's budget), and prints nothing when there are none. A worktree's
`session.md` is not format-checked; its `ASK OLA:` lines are still listed
(§3.2). The 300-character limit is Ola's
ruling 2 (§7), made after a 971-character `NOW:` line.

### 3.4 Fallbacks and the longest stretch without a commit (item 1)

The rule is text only (§4). Item 17, if ruled, would extend `session_format`:
a `FALLBACK:` item inside the `QUEUE:` line names the paths it writes, and the
check flags any governed one, using `GOVERNED` from
`.claude/hooks/guard_governance.py` (or `tools/governed.py` once h5 lands).

In `away.py`:

- `longest_quiet(since: datetime, end: datetime, commits: list[tuple[datetime, str]]) -> Quiet`,
  pure. `Quiet` is a frozen dataclass `(start, end, opener: str | None,
  count: int)`: the longest interval between consecutive points of
  `since`, the commit times inside `[since, end]` sorted, and `end`. `opener`
  is the short hash of the commit that opens the interval, `None` when it
  opens at `since`. `count` is the number of commits inside the window.
  Commits outside the window are ignored; input order does not matter; on a
  tie the earliest interval wins.
- `commits_between(root, since, end) -> list[tuple[datetime, str]] | None`:
  `git -C root log --all --format='%h %cI'`, parsed, and filtered to
  `[since, end]` in Python; `None` on git failure. No `--since`/`--until`:
  git stops a date-limited walk early when commit dates are out of order,
  and the whole history (1,355 commits on 2026-10-03, `git rev-list --all
  --count`) is cheap to read. `--all` because work lands on worktree
  branches.
- In `back()`, after the `Back. Unattended since ...` line, when the flag was
  present: the window is `[since, min(now, until)]` (`until` = the flag's
  `until`, or `now` when unreadable), and it prints exactly

  ```
  Longest stretch without a commit (any ref): {h} h {mm:02d} min, {start:%Y-%m-%d %H:%M} to {end:%Y-%m-%d %H:%M} UTC, {opener}; {count} commit(s) in the window.
  ```

  `{h}`/`{mm}` are the whole minutes of the interval, floored; times in UTC;
  `{opener}` is `after <hash>` or `from the window's start`. When the flag's
  `since` cannot be read (a `{}` flag), it prints
  `Longest stretch without a commit: unknown (the flag has no readable since).`;
  when git fails, `... unknown (git log failed).` With no flag, nothing.

The figure is not written to `windows.jsonl` (§7, ruling 4).

### 3.5 Running background jobs (item 5)

A Bash-tool command in Claude Code runs as a shell whose command line sources
a snapshot under `~/.claude/shell-snapshots/` and `eval`s the command, and a
background one keeps running under the Claude Code daemon after a restart
(evidence §2). Observed on this Mac on 2026-10-03 with `ps -e -ww -o
pid=,ppid=,etime=,command=`:

```
/bin/zsh -c source /Users/.../.claude/shell-snapshots/snapshot-zsh-....sh 2>/dev/null || true && ... && eval '<the command>' < /dev/null && pwd -P >| /tmp/claude-...-cwd
```

The format is Claude Code's, undocumented, and may change; the parse fails
towards printing nothing, so the section's heading says what it looks for.

- `background_jobs(ps_text: str, exclude: set[int]) -> list[str]`, pure: each
  line whose command contains `/.claude/shell-snapshots/` and `eval '`, not in
  `exclude`, as `pid <pid>, running <etime>: <command>`, where `<command>` is
  the text between `eval '` and the next `' < /dev/null` (or the end of the
  line), whitespace collapsed, cut to 100 characters. At most 4, then
  `... <n> more` (§3.6's budget).
- The caller runs `ps` as above (both macOS and procps accept it) and passes
  the PIDs of its own ancestors (walked through the `ppid` column from
  `os.getpid()`), so a recap run by hand does not list its own shell.
- **Format self-check**, so drift in Claude Code's wrapper shows rather than
  reading as "nothing running". `wrapper_check(ancestors: list[str]) ->
  Literal["confirmed", "drift", "unknown"]`, pure, takes the command lines
  of the ancestors from the parent upward and stops at the first Claude Code
  process (first token's basename `claude` or `claude.exe`):
  - an ancestor recognised as a wrapper (the test `background_jobs` uses)
    before that: `confirmed`;
  - otherwise, an ancestor before it that is a shell run with `-c`
    (basename `sh`, `bash` or `zsh`, second token `-c`) and is not the hook
    launcher: `drift`. The hook launcher is a `-c` shell whose argument
    starts with `python` and names `tools/session_state.py`: the
    `SessionStart` command in `.claude/settings.json`;
  - otherwise (no Claude Code ancestor, as in a terminal or a test; or only
    the hook launcher): `unknown`.

  On `drift` the recap prints `(wrapper format not recognised)` in place of
  `(none)`, and after the job lines when there are any. The hook-launcher
  shape is an assumption: Claude Code's hook process tree is not documented
  and was not observed for this design. The first `SessionStart` after green
  is the probe: the main session checks that the recap does not print the
  drift line, and reports it if it does.
- The recap prints, after "In flight":
  `Running Bash-tool jobs, any session (Claude Code shell wrappers):` then the
  lines, `(none)`, `(wrapper format not recognised)`, or
  `(could not list: <error>)`.

Foreground calls of other sessions show up too. That is wanted: it is what
the "agents in pairs" note needs before a second C++ build.

### 3.6 The size table (item 7)

New `tools/rule_sizes.py`, governed (§5):

- `RULE_FILES = ("CLAUDE.md", ".claude/REQUIRED-READING.md",
  "docs/increments/README.md", "docs/PRINCIPLES.md")`, then
  `.claude/agents/*.md` and `.claude/skills/**/*.md`, each group sorted.
  `docs/PRINCIPLES.md` is the fourth governed rule file, in by §7, ruling 3.
- `words(text: str) -> int`: `len(text.split())`, which is what `wc -w`
  counts for ASCII text.
- **The reference is derived, not stored**: the newest commit reachable from
  `HEAD` that added a file matching `docs/retrospectives/2???-??-??-*.md`
  (`git log -1 --diff-filter=A --format=%h --name-only -- <pattern>`), and
  the counts are the files matching the same set at that commit, listed with
  `git ls-tree -r --name-only <rev> -- <paths>` (so a file removed since is
  found) and read with `git show <rev>:<path>` or one `git cat-file --batch`. @orchestrator writes a dated
  file at each retrospective, so the reference moves with no extra step and
  there is no reference file to keep, forget or reset (§7, ruling 1).
- `table(now: dict[str, int], then: dict[str, int] | None, label: str) -> list[str]`,
  pure. First line `Rule text in words, change since <label>:` where label is
  `<retrospective file name> (<hash>)`; then one line per file,
  `  <words right-aligned to 5> <change> <path>`, change as `+12`, `-3`, `0`,
  `new`, or for a file present only at the reference `removed (was <n>)`;
  then `  <total> <total change> total`. With no reference:
  `Rule text in words (no retrospective found):` and no change column.
- `main()` prints the table, so `python3 tools/rule_sizes.py` works on its own.

The recap prints the table last in `== recap ==`, after "Next on ROADMAP.md".
At today's 14 files it is about 16 lines and 600 characters.

**Budget against the 10,000-character cap.** The recap printed 6,025
characters on 2026-10-03 (round 1 review, `python3 tools/session_state.py |
wc -m` in the main checkout), leaving about 3,975. That headroom moves with
the length of the predecessor turns, so the four sections h8 adds or widens
get a fixed budget of **3,000 characters together**, newlines and headings
included, which leaves about 1,000 of margin at today's size. Each section's
cap, with its worst case:

| Section | Cap | Worst case, characters |
|---|---|---|
| Running Bash-tool jobs (§3.5) | 4 lines of at most ~136 (`  pid`, etime, 100-character command), then `... <n> more` | ~630 |
| Waiting on Ola (§3.2) | 5 lines, each cut to 160 characters for display, then `... <n> more` | ~850 |
| `session.md format:` (§3.3) | 4 lines of at most ~112, then `... <n> more` | ~490 |
| Size table (§3.6) | one line per file, `  <words:5> <change:>6> <path>`; 14 files today, longest path 48 characters, so lines of at most 63 | ~980 |
| **Total** | | **~2,950** |

The 160-character cut on an `ASK OLA:` line is a display cut in the recap
only, marked with `...`; `session.md` lines may be up to 300 characters
(ruling 2), and the main checkout's `session.md` is printed in full further
down by `print_current_task`, so a cut or dropped line there is still seen.
The size table grows by up to 63 characters per added rule file, so about
one more file uses up the budget's slack; test 11 then fails, and the fix is
to raise the budget against the measured recap or to lower a cap.

### 3.7 Recap order

Header (unchanged); `== recap ==`; Last landed; In flight; Running Bash-tool
jobs (§3.5); Waiting on Ola (§3.1-3.2); `session.md format:` (§3.3, only when
warnings); the harness queue (unchanged); Next on ROADMAP.md; the size table
(§3.6). Then the current-task print (from the main checkout) and the
predecessor turns.

## 4. Rule text, and what it replaces

All in `.claude/REQUIRED-READING.md`. The word counts are `len(text.split())`
of the passages at this design's base (595c56a) and of the new text below;
the current-task row counts the closing "The spawner" on both sides.

| Passage (lines at 595c56a) | Change | Removed | Added | Net |
|---|---|---|---|---|
| Cold start, steps 1-2 (12-20) | step 1 lists "running background jobs"; step 2 gains the item 5 sentence | 95 | 108 | +13 |
| `.claude/current-task/` bullet (30-37) | the three-line rule moves to the subagent files; `session.md` gets item 4's format and ruling 2's limit | 89 | 104 | +15 |
| new bullet after "One session per working tree" | item 6 | 0 | 27 | +27 |
| Unattended, last sentence (158-159) | item 1 | 19 | 48 | +29 |
| `SessionStart` paragraph, tail (168-175) | the cap sentence shortened; "thin recap in a worktree" replaced by what item 2 shares | 89 | 75 | -14 |
| **Total** | | 292 | 362 | **+70** |

`CLAUDE.md`, `docs/increments/README.md`, the persona and skill files:
unchanged (0 added, 0 removed); ruling 1 keeps `orchestrator.md` out of it.

The text, verbatim (Markdown as it goes in):

Steps 1-2:

```
1. The recap from `tools/session_state.py`: the round recap (last landed, in
   flight, running background jobs, decisions waiting on Ola, next ROADMAP
   items), `.claude/current-task/` and the predecessor's last human turns,
   **including prompts Ola queued and the harness absorbed mid-turn** (they
   never appear as normal turns). The `SessionStart` hook puts it in context
   on startup, resume, `/clear`, compaction and fork; if it is missing, run
   `python3 tools/session_state.py` yourself.
2. Judge each surfaced turn against the tree: a turn with no answering commit,
   file or PR is still pending, and an earlier pending turn outranks your
   reconstruction of "what comes next". Start no background job before the
   running ones are listed.
```

The current-task bullet (the rest of the subagent sub-bullet, from "The
spawner names the path", is unchanged):

```
- **The current ask lives in `.claude/current-task/`**, untracked and
  gitignored. **One file per writer; the spawner deletes it.**
  - `session.md` is the main session's, and only the main session writes it,
    in place; it is deleted when the round lands. Each line starts with
    `NOW:` (one), `QUEUE:` (one) or `ASK OLA:` (one per decision waiting on
    Ola), at most 300 characters each; no rulings, no history. The recap
    lists every worktree's `ASK OLA:` lines and warns on any other line.
  - Every other file is one subagent's, `<persona>-<HHMMSS>.md`, three lines
    or fewer: the ask, the persona and the file it will produce. The spawner
```

New bullet, after "One session per working tree":

```
- **After merging a change to `CLAUDE.md` or `.claude/agents/`, restart
  before spawning the changed persona**; until then its brief says to read
  the persona file from disk.
```

Unattended, replacing its last sentence:

```
When Ola says he is leaving, ask him how long, and ask him to run
`away.py` with that duration. Before he leaves, the `QUEUE:` line names at
least one fallback that needs no ruling and writes no governed path; idle
is accepted only when no such item exists.
```

`SessionStart` paragraph, from "spawns with the hook live" to its end:

```
spawns with the hook live for a recap it should not have. Claude Code caps
the hook's stdout at 10,000 characters (`python3 tools/session_state.py | wc
-m` measures it) and passes only a 2,000-character preview past the cap, so
keep `session.md` and the subagent files short. The recap finds the main
checkout from the repository's common git dir, so a session launched inside
a worktree gets the main checkout's `session.md` and every worktree's
`ASK OLA:` lines.
```

Item 2's sentence is true only once §3.1 is green; it lands in the green
commit's PR, not before. Two memory notes become redundant once this merges
and can be cut: `fill-the-unattended-window.md` (now a rule) and
`background-jobs-survive-restart.md` (now in the recap and step 2).

## 5. Files

| File | Change | Production lines (CLAUDE.md §2) |
|---|---|---|
| `tools/session_state.py` | §3.1, 3.2, 3.3 (with the line limit), 3.5, 3.7; docstring | ~98 |
| `tools/rule_sizes.py` | new, §3.6 | ~50 |
| `tools/away.py` | §3.4, decisions from every worktree | ~35 |
| `.claude/hooks/guard_governance.py` | `"tools/rule_sizes.py"` in `GOVERNED`: it runs inside the `SessionStart` hook, like the rest of the self-protecting set | 1 |
| `.claude/REQUIRED-READING.md` | §4 | prose |
| `tests/python/harness_fixtures.py` | `tools/rule_sizes.py` in `COPIED` | test |
| **Total** | | **~185**, under the 700 ceiling |

Overlap: h5 (red step on `worktree-h5-state-check`, not merged) also edits
`tools/away.py` (`--back` deletes `tampered.json`; `enter()` installs its git
shims), appends tests at the end of `test_away.py`, edits
`harness_fixtures.py`, and moves `GOVERNED` into `tools/governed.py`. The
conflicts are a few lines in `back()` and `enter()`, the appended tests and
one tuple; whichever merges second resolves them. Nothing in h8 depends on h5.

No suite here is invariant-critical (README, *Cost constraints*): no mutation
round. No refine or mesh code: no `@perf` acceptance.

## 6. Test plan for @tester (red, before any code)

`tests/python/test_session_state.py`, `tests/python/test_away.py`, and a new
`tests/python/test_rule_sizes.py`. Fixture repositories as
`harness_fixtures.make_repo`; worktrees with `git worktree add`; commit times
set with `GIT_COMMITTER_DATE`. Pure functions are tested on strings and lists;
one end-to-end run per section through `run_script`.

Amended (a ruling changed them, so the amendment says so in its commit):
`test_pending_decisions_collects_ask_ola_lines_from_every_task_file` (the
lowercase line no longer counts), and
`test_back_survives_a_broken_session_state` and
`test_back_prints_the_queue_by_branch_and_archives_it`
(`tests/python/test_away.py:396` asserts the old heading; new heading, §3.1).

1. **Matching (§3.2).** Counted: `ASK OLA: x`, `- ASK OLA: x`, `* ASK OLA:
   x`, `  + ASK OLA: x`. Not counted: `ask ola: x`, `Note: ASK OLA: x`,
   `**ASK OLA:** x`, `1. ASK OLA: x`, `ASK OLA x`. Empty: `ASK OLA:` and
   `- ASK OLA:   ` give `<file>:<lineno>: WARNING: empty ASK OLA line`.
2. **Checkouts (§3.1).** From the main checkout and from inside a worktree,
   `main_checkout` is the main path (compare `resolve()`d); `checkouts` lists
   main first, then both worktrees; a worktree whose directory was deleted is
   skipped; outside a repository both fall back to the given path.
3. **Every worktree's decisions.** Lines in main and two worktrees, listed
   with the §3.2 prefixes; the same list whichever checkout it starts from.
   End to end: the recap run from a worktree copy prints the main checkout's
   `session.md` and the worktrees' lines.
4. **Format (§3.3).** No warnings for NOW, QUEUE, two ASK OLA lines and a
   blank line, with and without bullets. Each warning text once: no NOW; two
   QUEUE; empty NOW; `Rulings: ...` on line 4 (the warning names line 4,
   §3.3's text); `now: x` (case) is "not a" line. Order: a text with every
   kind of fault gives the warnings in §3.3's order. Line length: a `NOW:`
   line of exactly 300 characters passes, a length over 300 warns with its
   line number and length, and a long `ASK OLA:` line warns too; trailing
   spaces do not count. Absent file: none. End to end: the recap prints the
   `session.md format:` block for a faulty main `session.md`, nothing for a
   good one, and nothing for a faulty `session.md` in a worktree (whose
   `ASK OLA:` lines are still listed).
5. **Background jobs (§3.5).** On a `ps` text built from the observed format:
   wrapper lines listed with pid, etime and the eval'd command; non-wrapper
   lines and excluded pids not; a command past 100 characters cut; 6
   wrappers give 4 lines and `... 2 more`. End to end: start
   `/bin/sh -c ": /x/.claude/shell-snapshots/s.sh && eval 'sleep 60' < /dev/null"`,
   run the recap, find its pid listed; kill it. With `ps` unavailable (a seam
   or an empty `PATH` for that call), the recap prints `(could not list: ...)`
   and the rest of the recap. `wrapper_check`, pure: a recognised wrapper
   below `claude` gives `confirmed`; a `/bin/zsh -c` ancestor of another
   shape below `claude` (drift) gives `drift`, and the recap then prints
   `(wrapper format not recognised)` instead of `(none)`; the hook launcher
   (`/bin/sh -c python3 "$CLAUDE_PROJECT_DIR/tools/session_state.py"`)
   below `claude` gives `unknown`; no `claude` ancestor gives `unknown`.
6. **Longest quiet (§3.4), pure.** No commits: the whole window, opener
   `None`, count 0. Commits at +1 h and +2 h in a 9 h window: 7 h after the
   second. Commits before `since` or after `end` ignored; unsorted input; a
   tie takes the earliest; the gap from `since` to the first commit wins when
   longest.
7. **`--back` prints it.** Flag window from 3 h ago to 1 h ahead, `now`
   passed in; commits on two branches at chosen times (`--all`); the printed
   line matches §3.4's template exactly, including `after <hash>`. The window
   ends at `now`, not the flag's `until`; with an expired flag, at `until`. A
   `{}` flag prints the `unknown` line; no flag, no line. With git failing
   (`commits_between` returning `None`, through a seam or a `root` that is
   not a repository), `--back` prints `... unknown (git log failed).` and
   still deletes the flag and archives the queue.
8. **`--back` from a worktree** lists the main checkout's `session.md`
   decisions: the reproduction of evidence §1b.
9. **Sizes (§3.6).** `words` equals `wc -w` on an ASCII fixture. In a
   fixture repo: rule files committed, a retrospective file added in a
   commit, then one file grown by 12 words, one cut by 3, one added and one
   deleted: the table shows `+12`, `-3`, `new`, `removed (was n)` and the
   total change. No retrospective: the no-reference heading. A later
   retrospective moves the reference. The recap prints the table, and with
   `tools/rule_sizes.py` broken it prints one `(size table unavailable: ...)`
   line and the rest.
10. **Governed.** A write to `tools/rule_sizes.py` gets `ask` from
    `guard_governance.py` (add it to the existing list in
    `test_guard_governance.py`).
11. **Budget (§3.6).** With 15 wrapper processes whose commands are 300
    characters long, 30 `ASK OLA:` lines of 300 characters across three
    worktrees, 15 format faults in the main `session.md`, and the size table
    over the real file set, each of the four sections prints its cap and
    `... <n> more`, and the four sections' blocks in the recap output (from
    each section's heading to the line before the next section's heading,
    newlines included) sum to at most 3,000 characters. It measures those
    blocks directly, not the recap with and without them, so the result does
    not depend on the fixture's transcripts or `ROADMAP.md`. An `ASK OLA:`
    line is shown cut to 160 characters ending in `...`.

## 7. Rulings (Ola, 2026-10-03)

The four questions of the first draft, answered:

1. **Where the reference word counts live:** derived from git, the commit
   that added the newest dated retrospective file (§3.6). Nothing stored.
2. **A length limit on `session.md`:** yes, 300 characters per line; the
   recap warns beyond it (§3.3, test 4). Ola changed his mind on seeing a
   971-character `NOW:` line.
3. **`docs/PRINCIPLES.md` in the size table:** yes (§3.6).
4. **The longest quiet stretch in `windows.jsonl`:** no; `--back` prints it
   only (§3.4).

## Review

**Design review, round 1, 2026-10-03.** Range `5e520fe..e5d308f` (merge ea813fc). Verdict: CHANGES REQUESTED. LOC: 0 (design only); the ~185 estimate is plausible. Rulings 1, 3, 5, 6 and Ola's four answers implemented as ruled; rule text word for word; code sites checked. Blocking: (1) word counts off by two (292 removed, 355 added; net +63 right); (2) reading predecessor transcripts from the main and running checkouts was not ruled, is untested, and would print the running main session's turns as pending; (3) format warnings for every checkout's session.md go beyond item 4 (about 25 warning lines from four stale files); (4) the background-job listing prints "(none)" if the wrapper format drifts: add a recognised-format self-check; (5) `test_away.py:396` also asserts the old heading; (6) "every file it changes is governed" should say production files; the 961 vs 971 figure. Not pushed; no CI.

**Design review, round 2, 2026-10-03.** Range `2d5ac21..7523c8d`. Verdict: CHANGES REQUESTED. LOC: 0 (design only). Round-1 findings fixed; word counts recounted: 292 removed, 362 added, net +70. Blocking: §3.6's budget does not fit its own caps (at the caps the new sections add about 6,000 characters to a 6,025-character recap, over the 10,000 limit); fit the caps to the headroom and say what test 11 measures against; state the 160-character ASK line cut as a display cut, not a second limit. Not pushed; no CI.

**Design review, round 3, 2026-10-03.** Range `352de1e..d636eea`. Verdict: APPROVED. LOC: 0 (design only); about 185 estimated. The budget fits its cap: recomputed worst cases 666 (jobs), 845 (Waiting on Ola), 481 (format warnings), 726 (size table, actual files), 2,718 in all (2,972 with the design's cautious table figure), under 3,000. The display cut, the full session.md print, test 5 and test 11 are consistent. Suggestions: say ~670 for the jobs row; `removed (was n)` lines run past 63 characters. Not pushed; no CI.
