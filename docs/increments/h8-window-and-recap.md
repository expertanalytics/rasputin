# Harness h8: the unattended window, the recap and the size table

Status: **design**, @architect, 2026-10-03. One PR, by day: every file it
changes is governed. Implements Ola's rulings of 2026-10-03 on items 1 to 6 of
`docs/retrospectives/next.md`, section "The window of 2026-10-02, the restart
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
| `session.md` format warnings | each checkout's `session.md` that exists (§3.3) |
| Predecessor turns | the transcript folders of both the main checkout and the running checkout (`~/.claude/projects/<slug>`, slug as now), merged by time; one folder when they are the same |

`TRANSCRIPTS` and the `.claude/current-task` path stop being import-time
constants derived from `__file__`; they are computed in `main()` so tests can
point them at a fixture repository.

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
  line. Others are listed as now: `<file>: <line stripped>`.
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
  - `<name>:<lineno>: not a NOW, QUEUE or ASK OLA line: <first 60 characters>`.
- An absent `session.md` gives no warning (it is deleted when a round lands).

The recap prints the warnings of every checkout's `session.md` under
`session.md format:`, prefixed as in §3.2, directly after "Waiting on Ola",
and prints nothing when there are none. There is no length check: the ruling
was the format (see §6, question 2).

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
  `git -C root log --all --format='%h %cI' --since=<since> --until=<end>`,
  parsed; `None` on git failure. `--all` because work lands on worktree
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

The figure is not written to `windows.jsonl` (§6, question 4).

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
  line), whitespace collapsed, cut to 100 characters. At most 10, then
  `... <n> more`.
- The caller runs `ps` as above (both macOS and procps accept it) and passes
  the PIDs of its own ancestors (walked through the `ppid` column from
  `os.getpid()`), so a recap run by hand does not list its own shell.
- The recap prints, after "In flight":
  `Running Bash-tool jobs, any session (Claude Code shell wrappers):` then the
  lines, `(none)`, or `(could not list: <error>)`.

Foreground calls of other sessions show up too. That is wanted: it is what
the "agents in pairs" note needs before a second C++ build.

### 3.6 The size table (item 7)

New `tools/rule_sizes.py`, governed (§5):

- `RULE_FILES = ("CLAUDE.md", ".claude/REQUIRED-READING.md",
  "docs/increments/README.md", "docs/PRINCIPLES.md")`, then
  `.claude/agents/*.md` and `.claude/skills/**/*.md`, each group sorted.
  `docs/PRINCIPLES.md` is the fourth governed rule file; Ola named three
  (§6, question 3).
- `words(text: str) -> int`: `len(text.split())`, which is what `wc -w`
  counts for ASCII text.
- **The reference is derived, not stored**: the newest commit reachable from
  `HEAD` that added a file matching `docs/retrospectives/2???-??-??-*.md`
  (`git log -1 --diff-filter=A --format=%h --name-only -- <pattern>`), and
  the counts are the same files read at that commit (`git show
  <rev>:<path>`, or one `git cat-file --batch`). @orchestrator writes a dated
  file at each retrospective, so the reference moves with no extra step and
  there is no reference file to keep, forget or reset. A stored JSON is the
  alternative (§6, question 1).
- `table(now: dict[str, int], then: dict[str, int] | None, label: str) -> list[str]`,
  pure. First line `Rule text in words, change since <label>:` where label is
  `<retrospective file name> (<hash>)`; then one line per file,
  `  <words right-aligned to 5> <change> <path>`, change as `+12`, `-3`, `0`,
  `new`, or for a file present only at the reference `removed (was <n>)`;
  then `  <total> <total change> total`. With no reference:
  `Rule text in words (no retrospective found):` and no change column.
- `main()` prints the table, so `python3 tools/rule_sizes.py` works on its own.

The recap prints the table last in `== recap ==`, after "Next on ROADMAP.md".
At today's 14 files it is about 16 lines and 600 characters, against the
10,000-character cap on the hook's output.

### 3.7 Recap order

Header (unchanged); `== recap ==`; Last landed; In flight; Running Bash-tool
jobs (§3.5); Waiting on Ola (§3.1-3.2); `session.md format:` (§3.3, only when
warnings); the harness queue (unchanged); Next on ROADMAP.md; the size table
(§3.6). Then the current-task print (from the main checkout) and the
predecessor turns.

## 4. Rule text, and what it replaces

All in `.claude/REQUIRED-READING.md`. The word counts are `len(text.split())`
of the passages at this design's base (595c56a) and of the new text below.

| Passage (lines at 595c56a) | Change | Removed | Added | Net |
|---|---|---|---|---|
| Cold start, steps 1-2 (12-20) | step 1 lists "running background jobs"; step 2 gains the item 5 sentence | 95 | 108 | +13 |
| `.claude/current-task/` bullet (30-37) | the three-line rule moves to the subagent files; `session.md` gets item 4's format | 87 | 97 | +10 |
| new bullet after "One session per working tree" | item 6 | 0 | 27 | +27 |
| Unattended, last sentence (158-159) | item 1 | 19 | 48 | +29 |
| `SessionStart` paragraph, tail (168-175) | the cap sentence shortened; "thin recap in a worktree" replaced by item 2 | 89 | 68 | -21 |
| **Total** | | 290 | 348 | **+58** |

`CLAUDE.md`, `docs/increments/README.md`, the persona and skill files:
unchanged (0 added, 0 removed). With question 1 answered "stored JSON", add
about 25 words to `.claude/agents/orchestrator.md`.

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
    Ola); no rulings, no history. The recap lists every worktree's `ASK OLA:`
    lines and warns on any other line.
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
a worktree gets the same recap.
```

Item 2's sentence is true only once §3.1 is green; it lands in the green
commit's PR, not before. Two memory notes become redundant once this merges
and can be cut: "Fill the unattended window" (now a rule) and the "check
`ps` before resuming" note (now in the recap).

## 5. Files

| File | Change | Production lines (CLAUDE.md §2) |
|---|---|---|
| `tools/session_state.py` | §3.1, 3.2, 3.3, 3.5, 3.7; docstring | ~95 |
| `tools/rule_sizes.py` | new, §3.6 | ~50 |
| `tools/away.py` | §3.4, decisions from every worktree | ~35 |
| `.claude/hooks/guard_governance.py` | `"tools/rule_sizes.py"` in `GOVERNED`: it runs inside the `SessionStart` hook, like the rest of the self-protecting set | 1 |
| `.claude/REQUIRED-READING.md` | §4 | prose |
| `tests/python/harness_fixtures.py` | `tools/rule_sizes.py` in `COPIED` | test |
| **Total** | | **~180**, under the 700 ceiling |

Overlap: h5 (red step on `worktree-h5-state-check`, not merged) also edits
`tools/away.py` (`--back` deletes `tampered.json`), `test_away.py` and
`harness_fixtures.py`, and moves `GOVERNED` into `tools/governed.py`. The
conflicts are a few lines in `back()` and one tuple; whichever merges second
resolves them. Nothing in h8 depends on h5.

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
lowercase line no longer counts) and `test_back_survives_a_broken_session_state`
(new heading, §3.1).

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
   QUEUE; empty NOW; `Rulings: ...` on line 4 (the warning names line 4, §3.3's text);
   `now: x` (case) is "not a" line. Absent file: none. The recap prints the
   `session.md format:` block only when there are warnings.
5. **Background jobs (§3.5).** On a `ps` text built from the observed format:
   wrapper lines listed with pid, etime and the eval'd command; non-wrapper
   lines and excluded pids not; a command past 100 characters cut; 12
   wrappers give 10 lines and `... 2 more`. End to end: start
   `/bin/sh -c ": /x/.claude/shell-snapshots/s.sh && eval 'sleep 60' < /dev/null"`,
   run the recap, find its pid listed; kill it. With `ps` unavailable (a seam
   or an empty `PATH` for that call), the recap prints `(could not list: ...)`
   and the rest of the recap.
6. **Longest quiet (§3.4), pure.** No commits: the whole window, opener
   `None`, count 0. Commits at +1 h and +2 h in a 9 h window: 7 h after the
   second. Commits before `since` or after `end` ignored; unsorted input; a
   tie takes the earliest; the gap from `since` to the first commit wins when
   longest.
7. **`--back` prints it.** Flag window from 3 h ago to 1 h ahead, `now`
   passed in; commits on two branches at chosen times (`--all`); the printed
   line matches §3.4's template exactly, including `after <hash>`. The window
   ends at `now`, not the flag's `until`; with an expired flag, at `until`. A
   `{}` flag prints the `unknown` line; no flag, no line.
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

## 7. Questions for Ola (each with a default)

1. **Where the reference counts live.** Default: derived from git, the
   commit that added the newest dated retrospective file (§3.6); nothing to
   store or update. Alternative: a committed JSON that @orchestrator updates
   at each retrospective (about 25 more words in `orchestrator.md`, and a
   reference that can be reset by hand).
2. **A size bound on `session.md`.** The ruling enforces the kinds of line,
   not their length. Measured 2026-10-03: the main checkout's `session.md`
   is 10 lines and 2,468 characters (`wc -lc`), nearly all of it in the
   `NOW:` and `QUEUE:` lines (961 and 943 characters), and it passes §3.3
   as designed (its two blank lines are ignored). Default: no bound in h8; if wanted, a warning past
   1,500 characters is three lines.
3. **`docs/PRINCIPLES.md` in the size table** (1,421 words today). Default:
   in, as the fourth rule file `guard_governance.py` governs.
4. **Record the longest quiet stretch in `windows.jsonl`** as well as
   printing it, so the nightly trials have a series. Default: no; the print
   is what was ruled, and one field can follow.
