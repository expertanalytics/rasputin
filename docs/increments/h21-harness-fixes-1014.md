# Harness h21: nine fixes from the 2026-10-14 morning check

Status: **design**, 2026-10-14, on `worktree-harness-1014` off master
`55b9f564`. Two PRs, by day (section 11): every production file is governed.
Implements nine changes Ola approved on 2026-10-14 ("yes to all, push and
enqueue"): P1 to P5 of `docs/retrospectives/2026-10-14-night-snow-46.md`
section 6; the three guard false positives that retrospective lists as
waiting on Ola in its section 2 (the GraphQL read, the scratchpad clone, the
after-commit hook on a conflicted merge); and P2 and P4 of
`docs/retrospectives/2026-10-13-merge-43.md` section 7 (P2: a read-only
`gh api graphql` query is not a forge write; P4: the session recap lists
every open PR with its age and state). `session.md`'s `ASK OLA` lines, and
the 2026-10-14 retrospective's section 2, call those two "Q1" and "Q2".

Acceptance (PR A, prose): `@reviewer` checks each rule's incident commit; no run beyond the governance gates.
Acceptance (PR B, code): the harness suites (`pytest -m harness`) and the three live checks of section 9, recorded in this file; no `@perf` run (no refine or mesh code).

Words used here. **Pre hook / post hook**: a Claude Code hook that runs
before (`PreToolUse`) or after (`PostToolUse`) a tool call. **Scratchpad**: a
session's temporary directory, `/private/tmp/claude-<uid>/<project>/<session>/scratchpad`
(`tools/scratchpad.py`). **Governed**: a file `guard_governance.py` asks
about; a change to it needs Ola's prompt, so it is made by day.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature.*

- GraphQL, *Specification*, October 2021 edition, section 2.3 *Operations*
  (spec.graphql.org/October2021, read 2026-10-14): "query – a read-only
  fetch", "mutation – a write followed by a fetch", "subscription – a
  long-lived request ..."; the *Query shorthand* omits "the query keyword and
  operation name" and is a query. An operation type is one token, so a
  mutation's document holds the word `mutation` as written (GraphQL names
  are ASCII, `[_A-Za-z][_0-9A-Za-z]*`, and no escape forms a keyword). Item 1
  reads an operation as a read on that basis. Departure: it does not parse
  GraphQL; it refuses a line that contains the word `mutation` or
  `subscription` anywhere, comments and strings included, which
  over-refuses (a read whose text names the word). Its soundness rests on
  two premises, each enforced or stated in section 4.1: the text the guard
  reads is the text `gh` sends (no shell expansion; `-f` sends "static
  string parameters", where `-F` fills `{owner}`, `{repo}` and `{branch}` and
  reads `@file`, per `gh api --help`, gh 2.101.0), and nothing else on the
  line changes what `gh` runs.
- Claude Code hooks reference (code.claude.com/docs/en/hooks, read
  2026-10-14): every event's `cwd` is the "Current working directory when the
  hook is invoked", and "the new directory after Claude runs `cd`"; both tool
  events carry `tool_use_id`; `PostToolUse` fires after a tool call
  succeeds, a failed one goes to `PostToolUseFailure`; Pre events carry
  `scratchpad_dir`. Item 4 rests on the first two: a post hook's `cwd` is the
  directory *after* the command, so a line that ends in a `cd` moves the
  directory the after-commit hook resolves from (section 5).
- git 2.55: a fast-forward writes a HEAD reflog entry whose subject ends
  `: Fast-forward` (`merge side: Fast-forward`, `pull -q --ff-only:
  Fast-forward`); "Already up to date" writes none; a merge stopped on a
  conflict leaves `MERGE_HEAD` and HEAD unchanged; a real merge writes
  `merge main: Merge made by the 'ort' strategy.`. Measured by a probe script
  in the session scratchpad on a throwaway repository, 2026-10-14. Not
  probed, taken from git's documentation of `git log -g` and pinned by
  tests 12 and 17: a reflog entry's `%H` is the commit HEAD moved to, and
  `rev-parse --git-common-dir` in a worktree or behind a `.git` gitfile names
  the other repository's directory. A second probe was refused by the push
  guard (`git reset --hard` in a scratchpad repository, the false positive
  item 5 removes).
- Merge-queue practice (GitHub's merge queue, already in use here per
  `CLAUDE.md` §4) tests each PR on top of those ahead of it; item 6 does the
  same check by hand before the push question, with `git merge-tree`.

*Legacy.* Nothing: the legacy tree has no harness.

```
$ git grep -l -i -E 'graphql|merge-tree|reflog|MERGE_HEAD|gh pr list' legacy-archive -- legacy | wc -l
       0
```

## 2. The ideas file

Taken up or overlapped:

- *Recap: items in flight and open ideas*: item 2 takes the open-PR list
  with ages; the "open ideas: N" count and "in flight: N of 3" are left out.
- *After-commit hook: keep the cause when two faults meet*: overlaps item 4
  (both change `committed_trees`); left out, its case (`shell_scan` missing
  and the directory in no work tree) is unchanged here.

Left out, handed back unchanged: *Lead time ... in every recap*, *Monitor in
every persona's tool list*, *Governed-file asks answered without a prompt in
auto mode*, *Cap on runner and shell words per line in the push guard*.

## 3. Scope

| # | What | Production files | Rule files |
|---|---|---|---|
| 1 | `gh api graphql` with a read-only query is not a forge write | `.claude/hooks/guard_push.py` | none |
| 2 | The recap lists every open PR with its age | `tools/session_state.py` | none |
| 3 | The after-commit hook says "new commit", not "amend" | `.claude/hooks/gates_after_commit.py` | none |
| 4 | The after-commit hook is quiet when no commit was made, and resolves directories from where the line started | `.claude/hooks/record_bash_start.py` (new), `.claude/hooks/gates_after_commit.py`, `tools/shell_scan.py`, `.claude/settings.json` | `.claude/REQUIRED-READING.md` |
| 5 | git's local writes, and file writes, in a scratchpad repository pass the guards | `.claude/hooks/guard_push.py`, `.claude/hooks/guard_governance.py`, `tools/shell_scan.py` | `.claude/REQUIRED-READING.md` |
| 6 | Before a push question, check conflicts and merge master in | none | `CLAUDE.md` §3 |
| 7 | Briefs quote the increment's own `Acceptance` lines | `tools/brief.py` | `docs/increments/README.md` step 1 |
| 8 | Briefs name one log folder per worktree | `tools/brief.py` | `.claude/briefs/common.md` |
| 9 | `reviewer.md` points to step 4 instead of restating it | none | `.claude/agents/reviewer.md`, `docs/increments/README.md` step 4 |

Every file in both columns is governed (`guard_governance.py`'s `GOVERNED`,
`GOVERNED_PREFIXES` and `GOVERNED_GLOBS`), so each write needs Ola's prompt,
and none can be made while unattended mode is on. `.claude/settings.json`
(item 4) is also a settings change: Ola's own yes, by day.

## 4. Items 1 and 5: the push guard

### 4.1 A GraphQL read (item 1)

`segment_why` today calls any `gh api` with a field and no method a write,
because `gh` then POSTs, and every GraphQL call is a POST. New pure function
in `guard_push.py`, `graphql_read(rest: list[str], line: str) -> bool`, called
where the `gh api` write is decided, with the simple command's argv after
`api` and the whole command line as written; `True` means "not a write". It
is called only on the `shell_scan` path: a line `shell_scan.parse` cannot read
is judged by today's text rules, which ask. It is `True` only when all hold:

1. The endpoint (the first word after `api` that is not an option or an
   option's value) is exactly `graphql`.
2. Every option is one of: a field (`-f`, `-F`, `--field`, `--raw-field`,
   joined or separate, `--raw-field=k=v` included), `-H`/`--header`,
   `-q`/`--jq`, `-t`/`--template`, `--cache`, `--hostname`, `-p`/`--preview`
   (each with its value), `--paginate`, `--slurp`, `-i`/`--include`,
   `--silent`, `--verbose`; or `-X`/`--method` with `POST` or `GET`. Anything
   else (`--input` included) is `False`.
3. Exactly one field has the key `query`, and it is given by `-f` or
   `--raw-field` (sent as written). A `query` given by `-F`/`--field` is
   `False`: `gh` fills `{owner}`, `{repo}` and `{branch}` in it and reads
   `@file`, so the guard does not see the text sent.
4. That value, stripped of leading whitespace, starts with `{`, `query` or
   `fragment`.
5. `literal(line)`: the line, as written, has no shell expansion outside
   single quotes. A new pure function in `tools/shell_scan.py`, a left-to-right
   scan with three states (outside quotes, in `'...'`, in `"..."`; a
   backslash outside single quotes takes the next character): `False` on `$`
   or a backtick outside single quotes, and outside any quotes on `{`, `}`,
   `*`, `?`, `[`, `~` or a backslash (brace, glob and tilde expansion, and
   escapes such as `mu\tation`); also `False` on an unclosed quote. So
   `$M`, `"$M"`, `$'\x6d'`, `$(printf mut)`, `` `…` `` and `mu{tation,}`
   are refused, and a query written in single quotes passes.
6. The whole line, as written, contains neither `mutation` nor
   `subscription` as a word (`\b(mutation|subscription)\b`): every field,
   every other command on the line, comments and strings included.

Other fields (GraphQL variables such as `-F number=5`, `-f owner=x`) are
allowed: a variable is a value, never an operation type.

What it does not claim. Given 5, the query text `gh` sends is the text the
guard reads, and given 6 a mutation's keyword is on the line, so a mutation
or subscription sent through `gh api graphql` by a line that passes 1 to 6
is not known to exist. It does not hold against a `gh` that is not GitHub's
(a wrapper earlier on `PATH`, `GH_*` settings pointing elsewhere); the
guards are tripwires (`.claude/REQUIRED-READING.md`, *The harness*). Rule 5
over-refuses any line with a `$` outside single quotes (`cd "$W" && gh api
graphql …` asks); write such reads as one plain line.

### 4.2 Git in a scratchpad repository (item 5, push guard)

`guard_push.py` today asks for git's *local* writes (`config`, `remote`,
`update-ref`, `symbolic-ref`, `replace`, a fetch into a named ref, history
rewrites, `--no-verify`) wherever they run. `git clone` and `git commit` pass
already (probe of 2026-10-14 against `guard_push.py@55b9f564`; pinned below as
regression tests). New: such a write passes when its repository is a
scratchpad's own, by all of the checks below; a failed check means today's
ask. Never passes, whatever the directory: `git push` (`PUSH`, `FORCE`),
`gh` of any kind, `curl`.

1. The line is one `shell_scan` reads, and the git command's directory is a
   `Path` from `directories(..., chain_only=True)` (below) under a
   scratchpad (`scratchpad.under`, by real path).
2. Its repository is the scratchpad's: `git -C <dir> rev-parse
   --path-format=absolute --git-common-dir`, run by the guard (read-only, 5 s
   timeout, assuming a local rev-parse of milliseconds, one per git write
   on the line; `GIT_*` removed from its environment), exits 0 and prints a path
   whose real path is under a scratchpad. A worktree or a `.git` gitfile in
   the scratchpad whose common directory is the real repository's (shared
   refs and config) fails here, as does a directory that does not exist yet
   when the guard runs (`git clone <url> <pad>/c && git -C <pad>/c config …`
   on one line asks; run the clone first).
3. The git command names no file outside its repository: `config` with
   `--global`, `--system` or `--file`/`-f` (joined or separate) fails, as
   do `--git-dir`/`--work-tree` and `-c include.path=…` anywhere in its argv.
4. The line's text holds no `GIT_[A-Z_]*` word (`GIT_DIR`, `GIT_CONFIG_GLOBAL`,
   `GIT_CONFIG_PARAMETERS`, …): the environment could move any of 1 to 3.
5. Nothing else on the line can change which repository the directory is
   between the guard's check and the write: every other simple command is
   `cd`/`pushd`/`popd`, or a git command whose subcommand is not `init`,
   `clone`, `worktree` or `submodule`, or a command with no `writes`, not
   `unknown`, and no interpreter `program` or `script`.

The directory of each git command comes from new pure functions in
`tools/shell_scan.py`, moved out of `gates_after_commit.py` so both hooks and
the new recorder share one reading:

```python
def directories(simples: list[Simple], start: Path, *, chain_only: bool) -> list[Path | str]:
    """The directory each simple command runs in: `cd`/`pushd`/`popd` followed
    from `start`; a str (the word as written) once a directory is unknown.
    With chain_only, a command's directory is a Path only while every
    top-level command up to it is `chained`, and a `nested` command's is
    always a str; the guards use this, the after-commit hook does not."""

def git_directory(argv: list[str], shell_dir: Path | str) -> Path | str:
    """The directory git works in: shell_dir, moved by `-C`; a str when
    `--git-dir` or `--work-tree` is given (h19 §3.2's rule, unchanged)."""

def literal(line: str) -> bool:
    """True when the line, as written, has no shell expansion outside single
    quotes (section 4.1, rule 5)."""
```

`Simple` gains two fields the scanner sets (item 5, review point 3):

- `nested: bool`: it runs inside `( … )`, `$( … )`, backticks, `<( … )`, or
  a nested parse (`sh -c`, a shell's heredoc, `eval`), so its `cd` moves no
  top-level command.
- `chained: bool`: it is a top-level command, no reserved word (`!`, `if`,
  `{`, `while`, …) was stripped from its front, and the only token between it
  and the previous top-level command is one `&&` (the line's first command:
  none). `;`, a newline, `||`, `|`, `&`, `(`, `)` or a reserved word there
  makes it `False`.

So with `chain_only`, `cd <pad> && git config a.b c` passes; `(cd <pad>) ;
…`, `false && cd <pad>; …`, `cd <pad> || true && …`, `! cd <pad> && …` and
`cd <pad> | …` give every later command a `str` directory. A `cd` that fails
inside an `&&` chain stops the line, so no later write runs in the wrong
place. `gates_after_commit.git_dirs` keeps its signature (its tests load it
from the hook) and becomes a composition of `directories(...,
chain_only=False)` and `git_directory`, following `cd` in a straight line as
today: there a wrong guess gates a wrong tree, it does not pass a write.
`start` is the Pre event's `cwd` (before the command, section 5).

## 5. Items 3 and 4: the after-commit hook

### 5.1 What goes wrong (item 4)

Three cases, all seen 2026-10-13 and 2026-10-14 (transcript of session
`63d5ece9`, and the `ASK OLA` lines in `session.md`):

- A `git merge` that stopped on a conflict, or found nothing to merge, made
  no commit, but the hook gated the tree and said "The commit is in".
- A `git merge --ff-only origin/master` fast-forward in the main checkout
  gated master's content as if it were this session's commit.
- `git merge -q --ff-only origin/master && ... && cd ../rasputin_scratch/ideas
  && ...` printed `NOT CHECKED: git ran in .../rasputin_scratch/ideas`. The
  merge ran in the main checkout; the hook resolved it from the post event's
  `cwd`, which is the directory after the line's `cd` (section 1). The doubled
  paths `.../worktrees/snow-anim/.claude/worktrees/snow-anim` in the same
  transcript are the same fault. So a line that commits and then `cd`s gates
  the wrong tree or none.

### 5.2 The start record: `record_bash_start.py` (new pre hook)

A pre hook on `Bash`, after the guards. When the command passes
`gates_after_commit`'s trigger (the same `COMMITTING` set; the text test
`TEXT_TEST` when `shell_scan` cannot read the line), it writes
`<git-common-dir>/harness/bash-start/<tool_use_id>.json`:

```json
{"cwd": "<event cwd>", "heads": {"<directory as resolved>": "<HEAD sha>" | null}}
```

one entry per committing git command's directory (section 4.2's functions,
from the event `cwd`), `null` where the directory is not in a work tree or
does not exist yet. It also deletes files in `bash-start/` older than 24 h.
It prints nothing and exits 0 in every case, a failure included: a missing
record only returns the post hook to today's behaviour. `tool_use_id`
missing or not `[A-Za-z0-9_]+`: no record.

The common dir is the hook's own repository's (`HOOK_ROOT`), as for the
rest of the harness state. Constant: 24 h for the sweep; assumes a call's
post hook runs within a day of its pre hook; a call whose post event never
comes (`PostToolUseFailure`, an interrupt) leaves one file of under 1 KB
until the sweep.

### 5.3 The post hook's rule

`gates_after_commit.py` reads and deletes the record for the event's
`tool_use_id`. **Without a record, it does exactly what it does today**
(every existing test in `tests/python/test_gates_after_commit.py` runs
without one and stays as it is). A record it cannot read counts as none: a
file that is not JSON, not an object, or whose `cwd` is not an absolute
path or `heads` not an object of 40-hex strings or `null`; and a directory
missing from `heads` is judged by today's rule (gated). With a record:

1. Directories are resolved from the record's `cwd`, not the event's.
2. For each resolved directory, in order, its work-tree top and HEAD now:
   - the directory exists and is in no work tree, now and in the record
     (`null` both times): skip, silently (git refused there; no commit);
   - HEAD equals the recorded HEAD and no operation is in progress: skip,
     silently (up to date, nothing to commit, a failed commit);
   - an operation is in progress (`MERGE_HEAD`, `CHERRY_PICK_HEAD`,
     `REVERT_HEAD`, or the `rebase-merge` / `rebase-apply` git path exists):
     if HEAD is unchanged, skip silently; if HEAD moved (a commit, then a
     merge that stopped), one line: `NOT CHECKED: <top> has a <merge|...> in
     progress; run the gates after the commit that ends it.`;
   - HEAD moved and **every** HEAD reflog entry since the recorded HEAD is a
     fast-forward: skip, silently. Read with `git log -g -n 50
     --format=%H%x09%gs HEAD`, newest first; each entry's `%H` is the commit
     HEAD moved to, so the entries since the recorded HEAD are those before
     the first whose `%H` is the recorded sha. Skip only when that sha is
     found, at least one entry precedes it, and each such entry's subject
     ends with `: Fast-forward`. The sha not among the 50, git failing, or
     one entry of another kind (`commit: …`, `merge …: Merge made by …`,
     `reset: …`): gate it. So `git commit -m x && git merge --ff-only y` is
     gated, and a fast-forward alone is not. Constant: 50 entries, assuming
     one Bash line moves HEAD fewer than 50 times; more than 50 gates, the
     safe side;
   - otherwise gate it, as today.
3. A `str` directory (unknown) stays `NOT CHECKED`, as today.

Nothing printed means exit 0.

### 5.4 The wording (item 3)

The opening sentence becomes: "The commit is in, and below are the gates
that are red or the trees this hook did not check. Fix it in a new commit
(an amend is a history rewrite, which the push guard asks about) -- do not
push or open a PR on this." The module docstring's "A commit is cheap to
amend" becomes "A fix is one more commit".

## 6. Item 2: open PRs in the recap

`tools/session_state.py` gains `open_prs(now) -> list[str]`, printed under
`In flight:` as `Open PRs (N):`, one line per PR, oldest first:

```
  #251  12h  CONFLICTING  worktree-example  An example title
  #255  2d 3h  worktree-catalogue-data  Increment 44: land-cover catalogue as data
```

From one call, `gh pr list --state open --limit 100 --json
number,title,headRefName,createdAt,mergeable`, run with a 10 s timeout from
the main checkout. The state word after the age: `CONFLICTING` when
`mergeable` is `CONFLICTING`; `mergeable unknown` when it is `UNKNOWN`;
nothing when `MERGEABLE`; any other value is printed as given. `UNKNOWN` is
common, not an error: GitHub computes the value on demand, and the first
call of 2026-10-14 (08:00 UTC) gave `UNKNOWN` for all three open PRs, the
next, seconds later, `MERGEABLE` for #258 and `CONFLICTING` for #251 and
#252. The recap does not wait or retry for it.

Left out of P4's list, with the reason: *reviewed* and *approved* (GitHub's
`reviewDecision` is `""` on every open PR and on the five newest merged
ones, 2026-10-14: reviews here are `## Review` lines in increment files, not
GitHub reviews), and *not enqueued* (`gh pr list --json` has no merge-queue
field in gh 2.101.0; it needs a second, GraphQL, call at every session
start). Age is from `createdAt` (opened, not last updated): under 1 h
`Nm`, under a day `Nh`, else `Nd Nh`. The title is cut at 50 characters with
`...`. No PRs: `Open PRs: none`. `gh` missing, exit not 0, timeout, or JSON
that is not a list of such objects: one line `Open PRs: unavailable
(<reason, at most 80 characters>)`, and the recap goes on.

Constants: 10 s timeout, assumes one GitHub API round trip (typically under
2 s); 100 PRs, about 9,000 characters at most, which with the rest could
pass the recap's 10,000-character cap, so the list stops at 30 lines with
`... N more (gh pr list)`. Every open PR is listed while there are 30 or
fewer; `gh pr list --state open` gave 3 on 2026-10-14 at 07:40 UTC.

## 7. Items 7 and 8: briefs

### 7.1 The acceptance lines (item 7)

`docs/increments/README.md` step 1 gains the rule (section 8); `brief.py`'s
`_from_increment` gains, for `@developer` and `@reviewer` only: every line
outside `## Review` that matches `^Acceptance( \([^)]*\))?:`, verbatim, cut
at 300 characters with `...`, printed right after the Status line and
outside the 8-line `MAX_QUOTED` cap. None found: `Acceptance: none stated in
<file>; read its sections on acceptance before you start.` The existing
`QUOTED` lines are unchanged.

Why a line and not the regex: increment 46's acceptance for PR 3 sat at line
780, past eight earlier `QUOTED` hits, and the one hit nearer the top (line
779) reads "the acceptance rule does not apply"; the brief for PR 3's review
said "no @perf" (retrospective section 3, item 1).

### 7.2 The log folder (item 8)

`common.md` gains one line, with a new template variable `$logs`:

"Write logs (suite, probe and build output) only under `$logs`; create it
first, and never write into another worktree's."

`brief.py` computes `$logs` as `<scratchpad>/logs/<worktree directory
name>`, where `<scratchpad>` is
`/tmp/claude-<uid>/<slug>/<CLAUDE_CODE_SESSION_ID>/scratchpad` resolved
(`<slug>`: the main checkout's path with every character outside
`[A-Za-z0-9]` replaced by `-`, as `session_state.SLUG`), when that directory
exists; otherwise `<tempfile.gettempdir()>/rasputin-logs/<worktree
directory name>`. The root `/tmp` is one module literal, `SCRATCH_ROOT`, so
the tests rewrite it as `point_scratchpad` does. One folder per worktree is
enough: `brief.py` already refuses two agents in one worktree with a writer
among them.

## 8. Rule text (lands with Ola's prompts, by day)

No dates, incidents or quotations: rules and checks only.

**Item 6, `CLAUDE.md` §3, a new bullet after "No asking between internal
steps, but stop at the remote"** (about +55 words):

> * **Conflicts before a push question.** Before asking Ola to push a branch,
>   run `git fetch` and `git merge-tree --write-tree origin/master <branch>`,
>   and the same against each branch queued or asked for ahead of it. When
>   master conflicts, a writer merges it in first (`tools/merge_master.py`);
>   when a branch ahead conflicts, the question names the files.

**Item 7, `docs/increments/README.md` step 1, after "LOC estimate. No
production code."** (about +40 words):

> The file states its acceptance in lines that start `Acceptance:` (or
> `Acceptance (PR n):` when it splits), each naming the run beyond the suites
> and its owner, or "none beyond the suites"; `tools/brief.py` quotes them to
> `@developer` and `@reviewer`.

**Item 9, `.claude/agents/reviewer.md` "Feedback format" item 1** becomes
(about -35 words):

> 1. **Verdict:** one of step 4's three (`docs/increments/README.md`), and
>    the handback ends with step 4's `Review:` line.

and `docs/increments/README.md` step 4's "`@reviewer`'s handback ends with
one line, `Review: <verdict> <sha>`" becomes "ends with exactly one line
`Review: <verdict> <sha>`, and no other line starts with `Review:`" (+8
words; the rule now lives only there).

**Items 4 and 5, `.claude/REQUIRED-READING.md` *The harness*:**

- "`gates_after_commit.py` puts the gates' own output in the transcript
  after a commit or merge" becomes "... after a git command that made a
  commit (not a fast-forward, nor a merge that stopped); `record_bash_start.py`
  notes where the line started, so a `cd` later on the line does not move
  it".
- "except an absolute path in a session scratchpad, written by Edit or Write
  or by a shell line that is one plain command (no `&&`, `;`, pipe,
  substitution or interpreter program). `guard_push.py` passes no git write
  in a scratchpad repository." becomes "except a path in a session
  scratchpad: written by Edit or Write, or by a shell line the guard can read
  whose directories are literal and reached only through `cd`s joined by
  `&&`, with no interpreter program and no link made on the line.
  `guard_push.py` passes git's local writes in a repository that is the
  scratchpad's own (not a worktree or gitfile of another), never a push; a
  `gh api graphql` query, written in single quotes and sent by `-f`, is a
  read."

**Item 8, `.claude/briefs/common.md`:** the line of section 7.2.

Net, by `tools/rule_sizes.py`: about +160 words (+30 at design review
round 1, for the `&&` chain, the scratchpad's own repository and the
GraphQL read).

## 9. Item 5, the governance guard, and the live checks

### 9.1 Scratchpad writes on a compound line (item 5)

`guard_governance.judge_bash` today exempts a scratchpad target only on a
line that is one plain command. New: on a line `shell_scan` can read, a
target is exempt when it, joined to its simple command's directory
(`shell_scan.directories(..., chain_only=True)` from the Pre event's `cwd`,
section 4.2: only through `cd`s joined by `&&` at top level), lies under a
scratchpad, provided no simple command on the line:

- is `unknown`, or has an interpreter `program` or `script`;
- makes a link: `ln`, `link`, or `cp` with `-s`, `-l`, `--symbolic-link`,
  `--link`, `-R`, `-r`, `-a`, `-P` (a recursive copy copies links as links on
  macOS).

Its directory must be a `Path` (a `str` directory exempts nothing). A line
`shell_scan` cannot read keeps today's text rule: the incident line, with
`for f in $(find …)`, stays refused; the brief line of section 7.2 and the
rule text say to write scratchpad work as a readable line or as a script file
in the scratchpad run by one plain command.

Residual, stated: a link inside a repository cloned on the same line is not
seen (`git clone` is not a linker here). The guards are tripwires
(`.claude/REQUIRED-READING.md`, *The harness*).

### 9.2 Live checks (after the merge, on master, since hooks run from the settings a session starts with)

1. A real `gh api graphql -f query='{ viewer { login } }'` runs without a
   prompt; recorded in this file.
2. The SessionStart recap shows the open PRs.
3. One real `git merge --ff-only origin/master` in a worktree, then one
   commit: the first is silent, the second gated; and both events carried
   `tool_use_id` (one record file written and consumed per call).

## 10. Tests `@tester` writes first

All in the harness suites (`-m harness`), on `harness_fixtures.make_repo`
copies; `COPIED` gains `.claude/hooks/record_bash_start.py`. Real git, no
network; `gh` is a fake on `PATH` where used.

### 10.1 GraphQL (test_guard_push.py; `graphql_read` pure, and end to end)

1. Pass (no decision): `-f query='{ repository(owner:"a", name:"b") {
   pullRequest(number:1) { isInMergeQueue } } }'`; `-f query='query Q {
   viewer { login } }'`; `-fquery=...`; `--raw-field=query=...`; `-F
   number=5 -f query='query($number:Int!){...}'`; with `--jq .data`,
   `--paginate`, `-H 'Accept: x'`, `-X POST`.
2. Ask, operation: `-f query='mutation { ... }'`; `query Q {...} mutation M
   {...}`; `subscription {...}`; a read whose text contains `# mutation`; a
   read followed by `&& echo mutation` (the word anywhere on the line).
3. Ask, shell expansion (review point 1): `-f query="$M"`; `-f
   query=$M`; `-f query=$'\x6d'utation...`; `-f query="$(printf
   mut)ation {...}"`; a backtick form; `-f query=mu{tation,}...` unquoted;
   `-f query=mu\tation...` unquoted; an unclosed quote (the text rule asks);
   and `literal` pure: `'$x'` is `True`, `"$x"`, `$x`, `a{b,c}`, `a*` and
   `~/x` are `False`.
4. Ask, fields and options: `-F query='{ viewer { login } }'` (query by
   `-F`, review point 1); `-F query=@q.graphql`; `--input q.json`; two
   `query` fields; `-X PUT`; an unknown option `--foo`; endpoint `/graphql`
   or `https://api.github.com/graphql`; `repos/a/b/issues -f title=x`
   (unchanged).

### 10.2 Scratchpad (test_guard_push.py, test_guard_governance.py, test_shell_scan.py, with `point_scratchpad`)

Fixture: `<pad>/c` made by `git init -q` (its own common directory);
`<pad>/w` a worktree of the test repository (`git -C <repo> worktree add`);
`<pad>/g` a directory whose `.git` file reads `gitdir: <repo>/.git`.

5. Push guard, pass: `git -C <pad>/c config user.name x`; `cd <pad>/c && git
   update-ref refs/x HEAD`; `git -C <pad>/c commit --amend -m m`; and, as
   regression pins, `git clone <url> <pad>/d && git -C <pad>/d commit -m m`.
6. Push guard, ask: `git -C <pad>/c push origin x`; `cd <pad>/c && cd <repo>
   && git config a.b c`; `git -C <pad>/c --git-dir=<repo>/.git config a.b
   c`; `GIT_DIR=<repo>/.git git -C <pad>/c config a.b c`;
   `GIT_CONFIG_GLOBAL=x git -C <pad>/c config a.b c`; `git -C "$S"
   config a.b c`; `git -C <pad>/link config a.b c` where `link` points to the
   repository.
7. Push guard, ask, common directory (review point 2): `git -C <pad>/w
   config a.b c`; `git -C <pad>/w update-ref refs/x HEAD`; `git -C <pad>/g
   config a.b c`; `git clone <url> <pad>/e && git -C <pad>/e config a.b c`
   (no repository yet); `git -C <pad>/c config --global a.b c`, `--system`,
   `--file <repo>/.git/config`, `-f <repo>/.git/config`; `cd <pad>/c && rm
   -rf .git && git config a.b c` (rule 5: a writer on the line).
8. Directory chain (review point 3), `directories(chain_only=True)` pure and
   both guards end to end: `(cd <pad>/c) ; git config a.b c`; `(cd <pad>/c)
   && git config a.b c` (the cd moves nothing: the repo is judged); `false &&
   cd <pad>/c; git config a.b c`; `cd <pad>/c || true && git config a.b c`;
   `! cd <pad>/c && git config a.b c`; `cd <pad>/c | git config a.b c`; and
   the governance forms `(cd <pad>) ; printf x >> CLAUDE.md`, `false && cd
   <pad>; printf x >> CLAUDE.md`: each asks; `chained` and `nested` set as
   section 4.2 says on each (`test_shell_scan.py`). With `chain_only=False`,
   `cd <wt>; git commit` resolves to `<wt>` (today's reading kept).
9. Governance, pass: `cd <pad>/c && printf x >> CLAUDE.md && git -C <pad>/c
   commit -qam m`; `cd <pad> && mkdir -p c/.claude/hooks && cp a.py
   c/.claude/hooks/`.
10. Governance, ask: the same with `cd <repo>`; `cd <pad> && ln -s
    <repo>/CLAUDE.md x && printf x >> x`; `cd <pad> && cp -R <repo>/.claude d
    && printf x >> d/agents/a.md`; `cd $S && printf x >> CLAUDE.md`; `cd
    <pad> && python3 -c '...' && printf x >> CLAUDE.md`. Existing one-plain-
    command cases unchanged.

### 10.3 The after-commit hook (test_gates_after_commit.py, new section, real git)

Each case runs the recorder copy on a Pre event, then really runs the git
command, then the hook copy on the Post event with the same `tool_use_id`
(the fixture's `Repos`, red detria stub in the worktree).

11. A conflicting `git merge`: exit 0, nothing printed, no gate ran, record
    file gone.
12. `git merge` already up to date; `git merge --ff-only` fast-forward;
    `git commit` with nothing to commit: each exit 0, silent, no gate.
13. A real commit in the worktree: exit 2, `RED IN <wt>`.
14. The `cd` case: line `git -C <wt> commit -qm x && cd <plain dir>`, Pre
    `cwd` the main repository, Post `cwd` the plain directory: the worktree is
    gated (`RED IN <wt>`), no `NOT CHECKED`.
15. `cd <plain dir> && git commit -m x` (git refuses there): silent.
16. `git commit -qm a && git merge <conflicting>` in the worktree: one `NOT
    CHECKED ... in progress` line, exit 2, no gate ran.
17. Fast-forward after a commit (review point 4): `git commit -qm a && git
    merge -q --ff-only <side>` in the worktree, `<side>` one commit ahead of
    the new HEAD: exit 2, `RED IN <wt>` (the newest entry is a fast-forward,
    the one before it a commit). Two fast-forwards on one line, `git merge
    --ff-only <s1> && git merge --ff-only <s2>`: silent. A record whose
    `heads` sha is not in the reflog (a made-up 40-hex sha): gated.
18. No record (Post event with an unknown `tool_use_id`): today's behaviour,
    `RED IN <wt>` for `git -C <wt> commit`. An unreadable record (review
    point 4), each written by hand under the call's `tool_use_id`: `not
    json`, `[]`, `{"cwd": "rel", "heads": {}}`, `{"cwd": "<abs>", "heads":
    {"<wt>": "xyz"}}`: each gives today's behaviour, `RED IN <wt>`, and the
    file is deleted.
19. The recorder writes nothing for `git status`, a missing `tool_use_id`,
    or `tool_use_id` `../x`; deletes a 25-hour-old file in `bash-start/` and
    keeps a 23-hour-old one (mtimes set with `os.utime`); exits 0 and prints
    nothing when `bash-start/` cannot be made (a file at its path).
20. Item 3: a red gate's stderr contains "new commit" and not "Fix and amend".

### 10.4 The recap (test_session_state.py, a fake `gh` like `_fake_ps`)

21. Two PRs, opened 3 days 4 hours and 50 minutes before the test's `now`:
    `Open PRs (2):`, the older first, ages `3d 4h` and `50m`, a 70-character
    title cut to 50 with `...`.
22. State (review point 5): three PRs with `mergeable` `CONFLICTING`,
    `UNKNOWN` and `MERGEABLE`: their lines carry `CONFLICTING`, `mergeable
    unknown` and no state word; the fake `gh` is called with `mergeable` in
    its `--json` list.
23. `[]`: `Open PRs: none`. 35 PRs: 30 lines then `... 5 more (gh pr list)`.
24. Unavailable, each with the rest of the recap still printed: no `gh` on
    `PATH`; `gh` exits 1 (reason from its first stderr line); `gh` sleeps
    30 s (the recap ends within 15 s); `gh` prints `{}` or `not json`.

### 10.5 Briefs (test_brief.py)

25. An increment with nine `QUOTED` lines, then `Acceptance (PR 3): @perf
    times the default spring on the 1 m mesh`: the `@reviewer` and
    `@developer` briefs carry the line verbatim after the Status line; the
    `@tester` brief does not; an `Acceptance:` line inside `## Review` is not
    quoted; a 400-character line is cut to 300 with `...`.
26. No `Acceptance` line: the `none stated in <file>` line.
27. With `SCRATCH_ROOT` pointed at a temporary base and
    `CLAUDE_CODE_SESSION_ID` set and its scratchpad made: the block names
    `<pad>/logs/<worktree name>`; without the session's scratchpad: the
    `rasputin-logs/<worktree name>` fallback under the temp directory; two
    worktrees get two folders.

**Mutation targets** (this increment's invariant-critical suites are 10.1,
10.2 and 10.3, the guards and the hook; the kill record covers each; the
number is the killing test):

- GraphQL: the `mutation` word test dropped (2); the word test applied to
  the query value only, not the line (2, `&& echo mutation`); `literal`
  not called (3, `$M`); `literal` ignoring double quotes (3, `"$M"`);
  `literal` without the brace check (3, `mu{tation,}`); `-F` accepted for
  `query` (4); the unknown-option refusal dropped (4, `--foo`).
- Scratchpad, push guard: the scratchpad test applied to a `str` directory
  (6, `"$S"`); the `--git-dir` refusal dropped (6); the `GIT_` text test
  narrowed to the old four names (6, `GIT_CONFIG_GLOBAL`); the push excluded
  from the pass dropped (6); the common-directory check dropped (7, `<pad>/w`);
  the common-directory check replaced by `--git-dir` (7, `<pad>/g`); the
  `--global`/`--system`/`--file` refusal dropped (7); rule 5 dropped (7,
  `rm -rf .git`).
- Directory chain: `chained` ignored (8, `false && cd <pad>/c; …`); `nested`
  ignored, a subshell's `cd` followed (8, `(cd <pad>/c) && …`); `||` read as
  `&&` (8); a stripped `!` not breaking the chain (8).
- Governance: the linker test dropped (10, `ln`).
- After-commit hook: resolving from the event's `cwd` (14); the
  HEAD-unchanged skip dropped (12); the in-progress test dropped (11 or 16);
  the fast-forward skip dropped (12); the fast-forward test reading only the
  newest entry (17, commit then fast-forward); the no-record fallback
  replaced by "skip" (18); an unreadable record read as "skip" (18).

## 11. Size and the split

Net production lines (`CLAUDE.md` §2), estimates:

| File | Item | Estimate |
|---|---|---|
| `.claude/hooks/guard_push.py` (`graphql_read`; scratchpad checks 1 to 5 with the `rev-parse` call) | 1, 5 | about 75 |
| `.claude/hooks/guard_governance.py` | 5 | about 20 |
| `tools/shell_scan.py` (`directories`, `git_directory` moved; `literal`; `chained`, `nested`) | 1, 4, 5 | about 65 |
| `.claude/hooks/gates_after_commit.py` (less the moved code; reflog walk; record check) | 3, 4 | about 35 |
| `.claude/hooks/record_bash_start.py` (new) | 4 | about 40 |
| `tools/session_state.py` (with the state word) | 2 | about 40 |
| `tools/brief.py` | 7, 8 | about 25 |
| Total | | about 300 |

Up from about 215 at design review round 1: the five review points add
about 85, under the 700-line cap with room.

`.claude/settings.json`'s new pre-hook entry is JSON (not counted):

```json
{ "matcher": "Bash", "hooks": [
  { "type": "command", "command": "$CLAUDE_PROJECT_DIR/.claude/hooks/record_bash_start.py" } ] }
```

**Two PRs.** PR A, prose only (items 6 and 9, and section 8's other text
held for PR B): `CLAUDE.md`, `reviewer.md`, README step 4; CI runs only the
governance gates, one review, and it can merge the same day. PR B, code
(items 1 to 5, 7 and 8, with their rule text): full path, red then green,
one code review, Ola's settings yes for the entry above, the live checks
of 9.2 after the merge. Both count as one item toward the three in flight.

## 12. Order

PR A: Ola's prompts for the three rule files; `@reviewer`; push on Ola's yes.
PR B: design review; red (section 10, all failing); green; code review;
Ola's yes for the settings entry and the rule text; push on his yes; live
checks after the merge, recorded here. All by day: every file is governed.

## Review

**Design review, round 1, 2026-10-14.** Commit `9b123d8e`, about 215 net lines estimated, @reviewer: CHANGES REQUESTED. (1) The read-only GraphQL rule passes a mutation through shell expansion (`$M`, `$'\x6d'`, `$(printf mut)`) or gh's `{branch}` fill in `-F`: pass only lines `shell_scan` reads, refuse expansion outside single quotes, require `query` through `-f`/`--raw-field`, test the words over the whole line. (2) A scratchpad worktree or gitfile shares the real repository's refs and config: require `git rev-parse --git-common-dir` under a scratchpad; refuse `config --global`, `--system`, `--file`. (3) The directory tracking follows `cd` in a straight line, so `(cd <pad>) ; …`, a failed `cd` and `false && cd <pad>; …` exempt a write to the real tree: pass only through `cd`s joined by `&&` at top level. (4) The fast-forward skip hides `git commit && git merge --ff-only x`: skip only when every reflog entry since the recorded HEAD is a fast-forward; an unreadable record counts as none. (5) Item 2 drops the open PR's state asked for with it: add `mergeable` or say why not. A mutation target and tests for each. Held: items 3, 6, 7, 8 and 9 match their sources; no history in rule text; citations resolve.

**Design review, round 2, 2026-10-14.** Range `8f3785ca..d42f04ee`, about 300 net lines estimated, @reviewer: CHANGES REQUESTED. Points 4 and 5 closed. Still open, each broken by a probe: (1) quotes split the word, `-f query='query Q { a } mut''ation M { b }' -f operationName=M` (or `"mut"ation`) passes rule 6 and gh sends `mutation M`: also test every argument after quote removal; (2) `git -C <pad>/c config core.worktree <repo>` passes rules 1-5 and a later `reset --hard` there would overwrite the real tree: require `--show-toplevel` under a scratchpad by real path, or refuse `core.worktree`; (3) zsh's `cd OLD NEW`, `pushd +1` and a `CDPATH=…` prefix move elsewhere than the guards compute: give a path only for a `cd`/`pushd` of one word that is absolute or starts `./` or `../`, with no assignment prefix. A test and a mutation target for each.
