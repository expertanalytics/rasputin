# Harness h9: briefs come from files, and spawns are checked

Status: **design**, @architect, 2026-10-04, at master 2060f14. One PR, by
day: it edits `.claude/settings.json` and other governed files. Implements
stage 1 of the plan in `docs/retrospectives/2026-10-03-dispatcher-control.md`
§5, which Ola approved on 2026-10-04: "1: yes to stage 1, then stages 2 and
3". His approval of stage 1 includes the hook, so the settings change in §6
is approved as written there and nothing more. Evidence: rows 1 to 8 of that
file's §2 (errors in writing a brief).

Words used below: a **brief** is the prompt the main session writes when it
starts a persona (`Agent`) or sends a follow-up to a running one
(`SendMessage`); the **block** is the fixed part of a brief that
`tools/brief.py` prints; a **resume** is a `SendMessage` to a persona that is
already running or has stopped.

## 1. Prior art

Tooling; no novelty claimed.

*Literature.* Anthropic's account of its research system says what a
delegation must carry: "Each subagent needs an objective, an output format,
guidance on the tools and sources to use, and clear task boundaries", and
"Without detailed task descriptions, agents duplicate work, leave gaps, or
fail to find necessary information" ("How we built our multi-agent research
system", 2025-06-13, as quoted in the plan's §4, item 4). The block is that
list, made fixed: the files to read (guidance), the worktree and write limit
(boundaries), the handback headings (output format). The task stays free
text. Claude Code's documentation gives the mechanism and the reason for it:
memory and `CLAUDE.md` are "context, not enforced configuration. To block an
action regardless of what Claude decides, use a PreToolUse hook"
(code.claude.com/docs/en/memory); a `PreToolUse` hook may answer
`permissionDecision: "deny"` with a reason; and of the working directory,
"`${CLAUDE_PROJECT_DIR}` stays put: it still points at the project root where
the session started", while "`cwd` follows Claude: the `cwd` field in the
hook's input JSON is the worktree root after Claude enters a worktree, and the
new directory after Claude runs `cd`" (code.claude.com/docs/en/hooks, read
2026-10-04). The block's check is a content hash, as git names objects by
theirs: it shows the text was not changed after it was made, and claims
nothing more (§3.4).

*Legacy.* Nothing; the legacy tree has no harness.

```
$ git grep -l -iE 'brief|PreToolUse|subagent|CLAUDE_PROJECT_DIR' legacy-archive -- legacy
(no output, exit 1)
```

## 2. Scope

| # | What | Where |
|---|---|---|
| 1 | Print the block for a persona: files to read, worktree, note file, write limit, concurrency, the increment file's own lines on required steps and its review rounds, Ola's words checked against the transcript | new `tools/brief.py` |
| 2 | The fixed text: one shared template and one per persona | new `.claude/briefs/common.md` and six `.claude/briefs/<persona>.md` |
| 3 | Refuse a persona spawn without an intact, current block; refuse any spawn or resume while the session is not in the directory it started in | new `.claude/hooks/guard_spawn.py`, `.claude/settings.json` |
| 4 | Delete what the templates replace | `CLAUDE.md` §3 "Briefs" bullet; memory note `lean-agent-briefs.md` |
| 5 | Keep the new rule text governed and measured | `guard_governance.py`, `tools/rule_sizes.py` |

Out: the next free increment id (row 6 of the evidence) is stage 2's
`tools/pipeline.py`; concurrency counting by hook is stage 4. Stage 1 stops
rows 1 to 5, 7 and 8, not row 6.

## 3. Design

### 3.1 `tools/brief.py`

Run by the main session, from the main checkout (the hook refuses spawns from
anywhere else, §3.5), so the main checkout's copy is the one both sides use.

```
python3 tools/brief.py <persona> --worktree <path> --beside <none|persona>
                       [--increment <file>] [--no-build] [--ola <text>]...
```

- `<persona>`: one of `PERSONAS`, the template names in `.claude/briefs/`
  other than `common` (architect, developer, orchestrator, perf, reviewer,
  tester). Anything else: exit 2.
- `--worktree`: an existing checkout of this repository
  (`git -C <path> rev-parse --show-toplevel` resolves to the path; the main
  checkout is allowed). Printed absolute and resolved. A path containing
  whitespace, or not a checkout: exit 2.
- `--increment`: required for tester, developer and reviewer, and must exist
  (exit 2 otherwise). For architect it may name a file that does not exist
  yet: the block then says `<file> (new: you create it)`. Optional for perf
  and orchestrator.
- `--beside`: required, so the concurrency choice is made every time
  (memory note "agents in pairs"). Exit 2 when the persona is `perf` and
  `--beside` is not `none`, or when `--beside perf`: nothing runs beside a
  timing run. `--no-build` adds that this persona may not build C++ in this
  run (only one C++ build at a time).
- `--ola <text>`, repeatable: a quotation of Ola. Each is checked against the
  human turns of the running session's transcript,
  `session_state.TRANSCRIPTS / f"{CLAUDE_CODE_SESSION_ID}.jsonl"`, read with
  `session_state.human_turns` (which already includes prompts absorbed
  mid-turn). It is found when its whitespace-collapsed text is a substring of
  one turn's text. Printed as given (whitespace collapsed), in double quotes,
  with that turn's timestamp. Not found, no session id, or no transcript:
  exit 2, naming the quotation. This turns row 1 (a paraphrase printed as his
  words) into a refusal at brief time.
- The note file: `<main checkout>/.claude/current-task/<persona>-<HHMMSS>.md`,
  local time, main checkout from `session_state.main_checkout`. `brief.py`
  names it and does not create it (`REQUIRED-READING.md`: the subagent writes
  it, the spawner deletes it).

`WRITES`, a dict in `brief.py`, is the write limit per persona, taken from
h6's table (`docs/increments/h6-role-limits.md` §3), plus "your note file":

| persona | may write |
|---|---|
| tester | `tests/` |
| developer | `src_python/`, `include/`, `src/`, `bindings/`, `tools/`, `.claude/hooks/`, `.github/`, `CMakeLists.txt`, `pyproject.toml` |
| perf | `docs/benchmarks/` |
| architect | `docs/` except `docs/retrospectives/`, `ROADMAP.md`, `CLAUDE.md`, `.claude/` files ending in `.md` |
| orchestrator | `docs/retrospectives/` |
| reviewer | nothing |

When h6 lands (it waits on h5), `WRITES` is replaced by an import of h6's
`ROLES`, so the limit has one statement; whichever of h6 and h9 merges second
makes that change.

**From the increment file**, printed as the file's own lines, never
summarised (rows 3 and 5: required mutation tests dropped from a brief):

- the first line starting `Status:`;
- outside the `## Review` section, every line matching
  `(?i)invariant-critical|mutation|@perf|acceptance run|tools/bench\.py`, as
  `  <lineno>: <line stripped>`, at most 8, then
  `  ... <n> more; read them in the file`;
- `Review rounds recorded: <n>`, where n counts the lines inside
  `## Review` (to the next `## ` heading or the end) that contain `APPROVED`
  or `CHANGES REQUESTED`, and, when n > 0, `last: <that line, cut to 200
  characters>`; with no `## Review` section, `Review rounds recorded: none`.
  The entries' formats vary (h8 and the 2x files differ); stage 2 fixes the
  format, and this count is a pointer, not a ledger.

### 3.2 The block

```
<<<BRIEF persona=<p> worktree=<abs path> head=<40-hex> hash=<12-hex>>>>
<common.md, filled in>
<p>.md, filled in>
Concurrency: <generated, §3.1>
Write limit: <WRITES[p]>, and your note file.
From <increment>, its own lines:
  ...
Review rounds recorded: ...
Ola, verbatim, checked against this session's transcript:
  "<text>" (<timestamp>)
<<<END BRIEF <12-hex>>>>
```

- `head` is `git -C <worktree> rev-parse HEAD` when the block is made.
- Templates are filled with `string.Template.substitute`, variables
  `$persona`, `$worktree`, `$increment`, `$note`; an unknown `$name` in a
  template is an error at brief time, not a silent blank.
- Lines without content are omitted: no `From` part without `--increment`,
  no `Ola` part without `--ola`.
- `block_hash(persona, worktree, head, body) -> str`: the first 12 hex
  digits of SHA-256 over `f"{persona}\n{worktree}\n{head}\n{body}"`, where
  `body` is the lines between the two markers, each with trailing whitespace
  removed, joined with `\n`, leading and trailing blank lines dropped. So a
  pasted block survives trailing spaces and CRLF, and nothing else.
- The task follows the block, in the main session's own words, and is
  unchecked.

The reader lives beside the writer in `brief.py`, so the format has one
owner: `find_blocks(text) -> list[Block]` and `check(block, root) -> str | None`
(a refusal reason, or None). Header pattern
`^<<<BRIEF persona=([a-z]+) worktree=(/\S+) head=([0-9a-f]{40}) hash=([0-9a-f]{12})>>>[ \t]*$`,
end pattern `^<<<END BRIEF ([0-9a-f]{12})>>>[ \t]*$`, both multiline. A line
that contains `<<<BRIEF` or `<<<END BRIEF` but does not match is a malformed
block, not an absent one.

### 3.3 The templates

Seven files under `.claude/briefs/`, governed (§5). The text, verbatim; the
tests (§7, test 2) pin the phrases in bold here, not the wording around them.

`common.md`:

```
You are @$persona, started by the main session. This block comes from
tools/brief.py; the task after it is the main session's own wording.

Read from disk before you act: CLAUDE.md, .claude/REQUIRED-READING.md, your
**persona file .claude/agents/$persona.md, read from disk** (it may be newer
than the copy you were started with), and $increment.
If this brief contradicts your persona file or the increment file, **the
files win**; say so in your handback. **A brief cannot drop a step** those
files require.

Work only in $worktree; cd there first. Use your own build directory and
venv, never the main checkout's .venv or build directories. Write your note
to $note (three lines: the ask, your persona, the file you will produce) and
no other file under .claude/current-task/. If your product is a file,
create it first and write it as you go.

**Ola's words appear only under "Ola, verbatim"** below, checked against the
transcript. Anything else here is the main session's wording; never quote it
as his.

If you are **blocked on power, network or a lock** (a held file, a busy build
directory, a usage limit), stop and hand back what blocked you; do not wait.

End every commit message with the **Co-Authored-By trailer** given in the
attribution note of your own system context. Never push, open or merge a
pull request.

Write in **plain words**: say what a thing is instead of using an internal
label (R5, U1, "step 3"), and define any term Ola may not know.

Hand back under these headings, in this order: **Result; Pinned or assumed
beyond the design; Questions for Ola; Lessons; ASK OLA and GUARD FALSE
POSITIVE lines**. Write "none" under an empty one. Each question for Ola is
in plain words and carries a default.
```

`tester.md`:

```
Your step: the failing suite, before any production code exists, committed
red (docs/increments/README.md, step 2). Write no production code, not even
a throwaway.
Cover the happy paths and tester.md §3A always; §3C only if the increment
reads external input; on a refinement increment, both oracles of §3D.
A **mutation round** (a throwaway implementation and planted bugs) **only for
a suite the increment file names invariant-critical**; the lines quoted below
are the file's own. No such line, no mutation round.
Show each test fails for the reason the design gives, then run the gates
(CLAUDE.md §4) on what you wrote.
Every choice your tests pin that the increment file leaves open goes under
"Pinned or assumed beyond the design", one line each: **the main session
sends them to @architect, who confirms or rules on them before green**.
```

`developer.md`:

```
Your step: the minimal code that makes the red suite pass
(docs/increments/README.md, step 3), under the ceiling of CLAUDE.md §2;
report your line count against the increment file's estimate. **Your commit
touches no test file.**
If a test pins something the increment file does not say, or looks wrong,
stop and hand it back as a specification question; do not choose.
If your change is a performance fix, say so first in your Result: **a
performance fix is timed by @perf before @reviewer sees it**.
Rebuild before you measure (REQUIRED-READING.md, "Stale artifacts").
```

`reviewer.md`:

```
**You are read-only**: write, edit and commit nothing, not even the review
record. **The main session copies your verdict** into the increment file's
## Review section.
This is a request for your review; green CI is not done. Make **the three
checks of reviewer.md §5**, not what the gates cover; on a branch without
production code, the scope in REQUIRED-READING.md, "Before you publish".
Your Result, in this order: Verdict; the commit range reviewed; the
production line count by CLAUDE.md §2's rule, against the estimate; Blocking
issues; Suggestions.
If the increment touches refine or mesh code, say whether @perf's acceptance
run is recorded (docs/increments/README.md, "Acceptance").
```

`architect.md`:

```
Your step: the design in $increment, before @tester writes anything
(docs/increments/README.md, step 1): types, invariants, exclusions,
degeneracy policy, the LOC estimate, and the tests @tester writes red. No
production code.
Write the **Prior art** section first, with the legacy grep pasted with what
it returned.
When settling choices @tester pinned beyond the design, confirm or rule on
each in the increment file, and list the tests that change.
Questions for Ola only if unavoidable.
```

`perf.md`:

```
Your step: as perf.md says: the acceptance run, or the timing asked below,
from tools/bench.py, with the power state recorded and like compared with
like (docs/increments/README.md, "Acceptance"). Figures only from finished
runs.
**Nothing runs beside a timing run**; the concurrency line below says so.
A change to tools/bench.py follows the test-first loop (perf.md §1): report
the change it needs; do not make it.
**A performance fix by @developer is timed by you before @reviewer sees it.**
```

`orchestrator.md`:

```
Your step: the check, retrospective or research asked below (orchestrator.md
§1 to §3). A change to a rule or to the harness goes in your report as a
proposal, with its evidence.
Measure the rule text (python3 tools/rule_sizes.py) and propose a cut.
Cite transcripts as session and line; **quote Ola only from the transcript**.
```

The concurrency line `brief.py` generates: `Concurrency: you run alone.`;
with `--beside X`, `Concurrency: @X runs at the same time, in its own
worktree. You may build C++ in this run.` (or `You may not build C++ in this
run.` with `--no-build`). For perf, `Concurrency: you run alone; nothing runs
beside a timing run.`

What the templates carry, against the recurring lines of the brief for this
increment:

| Recurring line | Where |
|---|---|
| pairs; alone beside @perf | `--beside`, `--no-build`, the generated concurrency line, the exit-2 rule |
| stop and report if blocked on power, network or a lock | `common.md` |
| the attribution trailer | `common.md`, by reference to the persona's own system context (this @architect run received the trailer that way), so a model change does not leave a stale literal in a governed file |
| Ola's words pasted verbatim | `--ola`, checked against the transcript; `common.md` |
| red-step choices go to @architect before green | `tester.md`; the handback heading in `common.md` |
| @reviewer is read-only; the spawner records | `reviewer.md` |
| a perf fix is timed by @perf before review | `developer.md`, `perf.md` |
| plain language | `common.md` |
| the "Briefs" bullet of `CLAUDE.md` §3 | `tester.md` (§3A, §3C, §3D), `developer.md` (minimal code, ceiling), `reviewer.md` (§5's three checks, explicit request) |
| the "lean briefs" note | `tester.md` (mutation only where the increment file names it, quoted) |

### 3.4 What the hash proves

That the block between the markers is what `brief.py` printed for that
persona, worktree and commit. It is not a secret: anything that can run
Python can compute it, and the main session can. It is a tripwire for the
errors of rows 2 to 7, which were made by editing or not including the fixed
text, not by forging it. A keyed hash (HMAC) would add a key the main session
can read, so it would prove nothing more; rejected.

### 3.5 `.claude/hooks/guard_spawn.py`

`PreToolUse`, matcher `Agent|SendMessage` (§6). Reads the event from stdin;
always exits 0; a refusal is
`{"hookSpecificOutput": {"hookEventName": "PreToolUse", "permissionDecision": "deny", "permissionDecisionReason": <reason>}}`,
as `guard_unattended.py` builds it. Input fields used: `tool_name`, `cwd`,
`agent_id`, and `tool_input.prompt`, `tool_input.subagent_type` and
`tool_input.isolation` (Agent) or `tool_input.message` (SendMessage). The
input field names are those of 527 `Agent` and 338 `SendMessage` calls in
this project's transcripts (2026-10-04); the hook input itself is not
documented per tool, so §9's probe confirms it.

In this order; the first that applies decides:

| # | Condition | Verdict | Reason starts with |
|---|---|---|---|
| 1 | stdin not a JSON object; `tool_name` not `Agent` or `SendMessage`; or `agent_id` present (the event comes from inside a subagent) | pass, no output | |
| 2 | `cwd` is not the same directory as `$CLAUDE_PROJECT_DIR` (§3.6) | deny | `cwd:` and both paths, and `run: cd <project dir>` |
| 3 | `brief` cannot be imported from the hook's own `tools/` | pass, with `additionalContext` (§3.7) | |
| 4 | the text holds a malformed marker, or more than one block | deny | `brief: malformed` / `brief: two blocks` |
| 5 | Agent, `subagent_type` in `PERSONAS`, no block | deny | `brief: no block` and the command to print one |
| 6 | a block whose hash does not match | deny | `brief: edited` and "paste the output of tools/brief.py unchanged" |
| 7 | Agent, block persona differs from `subagent_type` | deny | `brief: persona` and both names |
| 8 | the block's worktree is not a checkout, or its `HEAD` is not the block's `head` | deny | `brief: worktree` / `brief: stale` and "run brief.py again" |
| 9 | Agent with a block and a non-empty `isolation` | deny | `brief: isolation`: the persona works in the block's worktree, not a new one |
| 10 | otherwise | pass, no output | |

Rule 1 passes subagent events because subagents do not start personas here,
and a subagent's `cwd` is its worktree by design. Rule 5 asks no block of
other subagent types (`general-purpose`, `Explore`, forks): they have no
template, and h6 gives any name outside the persona table no write rights, so
they cannot do a persona's work once h6 lands. Until then this is a gap,
named in §8. Rule 8 catches a block reused from an earlier step: any commit
in the worktree makes it stale, so a brief is made just before its spawn.

Unattended mode does not change any verdict, and nothing is queued: every
refusal is one the main session fixes itself (a `cd`, a run of `brief.py`),
so none waits for Ola. `harness_mode` is not imported.

### 3.6 The working-directory comparison, settled

`os.path.samefile(event["cwd"], os.environ["CLAUDE_PROJECT_DIR"])`, with
`OSError` (a directory deleted since) counting as different. Same file, not
same string: it ignores trailing slashes, symlinks (`/tmp` and
`/private/tmp`) and the case of a path on a case-insensitive volume. Any other
directory is refused, a subdirectory of the project included: a persona
inherits the main session's directory, and the briefs' relative paths assume
the project root. A session Ola starts inside a worktree has that worktree as
its `$CLAUDE_PROJECT_DIR`, so it is not refused (the plan's requirement). A
session that entered a worktree with `EnterWorktree` is refused until it
leaves it (`ExitWorktree`), which is the same hazard as row 8 by a different
route. With `cwd` missing from the event or `$CLAUDE_PROJECT_DIR` unset (a
hand run), the comparison is skipped and said so in `additionalContext`.

### 3.7 When a part is unavailable

- **`tools/brief.py` missing or failing to import, or the hook failing after
  it has read the event.** The spawn goes through, with `additionalContext`:
  `guard_spawn: brief check unavailable (<error>); this spawn was not
  checked. Tell Ola; repairing tools/brief.py is a @developer task.` The
  working-directory check (rule 2) needs no `brief.py` and still applies. This
  is the opposite of `guard_unattended.py`, which refuses on a crash, and the
  reason is the cost of each mistake: there, letting a prompt through stalls
  a night; here, refusing every spawn would also refuse the @developer that
  repairs `brief.py`, and the main session does not write code.
- **The hook file missing or not executable.** Claude Code then runs the
  tool and shows a non-blocking hook error (code.claude.com/docs/en/hooks:
  "When the script path doesn't exist or isn't executable ... For most hook
  events, the action proceeds"). The wiring test (§7, test 24) checks that
  the file exists and is executable in git's index.
- **`brief.py` itself refuses** (exit 2): it prints the reason on stderr and
  no block; the main session fixes its arguments.

### 3.8 Resumes

A `SendMessage` needs no block: the persona was started with one and keeps
it in its context. The working-directory rule applies to every resume
(row 8 was a resume). A message that holds a block gets rules 4, 6 and 8
(rule 7 cannot apply: `to` is an agent id, not a persona name). So a resume
that starts a new step, such as a @tester amendment after a ruling, may carry
a fresh block and has it checked; the rule text (§4) says it should. That a
new step can be sent without one is named in §8.

### 3.9 Considered and rejected

**The hook writes the block in.** The main session would write one line,
`BRIEF tester --worktree ...`, and the hook would replace it with the block
through `updatedInput`, so nothing is copied by hand. Rejected for now:
whether `updatedInput` applies to `Agent` is not documented, and if it did
not, the persona would start with one line and no block, which nobody would
see. Worth a probe in stage 4; it would make the hash unnecessary.

## 4. Rule text, and what it replaces

`CLAUDE.md` §3, "The main session dispatches": the whole "Briefs" bullet is
deleted (55 words, `len(text.split())` at 2060f14), from `* **Briefs.**` to
`since green CI is not done.` No replacement in `CLAUDE.md`.

`.claude/REQUIRED-READING.md`, *The harness*, first paragraph: after
"configuration changes while unattended mode is on;" insert (47 words):

```
`guard_spawn.py` refuses a persona spawn whose prompt lacks an unchanged,
current block from `python3 tools/brief.py`, and any spawn or resume while
the session is not in the directory it started in; a resume needs no block,
but one that starts a new step carries a fresh one;
```

Net rule text: -8 words, plus the seven templates (826 words, in the
size table's new row, §5).

**Outside the repository:** the memory note
`~/.claude/projects/-Users-skavhaug-projects-rasputin/memory/lean-agent-briefs.md`
and its line in that folder's `MEMORY.md`. The main session removes both,
after the merge and after the hook's first live refusal (§9), since memory is
its own and no persona writes there. Its 2026-09-26 postscript ("no PNGs
unasked") is already `CLAUDE.md` §3's "No unasked images"; its "check a step
that runs past about 20 minutes" is the "agents report blocks" note. Nothing
is lost. The "agents in pairs" note stays until stage 4's hook counts running
personas; it then goes too.

## 5. Files

| File | Change | Production lines (`CLAUDE.md` §2) |
|---|---|---|
| `tools/brief.py` | new, §3.1, §3.2 | ~130 |
| `.claude/hooks/guard_spawn.py` | new, §3.5 to §3.8 | ~50 |
| `.claude/settings.json` | §6 | 6 |
| `.claude/hooks/guard_governance.py` | `"tools/brief.py"` in `GOVERNED` (the hook imports it, so it is live before review, like the self-protecting set); `".claude/briefs/"` in `GOVERNED_PREFIXES` | 2 |
| `tools/rule_sizes.py` | `.claude/briefs/*.md` as one row, words summed, in the current and the reference counts; ordered after the skills | ~10 |
| `.claude/briefs/*.md` | new, §3.3 | prose |
| `CLAUDE.md`, `.claude/REQUIRED-READING.md` | §4 | prose |
| `tests/python/harness_fixtures.py` | `tools/brief.py`, `.claude/hooks/guard_spawn.py` and the seven templates in `COPIED` | test |
| **Total** | | **~200**, under the 700 ceiling |

The plan said about 140; the difference is the block reader beside its
writer, the transcript check of `--ola`, and the concurrency arguments.

One row, not seven, in the size table: h8's test 11 holds the recap's new
sections to 3,000 characters, and its design review measured about 280 of
slack, less than seven rows of about 42 characters.

No suite here is invariant-critical: no mutation round. No refine or mesh
code: no `@perf` acceptance.

**`ROADMAP.md`:** no row. `docs/increments/README.md` asks a merge to update
"`ROADMAP.md`'s row for that increment", and h2 to h8 have none; whether
harness increments get rows is item 7 of `docs/retrospectives/next.md`,
"The window of 2026-10-02, the restart of 2026-10-03, h7", still open. Its
ruling covers h9 with the rest.

**Overlap.** h5 (`worktree-h5-state-check`) moves `GOVERNED` into
`tools/governed.py` and edits `harness_fixtures.py`; h6 replaces `WRITES`
(§3.1). Whichever merges second resolves a few lines.

## 6. The settings change

`.claude/settings.json` at 2060f14, lines 23-28 are the `AskUserQuestion`
entry. Line 28, `      }`, becomes `      },`, and these six lines follow it,
before line 29 (`    ],`):

```json
      {
        "matcher": "Agent|SendMessage",
        "hooks": [
          { "type": "command", "command": "$CLAUDE_PROJECT_DIR/.claude/hooks/guard_spawn.py" }
        ]
      }
```

Nothing else in the file changes. @developer makes the edit in the green
commit, by day; `guard_governance.py` asks Ola at the write, and that prompt
is the edit's own yes. It goes live in Ola's sessions only once merged into
the main checkout.

## 7. Tests for @tester (red, before any code)

New `tests/python/test_brief.py` and `tests/python/test_guard_spawn.py`;
additions to `test_settings_wiring.py`, `test_guard_governance.py` and
`test_rule_sizes.py`. Fixture repositories from `harness_fixtures.make_repo`,
worktrees with `git worktree add`, the hook run by path with
`$CLAUDE_PROJECT_DIR` set to the fixture, as `test_settings_wiring.py` does.
Transcripts are fixture `.jsonl` files under a temporary `HOME`, with
`CLAUDE_CODE_SESSION_ID` set.

`brief.py`:

1. **Shape.** For each persona, the output's first line matches §3.2's header
   pattern and its last the end pattern with the same hash; `head` is the
   worktree's `HEAD`; `worktree` is absolute and resolved; `block_hash` over
   the body reproduces it. Two runs in the same second give the same block.
2. **Templates.** Every block holds the bold phrases of `common.md` (§3.3);
   each persona's block holds the bold phrases of its own template and not
   another persona's (tester's mutation sentence is not in developer's
   block). An unknown `$name` in a fixture template: non-zero exit.
3. **Increment lines.** A fixture increment with a `Status:` line, two
   matching lines in the body, one matching line inside `## Review`, and a
   `## Review` with two verdict lines: the block quotes the status line and
   the two body lines with their numbers, not the review line, and says
   `Review rounds recorded: 2` with the last verdict line. Ten matching lines:
   8 and `... 2 more`. No `## Review`: `none`. A verdict line over 200
   characters is cut.
4. **Increment required.** tester, developer, reviewer without `--increment`,
   or with a missing file: exit 2. architect with a missing file: exit 0 and
   `(new: you create it)`. perf and orchestrator without one: exit 0, no
   `From` part.
5. **Worktree.** Not a checkout; a path with a space: exit 2. The main
   checkout: accepted.
6. **Concurrency.** `--beside none`: "you run alone"; `developer --beside
   tester`: names @tester and "may build"; with `--no-build`: "may not build";
   `perf --beside tester` and `tester --beside perf`: exit 2. No `--beside`:
   exit 2 (argparse).
7. **Ola.** A quotation present in a fixture human turn, with different
   spacing and a line break: printed collapsed, quoted, with the turn's
   timestamp. Present only in an assistant message or a tool result: exit 2.
   Present in an absorbed queued prompt: found. No session id, or no
   transcript: exit 2.
8. **Note file and write limit.** The note path matches
   `<main checkout>/.claude/current-task/<persona>-\d{6}\.md`, also when run
   with a worktree other than the main checkout; the file is not created.
   The write-limit line equals `WRITES[persona]`; reviewer's says nothing but
   the note file.
9. **Hash.** Changing one character of the body, the persona, the worktree
   or the head changes the hash; trailing spaces and CRLF line ends do not.

`guard_spawn.py`, each verdict of §3.5:

10. A fresh block for tester, `subagent_type: tester`, `cwd` the project:
    no output, exit 0.
11. Rule 5: tester without a block: deny, reason starts `brief: no block` and
    names `python3 tools/brief.py`. `general-purpose`, `Explore` and a
    missing `subagent_type` without a block: no output.
12. Rule 6: one body character changed; a body line deleted: deny, `brief:
    edited`. Trailing spaces added to every line: passes.
13. Rule 4: two blocks; a header with a 39-digit head; an END line missing;
    an END hash that differs from the header's: deny, `brief: malformed` or
    `brief: two blocks`.
14. Rule 7: a tester block on a developer spawn: deny naming both.
15. Rule 8: a commit in the worktree after the block was made: deny, `brief:
    stale`. The worktree removed: `brief: worktree`.
16. Rule 9: a valid block with `isolation: "worktree"`: deny.
17. Rule 2: `cwd` a worktree of the fixture, a subdirectory of the project,
    or a deleted directory: deny for Agent and for SendMessage, reason starts
    `cwd:` and names both paths. `cwd` a symlink to the project: passes. Rule
    2 is checked before the block: a valid block from the wrong directory is
    refused with `cwd:`.
18. Resumes: SendMessage with no block from the project: no output; with a
    fresh block: no output; with an edited or stale block: deny.
19. Rule 1: `agent_id` present with a wrong `cwd` and no block: no output.
    `tool_name: Bash`: no output. Stdin not JSON, or a JSON list: no output,
    exit 0.
20. Unavailable (§3.7): `tools/brief.py` removed from the fixture, or made to
    raise on import: a blockless tester spawn passes with `additionalContext`
    containing `brief check unavailable`; a wrong `cwd` is still denied.
21. Unattended: with the fixture's unattended flag on, test 11's verdicts are
    the same, and no queue file appears under the fixture's harness
    directory.
22. Missing `cwd` in the event, or `CLAUDE_PROJECT_DIR` unset: the
    comparison is skipped, said in `additionalContext`; the block checks
    still run.

Wiring and governance:

23. `test_guard_governance.py`: writes to `tools/brief.py`,
    `.claude/briefs/common.md` and `.claude/briefs/tester.md` get `ask`.
24. `test_settings_wiring.py`: a `PreToolUse` entry with matcher
    `Agent|SendMessage` and command
    `$CLAUDE_PROJECT_DIR/.claude/hooks/guard_spawn.py`; the file is
    executable on disk and in git's index (mode 100755); run by path from a
    fixture copy it denies a blockless tester spawn. The existing entries
    unchanged.
25. `test_rule_sizes.py`: with three files under `.claude/briefs/`, the table
    has one row `.claude/briefs/*.md` with their summed words, after the
    skills; at a reference without the folder the row says `new`.
26. Rule text: `CLAUDE.md` has no line starting `* **Briefs.**`;
    `REQUIRED-READING.md` names `guard_spawn.py` and `tools/brief.py` in *The
    harness*.

## 8. Not covered, named

- A resume that starts a new step without a block (§3.8).
- Spawns of types without a template (`general-purpose`, `Explore`, forks)
  until h6's role table gives them no write rights.
- The task text after the block: it can still paraphrase, contradict or
  drop things; the block's sentence "the files win" is the persona's
  defence, not a check.
- A block computed by hand rather than printed (§3.4).
- How many personas run at once (stage 4).
- The next free increment id (stage 2).

## 9. The first live firing

After the merge, by day, with Ola present, the main session checks the hook
does what the tests say, against the real tool input: (a) a @tester spawn
with no block is refused with `brief: no block`; (b) after `cd` into a
worktree, a `SendMessage` is refused with `cwd:`; then `cd` back; (c) a spawn
with a fresh block goes through, and the persona's handback has the five
headings. If (a) or (b) passes through, the tool input fields of §3.5 differ
from the transcripts', and the hook is reported to Ola before any other
spawn. Then the main session removes the memory note (§4).

## 10. Questions for Ola

None. Two choices made here that Ola may overturn: built-in helper types
need no block (§3.5, rule 5), and a broken `brief.py` lets spawns through
with a notice rather than stopping them (§3.7).
