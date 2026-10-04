# Harness h9: briefs come from files, and spawns are checked

Status: **design, round 2** (after design review round 1, below), @architect,
2026-10-04, at master 2060f14. One PR, by
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
| 2 | The fixed text: one template, holding only what the spawn alone knows; persona rules stay in the persona files (§3.3) | new `.claude/briefs/common.md`; one line each in `tester.md` and `reviewer.md`, one clause in `CLAUDE.md` |
| 3 | Refuse a persona spawn without an intact, current block; refuse any spawn or resume while the session is not in the directory it started in | new `.claude/hooks/guard_spawn.py`, `.claude/settings.json` |
| 4 | Delete what the block replaces | `CLAUDE.md` §3 "Briefs" bullet, two clauses of `REQUIRED-READING.md`; memory notes `lean-agent-briefs.md` and `agent-concurrency-pairs.md` (§4) |
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

- `<persona>`: one of `PERSONAS`, the keys of `WRITES` below (architect,
  developer, orchestrator, perf, reviewer, tester). Anything else: exit 2.
- `--worktree`: an existing checkout of this repository
  (`git -C <path> rev-parse --show-toplevel` resolves to the path; the main
  checkout is allowed). Printed absolute and resolved. A path containing
  whitespace, or not a checkout: exit 2.
- `--increment`: required for tester, developer and reviewer, and must exist
  (exit 2 otherwise). For architect it may name a file that does not exist
  yet: the block then says `<file> (new: you create it)`. Optional for perf
  and orchestrator.
- `--beside`: required, so the concurrency choice is made every time. Either
  `--beside none` alone, or one or more `--beside <persona>:<worktree>`, one
  per persona already running (worktree resolved as for `--worktree`, and it
  must exist). §3.1a gives the rules. `--no-build` adds that this persona may
  not build C++ in this run (only one C++ build at a time).
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

### 3.1a Concurrency

Ola's ruling of 2026-10-04, verbatim: "Read only agents should be allowed
even when two agents are working, as long as they don't require the same
files, ie that the files the read-only agent reads could change." It refines
the memory note "agents in pairs" (at most two at once, each in its own
worktree; nothing beside a `@perf` timing run).

A persona is **read-only** when its `WRITES` entry is empty (today only
`reviewer`); every other persona is a **writer**. Derived, not flagged, so
the two cannot disagree. Over the persona being briefed and every `--beside`
entry, `brief.py` exits 2, naming the rule, when:

1. `perf` is among them with anyone else: nothing runs beside a timing run
   (read-only agents included: a timing run needs a quiet machine);
2. there are more than two writers;
3. two writers share a worktree;
4. a read-only persona shares a worktree with a writer: the files it reads
   could change under it. Worktrees are the measure of "the same files"; a
   reader that reads another worktree's files is not seen (§8).

The generated line:

- `--beside none`: `Concurrency: you run alone.`
- perf: `Concurrency: you run alone; nothing runs beside a timing run.`
- otherwise: `Concurrency: also running: @tester (writer, in <path>), @reviewer
  (read-only, in <path>). At most two writers run at once; read-only agents
  do not count toward the two, as long as no writer changes the files they
  read.` then `You may build C++ in this run.` or, with `--no-build`, `You
  may not build C++ in this run.`

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
Concurrency: <generated, §3.1a>
Write limit: <WRITES[p]>, and your note file.
From <increment>, its own lines:
  ...
Review rounds recorded: ...
Ola, verbatim, checked against this session's transcript:
  "<text>" (<timestamp>)
<<<END BRIEF <12-hex>>>>
```

- `head` is `git -C <worktree> rev-parse HEAD` when the block is made.
- The template is filled with `string.Template.substitute`, variables
  `$persona`, `$worktree`, `$increment` (the path, or `none named` without
  `--increment`), `$note`; an unknown `$name` is an error at brief time, not
  a silent blank.
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

### 3.3 The template, and where the persona lines went

One file, `.claude/briefs/common.md`, governed (§5). It holds only what the
spawn alone knows or what has no other home: which files, which worktree,
which note file, whose words are Ola's, and the handback's shape. Concurrency,
the write limit, the increment's own lines and Ola's quotations are generated
(§3.1, §3.1a). Rules a persona file already states are not repeated, so a
persona-specific template would be empty; there are none. This departs from
the plan's "six templates" in the direction Ola asked for: less rule text.

The text, verbatim; tests pin the phrases in bold, not the wording around
them:

```
You are @$persona. This block comes from tools/brief.py; the task after it
is the main session's wording.
Read **.claude/agents/$persona.md from disk** and the increment file,
$increment. Where the task contradicts either, **the files win**, and you
say so; **a brief cannot drop a step** they require.
Work only in $worktree (cd there first), with your own build directory and
venv. Your note file is $note; write **no other file under
.claude/current-task/**.
**Ola's words appear only under "Ola, verbatim"** below.
If **blocked on power, network or a lock**, stop and hand back.
End each commit message with the **Co-Authored-By trailer** from your
system context.
Write in **plain words**: say what any internal label means.
Hand back under: **Result; Pinned or assumed beyond the design; Questions
for Ola; Lessons; ASK OLA and GUARD FALSE POSITIVE lines** ("none" under an
empty one). Each question for Ola is in plain words, with a default.
```

The trailer is carried by reference to the persona's own system context
(this @architect run received it that way), so a model change leaves no
stale literal in a governed file.

Where each recurring brief line now lives:

| Recurring line | Where | Already there, or moved by h9 |
|---|---|---|
| at most two writers; read-only agents beside them; nothing beside @perf | `brief.py`, §3.1a | new, by Ola's ruling |
| stop and report if blocked on power, network or a lock | `common.md` | from the memory note "agents report blocks" |
| the attribution trailer | `common.md` | by reference |
| Ola's words pasted verbatim | `--ola` and `common.md` | new |
| red-step choices go to @architect before green | `tester.md` §1, new bullet (below); the handback heading | moved |
| @reviewer is read-only; the spawner records | `reviewer.md` §5 | already there (8b8eb94) |
| @reviewer reports the commit range | `reviewer.md` §5, feedback item 2 (below) | moved |
| a performance fix is timed by @perf before review | `CLAUDE.md` §3, "Step order" (below) | moved: it is the dispatcher's order, not a persona's |
| plain language | `common.md` | from the memory note "no implicit labels", for personas |
| "Briefs" bullet, @tester part (§3A, §3C, §3D) | `tester.md` §3A ("every geometry test suite"), §3C ("Only for an increment that reads external input"), §3D ("every property test of refinement output") | already there; dropped |
| "Briefs" bullet, @developer part (minimal code, ceiling) | `developer.md`, first and fourth bullets | already there; dropped |
| "Briefs" bullet, @reviewer part (§5's checks, green CI is not done) | `reviewer.md` §5 | already there; dropped |
| mutation only for an invariant-critical suite ("lean briefs") | `docs/increments/README.md`, *Cost constraints*; `reviewer.md` §5 check 4 (8b8eb94); the quoted increment lines | already there; dropped |
| create the product file first | `REQUIRED-READING.md` | already there; dropped |
| never push | `REQUIRED-READING.md`, *Before you publish* | already there; dropped |

The persona-file lines, verbatim:

`tester.md` §1, a new last bullet:

```
* **Choices beyond the design:** list every choice your tests pin that the
  increment file leaves open, under the handback heading "Pinned or assumed
  beyond the design"; `@architect` confirms or rules on each before green.
```

`reviewer.md` §5, feedback item 2, from `**Size Metrics:** Confirm total LOC
and focus area.` to:

```
2. **Size Metrics:** The commit range reviewed, total LOC and focus area.
```

`CLAUDE.md` §3, "Step order", after "when refine or mesh code is touched
(`docs/increments/README.md`)." insert:

```
A performance fix is timed by `@perf` before review.
```

These edits are written against `worktree-agent-lines` (8b8eb94, which
changes `reviewer.md` §5 and `tester.md` §1). h9's PR is opened after that
branch merges, so the moved lines land on its text.

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
other subagent types (`general-purpose`, `Explore`, forks): they are not
personas, and h6 gives any name outside the persona table no write rights, so
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

Word counts are `len(text.split())`, Markdown emphasis removed, at 2060f14
for the repository and on 2026-10-04 for the memory notes.

| Passage | Change | Removed | Added |
|---|---|---|---|
| `CLAUDE.md` §3, "Briefs" bullet, from `* **Briefs.**` to `since green CI is not done.` | deleted | 55 | 0 |
| `CLAUDE.md` §3, "Step order" | the performance-fix sentence (§3.3) | 0 | 9 |
| `REQUIRED-READING.md`, *The harness*, after "configuration changes while unattended mode is on;" | the `guard_spawn.py` clause below | 0 | 37 |
| `REQUIRED-READING.md`, current-task bullet: "The spawner names the path in the prompt, and the subagent writes that path and no other." | becomes "`tools/brief.py` names the path, and the subagent writes no other." | 17 | 10 |
| `REQUIRED-READING.md`, restart bullet: "; until then its brief says to read the persona file from disk" | deleted; the block always says so (the bullet ends "...the changed persona**.") | 12 | 0 |
| `tester.md` §1 | the "Choices beyond the design" bullet (§3.3) | 0 | 35 |
| `reviewer.md` §5, item 2 | "Confirm" becomes "The commit range reviewed," | 1 | 4 |
| `.claude/briefs/common.md` | new (§3.3) | 0 | 150 |
| **Repository total** | | **85** | **245** |

**Net in the repository: +160 words.** The whole fixed part of a brief is
150 of them; nothing else grows by more than a sentence.

The `REQUIRED-READING.md` clause:

```
`guard_spawn.py` refuses a persona spawn without an unchanged, current block
from `python3 tools/brief.py`, and any spawn or resume made outside the
directory the session started in; a resume that starts a new step carries a
fresh block;
```

**Outside the repository**, loaded into every main session (the memory
folder `~/.claude/projects/-Users-skavhaug-projects-rasputin/memory/`):

| Note | Words | Why it can go |
|---|---|---|
| `lean-agent-briefs.md` | 185 | mutation only where the increment file names it: README, *Cost constraints*, quoted into each block; its postscript ("no PNGs unasked") is `CLAUDE.md` §3's "No unasked images"; "check a step past about 20 minutes" is the note "agents report blocks" |
| `agent-concurrency-pairs.md` | 296 | `brief.py`'s rules (§3.1a) and the template's "own build directory and venv". Its one sentence with no other home, "If a 'Usage limit reached' message appears again, go back to one agent at a time", moves into the note `agents-report-blocks.md` |

With both removed, what a main session loads shrinks by 481 words less one
moved sentence, so the rule text it carries goes down by about 300 words
overall. The main session removes both notes and their `MEMORY.md` lines,
after the merge and the first live check (§9): memory is its own, and no
persona writes there.

## 5. Files

| File | Change | Production lines (`CLAUDE.md` §2) |
|---|---|---|
| `tools/brief.py` | new, §3.1, §3.1a, §3.2 | ~140 |
| `.claude/hooks/guard_spawn.py` | new, §3.5 to §3.8 | ~50 |
| `.claude/settings.json` | §6 | 6 |
| `.claude/hooks/guard_governance.py` | `"tools/brief.py"` in `GOVERNED` (the hook imports it, so it is live before review, like the self-protecting set); `".claude/briefs/"` in `GOVERNED_PREFIXES` | 2 |
| `tools/rule_sizes.py` | `".claude/briefs/common.md"` appended to `RULE_FILES` | 1 |
| `.claude/briefs/common.md` | new, §3.3 | prose |
| `CLAUDE.md`, `.claude/REQUIRED-READING.md`, `.claude/agents/tester.md`, `.claude/agents/reviewer.md` | §3.3, §4 | prose |
| `tests/python/harness_fixtures.py` | `tools/brief.py` and `.claude/hooks/guard_spawn.py` in `COPIED` (the template is copied by `brief_fixtures.make_brief_repo`, §11 point 2) | test |
| **Total** | | **~200**, under the 700 ceiling |

The plan said about 140; the difference is the block reader beside its
writer, the transcript check of `--ola`, and the concurrency rules.

One template means one more row in the size table, about 50 characters;
h8's test 11 holds the recap's new sections to 3,000 characters, and its
design review measured about 280 of slack.

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
2. **Template.** Every persona's block holds the bold phrases of
   `common.md` (§3.3) with `$persona`, `$worktree`, `$increment` and `$note`
   filled in; without `--increment`, `none named`. An unknown `$name` in a
   fixture template: non-zero exit, no block.
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
6. **Concurrency (§3.1a).** Each with its own worktrees unless said:
   - `--beside none`: "you run alone"; no `--beside`: exit 2; `none` with
     another entry: exit 2; an entry without `:<worktree>`, or with a
     worktree that does not exist: exit 2.
   - `developer --beside tester:A`: names @tester as writer in A, "may
     build"; with `--no-build`, "may not build".
   - Read-only beyond two writers: `reviewer --beside tester:A --beside
     developer:B`, and `developer --beside tester:A --beside reviewer:C`:
     accepted; the line names @reviewer read-only and says read-only agents
     do not count toward the two. Two read-only entries beside two writers:
     accepted.
   - Rule 1: `perf --beside tester:A`, `tester --beside perf:A`, `perf
     --beside reviewer:C`, `reviewer --beside perf:A`: exit 2 naming the
     timing run.
   - Rule 2: `architect --beside tester:A --beside developer:B`: exit 2,
     "at most two writers".
   - Rule 3: `developer` in A `--beside tester:A`: exit 2.
   - Rule 4: `reviewer` in A `--beside developer:A`, and `developer` in A
     `--beside reviewer:A`: exit 2 naming the shared worktree.
   - Read-only is derived: with a fixture `WRITES` where `orchestrator` is
     empty, `orchestrator` counts as read-only.
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
    fresh block: no output; with an edited or stale block: deny; with a
    malformed marker (a `<<<BRIEF` line that does not match the header
    pattern, or an END line with no header): deny, `brief: malformed`.
19. Rule 1: `agent_id` present with a wrong `cwd` and no block: no output.
    `tool_name: Bash`: no output. Stdin not JSON, or a JSON list: no output,
    exit 0.
20. Unavailable (§3.7): `tools/brief.py` removed from the fixture, or made to
    raise on import: a blockless tester spawn passes with `additionalContext`
    containing `brief check unavailable`; a wrong `cwd` is still denied.
    The same when `brief.check` raises (a fixture `brief.py` whose `check`
    raises `RuntimeError`): pass, `additionalContext` names the error, exit
    0.
21. Unattended: with the fixture's unattended flag on, test 11's verdicts are
    the same, and no queue file appears under the fixture's harness
    directory.
22. Missing `cwd` in the event, or `CLAUDE_PROJECT_DIR` unset: the
    comparison is skipped, said in `additionalContext`; the block checks
    still run.

Wiring and governance:

23. `test_guard_governance.py`: writes to `tools/brief.py`,
    `.claude/briefs/common.md` and any new `.claude/briefs/x.md` get `ask`.
24. `test_settings_wiring.py`: a `PreToolUse` entry with matcher
    `Agent|SendMessage` and command
    `$CLAUDE_PROJECT_DIR/.claude/hooks/guard_spawn.py`; the file is
    executable on disk and in git's index (mode 100755); run by path from a
    fixture copy it denies a blockless tester spawn. The existing entries
    unchanged.
25. `test_rule_sizes.py`: the table has a `.claude/briefs/common.md` row
    after `docs/PRINCIPLES.md`; at a reference without it the row says `new`.
26. Rule text: `CLAUDE.md` has no line starting `* **Briefs.**`, and its
    "Step order" bullet holds "timed by `@perf` before review";
    `REQUIRED-READING.md` names `guard_spawn.py` and `tools/brief.py` in *The
    harness*, and no longer holds "until then its brief says";
    `tester.md` holds "Choices beyond the design"; `reviewer.md` holds "The
    commit range reviewed".

## 8. Not covered, named

- A resume that starts a new step without a block (§3.8).
- Spawns of types that are not personas (`general-purpose`, `Explore`,
  forks) until h6's role table gives them no write rights.
- The concurrency rules trust the `--beside` list: a persona left off it is
  not counted, and nothing checks that the listed ones are still running
  (stage 4 records starts and stops).
- A read-only persona that reads files outside its own worktree (master, a
  writer's worktree): "the same files" is measured by worktree only.
- The task text after the block: it can still paraphrase, contradict or
  drop things; the block's sentence "the files win" is the persona's
  defence, not a check.
- A block computed by hand rather than printed (§3.4).
- The next free increment id (stage 2).

## 9. The first live firing

After the merge, by day, with Ola present, the main session checks the hook
does what the tests say, against the real tool input: (a) a @tester spawn
with no block is refused with `brief: no block`; (b) after `cd` into a
worktree, a `SendMessage` is refused with `cwd:`; then `cd` back; (c) a spawn
with a fresh block goes through, and the persona's handback has the five
headings. If (a) or (b) passes through, the tool input fields of §3.5 differ
from the transcripts', and the hook is reported to Ola before any other
spawn. Then the main session removes the two memory notes (§4).

## 10. Questions for Ola

None. Two choices made here, which Ola may overturn:

1. **Claude's own helper agents need no brief block.** Claude Code has
   built-in helpers besides our six personas: a general one, and a
   read-only "Explore" one for searching code. The check applies only to
   the six personas. Consequence: the main session can still hand work to
   a helper with a hand-written prompt and nothing checks it. Once h6 lands,
   such a helper cannot write files in the project, so it can search but
   not do a persona's work. The other choice, requiring a block for helpers
   too, would need a template for them, and a quick search would cost a
   brief.
2. **If `tools/brief.py` breaks, spawns go through unchecked, with a
   notice.** Consequence: until it is repaired, a brief can again be
   missing its fixed part, and the main session is told so on every spawn.
   The other choice, refusing every spawn, would also refuse the
   `@developer` who repairs `brief.py`, and only Ola could unblock it by
   hand. The check on the working directory does not depend on `brief.py`
   and keeps working either way.

## 11. Settled after h9's red step (84fbde6)

@tester pinned ten choices the design left open. Rule numbers 1 to 4 in
points 5 and 8 are §3.1a's concurrency rules, not §3.5's hook rules.

1. **`brief.main(argv: list[str] | None = None) -> int`.** Confirmed; it
   returns the exit status (0, or 2 on a refusal) and does not call
   `sys.exit`, so test 6 can call it after emptying `WRITES` in the fixture
   copy.
2. **The template is not in `COPIED`.** Confirmed, and §5 is overruled on
   this point: `COPIED` gets `tools/brief.py` and
   `.claude/hooks/guard_spawn.py` only. Copying `common.md` into every
   harness fixture would move `test_rule_sizes.py`'s pinned totals;
   `brief_fixtures.make_brief_repo` copies it where it is needed.
3. **A pass with a notice.** Ruled: **no `permissionDecision` key at all**,
   only `hookSpecificOutput.additionalContext` (with `hookEventName`). Not
   `"allow"`: this hook only refuses or stays silent, and never grants past
   the permission system. Rule 22's notice mentions `cwd` or `directory`:
   confirmed.
4. **The rule-2 refusal** contains the event's `cwd`, the project path, and
   exactly `run: cd <project>`. Confirmed.
5. **Wordings pinned only where the design gives them**: concurrency rule 1
   "timing run", rule 2 "at most two writers", rule 4 the shared worktree's
   path, a missing `--ola` quotation its own text; rule 3 by exit 2 only.
   Confirmed.
6. **Test 8.** The write-limit line ends `and your note file.`; `WRITES` is
   compared as text; reviewer's line names no other persona's path.
   Confirmed.
7. **Test 3.** The `last:` verdict is cut to exactly 200 characters, a
   prefix of the verdict line, with no marker. Confirmed (h8's display cut
   adds `...`; this one does not, because the persona reads the file).
8. **Also refused (exit 2):** `--beside nobody:<path>` (not a persona), and
   a `--worktree` that is a subdirectory of a checkout (§3.1 asks that
   `--show-toplevel` resolves to the path itself). Confirmed.
9. **Test 21.** After the hook, the fixture's harness directory holds only
   `unattended.json`. Confirmed: nothing queued (§3.5).
10. **Test 24.** The new entry directly after the `AskUserQuestion` entry.
    Confirmed: that is §6.

**What @tester amends before green** (one commit, reason in its message):
point 3, `tests/python/test_guard_spawn.py` lines 342 and 380,
`assert decision in (None, "allow")` becomes `assert decision is None`.

## Review


**Design review, round 1, 2026-10-04.** Range `2060f14..7cae124`. Verdict: CHANGES REQUESTED. LOC: 0 (design); estimate about 200, plausible. Every factual claim checked holds (`settings.json` lines 23-29, `guard_unattended.py`'s deny format and crash refusal, `session_state.py`'s helpers, the 55-word bullet, `WRITES` against h6's table, the tool-input fields, the quoted Claude Code docs); each hook rule has a test. Blocking: (1) concurrency cannot express Ola's ruling of 2026-10-04 (read-only agents may run beyond two writers if nothing they read is being changed): `--beside` takes one name; make it repeatable, mark read-only (or derive it from an empty `WRITES`), say so in the generated line, keep "nothing beside @perf" and "at most two writers", extend test 6; (2) the increment adds about 600 words of rule text (826 in templates plus 47) against 55 deleted, and §4's net is misstated; cut the templates to what only the spawn knows, or move persona lines into the persona files and delete what they repeat, and state the real net (candidates: the pairs note once (1) lands, REQUIRED-READING's current-task paragraph). Suggestions: state the two defaults for Ola in plain words with their consequences; a test where `check()` raises; a malformed block in a `SendMessage`. Not pushed; no CI.

**Design review, round 2, 2026-10-04.** Range `1c95f7a..9740246`. Verdict: APPROVED. LOC: 0 (design); estimate about 200. Both blockers and the three suggestions closed: §3.1a quotes Ola's ruling and has four exit-2 rules, each tested; the word table recounts exactly (+160 in the repository, about −320 with the two memory notes); every "already stated elsewhere" claim checks; the PR opens after `worktree-agent-lines`. Suggestions: three dropped template lines have no home (@orchestrator measures rule text and proposes a cut; @orchestrator quotes Ola only from the transcript; @reviewer says whether @perf's acceptance is recorded); list the usage-limit sentence's move in §9. Not pushed; no CI.

**Code review, round 1, 2026-10-04.** Range `2060f14..661ef55` (design 7cae124 and 9740246, red 84fbde6, rulings fac025a, amendment 1922d72, green 661ef55). Verdict: CHANGES REQUESTED. LOC: 300 (`brief.py` 218, `guard_spawn.py` 74, `settings.json` 6, `guard_governance.py` 1, `rule_sizes.py` 1) against about 200; the overrun is unpriced code (`Block`, `MalformedError`, the argparse subclass, review-section parsing, §3.1a's rules), not formatting; ten `# fmt: skip` regions stay readable. `settings.json` changes exactly as §6; the rule text matches §3.3 and §4 word for word; the refusal order matches §3.5; fail-open matches §10. Blocking: (1) the status line still says "design, round 2", and the measured 300 lines are unrecorded; (2) `--increment` is resolved against the checkout running `brief.py` (`brief.py:241`), so a branch-only increment file is refused and a branch with newer rulings gets master's stale quotes; resolve a relative path against `--worktree`, ruled in §3.1, with a test whose worktree copy differs. Suggestions: `brief.py:181` cites §3.1a for §3.1; align §5's `COPIED` row with §11 point 2; decide the three homeless template lines. Not pushed; no CI.
