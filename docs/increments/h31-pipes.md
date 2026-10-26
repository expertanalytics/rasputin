# Harness h31: development pipes

Status: **design to build, first slice**, ruled by Ola 2026-10-25 ("yes, build
the first slice"): rules option A, the rest as recommended. Section 8 holds
the rulings. Design round 3 by `@architect`, 2026-10-25, on
`worktree-pipes-design`; round 1's seven review findings are closed in
sections 3 to 9, each marked *(review 1, item N)*; round 2's three in
sections 4 and 11, marked *(review 2, item N)*, its wordings in 4, 7, 9, 11.
Amended after the red step (`5fc19f4a`), 2026-10-26: the tester's pins ruled
(section 11) and the recap block made one table (section 7, Ola: "We need a
nice way to look at the pipe, probably in a table format.").

Ola, 2026-10-25: "When I'm back, we should make a few development pipes. Then
we add tasks to these pipes. We are limited by the number of agents, but this
way we don't get things left behind." And: "One thing - when stuff is stuck
behind for a while, they should resurface." Settled the same day: a plugin
system, with a core pipe in this repository and further pipes added by the
private repository; priorities set by Ola as he goes, with a fallback so no
agent sits idle; rules stay global for now (section 6, option A).

Words used here. **Pipe**: a named, ordered list of items. **Item**: one
PR-sized piece of work with a state. **Writer**: an agent that changes files
(at most two run at once, `CLAUDE.md` §3). **Runnable**: an item whose next
step needs no answer from Ola and, while away, writes no governed path.
**Stuck**: an item that has not changed state for more than three days
(section 9). **In flight**: an item in state `designing`, `red`, `green` or
`review`; `queued` and `waiting` items are not in flight.

## 1. The problem, measured

The dispatcher's queue is `.claude/current-task/session.md`. Its format rule
(`.claude/REQUIRED-READING.md`, *While you act: in-flight state goes on disk*)
allows one
`NOW:` line, one `QUEUE:` line and one `ASK OLA:` line per decision, each
under 300 characters. Measured in the main checkout at 18:3x UTC today:

| Measure | Rule | Today |
|---|---|---|
| Lines, not blank | no number stated | 71 |
| `QUEUE:` lines | 1 | 19 |
| `ASK OLA:` lines | one per open decision | 44, some answered on the `RULINGS` line |
| Lines of other kinds (`RULINGS`, `PUSHED`, `FALLBACK`, `QUEUED`) | 0 | 6 |
| Lines over 300 characters | 0 | 22 |
| Recap size (`python3 tools/session_state.py \| wc -m`) | 10,000 cap | 33,238 |

The queue is a prose log, so an item is "left behind" when its line scrolls
out of the main session's attention; whether it is behind is not visible, and
"which pipe gets the next free writer" is answered from memory each time.
Age is invisible too: the same file at 19:2x UTC has 21 `QUEUE:` lines whose
dated mentions (`grep -o "10-[0-9][0-9]"`) run from 10-06 to 10-25, so items
from nine, five and four days ago sit beside today's with nothing marking them.

## 2. Prior art: legacy and literature

*Legacy.* `git grep -n -i "pipe\b\|lane" legacy-archive -- legacy` returns two
lines, both unrelated (`compute_slopes`, triangle orientation); nothing is
carried across. The harness's own legacy is the dispatcher-control plan,
`docs/research/2026-10-03-dispatcher-control.md`, stage 2 (lines 312-322): a
`tools/session.py` writing the queue in a fixed format and a
`tools/pipeline.py` printing "one line per open branch saying which step is
next". Ola ruled yes to stage 2 on 2026-10-04 (`docs/retrospectives/next.md`,
*Taking recurring decisions off the main session*). Neither file exists
(`ls tools/pipeline.py tools/session.py`: no such file). This note is stage 2
with the parts Ola asked for today added: several lists, priorities between
them, a fallback, plugin pipes, and resurfacing.

*Literature.* Tooling; no novelty claimed. Four borrowed ideas: Kanban's
work-in-progress limit and swimlanes (Anderson, *Kanban*, 2010: a lane per
class of work, a cap on items in progress; the lean slogan is to stop starting
and start finishing),
which is what `CLAUDE.md` §3 *Finish before starting* already states for one
lane; a work-conserving scheduler, which never leaves a resource idle while
any queue has ready work (Kleinrock's term, *Queueing Systems* vol. 2, 1976),
which is Ola's fallback; strict priority with a floor against starvation,
as in Linux's scheduling classes, so the lowest pipe still moves when the
higher ones are blocked; and ageing, the classic cure for starvation under
strict priority (Silberschatz, Galvin and Gagne, *Operating System Concepts*,
the *Priority Scheduling* section: a waiting process's priority rises with its
wait), which is Ola's resurfacing, done here as one threshold rather than a
continuous ramp because a recap reader needs a short list, not a reordering
he has to notice. The plugin side copies this repository's own
`docs/plugins.md`: a plugin adds commands by declaring a file the core
discovers; the core never names the plugin.

Ideas-file entries this note takes up or overlaps (`../rasputin_scratch/ideas/ideas.md`):
*Recap: items in flight and open ideas* (the recap shows only the current
branch; this slice prints every pipe's items with their states, which is that
entry's "In flight: N of 3"); *Lead time from design to merge, per increment*
(item states with times give it without a second tool; the report itself is
not in this slice); *Morning check proposals*, P8
(`session.md`'s format enforced: the format shrinks instead, section 7).

## 3. What a pipe is

A pipe is two things, kept apart because they change at different speeds:

- Its **definition**, tracked: `.claude/pipes/<name>.toml`, a few lines: the
  name and one sentence of purpose. It changes through a PR, like a rule. The
  second slice may add checks (section 6, option B); the first slice reads
  only `name` and `purpose` and warns on any other key, so a definition
  written for the second slice is not silently half-read.
- Its **live list**, untracked: `<git-common-dir>/pipes/<name>.jsonl`
  *(review 1, item 1: not under `harness/`, which no agent may write;
  `.claude/REQUIRED-READING.md`, *The harness*)*. One JSON object per line,
  the newest line for an item id being its current record: id, title, state,
  pipe, `next` (the next step, free text), branch, PR number, the open
  question and its default if any, `plain` (writes no governed path, section
  4), `fallback` (section 4), `since` (ISO UTC time of the last state change)
  and `priority_words` on the pipe's own record. Written only by
  `tools/work_pipes.py`, appended, never rewritten, so a crash leaves a readable
  file and `git` is not involved. Checked today: `guard_governance.judge_bash`
  passes `python3 tools/work_pipes.py set h31 red`, `... add core "x" --plain`,
  `... ask h31 "q?" --default "yes"` and a shell append to
  `.git/pipes/core.jsonl` (returns `None` for each), and denies an append to
  `.git/harness/pipes/core.jsonl`; so no guard entry is needed, and the
  harness directory stays the hooks' and `away.py`'s alone.

States: `queued`, `designing`, `red`, `green`, `review`, `waiting` (on Ola),
`merged`, plus `parked` and `dropped` (Ola's words, kept so the item is not
re-proposed). Each change is one command, `python3 tools/work_pipes.py set <item>
<state> [--branch B] [--pr N] [--next "<step>"]`, by the main session; the
recap prints one table, one row per item (section 7). Where
git can tell the state (a `## Review` record, a merge on master), the tool
prints a warning when the file disagrees, and nothing more in the first slice.

**How the private repository registers pipes without this repository naming
them** *(review 1, item 3)*. The private remote is a second remote of the
same checkout, and its `master` holds commits this repository's `master`
does not (`tools/push_grant.py`, `home`). A private pipe definition is a file
`.claude/pipes/<name>.toml` in one of those commits. `work_pipes.py` reads
definitions from the working tree, then for every remote `git remote` lists
other than `origin`:

- `git rev-parse --verify -q <remote>/master`; when that fails (exit 128,
  checked today with `git ls-tree nosuch/master`), it prints one warning,
  `pipes: no <remote>/master locally; git fetch <remote> to read its pipes`,
  and goes on. **No fetch**: the recap runs inside a hook, and a hook must
  not touch the network.
- `git ls-tree --name-only <remote>/master .claude/pipes/`, then `git show
  <remote>/master:<path>` per file. A name already read from the working
  tree is not read again (the private `master` contains the public one, so
  `core.toml` is on both). The block header prints what was read:
  `pipes read from private/master 79848102, committed 2026-10-25` (the sha
  from `git rev-parse --short`, the date from `git log -1 --format=%cs`).
  Checked today in the main checkout: `private/master` is `79848102`,
  committed 2026-10-25, and `ls-tree` of `.claude/pipes/` there prints
  nothing with exit 0, which the tool reports as `no pipe definitions`.

Live lists live in the git common directory, outside both repositories'
trees, so a private item's title is on nobody's `master`.

*Alternatives* (unchanged from round 1): (a) one tracked file for all pipes:
every state change is a commit, conflicts between worktrees, and private names
would have to be in it. (b) Live lists tracked too: the same conflicts; the
private list would need a private branch per change. (c) No files, states
derived from git and increment files (stage 2's `pipeline.py` alone): clean
for items in flight, but a queued item has no git trace, which is exactly the
item that gets left behind. (d) Keep `session.md`, enforce its format: the
format was ruled 2026-10-03 and today's file breaks it on 37 of 71 lines; a
format nothing reads is not enforced by a warning.

## 4. Priorities, the writer rule, and the fallback

**Setting priorities.** Ola says it in chat, in plain words ("core first, then
the private work"; "put the showcase on top"). The main session turns it into
one command, `python3 tools/work_pipes.py priority core=1 <name>=2 ... --words
"<Ola's words>"`, which stores an integer per pipe (1 highest) with the time
and Ola's words; `work_pipes.py top <item>` moves an item to the head of its pipe.
`work_pipes.py next` prints which item would get the next free writer and why, so
Ola sees the effect before an agent is spawned. No file edits by hand, no
numbers Ola must remember: the recap shows the current order.

**Adding an item.** `python3 tools/work_pipes.py add <pipe> <id> "<title>"
[--plain] [--fallback] [--next "<step>"]`. The id is a short slug, unique
across pipes (the tool refuses a duplicate). `--plain` records that the item
writes no governed path *(review 1, item 4: the field the away rule needs,
set when the item is added, by the main session's judgement of the item's
files; `work_pipes.py set` may change it later with `--plain` / `--governed`)*.
`--fallback` marks the pipe's standing fallback (below).

**What a free writer takes.** `work_pipes.py next` applies, in this order, the
first that holds; "runnable" excludes `waiting` items always, and while
unattended (`harness_mode.read_mode`) also items without `plain`. An
in-flight item is **busy** while `.claude/current-task/` holds a note `<persona>-<name>-<HHMMSS>.md`, the name being everything between the first hyphen and the six-digit time, for its branch's worktree, the branch being
`worktree-<name>` (`tools/brief.py` writes the note, `tools/new_worktree.py`
names the branch): a writer already has it *(review 2, item 1)*.

1. Across all pipes, the runnable in-flight item that is not busy, nearest to
   merging (`review` before `green` before `red` before `designing`); this
   keeps `CLAUDE.md` §3's rule and Kanban's "finish first". Ties go to the
   higher-priority pipe, then to the older `since`.
2. If the in-flight count is under three: a stuck runnable `queued` item,
   oldest first, from any pipe *(section 9)*; otherwise the highest-priority
   pipe's first runnable `queued` item. Steps 2 and 3 skip fallback items
   *(review 2, item 2)*.
3. The next pipe by priority with a runnable `queued` item, under the same
   count.
4. Under the same count *(review 2, item 3)*, any pipe's fallback item,
   highest priority first, when it is runnable and neither in flight nor
   busy.
5. Idle, printed as `IDLE: no runnable item in any pipe`, which the morning
   check counts.

The output is one line, `next: <pipe>/<id> <state> (<reason>)`, where the
reason names the step that chose it ("nearest to merging", "stuck 4 d",
"core is first", "fallback"). `next` also prints one `busy: <pipe>/<id> (<note file>)` line for each in-flight item it skipped as busy.

**Counting.** The cap stays three items in flight across all pipes
(`CLAUDE.md` §3). A `waiting` item is not in flight *(review 1, item 4; this
is `CLAUDE.md` §3's own "except while it waits only on Ola")*, nor is a
`queued` one. To stop one pipe from holding all three while another has
ready work, a pipe holds at most two in flight unless it is the only pipe
with runnable items. *Alternative:* weighted shares (pipe A gets two writers
of three); more rule text, and priorities already express it.

**The fallback** *(review 1, item 5)*. A fallback is an item in a pipe's live
list marked `--fallback` at `add`: work that needs no ruling and writes no
governed path (so `--fallback` implies `--plain`, and the tool refuses
`--fallback --governed`). Step 4 takes it only when steps 1 to 3 find
nothing, which is the case Ola described: every top item waits on him. It is
a live-list fact, not a tracked field, so using it or naming a new one is a
command, never a commit. Today's `FALLBACK:` line, the basin-scale pieces, is
the core pipe's first fallback item. The sentence in
`.claude/REQUIRED-READING.md` that asks the `QUEUE:` line to name a fallback
before Ola leaves **stays in this slice**; the second slice makes `away.py`'s
enter step read `work_pipes.py next` instead and refuse, with the reason, when no
pipe has a fallback, and only then is that sentence replaced and its refusal
test written.

## 5. Questions to Ola, attached to items

`python3 tools/work_pipes.py ask <item> "<question>" --default "<what is taken if
silent>"` puts the question on the item and sets it `waiting`; only that item
stops, and the rule above picks the next one. `work_pipes.py answer <item>
"<Ola's words>"` records the ruling on the item's record and returns it to the
state it had before `ask`. The recap lists open item questions grouped by
pipe, numbered in one sequence, each with its default, which is the form Ola
asked for; stuck ones come first (section 9). An answered question is no
longer printed (it stays in the appended history).

**What stays an `ASK OLA:` line** *(review 1, item 2)*. The `ASK OLA:` lines
in `session.md` are read by `tools/harness_mode.py` (the refusal text tells
the agent to write one), `.claude/hooks/guard_unattended.py` (the same),
`tools/away.py --back` (prints every worktree's) and
`tools/session_state.py` (*Waiting on Ola*), checked with
`grep -n "ASK OLA" tools/*.py .claude/hooks/*.py`. So: **refusals and
defaults taken while away stay `ASK OLA:` lines in `session.md`, in the form
`.claude/REQUIRED-READING.md` gives; only questions about an item move to
`work_pipes.py ask`.** A question about the harness or process, not about an item,
goes on a pipe instead: `work_pipes.py ask core "<question>" --default "..."`
puts it on the pipe's own record (the one that holds its priority words),
where it is listed with the item questions and answered the same way; no
item changes state. Today's 44 `ASK OLA:` lines become: refusals and
defaults as they are, and at most one open question per waiting item or pipe.

*Alternative:* keep every question as an `ASK OLA:` line prefixed with the
item id. Cheaper (no tool), but it is what exists today, and section 1 shows
the result.

## 6. Rules in pipes

Ola: "I think rules could be a useful concept in pipes. But we need to discuss
it more." Ruled 2026-10-25 with the first slice: **option A now**; B after a
week of use; C is not offered. The cost is measured in rule words, which
`CLAUDE.md` §2 says must be paid for.

**A. Rules stay global; pipes carry tasks only** (this slice). A pipe
definition has no rule field. Cost: 0 words. The private repository keeps its
own rule copies as today, and the drift between them (12,080 words against
9,494 public, per an `ASK OLA:` line in today's `session.md`) stays a
separate question.

**B. A pipe adds checks from a fixed menu** (second slice, after a week on
A). The definition may name: extra skills every persona loads; an
acceptance command run before review (as `tools/bench.py` is for refine and
mesh code); a reviewer focus line; a size budget below 700 lines. `brief.py
--pipe <name>` appends these to the brief as lines the hash covers. The menu
is closed, so a pipe can tighten but never loosen a rule. Cost: about 50
words in `docs/increments/README.md` (the menu and that a pipe only
tightens) and about 30 lines in `brief.py`; pipe definitions are then counted
by `tools/rule_sizes.py` so growth there shows in every recap. This is the
plugin hook Ola asked for, in the same shape as `docs/plugins.md`: a file the
core discovers, a closed interface.

**C. Free-text rules per pipe,** a Markdown file the brief includes. The most
flexible, and the thing §2 fights: unbounded, unreviewed by this repository's
`@reviewer`, and a second place for every persona to read before acting.
Rejected.

## 7. Blueprint: modules, data flow, and what the recap prints

`tools/work_pipes.py` is standard-library only (`json`, `tomllib`, `argparse`,
`dataclasses`, `datetime`, `subprocess` for the five git calls), like the
other harness tools, because the recap runs from a hook on the system
`python3` with no venv. Four layers, each testable without the one below:

1. **Records.** A frozen `Item` dataclass (the fields of section 3) and
   `Pipe` (name, purpose, priority, source: `tree` or `<remote>/master
   <sha7>`). Pure: `load_items(path) -> dict[str, Item]` folds the appended
   lines to the newest per id; `fold` is a pure function over a list of
   dicts, tested on strings.
2. **Definitions.** `definitions(repo) -> tuple[list[Pipe], list[str]]`: the
   working tree's `.claude/pipes/*.toml`, then each other remote as section 3
   says; the second value is the warnings. Git is reached through one
   injectable `run(args) -> (code, stdout)` so tests pass a fake.
3. **Rules.** `order(pipes, items, now, unattended, busy) -> Choice`
   (section 4's five steps, returning the item and the reason; `busy` is the
   set of branches `busy_branches(tasks) -> set[str]` lists from the note
   files in `.claude/current-task/`), `stuck(items, now) -> list[Item]`
   (section 9), `in_flight(items) -> int`. Pure functions of their
   arguments; `now`, `unattended` and `busy` are parameters, never read
   inside.
4. **Commands.** `add`, `set`, `ask`, `answer`, `priority`, `top`, `next`,
   `show`; each parses, reads, applies one pure function, appends one line.
   `show` prints the block below; `session_state.print_recap` calls
   `pipes.show_lines(repo, now)` and prints them after *Waiting on Ola*,
   before the `session.md format:` warnings.

**The recap block is one table**, the same text as `work_pipes.py show` and
`python3 tools/work_pipes.py` with no arguments; there is no second rendering.
Capped at 40 lines then `... N more`, each line at most 120 characters, so the
block is at most 4,800 of the recap's 10,000-character cap:

```
Pipes (in flight 2 of 3; pipes read from private/master 79848102, committed 2026-10-25)
     PIPE    PRI ITEM                          STEP      WAITING ON  SINCE      NEXT
     core    1   h30 governed files            waiting   Ola         4 d stuck  question 1 below
     core    1   basin basin-scale pieces      queued    free agent  9 d stuck  fallback: design
next core    1   tri-budget triangle budget    review    free agent  2 h        reviewer round 2 #351
     core    1   h31 development pipes         red       its agent   5 h        green
     core    1   fable-june Fable June data    queued    free agent  1 d        design
     lsm     2   snow-shelter wind shelter     queued    Ola present 2 d        design
Questions for Ola:
  1. core/h30: finish the 4 governed files with you at the prompt?  default: waits
next: core/tri-budget review (nearest to merging)
```

One row per item that is not `merged`, `parked` or `dropped`. Columns, fixed
widths, a cell cut with `...` when longer: the mark (5: `next` on the row
`order` picks, else blank), PIPE (8), PRI (4, the pipe's priority), ITEM (30,
id then title), STEP (10, the state), WAITING ON (12: `Ola` when `waiting`;
`Ola present` when unattended and not `plain`; `its agent` when busy; `free
agent` otherwise), SINCE (11: whole days, or hours under one day, then
`stuck` past section 9's threshold), NEXT (the rest, to 120: `next`, then
`#<pr>` if set; a `waiting` row says `question N below`). Rows: stuck first,
oldest first; then by pipe priority; within a pipe, in flight nearest to
merging first, then `waiting`, then `queued` in the pipe's order (`top`).
Plain ASCII, no colour. The header line carries the unreadable count when
not zero (`1 unreadable`) and the definition warnings follow it.

`session_format` *(review 1, item 2)* stops requiring exactly one `QUEUE:`
line: zero or one is accepted (one, because the fallback sentence of section
4 still names it in this slice), more than one is warned as today. `NOW:`
stays exactly one; `ASK OLA:` lines stay as they are.

## 8. Rulings (Ola, 2026-10-25, "yes, build the first slice")

The eight questions of round 1 and what the yes settles, each at its default:

1. Two files per pipe, a tracked definition and an untracked live list: yes.
2. Priorities as one integer per pipe set by `work_pipes.py priority`, from Ola's
   words in chat: yes; no hand edits.
3. The writer order of section 4 with at most two in flight per pipe: yes,
   both; stuck items ranked as section 9 says (Ola's resurfacing).
4. Questions attached to items, one `waiting` state, answered questions no
   longer printed: yes.
5. Rules: option A now, option B after a week of use.
6. Pipes to start with: `core` here; one private pipe, defined on the private
   `master` by a private PR after this one merges (not this slice's work).
7. Today's `session.md` moved into the pipes by hand on the first day,
   dropping the stale asks with Ola at the prompt: yes, the tool's first use.
8. `@orchestrator`'s post-merge check reading `work_pipes.py show` for lead time:
   later, not in this increment.

Open for a later ruling, not asked now: whether the stuck threshold of three
days (section 9) should differ between `queued` and `waiting`.

## 9. Resurfacing: stuck items *(review 1, item 6; Ola's requirement)*

An item that is `queued` or `waiting` and whose `since` is more than **three
days** before `now` is **stuck**. The constant assumes Ola's check cycle of
about one day (he drops in daily) and allows a weekend; it was checked
against today's `session.md` (section 1), where
items dated 10-16, 10-20 and 10-21 would be stuck and the 16 dated 10-25
would not. It is one module constant, `STUCK_AFTER = timedelta(days=3)`.

What resurfacing does, in both places an item can be overlooked:

- **The recap.** Stuck items are the table's first rows, oldest first, SINCE
  reading `4 d stuck` (section 7), and a stuck `waiting` item's question is
  numbered first under *Questions for Ola*.
- **The order.** In `next`, step 2 prefers a stuck item to a fresher one,
  oldest first, across pipes: a stuck `queued`
  item in the lowest pipe is taken before a fresh one in the highest. A stuck
  `waiting` item is not runnable, so it only surfaces; its ruling is Ola's.

Not done: an in-flight item (`designing` to `review`) unchanged for three days
is a stalled writer, not a forgotten item; step 1 already takes it first,
and the morning check reports idle writers. An item that resurfaces and is
still not wanted is `parked` or `dropped` by Ola's word, which ends it.

## 10. What changes, and the first slice

| Where | Change | Size |
|---|---|---|
| `tools/work_pipes.py` (new) | records, definitions (tree and other remotes, no fetch), the order rule, stuck, busy, `add`, `set`, `ask`, `answer`, `priority`, `top`, `next`, `show` | about 260 lines (220 + about 15 for stuck + about 10 for busy + about 15 for the table: about 25 for columns, row order and cell cutting, less about 10 for the per-pipe blocks and `Stuck:` section it replaces) |
| `tools/session_state.py` | print `pipes.show_lines` where `QUEUE:` was read; `session_format` accepts zero or one `QUEUE:` line | about +15, -20 |
| `.claude/hooks/guard_governance.py` | `"tools/work_pipes.py"` added to `GOVERNED`, because `session_state.py` imports it | +1 |
| `.claude/REQUIRED-READING.md` | the `session.md` format line: `NOW:` (one), `ASK OLA:` lines, `QUEUE:` (at most one); item questions go to `work_pipes.py ask`; the fallback sentence kept | net about -15 words |
| `CLAUDE.md` §3 | *Finish before starting*: "a free writer goes first ... process work" becomes "`python3 tools/work_pipes.py next` says which item" | net about -20 words |
| `.claude/pipes/core.toml` | `name = "core"`, `purpose = "..."` | 2 lines, not counted |
| `tests/python/harness_fixtures.py` | `tools/work_pipes.py` added to `COPIED` | uncounted |
| `tests/python/test_work_pipes.py` | section 11 | uncounted |

Net about 255 lines of counted code *(review 1: the round-1 table summed to
220, not 250; the table view adds about 15)*, one PR under the 700 cap. **Governed paths throughout**
(`session_state.py`, `guard_governance.py`, both rule files), so built by
day, with Ola at the prompt. **The second slice**: `away.py`'s enter step
reads `next` and refuses with no fallback (then the REQUIRED-READING fallback
sentence is replaced and its refusal test written), `brief.py --pipe`, and
option B.

Governed files after this slice: `tools/work_pipes.py` is in `guard_governance.py`'s
`GOVERNED` tuple, beside `tools/session_state.py`, for the reason the
`brief.py` entry gives there: a file the harness imports is harness code.
`tools/check_citations.py`'s own `GOVERNED` (rule files for the citation
check) is unchanged: `work_pipes.py` is code, not a rule file.

## 11. Tests (`tests/python/test_work_pipes.py`)

The module is loaded from its path as `test_session_state.py` loads
`session_state.py`; the pure functions get their `now` as an argument, so no
test sleeps or patches the clock. The fixture `now` is
`2026-10-25T19:00:00Z`; ages below are measured from it. The temporary
repository is `harness_fixtures.make_repo`, which copies `tools/work_pipes.py`
with the rest.

1. **Fold.** Three appended lines for `h31` (`queued` at 19:00 minus 2 d,
   `designing` at minus 1 d, `red` at minus 1 h) fold to one record in state
   `red` with `since` one hour before `now`; an unreadable line is counted,
   not raised.
2. **Order, step 1.** Two pipes, `core` priority 1 and `lsm` priority 2;
   items `a` `review` in `lsm`, `b` `green` in `core`, `c` `queued` in `core`:
   `next` is `lsm/a`, reason "nearest to merging".
3. **Order, steps 2 and 3, and the caps.** `core` has two in flight (`red`,
   `green`), both busy (the fixture's `.claude/current-task/` holds one
   `developer-<name>-190000.md` per branch, plus a `session.md` that
   `busy_branches` ignores), and one `queued`; `lsm` has none in flight and
   one `queued`: in flight is 2 of 3, `core` is at its per-pipe two, so
   `next` is `lsm`'s queued item. With the `red` item's note removed, `next`
   is it, reason "nearest to merging". With three in flight across the
   pipes, all busy, `next` is `IDLE`, never a `queued` or fallback item
   *(review 2, items 1 and 3)*. A note `developer-<name>-x-190000.md` does not make `<name>` busy.
4. **Waiting is not in flight.** Five items: three `waiting`, one `red` and
   busy, one `queued`; `in_flight` is 1, and `next` is the `queued` item
   *(review 2, item 1)*.
5. **Unattended.** The same items with `unattended=True`: a `queued` item
   without `plain` is skipped; with none `plain`, `next` prints
   `IDLE: no runnable item in any pipe`.
6. **Fallback.** Every non-fallback item `waiting`; one `queued --fallback`
   item in `lsm`: `next` is it, reason "fallback"; with a plain `queued`
   item added to `core`, `next` is that one, not the fallback *(review 2,
   item 2)*. `add --fallback --governed` is refused with a message naming
   both flags.
7. **Private definitions.** A second bare repository, added as remote
   `private`, with one commit holding `.claude/pipes/lsm.toml`; after
   `git fetch private` (in the test, not in the tool), `definitions` returns
   `core` from the tree and `lsm` from `private/master <sha7>`, with the
   header line naming that sha and its commit date. With the remote added but
   never fetched, the result is `core` alone plus the one warning
   `no private/master locally`, exit 0. The fake `run` records that no
   argument list starts with `fetch`.
8. **Stuck.** Items `old` `queued` with `since` 4 d before `now`, `fresh`
   `queued` 1 d before, `q` `waiting` 5 d before; `stuck` returns `q` then
   `old` (oldest first); in `show_lines` the first two table rows are `core`
   `q` (SINCE `5 d stuck`) and `lsm` `old` (`4 d stuck`, marked `next`), and
   the `fresh` row (`1 d`, no `stuck`) comes after both *(amended after red:
   no `Stuck:` heading, no per-pipe headings)*; `next` with `fresh` in the priority-1 pipe and
   `old` in the priority-2 pipe is `old`, reason "stuck 4 d". An item
   `queued` 3 d minus one second before `now` is not stuck (the threshold is
   strict).
9. **Questions.** Two `waiting` items, one stuck; *Questions for Ola* numbers
   the stuck one `1.`, each line ends with its default.
10. **Recap integration.** `harness_fixtures.run_script(repo,
    "tools/session_state.py", "")` with one pipe and two items, one `red`,
    one `queued`, prints the
    `Pipes (in flight 1 of 3` header after *Waiting on Ola*; `session_format`
    on a file with one `NOW:` and no `QUEUE:` line returns no warning, with
    two `QUEUE:` lines returns one.
11. **Guard.** `guard_governance.judge_bash("python3 tools/work_pipes.py set h31
    red")` is `None`; `governed("tools/work_pipes.py")` is `True` after the
    `GOVERNED` entry (it is `False` today, checked).
12. **The table** *(amended after red; `@tester` adds these)*. In the fixture
    repository with `core` priority 1 and `lsm` priority 2: items `r` `red`
    (branch `worktree-r`, busy by a note), `g` `green` not busy, `w`
    `waiting` with a question, `q` `queued`, all in `core`, `l` `queued` in
    `lsm`, and `m` `merged`: `python3 tools/work_pipes.py` with no arguments
    prints the same text as `show`; under the column heading row (it holds
    `PIPE`, `PRI`, `ITEM`, `STEP`, `WAITING ON`, `SINCE`, `NEXT`) there is
    exactly one row per item but `m` (five rows); only the `g` row starts
    with `next`, matching the `next` command's choice; WAITING ON reads `its
    agent` on `r`, `Ola` on `w`, `free agent` on `q`; and every row's
    columns start at the heading's offsets. Test 8's set covers the `stuck`
    mark and its order; test cap (40 lines of 120) stays as it is, the
    200-character titles now cut at the ITEM width.

The red step's pins (`5fc19f4a`), ruled after it: all seven accepted as
pinned. (1) Item record keys `id`, `title`, `state`, `pipe`, `next`,
`branch`, `pr`, `question`, `default`, `plain`, `fallback`, `since`; the
pipe's own record also holds `priority`, `priority_words` and its open
`question`/`default` (section 5). (2) `Pipe(name, purpose, priority,
source)`, keyword arguments. (3) `order(...)` returns a value with `.item`
(`None` when idle) and `.reason`. (4) `fold` returns `dict[str, Item]` by id,
the type `order` and `stuck` take; `load_items` is `fold` over a file. (5)
`definitions(repo, run=...)`, `run` defaulting to the real git call. (6)
`show` says `1 unreadable`, on its header line (section 7). (7) More than one
`QUEUE:` line warns `N QUEUE lines; at most one expected`.

No mutation round: this suite is not invariant-critical (`docs/increments/
README.md`, *Cost constraints*).

## As built (first slice, green at `c66ae776`)

- **Size.** 340 net counted lines (`python3 tools/count_loc.py master
  c66ae776`: `tools/pipes.py` 333, `tools/session_state.py` 6,
  `guard_governance.py` 1), against about 255 in section 10.
- **Renamed** `tools/pipes.py` to `tools/work_pipes.py` (`b05d7081`): `pipes`
  shadowed Python 3.12's standard-library module.
- **Beyond the design.** Two record keys: `rank` (order inside a pipe; `top`
  sets it below the pipe's lowest) and `before` (the state `ask` left, which
  `answer` restores). A pipe's own record has state `pipe`, listed nowhere as
  an item. `add` refuses a pipe with no definition (`pipes: no pipe <name>
  defined`), so the tool takes no item until `core.toml` exists. An
  unreadable unattended flag counts as unattended (`read_mode` state not
  `off`): the safe side, plain items only.
- **Rule files.** `CLAUDE.md` §3 *Finish before starting*: the sentence on
  where a free writer goes becomes "`python3 tools/work_pipes.py next` names the
  item a free writer takes" (-10 words). `.claude/REQUIRED-READING.md`, the
  `session.md` line: `QUEUE:` at most one; `ASK OLA:` for a refusal or a
  default taken while away; items and their questions go in `tools/work_pipes.py`;
  "only the main session writes it" shortened (+8 words). Net -2 words, not
  the -35 section 10 estimated; the fallback sentence is unchanged.
- **`.claude/pipes/core.toml`, not yet written** (outside both writers'
  limits; who adds it is Ola's call). Exact content, checked with `tomllib`:

  ```toml
  name = "core"
  purpose = "Public rasputin work: the terrain mesh product, its tools, and the harness and process that build them."
  ```

## Review

- **Design review, round 1, 2026-10-25.** `@reviewer`, head `c2b2dc69fc15f9ddeb38f1ecad25a0db2a1053b7`, range `d886847d..c2b2dc69`. Verdict: CHANGES REQUESTED. LOC: 0 (design); table sums to about 220, not 250. Blocking: (1) live lists go in `<git-common-dir>/pipes/`, not `harness/` (REQUIRED-READING 163-164: no agent writes harness state); no guard entry needed (`judge_bash` already passes `python3 tools/pipes.py set`); (2) `ASK OLA:` lines are read by `harness_mode.py`, `guard_unattended.py`, `away.py --back` and `session_state.py`: "Refusals and defaults taken while away stay `ASK OLA:` lines in `session.md`; only questions about an item move to `pipes.py ask`"; `session_format` stops requiring one `QUEUE:` line; (3) private definitions read with `git ls-tree --name-only <remote>/master .claude/pipes/` and `git show`, no fetch (the recap runs in a hook), printing "read from <remote>/master <sha7>, committed <date>"; a missing ref warns, does not fail; (4) `waiting` items do not count toward the cap of three; "writes no governed path" needs an item field set at `add`; (5) a tracked fallback "re-filled when used" needs a commit per use: keep the REQUIRED-READING fallback sentence until `away.py` checks it (second slice), the refusal test moves there too; (6) resurfacing (Ola: "when stuff is stuck behind for a while, they should resurface"): "An item `queued` or `waiting` longer than three days is listed first in the recap under `Stuck:` with its age, and `next` ranks it ahead of fresher items of equal step", about 15 lines plus a test; (7) Status: "Status: **design to build, first slice**, ruled by Ola 2026-10-25 (\"yes, build the first slice\"): rules option A, the rest as recommended. Section 8 holds the rulings." If `session_state.py` imports `pipes.py`, `pipes.py` joins GOVERNED. Not pushed; no CI.
- **Design review, round 2, 2026-10-25.** `@reviewer`, head `bbc7dfdc1ff041527c46f581b78406b09b54a1b0`, range `82d96911..bbc7dfdc`. Verdict: CHANGES REQUESTED. LOC: 0 (design); table about 231. Round 1's seven closed and checked (`judge_bash` passes `pipes.py set`/`ask` with `tools/pipes.py` in GOVERNED; `private/master` 79848102; fallback sentence kept). Blocking, writer order (section 4): (1) busy items: "An in-flight item is busy while `.claude/current-task/` holds a `<persona>-<worktree>-*.md` file for its branch's worktree; steps 1 and 4 skip busy items", tests 3 and 4 mark their in-flight items busy and test 4 gains a `queued` item that `next` picks; (2) "Steps 2 and 3 skip fallback items."; (3) step 4 starts a fourth item unattended: add "under the same count". Wordings: drop "Ties go to the stuck item (section 9), then" from step 1 and "step 1's ties and" from section 9; section 7: "after *Waiting on Ola*, before the `session.md format:` warnings" (`tools/session_state.py` 399-408); test 10: "one `red`, one `queued`". Not pushed; no CI.
- **Design review, round 3, 2026-10-25.** `@reviewer`, head `3b0d2d981dfdc348ddc6b37d6fcfa2e5743ebecd`, range `51b3dc93..3b0d2d98`. Verdict: APPROVED WITH WORDINGS. LOC: 0 (design); about 240. Round 2 closed as given. Applied in this commit: the busy note name is `<persona>-<name>-<HHMMSS>.md` (a longer worktree name no longer collides, e.g. `tri-budget-rms`), with a test 3 case; `next` prints a `busy:` line per skipped item, so a stale note leaves a trace. Not pushed; no CI.
- **Code review, round 1, 2026-10-26.** `@reviewer`, head `8691037448c70ea461c595d9a5c9f7015f4d48a3`. Verdict: CHANGES REQUESTED. Blocking: `answer` dropped the ruling; `answer` without an open question changed state.
- **Code review, round 2, 2026-10-26.** `@reviewer`, head `316adc0a7babc3ad75ad60cdb269c5c78ae7f038`. Verdict: CHANGES REQUESTED. Blocking: CI red on Python 3.12, `tools/pipes.py` shadowed the standard-library `pipes` module; renamed to `tools/work_pipes.py`.
