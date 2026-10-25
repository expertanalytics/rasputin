# Harness h31: development pipes

Status: **design note for discussion**, 2026-10-25, on `worktree-pipes-design`.
Nothing is built; nothing here is a build plan until Ola has ruled on
section 8. Written by `@architect`.

Ola, 2026-10-25: "When I'm back, we should make a few development pipes. Then
we add tasks to these pipes. We are limited by the number of agents, but this
way we don't get things left behind." Settled the same day: a plugin system,
with a core pipe in this repository and further pipes added by the private
repository; priorities set by Ola as he goes, with a fallback so no agent sits
idle. Open: rules inside pipes (section 5).

Words used here. **Pipe**: a named, ordered list of items. **Item**: one
PR-sized piece of work with a state. **Writer**: an agent that changes files
(at most two run at once, `CLAUDE.md` §3). **Runnable**: an item whose next
step needs no answer from Ola and, while away, writes no governed path.

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
them, a fallback, plugin pipes.

*Literature.* Tooling; no novelty claimed. Three borrowed ideas: Kanban's
work-in-progress limit and swimlanes (Anderson, *Kanban*, 2010: a lane per
class of work, a cap on items in progress; the lean slogan is to stop starting
and start finishing),
which is what `CLAUDE.md` §3 *Finish before starting* already states for one
lane; a work-conserving scheduler, which never leaves a resource idle while
any queue has ready work (Kleinrock's term, *Queueing Systems* vol. 2, 1976),
which is Ola's fallback; and strict priority with a floor against starvation,
as in Linux's scheduling classes, so the lowest pipe still moves when the
higher ones are blocked. The plugin side copies this repository's own
`docs/plugins.md`: a plugin adds commands by declaring a file the core
discovers; the core never names the plugin.

Ideas-file entries this note takes up or overlaps (`../rasputin_scratch/ideas/ideas.md`):
*Recap: items in flight and open ideas* (the recap shows only the current
branch); *Lead time from design to merge, per increment* (item states with
times give it without a second tool); *Morning check proposals*, P8
(`session.md`'s format enforced).

## 3. What a pipe is

**Recommendation.** A pipe is two things, kept apart because they change at
different speeds:

- Its **definition**, tracked: `.claude/pipes/<name>.toml`, a few lines: the
  name, one sentence of purpose, its standing fallback (section 4), and later
  its checks (section 5). It changes through a PR, like a rule.
- Its **live list**, untracked: `<git-common-dir>/harness/pipes/<name>.jsonl`,
  beside the away flag and the refusal queue, written only by `tools/pipes.py`.
  One JSON object per item: id, title, state, next step, branch, PR number,
  the open question if any, size estimate, and the time of each state change.

States: `queued`, `designing`, `red`, `green`, `review`, `waiting` (on Ola),
`merged`, plus `parked` and `dropped` (Ola's words, kept so the item is not
re-proposed). Each change is one command, `python3 tools/pipes.py set <item>
<state>`, by the main session; the recap prints every pipe as one block, each
item one line. Where git can tell the state (a red commit, a `## Review`
record, a merge on master), the tool prints a warning when the file disagrees,
and nothing more in the first slice.

**How the private repository registers pipes without this repository naming
them.** The private remote is a second remote of the same checkout, and its
`master` holds commits this repository's `master` does not
(`tools/push_grant.py`, `home`). A private pipe definition is a file
`.claude/pipes/<name>.toml` in one of those commits. `pipes.py` reads
definitions from the working tree and from `git show <remote>/master:.claude/pipes/`
for every remote other than `origin`; the public code says "every other
remote", never a name. Live lists live in the harness directory, outside both
repositories, so a private item's title is on nobody's `master`. The guard
that refuses writes to harness state (`guard_governance.py`) must pass
`pipes.py`'s writes, as it passes `away.py`'s; that is one entry.

*Alternatives.* (a) One tracked file for all pipes: every state change is a
commit, conflicts between worktrees, and private names would have to be in
it. (b) Live lists tracked too: the same conflicts; the private list would
need a private branch per change. (c) No files, states derived from git and
increment files (stage 2's `pipeline.py` alone): clean for items in flight,
but a queued item has no git trace, which is exactly the item that gets left
behind. (d) Keep `session.md`, enforce its format: the format was ruled
2026-10-03 and today's file breaks it on 37 of 71 lines; a format nothing reads is
not enforced by a warning.

## 4. Priorities, the writer rule, and the fallback

**Setting priorities.** Ola says it in chat, in plain words ("core first, then
the private work"; "put the showcase on top"). The main session turns it into
one command, `python3 tools/pipes.py priority core=1 <name>=2 ...`, which
stores an integer per pipe (1 highest) with the time and Ola's words; `pipes.py
top <item>` moves an item to the head of its pipe. `pipes.py next` prints
which item would get the next free writer and why, so Ola sees the effect
before an agent is spawned. No file edits by hand, no numbers Ola must
remember: the recap shows the current order.

**What a free writer takes.** In this order, the first that applies:

1. Across all pipes, the runnable item nearest to merging (`review` before
   `green` before `red` before `designing`); this keeps `CLAUDE.md` §3's rule
   and Kanban's "finish first". Ties go to the higher-priority pipe.
2. The highest-priority pipe's first runnable `queued` item, if the in-flight
   count is under three.
3. The next pipe by priority with a runnable item.
4. Any pipe's standing fallback (below), highest priority first.
5. Idle, printed as `IDLE: no runnable item in any pipe`, which the morning
   check counts.

The cap stays three items in flight across all pipes (`CLAUDE.md` §3). To
stop one pipe from holding all three while another has ready work, a pipe
holds at most two in flight unless it is the only pipe with runnable items.
*Alternative:* weighted shares (pipe A gets two writers of three); more rule
text, and priorities already express it.

**The fallback.** Each pipe's definition names a standing fallback: an item
that needs no ruling and writes no governed path, re-filled when used (today's
`FALLBACK:` line, the basin-scale pieces, is the core pipe's). `away.py`'s
enter step refuses, with the reason, when no pipe has one, replacing the
sentence in `.claude/REQUIRED-READING.md` that asks the `QUEUE:` line to name
one. Step 4 above is what keeps a writer busy when every top item waits on
Ola, which is the case Ola described.

## 5. Questions to Ola, attached to items

**Recommendation.** `python3 tools/pipes.py ask <item> "<question>" --default
"<what is taken if silent>"` puts the question on the item and sets it
`waiting`; only that item stops, and the rule above picks the next one.
`pipes.py answer <item> "<Ola's words>"` records the ruling on the item and
returns it to its previous state. The recap lists open questions grouped by
pipe, numbered in one sequence, each with its default, which is the form Ola
asked for. A question about the harness or process, not about an item, goes
on the core pipe's standing process item. The 44 `ASK OLA:` lines of today
become at most one per waiting item, and an answered one disappears instead
of staying in the file.

*Alternative:* keep `ASK OLA:` lines in `session.md`, prefixed with the item
id. Cheaper (no tool), but it is what exists today, and the file shows the
result.

## 6. Rules in pipes

Ola: "I think rules could be a useful concept in pipes. But we need to discuss
it more." Three options; the cost is measured in rule words, which
`CLAUDE.md` §2 says must be paid for.

**A. Rules stay global; pipes carry tasks only.** A pipe definition has no
rule field. Cost: 0 words. The private repository keeps its own rule copies
as today, and the drift between them (12,080 words against 9,494 public, per
an `ASK OLA:` line in today's `session.md`) stays a separate question.

**B. A pipe adds checks from a fixed menu (recommended, after a first week
on A).** The definition may name: extra skills every persona loads; an
acceptance command run before review (as `tools/bench.py` is for refine and
mesh code); a reviewer focus line; a size budget below 700 lines. `brief.py
--pipe <name>` appends these to the brief as lines the hash covers. The menu
is closed, so a pipe can tighten but never loosen a rule. Cost: about 50
words in `docs/increments/README.md` (the menu and that a pipe only
tightens) and about 30 lines in `brief.py`; pipe definitions are counted by
`tools/rule_sizes.py` so growth there shows in every recap. This is the
plugin hook Ola asked for, in the same shape as `docs/plugins.md`: a file the
core discovers, a closed interface.

**C. Free-text rules per pipe,** a Markdown file the brief includes. The most
flexible, and the thing §2 fights: unbounded, unreviewed by this repository's
`@reviewer`, and a second place for every persona to read before acting.
Rejected.

## 7. What changes, and the first slice

| Where | Change | Size |
|---|---|---|
| `tools/pipes.py` (new) | read definitions (tree and other remotes), live lists, `set`, `add`, `ask`, `answer`, `priority`, `top`, `next`, `show` | about 220 lines |
| `tools/session_state.py` | print `pipes.py show` instead of `session.md`'s `QUEUE:` lines; keep `NOW:` | about +15, -20 |
| `tools/away.py` | refuse to enter when no pipe has a fallback | about 10 |
| `.claude/hooks/guard_governance.py` | pass `pipes.py`'s writes under `harness/pipes/` | about 3 |
| `.claude/REQUIRED-READING.md` | `session.md` keeps `NOW:` only; questions and queue move to `pipes.py`; fallback sentence replaced | net about -40 words |
| `CLAUDE.md` §3 | *Finish before starting*: "a free writer goes first ..." becomes "`python3 tools/pipes.py next` says which item" | net about -20 words |
| `.claude/pipes/core.toml` | the core pipe, with the basin pieces as fallback | 6 lines, not counted |
| `tests/python/test_pipes.py` | the order rule (steps 1-5), the per-pipe cap, private definitions from a second remote in a temporary repository, a `waiting` item skipped, the refusal with no fallback | uncounted |

Governed paths throughout, so built by day, with Ola at the prompt. **First
slice:** `pipes.py` with definitions, live lists, `set`/`add`/`ask`/`answer`/
`priority`/`next`/`show`, the recap change and the two rule-file cuts: about
250 net lines, one PR under the 700 cap, option A for rules. The `away.py`
check, `brief.py --pipe` and option B are the second slice. Today's
`session.md` is moved into the live lists by hand once, item by item, and its
stale `ASK OLA:` lines are dropped with Ola, which is the first use of the tool.

## 8. Open questions for Ola

1. Two files per pipe (a tracked definition, an untracked live list) as
   section 3 recommends? Default: yes.
2. Priorities as one integer per pipe set by `pipes.py priority`, from your
   words in chat? Default: yes; no hand edits.
3. The writer order of section 4 (nearest to merging first, then priority,
   then fallback), with at most two in flight per pipe? Default: yes, both.
4. Questions attached to items, one `waiting` state, answered questions
   deleted from the list? Default: yes.
5. Rules: option A now, option B after a week of use? Default: yes; C is not
   offered.
6. Which pipes to start with? Default: two: `core` (this repository: product,
   harness, process) and one private pipe, defined on the private `master`.
7. Move today's 72-line `session.md` into the pipes by hand on the first day,
   dropping the stale asks with you? Default: yes, about 20 minutes with you
   at the prompt.
8. `@orchestrator` owns the process side (your 2026-10-25 note): its
   post-merge check reads `pipes.py show` for lead time per item, replacing
   `docs/retrospectives/lead_time.py` later? Default: later, not in this
   increment.
