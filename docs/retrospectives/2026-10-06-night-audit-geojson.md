# 2026-10-06 night: audit PR C (one GeoJSON reader) and the first three hours away

Recorded by `@orchestrator` on 2026-10-06 at about 02:25 Oslo time, on
`worktree-retro-1006` off master `44fa7f5`. Times are Oslo time (UTC+2).

**This covers only part of the window.** Ola entered unattended mode at
23:23:19 for a window that ends at 08:59:19
(`<git-common-dir>/harness/unattended.json`, `windows.jsonl`). This file
covers 23:23 to 02:22, the first 179 minutes. The rest, 6 h 37 min, is for the
morning report or an addendum.

Sources: the main session's transcript,
`~/.claude/projects/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a.jsonl`
(cited as "line N"), and its subagent transcripts in the `subagents/`
directory beside it; the harness queue; the branches `worktree-audit-geojson`
(`32b5092..b26beb8`) and `worktree-geotiff-param-crs` (`80a48a2`); the main
session's lessons list, as passed to me. Handback quotes are the personas'
words, not Ola's. Rule changes below are proposals. Ola decides. The
2026-10-05 retrospective (`d1a2a6c`, unpushed branch `worktree-retro-1005b`)
is not repeated.

Labels: **PR C** is the third Python-audit PR. It merges the GeoJSON readers
into one reading path, `io/geojson.py` (section 10 of
`docs/increments/python-audit.md`). **PR B** is the CRS-helpers audit PR it
stacks on. **The GeoTIFF branch** reads a GeoTIFF whose CRS is given only by
parameters, when the parameters match an EPSG code (Ola's ruling D9: a).
A **red step** commits failing tests and a **green step** commits the code
that passes them. A **pinned list** is the tester's "Pinned or assumed beyond
the design" handback section.

## 1. What happened

- 23:23 away mode on. 23:32, line 5689, Ola: "We're runnig in away mode. Good
  night!"
- GeoTIFF branch, 23:24 to 00:35: design (`fbbfb56`), two design rounds,
  red `0f6c271`, green `c8984e8`, code review round 1, red `5270a64`, green
  `35c44c8`, code review round 2 approved, recorded at `80a48a2` (+159 net
  production lines). Clean: every step in order, every commit in its
  persona's files.
- PR C, 23:14 to 02:20: design and three design rounds (to `32b5092`), red
  `1928569`, green `3b739ac`, then five code review rounds. Approved at
  `efe0eea` (+29), recorded at `b26beb8`. Details in section 3.
- By 00:11 CI had come back after GitHub's outage. #188 and #190 had merged
  through the queue; #189 and #191 had green CI but no auto-merge.
- 00:36: context compaction (line 6509). The `SessionStart` hook printed the
  recap (line 6522), and the main session read `session.md` before its next
  spawn.
- 02:21: this retrospective and a `ROADMAP.md` note on GeoPackage colours
  (`bf2d42f`, `worktree-gpkg-colour`) started. Nothing was pushed.

## 2. Idle time and capacity, 23:23:19 to 02:22:30

Measured from the spawn and handback times in the transcript. For the three
agents whose completion did not arrive as a handback message, the times come
from their result lines (5882, 6031 and 6412).

| Agents running | Minutes |
|---|---|
| none | 7.7 |
| one | 105.8 |
| two | 63.0 |
| three (one writer and two read-only reviewers) | 2.7 |

- **No idle time.** The 7.7 minutes with no agent running are hand-off
  gaps. The longest is 41 s, at 02:20:23.
- **One writer slot empty for 1 h 47 min**, from 00:35:17 (the GeoTIFF
  recording handback, line 6487) to 02:21:52. All of PR C's code review ran
  one agent at a time. Two items that needed no ruling and wrote no governed
  path were queued all along: this retrospective and the `ROADMAP.md` note.
  The `QUEUE:` line rewritten at 00:36 (line 6552) put both "after C review".
  This breaks no rule. Ola's standing note asks for a deep queue, and the
  idle rule counts only a window with no work running. But half the night's
  capacity went unused while work was waiting.
- **No sleep.** `caffeinate` held. No 600 s watchdog stall appears in the
  transcript.

## 3. PR C: five code review rounds, and the pattern

| Round | Ended | Blocking (all prose) | Found besides, by `@architect`'s probe while fixing the prose | Fix |
|---|---|---|---|---|
| 1 | 00:41 | stale `project_structure.md` rows; wording table missing rows; status | this PR made the station readers say `None is a None, not a Point` for a feature with no geometry (line 6676) | red `28fe0cf`, green `693f560` |
| 2 | 01:10 | wording rows missing shapes; known gap 1 false | this PR turned a clean refusal of an empty non-null `geometry` into a traceback in `--features` and `catchment --lakes` (line 6932) | red `af5cb4e`, green `0a11a34` |
| 3 | 01:39 | table lacks `--features` and `--lakes` rows; status | this PR turned the refusal of `"geometry": 7` into a traceback in four readers (line 7126) | probe committed as the record (`9e002e2`); red `5a02bc9`, green `21f49d6` |
| 4 | 02:14 | one wrong summary row | none | prose `efe0eea` |
| 5 | 02:19 | none: approved | | record `b26beb8` |

Every reviewer verdict said "CHANGES REQUESTED, prose only". **All three
defects this PR introduced were found by the design's author probing in
order to write the prose**, not by the red step's tests, not by `@developer`,
and not by `@reviewer`. The brief I was given says two regressions; the
transcript shows three. The first is a wording defect, not a crash, and the
main session called it "introduced" by this PR at line 6676. Code review
took 1 h 46 min (00:33 to 02:19) for a +29-line PR. Three of its five rounds
existed because of the hand-written table.

**Why the reviewer missed them.** Round 1's record says "Behaviour matrix
re-run at base and head ... as designed". That matrix came from the design,
so it tested the cases the design had foreseen. The `@architect` lesson
(rounds 2 and 3) names the cause: "a probe built from the bug as found
misses shapes only the fix makes reachable, and callers that share the bad
input but not the fix; build it from the code as fixed with every empty value
(None, "", 0, [], {}) per field tested." `REQUIRED-READING.md` already says
"Derive the probe set from the code as fixed"
(`.claude/REQUIRED-READING.md@44fa7f5:68`). The principle was in force. What
was missing was a way to apply it. Each round's probe grew: 22 shapes × 7
readers in round 2, 56 × 8 in round 3, 68 × 8 committed. The record then
stopped being hand-written. MAST names this kind of failure "No or incomplete
verification" (FM-3.2), and notes that verifiers "perform only superficial
checks, despite being prompted to perform thorough verification" (Cemri et
al. 2025, arXiv 2503.13657; this quotation and the Wikipedia one below are
as returned by a fetch of the page on 2026-10-06, not checked against the
PDF).

**What fixed it.** At 01:45 (`9e002e2`) the probe,
`docs/increments/python-audit-probes/geojson_wordings.py` (202 lines),
became the record, with its output committed at the base and at the head.
The diff of the two outputs is the full list of changed wordings, and the
doc's table became a 9-row summary. Round 4 re-ran the probe "byte for byte"
and found one wrong row. Round 5 approved. UNCAUGHT lines (tracebacks) fell
from 50 at the base to 2, both a known gap that was already there at the
base.

**Is it worth a rule? Yes, but a narrow one, placed at design time.** The
method is old. Feathers calls it a characterization test, or golden master:
it "describe[s] the actual behavior of an existing piece of software, and
therefore protect[s] existing behavior ... against unintended changes"
(en.wikipedia.org/wiki/Characterization_test, read 2026-10-06). Running two
versions on the same inputs and diffing the outputs is differential testing
(McKeeman, "Differential Testing for Software", Digital Technical Journal 10,
1998; abstract read, paper not read). Anthropic's Claude Code guide lists "a
script that diffs output against a fixture" as a check an agent can run
(code.claude.com/docs/en/best-practices, read 2026-10-06). Tonight the probe
came after the green step, so each defect it found cost another
red/green/review cycle. Had it been run at the base and committed with the
design, the red step could have taken its tests from the base output, and
the green step's diff would have shown all three defects before review. See
proposal P1.

## 4. Deviations

**D1. Six tester pinned lists went to green with no `@architect` ruling.**
`.claude/agents/tester.md@44fa7f5:27-29`: "`@architect` confirms or rules on
each before green." Each red handback went straight to a `@developer` spawn
16 to 26 s later: lines 5852→5873, 6275→6296, 6329→6351, 6682→6703,
6938→6958 and 7138→7158. PR C's first red step pinned 15 items. The
developer's brief did not contain them. The developer reported: "**The
tester's pinned list (15 items) was not in my brief.** I worked from the
tests". Another tester (line 6682's handback) asked: "`@architect` should
confirm." That pin was the `feature 0` index in the refusal wording, which
went into the code unruled. Pin 7 of the first list ("A bare geometry is
wrapped with `geometry == the document itself`") is the behaviour behind
round 1's `None is a None` wording. An architect reading it against the
readers might have caught that before review. That is a guess, not checked.
MAST calls this "Information withholding" (FM-2.4).

**D2. A second enqueue tried on a grant for one.** At 00:11:16 (line 6202)
the main session ran `gh pr merge 189` and `gh pr merge 191`, because their
auto-merge had dropped. The rule (`.claude/REQUIRED-READING.md`, *A grant
covers one occurrence*): "one yes covers one enqueue". In away mode
`guard_push.py` refused it and queued it (the only queue entry in the
window). The session did not retry, and it wrote the `ASK OLA:` line. The
guard did its job. The `ASK OLA:` line still says "(you approved pushing
1-4)", which reads as if the old grant covered the re-enqueue. It does not;
Ola needs to give a fresh yes.

**D3. Two unchecked facts in the main session's briefs.** At about 02:14 it
told `@architect` that `--features` reads a bare geometry under the CORINE
maps. Two of the three maps refuse it (`efe0eea`, lines 7362 and 7417). At
02:21 it said the Graz GeoPackage uses the official CORINE legend colours.
The file uses rasputin's own palette (line 7584). The main session caught
both itself, and `REQUIRED-READING.md` *Claims* already covers them. No new
rule. As a practice, pass on another agent's claim as "the reviewer says X
(unchecked)".

**D4. A nominal increment file for a prose branch.** `tools/brief.py` refused
a reviewer brief for the `ROADMAP.md` branch ("@reviewer needs --increment",
line 7597). The session then passed `--increment ROADMAP.md` and finally
`--increment docs/increments/README.md` (line 7640). That reviewer's brief
now names README as "the increment file". No guard was worked around, but
the tool forced a false statement into the brief. See P4.

**D5. Note files.** Of 25 writing agents in the window, 20 wrote their note
file in their first seven tool calls. Two wrote it near the end: the PR C
round-3 prose `@architect` (call 48 of 50) and the GeoTIFF design `@architect`
(call 67 of 82). Three wrote none: `@tester` `a08253b7`, `@tester` `a81f5455`
and `@architect` `a29f4270`. No `@reviewer` can write one: `reviewer.md` lists
no Write or Edit, and shell writes are out of its role. Yet `tools/brief.py`
gives every persona "Write limit: ..., and your note file."
(`tools/brief.py@44fa7f5:265`).

**D6. Increment file against persona file on "throwaway".** PR C's red-step
section says "Lean: no throwaway implementation, no mutation round."
(`docs/increments/python-audit.md@b26beb8:1730`). `tester.md` requires the
named fix to be tried in a scratch copy
(`.claude/agents/tester.md@44fa7f5:76`). The tester followed the persona file
and reported that the run "found the one older test the fix breaks". In
`docs/increments/README.md@44fa7f5:114`, "throwaway implementation" means the
mutation-testing stand-in, so the two texts mean different things by the
same word. Fix the phrase, not the rule (P5).

**Not deviations.** Every commit in `32b5092..b26beb8` touches only its
persona's files (`@tester`: `tests/`; `@developer`: `src_python/`;
`@architect`: `docs/` and `project_structure.md`). Every fix went red, green,
then review. `guard_spawn.py` refused an edited brief at 02:24:00 (line 7620),
and the session respawned with a fresh block 22 s later. Neither branch
touches refine or mesh code, so no `@perf` run was due.

## 5. Guards: refusals and false positives

- **Window queue:** one entry, the D2 merge. A correct refusal.
- **My own run, two refusals.**
  - A read-only Python one-liner that opened the queue files was refused as
    a harness write ("Only Ola enters or leaves unattended mode, and only
    hooks and away.py write the harness state"). The guard judged a read
    from its text. False positive.
  - A planted-import check in a scratch copy under the session scratchpad
    was refused as a write to `tools/rule_sizes.py`. The shell had `cd`'d
    into the scratch copy, so the relative path pointed into the copy, not
    the repository. False positive in substance. Because the refusal names
    a governed path, I did not redo it by any route. So I have not checked,
    by planting, the `@architect`'s finding that `check_prohibited_deps.py`
    misses imports under `docs/`. I read the mechanism instead:
    `SOURCE_DIRS` omits `docs`
    (`tools/check_prohibited_deps.py@44fa7f5:27`).

## 6. Lessons the personas reported, grouped

**Probes (5).** Build the probe from the code as fixed, with every empty
value (`None`, `""`, `0`, `false`, `[]`, `{}`) for each field, and across
every caller that shares the input (`@architect`; `@reviewer` rounds 2 and
3). A probe's UNCAUGHT can be the probe's own bug (`@architect`). Stray
stderr broke a parse of `tools/scratch_copy.py` output; read stdout only
(`@reviewer`). A scratch copy of `src_python` is shadowed by the editable
install's import hook, and `PYTHONPATH` does not beat it: drop the
`editable` finder from `sys.meta_path` and assert `tin_engine.__file__`
(`@reviewer` round 1, `@architect`). Probe `.py` files under `docs/` escape
the prohibited-dependency gate (`@architect`; P2).

**Hand-offs (3).** The pinned list must reach the developer through the
brief or a tracked file (`@developer`; D1). A status paragraph that says
what the next round will check goes stale when the fix does more; say what
the fix changed (`@architect`). A claim in a brief needs the same check as
any other claim (main session; D3).

**zsh (2), recurring.** `"$T:path"` applies a history modifier; write
`"${T}:path"` (`@reviewer` round 1). Yesterday's retrospective recorded the
same trap with `"$rev:..."`. `echo ====` fails, because zsh treats a leading
`=` as a command lookup (`@architect`). This is the second night running,
so it now earns one line in `.claude/briefs/common.md` (about 20 words):
"zsh: write `${v}:path`, never `$v:path`; do not start an argument with
`=`."

**The stderr noise has a root cause on this Mac (checked).** Every bare
`python3` prints two `RuntimeWarning: Unexpected value in sys.prefix` lines.
`/opt/homebrew/bin/pyvenv.cfg`, dated 2023-09-12, describes a Python 3.10
virtualenv. Python 3.14 finds it beside the executable and takes
`/opt/homebrew/bin` as its prefix. `/opt/homebrew/bin/python3` prints prefix
`/opt/homebrew/bin` with the warnings. The same binary called by its Cellar
path prints the right prefix and no warnings. The main session appended
`2>/dev/null` to 69 of its 80 shell commands in this window, and one of them hid a failed
`session.md` edit (line 6208: the replace matched nothing, and `grep -c`
printed 0). The `@reviewer`'s polluted-parse lesson is the same noise. The
fix is outside the repository, so it is Ola's (question 2).

## 7. Proposals (Ola decides)

**P1. A differential probe is the record of changed behaviour (rule, about
70 words).** In `docs/increments/README.md`, design step: "A PR that changes
what existing entry points accept, or how they refuse, commits with its
design a probe that feeds every entry point every input shape (each field
the code reads, with each empty value and a wrong type), and the probe's
output at the base. The green step commits the output at the head. The diff
is the record; the design's table only summarises it, and `@reviewer`
re-runs the probe." Incidents: `591c026`, `9e002e2` (missing rows, rounds 1
to 3), and the fixes `693f560`, `0a11a34`, `21f49d6` (defects found after
green). Owner: whoever is briefed for `docs/increments/README.md`. It is a
rule file, so it is done by day. Cost: about 70 words, and the probe is
written at design time instead of late. Tonight's was 202 lines. Expected
saving on a PR like C: three red/green/review cycles, about 1 h 30 min.

**P2. Probe scripts go through the gates.** Either probes live under
`tools/probes/`, or `SOURCE_DIRS` gains `docs/increments`
(`tools/check_prohibited_deps.py@44fa7f5:27`); ruff and mypy exclude `docs/`
too. Incident: the `@architect`'s planted `import rasterio`, not caught (its
handback; my own plant was refused, section 5). Owner: `@developer` through
the pipeline. It writes `tools/`, so by day. Cost: one line plus a test. Do
it together with P1, or P1 adds ungated code.

**P3. The pinned list goes into the increment file in the red commit.**
Change `.claude/agents/tester.md@44fa7f5:27-29` from "under the handback heading" to "in the
increment file's red-step section, in the red commit". `@architect` then
rules in the file, and `@developer` reads it there. Incident: D1, six
hand-offs. Owner: `@tester`'s file (governed, so by day). Cost: about 10
words. It also suits stage 2 of `docs/research/2026-10-03-dispatcher-control.md`
(a script that says what step comes next on each branch): "pins unruled"
becomes a state that script can see.

**P4. `tools/brief.py --increment none` for prose-only branches.** The brief
would then say "no increment file (prose branch)" instead of naming
README. Incident: D4. Owner: `@developer`. Cost: about 5 lines. Low priority.

**P5. Say "no mutation round", not "no throwaway".** In the main session's
memory note on lean briefs and in increment red-step sections, use "no
mutation round". The scratch-copy try of a named fix (`.claude/agents/tester.md@44fa7f5:76`) stays.
Incident: D6. Cost: none; no rule file changes.

**Also from tonight, not a proposal:** the empty writer slot (section 2).
The main session can start a queued item that needs no ruling as soon as a
slot is free, instead of queueing it "after" a running item. This is already
the intent of Ola's note on filling the window.

## 8. Rule text: size and one cut

`python3 tools/rule_sizes.py` at `44fa7f5`: **10,564 words**, the same as
the 2026-10-05 measure. No rule file has changed since then (on
`origin/master` either, `git diff 44fa7f5 origin/master` over the rule
files is empty). Yesterday's cut, a one-line review record in
`docs/increments/README.md` step 4, has not landed. Tonight's PR C code
records are each one line of Markdown, but run 78 to 211 words (rounds 1 to
5: 125, 165, 211, 108, 78), against the 45-word format proposed.

**Tonight's cut: the first paragraph of `.claude/REQUIRED-READING.md`'s "The
harness"** (lines 137 to 155 at `44fa7f5`, 196 words). It lists each hook's
commands, and each hook's docstring already holds that list. Agents learn of
a guard when it refuses, with the reason, and tonight's one overnight
refusal came on an act the list names (`gh pr merge`), so the list did not
prevent it. Replace it with:

> **The harness.** Hooks in `.claude/settings.json` ask before a publishing
> act, a write to a git ref, remote or config, or a write to a rule file;
> while Ola is away they refuse these, and questions to him and
> configuration changes, instead. They refuse a spawn without an unchanged
> `tools/brief.py` block (a resume that starts a new step carries a fresh
> one), and run the gates after a commit. Each hook's docstring lists what
> it guards. They are tripwires: the boundary is still yours, and auto mode
> can let an unapproved push through.

That is about 90 words, so about 105 fewer. Risk: a reader no longer sees in
one place that `git config` writes and `gh api` writes are guarded. The
refusal message names them, and so does `guard_push.py`'s docstring.

## 9. Questions for Ola

1. **P1, the committed differential probe as the record, written into the
   design step?** Default: yes, together with P2.
2. **The stray `/opt/homebrew/bin/pyvenv.cfg` (2023, Python 3.10) makes
   every `python3` print two warnings, and agents hide them with
   `2>/dev/null`, which also hides real errors. Rename it to
   `pyvenv.cfg.bak`?** This is your machine, outside the repository; no agent
   touches it. Default: yes, you rename it.
3. **The "The harness" cut in section 8?** Default: yes, in the same PR as
   the cuts you have already approved.
4. **P3, the pinned list goes into the increment file?** Default: yes.
5. **Re-enqueueing #189 and #191** needs your fresh yes; the overnight
   attempt was refused (D2). Default: yes, once you have looked at them.
