# Harness h35: the question tool

Status: **designed**, 2026-10-27, on `worktree-question-tool` off master
`7c2aa7c0`. Two PRs on this one design: the first is the tool alone, about
180 counted lines, no governed file, so it is built and pushed while Ola is
away; the second wires it into the recap and the rules, writes governed files
(`tools/session_state.py`, `tools/work_pipes.py`, `tools/harness_mode.py`,
`.claude/REQUIRED-READING.md`, `CLAUDE.md`), and is built with Ola at the
guard's prompt. No C++, no refine or mesh code, no `@perf` run.

Source: Ola, 2026-10-27: "Like, a question tool, so I don't have to repeat how
you should pose questions all the time :)". The rules it checks are the ones he
has given across many corrections (the main session's memory notes on
complete, distinct-answer, few and short questions, plain words, no labels).

## Prior art: legacy and literature

*Literature.* Survey methodology already names two of the faults checked here:
the double-barrelled question (two decisions in one) and the leading or
unanswerable question; Dillman, Smyth and Christian, *Internet, Phone, Mail,
and Mixed-Mode Surveys* (4th ed., Wiley, 2014), ch. 4, advise one concept per
question and a stated answer set. The tool follows that and adds what a survey
lacks: a default for silence and one running number across askers. No
novelty is claimed.

*Legacy.* `git grep -n -i -E 'ask_ola|questions for|def ask' legacy-archive -- legacy`
returns nothing; nothing is carried across.

*Ideas file.* No entry in `ideas.md` is taken up or overlaps.

*In the tree.* `tools/work_pipes.py` `ask`/`answer` hold one question per
item, with a default, and `show` numbers them; `.claude/current-task/session.md`
holds `ASK OLA:` lines that the recap lists unnumbered. The tool replaces both
(second PR).

## 1. What it is

`tools/questions.py`, standard library only, not governed in the first PR.
Its store is `<git-common-dir>/questions.jsonl`, untracked, appended only, one
JSON record per line, the newest record per id wins (the pipes' pattern), so
every worktree sees one list and nothing about private work reaches git. Only
the main session runs it: an agent's questions arrive in its handback, as
today (`.claude/REQUIRED-READING.md`, *Data, scratch and temp folders are not a
channel*).

A record: `id` (`q1`, `q2`, ... never reused), `kind` (`question`, `act` for a
publishing act awaiting a yes, `taken` for a default taken), `text`,
`default`, `choices` (empty means yes/no), `item` (optional pipes item id,
never printed), `source` (who asked: persona or `main`), `state` (`open`,
`held`, `answered`, `dropped`), `words` (Ola's answer, verbatim), `at`.

Commands:

| Command | Effect |
|---|---|
| `add TEXT --default D [--choices A\|B\|C] [--act] [--item ID] [--source WHO]` | runs the checks of §2; on a refusal prints the reason, exit 1, stores nothing |
| `taken TEXT --default D [--source WHO]` | records a small item decided on its default (checks C1, C4, C5, C6) |
| `list` | prints the block below and stores which number meant which id |
| `answer N WORDS` | closes the question numbered N in the last `list`, with Ola's words verbatim; prints the ruling line and any question it releases from `held` |
| `drop N REASON` | closes without an answer (superseded, no longer needed) |
| `rulings [--since DATE]` | the answered records, newest first |

`list` prints exactly what the main session copies to Ola, one running
sequence, questions first, then acts, then the defaults taken since the last
`list`, on one vetoable line each:

```
1. Should the land-cover repair run once per area and be cached? Default: no, repair on every run.
2. Push the question tool and open its pull request? Default: wait for your yes.
3. Taken: the cache lives in the user's cache folder. Veto by number.
```

Numbers are assigned at `list` time; `answer 2 "yes"` resolves through the
stored mapping, so a question added after the list was shown cannot shift
what Ola's "2" meant. An `answer` naming a number the last `list` did not
print is refused.

`answer` prints `Ola, 2026-10-27: "<words>"`, which the main session pastes
into the increment file or ROADMAP row the question was about, as rulings are
recorded today. The store is the full ruling record; the tracked record stays
where it is now, in the increment files, and only for public work.

## 2. The checks (each refuses with one plain sentence; none edits the text)

Constants are chosen, not derived. The probe below measured them against the
139 numbered items that carry a `?` and a default in `docs/increments/` and
`docs/retrospectives/` at `7c2aa7c0` (paragraphs joined; script removed after
the run): median 220 characters, 90th percentile 628, longest 1297; 77 carry
an internal label, 10 a relative path, 12 a reference to another item.

* **C1 default.** `--default` is empty or one of `none`, `tbd`, `?`, `-`
  (case-folded): refused, "a question needs a default: what happens if Ola is
  silent".
* **C2 one decision.** The text has more than one `?`, or does not end with
  `?` (`question` and `act` kinds): refused.
* **C3 distinct answers.** A new yes/no question while another yes/no question
  or act is open, or a new `--choices` question whose choices share a word
  with an open question's choices: stored as `held`, not refused, with the
  printed reason "held until question qN is answered: the same short answer
  would fit both". `list` does not show held questions; `answer` releases the
  oldest held one. (The brief said refuse; holding keeps a question from an
  agent's handback from having no home, and it is printed, so not silent.)
* **C4 self-contained.** Refused on `see above`, `above`, `below`, `as
  before`, `previous question`, `question <n>`, `same as`: "each question must
  read alone".
* **C5 plain words.** Refused, naming the token, on: a design or harness label
  (`h<n>`, a capital letter and one or two digits such as `P4` or `Q11`,
  `increment <n>`, `row <n>`); a pull request number `#<n>`; a commit hash (7
  to 40 hex characters with a digit); a persona handle `@word`; a branch name
  `worktree-...` (allowed with `--act`, where Ola needs it); a path with a `/`
  that does not start with `/` or `~` and is not a URL; a capitalised
  abbreviation of 2 to 6 letters not in the tool's short list of known ones
  (`CI`, `PR`, `DEM`, `CRS`, `GIS`, `MCP`, `CPU`, `UTC`, `CLI`, `API`, `URL`,
  `CSV`, `JSON`, `EPSG`, `QGIS`, `NVE`, ...) and not spelled out in the same
  text as `words (ABC)`.
* **C6 short.** Text plus default over 240 characters, two lines of the
  recap's 120: refused, "at most 240 characters; say the rest in the
  increment file". 78 of the 139 probed items fit.
* **C7 not re-asked.** The text, case-folded with whitespace collapsed, equals
  an open, held or answered record's: refused, quoting the earlier answer and
  its date when there is one.

What the tool cannot judge (whether a question changes direction, whether a
plain sentence is clear) stays the main session's; `taken` is the route for
the rest.

## 3. Handbacks and away mode

A persona's handback keeps its *Questions for Ola* heading, each with a
default. The main session runs `questions add ... --source <persona>` for
each; on a refusal it rewrites the question itself and adds again, or, when
the default is clearly fine, records it with `taken`. While Ola is away a
refused guarded act becomes `add --act`, and a default taken becomes `taken`;
the first report after the window is `list`.

## 4. Second PR: replacing ASK OLA lines

* `tools/session_state.py`: the recap's *Waiting on Ola* block prints
  `questions.list_lines()` instead of reading `ASK OLA:` lines; `session.md`
  holds `NOW:` and `QUEUE:` only, and an `ASK OLA:` line is warned as an
  unknown line. `tools/questions.py` joins the governed list, since a hook then
  imports it.
* `tools/work_pipes.py`: `ask` and `answer` are cut (no alias); `show` marks an
  item `waiting` when it has an open question in the store (`item` field), and
  prints no question block of its own.
* `tools/harness_mode.py` and `.claude/hooks/guard_unattended.py`: the refusal
  text says `questions add --act` (subagent: in its handback) instead of the
  `ASK OLA:` line.
* Rule text, net words cut: `CLAUDE.md`'s "Two or more questions or acts
  awaiting Ola are numbered in one sequence, each explained, with its default"
  becomes "Questions to Ola are the output of `python3 tools/questions.py
  list`"; `.claude/REQUIRED-READING.md`'s `ASK OLA:` paragraphs shrink to a
  pointer. Counted with `tools/rule_sizes.py` in that PR.

About 40 counted lines added and 45 removed.

## 5. Tests (first PR), `tests/python/test_questions.py`

Each runs `command()` on a store in `tmp_path`; inputs are literal strings.

1. C1: `--default ""` and `--default "TBD"` refused; `"no"` accepted.
2. C2: `"Push it? And merge?"` refused; `"Push the tool."` refused; `"Push the tool?"` accepted.
3. C3: two yes/no adds: the second is `held`, `list` prints one line; `answer 1 yes` prints the released one; `--choices "metre|foot"` beside an open yes/no is open.
4. C4: `"Same as above?"` refused naming `above`.
5. C5: one refusal each for `h33`, `P4`, `#381`, `7c2aa7c0`, `@tester`, `worktree-x` (accepted with `--act`), `tools/questions.py` (accepted as `/Users/x/tools/questions.py` and `https://a.b/c`), `ABCD`; `"Large Shapes (LS)" ... LS` accepted.
6. C6: 240 characters accepted, 241 refused.
7. C7: re-adding an answered question is refused and the refusal quotes Ola's words.
8. `list`: questions, acts, then taken defaults, numbered 1..n in one sequence; a `taken` line appears in one `list` only.
9. `answer`: resolves through the stored mapping after a later `add`; a number not in the last list is refused; the words are stored verbatim (quotes, newlines).
10. `drop` closes; `rulings` lists answered records newest first.
11. The store: a malformed line is skipped and counted, as in the pipes.

No mutation round; this is not an invariant-critical suite.

## 6. Size

First PR about 180 counted lines: store and records about 35, checks about 60,
`list`/`answer` and the number mapping about 40, the command line about 45.
Second PR: about 40 added and 45 removed, governed.
