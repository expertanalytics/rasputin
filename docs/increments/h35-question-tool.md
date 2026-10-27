# Harness h35: the question tool

Status: **designed**, 2026-10-27, on `worktree-question-tool` off master
`7c2aa7c0`. One PR: the tool alone, about 185 counted lines, no governed
file. Wiring it into the recap and the rules in place of the `ASK OLA:` lines
writes governed files and is its own increment,
`docs/increments/h36-questions-in-recap.md`, built after this one with Ola at
the guard's prompt. No C++, no refine or mesh code, no `@perf` run.

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
(`h36-questions-in-recap.md`).

## 1. What it is

`tools/questions.py`, standard library only, not governed.
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
| `list` | prints the block below, stores which number meant which id, and marks the `taken` lines it printed as shown |
| `answer N WORDS` | closes the question numbered N in the last `list`, with Ola's words verbatim; prints the ruling line and any question it releases from `held` |
| `drop N REASON` | closes without an answer (superseded, no longer needed) |

`list` prints exactly what the main session copies to Ola, one running
sequence, questions first, then acts, then the defaults taken since the last
`list`, on one vetoable line each:

```
1. Should the land-cover repair run once per area and be cached? Default: no, repair on every run.
2. Push the question tool and open its pull request? Default: wait for your yes.
3. Taken: the cache lives in the user's cache folder. Veto by number.
```

When any question is held (C3), the block ends with one line, `and N more
yes/no questions held until one above is answered`.

The block is built by a pure function, `list_lines(records) -> list[str]`,
that reads the store and writes nothing; the recap (`h36`) calls only that.
Only the main session's explicit `list` command stores the number mapping and
marks `taken` lines shown, so a recap printed between a `list` and Ola's
answer cannot renumber what he saw or swallow a default taken before he saw
it.

Numbers are assigned at `list` time; `answer 2 "yes"` resolves through the
stored mapping, so a question added after the list was shown cannot shift
what Ola's "2" meant. An `answer` naming a number the last `list` did not
print is refused.

`answer` prints `Ola, 2026-10-27: "<words>"`, which the main session pastes
into the increment file or ROADMAP row the question was about, as rulings are
recorded today. The store keeps every answer; the tracked record stays where
it is now, in the increment files, and only for public work. No command prints
past answers: Ola asked for a tool to pose questions, and C7 already quotes an
earlier answer where it matters.

## 2. The checks (each refuses with one plain sentence; none edits the text)

Constants are chosen, not derived. A probe (script removed after the run)
measured them against the 142 numbered items that carry a `?` and a default
in `docs/increments/` and `docs/retrospectives/` at `7c2aa7c0`, paragraphs
joined: median 220 characters, 90th percentile 628, longest 1297.

*False refusals.* As the sample of good past questions, the probe took the 78
items of 240 characters or less (this file's own excluded): questions Ola
answered in the length the tool allows. The checks below refuse 49. 46 are
true refusals: the item numbers itself (`Q3.`, `P5,`: the tool numbers), or
carries a design label, a commit hash, a branch, or a path inside the
repository. 3 are false: `Bane NOR`, a company's name, read as an
abbreviation; `--term NOR` accepts it. A first version of the checks refused 6
more, now accepted: the compiler flag `-O1` as a label, four persona handles,
`$RASPUTIN_DATA/sources.toml` as a relative path, and `CLAUDE` in `CLAUDE.md`
as an abbreviation.

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
* **C4 self-contained.** Refused, case-folded, only on phrases that point
  elsewhere: `see above`, `as above`, `as before`, `same as above`, `previous
  question`, `question <n>`: "each question must read alone". Bare `above`,
  `below` and `same as` pass (`above 1000 m`, `the same as master`).
* **C5 plain words.** Refused, naming the token, on:
  * a design or harness label: `h<n>`, a capital letter and one or two digits
    such as `P4` or `Q11` (not after `-`, so `-O1` passes), `increment <n>`,
    `row <n>`;
  * a pull request number `#<n>`;
  * a commit hash: 7 to 40 hex characters with at least one digit and one
    letter `a` to `f` (`1400000` passes);
  * a handle `@word` that is not a persona in `.claude/agents/` (`@perf`
    passes; Ola knows the personas);
  * a branch name `worktree-...`, allowed with `--act`, where Ola needs it;
  * a path inside the repository: a token not starting with `/`, `~` or `$`,
    not a URL, that either begins with a tree directory (`tools/`, `docs/`,
    `src_python/`, `tests/`, `include/`, `.claude/`, `.github/`) or has
    directories and ends in a file extension (`a/b.py`); `m/s`, `yes/no`,
    `W/m²` and `on/off` pass;
  * a capitalised abbreviation of 2 to 6 letters, not part of a file name
    (`CLAUDE.md` passes), not in the tool's short list of known ones (`CI`,
    `PR`, `DEM`, `CRS`, `GIS`, `MCP`, `CPU`, `UTC`, `CLI`, `API`, `URL`, `CSV`,
    `JSON`, `EPSG`, `QGIS`, `NVE`, `UTM`, `AC`, `CORINE`, ...), not spelled out
    in the same text as `words (ABC)`, and not named with `--term ABC` (a
    proper name such as `Bane NOR`).
* **C6 short.** Text plus default over 240 characters, two lines of the
  recap's 120: refused, "at most 240 characters; say the rest in the
  increment file". 78 of the probed items fit.
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

## 4. Tests, `tests/python/test_questions.py`

Each runs `command()` on a store in `tmp_path`; inputs are literal strings.

1. C1: `--default ""` and `--default "TBD"` refused; `"no"` accepted.
2. C2: `"Push it? And merge?"` refused; `"Push the tool."` refused; `"Push the tool?"` accepted.
3. C3: two yes/no adds: the second is `held`, `list` prints one line; `answer 1 yes` prints the released one; `--choices "metre|foot"` beside an open yes/no is open.
4. C4: `"Same as above?"` and `"As in question 2?"` refused; `"Keep cells above 1000 m?"` and `"Is it the same as master?"` accepted.
5. C5 refusals, one each: `h33`, `P4`, `#381`, `7c2aa7c0`, `@someone`, `worktree-x` (accepted with `--act`), `tools/questions.py`, `docs/`, `a/b.py`, `ABCD`.
6. C5 acceptances: `-O1`, `1400000`, `@perf`, `/Users/x/tools/questions.py`, `~/x/y.toml`, `$RASPUTIN_DATA/sources.toml`, `https://a.b/c.py`, `m/s`, `yes/no`, `W/m²`, `on/off`, `CLAUDE.md`, `"Large Shapes (LS)" ... LS`, `Bane NOR` with `--term NOR` (refused without it).
7. C6: 240 characters accepted, 241 refused.
8. C7: re-adding an answered question is refused and the refusal quotes Ola's words.
9. `list`: questions, acts, then taken defaults, numbered 1..n in one sequence; a `taken` line appears in one `list` only; with two held questions the last line is `and 2 more yes/no questions held until one above is answered`.
10. `list_lines` is read-only: called between a `list` and an `answer` it leaves the store byte-identical, `answer 1` still means the first `list`'s question 1, and a `taken` line it printed still appears in the next `list`.
11. `answer`: resolves through the stored mapping after a later `add`; a number not in the last list is refused; the words are stored verbatim (quotes, newlines).
12. `drop` closes the numbered question; it no longer appears in `list`.
13. The store: a malformed line is skipped and counted, as in the pipes.

No mutation round; this is not an invariant-critical suite.

## 5. Size

About 185 counted lines: store and records about 35, checks about 65,
`list_lines`, `list`, `answer` and the number mapping about 40, the command
line about 45.
