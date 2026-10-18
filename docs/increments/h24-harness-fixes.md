# Harness h24: fixes from the night of 2026-10-17/18 and the morning check

Status: **designed**, 2026-10-18, on `worktree-h24-harness-fixes` off master
`1806288c`. Five PRs (section 4). Ola, 2026-10-18 07:04 UTC: "yes, start the
harness design". Items 2 and 3: Ola, 2026-10-18, "yes to both"; everything
else has his yes or needs none. Revised for design review round 1 (items 1,
5, 6, 13 and sections 4 and 6), and for PR A's code review round 1 on
`188089f1` (items 2 and 5, new item 14, sections 4 and 5), and for design
review round 3 on `a4c19f4e` (items 5 and 14).

Sources: `docs/retrospectives/2026-10-18-merges-284-285.md` (sections 3, 5
and 7), `docs/retrospectives/2026-10-17-overdue-check-night-merges-audit.md`
section 4, and the `FINDING`, `LESSON`, `GUARD FALSE POSITIVE` and `IDEA`
lines of the main checkout's `.claude/current-task/session.md`.

Acceptance: the harness suites (`pytest -m harness --no-cov`, whose pass count
is read, not the coverage floor) and the live checks of section 5, recorded
in this file. No `@perf` run (no refine or mesh code).

Words used here. **Grant**: a branch Ola named with `away.py --grant`; the
guard lets its push, pull request and enqueue through while he is away
(`tools/push_grant.py`). **Governed**: a file `guard_governance.py` asks
about; a change to it needs Ola at the prompt, so its green step runs by day.
**Bare command**: a line holding one simple command with every word written
out (no variable, no `&&`, `;`, pipe, redirect or heredoc). **Panel entry**:
a line in Ola's agent list in the Claude Code window.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature* (each read 2026-10-18).

- GraphQL specification, October 2021, section 2.3 *Operations* (item 5):
  "There are three types of operations that GraphQL models: query – a
  read-only fetch. mutation – a write followed by a fetch. subscription – a
  long-lived request that fetches data in response to source events." And
  the shorthand: such an operation "may be represented in a short-hand form
  which omits the query keyword and operation name". So a document whose
  operations all start with `{` or `query` reads; the guard passes only
  that. Departure: the spec allows `mutation` as a field name or inside a
  string; the guard asks on the letters anywhere in the text, comments and
  strings included, which fails closed.
- `gh api` manual (cli.github.com/manual/gh_api): "The default HTTP request
  method is `GET` normally and `POST` if any parameters were added"; "`-f/
  --raw-field` … static string parameters"; for `-F/--field`, "if the value
  starts with `@`, the rest of the value is interpreted as a filename to read
  the value from"; "For GraphQL requests, all fields other than `query` and
  `operationName` are interpreted as GraphQL variables"; a body "may be read
  from file specified by `--input`". So a GraphQL read is a POST, which is
  why today's rule (a field means a write) refuses it; and a query from a
  file or `--input` cannot be read by the guard, so those still ask.
- Python, `open()`: mode defaults to `'r'`, "open for reading"; and
  "A custom opener can be used by passing a callable as *opener*" (`help(open)`),
  so it runs any code, with the file name as its argument. Item 5's
  Python exemption rests on the first and refuses the second.

*Legacy.* Nothing: the legacy tree is the old meshing product, with no harness.

```
$ git grep -l -i -e 'push_grant' -e 'session_state' -e 'scratch_copy' -e 'away.py' -e 'graphql' legacy-archive -- legacy
(no output, exit 1)
```

*Earlier increments.* h23 (`tools/duties.py`, branch `worktree-h23-duties`,
design approved, not built) overlaps items 8 and 9; section 6 says how they fit.

## 2. The ideas file

`../rasputin_scratch/ideas/ideas.md`, read 2026-10-18.

- Taken up in part: *Monitor in every persona's tool list* (its second
  option, "make `.claude/briefs/common.md` name the fallback", is item 6;
  the first, adding `Monitor` to six `tools:` lines, is left out: whether a
  subagent gets the tool at all is unchecked, and item 6 needs no new tool).
- Overlapping, left out: *Morning check proposals P1-P12* (its P8,
  `session.md`'s format enforced, is items 4 and 8 here; the rest are not
  checked); *Recap: items in flight and open ideas* (shares the recap's
  10,000-character cap with items 4 and 8); *Clear out stale worktrees*
  (item 11 is about processes and panel entries, not worktrees).

## 3. The items

Each: the change, the test that pins it (with its fixture's measured
figures), whether the files are governed, and the size in counted lines
(`CLAUDE.md` §2; tests and docstrings do not count).

### Item 1. Post-review `@perf` evidence counts like prose

*Why.* `docs/increments/README.md` runs `@perf` after the code review, so
its evidence commit always follows the verdict, and `changes`
(`tools/push_grant.py@1806288c:227`) allows only prose after a review. Two
increments in one day could not be pushed under a grant: 66 (`5eee02e8`)
and 67 part 2 (`d1188d5c`). Each commit adds 13 files and changes nothing
else, `5eee02e8` under `docs/benchmarks/2026-10-17/`, `d1188d5c` under
`docs/benchmarks/2026-10-18/` (`git show --name-status --format= <sha>`: 13
lines each, all `A`; 4 of the 26 are `run.json`).

*Change.* In `changes`, a path also passes when it matches
`BENCH = re.compile(r"docs/benchmarks/\d{4}-\d{2}-\d{2}/.+")` and the commit
adds it (`git diff --name-only --no-renames --diff-filter=A -z base commit`,
one more git call). Added only, so earlier evidence cannot be edited; dated
only, so `docs/benchmarks/quick/` (`tools/bench_quick.py@1806288c:32`) stays
code. Nothing runs the dated folders' scripts: `testpaths =
["tests/python"]` and ruff excludes `docs/benchmarks` (`pyproject.toml@1806288c:109`,
`:134`), and CI still runs on them.

*What `tools/bench.py` reads there.* Two things. `DEFAULT_QUARTER`
(`tools/bench.py@1806288c:46`) is an existing file under
`docs/benchmarks/2026-09-26/`; add-only means no commit under the grant can
change it. And `find_baseline` (`tools/bench.py@1806288c:410-428`) takes the
newest `run.json` at any depth under `docs/benchmarks` (`rglob("run.json")`)
that is older than the new run, comparable with it, and whose commit is an
ancestor of the new run's. So an added `run.json` becomes the baseline of
later `bench.py` runs on commits that descend from it: a forged one could
hide a regression in a later increment. Item 1 accepts this and does not
exclude `run.json`: every `@perf` evidence commit carries one (4 of the 26
files above), so excluding it would leave item 1 useless; the same commit's
`README.md` carries the same figures and already passes as prose today; and
add-only keeps every existing baseline as it was. The risk is to a later
measurement, which `@perf` reads and `@reviewer` checks then, not to the
remote or to a rule.

*Tests* (`tests/python/test_push_grant.py`, in the existing `built` repo
fixture style): a commit after the approved review that adds three files
(`.tsv`, `.json`, `.py`) under `docs/benchmarks/2026-10-17/x/` passes; one
that adds a `run.json` there passes (the accepted risk above); one
that edits an existing file there is refused naming it; one that adds
`docs/benchmarks/quick/new.json` is refused; one that adds a dated file and
also changes `tools/x.py` is refused naming `tools/x.py`.

*Governed* (`tools/push_grant.py`). About 5 lines.

### Item 2. A grant refusal for a line that is not one bare command says so

*Decided*: yes (Ola, 2026-10-18, "yes to both"; `@orchestrator`'s question,
retrospective 2026-10-18 section 7.3).

*Why.* At 00:22:54 the main session wrote a PR body and ran `gh pr create`
for granted branch `worktree-pytest-home` in one line (`S=… && cat >
$S/pr68.md <<'EOF' … EOF` then `gh pr create … --body-file $S/pr68.md`).
`granted_act` takes only a bare command, so the line fell through to the
plain refusal ("gh pr changes a pull request"), which does not say why; the
PR waited 6 h 17 min.

*Change.* In `judge`, when `granted_act(command)` is None: for each simple
command of `shell_scan.parse(command)` that is `git push`, `gh pr create` or
`gh pr merge` and names a granted branch (a word equal to it, its `--head`
value, or a `<sha>:refs/heads/<branch>` ref), return
`Grant(False, branch, kind, None, NOT_BARE)` with
`NOT_BARE = "the grant takes only one bare command, every word written out
(no variable, no &&, ;, pipe, redirect or heredoc); run it alone"`. The guard
already prints "granted branch X, but <why>" and queues the refusal. The
detection is loose on purpose: it only chooses the refusal's words, never
lets anything through; the rerun is judged by `granted_act` as today.

*Rule* (`.claude/REQUIRED-READING.md`, *The harness*, its Unattended mode paragraph, after "do not retry
it or work around it"): "except a refusal that says the grant takes only one
bare command: run that command once more, alone and written out" (about +20
words). **The sentence goes in PR C, the rule-text PR, not PR A.**
`@developer`, who writes PR A's green step, cannot write a rule file
(`tools/brief.py@1806288c:39-40` gives it code folders and build files,
no Markdown; its green commit `188089f1` says so), and PR C already edits that file, by a writer who can.
No order between A and C is needed: until C lands, an agent that meets the
`NOT_BARE` refusal follows today's rule (no retry; an `ASK OLA:` line), which
is today's behaviour; and an agent that does retry gains nothing unsafe,
because the rerun is judged afresh by `granted_act`.

*Tests:* the 00:22:54 line (body cut to one line), branch granted, mode on:
refused, reason contains `NOT_BARE`; the same line with the branch not
granted: the plain refusal, no `NOT_BARE`; `git push origin $SHA:refs/heads/<granted>`
(a variable): `NOT_BARE`; the bare `gh pr create --head <granted> --title 'T'
--body-file /abs/pr.md` with an approved head: passes (regression).

*Governed* (`tools/push_grant.py`; the sentence, `.claude/REQUIRED-READING.md`, in PR C). About 12 lines.

### Item 3. A fix that lives under `tests/` is `@tester`'s green step

*Decided*: yes (Ola, 2026-10-18, "yes to both"; retrospective 2026-10-18
section 7.4).

*Why.* 68's whole fix was `tests/python/conftest.py`; `tools/brief.py` gives
`tests/` to `@tester` alone, so the design's "green by `@developer`" had no
writer, and the green commit broke "The green commit touches **no test
file**" (`docs/increments/README.md@1806288c:45-46`) while keeping its purpose.

*Rule* (`docs/increments/README.md`, step 3, after the sentence above): "A fix
that lives under `tests/` (a conftest hook, a shared fixture) is `@tester`'s
green step, in its own commit after the red one, leaving the red commit's
files untouched; `@reviewer` checks both." The short path's "never dropped"
list changes with item 10's cut C6. About +30 words. No test (prose); the
check is `@reviewer`'s, by `git show --stat` of the two commits.

*Governed* (`docs/increments/README.md`). 0 lines.

### Item 4. The recap warns, first, when it passes its cap

*Why.* `python3 tools/session_state.py | wc -m`: 28,331 characters on
2026-10-17 (retrospective section 5 item 5), 38,650 at 07:15 UTC today; the
cap is 10,000 (`.claude/REQUIRED-READING.md`, end of *The harness*). A warning
at the end would be cut off: the `SessionStart` hook's output is capped at
10,000 characters (Claude Code hooks guide, as h23 section 1 quotes it).

*Change.* `main` renders the whole output into a buffer
(`contextlib.redirect_stdout`), then prints, before it, one line when the
length passes `CAP = 10_000`: `RECAP OVER CAP: <n> characters, cap 10,000;
session.md is <m> of them; prune it to its NOW, QUEUE and ASK OLA lines.`
A pure `cap_warning(length: int, session_chars: int) -> str | None` makes
the line. `CAP` assumes the hook's limit above; checked at 38,650.

*Tests* (`tests/python/test_session_state.py`): `cap_warning(10_000, 0)` is
None; `cap_warning(28_331, 25_000)` starts `RECAP OVER CAP: 28,331` and names
25,000; and one run of `main` with a 12,000-character `session.md` in a tmp
tasks folder (the existing monkeypatch style) prints that line first.

*Governed* (`tools/session_state.py`). About 10 lines.

### Item 5. Two guard false positives: GraphQL reads, Python reads

Six refusals in two days (`session.md` lines 47, 74, 75, 83, 112, 114).
Reproduced on `1806288c`:
`judge_bash("python3 -c 'import json; print(json.load(open(\"/x/.git/harness/unattended.json\")))'")`
gives `('deny', [], '')`; with `.claude/settings.json` it gives an ask; and
`segment_why(["gh","api","graphql","-f","query={ viewer { login } }"])` gives
"gh api with a writing method changes the forge".

**5a. `gh api graphql` reads** (`.claude/hooks/guard_push.py`, `segment_why`).
The guard reads the words before the shell expands them, so a `$` or a
backtick hides what gh will receive. A `gh api` segment passes as a read
only when all of these hold; otherwise the present rule applies (a field
means a write, so it asks):

- no word of the segment holds `$` or a backtick (so `query($n:Int!)` asks:
  a query with variables cannot be told from one with an expansion);
- the first word after `api` is `graphql`, and every later word is a flag
  from `-f`, `--raw-field`, `-F`, `--field`, `-X`, `--method`, `-q`,
  `--jq`, each followed by its value as a separate word. Any other word asks:
  a second plain word, `--input`, `-H`, `--hostname`, `--paginate`, a flag
  with its value glued on;
- the method, if given, is `POST`;
- the fields are exactly one `query` and at most one `operationName`, no
  others (no variables); a `-F` value does not start with `@`;
- `read_query(text)` holds for the `query` value: after leading whitespace
  it starts with `{` or the word `query`, and the letters `mutation` and
  `subscription` appear nowhere in the text as given. Nothing is stripped
  first: a `#` inside a string is not a comment, so a document with a
  leading comment asks.

Tests (`tests/python/test_guard_push.py`, `segment_why` on split words).
Passes: the 23:03:56 read
(`query='{ repository(owner:"expertanalytics", name:"rasputin") { pullRequest(number:285) { userContentEdits(first:5) { nodes { editedAt } } } } }'`);
`query='query A { viewer { login } }' -f operationName=A`; `-X POST` with
`{ viewer { login } }`. Asks, each a write or unreadable:

- the review's three: `query='query A { a(x:"#") } mutation B { addComment(…) }' -f operationName=B`;
  `-f query="query A {a} $M" -f operationName=B` (the word holds `$M`),
  and `gh api graphql $X -f query='{ a }'`; `gh api graphql repos/a/b/issues -f query='{ a }'`
  (a second plain word);
- `query='query($n:Int!){ a }' -F n=285` (a `$`), `mutation { addComment(…) }`,
  `query { a } mutation { b }`, `subscription { a }`, `-F query=@q.graphql`,
  `--input q.json`, `-X PUT`, a second `-f query=`, `-f n=1` beside a
  query, `-H 'X-HTTP-Method-Override: PUT'`, and REST `gh api repos/a/b/issues
  -f title=x` (unchanged).

**5b. Python reads** (`tools/shell_scan.py`, `candidates`). A string
literal is not a write candidate when it sits in one of three exact shapes
and the whole program meets four conditions. `OPEN` below is `open(LIT)` or
`open(LIT, M)`, `LIT` a string literal, `M` a string literal in `r`, `rb`,
`rt`, and no keywords but `encoding` and `errors`, each a string literal
(so no `opener`, and no `mode=`).

- Shapes: `json.load(OPEN)` with no other argument; `OPEN.read()` or
  `OPEN.readlines()` with no argument; `Path(LIT).read_text()` or
  `.read_bytes()` (or `pathlib.Path(LIT)…`), `Path` with that one argument,
  the read with no arguments but literal `encoding` and `errors`. The file
  or path object is never named, so nothing later can reach it
  (`open(LIT).name` reused by a later `open(…, "w")` stays a candidate).
- Imports: only `import json`, `import pathlib` and `from pathlib import
  Path`, none with `as`; no other import and no `import *`. One allowlist
  shuts out `sys`, `inspect`, `gc`, `ctypes`, `importlib` and `builtins`
  together, where a list of forbidden modules would keep growing.
- No rebinding: the names `open`, `json`, `pathlib` and `Path` are bound
  nowhere else: no assignment or `del` target (plain, augmented, annotated,
  walrus, `for`, `with … as`, `except … as`, comprehension), no function
  parameter, `def`, `class`, `global` or `nonlocal` of them.
- No patching or reflection: no assignment, augmented assignment or `del`
  whose target is an attribute or a subscript (`builtins.open = …`,
  `Path.read_text = …`); none of the names `setattr`, `delattr`, `getattr`,
  `vars`, `locals`, `globals`, `exec`, `eval`, `compile`, `breakpoint`,
  `__import__`, `__builtins__`, `builtins`; and no name or attribute that
  starts with `_` (the routes to the literal object through `__code__`,
  frames or `__class__`).

If any condition fails, every literal is a candidate, as today. The literal
then flows only into the read and is gone: no argument runs code, nothing
replaces `open`, `json.load` or `Path.read_text`, and nothing reaches the
literal through code objects. The scanner stays a tripwire: a path built
from pieces is no candidate today, and the exemption leaves that as it is.
Reads written as `for line in open(LIT)` or `with open(LIT) as f` stay
candidates; write `open(LIT).readlines()`.

Tests (`tests/python/test_guard_targets.py` for `candidates`,
`tests/python/test_guard_governance.py` for the governance guard; `H` is
the absolute path of the harness folder's `unattended.json`, `W` any
function). Not candidates:
`import json; print(json.load(open(H)))` by `-c` and by heredoc
(`python3 - <<'X'`); `from pathlib import Path; print(Path(H).read_text())`;
`print(len(open(H, "r", encoding="utf-8").readlines()))`;
`open(".claude/settings.json").read()` passes the governance guard.
Still candidates, each its own test:

- the review's four: `n = open(H).name; open(n, "w").write("x")`;
  `open(H, opener=W).read()`; `import builtins; builtins.open = W;
  print(open(H).read())`, and `from pathlib import Path; Path.read_text = W;
  Path(H).read_text()`, and `setattr(Path, "read_text", W)` before the
  read; `import shutil as json; json.load(open(H))`, `from x import Path;
  Path(H).read_text()`, `from os import *; print(open(H).read())`, and
  `json = W; json.load(open(H))`;
- the reflection route: `print(open(H).read()); (lambda: 0).__code__`, and
  `import sys; print(open(H).read())`;
- unchanged: `open(p, "w")`, `open(p, "a")`, `open(p, mode)`,
  `open(H, mode="r").read()`, `Path(p).write_text("x")`, `Path(p).open("w")`,
  `p = Path(LIT); p.read_text()`, `f = open(H); f.read()`,
  `for l in open(H): print(l)`, the harness folder's path joined with `/`,
  `shutil.copy(LIT, LIT)`, and a program with `def open(p): …` before
  `open(H).read()`.

**Code review round 1 on `188089f1`: what the shell expands.** Both
exemptions judge words the shell has not yet expanded, and the `$`/backtick
check of 5a covers only two of the shell's expansions. The review found
writes passing both; `@tester` pinned them in `3c227e08` (9 tests, all
through the whole hook, all red on `3c227e08`). Two more conditions:

- **5a: the line expands nothing.** `graphql_read` gains a third argument,
  `raw`, the raw command text, and passes only when `raw` matches, whole,
  `EXPANDS_NOTHING = re.compile(r"(?:[A-Za-z0-9_./= \t|-]|'[^']*')*")`, and
  no word of the gh segment starts with `=` (`any(w.startswith("=") for w in
  words)`; zsh expands `=cmd` to a path, also when written `''=ls`, whose
  word is `=ls`). So every character outside single quotes is a letter,
  digit, `_ . / = -`, a blank or a pipe: no brace, glob, tilde, `$`,
  backtick, double quote, backslash, redirect, `;` or `&` anywhere in the
  line, and the scanner's words are then exactly what gh receives. This is
  the review's rule (plain `[A-Za-z0-9_.=-]` or single-quoted), widened by
  `/` (a path to `gh` or in a pipe) and the pipe (`| jq …`, `| head`),
  neither of which the shell expands. *Linear time:* each step of the
  pattern takes one character or one whole quoted string, and `'` is not in
  the character class, so there is one way to match any prefix and no
  backtracking. (Design review round 3: `a4c19f4e`'s pattern nested two
  repeats; on a failing line its time grew about fourfold per two
  characters, 5.8 s at 36 characters, and a hook that times out lets the
  call through.) It reads the raw text the hook already has:
  `segment_why(words, text)` receives the whole line from `main`
  (`.claude/hooks/guard_push.py@3c227e08:360`), and the call at `:261`
  passes `text or shlex.join(words)`. The fallback keeps the existing
  word-level tests (`segment_why(words)`, `text` empty): `shlex.join`
  single-quotes every word that holds a character outside its safe set, so
  a split read like `query={ viewer { login } }` rejoins as
  `'query={ viewer { login } }'` and still passes; in production `main`
  always passes the line, never an empty one. The check is on the whole
  line, not the gh segment: the scanner keeps no word spans, and matching a
  word to a quoted copy of it elsewhere on the line would let
  `echo '{x,-X,PUT}'; gh api graphql … -q {x,-X,PUT}` through. Cost: a read
  behind `cd /x &&`, `;`, a redirect or in double quotes asks, as on master.
- **5b: no exemption when the shell rewrites the program.** In `candidates`
  (`tools/shell_scan.py@3c227e08:454`), `exempt = set()` when the simple
  command is `unknown`, which item 14 sets exactly when the shell may
  rewrite the program: a `-c`/`-e` program holding `$` or a backtick (new
  two-line helper `substituted(text) -> bool`), or a heredoc program whose
  delimiter is unquoted and whose body holds `$`, a backtick or a backslash.
  The program text needs no new reading of the raw line for `$` and
  backticks: it keeps them as written, since `_Scan.dollar` returns
  `text[start : self.at]` (`:233`) and a heredoc body is kept line for
  line. The backslash matters only in a heredoc: in an unquoted one the
  shell turns `\\` into `\` before Python reads it, while the scanner keeps
  the body as written, so `x = 'a\\'; y = 1 #'; open(…,'w').write('X')`
  shows the scanner a comment where Python runs a write (design review
  round 3; also on master). In a double-quoted `-c` the scanner already
  removes the escapes as the shell does (`_Scan.quoted`). A `-c` program is
  not told single- from double-quoted: one holding `$` loses the exemption
  either way, master's behaviour for every literal. A quoted heredoc
  (`<<'EOF'`, `<<"EOF"`, `<<\EOF`) expands nothing and keeps the exemption,
  `$` and backslashes included.

Tests added by round 3 (`tests/python/test_guard_push.py`): through
`segment_why` with the line as text, `gh api graphql -f query='{ a }' -q
''=ls` names the forge; and a failing line of at least 60 characters,
`gh api graphql -f query='{ a }' --jq .data.repository.id > /tmp/out.json`
(72), names the forge in under 1 s, timed around the call (on `a4c19f4e`'s
pattern a 36-character failing line took 5.8 s; on the new one, under
0.001 s). 5b's backslash tests are item 14's.

Measured, both fixes and item 14 as revised here, applied to a scratch copy
of `3c227e08` (`tools/scratch_copy.py`): `pytest -m harness --no-cov`
2530 passed, 14 skipped (the scratch copy is no git work tree), 0 failed:
the 9 red tests pass and every earlier test, the word-level 5a and 5b ones
included, still does. Extra probes through `segment_why` with the line:
`gh api graphql -f query='{ a }' | head -5` and `… --jq .data | jq -r '.x'`
pass; `echo '{x,-X,PUT}'; gh api graphql -f query='{ a }' -q {x,-X,PUT}`,
`-q =x`, `-q ''=ls`, `-q ~`, `query="{ a }"`, `cd /x && gh api graphql -f
query='{ a }'` and the 72-character redirect line ask. Through `judge_bash`,
`python3 -c 'print(open(".claude/settings.json").read().count("$"))'` asks
(silent on `3c227e08`, asks on master): the accepted cost above; the same
read by quoted heredoc stays silent.

*Governed* (both files). Green `188089f1`: 98 counted lines
(`python3 tools/count_loc.py 61d3ce4b 188089f1`: `guard_push.py` 21,
`push_grant.py` 20, `shell_scan.py` 57), against the 72 estimated. The round-1
and round-3 fix adds about 3 (5a) and 0 (5b, one changed line).

### Item 6. Briefs name a wait that personas can run

*Why.* `.claude/briefs/common.md@1806288c:12` says "Wait for a background run with the
Monitor tool, never `sleep`"; no persona's `tools:` line has Monitor
(`grep -n '^tools:' .claude/agents/*.md`). Personas waited with `tail -f`,
which never exits (8 left over, up to 2 days old, on 2026-10-17), and
`@perf`'s `pgrep`/`kill -0` loops reported done while `bench.py` still ran.

*Change.* The sentence becomes: "Wait for a background run by its exit
file: start it with Bash `run_in_background` as `<cmd> > <log> 2>&1; echo $?
> <log>.exit`, then wait in one Bash call with `timeout: 600000`, `for i in
$(seq 55); do test -e <log>.exit && break; sleep 10; done`, repeated as
needed; never `tail -f`, never a wait on a process id." The timeout is
needed: the Bash tool's default is 120 s, and the loop runs up to 550 s,
under the tool's 600 s maximum. The exit file is
written only when the command ends, so it cannot report done early. A
bounded loop with `sleep` runs in a persona's Bash (checked by this design's
own run). About +40 words; `.claude/briefs/common.md` is counted by
`tools/rule_sizes.py`.

*Tests:* `tests/python/test_brief.py`, one assertion that the rendered
block contains `.exit` and not `Monitor`.

*Governed* (`.claude/briefs/`). 0 lines.

### Item 7. `tools/scratch_copy.py` refuses a worktree with no venv

*Why.* With no `<worktree>/.venv`, it warns "no built extension" and still
prints a command naming `<worktree>/.venv/bin/python`, which does not exist
(`@tester`, 68).

*Change.* Before `git archive`, when `ROOT/.venv/bin/python` is missing,
refuse (exit 2, nothing extracted): `scratch_copy: <ROOT> has no .venv; make
one with this worktree's install before copying`.

*Tests* (`tests/python/test_scratch_copy.py`): ROOT pointed at a tmp repo
with no `.venv`: exit 2, the stderr line, the target not created.

*Not governed.* About 3 lines.

### Item 8. The queue carries a status per item, checked before a night

*Why.* At 00:33 the main session found no unblocked work, but item 1 of
`session.md` line 24 (the private routing work) had been unblocked since
#279 merged; no item on that line was marked done, parked or blocked
(retrospective section 7.2). The window idled 3 h 15 min.

*Format* (`.claude/REQUIRED-READING.md`, the `session.md` bullet): "`QUEUE:`
(one) lists items split by `; `, each opening with its status in brackets:
`next`, `running`, `waits-ola`, `blocked`, `parked` or `done`, and
`fallback` after it for an item that needs no ruling and writes no governed
path (`[next fallback]`)." The *Unattended mode* sentence on the fallback
gains "`away.py` warns when none is open". About +45 words.

*Change.* `session_state.py`: `queue_items(text) -> list[tuple[str, str]]`
(status, text) from the `QUEUE:` line, pure; `session_format` adds a warning
per item with no or an unknown status, and one line `<k> done items: prune
them`. The recap prints `Queue:` and each item, cut to 100 characters, after
the `session.md` format lines. `away.py <duration>` (entering) prints the
queue the same way and, when no item is `[next fallback]` or `[running
fallback]`, prints `No open fallback item: the window may idle.` It warns
only; Ola decides.

*Tests:* `queue_items` on `QUEUE: [next] a; [waits-ola] b; [next fallback] c`
gives three pairs; `session_format` warns on `QUEUE: a; [done] b` (one
missing status, one done item); `away.py` entering with a queue holding no
open fallback prints the warning, and with one prints none (the existing
`test_away.py` tmp-state fixture, `session.md` in a tmp tasks folder).

*Governed* (`tools/session_state.py`, `tools/away.py`, `.claude/REQUIRED-READING.md`).
About 25 lines.

### Item 9. `away.py --back` prints the window's timeline

*Why.* Each morning's idle account so far came from the main session's or
`@orchestrator`'s reading of transcripts; Ola: "You seem to forget these
after a while". `--back` prints only the longest stretch without a commit.

*Change.* After that line, `--back` prints:

```
Timeline, UTC (idle: no commit and no session or agent event for 15 min or more):
  22:59-00:32  working  1 h 33 min
  00:32-03:47  idle     3 h 15 min
Idle: 3 h 15 min of 4 h 48 min (68 %). Agent runs: 9 (architect 2, developer 2, perf 1, reviewer 3, tester 1).
```

- *Events*: every commit time inside the window (`commits_between`, already
  there) and every `"timestamp"` inside the window in a transcript file under
  `session_state.TRANSCRIPTS` (main sessions and their `subagents/`), read
  only from files modified after the window's start, by a regex on each
  line, never a full JSON parse.
- *Agent runs*: subagent transcripts with at least one event inside the
  window, counted by `agentType` from the `.meta.json` beside each.
- Pure `timeline(since, end, moments, gap) -> list[Stretch]` with
  `Stretch(start, end, idle)`: a gap of at least `GAP` between neighbouring
  events (the window's start and end included) is idle, the rest working.
  At most 12 lines, then `… <n> more`; a fault prints `Timeline: unknown
  (<reason>)` and `--back` goes on, as the stretch line does.
- `GAP = 15 min`. Assumes no tool call runs longer than 10 minutes (the
  Bash tool's longest timeout), so a working agent writes an event at least
  every 10 minutes; item 6's wait loop keeps to that. A background job that
  nobody waits on writes no event and is counted idle; that is a limit, said
  in the docstring. Checked on one window, in section 5.

*Tests* (`tests/python/test_away.py`): `timeline` on a 4 h 48 min window
(22:59:38 to 03:47:38) with events every 5 min until 00:32 and none after
gives two stretches, the idle one 3 h 15 min; events 14 min apart give one
working stretch; no events gives one idle stretch; and `--back` with a tmp
transcripts folder holding one subagent file and its `.meta.json` prints the
`Agent runs: 1 (tester 1)` line.

*Governed* (`tools/away.py`). About 45 lines.

### Item 10. The rule-text cut

*What is still owed.* The "−347 words" carried since 2026-10-14 (its P9) is
five proposals; one is done and three describe the same sentence:

| Proposal | State on `1806288c` |
|---|---|
| −35, `reviewer.md`'s verdict item | done: `.claude/agents/reviewer.md@1806288c:55` is the short form |
| −155 of 2026-10-12: pipeline twice in `CLAUDE.md` (−40), branch-protection history (−45), governance examples (−70) | not done |
| −67 of 2026-10-13: the same governance sentence | same sentence as the −70 |
| −30, the gate-planting recipe | not done |
| −60, the docstring rule's examples | not done |

The 2026-10-17 audit's −90 (the hook-by-hook account) and the 2026-10-18
retrospective's −45 also cover that sentence. So the distinct cuts are C1 to
C6 below, measured with a scratch script (old span and new text, words by
whitespace as `tools/rule_sizes.py` counts): **−229 words**, not −437.

| Cut | File | Old → new | Words |
|---|---|---|---|
| C1 | `.claude/REQUIRED-READING.md`, *The harness*, "Active in …" to "unapproved push through." | below | 267 → 156, −111 |
| C2 | same, *Claims*: "Verify a gate by planting … restore the file." | "Verify a gate by planting what it forbids." | −17 |
| C3 | `CLAUDE.md` §3, "Every code alteration must execute … before merge readiness." (the numbered list) | "Every code alteration follows the step order below, dispatched by the main session." | −44 |
| C4 | `CLAUDE.md` §4 *CI*, "prose-only change; it is meant … never report." | "prose-only change, and it is the one required check on `master`." | −28 |
| C5 | `.claude/skills/python-development/SKILL.md`, "* They never carry history: … keep that." | "* They never carry history: no increment, design-file, review or PR labels, no record of who decided, no commit hashes, no account of what the code used to do. Git and `docs/increments/` keep that." | −21 |
| C6 | `docs/increments/README.md`, short path, "Never dropped: … governed file." | "Never dropped: the red commit ahead of the green one, a green commit that leaves the red suite untouched, and the approval boundary of `.claude/REQUIRED-READING.md`." | −8 |

C4 also corrects a stale claim: `gh api
repos/expertanalytics/rasputin/branches/master/protection/required_status_checks
--jq .contexts` gives `["CI result"]`, so it is already the one required check.
C5 keeps the tool-directive bullet (the 2026-10-13 proposal cut it): a
writer who takes `# type: ignore[...]` for a comment would delete it.

C1's new text:

> Active in `.claude/settings.json`; each hook's docstring lists the forms it
> judges. `guard_push.py` asks before any act that publishes, writes to the
> forge, rewrites history or writes refs, remotes or git config, and before
> `--no-verify`; `guard_governance.py` asks before any write to a file that
> states rules or to a `tools/` file named after a standard-library module,
> however the path is spelled, except a scratchpad path written by Edit,
> Write or one plain command; `guard_unattended.py` refuses
> `AskUserQuestion`, permission prompts and configuration changes while
> unattended mode is on; `guard_spawn.py` refuses a persona spawn without an
> unchanged, current block from `python3 tools/brief.py`, or from outside the
> session's starting directory, and a resume that starts a new step carries a
> fresh block; `gates_after_commit.py` puts the gates' output in the
> transcript after a commit or merge, and exits 2 when one is red. These are
> tripwires: the boundary is still yours to keep, and auto mode can let an
> unapproved push through.

What C1 moves into the hook docstrings (no history labels:
`tools/check_code_history.py` checks them): to `guard_push.py`, "not `gh pr
close` or `gh pr comment`", "no git write in a scratchpad repository", and
the parsed-or-text sentence (already there); to `guard_governance.py`, the
forms "a copy or move into a governed directory (`cp json.py tools`, with or
without the trailing `/`)", "a path judged with `.`, `..` and a doubled `/`
resolved (`tools/./json.py`)", "not a write to a whole governed directory
(`rm -r .claude/hooks`)", and "one plain command: no `&&`, `;`, pipe,
substitution or interpreter program". Docstrings are not counted lines.

With items 2, 3, 6, 8, 11 and 12 (+20, +30, +40, +45, +30, +10) and section
4's review-line wording (+12), the net is about **−43 words**. Measure on the branch with `python3 tools/rule_sizes.py`
and record the figure in this file.

*Governed* (all four rule files and both hooks). 0 lines.

### Item 11. Leftover processes and stale panel entries

*Where the rule lives.* Nowhere yet (`grep -n -i 'panel\|leftover\|TaskStop'
CLAUDE.md .claude/REQUIRED-READING.md` finds nothing). New bullet in
`CLAUDE.md` §3 *The main session dispatches*: "**After each handback.** Stop
the agent's panel entry (TaskStop) and any shell it left running: the
recap's *Running Bash-tool jobs* lists them with their age." About +30 words.

*Can a hook do it?* Not the panel entry: TaskStop is a tool the model calls,
and a hook only runs a command. Not the shells, safely: a `SubagentStop`
hook gets no list of the shells that agent started, and the recap's job list
(`session_state.py`, `running_jobs`) cannot tell whose a shell is, so a hook
would have to guess what to kill. The cause goes with item 6: a wait that
ends on its own leaves no shell. If leftovers persist after item 6, the next
retrospective counts them and proposes a hook with evidence. As the main
session told Ola: the three-item cap counts work items, not panel entries,
so stale entries cost clarity and leftover shells, not spawn slots.

*Governed* (`CLAUDE.md`). 0 lines; no test (prose).

### Item 12. Reviewer wordings: the split Ola kept

Ola, 2026-10-18 06:44 UTC: "yes, keep the split". No rule says it:
`docs/increments/README.md@1806288c:61` says only "The main session applies
them exactly", and `.claude/agents/orchestrator.md` §3 says the
retrospectives are `@orchestrator`'s. Add to README step 4, after that
sentence: "On a retrospective, `@orchestrator` applies them instead." About
+10 words.

*Governed* (`docs/increments/README.md`). 0 lines.

### Item 13. A recorded signature for every name `plugin_api` promises

Ola, 2026-10-18: "add the plugin_api signature test to h24".

*Why.* `src_python/tin_engine/plugin_api.py`'s docstring promises each name
is "removed or changed only in a release that says so", and
`tests/python/test_plugin_api.py` checks only the list and that each name is
the home module's object. Two promised functions changed shape with no
public test failing: `destination(out, out_parent, name, hint)` became
`destination(out, hint)` (`7d6e9977`, increment 57), and `repository_for`
returns the repository alone, not a pair (`81838dcf`, increment 60). Moving
rasputin-lsm's pin past them broke 28 private tests.

*Where the record lives.* In `tests/python/test_plugin_api.py` itself, as a
literal dict `SIGNATURES: dict[str, str]` beside `PROMISED`, written out by
hand, never generated by a flag: a promise change is then a visible edit of
a test file, which only `@tester` writes. One test,
`test_each_promised_name_keeps_its_recorded_signature`, parametrized over
`__all__`, compares the record with:

- a function: `str(inspect.signature(obj))`, which holds the parameters,
  their kinds and defaults, and the return annotation. The modules use
  `from __future__ import annotations`, so annotations are the source text
  (`-> 'Path'`). The return shape is part of the promise this way:
  `repository_for`'s `-> 'tuple[…]'` against `-> 'TiffDemRepository …'`
  fails. Same text on CI's three Pythons: the 15 functions' `def` lines,
  bodies cut, run through `str(inspect.signature(…))` under 3.12.12,
  3.13.15 and 3.14.7 (`uv run --no-project --python <v>`), print the same
  15 lines as the main checkout's `.venv`, byte for byte (`diff`, empty).
- a class: all five (`Bounds`, `DemRequest`, `DemTile`, `RasterMeta`,
  `Resources`) are Pydantic models, and their `inspect.signature` is not
  source text: it prints `tin_engine.io.models.Bounds` and
  `_PydanticGeneralMetadata(...)`, which Pydantic may change. So a class is
  recorded as its fields from `model_fields`, in order, joined by `, `: the
  name alone when required, `name=<repr of the default>` otherwise (no
  field uses a default factory), e.g. `RasterMeta`'s `…, epsg, crs='',
  geographic=False, …, z_unit='metre', file_crs=None`. A renamed, added,
  dropped or newly required field, or a changed default, fails; a changed
  field type does not (a limit, said in the test's docstring).
- `THREADS`, an `Annotated[...]` alias (its `inspect.signature` is a
  meaningless `(*args, **kwargs)`): recorded as its type, `str` of the first
  `Annotated` argument (`int | None`, the same under all three Pythons), and
  the Typer option's flag, `min` and `metavar`: `int | None, --threads,
  min=0, metavar=N`. Typer stores the flag as the option's `default`, and
  its `param_decls` is empty (`()`), so `param_decls` is not the field to read.

A missing or extra key in `SIGNATURES` fails too, so a name added to
`__all__` needs its record in the same commit. The failure message gives the
name, the recorded and the current text, and the way out: "a promised
signature changed: plugins (rasputin-lsm among them) break. If on purpose,
update SIGNATURES in its own commit and say so in the PR body."

*A deliberate change.* `@tester` updates the entry in its own commit (never
inside `@developer`'s green commit; `docs/increments/README.md` step 3
already forbids that), and the PR body names the changed signature, so the
private repository's pin move knows what to adapt. The increment that
changes the signature says so in its design.

*Fixture, measured* on the main checkout's `.venv` at 07:20 UTC (master
`91ce6262`): 21 names in `__all__`: 15 functions, 5 classes, and `THREADS`;
`destination` reads
`(out: 'Path', hint: 'str' = '--out') -> 'Path'`. The record type is `SIGNATURES: dict[str, str]`
for all 21. *Probe able to fail:* the
suite run on a scratch copy of `cea32bad` (`tools/scratch_copy.py
cea32bad <dir>`, the old public pin) fails on exactly `destination` and
`repository_for`, recorded in this file.

*Not governed* (`tests/` only). 0 counted lines; about 60 test lines. No
red/green pair: the test pins today's behaviour and passes on arrival, so
the probe above stands in for the red step.

### Item 14. A Python program the shell rewrites is not readable

*Why.* Found by `@tester` beside PR A's review (`3c227e08`, its "Left out"),
on master too: a literal holding a space is skipped as prose
(`tools/shell_scan.py@3c227e08:461`), so in a double-quoted `-c` or an
unquoted heredoc, ``open('<harness path>`cat q`')`` or
`open('.claude/settings.json$(cat q)')` gives no candidate, and the shell may
splice any code into the program before Python reads it (q holding
`','w').write('x'))#` turns the read into a write). Measured through
`judge_bash` on a scratch copy of `3c227e08`: the backtick read of the
harness state by `-c` and by unquoted heredoc, the backtick and the
`$(cat q)` read of `.claude/settings.json`, and
`python3 -c "open('$(echo CLAUDE.md)','w').write('x')"` all give None
(silent).

*Change* (`tools/shell_scan.py` only; the guards are unchanged). Where the
program text holds `$` or a backtick (item 5's `substituted`), what the
interpreter runs is decided by the shell at run time, so the scanner treats
it as unknowable, in two places:

- In `candidates`, a Python string literal for which `substituted` holds is a
  candidate whatever spaces it holds (spaces mark prose only where nothing is
  substituted). The guard already judges a candidate by its static form
  (`static(t)`, and `static(t, ".git")` for the harness state), so
  ``'<H>`cat q`'`` is judged as `<H>` and denied, and
  `'.claude/settings.json$(cat q)'` asks naming it.
- The simple command is marked `unknown` ("writes files it does not name")
  when it is an interpreter, not a shell (a shell's program is parsed as
  commands already), and the shell may rewrite its program. For a `-c`/`-e`
  program: `substituted(program)`, set in `classify` after `interpreter(…)`
  (`:323`). For a heredoc program: its delimiter is unquoted and the body
  holds `$`, a backtick or a backslash (item 5b says why the backslash), set
  in `finish_heredoc` after `simple.program = simple.heredoc` (`:303`).
  Whether the delimiter is quoted is read from the raw text in `redirect`:
  `start, target = self.at, self.word()`, and the delimiter expands when
  `self.text[start : self.at]` holds no `'`, `"` or `\` (the shell's rule:
  any quoting of the delimiter turns expansion off). That flag rides as a
  third element of the `heredocs` and `pending` tuples into
  `finish_heredoc(simple, expands)`. The governance guard then does
  what it does for any unreadable target
  (`.claude/hooks/guard_governance.py@3c227e08:158`, `:168`): denies when the
  line names `.git/harness`, and otherwise, when no candidate matched, judges
  the line by its text rule, which asks when the line names a governed
  file. So `open('$(echo CLAUDE.md)','w')`, whose static form is empty, asks
  naming `CLAUDE.md`. A splice that reads its path from a file and names no
  governed file in the line still passes: the guard is a tripwire, as for a
  write target that is wholly a variable.

Ask or refuse: as for every other unreadable line, ask (deny for the harness
state); refusing every program with a `$` would ask on each regex ending in
`$`. The cost: an interpreter line whose program holds `$` or a backtick and
whose text names a governed file asks even when that name is only prose or a
comment. On master such a line passes; it is rare, and the ask names the
file.

*Tests* (`tests/python/test_guard_governance.py`, through the whole hook as in
`3c227e08`; `H` the harness state's absolute path, q holding item 5b's
`SPLICE`): ``python3 -c "print(open('<H>`cat q`').read())"`` and the same by
unquoted heredoc, modes off and on: denied with the always-denied text,
nothing queued; ``python3 -c "print(open('.claude/settings.json`cat q`').read())"``
and the `$(cat q)` form: ask, naming `.claude/settings.json`;
`python3 -c "open('$(echo CLAUDE.md)','w').write('x')"`: asks naming
`CLAUDE.md`. The backslash (design review round 3; `B(t)` is the body
`x = 'a\\'; y = 1 #'; open('t','w').write('X')`): `python3 - <<EOF`
with `B(CLAUDE.md)` asks naming `CLAUDE.md`; with `B(<H>)`, modes off and
on, denied, nothing queued; `python3 - <<'EOF'` with `B(CLAUDE.md)` stays
silent (quoted: Python too reads the write as a comment). Still silent:
`python3 -c 'print("a $ b")'`,
`python3 -c 'print("see CLAUDE.md for rules")'` (a prose literal, nothing
substituted), `python3 -c 'print(open("/x/notes.txt").read().count("$"))'`,
and by quoted heredoc `import json; print("\n".join(json.load(open(H))))`
and `print(open(H).read().count("$"))` (the read exemption holds).
In `tests/python/test_guard_targets.py`: `parse` of a double-quoted `-c`
program holding `$(cat q)` gives the `python3` simple with `unknown` True
(False on `3c227e08`); of `bash -c "echo $(cat q)"`, four simples, none
`unknown`; of `python3 - <<EOF` with `B(notes/plain.txt)`, `unknown` True,
and of the same with `<<'EOF'`, `<<"EOF"` and `<<\EOF`, False. Each line
above was run, through `judge_bash` or `parse`, on the patched scratch copy,
with the result stated; all of the round-3 heredoc lines are silent on
`3c227e08`.

*Governed* (`tools/shell_scan.py`). About 5 counted lines (the 2-line helper
is shared with 5b; 2 lines set `unknown`; 1 in `redirect`), and about 8
changed lines (`candidates`, and the heredoc tuples' types and unpacking).
In PR A, beside item 5b: same file, same green step.

## 4. The PRs

Total about 190 counted lines and about 60 test lines for item 13, far under 700; split by risk and by night
eligibility, not by size.

| PR | Items | Files | Counted lines | Risk | When |
|---|---|---|---|---|---|
| A, guards | 1, 2, 5, 14 | `tools/push_grant.py`, `.claude/hooks/guard_push.py`, `tools/shell_scan.py` | about 106 (98 in green `188089f1`, about 8 more for the review fixes and item 14) | a guard that lets a write through; review probes it adversarially | by day; red can be written any time, green needs Ola at the prompt |
| B, reports | 4, 8, 9 | `tools/session_state.py`, `tools/away.py`, item 8's sentences | about 80 | none to the boundary: they print | by day, as A |
| C, rule text | 2's sentence, 3, 6, 10, 11, 12, the review-line wording below | four rule files, `.claude/briefs/common.md`, two hook docstrings | 0 | rule wording only | by day (governed); no red or green, `@reviewer` once |
| D, scratch copy | 7 | `tools/scratch_copy.py` | about 3 | none | short path; not governed, so a night fallback |
| E, plugin signatures | 13 | `tests/python/test_plugin_api.py` | 0 | none | `@tester` alone, then `@reviewer`; not governed, so a night fallback |

Order: none needed among A, B and C (item 2 says why A does not wait for
its sentence in C). No added sentence falls inside a span C cuts:
item 2's sentence goes in the *Unattended mode* paragraph and item 8's in
the `session.md` bullet and that paragraph, while C1 cuts the "Active in"
paragraph, C2 the *Claims* gate recipe, and C3 to C6 other files. B and C
edit `.claude/REQUIRED-READING.md`, in different paragraphs; A no longer
does. Whichever lands later merges master first (`tools/merge_master.py`).
A and B are independent and can run as the two writers; C needs no writer
beyond its prose. D and E go any time.

**The review line.** Every `@reviewer` handback on these PRs, design and code
rounds alike, ends on one line written exactly
`Review: <VERDICT> <full sha>`: `<VERDICT>` one of `APPROVED`,
`APPROVED WITH WORDINGS`, `CHANGES REQUESTED`; the sha all 40 characters;
nothing else on the line (no bold, no backticks, no list dash); and no
other line starting `Review:`. The hook `record_review.py` keeps only such
a line (`REVIEW`, `tools/push_grant.py@1806288c:39-41`), and the away-mode
push grant reads only what it kept. The increment file's one-line record
(`- **Code review, round 1, …** On <sha>, @reviewer: …`) is the main
session's copy, written afterwards; a handback in that form is not recorded.
On 2026-10-18 the main session counted 4 of 19 verdicts unrecorded for that
reason (its count, not re-checked here). The rule already says this
(`docs/increments/README.md` step 4, "ends with exactly one line
`Review: <verdict> <sha>`"); PR C makes the wording exact there: `<sha>`
becomes "the full 40-character sha", and "written plain, with nothing else
on the line" is added (about +12 words). Until C lands, the main session
quotes the format in each `@reviewer` brief.

## 5. Live checks, after green

1. Item 1: with a grant for a test branch whose only post-review commit is
   `5eee02e8`'s kind, `python3 -c` of `push_grant.reviewed(...)` on a local
   copy of 66's branch (`worktree-coincident-vertices`, head `a5ae6356`)
   returns None.
2. Item 4: `python3 tools/session_state.py | head -1` on the main checkout
   starts `RECAP OVER CAP` while `session.md` is unpruned.
3. Item 5: the two reproductions of section 3 item 5, rerun through
   `judge_bash` and `segment_why`, pass.
4. Item 9: `timeline` over the window 2026-10-17 22:59:38 to 2026-10-18
   03:47:38 with this machine's transcripts gives an idle stretch from about
   00:32 to 03:47 (the retrospective's 3 h 15 min, section 7.2), and the
   call takes under 5 s (time it).
5. Item 5's round-1 fix and item 14: on the main checkout after PR A merges,
   `judge_bash` of ``python3 -c "print(open('<H>`cat q`').read())"`` (`H` the
   real harness state's path) is `('deny', [], '')`, and `segment_why` of
   `gh api graphql -f query='{ a }' -q {x,-X,PUT}`, given that line as its
   text, names the forge.

## 6. How h23 fits

h23 (`tools/duties.py`) says which `@orchestrator` runs are owed: a morning
check after an away window of 4 h or more, among others. Item 9 gives that
check its facts (the timeline); it does not decide whether a check is owed,
and h23 does not read transcripts. They share no code: h23 reads
`windows.jsonl`, item 9 the flag being ended. Both add lines to the recap
(h23 one line after *Last landed*; item 8 the queue), so both count against
the 10,000-character cap item 4 warns about. Both change `tools/session_state.py`
in different functions; whichever lands second merges master first
(`tools/merge_master.py`). Item 8's fallback warning is a check before a
window; h23's morning-check line is a check after it.

They also meet in two rule files (h23's design, section 6, "Prose in the
same PR", on `worktree-h23-duties` at `2843188f`). h23 names its hook and
status line in `.claude/REQUIRED-READING.md` *The harness*, the "Active in"
paragraph that cut C1 rewrites; and it adds a sentence to `CLAUDE.md` §3's
*When to spawn `@orchestrator`* bullet, in the same section where cut C3
shortens the step list and item 11 adds a bullet under *The main session
dispatches*. Whichever of h23 and PR C lands second merges master and
writes its sentence onto the other's text: if C is second, its C1 text
keeps h23's hook clause in the same short form as the others ("…;
`<h23's hook>` …"); if h23 is second, it adds its clause to C1's text.
Both are prose, so the merge is a wording step, not a design change.

## Review

- **Design review, round 1, 2026-10-18.** On `808782fd`, @reviewer: CHANGES REQUESTED. The GraphQL read exemption let writes through (a `#` inside a string, shell `$` expansion, extra words); the Python read exemption let writes through (reusing an opened file's name, `opener=`, patching `open`, rebinding `json` or `Path`); item 13's class signatures were not stable text; four wordings. Revised in `92247dec` by a fresh `@architect`.
- **Design review, round 2, 2026-10-18.** On `92247dec`, @reviewer: APPROVED. Every round-1 write is pinned by a test that it still asks; no new write passes. Items 2 and 3: Ola said yes the same day.
