# Harness h25: programs a guard cannot see because they arrive on stdin

Status: **designed**, 2026-10-19, on `worktree-h25-stdin-routes` off master
`c5da0385`. One PR. Builds on h24 PR A (items 14, 15 and 16 of
`docs/increments/h24-harness-fixes.md`, on `worktree-h24a-guards` at
`be2ab149`, not yet merged): its red and green steps start after PR A merges,
on a branch that has merged master (`tools/merge_master.py`).

Source: `@reviewer`'s findings in h24's design review round 5 (h24's
`## Review`, on `worktree-h24a-guards`): "`/dev/stdin` or `/dev/fd/0` as a
script, a pipe into an interpreter or shell, `<(…)` as a script, `.
/dev/stdin`, `python3 -m code` and `uv run python3` reading stdin; `python3 <
file` stays a known blind spot."

Acceptance: the harness suites (`pytest -m harness --no-cov`, pass count
read), the replay of section 5 rerun on the green commit, and `@reviewer`.
No `@perf` run (no refine or mesh code).

Words used here. **Scanner**: `tools/shell_scan.py`, which tells the guards
what a shell line runs and writes. **Interpreter**: a name the scanner reads a
program for (`python`, `python3`, `python3.N`, `perl`, `ruby`, `node`, and the
twelve shells in `SHELLS`). **Unknown**: the scanner's mark on a command that
writes files it does not name (`Simple.unknown`); with h24 item 15, the
governance guard then judges the whole line by its text, asking when the line
names a governed file. **Governed**: a file `guard_governance.py` asks about;
`tools/shell_scan.py` and `.claude/hooks/guard_push.py` are both, so the green
step needs Ola at the prompt. **Read** a program: give the guards its text,
as the scanner does for `python3 -c`, a heredoc and (item 15) a here-string.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature* (each read 2026-10-19).

- POSIX.1-2024, Shell Command Language, 2.9.2 *Pipelines*: "The standard
  output of command1 shall be connected to the standard input of command2."
  So a pipe hands the next command whatever the previous one prints; the
  scanner sees the producer's words, never its output.
- Python, *Command line and environment* (docs.python.org/3/using/cmdline.html):
  `-` is "Read commands from standard input (sys.stdin)"; with no interface
  option "-i is implied"; `-i` is "Enter interactive mode after execution",
  also "When a script is passed as first argument", and "Interactive mode
  will start even when sys.stdin does not appear to be a terminal". So
  `python3 -i x.py <<< S` runs `x.py`, then `S`.
- `uv help run` (uv 0.12.15): "When used with `-`, the input will be read
  from stdin, and treated as a Python script." `uv run [OPTIONS] [COMMAND]`
  runs `COMMAND`, so `uv run python -c S` is `python -c S`.
- zsh's `MULTIOS` (on by default, as h24 item 15 found for two here-strings):
  a pipe and a here-string both feed stdin. Checked 2026-10-19:
  `zsh -c "echo 'print(1)' | python3 - <<< 'print(2)'"` prints 1 and 2; bash
  prints 2.

Every route was run under zsh and bash on this machine on 2026-10-19 with a
harmless `print` or `echo` and runs the text it is given: `python3
/dev/stdin <<<`, `python3 /dev/fd/0 <<<`, `bash /dev/stdin <<<`, `printf … |
python3`, `echo … | sh`, `python3 <(echo …)`, `. /dev/stdin <<<`, `echo … |
python3 -m code`, `python3 -qi x.py <<<`, `python3 < <(echo …)`, `echo … |
tee >(sh)`, `echo … | uv run --no-project -`; and `source <(echo …)` in zsh
(macOS's bash 3.2 runs nothing for it).

*Legacy.* Nothing: the legacy tree is the old meshing product, with no harness.

```
$ git grep -l -i -e 'shell_scan' -e '/dev/stdin' -e 'guard_governance' -e 'here-string' legacy-archive -- legacy
(no output, exit 1)
```

*Earlier increments.* h4 (the scanner), h16 (wrappers and shells), h24 items
14 and 15 (a program the shell rewrites is unknown; a here-string is read).
This increment uses item 15's `finish_heredoc(simple, body, rewritten)` and
its guard line ("a line the scanner marks unknown is judged by its text").

## 2. The ideas file

`../rasputin_scratch/ideas/ideas.md`, read 2026-10-19. Overlapping, left out:
*Push guard: `watch git myalias` skips the unknown-command ask* (another
route by which a program reaches a shell, through a string runner, not
stdin); *Cap on runner and shell words per line in the push guard* (this
design adds no shell word to `runs`). Dropped already: *PowerShell flags
beyond -c*. Nothing taken up.

## 3. The design

### 3.1 One rule, three sources of stdin

A command **runs its stdin** when it is one of:

- an interpreter with no `-c`/`-e` program and a script operand that is
  absent, `-`, `/dev/stdin` or `/dev/fd/0`;
- Python running a console (`-m code`, `-m asyncio`, `-m pdb`) or with `-i`
  in a leading option cluster (`-i`, `-qi`; Python's one-letter flags only,
  so `-Wignore` is not `-i`);
- `.` or `source` whose operand is `/dev/stdin` or `/dev/fd/0`.

Where its stdin comes from decides what the scanner does:

| Stdin is | Example | The scanner |
|---|---|---|
| a heredoc or here-string on the line | `python3 /dev/stdin <<'X'`, `. /dev/stdin <<< S` | **reads** it: the program of an interpreter, parsed as commands for a shell or `.` (item 15's path) |
| another command's output | `printf S \| python3`, `echo S \| sh`, `python3 < <(echo S)`, `tee >(sh)`, `echo S \|& python3` | marks it **unknown** |
| a file | `python3 < prog.py` | nothing (section 3.4) |

And whatever stdin is, a script or `.`/`source` operand that is a process
substitution (`python3 <(…)`, `source <(…)`) or `/dev/fd/N` with N not 0
(`python3 /dev/fd/3 3< <(…)`) marks the command unknown: another command
writes the program.

Why unknown, not read, for another command's output: the producer's bytes
are not its words. zsh's `echo` turns `\n` into a newline and bash's does
not; `printf` has its formats; `cat f` prints a file. Reading them would need
each producer's rules and still miss `cat` and `curl`. The text rule sees the
same words the producer prints when they are on the line, which is the case
the guards can judge at all. Both are tripwires, as h24 item 14 says.

### 3.2 Blueprint (`tools/shell_scan.py`)

Data flow is unchanged: `parse` → `_Scan.run` → `finish` → `classify`, with
heredoc bodies through `read_heredocs` → `finish_heredoc`. New, all pure:

- Constants: `STDIN_SCRIPTS = {None, "-", "/dev/stdin", "/dev/fd/0"}`;
  `CONSOLES = {"code.py", "asyncio.py", "pdb.py"}` (as `interpreter` names a
  `-m` module); `SOURCES = {".", "source"}`; `FED_SCRIPT`, a pattern for
  `<(` or `/dev/fd/N` with N not 0.
- `stdin_program(simple) -> bool`: the rule of 3.1. For an interpreter it
  first checks `interpreter(name, argv[1:])` for an inline program (as item
  15's inline check does), so a shell's `-c`, already parsed and cleared,
  keeps stdin as data.
- `interactive(args) -> bool`: `-i` in a leading Python option cluster,
  stepping over `-W` and `-X` values.
- `fed(simple, piped) -> bool`: the script or `.` operand matches
  `FED_SCRIPT`, or `piped and stdin_program(simple)`.
- `_Scan.run(close, piped=False)`: a local `piped`, true for the simple after
  `|` or `|&` (not `||`), after a `<` redirect from `<(…)`, and for the first
  simple inside `>(…)` (`word` calls `run(")", True)` for `>(` only). Reset
  after `;`, `&&`, `||`, newline. A subshell `( … )` after a pipe gets it
  (`echo S | (python3)`), and `{` is already dropped as a reserved word.
- `_Scan.redirect` returns whether it fed stdin from a command; it accepts a
  process substitution as a target (today `x < <(…)` is unreadable, judged as
  text) and never records `>(…)` as a written file.
- `finish(…, piped)`: `simple.unknown |= fed(simple, piped)` after `classify`.
- `finish_heredoc`: its gate becomes `stdin_program(simple)`, and `.` and
  `source` parse the body as a shell does. A piped command with a heredoc
  gets both (zsh feeds both).
- `uv run` is stripped like a wrapper, in `unwrap`: `uv_run(words after
  "run")` skips options and their values up to the first interpreter word
  (not the value of `-p`/`--python`), and maps uv's own `-` and `-m`/
  `--module` to `python3 -` and `python3 -m`. A plain word that follows no
  option ends the search: `uv run pytest`, `uv run ruff check python` keep
  `uv` as argv[0], as today. This also reads `uv run python -c S`, silent
  today (section 5).

### 3.3 The guards

- `guard_governance.py`: no change. Item 15's guard line already judges by
  text any line with an unknown command, and the `.git/harness` backstop
  denies a line that names the harness state and has an unknown command.
- `guard_push.py` (`main`, about 3 lines): the push guard reads a shell's
  program only through the scanner's nested simples, and never uses
  `unknown`. When any simple is a shell or `.`/`source` marked unknown, its
  reasons become the text rules' (`PUBLISHES` over the whole line), as for a
  line the scanner cannot read: `echo 'git push origin x' | sh` asks. An
  interpreter (`… | python3`) adds nothing to the push guard: it never reads
  a Python program for `git push` (`python3 -c "os.system('git push')"`
  passes today), and this design keeps that.

### 3.4 `python3 < file`, and what stays out

`python3 < prog.py` is `python3 prog.py` by another spelling: the program is
a file, which neither guard reads for any script. It is bounded the same
way: the file must exist on disk before the line runs, and was written either
earlier (a write the guards judged then: an Edit, a Write or a shell write,
judged by its target, not its content) or in the same line (`echo S >
p.py; python3 < p.py`, where the guards judge the write of `p.py` and not
what it holds, as for `python3 p.py`). A guard that opened the file would read
it before the line runs, not when Python does. So no change; the gap is the
script gap, named in the docstring.

Also not covered, each as today: `cat f | sh` and `curl … | sh` (unknown,
but the text rule finds no name in a file or a download); an interpreter or
module that `exec`s stdin of its own accord (`python3 x.py`,
`python3 -m mymod`): only the three known consoles count; `uv run x.py`,
`uv run https://…`, `uvx`, `poetry run`, `pipenv run` (none in the replay's
`run` commands but `uv`; 0 uses of the last two); `bash <&3` and other
descriptor copies; `eval "$(…)"` and `bash -c "$(cat f)"` (not stdin; see
*Ideas* in the handback).

### 3.5 Cost in false asks

A new ask needs both: a command that runs its stdin from another command (or a
fed script), and a governed name in the line's text, or a `git push`-shaped
text beside a fed shell. The replay (section 5) found 0 verdict changes over
67,534 past commands; agents write `python3 - <<'EOF'`, not pipes into
interpreters. The likely false ask: `cat CLAUDE.md | python3 -` reading a
rule file on stdin as a program. `cat CLAUDE.md | python3 -c S` (stdin as
data) stays silent.

## 4. Tests

`H` is the harness state's absolute path; `W` is
`open('CLAUDE.md','w').write('x')`. Every verdict below was measured on
2026-10-19 on two scratch copies of `be2ab149` with h24 items 14 and 15
applied as designed (`old`), and the same plus this design (`new`); "old"
results are all silent unless stated.

Through the whole governance hook (`tests/python/test_guard_governance.py`,
`repo` fixture), mode off unless named:

- ask naming `CLAUDE.md` (read, so the reason is "it writes CLAUDE.md"):
  `bash /dev/stdin <<< 'echo x > CLAUDE.md'`; `python3 /dev/stdin <<< "W"`;
  `python3 /dev/fd/0 a b <<'X'` with body `W`; `. /dev/stdin <<< 'echo x >
  CLAUDE.md'`; `python3 -m code <<< "W"`; `python3 -qi x.py <<< "W"`;
  `uv run python3 - <<< "W"`; `uv run python -c "W"`; `uv run --with vtk
  python -c "W"`; `uv run -p 3.12 python3 -c "W"`; `uv run -m code <<< "W"`;
- ask naming `CLAUDE.md` (unknown, the text rule): `printf "W" | python3`;
  `echo 'echo x > CLAUDE.md' | sh`; `cat <<'X' | python3 -` with body `W`;
  `python3 <(echo "W")`; `source <(echo 'echo x > CLAUDE.md')`; `echo "W" |
  python3 -m code`; `echo "W" | uv run -`; `python3 < <(echo "W")` (on old
  unreadable and silent); `echo "W" | { python3; }`; `echo "W" | (python3)`;
  `echo "W" |& python3`; `echo "W" | env python3`; `echo "W" | perl`; `echo
  "W" | python3 -m pdb x.py`; `python3 /dev/fd/3 3< <(echo "W")`;
- denied with the always-denied text, nothing queued, modes off and on:
  `echo "open('H','w').write('x')" | python3` and `bash /dev/stdin <<< 'echo
  x > H'`;
- silent, on both copies: `python3 < prog.py`; `echo "W" || python3`; `echo
  "W" && python3`; `echo "W" | python3 x.py`; `git diff | python3
  tools/count_loc.py`; `echo '{}' | python3 -m json.tool`; `cat CLAUDE.md |
  python3 -c 'import sys; print(len(sys.stdin.read()))'`; `source
  .venv/bin/activate`; `uv run pytest tests/python`; `uv run ruff check
  python`; `printf 'print(1)' | python3`; `python3 -Wignore x.py <<< "W"`;
  `while read l; do echo $l; done < <(ls CLAUDE.md)`.

Through the whole push hook (`tests/python/test_guard_push.py`, `ASKED` and
`PASSED`): ask: `bash /dev/stdin <<< 'git push origin x'`, `source
/dev/stdin <<< 'git push origin x'` (read), `echo 'git push origin x' | sh`,
`bash <(echo 'git push origin x')`, `echo 'git push origin x' | tee >(sh)`
(unknown, text rules); unchanged asks: `echo 'git push origin x' > >(sh)`,
`bash < <(echo 'git push origin x')`, `bash -s < <(echo 'git push origin x')`
(old: unreadable, text rules; new: read as fed). Passes: `echo 'git push
origin x' | cat`, `echo 'git push origin x' | python3` (interpreter, 3.3).

In `tests/python/test_guard_targets.py`, `parse`: `printf S | python3` gives
`python3` with `unknown` True and `echo S | python3 x.py` False; `python3
/dev/stdin <<< "W"` gives program `W\n` and candidate `CLAUDE.md`; `uv run
--with vtk python -c S` gives argv[0] `python`; `uv run pytest` keeps argv[0]
`uv`; `echo x > >(sh)` has no write.

Existing suites stay green on the green commit (section 5 measured the
prototype: no test changes result).

## 5. Measured

*Replay.* Every distinct Bash command in this machine's Claude Code
transcripts (`~/.claude/projects/**/*.jsonl`, read with `json`, commands
never printed except those matching the routes, a line with `key`, `token`,
`secret`, `passw`, `authoriz` or `bearer` redacted): 67,534. Through both
copies' `judge_bash` and the push hook's `main`: **0 verdicts changed**
(governance asks 432, push asks 748, on both). The probe can fail: the route
lines of section 4 change verdict on the same harness. `parse` output
changed for 168 lines: 148 newly read programs (mostly `uv run python -c`,
and `python3 /dev/stdin args <<'X'`, which agents use to pass arguments with
a heredoc program), 16 newly unknown (all `uv run python -c "…$S…"`, item 14's
rule reached through the new wrapper), 4 lines newly readable (`… < <(…)`).
Each took about 13 s for the whole corpus per copy.

*Suites.* `pytest -m harness --no-cov` (the main checkout's venv, run from
each copy): `old` 2626 passed, 9 failed, 14 skipped; `new` the same, with the
same 9 failures. They are PR A red tests the prototype does not meet (5a's
review-round fixes, and item 14's judging by static form): PR A's green
step, not this design. The 14 skips: a copy is no git work tree.

## 6. Size and order

`tools/shell_scan.py` about 40 counted lines (constants 4, `stdin_program` 9,
`interactive` 6, `fed` 5, `uv_run` 9, `run`/`redirect`/`finish`/`unwrap`
about 8, `finish_heredoc` net -2), `.claude/hooks/guard_push.py` about 3;
total about 45, measured on the prototype diff (+55 -17 counted in the
scanner, +3 in the push guard). Tests about 120 lines.

Path: the full path despite the size (under the short path's 50): it is a
guard, where review probes routes adversarially, and a design round costs
less than a governed green step redone. Order: design review now; red after
h24 PR A merges and this branch merges master; green by day with Ola at the
prompt (both files governed). Not a night fallback.
