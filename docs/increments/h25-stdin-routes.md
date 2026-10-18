# Harness h25: programs a guard cannot see because they arrive on stdin

Status: **designed**, revised 2026-10-19 for design review round 1, on
`worktree-h25-stdin-routes` off master `c5da0385`. One PR. Builds on h24 PR A
(items 14, 15 and 16 of `docs/increments/h24-harness-fixes.md`, on
`worktree-h24a-guards` at `be2ab149`, not yet merged): its red and green
steps start after PR A merges, on a branch that has merged master
(`tools/merge_master.py`).

Source: `@reviewer`'s findings in h24's design review round 5 (h24's
`## Review`, on `worktree-h24a-guards`): "`/dev/stdin` or `/dev/fd/0` as a
script, a pipe into an interpreter or shell, `<(…)` as a script, `.
/dev/stdin`, `python3 -m code` and `uv run python3` reading stdin; `python3 <
file` stays a known blind spot." Then this design's review round 1
(CHANGES REQUESTED): the first design listed the routes one by one, and the
review found more than a dozen it missed (section 3.1 lists them); and the
push guard reads no interpreter program at all.

Acceptance: the harness suites (`pytest -m harness --no-cov`, pass count
read), the replay of section 5 rerun on the green commit, and `@reviewer`.
No `@perf` run (no refine or mesh code).

Words used here. **Scanner**: `tools/shell_scan.py`, which tells the guards
what a shell line runs and writes. **Simple**: one simple command as the
scanner returns it. **Interpreter**: a name the scanner reads a program for
(`python`, `python3`, `python3.N`, `perl`, `ruby`, `node`), and **shell**:
one of the twelve names in `SHELLS`. **Unknown**: the scanner's mark on a
command that writes files it does not name (`Simple.unknown`); with h24 item
15, the governance guard then judges the whole line by its text, asking when
the line names a governed file. **Governed**: a file `guard_governance.py`
asks about; `tools/shell_scan.py` and `.claude/hooks/guard_push.py` are both,
so the green step needs Ola at the prompt. **Read** a program: give the guards
its text, as the scanner does for `python3 -c`, a heredoc and (item 15) a
here-string. **Text rule**: judging a line by the names (governance) or the
`PUBLISHES` patterns (push) in its whole text.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature* (each read 2026-10-19).

- POSIX.1-2024, Shell Command Language, 2.9.2 *Pipelines*: "The standard
  output of command1 shall be connected to the standard input of command2."
  Its grammar (2.10.2) has `pipe_sequence '|' linebreak command`: a newline
  may follow `|` and the pipeline goes on. A pipeline element may be a
  compound command, whose every inner command reads the pipe.
- bash(1), macOS bash 3.2: "-s If the -s option is present, or if no
  arguments remain after option processing, then commands are read from the
  standard input." So `bash -s foo` runs stdin, and a shell's options and
  operands cannot be read safely to decide whether it reads stdin.
- zshexpn(1), *Process Substitution*: "If =(...) is used instead of <(...),
  then the file passed as an argument will be the name of a temporary file
  containing the output of the list process."
- Python, *Command line and environment* (docs.python.org/3/using/cmdline.html):
  `-` is "Read commands from standard input (sys.stdin)"; `-i` is "Enter
  interactive mode after execution", "Interactive mode will start even when
  sys.stdin does not appear to be a terminal". So `python3 -i -c pass <<< S`
  runs `S`.
- `uv help run` (uv 0.12.15): "When used with `-`, the input will be read
  from stdin, and treated as a Python script."
- zsh's `MULTIOS` (on by default): a pipe and a here-string both feed stdin.
  `zsh -c "echo 'print(1)' | python3 - <<< 'print(2)'"` prints 1 and 2; bash
  prints 2.

Run under zsh and bash on this machine on 2026-10-19 with a harmless
`print`, each running the text it is fed: a newline after `|`;
`| { true; python3; }`; `| bash -c python3`; `| bash -s foo` (bash parses
the Python text and reports a syntax error); `python3 -i -c pass <<< S`;
`exec 0< <(…); python3`; `f() { python3; }; … | f`;
`python3 3< <(…) <&3`; `| xargs -I{} python3 -c {}`; in zsh,
`python3 =(…)` and `coproc python3 -i` fed by `print -p`. macOS's bash 3.2
has no `coproc`. The first design's routes were run the same way (its
section 1, commit `d5e8c192`).

*Legacy.* Nothing: the legacy tree is the old meshing product, with no harness.

```
$ git grep -l -i -e 'shell_scan' -e '/dev/stdin' -e 'guard_governance' -e 'here-string' legacy-archive -- legacy
(no output, exit 1)
```

*Earlier increments.* h4 (the scanner), h16 (wrappers and shells), h24 items
14 and 15 (a program the shell rewrites is unknown; a here-string is read;
"a line the scanner marks unknown is judged by its text").

## 2. The ideas file

`../rasputin_scratch/ideas/ideas.md`, read 2026-10-19. Overlapping, left out:
*Push guard: `watch git myalias` skips the unknown-command ask* (a string
runner, not stdin; this design's runner rule does see `watch python3`, but
does not change how the push guard reads `watch`'s words); *Cap on runner and
shell words per line in the push guard* (this design adds no shell word to
`runs`). Nothing taken up.

## 3. The design

### 3.1 Why a broader rule, and why not the broadest

The first design listed the ways a program reaches stdin (`-`,
`/dev/stdin`, a pipe, `-m code`, `-i`, …) and reset its "piped" mark at every
`;` and newline. Review round 1 ran each of these under zsh and bash; each
runs its fed text and escaped that design:

1. a newline after a pipe (`echo S |⏎python3`);
2. a compound command after a pipe (`| { true; python3; }`, `| (true;
   python3)`, `| if true; then python3; fi`);
3. a shell's `-c` program inheriting fed stdin (`echo S | bash -c python3`,
   `bash -c python3 <<< S`);
4. words read as the script operand (`bash -s foo`, `bash -o errexit`,
   `ruby -r json`, `sudo -E python3 -`);
5. `python3 -i -c pass`;
6. zsh's `=(…)` as a script;
7. `python3 3< <(…) < /dev/fd/3` and `python3 3< <(…) <&3`;
8. `xargs -I{} python3 -c {}`, `env -S 'python3 -'`, `| $SHELL`,
   `//dev/stdin`, `exec 0< <(…); python3`;

and, lower priority, `python3 -c 'exec(input())'`, `python3 -m cProfile
/dev/stdin`, `arch -arm64 python3`. Its claim that `python3 < /dev/fd/N`
was "the same gap as `python3 file`" was false: the descriptor is fed from
the line.

So this design stops reading options to decide whether a command runs its
stdin. It asks instead whether the command **could** run text (a *runner*)
and whether its stdin or script **could** carry text the scanner did not
read (a *feed*). Both are judged broadly; the one narrowing is a Python data
reader (3.4), whose CLI is small and documented.

The broadest form the review suggested, "any interpreter or shell anywhere
on a line that also has a pipe, a redirect or a compound", was measured
first: on the replay of section 5 it adds **5,785** governance asks out of
67,715 past commands (8.5 %), because `python3 tools/check_citations.py |
tail -5` names a governed tool and holds a pipe. A guard that asks on one
command in twelve gets approved without being read. So the feed is judged on
the runner's own stdin, which keeps the pipe *after* a runner from counting,
and every other rule stays broad. With that, the added asks are 33
governance and 22 push (section 5).

### 3.2 The rule

A simple is **unread** (new field `Simple.unread: bool`) when it is a runner
and it is fed. An unread simple is also marked `unknown`, so the governance
guard judges the line by its text (item 15's guard line), and the push guard
applies its text rules to the whole line (3.5).

**Runners.** A simple is a runner when any of these holds:

- after wrappers, its command word is an interpreter or a shell, or `.`,
  `source`, `eval`, `exec`, `coproc`, `parallel` (with no command, GNU
  `parallel` runs its input lines);
- its command word as written holds `$` or a backtick (`$SHELL`,
  `$(which python3)`, `"$PY"`): it could be anything;
- any whitespace-separated token of any of its words is an interpreter or
  shell name, `uv` or `uvx`, unless the command word is one of `NO_RUN`
  (below). This catches `sudo -E python3 -`, `env -S 'python3 -'`,
  `arch -arm64 python3`, `uv run python3`, `xargs sh` and options whose
  values the wrapper stripping misreads, without listing any of them;
- it is `xargs` in replacement mode (a word holding `{}`, or starting `-I`,
  `-J`, `-i` or `--replace`) or with an expansion in a word: its input can
  become the command;
- its command word is a function defined on the line (`f() …`,
  `function f …`) or an alias defined on it (`alias f=…`).

`NO_RUN` is the commands that run no program whatever their arguments, so an
interpreter's name among their words is data (`grep -i python`): `grep`,
`egrep`, `fgrep`, `rg`, `head`, `tail`, `cat`, `wc`, `sort`, `uniq`, `cut`,
`tr`, `ls`, `echo`, `printf`, `diff`, `cmp`, `comm`, `basename`, `dirname`,
`realpath`, `test`, `[`, `which`, `type`, `file`, `stat`, `du`, `df`, `tee`,
`mkdir`, `rm`, `cp`, `mv`, `touch`. A list of what is safe, not of what is
dangerous: a command missing from it is a runner when an interpreter's name
is among its words, never the reverse.

**Feeds.** A runner is fed when any of these holds:

- **Piped.** It stands in a pipeline after `|` or `|&`. The scanner keeps a
  piped state per level: set after `|`; kept across a newline that directly
  follows the pipe (blanks, a comment or a backslash-newline between);
  inherited by everything inside a compound command (`{ }`, `if … fi`,
  `while`/`until`/`for`/`select … done`, `case … esac`), a subshell `( )`,
  `$(…)`, backticks, and the `-c` program of a piped shell; set inside
  `>(…)`. A `;`, `&`, `&&`, `||` or a newline resets it to the enclosing
  compound's state, not to "not piped". So `echo S | if true; then python3;
  fi` is fed, and `python3 x | tail; python3 y` is not.
- **On its own command**: an input redirect `<`, `<&` or `<>`, from any
  descriptor, whose target is not `/dev/null`; a process substitution as a
  word or a redirect target (`<(…)`, `>(…)`, zsh's `=(…)`).
- **On the line**, feeding every runner on it: an input redirect, a
  process substitution, a heredoc or a here-string on a compound command,
  a subshell or `exec` (`{ python3; } < <(…)`, `done <<< S`,
  `exec 0< <(…)`); `coproc` (bash, zsh); zsh's `cmd |&` with nothing after
  it (a coprocess); a heredoc or here-string on a runner whose command word
  is not literally an interpreter or shell (`$SHELL <<< S`, `sudo -E
  python3 - <<< S`); a heredoc or here-string on a shell that has any
  argument (`bash -c python3 <<< S`, `bash -s foo <<< S`: the inner commands
  or the shell's own reading take it), or whose body does not parse as
  commands.

**Reading a heredoc or here-string.** Item 15's gate ("an interpreter with
no inline program and a script that is absent or `-`") is replaced: a body
on an interpreter whose command word is literally its name is appended to
its program **whatever its options and operands**, unless the command is a
data reader (3.4). This reads `python3 -i -c pass <<< S`, `python3
//dev/stdin <<< S`, `python3 -m cProfile /dev/stdin <<< S` and `python3 -c
'exec(input())' <<< S` without knowing any of those options. A body on a
literal shell, `.` or `source` is parsed as commands, as today, and every
command in it counts as piped (it reads the rest of the body).

**What it covers.** Every escape in 3.1 asks (section 4), and so do the
first design's routes, `| uv run python3`, `uv run python3 - <<< S`,
`f() { python3; }; … | f`, `coproc python3; echo S >&…`, `| parallel`, and
`… | while read l; do eval "$l"; done`.

### 3.3 Blueprint (`tools/shell_scan.py`)

Data flow is unchanged: `parse` → `_Scan.run` → `finish` → `classify`, with
heredoc bodies through `read_heredocs` → `finish_heredoc`. All new code is
pure and stdlib.

- Constants: `RUNNER_FIRST` (`.`, `source`, `eval`, `exec`, `coproc`,
  `parallel`), `RUNNER_ANY` (the shells, `perl`, `ruby`, `node`, `uv`,
  `uvx`; Python by `PYTHON`), `NO_RUN`, `OPENERS` (`{ if while until for
  case select`), `CLOSERS` (`} fi done esac`), and for 3.4 `PLAIN_FLAGS`,
  `DATA_MODULES`, `RUNS_CODE`.
- `Simple` gains `unread: bool = False`; the scanner keeps per simple, as
  private state, its words as written (before wrappers are stripped), its
  piped state and whether its own command is fed.
- `runner(simple) -> bool`, `data_reader(simple) -> bool`,
  `data_program(text) -> bool`: the rules of 3.2 and 3.4.
- `parse(command, piped=False)`: after the scan, marks every simple that is
  a runner, fed (piped, on its own command, or the line), and not a data
  reader: `unread = unknown = True`. A nested parse (a shell's `-c` program,
  `eval`'s words, a heredoc body for a shell) passes the outer command's
  piped state (`True` for a heredoc body).
- `_Scan.run(close, inherit=False)`: the piped state and the compound stack
  of 3.2, keyed on `OPENERS`/`CLOSERS` met in command position; subshells,
  `$(…)`, backticks and process substitutions call it with the current
  state (`>(…)` with `True`).
- `_Scan.redirect`: accepts a process substitution as a target (today
  `x < <(…)` is unreadable and judged as text, which needs `WRITES` to
  match, so `python3 < <(…)` passes silently), never records `<(…)` or
  `>(…)` as a written file, and notes an input redirect.
- `finish`: notes a fed compound or `exec` for the whole line; records the
  names of functions and aliases defined.
- `finish_heredoc(simple, body, …)`: the reading rule of 3.2 in place of
  item 15's gate; sets the line feeds of 3.2.

### 3.4 The one narrowing: a Python data reader

Without it, 121 of the governance asks of the replay are `echo JSON |
python3 .claude/hooks/guard_x.py` or `… | python3 -c 'import json,sys; …'`
on a line that names a governed file: stdin is data there. A simple is a
**data reader** when all hold:

- its command word is literally `python`, `python3` or `python3.N`, with no
  wrapper before it;
- its options are only clusters of the no-value flags `-B -E -I -O -P -q -s
  -S -u -b -d -R -v -x`, or `-W`/`-X` with their value: never `-i`, `-m`,
  `-` or a long option;
- then either a script operand that is a plain path (no `$`, backtick, `<`,
  `>`, `=`, parenthesis; not under `/dev/`; not starting `//`), or `-c PROG`
  where `PROG` parses, imports only `DATA_MODULES` (`json sys re collections
  csv math itertools statistics textwrap datetime string functools operator
  unicodedata`), and uses no name or attribute in `RUNS_CODE` (the scanner's
  `REFLECTION`, plus `input open pickle marshal os subprocess runpy code pdb
  importlib ctypes shutil pathlib Path`) and no attribute starting `_`.

Why this holds: Python runs text as code only through `exec`, `eval`,
`compile`, an import of code (`importlib`, `runpy`, `pickle`, `code`,
`pdb`) or a child process (`os`, `subprocess`, `ctypes`); a `-c` program
that names none of them cannot run its stdin. A named script that runs its
stdin is the script gap of 3.6, as for `python3 x.py` today. A data reader's
heredoc stays data, as item 15 has it. `python3 -c 'exec(input())'`, `-i`,
`-m` and `/dev/stdin` are never data readers.

### 3.5 The push guard (`.claude/hooks/guard_push.py`, `main`)

Two additions, after `parse` and before `runs`:

- **An unread simple**: the reasons gain the text rules' (`PUBLISHES` over
  the whole line), as for a line the scanner cannot read. `echo 'git push
  origin x' | sh` and `… | { true; sh; }` ask. `unread`, not `unknown`:
  `git apply x.patch; echo 'git push origin x'` stays silent.
- **Every program the scanner reads.** Confirmed in review round 1 on h24
  PR A's code too: the push guard reads a shell's program through the
  nested simples but never an interpreter's, so `python3 -c
  "os.system('git push origin x')"` passes. Now the `PUBLISHES` patterns run
  over each simple's `program` (`-c`/`-e` text, a heredoc, a here-string)
  that **starts a process**: for Python, its syntax tree names one of
  `SPAWN_NAMES` (`subprocess system popen Popen run call check_call
  check_output getoutput getstatusoutput pty spawn create_subprocess_exec
  create_subprocess_shell __import__ getattr importlib ctypes runpy eval
  exec compile fork forkpty posix_spawn posix_spawnp startfile`, or a name
  starting `exec` or `spawn`); for Perl, Ruby, Node or Python that does not
  parse, the pattern `\b(subprocess|system|popen|exec\w*|spawn\w*|fork|pty|child_process|open3|qx)\b`
  or a backtick.

*A departure from the review's wording*, which asked for the text rules
over every program: that form adds 115 push asks on the replay; with "starts
a process" it adds 22, every one a harness probe. A program that names `git
push` without starting a process (an edit of a design file, a `print`)
publishes nothing.

### 3.6 What stays uncovered

Each is bounded the same way: the text a runner gets is not on the line, so
the text rule has nothing to find.

- **Text from a file or the network**: `cat f | sh`, `curl … | sh`, `python3
  < prog.py` (unread, but the file's names are not in the line), `eval
  "$(…)"`, `bash -c "$(cat f)"`, `source f`.
- **A named script that runs its stdin**: `echo S | python3 x.py` where
  `x.py` execs it (a data reader), as `python3 x.py` is read for no script.
- **Commands that run their input and are no runner here**: `awk
  '{system($0)}'`, GNU `sed`'s `e` command, `find -exec`, and any program
  not named an interpreter or shell (`osascript`, `lua`, `php`, `Rscript`).
  Adding a name to `RUNNER_ANY` is one word; none appears in the replay.
- **A program that builds the name** (`'CLA' + 'UDE.md'`): the text rule's
  bound, as in h24 item 14.
- **`uv run python -c S`** and other inline programs behind `uv run`: not a
  stdin route; the program is not read today and is not read here (the first
  design stripped `uv run`; dropped as unneeded for stdin, handed back as an
  idea).
- **The Bash tool's own stdin**: what a runner inherits with no feed on the
  line is not the agent's to fill.

### 3.7 Cost in false asks

A new governance ask needs an unread runner and a governed name in the
line's text, or a read heredoc naming a governed file. A new push ask needs
an unread runner and a `PUBLISHES` match in the line, or a process-starting
program that matches. Measured in section 5: 33 and 22 over 67,715 commands.
Inspected by hand: the push asks are all harness probes (Python programs
that build `git push` events and run the hooks). Of the governance asks,
about half are guard probes piping events into `.claude/hooks/*.py` through
a shell function, `$PY` or `sh -c`; the rest pipe a governed tool's output
(`tools/check_citations.py`, `tools/session_state.py`) or a rule file into a
Python program the data-reader rule does not exempt, or a heredoc that
names a governed file into a Python script whose path is a variable.

## 4. Tests

`H` is the harness state's absolute path; `W` is
`open('CLAUDE.md','w').write('x')`. Every verdict below was measured on
2026-10-19 on a prototype (section 5): `old` is `be2ab149`'s scanner and
guards with item 15 as designed (here-strings read; the guard line), `new`
the same plus this design. "old" results are all silent unless stated.

Through the whole governance hook (`tests/python/test_guard_governance.py`,
`repo` fixture), mode off unless named:

- ask naming `CLAUDE.md`, the review's escapes: `echo "W" |⏎python3`;
  `echo "W" | { true; python3; }`; `echo "W" | (true; python3)`; `echo "W" |
  if true; then python3; fi`; `echo "W" | bash -c python3`; `bash -c python3
  <<< "W"`; `echo 'echo x > CLAUDE.md' | bash -s foo`; `bash -s foo <<< 'echo
  x > CLAUDE.md'`; `echo 'echo x > CLAUDE.md' | bash -o errexit`; `echo "W" |
  sudo -E python3 -`; `python3 -i -c pass <<< "W"`; `echo "W" | python3 -i
  -c pass`; `python3 =(echo "W")`; `python3 3< <(echo "W") < /dev/fd/3`;
  `python3 3< <(echo "W") <&3`; `echo "W" | xargs -I{} python3 -c {}`;
  `echo "W" | env -S 'python3 -'`; `echo 'echo x > CLAUDE.md' | $SHELL`;
  `python3 //dev/stdin <<< "W"`; `exec 0< <(echo "W"); python3`; `echo "W" |
  python3 -c 'exec(input())'`; `python3 -c 'exec(input())' <<< "W"`; `echo
  "W" | python3 -m cProfile /dev/stdin`; `echo "W" | arch -arm64 python3`;
- ask naming `CLAUDE.md`, further routes: `echo "W" |  # c⏎python3`;
  `echo "W" | (python3)`; `echo "W" |& python3`; `echo "W" | python3 -m
  code`; `echo "W" | uv run -`; `echo "W" | uv run python3`; `uv run python3
  - <<< "W"`; `echo "W" | tee >(python3)`; `python3 <(echo "W")`; `source
  <(echo 'echo x > CLAUDE.md')`; `echo 'echo x > CLAUDE.md' | while read l;
  do eval "$l"; done`; `echo "W" | echo "$(python3)"`; `echo "W" | { cat; }
  | python3`; `f() { python3; }; echo "W" | f`; `function g { sh; }; echo
  'echo x > CLAUDE.md' | g`; `coproc python3; echo "W" >&"${COPROC[1]}"`;
  `$SHELL <<< 'echo x > CLAUDE.md'`; `bash -c 'echo "W" | python3'`; `echo
  "W" | sh -c 'cat | python3'`; `{ python3; } <<< "W"`; `{ python3; } <
  <(echo "W")`; `( python3 ) < <(echo "W")`; `while read l; do python3;
  done < <(echo "W")`; `echo 'echo x > CLAUDE.md' | parallel`; `cat
  CLAUDE.md | python3 -c "import os,sys; os.system(sys.stdin.read())"`;
  unchanged (asks on both): `echo "W" > >(python3)`;
- denied with the always-denied text, nothing queued, modes off and on:
  `echo "open('H','w').write('x')" | python3`; `bash -s foo <<< 'echo x >
  H'`; `{ python3; } < <(echo "open('H','w').write('x')")`;
- silent, on both: `python3 < prog.py`; `echo x || python3`; `echo x &&
  python3`; `echo "W" | python3 x.py` (a data reader); `git diff | python3
  tools/count_loc.py`; `cat CLAUDE.md | python3 -c 'import sys;
  print(len(sys.stdin.read()))'`; `python3 -c 'import json,sys;
  print(json.load(sys.stdin))' < CLAUDE.md`; `python3 tools/x.py <<'EOF'`
  with body `{"path": "CLAUDE.md"}`; `python3 x.py < /dev/null`; `diff
  <(echo "W") x; python3 x.py`; `echo '{}' | python3 -m json.tool`;
  `printf 'print(1)' | python3`; `source .venv/bin/activate`; `uv run pytest
  tests/python`; `while read l; do echo $l; done < <(ls CLAUDE.md)`.

Through the whole push hook (`tests/python/test_guard_push.py`, `ASKED` and
`PASSED`). Ask: `python3 -c "import os; os.system('git push origin x')"`;
`python3 -c "import subprocess; subprocess.run(['git','push','origin','x'])"`;
`perl -e 'system("git push origin x")'`; `bash -c python3 <<< "import os;
os.system('git push origin x')"`; `echo 'git push origin x' | sh`; `echo
'git push origin x' | bash -s`; `echo 'git push origin x' | { true; sh; }`;
`sh <(echo 'git push origin x')`. Pass: `echo 'git push origin x' | cat`;
`echo 'git push origin x' | python3 x.py`; `python3 -c "print('git push
origin x')"`; `git apply x.patch; echo 'git push origin x'`.

In `tests/python/test_guard_targets.py`, `parse`: `echo S |⏎python3` gives
`python3` with `unread` and `unknown` True; `python3 x | tail -1; python3 y`
gives both `python3` simples `unread` False; `echo S | if true; then
python3; fi` gives `python3` `unread` True; `grep -i python CLAUDE.md |
head` gives no `unread`; `python3 -i -c pass <<< "W"` gives program
`pass\nW\n` and candidate `CLAUDE.md`; `echo x > >(sh)` has no write.

Existing suites: the prototype was not run against them (section 5); the
green step reports any existing test whose verdict changes, with the line.

## 5. Measured

*Prototype.* A scratch copy of `be2ab149`'s `tools/` and `.claude/hooks/`,
with item 15 as designed (a here-string read, and the guard line) and this
design behind a switch, so `old` and `new` are the same files. Item 14 is
not in it: its asks are on lines with `$` in a program, which the replay's
`old` count therefore misses, so the added asks below are an upper bound.

*Replay.* Every distinct Bash command in this machine's Claude Code
transcripts (`~/.claude/projects/**/*.jsonl`, read with `json`): 67,715.
Through `judge_bash` and the push hook's `main` (its `emit` captured, no
grant module), `old` against `new`. Commands holding `key`, `token`,
`secret`, `passw`, `authoriz` or `bearer` were counted but never printed.

| Rule | governance asks added | push asks added |
|---|---|---|
| any runner on a line with any pipe, redirect or process substitution | 5,785 | 241 |
| fed on the runner's own stdin, no data reader | 121 | 60 |
| **this design** (3.2 to 3.5) | **33** (4 redacted) | **22** (2 redacted) |
| this design, push text rules over every program | 33 | 115 |

Governance asks: `old` 438, `new` 471; push asks: `old` 748, `new` 770.
No verdict went from ask to silent. The probe can fail: every route line of
section 4 changes verdict on the same harness, and each silent line was run
too. Each pass over the corpus took about 27 s.

*Suites.* Not run on the prototype: the copy of the test tree into the
scratchpad was refused by the governance guard in unattended mode (it
writes copies of governed files through a shell variable), so the suites'
figures come from the red and green steps.

## 6. Size and order

`tools/shell_scan.py` about 90 counted lines: constants about 14, `runner`
8, `data_reader` and `data_program` 25, the marking in `parse` 7, the piped
state and compound stack in `run` 15, `word`/`dollar`/`redirect` 8,
`finish` 6, `finish_heredoc` net 7. `.claude/hooks/guard_push.py` about 15
(`spawns` and its two constants 11, `main` 4). Total about 105, from the
prototype diff less its switches. Tests about 160 lines.

Path: full. Order: design review now; red after h24 PR A merges and this
branch merges master; green by day with Ola at the prompt (both files
governed). Not a night fallback.
