# Harness h37: governance guard follow-ups of h33

Status: **designed**, 2026-10-27, on `worktree-guard-zsh` off master
`aabac548`. Short path (`docs/increments/README.md`): about 40 counted lines,
no C++ file. Both code files are governed (`.claude/hooks/guard_governance.py`,
`tools/shell_scan.py`), so the green step runs with Ola at the guard's prompt,
by day. Source: `h33-harness-fixes-1027.md` §3, *Out of scope, follow-up*,
and its code review.

## Prior art: legacy and literature

*Literature.* None applies: repository tooling, no published method, no
novelty claimed. The shell forms are zsh's own (`man zshmisc`, *Redirection*).
*Legacy.* `git grep -n -i -E 'noclobber|guard_governance|shell_scan'
legacy-archive -- legacy` returns nothing. *Ideas file.* Overlaps "governance
ask reach Ola in auto mode" (not taken up: these fixes only matter if the ask
reaches him).

## 1. Probe of today's guard

`judge_bash` and the Write arm at `aabac548`, loaded from a scratch probe;
"refused" is a deny or an ask. zsh 5.9 with `noclobber` set confirmed, in a
scratch folder, that each redirect form below writes (the `2>` forms
truncate), and that a glob reaches an existing file in a redirect and in
`touch`.

| Case | Today |
|---|---|
| (1) `echo x OP CLAUDE.md`, OP one of `>!` `>>!` `>>|` `&>!` `>&!` `>&|` `>>&` `&>>!` `&>>|` `>>&!` `>>&|` `2>!` `2>>|`; also `>>|CLAUDE.md`, `>! CLAUDE.md` alone, `exec {fd}>! CLAUDE.md`, `>!` into `.claude/settings.json` and `tools/away.py` | **passes** (each) |
| (1) `>`, `>|`, `&>|`, `&>>`, `>!CLAUDE.md` (no space: target `!CLAUDE.md`, tail match), `exec {fd}> CLAUDE.md` | refused |
| (2) `echo x >` or `touch`/`rm` on `.git//harness/m`, `.git/./harness/m`, `./.git//harness/m`, `.git///harness/m`, `/r/.git//harness/m` | **passes** |
| (2) Write to `/r/.git//harness/m`, `/r/.git/./harness/m` | **passes** |
| (2) `.git/hooks/../harness/m` (asks: `.git/hooks/` is governed), `.git/harness//m`, `.git/harness/../harness/m`, `\.git/...`, `.g"it"/...`, `cd .git//harness && ...`, `cd .git/. && ...`; Write to `.git/hooks/../harness/m` | refused |
| (2) case: `> claude.md`, `> .Claude/agents/x.md`, `> Tools/Away.py`, `> docs/principles.md`, Write `/r/claude.md` (the Mac disk ignores case: each writes the governed file) | **passes** |
| (3) `cd $(git rev-parse --git-dir)/harness && touch x`; `rsync -a s/ $(git rev-parse --git-dir)/harness/` | **passes** |
| (3) `cp x` or `echo x >` to `$(git rev-parse --git-common-dir)/harness/m` | refused |
| (3) `cd $(git rev-parse --git-dir) && touch harness/x` | **passes** |
| (3) globs: `touch .git/harnes?/x`, `> .git/harnes[s]/x`, `touch .gi?/harness/x`, `> .*/harness/m` | **passes** |
| (3) symlink: `ln -s .git g && echo x > g/harness/m`; `ln -s $(git rev-parse --git-common-dir) g`; later line `echo x > g/harness/m` | **passes** |
| (3) script: `bash s.sh`, `python3 tools/x.py`, the script writing the folder | **passes** |

## 2. The rules

Two helpers first. **`known`**: the line's literal assignments,
`NAME=value` with a value free of quotes, `$`, backticks and separators
(`S=/private/tmp/.../scratchpad`); a `$NAME` or `${NAME}` it holds is
replaced by its value before any rule below reads a word. **`in_folder(p)`**:
`p` with each remaining expansion as a placeholder, lexically normalised
(`posixpath.normpath`), split into parts; true when two adjacent parts stand
for `.git` and `harness` (case folded) **and at least one of the two is
written out literally**. A part stands for its name when it is that name,
when it is an expansion, or, in the `globs` mode, when it is a glob that
matches the name by `fnmatch` (a `{...}` group read as `*`). So `$X/*`,
`$X/{a,b}` and `*/*` never match; `$G/harness`, `.git/harnes?` and
`.git//harness` do.

1. **zsh redirects** (`tools/shell_scan.py`). `REDIRECT` becomes
   `(\d*|&)(>>&?[|!]?|>&[|!]|>[|!]|>&|>|<<<|<<-|<<|<>|<&|<)`, and every
   operator starting `>` (or `<>`, or with prefix `&`) records its target as
   a write, as `>` does; `>&` with a digit or `-` stays a descriptor copy.
   The text rule `WRITES` (lines the parser cannot read) widens only by the
   new operators: `>>?\s` becomes `>>?&?[|!]?\s`, so `->` and `=>` read as
   they do today. A bash reading of `>!` (a file named `!`) is never safer
   than zsh's, so zsh's wins.
2. **Write targets.** In the target clause, in `only_reads` and in the
   Write arm, `in_folder` (globs mode) replaces each bare
   `HARNESS_PATH.search`; each target is first joined to the `cd`/`pushd`
   path so far. Closes all of (2), the glob rows of (3), and `cd $(git
   rev-parse --git-dir) && touch harness/x`. `cd $X && rsync s/ harness/`
   stays open: rsync names no target, and refusing it would refuse every
   `cd $X && <program>`.
3. **Words, per command.** A simple command one of whose own words passes
   `in_folder` (globs off: a quoted `find -path '*/.git/*'` is a pattern,
   not a path) must be on h33's short list (`allowed`) with its redirect
   targets plain and outside the folder; else the line is denied. The other
   commands on the line are judged as today. Closes `cd $(...)/harness &&
   touch x` (the `cd` walk of rule 2 then denies `touch`) and the `rsync` row.
4. **Links.** (a) An `ln` source (every non-option argument but the last)
   that, after `known`, is not under a scratchpad and ends in a literal
   `.git` part, holds `rev-parse ... --git-dir` (or `--git-common-dir`,
   `--absolute-git-dir`), or passes `in_folder`, is denied. The link name is
   not judged. (b) Each named write target and the Write arm's path are also
   judged by `os.path.realpath`, joined to the event's `cwd`, for the folder
   and for `governed()`: a link an earlier line made is followed.
   `judge_bash` takes `cwd` as a keyword, default the process's.
5. **Case** (main session's default, Ola's call open). `governed()`
   compares names case-folded (`GOVERNED`, prefixes, suffixes, globs); the
   standard-library-module rule stays exact (Python's import checks case on
   macOS). No tracked file becomes governed by this (`git ls-files | grep
   -i` over the names, minus the exact ones: none).

**Stays open.** A script, config or program whose own text builds the path
(`bash s.sh`, `python3 tools/x.py`, `make`): what a program writes is fixed
when it runs, and reading its file is defeated by `'har'+'ness'`. Also
`cd $X && <unknown writer> harness/...` (rule 2) and a path held in a
variable set to an expansion (`H=$G/harness; touch $H/x`). The one full
close is the operating system's: Claude Code's sandbox,
`sandbox.filesystem.denyWrite` on the folder, enforced by Seatbelt on child
processes too; hooks run outside it, so they still write their state. It
turns the sandbox on for every command (write access limited to the working
tree, network filtered), so it is a question for Ola, not part of this
change.

## 3. Cost over the transcripts

A scratch model of each rule beside today's `judge_bash`, run over the
92 968 distinct Bash commands in the session transcripts (every
`tool_use` named `Bash` under `~/.claude/projects/<this project>/`,
subagents included). Newly refused, per rule alone:

| Rule | Newly refused | What they are |
|---|---|---|
| 1 redirects and `WRITES` | 0 | |
| 2 and 3 folder | 8 | reads of the folder through `G=$(git rev-parse --git-common-dir); ... $G/harness/...` by `sed -n`, `python3 -I -` or a probe script: 3 shapes, each with and without a leading `cd`. Acceptable: the same reads with `.git/harness` written out are refused since h33; read with `cat`, `head`, `tail`, `rg` |
| 4a `ln` | 1 | `ln -sfn $R/.git $SP/lk/.git` with `R` the main checkout (a guard probe linking the real git folder into a scratchpad), asked today, denied now: the route this rule closes |
| 4b real path | 1 | `cd $S/t/tools && cp .../tools/check_citations.py .` in a scratchpad, now asked: the copy is named `tools/check_citations.py`, which asks today when written as one path in a multi-command line (h16: the scratchpad exemption holds for one plain command only) |
| 5 case | 0 | |

The first draft (expansions alone matching, globs on every word) refused
353 lines; the literal-part condition, `known`, per-command words and
`ln` sources only bring it to the 10 above.

## 4. Tests (red step)

Cases in `tests/python/test_guard_governance.py` and
`test_guard_targets.py`: every **passes** row of §1 is refused (asked for
(1) and the case rows, denied for the folder rows; the later-line link row
with a real link made in `tmp_path`), except the rows under *Stays open*,
pinned as passing with a comment naming this section; the model refuses all
38 such rows. Still passing: h33's 36-case pass list (not run through the model; the
suite reruns it); and, each passed by the model: `cd $CLAUDE_JOB_DIR/tmp/x && pytest`; `cd $(git rev-parse
--show-toplevel) && ls`; `echo $(git rev-parse --git-dir)/harness`; `make
2>&1 | grep x`; `echo x >&2`; `exec 3>&-`; `git grep -nE
'^(<<<<<<<|>>>>>>>)'`; `ln -s ../build b`; `ln -s $DATA data`; `rm $D/*`;
`find . -not -path '*/.git/*' -name x`; `S=/tmp/s; cd $S/harness && echo x
> y`; a Write to `/r/notes/claude-notes.md`. Every refused row of §1 stays
refused.

## 5. Size

About 40 counted lines: `shell_scan.py` 2 changed; the hook `known` 4,
`in_folder` 10, `cd` join 5, per-command words 5, links 9, case 3,
`WRITES` 0 net. No rule text changes.

## Review

- Design round 1: CHANGES REQUESTED, aabac548..650324b0, 0 net production lines (rules 2 and 3 refused 353 harmless transcript lines; three other new refusals uncosted).
- Design round 2: APPROVED, ..c04efa9e, 0.
