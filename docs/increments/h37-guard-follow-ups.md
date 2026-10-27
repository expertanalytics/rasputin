# Harness h37: governance guard follow-ups of h33

Status: **designed**, 2026-10-27, on `worktree-guard-zsh` off master
`aabac548`. Short path (`docs/increments/README.md`): about 30 counted lines,
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

1. **zsh redirects** (`tools/shell_scan.py`). `REDIRECT` becomes
   `(\d*|&)(>>&?[|!]?|>&[|!]|>[|!]|>&|>|<<<|<<-|<<|<>|<&|<)`, and every
   operator starting `>` (or `<>`, or with prefix `&`) records its target as
   a write, as `>` does; `>&` with a digit or `-` stays a descriptor copy.
   The text rule `WRITES` reads any `>` (today `>>?\s`). A bash reading of
   `>!` (a file named `!`) is never safer than zsh's, so zsh's wins. Cost:
   2 of the 690 594 Bash command lines in the transcripts (first 60
   characters) hold these forms,
   both inside quoted patterns, which the parser already reads as text.
2. **One folder test.** `in_folder(p)`: `p` with expansions read as `.git`,
   lexically normalised (`posixpath.normpath`), split into components; true
   when two adjacent components match `.git` and `harness`, case folded,
   each by `fnmatch` (so a glob or a brace counts as matching what it could
   match; a `{...}` group is read as `*`). Used by the target clause, by
   `only_reads`, by the `cd` path walk and by the Write arm, in place of each
   bare `HARNESS_PATH.search`. Closes all of (2) and the glob rows of (3).
3. **Words, not just text.** `names_folder` is also true when any word or
   write target of any simple command passes `in_folder`. Closes
   `cd $(...)/harness && touch x` and the `rsync` row (rsync is not on the
   short list).
4. **`cd` into an expansion.** In the target clause, targets are joined to
   the `cd`/`pushd` path so far (expansions read as `.git`), so `cd
   $(git rev-parse --git-dir) && touch harness/x` is denied. Only named
   targets: `cd $X && rsync s/ harness/` stays open, because closing it
   would refuse `cd $X && <any program>`, which agents run (268 of the 690 594
   transcript command lines hold a `cd $...`, mostly `cd
   $CLAUDE_JOB_DIR/tmp/...`).
5. **Links.** (a) An `ln` argument that, with expansions read as `.git` and
   normalised, ends in a `.git` component or passes `in_folder` is denied.
   (b) Each named write target and the Write arm's path are also judged by
   `os.path.realpath`, joined to the event's `cwd`: a link an earlier line
   made is followed, for the folder and for `governed()`. `judge_bash` takes
   `cwd` as a keyword, default the process's.
6. **Case.** `governed()` compares names case-folded (`GOVERNED`, prefixes,
   suffixes, globs); the standard-library-module rule stays exact (Python's
   import checks case on macOS). No tracked file becomes governed by this
   (`git ls-files | grep -i` over the names, minus the exact ones: none).

**Stays open.** A script, config or program whose own text builds the path
(`bash s.sh`, `python3 tools/x.py`, `make`): what a program writes is fixed
when it runs, and reading its file is defeated by `'har'+'ness'`. Also
`cd $X && <unknown writer> harness/...` (rule 4). The one full close is the
operating system's: Claude Code's sandbox, `sandbox.filesystem.denyWrite` on
the folder, enforced by Seatbelt on child processes too; hooks run outside
it, so they still write their state. It turns the sandbox on for every
command (write access limited to the working tree, network filtered), so it
is a question for Ola, not part of this change.

## 3. Tests (red step)

`test_shell_scan.py`-style cases in `tests/python/test_guard_governance.py`
and `test_guard_targets.py`: every **passes** row of §1 is refused (asked for
(1) and the case rows, denied for the folder rows), except the rows under
*Stays open*, which are pinned as passing with a comment naming this section.
Still passing: h33's 36-case pass list unchanged; `cd $CLAUDE_JOB_DIR/tmp/x
&& pytest`; `cd $(git rev-parse --show-toplevel) && ls`; `echo
$(git rev-parse --git-dir)/harness`; `2>&1 | grep`, `>&2`, `>&-`;
`git grep -nE '^(<<<<<<<|>>>>>>>)'`; `ln -s ../build b`; a Write to
`/r/notes/claude-notes.md`. Every refused row of §1 stays refused.

## 4. Size

About 30 counted lines: `shell_scan.py` 2 changed; the hook `in_folder` 8,
words 3, `cd` join 5, links 7, case 3, `WRITES` 1. No rule text changes.
