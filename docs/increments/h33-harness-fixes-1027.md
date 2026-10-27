# Harness h33: four harness fixes of 2026-10-27

Status: **designed**, 2026-10-27, on `worktree-harness-fixes-1027` off master
`29914d50`. Full path (`docs/increments/README.md`): about 80 counted lines,
over the short path's 50; no C++ file, no refine or mesh code. Every code file
it changes is governed (`tools/brief.py`, `tools/push_grant.py`,
`.claude/hooks/guard_governance.py`, `tools/shell_scan.py`,
`tools/check_citations.py`), so the green
step runs with Ola at the guard's prompt, by day.

Source: Ola's answers of 2026-10-27 to the numbered proposals: "1: yes 2: yes
3: ? 4: yes 5: yes 6: yes", and item 3 answered later in words (mentioning the
harness record folder is allowed; only writing there is refused). Items 1, 2,
3 and 6 are built here as fixes 1 to 4; item 4 is checked in §5. Open, defaults taken: the increment-file check covers the public repository only, and skips retrospective-only branches (§2).

## Prior art: legacy and literature

*Literature.* None applies: these are repository tooling rules, not a method
with a published form. No novelty is claimed.

*Legacy.* `git grep -n -i -E 'push_grant|check_citations|guard_governance|git/harness' legacy-archive -- legacy`
returns nothing; nothing is carried across.

*Ideas file.* No entry is taken up. "The citation gate skips absolute-path
cites" is related to fix 4 and already shipped (`ABSOLUTE` in
`tools/check_citations.py`).

## 1. An untimed @perf run may have a read-only agent beside it

Today `tools/brief.py@29914d50:175` refuses an untimed set unless it is
exactly `perf` plus one writer. New rule: with `--untimed`, the set holds
`perf`, at most one other writer, and any read-only personas (today only
`reviewer`); rules 2 to 4 of `h9-spawn-briefs.md` §3.1a still apply (a
read-only agent never shares a writer's worktree, `perf`'s included). The
refusal text for a second writer stays `an untimed @perf run goes beside one
writer at most`. A timing run (no flag) still runs alone. `h9`'s rule 1 text
is updated in the same PR. About 2 counted lines.

## 2. The push grant wants the review in the increment file

After `review_of` passes (`tools/push_grant.py@29914d50:386`), an act on
origin is also checked against the branch's increment file:

- **Finding the file.** The paths `git diff --name-only -z
  <remote>/master...<head>` lists (merge base, so master merged in does not
  count) that match `docs/increments/<name>.md`, one level deep,
  `README.md` excluded. Each is read at the head (`git show <head>:<path>`).
  Under the one-file-per-PR rule this is one file, but any of them may hold
  the record.
- **The record.** Its `## Review` section (from that heading to the next
  `## ` heading) holds a hexadecimal word of 7 to 40 digits that is a prefix
  of the recorded `Review:` sha. A sha in the Status line or elsewhere does
  not count.
- **Refusal**, one reason for every miss (no file, no section, no sha):
  `head <h7> is not reviewed: no increment file this branch changes names
  <r7> under ## Review; it may be redone once the review of <r7> is recorded
  there`. It is a "not reviewed" refusal, so the redo rule of
  `.claude/REQUIRED-READING.md`, *The harness*, already covers it.
- **Not checked** (the review log alone, as today): a branch whose every
  changed path is under `docs/retrospectives/` (a retrospective has no
  increment file), and acts on the private repository (8 of its last 10
  merges changed no file under `docs/increments/`; see Questions).

A wordings commit after an `APPROVED WITH WORDINGS` review still passes:
`review_of` passes it as prose, and the file at the head names the reviewed
sha. About 30 counted lines.

## 3. The guard judges writes to the record folder, not mentions

`.claude/hooks/guard_governance.py@29914d50:159` denies a parsed line that
names `.git/harness` anywhere and writes anything; that is the false refusal
(a heredoc into `session.md` quoting the path, a `tools/brief.py` task
quoting it). Four design rounds tried to list every way a program can
write; each review found another (path arguments, `cd` then any program,
`sed`/`awk` script words, `ed` reading its script on standard input). This
design stops listing writers: on a line that names the folder, only a short
list of programs that cannot write passes.

**When the rule applies.** The text names `.git/harness`, or a `cd`/`pushd`
argument has `.git` as a path component (`(^|/)\.git(/|$)`; this one
condition also closes `cd .git && touch harness/x`). Any other line is
judged exactly as today.

**Unchanged first.** The first two clauses of `judge_bash` stay (`.claude/hooks/guard_governance.py@29914d50:157-158`): a write target or
interpreter candidate matching `HARNESS_PATH` with expansions read as
`.git`, or a run of `away.py`, is denied on any line. The Edit/Write arm
stays. An unparseable line the rule applies to is denied (today it was
denied only when its text looked like a write).

**The rule.** A line the rule applies to passes only if all of these hold;
otherwise it is denied, as today:

1. **Every simple command is allowed**: its command word as written is the
   program (no wrapper such as `env`, `xargs`, `sudo`, `exec`), it has no
   leading assignment (`RIPGREP_CONFIG_PATH=`, `PYTHONPATH=`, `GIT_*=` change
   what a listed program runs), it is not `unknown`, and it is one of:
   - `echo`, `printf`, `cat`, `grep`, `head`, `tail`, `wc`, `ls`, `stat`,
     `diff`, `cut`: none has an option that writes a file or runs a program
     (checked on this machine: BSD `head`, `tail`, `wc`, `stat`, `cut`,
     `cat`, Apple `diff`; `echo` and `printf` are zsh builtins; `grep` is a
     zsh function that hands `--filter`, `--config`, `--save-config`,
     `--pager`, `--view` and `--format-open` to BSD `grep`);
   - `rg` without a word starting `--pre` or `--hostname-bin` (both run a
     program; ripgrep 14.1.1);
   - `git log`, `git show`, `git diff`, `git status`, the subcommand as the
     second word (no `-c`, `-C`, `--git-dir`), without a word starting
     `--output` (`--output=<file>` writes; git 2.55.0 refuses its
     abbreviations as ambiguous, probed);
   - `git commit` with a message word (`-m`, `-F`, `--message`, `--file`,
     or a short cluster holding `m` or `F`, such as `-qm`) and without
     `-e`/`--edit`: no editor runs;
   - `python3` or `python` with `tools/brief.py` as the second word, on a
     line with no `cd`/`pushd` (so the script is the governed one in the
     session's checkout); it writes only standard output and runs only
     `git rev-parse`;
   - `cd`/`pushd` with at most one argument, a plain path (below) not
     starting with `-`;
   - no command word at all (a bare redirect, `{ ...; } > f`).
   Left out, each because an option or its input writes or runs: `sort`
   (`-o`), `uniq` (a second operand), `tee`, `sed` (`w`), `awk`, `ed`, `ex`,
   `file` (`-C`), `jq` (not needed: `cat` and `rg` read), `find`, `rm`.
2. **Every write target is a plain path outside the folder**: each
   redirect target is plain (it matches `shell_scan.PLAIN_PATH` and holds
   no `!` or `|`; it does not start with `~` unless `~/`), and neither it
   nor it joined to the `cd`/`pushd` paths so far, lexically normalised,
   matches `HARNESS_PATH`. A glob target (`.git/harnes?/q`), an expansion
   (`"$N"`) and zsh's `>!` and `&>!` (read today as a target `!`) all fail
   it; zsh's `>>|` makes the line unparseable, so it is denied above.
3. **No `cd`/`pushd` leads into the folder**: the `cd`/`pushd` paths joined
   in order, normalised, do not match `HARNESS_PATH`. Reads after a `cd`
   into the folder are refused; read the file by its path instead.
4. **No expansion runs code the parser did not read**: the text holds no
   `${` other than `${NAME}` or `${N}`, no `$((` and no `$[` (zsh runs a
   command through the `(e)` parameter flag, probed); and no unquoted
   heredoc body holds `$(` or a backtick (the shell runs it; `shell_scan`
   reads the body as text). `$(...)` and backticks elsewhere are parsed,
   and their commands must pass 1 to 3.

`shell_scan` gains two fields on `Simple`, set where it already reads them:
`assigns` (the command had a leading `NAME=value`) and `runs_body` (it owns
an unquoted heredoc whose body holds `$(` or a backtick; `read_heredocs`
already knows the delimiter is unquoted). About 45 counted lines.

**Out of scope, follow-up.** Writes whose text names neither `.git/harness`
nor a `cd` into `.git` are judged as today, and some pass today and after
this change: a path built by an expansion and given to a program that is
not a listed writer (`cd $(git rev-parse --git-dir)/harness && touch x`,
`rsync -a s/ $(git rev-parse --git-dir)/harness/`); a glob that reaches the
folder (`touch .git/harnes?/x`); a symbolic link made by an earlier line; a
script or config file whose own text holds the path. Separately, zsh's
`>!` and `&>!` are misread by `shell_scan.REDIRECT` (target `!`), and
`>>|` leaves the line unparseable with no text-rule match (probed: zsh
writes the file each time), so `echo x >! CLAUDE.md` and `echo x >>|
CLAUDE.md` pass today without the governance prompt; a fix to `REDIRECT`
and to the text rule belongs to that follow-up.

## 4. Bare line citations in retrospectives are refused

`tools/check_citations.py@29914d50:441`: a bare line citation
(`path:line`, no `@commit`) whose source file is under
`docs/retrospectives/` is listed under `unpinned`, with the message
`pin it (path@<sha>:<n>)`; `legacy/` citations stay pinned to the tag. A
probe with that one condition on `29914d50` lists **43 bare citations in 24
files**. @orchestrator pins them in this PR, before the green commit: each
to the commit that last wrote its citing line (`git blame`), once the quoted
text is checked there; one it cannot pin so is reported. About 3 counted
lines.

## 5. Item 4: redo in the named form

Already complete: `.claude/REQUIRED-READING.md`, *The harness*, allows a
refused push, pull request or merge to be redone "in the form its reason
names", and `NOT_BARE` (`tools/push_grant.py@29914d50:56-58`) ends "it may be
redone in that form" for all three acts. Nothing more.

## 6. Tests (red step)

1. `test_brief.py`: perf(B) + developer(A) + reviewer(C) `--untimed` passes,
   on perf's and on reviewer's brief, the line listing both; perf + reviewer
   alone passes; reviewer sharing perf's worktree is refused (`is read-only
   and shares`); perf + two writers + reviewer is refused as today.
2. `test_push_grant.py`, a fixture branch off its origin/master: record at
   8 digits in `## Review` passes; 7 digits passes; 6 digits, the sha only
   in the Status line, no increment file changed, only `README.md` changed,
   and a file brought in only by a master merge are each refused with the
   exact reason; a wordings commit after the review passes; a
   retrospective-only branch and a private act pass on the log alone.
3. `test_guard_governance.py`. Pass list: a quoted and an unquoted
   heredoc into `session.md` naming `.git/harness/x` (the quoted one also
   holding `$(git rev-parse --git-dir)/harness`); `echo "see .git/harness" >
   note.md`; `printf '%s\n' .git/harness >> notes.md`; `git commit -m
   "guard: .git/harness"`, `-qm`, `-F msg.txt`; `python3 tools/brief.py a
   --task "the .git/harness/x file"`, also `> /tmp/brief.txt`; `cat`,
   `grep`, `rg` (also `> /tmp/hits`), `wc -l`, `cut`, `head | cut`, `tail`,
   `ls -la`, `stat`, `diff` of the folder or a file there; `git log
   --oneline`, `git show`, `git status`, `git diff --stat` beside a mention;
   `cd /abs/checkout && grep ...`; `cd .git && git status && echo
   .git/harness`; `echo $HOME .git/harness ${HOME} $1`; `cd $(git rev-parse
   --show-toplevel) && ls`. Deny list: the round-4 lines (`cd .git && touch
   harness/x`; after `cd .git/harness`: `python3 -c "open('x','w')"`, `rsync
   ... . > /tmp/log`, `tar xf ... > /tmp/log`, `vim x > /tmp/log`, `sed -n 'w
   x' in > /tmp/log`, `file -C -m m`; `cd .git && rsync -a /s/ harness/`;
   `python3 -c "...open(sys.argv[1],'w')" .git/harness/x`; `bash s.sh`,
   `python3 tools/x.py`, `python3 -m mod` on a path there; `rsync`, `tar
   -C`, `find -delete`, `ex -sc wq`, `vim`, `sort -o`, `uniq y x` there,
   each alone and with `> log`; the `sed ... w` pair and the `awk print >`
   line; `git apply p && cat .git/harness/x`; `N=.git/harness/x; echo x >
   "$N"`); the two `ed` lines (script on a pipe, script as a quoted
   heredoc); and this round's: `>!`, `>>!`, `&>!`, `>>|` into the folder;
   a glob target (`?`, `[s]`) beside a mention; `${(e):-'$(touch x)'}`; an
   unquoted heredoc body running `$(touch ...)`; a here-string `<<<
   "$(touch ...)"` and unquoted; the zsh glob qualifier `*(e:'touch ...':)`;
   `rg --pre`; `git log --output=`; `git -c ... log`; `RIPGREP_CONFIG_PATH=`
   and `PYTHONPATH=` prefixes; `cat <> .git/harness/x`; `cd .git && echo x
   > harness/y`; `cd .git/hooks && echo x > ../harness/y`; `cd /tmp/e &&
   python3 tools/brief.py ... .git/harness/x`; `env rg`; `| xargs touch`;
   `| sh`; `> ~/r/.git/harness/y`. The existing denials
   (`test_guard_targets.py` lines with `-harness`, the two harness lines of
   `DENIED_COMMANDS`) stay green unchanged.
   Measured at `f21119e7` (hook as at `29914d50`) with a scratch probe that
   loads the hook and models §3 on `shell_scan.parse` (the round-4 probe was
   not kept; its list is rebuilt from this file, 32 lines where round 4
   counted 36). Deny, 66 cases: today's `judge_bash` refuses 37 and passes
   29; the model refuses 65. The one it passes, `cd $(git rev-parse
   --git-dir)/harness && touch x`, names neither trigger and is in the
   follow-up. No case today refuses passes the model. Pass, 36 cases: today
   passes 29 (7 false refusals); the model passes 28. Former passes the
   model now refuses, by design: `sort` and `jq .` of a file there (not
   listed: use `cat`, `rg`), every read after `cd .git/harness` (rule 3:
   read by path), and `rm n.md && cat >> session.md` heredoc (refused today
   too: `rm` is not listed, as `rm -r` of an ancestor deletes the folder;
   run the `rm` as its own command). zsh confirmed, in a scratch folder,
   that the `(e)` flag, the glob qualifier, the unquoted heredoc body, `>!`,
   `&>!`, `>>|` and `<>` each run or write.
4. `test_check_citations.py`: a bare cite in a retrospective fails, a pinned
   one passes, a `legacy/` one passes, a bare cite elsewhere under `docs/`
   behaves as today.

## 7. Rule text and size

Words (`tools/rule_sizes.py`), applied by @architect after the design
review, with Ola's yes at the guard's prompt:

- `.claude/agents/perf.md` line 46: "may run beside one writer." becomes
  "may run beside one writer and read-only agents." (+3)
- `docs/increments/README.md` step 4: "a hook records it, and the away-mode
  push grant reads only that record." becomes "a hook records it, and the
  away-mode push grant also wants that sha in the `## Review` of the
  branch's increment file." (+9); cut "Check: `grep -n '^## Review'
  docs/increments/NN-name.md` on the branch." (-9; the grant now checks
  more than a heading) and "(CI is authoritative, `CLAUDE.md` §4)" (-5; it
  restates `CLAUDE.md`).

Net **-2** words. Counted lines: about 80 against the limit of 700.

## Review
