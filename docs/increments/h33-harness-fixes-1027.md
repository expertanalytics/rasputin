# Harness h33: four harness fixes of 2026-10-27

Status: **designed**, 2026-10-27, on `worktree-harness-fixes-1027` off master
`29914d50`. Full path (`docs/increments/README.md`): about 85 counted lines,
over the short path's 50; no C++ file, no refine or mesh code. Every code file
it changes is governed (`tools/brief.py`, `tools/push_grant.py`,
`.claude/hooks/guard_governance.py`, `tools/check_citations.py`), so the green
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
(a heredoc into `session.md` quoting the path). The new deny test for a
parsed line, replacing that clause:

- (a) unchanged: a write target or interpreter candidate, with expansions
  read as `.git`, matches `HARNESS_PATH`; or the line runs `away.py`;
- (b) traced through `cd`/`pushd`: each relative write target and each
  relative interpreter candidate is also judged joined to each directory
  argument the line changes to, and to all of them in order; a joined path
  is judged as (a) judges a path, each expansion read as `.git` (denied:
  `cd .git && touch harness/x`, `cd .git/harness && python3 -c
  "open('x','w')"`, `cd $(git rev-parse --git-dir)/harness && touch x`).
  When the directories changed to so far, joined in order with each
  expansion read as `.git`, match `HARNESS_PATH`, every later command not
  on the read-only list below is denied, whatever it names (`cd
  .git/harness && rsync -a /s/ . > /tmp/log`, `... && tar xf /t.tar >
  /tmp/log`, `... && vim x > /tmp/log`: no path of theirs matches alone);
- (c) arguments: in every simple command except a program that only reads,
  a word that is one path (no whitespace) matching `HARNESS_PATH`, as
  written or joined to the directories as in (b), counts as a write target
  (`cd .git && rsync -a /s/ harness/`). The programs that only read:
  `shell_scan.NO_RUN` less its writers (`tee`, `mkdir`, `rm`, `cp`, `mv`,
  `touch`) and less `file` (`file -C` writes), plus `jq`, `cd` and
  `pushd`; `sort` with `-o` or `--output` and `uniq` with two operands
  write, so they do not count as read-only. A program that takes a script
  of its own (`sed`, `awk`, `gawk`, `mawk`, `nawk`, `ed`, `ex`, `vi`,
  `vim`, `nvim`) has each word split on whitespace, quotes, braces,
  parentheses, `;`, `,`, `|`, `<` and `>`, and each piece is judged as a
  one-path word: a script's `w FILE` or `print > "FILE"` names its target
  inside a word with spaces. For an interpreter this covers the script's
  words, `-m`'s words and the words after a `-c` program (`python3 -c
  "…open(sys.argv[1],'w')" .git/harness/x`). Denied with or without
  another write on the line: `rsync -a s/ .git/harness/`, `tar xf t -C
  .git/harness`, `find .git/harness -delete`, `bash s.sh .git/harness/x`,
  `ex -sc wq .git/harness/x`, `vim .git/harness/x`, `python3 tools/x.py
  .git/harness/y`, `sort -o .git/harness/x y`, `uniq y .git/harness/x`;
  denied too: `sed -n 'w .git/harness/x' in > log`, `sed
  's/a/b/w .git/harness/x' in > log`, `awk '{print > ".git/harness/x"}' in
  > log`, and after `cd` into the folder `sed -n 'w x' in > /tmp/log` and
  `file -C -m m` (by (b)). Passing: `wc -l`, `sort`, `cut`, `rg` and `jq
  .` on a file there, also after `cd` into it; `sed` there is refused
  (`cat`, `head`, `tail` and `rg` read it);
- (d) untraceable, only when the text names `.git/harness`: a simple command
  is `unknown`; a write target holds an expansion (the `H=`, `export`,
  `read`, `declare` and `$(...)` routes); or, on a line that writes
  something, a `cd`/`pushd` argument holds an expansion.

A quoted mention holds spaces, so it is never one path and (c) passes it.
The unparseable-line rule (`:173`) and the Edit/Write arm stay. A Python
program naming the folder in a literal it does not open for reading stays
refused under (a): this design's own probe hit it (a `glob.glob` read);
`grep`, `cat` or the Read tool read the folder. About 40 counted lines.

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
3. `test_guard_governance.py`, a pass list: `rm n.md && cat >> session.md`
   heredoc naming `.git/harness/x`; `echo "see .git/harness" > note.md`;
   `git commit -m "guard: .git/harness"`; `python3 tools/brief.py a --task
   "the .git/harness/x file"`; `cat` and `grep` of a file there; `cd
   $(git rev-parse --show-toplevel) && ls`; every passing read of §3 (c),
   plain and after `cd`. A deny list: every example of §3 (b) and (c), the
   six non-interpreter (c) programs each alone and with `> log`; `python3
   -m mod .git/harness/x`; `git apply p` with the path in its text; `echo
   x > "$N"` with the path in its text. The existing denials
   (`test_guard_targets.py` lines with `-harness`) stay green unchanged.
   Measured at `29914d50` with a scratch probe that loads the hook and
   models §3 on `shell_scan.parse`: 36 deny cases, 20 refused by today's
   `judge_bash` and 16 passing today (the holes), all 36 denied by the
   model; 14 pass cases, 2 refused today (the false positives), all 14
   passed by the model. No mismatch. Without the script-word split, the
   two `sed ... w` lines and the `awk` line pass the model: taking `sed`
   off the read-only list alone does not close them.
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

Net **-2** words. Counted lines: about 60 to 70 against the limit of 700.

## Review
