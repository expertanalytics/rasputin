# Harness h4: guard fixes after the first unattended nights

Status: **design**, @architect, 2026-10-01. Ola asked for it on 2026-10-01
("let's do the guard fix increment now") and kept it small. Builds on
`docs/increments/h3-unattended-u1.md` (*U1*). Design only; one PR.

## 1. Prior art

Tooling; no novelty claimed. Codex splits a shell command into simple commands
with tree-sitter-bash only when it is plain words and operators, and treats
anything with redirection, substitution or expansion as one opaque command
(learn.chatgpt.com/docs/agent-configuration/rules). h4 parses further
(redirections, heredocs, `$(...)`) because the false positives below are
exactly those. bashlex (GPLv3) and tree-sitter are not options: the hooks run
under whichever `python3` is on `PATH`, so stdlib only.
*Legacy*: nothing.

```
$ git grep -l -iE 'shlex|bashlex|heredoc' legacy-archive -- legacy
(no output, exit 1)
```

## 2. The false positives, reproduced

Probe: `docs/increments/h4-probes/h4_probe.py`. It copies the hooks into a
temporary repository with the U1 fixture (`tests/python/harness_fixtures.py`),
feeds each guard the PreToolUse event, and prints the decision by day (flag off)
and at night (flag on). Given the night's queue file it replays those acts too:
`python3 docs/increments/h4-probes/h4_probe.py <git-common-dir>/harness/queue-2026-10-01.jsonl`.
Every case below gives `ask` by day and `deny` at night, except item 5, which is
`deny` in both modes.

| # | Case (probe label) | Cause in the current code |
|---|---|---|
| 1 | `printf ... > .claude/current-task/reviewer-roadmap.md; cat .claude/REQUIRED-READING.md` (queue line 1, `fp1-bash`) | Bash arm: some write construct anywhere (`> `) plus a governed name anywhere (the `cat`). The Write tool on the same path passes (`fp1-write`); the refused act was this Bash form |
| 2 | `cat > session.md <<'EOF' ... ASK OLA: git push ... EOF` (queue line 3, `fp2-*`) | `guard_push`'s `\bgit\b[^|;&]*\bpush\b` regex runs over heredoc bodies and quoted text |
| 3 | `sed -i` or a `python3 - <<EOF` editing `docs/increments/23-basin-scale.md` whose text names `CLAUDE.md`, `docs/increments/README.md` or runs `tools/check_citations.py` (queue lines 4, 7, 8; `fp3-*`) | as 1: a governed name in prose or in a command that only reads it |
| 4 | `git show b:.claude/settings.json \| python3 -c ...` (queue line 9); `python3 tools/check_citations.py > /dev/null; git diff > $CLAUDE_JOB_DIR/tmp/diff.txt; python3 - <<EOF` (queue lines 10, 11; `fp4-*`) | `python3 -c` and `> /dev/null` count as writes; the settings glob matches a `git show` argument |
| 5 | `grep -n 'def _open_tty\|/dev/tty' tools/away.py \| head`; `echo '... away.py + guards; ...' >> MEMORY.md`; `printf 'Ask: ... away.py ...; commit\n' > .claude/current-task/x.md` (`fp5-*`, from refused calls in the 2026-09-30/10-01 transcripts) | `SEGMENTS` splits on `\|`, `;` and newlines **inside quotes**; the fragment then fails `shlex`, falls back to `str.split`, and a fragment whose first word is not a reader "runs" `away.py` |

Queue lines 5 and 6 replay as `pass` only because the queue truncates `act` to
1000 characters and the governed name came later; they are the same case as 3.
The true positives in the same queue (line 0 `git rebase`, line 2 `git push`)
stay `ask`/`deny`, as do the probe's `tp-*` cases.

**One root cause.** Both guards read the command as a string: the governance
arm asks "is there a write construct?" and "is there a governed name?"
independently, and the segmenter does not know quoting or heredocs.

## 3. Design: judge targets, not words

**`tools/shell_scan.py`** (new, stdlib, pure, joins the self-protecting
governed set of U1 §3.6). `parse(command) -> list[Simple] | None`, where
`Simple(argv: list[str], writes: list[str], program: str | None)`:

- a quote-aware lexer (single, double, backslash) for `&&`, `||`, `;`, `|`,
  `&`, newline, `(`, `)`, `{`, `}`; leading reserved words (`if`, `then`, `do`,
  `else`, `!`, ...) and `VAR=x` assignments are skipped;
- redirections: `>`, `>>`, `>|`, `&>`, `N>` add their target to `writes`, except
  `/dev/null`, `/dev/std*` and `&N`; `<`, `<<<` read;
- heredocs (`<<`, `<<-`, quoted or not): the body is consumed, never lexed as
  commands; it becomes `program` when `argv` is an interpreter reading stdin;
- `$(...)`, backticks, `<(...)` are parsed recursively and their simple commands
  added; wrappers (`env`, `nohup`, `time`, `timeout N`, `nice`, `command`,
  `exec`, `xargs`, `script -q F`, `sudo`) are stripped; `sh|bash|zsh -c S`
  parses `S` recursively;
- writers by argv: `tee`, `cp`/`install`/`ln` (destination), `mv` (all
  operands), `rm`, `rmdir`, `unlink`, `touch`, `truncate`, `mkdir`, `sed -i`
  (files), `dd of=`, `git checkout|restore -- <paths>`; `git apply`, `patch`
  write unknown targets;
- interpreters (`python`, `python3`, `python3.N` by basename; `perl -e/-i`,
  `node -e`, `ruby -e`): `program` is the `-c`/`-e` text or the heredoc. For
  Python, the candidate targets are the string literals (`ast`) that are a
  single path, with no whitespace; prose literals are not paths. Other
  languages, or Python that does not parse: whitespace-separated tokens.
- `None` when it cannot parse: unbalanced quotes or `$(`, unterminated heredoc.

**Verdicts.**

- `guard_governance`, Bash: ask if any target of any simple command is governed
  (`governed()` unchanged). Always deny if any target matches `.git/harness`
  (an expansion counts as matching `.git`), or a simple command **runs**
  `away.py` (argv[0], or the script operand of an interpreter, has that
  basename; `-m away`). `cat`, `grep`, `git`, `pytest` naming it are reads by
  construction, so the `READERS` list goes.
- `guard_push`: the U1 rules and the six `PUBLISHES` regexes are applied to each
  simple command's `argv`, never to text: `git ... push`, `gh pr create|...`,
  `--no-verify` as a git argument, rebase/reset --hard/filter-branch/commit
  --amend, plus U1's `segment_why`. A `gh pr create` still asks; the words in
  its body or body file do not add a second reason.
- **Cannot judge.** The whole command unparseable, or a write target that is
  wholly an expansion (`> $OUT`, `rm $f`), or a writer with unknown targets
  (`git apply`): the guard falls back to today's text rule for that command
  (governed name plus write construct; the push regexes; `runs_away` on the
  text). That fallback is an `ask` by day and a refused, queued act at night,
  as U1. An expansion with a static tail (`$D/tmp/diff.txt`, `$D/CLAUDE.md`) is
  judged by the tail.
- The ask and the refusal name the targets: `why` becomes
  `it writes <governed targets>` (or `the guard cannot read this command's
  targets`), so the agent and Ola see what was judged.

**The rule on redoing a refused write** (rule text for
`.claude/REQUIRED-READING.md`, *The harness*): *A guard judges the file
written, whatever tool writes it, so the Edit and Write tools meet the same
guard as Bash. A refused write of a governed file is not retried by any route.
A refusal whose named targets are all ordinary files is a guard false positive:
redo it with Edit or Write, and put a `GUARD FALSE POSITIVE: <command>` line in
the handback (main session: in `session.md`).* The same paragraph's "all of them
read the command as text" becomes "they read the commands a shell line runs and
the files it writes; a line they cannot parse is judged as text".

Not in h4: the idle-time trace (deferred to its own increment); R-B; U2/U3.

## 4. Tests @tester writes red

`tests/python/test_guard_targets.py`, using the U1 fixture; each case flag off
and on. Commands are the probe's `CASES` strings, verbatim.

| Case | Day | Night |
|---|---|---|
| `fp1-write`, `fp1-bash` | pass | pass |
| `fp2-heredoc`, `fp2-body-file`, `fp2-commit-msg` | pass | pass |
| `fp3-sed`, `fp3-heredoc`, `fp3-commit-msg` | pass | pass |
| `fp4-git-show`, `fp4-diff-to-tmp` | pass | pass |
| `fp5-grep`, `fp5-echo-memory`, `fp5-printf-task` | pass | pass |
| `tp-push`, `tp-rebase`, `gh pr create --body 'x git push'` (one reason) | push ask | push deny, queued |
| `tp-claude-md-sed`, `tp-readme-heredoc`, `tp-settings-py`, `tp-check-cp`, `tp-edit-claude`, `rm CLAUDE.md`, `echo x > $D/CLAUDE.md`, `bash -c "sed -i s/a/b/ CLAUDE.md"` | governance ask | deny, queued |
| `tp-run-away`, `tp-harness-write`, `env python3 tools/away.py 8h`, `echo "$(python3 tools/away.py --back)"` | deny | deny, not queued |
| cannot judge: `echo 'unbalanced > CLAUDE.md`, `for f in CLAUDE.md; do sed -i s/a/b/ $f; done` | governance ask | deny, queued |
| cannot judge, no governed name: `rm $f`, `git apply x.patch` | pass | pass |

Plus unit tests of `shell_scan.parse` on the same strings (argv, writes, `None`
for the unparseable ones). U1's T4 and T7 stay green unchanged (T7's reader
rows now pass by construction). Not invariant-critical: no mutation round.

## 5. Estimate

One PR. Production lines (`CLAUDE.md` §2): `tools/shell_scan.py` ~170;
`guard_governance.py` +~25 net (Bash arm on targets, `READERS`/`runs_away`
replaced, fallback kept); `guard_push.py` +~20 net (argv checks, fallback);
the governed set +1. **~220**, far under 700. No `.claude/settings.json`
change. Rule text in `.claude/REQUIRED-READING.md` (prose, not counted).

## 6. Questions for Ola

1. **Cannot judge = today's text rule** (§3): an unparseable line asks or is
   refused only if it names a governed file, as now, rather than always.
   Recommendation: yes; always asking would turn every `for` loop over `$f`
   into a night refusal.
2. **The redo rule** (§3): an agent may redo a refused write of ordinary
   files with Edit or Write, reporting `GUARD FALSE POSITIVE:`.
   Recommendation: yes; the refusal names its targets, so the agent can tell.
