# Harness h29: the reviewer's verdict line, redoing a cured refusal, and amends after `cd`

Status: **designed**, 2026-10-21, on `worktree-h29-review-rules` off master
`5d8904ee`. Short path (`docs/increments/README.md`): about 30 counted lines,
no C++ file, no refine or mesh code. The code steps change governed files
(`tools/brief.py`, `tools/push_grant.py`, `tools/harness_mode.py`,
`.claude/hooks/guard_push.py`, `.claude/hooks/gates_after_commit.py`), so
they run with Ola at the prompt, by day.

Source: Ola's yes on 2026-10-21 to five proposals, taken at their defaults:
proposals 1 to 4 of the retrospective `docs/retrospectives/2026-10-21-merges-2.md`
(PR #325, branch `worktree-retro-1021c`), and a fifth, that the guard lets an
agent amend a commit that was never pushed.

**The picture.** Five small changes. (1) Every `@reviewer` brief demands the
one line the review hook records. (2) An away-mode grant refused only because
the head had no review may be redone once a review is recorded, and the
refusal says so. (3) `CLAUDE.md` stops restating who applies review wordings
and points to the one rule. (4) A design amendment recomputes the figures of
the tests it changes. (5) The amend exemption of h27 also covers the line
agents actually write, `cd <path> && git add … && git commit --amend …`, and
every message that tells an agent to amend names a form that passes.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature* (git 2.55.0, `git --version`).

- POSIX Shell Command Language, 2.9.3 *Lists*: in an AND-list
  `a && b`, `b` runs only when `a` exits 0. Probed here: `zsh -c 'cd
  /nonexistent-dir && echo RAN'` and the same in `bash` print the `cd` error,
  no `RAN`, exit 1. So a failed `cd` never lets the amend run in the wrong
  directory.
- `git help add`: "Add contents of new or changed files to the index". It
  writes the index, not HEAD. Probed in a scratch repository: HEAD's sha is
  the same before and after `git add f`.

*Legacy.* Nothing: the legacy tree is the old meshing product, with no
harness.

```
$ git grep -l -i -e 'amend' -e 'reviewer' -e 'verdict' legacy-archive -- legacy
(no output, exit 1)
```

*The ideas file* (`../rasputin_scratch/ideas/ideas.md`, read 2026-10-21): no
entry on review verdicts, refusal redos or amends; nothing taken up or
overlapped. The nearest, "Cap on runner and shell words per line in the push
guard", is about the guard's run time and is untouched.

## 2. What h27 already covers, and why amends were still refused

h27 (`docs/increments/h27-unpushed-amend.md`, built: `unpushed_amend` in
`.claude/hooks/guard_push.py`) passes an amend whose HEAD no remote-tracking
ref or tag holds, but only on a **bare one-command line**:
`git [-C <absolute path>] commit … --amend …`, nothing else on the line. It
asks for `cd <repo> && git commit --amend` on purpose (its test row
`9-cd-first`).

Today's three amends were all on commits never pushed, and all chained:

| When (UTC) | Who | Line, shortened | Outcome |
|---|---|---|---|
| 05:24 | `@orchestrator` (retro-1021) | `cd … && sed … && grep … && python3 tools/check_citations.py …; git commit -q --amend --no-edit -a && git log --oneline -1` | asked; Ola was at the prompt |
| 09:18:54 | `@orchestrator` (retro-1021c) | `cd … && python3 tools/check_citations.py … \| grep … \| tail -3 && git commit -q --amend -a --no-edit && git log … && git status …` | refused, away |
| 09:24:35 | `@architect` (lsm-ideas-warm-rain) | `cd … && git add docs/ideas/research/… && git commit -q --amend --no-edit && git log --oneline -1 && git status --short \| head` | refused, away |

Sources: `<git-common-dir>/harness/queue-2026-10-21.jsonl` (the two
refusals); the agents' transcripts for the 05:24 ask. Unpushed: the first
grant push of `worktree-retro-1021c` was at 09:23:16 and of
`lsm-ideas-warm-rain` at 09:31:24 (`grants-used.jsonl`), both after the
refusals. Probed with the hook itself on a scratch repository with one
unpushed commit: `git -C <abs> commit --amend --no-edit` passes (no output),
`cd <abs> && git commit --amend --no-edit` asks.

So the exemption works; the lines do not fit it. Two causes: the agents'
habit of starting every line with `cd` (their working directory resets
between calls), and the after-commit hook's "Fix and amend", which names no
form. What is missing is a wider accepted form that is still safe, and
messages that name it.

## 3. The changes

### 3.1 `tools/brief.py`: the verdict line (proposal 1)

For persona `reviewer` only, the body gains one line after the write limit:

```
End the handback with exactly one line `Review: <APPROVED | APPROVED WITH WORDINGS | CHANGES REQUESTED> <40-hex sha of the head you reviewed>`, and start no other line with `Review:`; without it no review is recorded (docs/increments/README.md step 4).
```

A module constant `REVIEW_LINE` holds it. The rule stays stated once, in
README step 4; the brief repeats the format because the hook
(`.claude/hooks/record_review.py`, regex `REVIEW` in `tools/push_grant.py`)
reads exactly that line. About 3 lines.

### 3.2 Redoing a refusal whose cure is done (proposal 2)

- `tools/push_grant.py`: the refusal for an unreviewed head (the
  `review_of` branch of `judge`, reason `head <sha7> is not reviewed: …`)
  appends `; it may be redone once @reviewer's verdict approving <sha7> is
  recorded`. Only that branch: the other refusals (already merged in this
  window, wrong repository, checks not green, …) are unchanged.
- `tools/harness_mode.py`: `QUEUE_REASON`'s "Do not retry it and do not
  work around it." becomes "Do not retry it or work around it unless its
  reason says it may be redone." The guard judges the redo again, so a redo
  before the review is recorded is refused again.

About 2 lines.

### 3.3 The amend exemption widened (item 5)

`unpushed_amend(command, cwd)` in `.claude/hooks/guard_push.py` accepts,
besides today's bare form, a line of simple commands **joined only by
`&&`**, in this order:

1. at most one `cd <absolute path>`, first, argv exactly `["cd", path]`;
2. any number of `git add …`, argv starting exactly `["git", "add"]`;
3. the amend, **last**, in today's bare form (`git [-C <absolute path>]
   commit … <AMENDS word> …`).

Each part, split at `&&` outside single quotes, must pass today's per-line
checks: `EXPANDS_NOTHING.fullmatch(part)`, and `shlex.split(part)` equals
the argv `shell_scan.parse` gives for that part, with as many parts as
simple commands (so a pipe, `;`, `||`, `&`, a redirect or a substitution
anywhere makes the counts or the match fail, and the line asks). The
repository judged is the amend's `-C` path if given, else the `cd` path,
else `cwd`. Then today's two `git` reads decide, unchanged.

Why it is safe: before the amend only `cd` (moves no ref) and `git add`
(writes the index, §1) run, and `&&` stops the line at the first failure
(§1), so HEAD when the guard looks is HEAD when the amend runs. The amend is
last because `reasons` is deduplicated: a later `git rebase` on the same line
would add the same history reason and hide behind the exemption. Anything
else before or after the amend asks, as today.

**The messages.** When the line's only reason is the history rewrite, it
holds an amend word, and the exemption is false, the guard's reason gains a
constant `AMEND_FORM`:

```
; an amend of a commit no remote holds passes when it comes last, after at most `cd <absolute path>` and `git add` commands joined by `&&`, or as `git -C <absolute path> commit --amend …` alone; it may be redone in that form
```

`.claude/hooks/gates_after_commit.py`: "Fix and amend -- do not push or open
a PR on this." becomes "Fix, then amend with `git -C <that tree> commit
--amend --no-edit -a` on a line of its own -- do not push or open a PR on
this." The module docstring of `guard_push.py` replaces its "on a bare
one-command line" sentence with the widened form.

About 20 lines.

### 3.4 Prose, in this design commit (proposals 2, 3 and 4)

- `.claude/REQUIRED-READING.md`, the unattended-mode paragraph: "A refused write to a
  rule file, or a refused push, is never redone by any route." becomes the
  rule for the two cured refusals: a refused rule-file write is never
  redone; a refused push, pull request or merge only when its sole reason
  is an unreviewed head, and only once an approving `@reviewer` verdict
  covering it is recorded; a refused amend in the form its
  reason names; the guard judges each redo again.
- `CLAUDE.md` §3, *Prose edits while writers are busy*, becomes *Prose
  edits*: a pointer to README step 4 for wordings and review records, the
  append-only merge kept under its two-busy-writers condition, and "composes
  no text in a persona's file and writes no code" kept.
- `.claude/agents/architect.md` §5, *Fixtures*: an amendment recomputes the
  figures of each existing test it changes.

## 4. Tests

Short path; no invariant-critical suite, so no mutation round. Fixtures are
built by the tests with `git` in `tmp_path`, as today.

**`tests/python/test_brief.py`**

1. The `reviewer` block contains `REVIEW_LINE` and the literal
   `` `Review: ``; the blocks of the five other personas do not.

**`tests/python/test_push_grant.py`** (fixture: today's granted branch with
one commit ahead of master, no record in `reviews.jsonl`)

2. The no-record refusal contains both `no review recorded for this branch`
   and `may be redone once @reviewer's verdict approving`.
3. With a `CHANGES REQUESTED` record on the head, the refusal carries the
   redo text too (it is a "not reviewed" refusal).
4. The "already merged in this window" refusal does not carry it.

**`tests/python/test_harness_mode.py`**

5. Existing test `test_decide_ask_while_unattended_denies_with_the_queue_reason`
   pins the phrase `"Do not retry"`. The new text still contains it, so the
   row holds; add `"unless its reason says it may be redone"` to its phrase
   list.

**`tests/python/test_guard_push.py`** (fixture: h27's — `make_repo`, a bare
`origin`, `git push origin master`, then `local_commit`: one pushed commit,
one unpushed HEAD)

Pass, silent in modes `off`, `on`, `expired`, nothing queued:

6. `cd <abs repo> && git commit --amend --no-edit`, event `cwd` the parent
   directory (another repository or none: the `cd` path is judged).
7. `cd <abs repo> && git add f && git add g && git commit -q --amend --no-edit -a`.
8. `git add f && git -C <abs repo> commit --amend --no-edit`, cwd the repo.
9. `cd <abs repo> && git commit --amend -m 'a && b'` (`&&` inside single
   quotes is not a split).

Ask by day, refused and queued while away, reason contains the history text
**and** `AMEND_FORM`:

10. Pushed HEAD: `cd <abs repo> && git commit --amend --no-edit`.
11. `cd rel && git commit --amend --no-edit` (relative `cd`).
12. `cd <abs> && cd <abs> && git commit --amend` (two `cd`).
13. `cd <abs repo>; git commit --amend`, `cd <abs repo> || git commit
    --amend`, `git add f | cat && git commit --amend` (other joins).
14. `cd <abs repo> && git checkout -q HEAD~0 && git commit --amend` (a
    command other than `git add` before the amend).
15. `cd <abs repo> && git commit --amend --no-edit && git log --oneline -1`
    (not last).
16. `cd <abs repo> && git commit --amend && git rebase master` (a second
    rewrite after the amend: the reason list is one entry, and it still asks).

17. Not an amend: `git rebase master` alone asks with the history reason and
    **without** `AMEND_FORM`.

**Existing rows that change** (recomputed, item 4): h27's `9-cd-first`
(`cd {repo} && git commit --amend --no-edit`, unpushed HEAD) moves from
`H27_ASKED` to the passed set; its pushed-HEAD twin is row 10. The other
`9-*` rows still ask; `@tester` checks each against §3.3 and says in the red
commit which moved.

**`tests/python/test_gates_after_commit.py`**

18. A red gate's text contains `commit --amend --no-edit -a` and `git -C`, and
    no longer `Fix and amend`.

## 5. Left out

- **Commands after the amend** (`&& git log`, `| head`). They cannot change
  what the amend replaced, but the guard's reasons are deduplicated, so a
  later rewrite on the same line would pass under the exemption. Allowing
  only named read-only tails is possible; not worth the parser surface. The
  message names the form instead.
- **Arbitrary commands before the amend** (`python3 tools/check_citations.py
  … &&`). The guard cannot know what they do to HEAD. Today's 05:24 and
  09:18 lines stay asked; the message tells the agent to amend on its own
  line.
- **`git -C <dir> add`** before the amend: asks; agents write `cd` then
  `git add`.
- **A check that the redo of a grant refusal happens only once.** The guard
  judges the redo on its own merits (review recorded, not used in this
  window), which is what the rule needs.

## 6. Estimate

About 30 counted lines (§3.1 about 3, §3.2 about 2, §3.3 about 20, the
after-commit text about 1), measured with `python3 tools/count_loc.py`.
Rule text (`python3 tools/rule_sizes.py`, against master `5d8904ee`): +44
words (`REQUIRED-READING.md` +37, `architect.md` +10, `CLAUDE.md` -3).
