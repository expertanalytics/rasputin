# 2026-10-04: the recovery round, lessons from the morning's handbacks

Recorded by `@orchestrator` on 2026-10-04 from five persona handbacks and one
main-session incident, on a worktree from origin/master d20126b. Quotes from
handbacks are the persona's words, not Ola's. Transcript times are UTC (local
is UTC+2); transcript lines are from
`~/.claude/projects/-Users-skavhaug-projects-rasputin/614a4cd1-2c58-450f-9f74-28e3ab0eaeaf.jsonl`.

Each lesson says whether a rule already covers it, and if not, the smallest
change. Rule changes are proposals for Ola; none is made here.

## 1. A fix that changes where a path resolves breaks older tests too (`@tester`, h9)

What happened: h9's design ruled that `tools/brief.py` resolves `--increment`
against `--worktree`. `@tester` applied the ruled one-line fix to a scratch
copy and ran the whole suite, not just the new red tests. afcceb9 (worktree-h9)
records the result: the shared `trees` fixture built worktrees without the
increment file, so "16 older tests (test_1, test_6, test_8, test_9) went red on
'does not exist'. Three of them (refused, no stderr check) would have passed
for that reason instead of the concurrency rule they name." The fixture was
fixed in the red commit, before `@developer`, who may not edit tests
(`docs/increments/README.md`, step 3), started.

Handback: when a fix changes which directory a path resolves against, apply it
to a scratch copy and run the whole suite.

Covered? In part.
- The refusal tests that check no reason: `tester.md` §4, *Explicit
  Assertions*, already asks for "the exact exception type and error message",
  but it is written for exceptions, not for a CLI refusal (exit status plus
  stderr).
- Running the whole suite against the fix: no rule. It also pulls against
  Ola's standing preference for lean red steps with no throwaway rounds (one
  red step with them took 78 minutes). This time it was cheap: the fix was one line
  the design named (it changes `path = ROOT / args.increment` at `tools/brief.py@afcceb9:241` to `worktree / args.increment`; 4ee0328), and the run was 151 tests.

Proposal (P1, `tester.md` §4, about 40 words): extend *Explicit Assertions*
to refusals ("a refusal test asserts the reason, on stderr or in the message,
not only the exit status"); and add "when the design names the fix as a code
location and it changes a shared default or path, apply it to a scratch copy
and run the whole suite once; an older test it turns red is fixed in the red
commit". This is limited to fixes the design pins, so it does not bring back
general throwaway implementations.

## 2. Merges bring in citations that the branch's own edits make false (`@reviewer`, 23b round 6; `@developer`, 23b)

What happened: 23b merged origin/master twice (e4f07c3, 24161fe). Master's
`27-node-sampling.md` and `25-plain-output.md` cite `refine.hpp` and
`scan.hpp` by line, numbered against master. 23b edits both files, so the
numbers were wrong on the merged branch. `@reviewer`'s round 6 (recorded at
d5aeb82, `docs/increments/23-basin-scale.md`, the round-6 entry) listed
`27-node-sampling.md:115,124,148,150,171,172,263` and `25-plain-output.md:115,120`.
`@developer`'s fix, e3a6add, also found two citations in `15f-edge-strip.md`
(`:1554`, `:1866`) that 23b's own edits had moved, after five review rounds
had approved.

Handbacks: `@reviewer`: "a master merge brings other increments' docs whose
line citations are numbered against master; the branch's own edits make them
false, and no step re-reads them after a merge." `@developer`:
"check_citations.py cannot tell a stale number from a moved one … On a merge
round, re-read every at-risk line whose cited file the branch itself edits."

Covered? The tool covers it. The step that uses the tool does not.
`tools/check_citations.py` computes "files this branch modifies" as
`git diff --name-only master...HEAD` (`changed_files`). After a merge of
master, the merge base is master's tip, so the set is the branch's own edits,
and a merged-in doc's citation into `refine.hpp` is listed as at-risk. (I
read this from the code; I did not run it on the 23b tree, where two reviewers
are working.) But `REQUIRED-READING.md`, *Before you publish*, asks for the
at-risk list to be re-read only "on a prose or tooling branch". `reviewer.md`
§5 check 2 covers "every prose claim the increment touched", and a doc that
came in with a merge was not touched by the increment. Every at-risk entry
also looks the same, whether the cited line is unchanged, moved or rewritten,
so a long list gets skimmed.

Proposals:
- P2a (rule, `reviewer.md` §5 check 2, one sentence): "After any merge into the
  branch, rerun `check_citations.py` and re-read its at-risk list as
  quotations, including docs the merge brought in."
- P2b (tool, `tools/check_citations.py`, about 40 lines plus tests,
  governed path, by day, through the pipeline): for each at-risk citation,
  compare the text at the cited line at the merge base with the text at HEAD.
  Label it `unchanged`, `moved to :N` (the old text is found at a new line) or
  `rewritten`. Only `moved` and `rewritten` need a human. This is
  `@developer`'s "stale versus moved" distinction, made by the tool.

## 3. "Already pushed" is not "CI ran" (`@reviewer`, 23b round 6)

What happened: PR #162 (23b) conflicted with master after #165 merged.
`gh pr checks 162` printed "no checks reported on the
'worktree-agent-a8b29c323dc2ad69d' branch" (transcript line 76, 06:32:40),
which looks the same as a run that has not started. The GCC 13 fix 4cd28dc
was therefore never built by CI. `gh pr view 162 --json mergeStateStatus`
still says `DIRTY` at the time of writing. `.github/workflows/main.yaml` runs
on `pull_request` and on pushes to master only. GitHub runs `pull_request`
workflows on the PR's merge ref, and it cannot make that ref while the PR
conflicts. (Source: GitHub community discussion 26304 and the Classmethod
article, from search results, not read in full; unchecked.)

Handback: "Check mergeStateStatus with gh pr checks."

Covered? No. `reviewer.md` §5, *Precondition: CI status*, treats red CI and
"a workflow that does not exercise the current build" as CHANGES REQUESTED,
but does not say what "no checks reported" means.

Proposals:
- P3a (rule, `reviewer.md` §5 precondition, one sentence): "'No checks
  reported' is not green: read `gh pr view <pr> --json mergeStateStatus`; on
  `DIRTY` CI cannot run until master is merged in."
- P3b (tool, `tools/session_state.py`, about 20 lines, governed): the recap
  lists each open PR with its `mergeStateStatus` and check summary, so a
  conflicting PR shows up on every recap. This could also be part of stage 3
  of the dispatcher plan (the merge runner), which Ola has already approved.

## 4. A "full pytest" that the worktree cannot run (`@developer`, h9)

What happened: the h9 brief asked for green on full `pytest` and forbade C++
builds. The h9 worktree's venv has no extension:
`.claude/worktrees/h9/.venv/bin/python -c 'import tin_engine._core'` raises
`ModuleNotFoundError: No module named 'tin_engine'` (checked now). The green
commit is 4ee0328.

Handback: "briefs should name the narrower check or provide a venv with the
extension."

Covered? No. `CLAUDE.md` §3, *Briefs*, says what each persona's brief
contains, not which checks the worktree can actually run.

Proposal (P4): put it in h9's `tools/brief.py` rather than in prose. When it
assembles a brief for a worktree whose `.venv` cannot import `tin_engine`, it
writes the narrower test command (the suites that do not import the
extension), or says that the extension must be copied in first
(`REQUIRED-READING.md`, *Stale artifacts*). About 10 lines on h9's own
branch. If Ola prefers prose: one sentence under *Briefs*: "a brief that
forbids a C++ build names the test command that can pass without one."

## 5. A merge go-ahead lost with the session (main session)

What happened, from the transcript (it differs in two places from the account
the main session gave Ola, and from the account in this task's brief):

- 06:33:06, line 109: after the restart ("continue, I killed the window by a
  mistake"), the main session ran `gh pr merge` on #167 and #168. It acted
  on the go-ahead recorded in `session.md`. `guard_push.py` asked (line 111,
  `"permissionDecision": "ask"`). The command ran at 06:38:09 (line 113), so
  it was approved at the prompt. **GitHub** refused both: the head branch was
  not up to date, and master's branch protection has
  `required_status_checks.strict: true`.
- 06:38:12 to 06:38:25, lines 123 to 135: the auto-mode classifier then denied
  **the read-only diagnostic** (`gh pr view … mergeStateStatus`, `gh api …
  branches/master/protection`) as "[Merge Without Review]", and after that a
  plain read of the h9 worktree. It did not deny the merge itself.
- 06:44:41, line 203: the main session told Ola the classifier "refused my
  `gh pr merge` as 'merge without review'". That is a cause it had not
  checked (`CLAUDE.md` §3, *Report only finished, verified results*).
- 06:41 to 06:47, lines 157 to 226: Ola ran the scripts himself. `--auto`
  failed because auto-merge is disabled on the repository
  (`allow_auto_merge: false`). Ola then said "merge #167 and #168 when
  green" in the new session, and a background job ran update-branch, waited
  for CI and merged, one PR after the other. #167 merged at 06:55:16.

Why the classifier stopped: Anthropic's engineering note on auto mode says
"The classifier sees only user messages and the agent's tool calls; we strip
out Claude's own messages and tool outputs", and "The prompt establishes what
is authorized; everything the agent chooses on its own is unauthorized until
the user says otherwise." A go-ahead that exists only in `session.md`, which
the session reads as tool output, can never count as user authority. That
holds even when the line is accurate.

Lesson candidate (main session): after a context loss, a recorded publishing
approval is asked again in the new session before anyone acts on it.

Covered? Not quite. `REQUIRED-READING.md`, *Before you publish*: "A grant
covers one occurrence", with the test "have I already performed this named
act once, and has the tree changed since?" A grant from a lost session passes
that test, because the act was never performed. The same section says
`session.md` holds "no rulings, no history", yet the live `session.md` QUEUE
line now reads "push (Ola yes)" twice. That is the same pattern, waiting to
happen again.

Proposals:
- P5a (rule, `REQUIRED-READING.md`, *Before you publish*, one sentence after
  "A grant covers one occurrence"): "A grant does not outlive the session it
  was given in: a yes recorded on disk (`session.md`, a handback) is asked
  again after a context loss."
- P5b (stage 3 of the dispatcher plan, the merge runner): because protection
  is strict and auto-merge is off, the runner does update-branch, waits for CI
  and merges, one PR at a time. It asks for one fresh yes per batch, at the
  start.

## 6. Other things the record shows

- **A review past two rounds.** 23b's code review reached round 6
  (`23-basin-scale.md`, Review section), and `session.md` now names round 7
  (67ca1ac). Rounds 3 to 6 were each started by new work (`@perf`'s
  regression verdict, two master merges, the N18 ruling), not by the same
  finding coming back. The "two-round cap" is quoted as Ola's in the 23a-2
  design entry of the same file, but no rule file states it: `grep -rn
  'two-round\|two round'` over `CLAUDE.md`, `.claude/` and
  `docs/increments/README.md` finds only `orchestrator.md`. Question for Ola,
  below.
- **Worktrees pile up.** `.claude/worktrees/` holds 57 entries. No finding
  yet, but nothing prunes them.

## 7. Rule-file size (Ola, 2026-10-03: every retrospective measures it)

`python3 tools/rule_sizes.py`, against the last retrospective's commit
(2d2d5a6): 9,890 words in all, +117. `REQUIRED-READING.md` grew by 70 words
(1,755; h8, b91f1be, the recap's running jobs and the `session.md` format) and
`perf.md` by 47 (620; increment 24, 23d4dad). No file shrank. Without
`docs/PRINCIPLES.md`, which the 2026-10-03 baseline did not count, the total
is 8,469, against 8,352 then. P1, P2a, P3a and P5a together would add about
90 words.

Proposed cut (C1, about 90 words, `REQUIRED-READING.md`, *The harness*): the
sentence that lists every command `guard_push.py` asks about (`git push`, `gh
pr create/merge/ready/edit`, `gh release`, … `curl` with a writing method)
repeats the hook's own docstring and patterns, and it goes stale whenever the
hook changes. Replace it with "`guard_push.py` asks before any act that
publishes, rewrites history or writes git configuration; its docstring is the
list." The cut is larger than the four proposals add.

## 8. Questions for Ola

- **Q1, the whole-suite run in the red step (P1).** Should `@tester` run the
  whole suite against a scratch copy of the fix, limited to fixes the design
  names by code location? Default: yes, limited that way.
- **Q2, the review-round cap.** Is the two-round cap counted per change set
  (each new piece of work starts at round 1 again) or per increment? If it is
  a rule, it belongs in `reviewer.md` or `docs/increments/README.md`. Default:
  per change set, written as one sentence in `docs/increments/README.md`.
- **Q3, P5a.** Should a recorded publishing yes lapse with the session that
  heard it? Default: yes. Until then, the "(Ola yes)" entries in the live
  `session.md` QUEUE line should be asked again before the push.

**Ruled 2026-10-04 (Ola, 07:48 and 07:50 UTC).** Q1: "Q1, yes": `@tester`
runs the whole suite against a scratch copy of a fix the design names by
code location. Q3: "Q3 yes": P5a as worded. Q2: Ola asked "when did I rule
this? Sounds like something you've come up with, but I might be wrong",
then: "the data is there, if we need to look. Duplicating this simply
doesn't seem important." No cap is written down, and `orchestrator.md`'s
watch-list line "a review that went past two rounds" is dropped, not moved.
No turn of Ola's in this project's transcripts states a two-round cap, so
the "(Ola's two-round cap)" in `docs/increments/23-basin-scale.md` (23a-2
design, round 2) is not his. P1's refusal half, P2a, P2b, P3a, P3b, P4, P5b
and C1 are still unruled.

## Review

**Round 1**, `@reviewer`, on dd48544 (base master d20126b): **CHANGES REQUESTED**. 0 production lines (two prose files). Not pushed, so no CI. `python3 tools/check_citations.py`: all resolve, none at risk. Checked and holding: the cited commits (afcceb9, 4ee0328, e4f07c3, 24161fe, d5aeb82, e3a6add, 4cd28dc), transcript lines 8 to 226 and Ola's two quotes, branch protection `strict` and `allow_auto_merge: false`, #167's merge time, #162 `DIRTY`, the rule-size figures against 2d2d5a6, and the two auto-mode quotes. Blocking: (1) Ola ruled Q1 to Q3 on 2026-10-04, so the rulings are recorded, including that no Ola turn states the two-round cap that `23-basin-scale.md` attributes to him; (2) `:35` showed the line before the fix as the fix.

**Round 2**, `@reviewer`, on 2313a5b (range dd48544..2313a5b): **APPROVED**, provided CI is green after the push. 0 production lines; one commit, two prose files, +18/-2, and nothing else changed. Both blocking items from round 1 are closed. (1) Ola's rulings on Q1 to Q3 are recorded after §8 and in `next.md`. Checked against the tree: `orchestrator.md`:23 still holds the "past two rounds" line the ruling drops, `23-basin-scale.md`:2582 holds the "(Ola's two-round cap)" attribution, and the list of unruled proposals matches the file. (2) `:35` now shows the fix: `tools/brief.py@afcceb9:241` is `path = ROOT / args.increment`, and 4ee0328 changes it to `worktree / args.increment`. `python3 tools/check_citations.py`: all resolve, none in a file this branch edits. Not blocking: the round-1 suggestions (`:177-179`, `:194-195`, `:202`, the URL at `:162`) still stand. The rule edit that drops `orchestrator.md`:23 is still to be made; the retrospective should cite that commit.
