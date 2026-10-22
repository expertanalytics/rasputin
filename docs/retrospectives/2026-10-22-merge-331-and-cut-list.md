# 2026-10-22 post-merge check of #331, and the rule cut list

Times are UTC. Sources: `gh pr view`, `git log`, `.git/harness/reviews.jsonl`,
and the main session's transcript (`3b8a12b2-…jsonl`, "line N"). What
`2026-10-22-merges-329.md` and `2026-10-22-night.md` record is not repeated.

## 0. Duties owed

`python3 tools/duties.py`: "post-merge check (1 merge since the last
retrospective, 2026-10-22 04:59 UTC: #331)". #331 merged 05:43:55; this check
was spawned 05:45:29: **2 min**. #332 and #333 (retrospectives, merged
05:19) and #334 (thread-pool design marked dropped, merged 05:45:38) are
covered in §1 below.

## 1. Post-merge check

**Lead time.** #331, the rule round: `@architect` spawned 04:49:08 (first
mention of `worktree-rules-1022` in the transcript), commit `1d778553` 05:01, review
05:06, push 05:16:57, merged 05:43:55: **55 min**. #334: one design file and
its ROADMAP row, head `0f641c6e` APPROVED in `reviews.jsonl` (05:38:57),
merged 05:45:38.

Findings:

1. **The main session edited a persona's rule file twice, and the recording
   commit carried a content change.** After the review (05:06:10) the main
   session put back `tester.md`'s sentence on the size ceiling (Edit, line
   26503) and committed it as `5ff8c3eb`, subject tagged `(@reviewer)`;
   `docs/increments/README.md` step 4 lets the recording commit carry only
   the reviewer's wordings, and `CLAUDE.md` §3 *Prose edits while writers are
   busy* lets the main session only apply them. It told Ola "a cut I didn't
   ask for" (05:06:43) without checking the cut against `CLAUDE.md` §2's
   "This is the only statement of the rule" — a cause not checked. Ola:
   "I don't like repeats, so why keep it in 2? Don't you obey CLAUDE.md
   §2?" (05:13:02). The main session then cut it again, `f0aea3bd`, with no
   persona tag. No harm reached the remote: `git diff 1d778553 f0aea3bd` is
   empty, so the pushed head is the reviewed tree.
2. **The rule round grew the rule text by 199 words and offset about 20.**
   `1d778553`'s message: "11,714 at base, 11,913 after (+199)". Its message
   also keeps a working note ("11,951 before edits on this base? No: …").
   Ola asked at 05:16:27 "But when did we cut last?"; part B answers.
3. **#331's rules are not in force for the main session.** The main checkout
   is at `2818da34`, 9 commits behind `origin/master` (`git rev-list --count
   HEAD..origin/master`), so `CLAUDE.md`'s *Worktrees* bullet and the
   unattended-mode default rule are not on disk there. Beyond that, the
   `CLAUDE.md` the session holds in context is older than its own disk: the
   copy passed to this spawn still lists the personas in §1 and keeps the
   *Prohibited dependencies* bullet, both cut on disk (`e074eff5`, #331).
   A fast-forward does not reload `CLAUDE.md`; only a restart does. The
   restart rule in `.claude/REQUIRED-READING.md` names the changed persona,
   not the main session. Proposal 1 of `2026-10-22-merges-329.md` (a duty
   for a stale main checkout) would show the first half; see proposal 2.
4. **A verdict given in a handback but not recorded.** One `@reviewer`
   judged two branches (line 26188): the hook recorded only the final
   `Review:` line (`f5c0889b`); the verdict on `1730e0fd` (#333) exists only
   as handback text, and `reviews.jsonl` holds no line for it. A push guard
   that reads the record would call #333 unreviewed. One reviewer, one
   branch.

#332's head `3575b3b8` is the recording commit after `f5c0889b`'s approval;
#334 in order.

## 2. Lessons, recorded

1. **(This check.) A cut of a restatement is in scope of any rule round.**
   Where a rule has one statement, a repeat elsewhere is a defect, not a
   change Ola must approve (finding 1).
2. **(This check.) The main session's `CLAUDE.md` is the one it started
   with.** A merge reaches it only on restart (finding 3).
3. **(This check.) One reviewer spawn judges one branch**, so each verdict
   is recorded (finding 4).

## 3. Proposals

1. **A rule change pays for itself** (the one-line rule Ola approved).
   Text: "* **Rule text pays for itself:** a change that adds words to the
   files `tools/rule_sizes.py` counts cuts as many in the same change,
   unless Ola explicitly allows the growth." (29 words.) Where: `CLAUDE.md`
   §2, after *Strict Size Limit*, as its one statement; `@reviewer` checks
   it under the existing prose-branch scope in `.claude/REQUIRED-READING.md`
   without a restatement. Its own 29 words are paid for by item 3.19 below.
   Cost: one rule line; a later `--base <rev>` option to `rule_sizes.py`
   would let the reviewer measure a branch against its base (Ideas).
2. **Restart after a rule merge reaches the main session too.** The
   restart bullet in `.claude/REQUIRED-READING.md` reads "After merging a
   change to `CLAUDE.md` or `.claude/agents/`, fast-forward the main
   checkout and restart before the next spawn" (about +5 words, paid by
   item 4.3). Evidence: finding 3.

## 4. The cut list (part B)

**Measured.** `python3 tools/rule_sizes.py` at `8d30ce71`: **11,914
words** (`CLAUDE.md` 1,255; `.claude/REQUIRED-READING.md` 2,014;
`docs/increments/README.md` 1,410; `docs/PRINCIPLES.md` 1,396;
`.claude/briefs/common.md` 236; agents: architect 541, developer 241,
orchestrator 482, perf 774, reviewer 604, tester 960; skills:
computational-geometry 545, geospatial-data-formats 392, modern-cxx 334,
python-development 730). The brief's 11,930 and its 10,935 for 2026-10-12
are not re-measured here.

Each saving is words removed minus words of the new text, counted with
`str.split()` over the exact passage (script kept in the session scratchpad,
not committed). "Ruled" marks an item that touches a rule Ola ruled on.

### 4.1 Repeats of a rule stated elsewhere

| # | File, passage | Stated also in | Change | Words |
|---|---|---|---|---|
| 1.1 | `docs/PRINCIPLES.md` A1-A5, B1, B2, all of C, D1, D4, E1, E3, E4 | `REQUIRED-READING.md` *Claims*, *Stale artifacts*, approval boundary; README *The loop*; `CLAUDE.md` §2-§3 | Cut; keep preface, A6, B3, B4, D2, D3. D1 "Dispatch serially" contradicts the two-writer limit. Proposed 2026-10-21-research §5, unruled. **Ruled**: D1 against Ola's two-writer rule; E1 is Ola's approval boundary, kept in `REQUIRED-READING.md` | −914 |
| 1.2 | `CLAUDE.md` §3, numbered steps 1-4 | Same section, *Step order*; README *The loop* | Replace with "Every code alteration follows the step order below, dispatched by the main session." | −44 |
| 1.3 | `REQUIRED-READING.md` "Also read `docs/increments/README.md` (…)" | `CLAUDE.md` §3 first paragraph | "Read the increment file for whatever you are working on." | −10 |
| 1.4 | README *Never dropped* "`@reviewer` before every push, Ola's yes for each push, and the prompt on a governed file" | `REQUIRED-READING.md` approval boundary and harness | "…and the approval boundary (`.claude/REQUIRED-READING.md`)". **Ruled**: push yes, kept by pointer | −11 |
| 1.5 | `developer.md` "`@reviewer` reconciles the two (`reviewer.md` §5, check 3)." | `reviewer.md` check 3 | Cut | −8 |
| 1.6 | `reviewer.md` §1 "; reconcile the count against the estimate (check 3 below)" | Same file, check 3 | Cut | −9 |
| 1.7-1.11 | Skills: computational-geometry *Dependencies* and *Change Limit*; geospatial *Prohibited* and *Change Limit*; python *Change Limit* | `CLAUDE.md` §2, which every session loads | Cut the five pointers. **Ruled** (§2 single statement, 2026-10-22): cutting a pointer keeps it | −49 |
| 1.12 | `architect.md` §3, item 3, "(Enforce the prohibited dependencies of `CLAUDE.md` §2 fiercely.)" | `CLAUDE.md` §2, gate | "Does a change add a heavy dependency?" | −7 |
| 1.13 | computational-geometry §3 *Test in the producer's relation* (15 lines) | `REQUIRED-READING.md` *Claims* object-identity question; says itself "not a rule of its own" | "**Oracles use the producer's predicate, never its records** (the general rule: `.claude/REQUIRED-READING.md`, *Claims*)." Also removes the "~93 %" anecdote and worked example | −170 |
| 1.14 | computational-geometry *Degeneracy Handling* | `tester.md` §3A | "handle the cases of `tester.md` §3A, and acute angles." | −9 |
| 1.15 | python-development §4 *Testability* | `tester.md` §2 tiers | Cut | −77 |
| 1.16 | `tester.md` "under the handback heading "Pinned or assumed beyond the design"" | `.claude/briefs/common.md` handback headings | "in the handback" | −7 |
| 1.17 | `perf.md` §1 "the 1 m benchmark and a thread-scaling sweep, compared with…" | README *Acceptance* (its one statement) | Cut, keep the pointer | −14 |
| 1.18 | `perf.md` §2 *Record the power state…* | README *Acceptance* | Cut | −30 |
| 1.19 | `architect.md` *Ideas* "read the private ideas file; … hand back each idea it leaves out" | `REQUIRED-READING.md` ideas paragraph | "a design names each ideas-file entry it takes up or overlaps." | −14 |
| 1.20 | `CLAUDE.md` "(recap, every prompt, Ola's status line)" | `REQUIRED-READING.md` *The harness* | Cut (carried from 2026-10-21-research) | −6 |
| 1.21 | `perf.md` §4 "No figure from a run still in progress, and no cause…" | `CLAUDE.md` §3 *Report only finished*; `REQUIRED-READING.md` *Claims* | "An unprofiled cause is "thought to be"." | −20 |

Subtotal **−1,399**.

### 4.2 Rules with no check, hook, test or incident in two weeks

Searched: every `docs/retrospectives/` file dated 2026-10-08 to 2026-10-22
(35), every `docs/increments/` file changed since 2026-10-08 (74), and the
tracked code (`src_python`, `include`, `src`, `tests`, `tools`,
`.claude/hooks`, `CMakeLists.txt`, `pyproject.toml`, `.github`: 556 files),
with one regular expression per rule, counting files that match. A rule
found only in its own file counts as unused.

| # | Passage | Search result | Change | Words |
|---|---|---|---|---|
| 2.1 | README *Template cross products* and *Match the model to the work* | 0 retrospectives; no agent runs a smaller model | Cut both | −96 |
| 2.2 | python-development *Settings Management* (`pydantic-settings`) | 0 in code | Cut | −10 |
| 2.3 | modern-cxx *Asynchrony & Streaming* (coroutines) | 0 `co_yield`/`std::generator` in code | Cut | −26 |
| 2.4 | geospatial §4 *Output Architecture* (three metadata tiers) | 0 anywhere | Cut | −81 |
| 2.5 | geospatial *Vector* "`geojson` or `ujson`" | 0 in code | "`shapely`." | −12 |
| 2.6 | geospatial *GeoJSON* streaming with `ijson` | 0 in code | Cut | −24 |
| 2.7 | `tester.md` *System & Integration Tests* ("multi-tiered metadata 2D surface") | 0 | "end to end, from input files through the C++ core to the written mesh." | −13 |
| 2.8 | `PRINCIPLES.md` preface, *Last exercised* and the three-increment retirement rule | every entry reads 5b, 5c or "not yet"; never run | Cut | −30 |

Subtotal **−292**. The skill items (2.2-2.6) come from the template
committed 2026-09-16 (`b92cdd6b`, Claude co-authored); Ola may have written
that template, so they are his to keep. Also unused and kept, no saving
claimed: "a subagent whose product is a file creates that file first"
(0 matches; preventive), reviewer check 6's ruling scope (merged today).

### 4.3 Shorter wording, same meaning

| # | Passage | Shorter text | Words |
|---|---|---|---|
| 3.1 | README step 1 *Legacy* | "*Legacy.* What the `legacy-archive` tag holds on the subject, with the command and the files it returned (`git grep -n <pattern> legacy-archive -- legacy`), and what is carried across or why nothing is. Where something is, `@architect` reads that source and states its intent before `@tester` starts." | −100 |
| 3.2 | README "Steps 2 and 3 are not strictly once each…" | "A ruling or a platform-retuned constant after the red commit is `@tester`'s amendment, in its own commit with the reason. `@developer` never edits tests, and no test change hides in an implementation commit." | −52 |
| 3.3 | README "These exist because a round that costs three times…" | Cut | −21 |
| 3.4 | README mutation paragraph, sentences 2-4 | Cut; the rule sentence stays. **Ruled**: mutation rounds kept (token-budget cuts) | −45 |
| 3.5 | README "A "documentation debt" section is a ledger…" | "No debt ledger." | −32 |
| 3.6 | `REQUIRED-READING.md` step 1's list of the recap's contents | "the round recap, `.claude/current-task/` and the predecessor's last human turns," **Ruled**: recap each round, kept | −14 |
| 3.7 | `tester.md` opening paragraph | "You write the failing suites and hold the coverage floor `testing.md` defines." | −38 |
| 3.8 | `architect.md` opening paragraph | "You rule on structure, boundaries and testability, before any production file exists." | −25 |
| 3.9 | `reviewer.md` "You are the Senior Code Reviewer…" | "Your mandate is the checks no gate makes:" | −13 |
| 3.10, 3.11 | `reviewer.md` "(There are no sections 2-4…)", `tester.md` "(There is no section B…)" | Cut | −23 |
| 3.12 | `perf.md` plot bullet | "A plot of observed against modelled series has three lines unless Ola asks for more: observed, the best so far and the run under investigation. Plain-word legends mark the best and give each run a reasonable explanation; when too long, they go in a caption at the base, like a figure legend." **Ruled** (2026-10-20, Ola's words only): keeps his phrases; skip if he prefers the text as is | −29 |
| 3.13 | `tester.md` *Coverage Floor* | "keep project-wide line coverage at or above `testing.md`'s floor (`--cov-fail-under`, Python only); a per-module floor is checked in review, and anything lower is justified in the PR. Every documented invariant and every specially-handled condition (NaN, NoData, empty or single-element input, boundary intersection) has a test that names it." | −31 |
| 3.14 | `tester.md` *Determinism*, *Numeric bounds*, *Zero Flakiness* | "**Determinism:** no races, no flaky test (a blocking bug); bitwise results or a stated tolerance. Every absolute bound states the scale it assumes and the largest input it was checked at." | −47 |
| 3.15 | `perf.md` *Sanitize a scratch build* | "a C++ patch for a timed or full-size run first runs once on a small case under `-fsanitize=address,undefined` (a C++ driver or test binary), or `-fsanitize=undefined` with libc++ extensive hardening through `_core`, where ASan cannot load. A crashed run has no numbers; a deliberate crash stays out of a script's default loop." | −49 |
| 3.16 | `developer.md` fp-contract explanation after the colon | Cut | −22 |
| 3.17 | `CLAUDE.md` *Milestone updates* | "After each milestone, one line to Ola on where the round stands." | −10 |
| 3.18 | `REQUIRED-READING.md` *The harness*, the two guards' first sentence | "`guard_push.py` asks before any act that publishes, writes to the forge, rewrites history or writes refs, remotes or git config, or uses `--no-verify` (its code lists them); `guard_governance.py` asks before any write to a file that states rules or to a `tools/` file named after a standard-library module (its code lists what it resolves), except an absolute scratchpad path written by Edit or Write or one plain shell command." (2026-10-21-research P5) | −70 |
| 3.19 | `CLAUDE.md` §4 *CI* "it is meant to be the one required check … never report." | Cut to "`CI result` is green only when every gating job passed, or was skipped on a prose-only change." The rest is stale: `gh api …/required_status_checks` returns `["CI result"]` today | −37 |

Subtotal **−658**.

### 4.4 History, dates, quotes and anecdotes in rule files

| # | Passage | Change | Words |
|---|---|---|---|
| 4.1 | README intro: the audit sentence and "`kernel-sufficiency-audit.md` is the first" | "An audit of shipped work may live here; its status line says so." | −35 |
| 4.2 | README "(for example an exact sup-norm guarantee on a DEM with constraints)" | Cut | −11 |
| 4.3 | `REQUIRED-READING.md` ""Ok, let's continue" is not an instruction. Neither is `git log`." | Cut (a quote) | −11 |
| 4.4 | `REQUIRED-READING.md` "A question costs a line; a guessed round costs ~150k tokens." | Cut | −11 |
| 4.5 | `REQUIRED-READING.md` the `CMakeLists.txt` planting example | "Verify a gate by planting what it forbids." | −17 |
| 4.6 | `REQUIRED-READING.md` "(no "current anchor" hash, no count of rounds or commits)" | Cut | −10 |
| 4.7 | `perf.md` "a crash in a Python process shows as a crash dialog on the owner's screen" | Cut | −15 |
| 4.8 | `PRINCIPLES.md` A6 Origin ", which measures 48 of 48 governance commits as agent-written" | Cut | −10 |
| 4.9 | `PRINCIPLES.md` "Last exercised: 5c." in the four kept entries | Cut | −12 |

Subtotal **−132**. The *Origin* commit pointers in the kept `PRINCIPLES.md`
entries are references by the file's own design, not anecdotes. The
example in `CLAUDE.md` §3 *Finish before starting* ("this makes 4 of 3") is the wording to use, not an anecdote.

### 4.5 Total

All items: **−2,481**, from 11,914 to about **9,433**. Under 11,000 needs
915: 1.1 alone gives exactly 11,000; 1.1 with 1.2 and 3.19 gives **10,919**.
Proposals 1 and 2 above add about 34 words, paid by 3.19 and 4.3.
Owner of every cut: `@architect`, after Ola sees this list.
