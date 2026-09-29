# Retrospective: orchestrator boundary & hook gates (2026-09-29)

Status: audit of work already shipped, dated to the commits that answered it.
Scope: increment 16e (PR #112, merged 2026-09-29T16:39Z, merge commit
`00c2c88`), the auto-catchment pair (PR #110 `a7e51d1`, PR #111 `7d9882d`,
both merged 2026-09-29), and the session workflow up to 2026-09-29. Requested
by Ola with role bleed as the named focus (`a612269`, agenda in
`docs/retrospectives/next.md`).

Every incident below cites the command, file or commit a sceptic would run to
check it. Where a claim could not be verified from disk, it says so rather than
asserting.

---

## 1. Orchestrator boundary & role bleed

### 1.1 Tool scoping is uniform: only `@reviewer` is read-only

`sed -n '1,12p' .claude/agents/<persona>.md` gives the `tools:` line for each
persona:

| persona | tools |
|---|---|
| orchestrator | Read, Grep, Glob, Bash, **Write, Edit**, Skill |
| architect | Read, Grep, Glob, Bash, **Write, Edit**, Skill |
| tester | Read, Grep, Glob, Bash, **Write, Edit**, Skill |
| developer | Read, Grep, Glob, Bash, **Write, Edit**, Skill |
| reviewer | Read, Grep, Glob, Bash, Skill (no Write/Edit — read-only by design) |
| perf | Read, Grep, Glob, Bash, **Write, Edit**, Skill |

Two structural facts fall out of this table:

1. **`@orchestrator` can author production and test files.** Its remit is
   coordination (`.claude/agents/orchestrator.md` description: "breaks them
   into sub-tasks, and drives the strict TDD loop"), which never needs Write or
   Edit on code. The grant is unused surface that only enables role bleed.
2. **Claude Code agent frontmatter scopes tool *names*, not *paths*.** So
   `@tester` and `@developer` each hold `Write, Edit` over the entire tree.
   Nothing in the frontmatter can express "`@tester` may not write
   `src_python/`" or "`@developer` may not edit `tests/`". That boundary can
   only be enforced by a `PreToolUse` hook — and no such hook exists
   (see 1.3 and §2).

### 1.2 Verified role-bleed incidents

These are drawn from the agenda in `docs/retrospectives/next.md` (built by Ola
and prior sessions) and re-checked against the tree. They are real and
disclosed; the pattern, not any single fix, is the finding.

- **Main session about to author a persona's deliverable (2026-09-29).** Two
  consecutive human turns, recovered by `tools/session_state.py` under "last
  human turns":
  `[2026-09-29T15:33:14] "Why would you write the @orchestator reflection doc?
  You're not the @.claude/agents/orchestrator.md"` and
  `[2026-09-29T15:34:35] "...You're not the @orchestrator"`. The main session
  is not the `@orchestrator` persona (MEMORY: `not-the-orchestrator-persona`);
  it was about to write an orchestrator/reflection document itself. Caught by
  Ola's eye, not by any mechanism — the file would have been a `docs/` write,
  which is ordinary tree work and trips no guard.
- **Redirect to the owning persona worked (2026-09-29).** Turn
  `[2026-09-29T15:48:53] "R1 left-align is fine, have @tester fix the three"`
  was honoured: commit `170bb2d` ("16e red amendment: three tests corrected to
  Ola's R1 left-align ruling") touches **only** `tests/python/`
  (`git show 170bb2d --stat`). The fix went to `@tester`, not the main session.
  This is the boundary working — but it worked because Ola named the persona.
- **`@tester` writing production, repeatedly** (from `next.md`, agenda §
  "Agents taking on each other's work"): ~560 lines of throwaway `src_python/`
  modules in the 16b-1/2 red step; a temporary edit-and-`git checkout` of
  `src_python/tin_engine/io/geopackage.py` on 2026-09-29. Brief-level prose did
  not stop it; `next.md` itself concludes "a hard limit would have".
- **`@architect` doing implementation in the design step (2026-09-29 night,
  increment 22):** an uncommitted Python prototype to check catchment area
  against NVE (304.91 vs 305.54 km²). Useful evidence, wrong step.
- **`@developer` writing behaviour before the test (2026-09-29 night,
  increment 22):** changed the window rule during the PR 1 green step, with the
  design amendment and pinning test landing *after* the code; and added
  user-visible lake-refusal reporting under `--lakes` in a review-fix commit
  with no failing test first (`@reviewer` caught it, `@tester` pinned it after).

### 1.3 The 16e trace itself is clean

Against the above, increment 16e's commit trace holds the red-before-green
discipline exactly (`git show --stat` per commit):

- `3599be1` design — only `docs/increments/16e-multi-features.md`
- `cb6f78b`, `170bb2d`, `3e2e084` red + amendments — only `tests/python/`
- `e3ec289` green — only `src_python/` (`cli.py`, `feature_input.py`); touches
  **zero** test files (`git show e3ec289 --stat | grep -c test` = `0`)
- `524b27e` — only `ROADMAP.md`

So at the *commit* granularity 16e is disciplined. The role bleed that occurred
this cycle was at the *session* granularity (1.2, first bullet) and left no
commit, which is exactly why it is hard to gate.

### 1.4 Verdict on restricting `@orchestrator`'s tools

Removing `Write, Edit` from `@orchestrator` is correct and costless — the
persona never needs them — but it does **not** solve the observed problem,
because:

- the 2026-09-29 near-miss was the **main session**, which does not run under
  `orchestrator.md`'s frontmatter and keeps the full default toolset; and
- the damaging bleed is `@tester`→`src_python/` and `@developer`→`tests/`,
  which frontmatter cannot express (path, not tool name).

The enforceable fix is a per-persona path guard as a `PreToolUse` hook (§2, §4).
Do the tool-list trim anyway as defence in depth.

---

## 2. Hook & gate blindspots

### 2.1 The three governance hooks are drafted but INACTIVE

This is the headline gate finding. The hooks exist in `.claude/hooks/`
(`guard_governance.py`, `guard_push.py`, `gates_after_commit.py`) but are wired
into **no active settings file**. Evidence:

- `ls .claude/settings*.json` returns only `settings.local.json`, which has a
  `permissions` block and **no `hooks` key**.
- The only file that wires the three hooks is
  `.claude/settings.json.pending-orchestrator-review`, added by commit
  `aeccda9` whose own subject reads **"Draft three Claude Code hooks --
  UNAPPROVED, INACTIVE, UNTESTED"**.
- The user-level `~/.claude/settings.json` *does* register `PreToolUse` and
  `PostToolUse`, but both point only at `~/.config/iterm2/cc-status` (a status
  line), not at any governance hook.

Consequence: for the entire 16e / auto-catchment period, `guard_governance.py`,
`guard_push.py` and `gates_after_commit.py` **never executed**. Every "the hook
would have caught this" in `next.md` describes a guard that was not running.
`.claude/REQUIRED-READING.md` (§ "The harness") is consistent with this — it
says these hooks must be *proposed* for Ola's approval and never self-added —
so the inertness is by policy, not by accident. But the policy leaves the
boundary entirely on the agent's habit, which is the condition every guard here
was written to remove.

### 2.2 `guard_governance.py` WRITES regex is bypassable by construction

Even when activated, the Bash arm is a substring tripwire, not a parser. The
file says so itself (`.claude/hooks/guard_governance.py:24-25`): "a path built
by variable expansion escapes it. It is a tripwire, not a sandbox", and
lines 132-134 repeat it in the user-facing reason. Concretely bypassable:

- **Variable-expanded path.** `WRITES` (`.hooks/guard_governance.py:56-59`)
  matches redirection/`tee`/`cp`/`sed -i`/`python3 -`, and `governed()` matches
  the path as a literal substring (`name.rstrip("/") in command`, line 104). A
  command like `f=CLAUDE.md; echo x >> "$f"` contains no literal `CLAUDE.md`
  token adjacent to the write and the path is only in `$f`, so the substring
  scan at line 104 misses it.
- **Subshell / indirect write.** `printf ... | tee "$(cfg_path)"` or a heredoc
  into a computed path likewise carries no literal governed token.
- **Write shapes not in the regex.** The regex enumerates specific verbs;
  anything writing by another means (e.g. a Python one-liner that is invoked as
  `python3 script.py` rather than `python3 -`) is outside it.

This is a known and accepted limitation (the hook is `ask`, not `deny`, and is
declared a tripwire). It is worth recording that the tripwire is only ever as
good as the literal path appearing in the command text.

### 2.3 `guard_push.py` is sound in intent but likewise inactive and text-only

`guard_push.py:20-30` matches `git push`, `gh pr create|merge|ready|edit`, `gh
release`, `--no-verify`, history rewrites and force pushes, returning `ask`.
The patterns are reasonable. But (a) it is inactive per 2.1, so it did not
guard any of the #110/#111/#112 pushes, and (b) like 2.2 it reads command text,
so a push issued through a wrapper or a variable (`$GIT push`) would not match
`\bgit\b`. The permission system is explicitly *not* the backstop
(`REQUIRED-READING.md` § "The harness": "auto mode has let unapproved pushes
through"), so with the hook inactive the push boundary rested entirely on the
agent asking — which is what §"Before you publish" of REQUIRED-READING already
mandates in prose.

### 2.4 `gates_after_commit.py` exit-2 path has no execution history

The task asks whether exit code 2 reliably alerted the agent loop to red gates
without human intervention. It cannot have, because the hook never ran (2.1):
it is registered only in the pending settings file. The code is correct in
design — it runs the four `tools/check_*` gates plus `ruff check`/`format`
after any `git commit`/`git merge`, prints their unsuppressed output to stderr,
and `return 2` "feeds stderr back to the agent"
(`gates_after_commit.py:72-86`). Its own docstring flags it as NOT
AUTHORITATIVE (mypy and ctest need a build; CI decides). The design is sound;
the *deployment* is the gap. There is no evidence — because there can be none —
that exit 2 ever surfaced a red gate to the loop.

---

## 3. Subagent skill & build-context leaks

### 3.1 Skill invocation cannot be verified from disk, and that is itself the risk

`.claude/REQUIRED-READING.md:49-54` requires each persona to invoke the
relevant `Skill` (modern-cxx / computational-geometry / python-development /
geospatial-data-formats) itself, because "the `skills:` frontmatter key does
not reliably preload them, and subagents do not inherit skills from the
caller." Whether a given subagent actually called `Skill` leaves no artifact in
the git tree or in `docs/`, so this audit cannot confirm compliance per
increment. The structural point stands: nothing on disk records a skipped
skill, so a persona that reasons from inherited caller context instead of its
own skill load fails silently. The only durable countermeasure is the brief
naming the skill explicitly and the persona reporting it invoked — neither of
which is currently machine-checked.

### 3.2 Stale `_core` risk was low this cycle; 16e touched no C++

`REQUIRED-READING.md:82-98` codifies the stale-extension trap: `pytest` does
not rebuild `_core`, and the rebuild+`cp`+`touch` steps are manual.

- **16e touched no C++** — `git show e3ec289 --stat` lists only
  `src_python/tin_engine/cli.py` and `feature_input.py`. No `_core` rebuild was
  needed, so 16e's Python tests could not have measured a stale extension.
- **Increment 22 did touch C++** (`include/terrain/hydrology/upstream.hpp`,
  `include/terrain/vector_simplify/area_collapse.hpp`, `bindings/core.cpp`).
  Its Bygdin acceptance run is recorded under
  `docs/benchmarks/2026-09-29/bygdin/` (README, `run.sh`, `logs/`, reduced
  GeoJSON present), which implies a build occurred. Whether the `cp` of the
  fresh `.so` into the venv preceded every test run is not recoverable from
  disk. No stale-`_core` incident is provable for this period; flagged as an
  unverifiable-from-disk gap, not a clean bill.

---

## 4. Process & governance compliance

### 4.1 LOC reconciliation

**Increment 16e** (`docs/increments/16e-multi-features.md:418-434`): estimate
~45 production lines, worst case ~70. The green commit `e3ec289` shows
`115 insertions(+), 61 deletions(-)` across `cli.py` and `feature_input.py`
(`git show e3ec289 --stat`), but that gross count includes the reformatting
churn of the `_feature_source`→`_feature_sources` rename and blank/comment
lines that CLAUDE.md §2 excludes. Net change is well within the worst-case ~70
and nowhere near the 700 ceiling. **No split seam was needed or missed.**

**Increment 22** (`docs/increments/22-auto-catchment.md:685-712`): the doc
carries a full LOC table with `@reviewer`'s as-built counts (PR 1 = 619 net
over `master..608e366`; PR 2 = 406 net). The per-file estimates were exceeded
in places — `catchment.py` estimated 170, built 244 in PR 1 (a ~44% overrun) —
but the increment was **pre-split into two PRs at defined seams** ("the fine
catchment" ~510 lines / "the reduction" ~350 lines,
`docs/increments/22-auto-catchment.md:715-725`), and both PRs landed under 700.
**No unchecked overrun; the split seam did its job.** This is a model of the
rule working: the estimate table names the seam before the code, and the
reviewer reconciles actuals against it in the same PR.

### 4.2 Was `@reviewer` run before every push?

- **Increment 22: yes, and it is verifiable.** The LOC table cites "`@reviewer`'s
  count at review (2026-09-29)" with the exact commit ranges it counted
  (`docs/increments/22-auto-catchment.md:707-712`). The review left a durable
  trace because its output was written into the increment doc.
- **Increment 16e: not verifiable from disk.** `@reviewer` is read-only
  (§1.1) and leaves no commit; `git log 3599be1~1..524b27e | grep -i review`
  returns nothing. The session's own `session.md` names `@reviewer` as a loop
  step, but there is no on-disk artifact proving it ran on 16e. This is a
  structural blind spot: unlike the red and green commits, the reviewer step
  has no mandated trace, so "was the reviewer run?" is answerable only when the
  reviewer happens to write into a tracked file (as in 22).
- **Documentation/tooling pushes:** `REQUIRED-READING.md:118-126` already
  requires `@reviewer` once on any branch before its first push "whether or not
  the branch contains production code", with a prose/tooling scope. Compliance
  on non-code branches is likewise unverifiable from disk for the same
  no-trace reason.

---

## 5. Rule updates (drafts only — not applied; do not push)

Per the task and `REQUIRED-READING.md:100-126`, these are drafted here for
Ola's decision. None is written into a governance file by this audit.

### R-A. Trim `@orchestrator` (and design/review-adjacent) tool grants

In `.claude/agents/orchestrator.md` frontmatter, change
`tools: Read, Grep, Glob, Bash, Write, Edit, Skill` to
`tools: Read, Grep, Glob, Bash, Skill`. Rationale: the orchestrator coordinates
and never authors code; removing Write/Edit is costless defence in depth
(§1.1). Machine-check: `grep '^tools:' .claude/agents/orchestrator.md` must not
contain `Write` or `Edit`. Caveat recorded in §1.4: this does not constrain the
main session, which is the persona that nearly bled on 2026-09-29.

### R-B. A per-persona path guard (the real fix), as a proposed `PreToolUse` hook

Frontmatter cannot scope paths (§1.1). Propose to Ola a new
`.claude/hooks/guard_persona_paths.py` (`ask`, never `deny`) that, on
`Write`/`Edit`/`NotebookEdit`, refuses:

- `@tester` writing under `src_python/` or `include/` or `bindings/` or `src/`;
- `@developer` (and `@perf`) editing under `tests/`;
- anyone but `@architect`/`@migration-expert` editing `docs/increments/**`
  (`@tester` and `@reviewer` amendments are the known exceptions — encode them
  explicitly).

Blocker to record honestly: it must be confirmed that the `PreToolUse` hook
payload carries the active subagent's identity; if Claude Code does not expose
which persona issued the tool call, this hook cannot distinguish them and the
boundary stays brief-level only. This uncertainty is why the rule is a *draft*
for Ola, not an applied change. Machine-check once built: plant a `@tester`
write to `src_python/x.py` and confirm the hook returns `ask` naming it, per
`REQUIRED-READING.md:68-70`'s "plant what it forbids" discipline.

### R-C. Give the `@reviewer` step a mandated on-disk trace

Add to `docs/increments/README.md` step 4 (or `REQUIRED-READING.md`'s
assessment section): the reviewer's verdict and LOC reconciliation are written
into the increment doc (as increment 22 already does,
`22-auto-catchment.md:707-712`), so "was the reviewer run?" is answerable after
the fact. Machine-check: an increment PR that touches `src_python/` or
`include/` must add or update a reviewer line in its `docs/increments/NN-*.md`.
Rationale: §4.2 — the reviewer is the only loop step with no required trace.

### R-D. Activate the drafted hooks, or record why not

`REQUIRED-READING.md` § "The harness" forbids self-adding hooks, correctly.
But the three hooks have sat inactive since `aeccda9` while `next.md` reasons as
if they were live. Draft for Ola: either approve wiring
`.claude/settings.json.pending-orchestrator-review` into
`.claude/settings.json` (a fresh-yes act per `REQUIRED-READING.md:110-112`), or
add a line to `next.md`/the harness section stating the hooks are deliberately
dormant so no future session cites them as protection. Do not wire them without
Ola's explicit yes.

---

## 6. Summary of verdicts

- **Role bleed is real and recurrent** (§1.2), and the enforceable boundary
  (per-persona paths) exists in neither frontmatter nor any active hook (§1.1,
  §2.1). The 16e *commit* trace is clean (§1.3); the bleed lives at session
  granularity and leaves no commit.
- **All three governance hooks are inactive** (§2.1). `gates_after_commit.py`'s
  exit-2 alerting has no execution history because it never ran (§2.4). The
  `WRITES` tripwire is bypassable by variable-expanded/subshell paths by its own
  admission (§2.2).
- **LOC governance held** (§4.1): 16e within worst case, 22 pre-split at named
  seams and reconciled by the reviewer; no unchecked overrun.
- **`@reviewer` provably ran on 22, unverifiable on 16e** (§4.2) — the reviewer
  step needs a mandated trace (draft R-C).
- **Skill invocation and stale-`_core` are unverifiable from disk** (§3) — no
  incident proven, no clean bill either.

Nothing in this audit was pushed. Rule changes in §5 are drafts for Ola.

---

## 7. Ola's rulings (2026-09-29) and the R-B analysis

**Applied on branch `worktree-harness-hooks-ra-rc`:**

- **R-A:** `@orchestrator`'s `tools:` line is now `Read, Grep, Glob, Bash,
  Skill`. Check: `grep '^tools:' .claude/agents/orchestrator.md`.
- **R-C:** `docs/increments/README.md` step 4 now requires a `## Review`
  section in the increment file, copied from `@reviewer`'s handback by its
  spawner and committed before the push.
- **R-D:** the three hooks are active. `.claude/settings.json.pending-orchestrator-review`
  was renamed to `.claude/settings.json`. Tested by feeding each hook planted
  input: the guards ask on a push, on an Edit of `CLAUDE.md` and on
  `sed -i .claude/settings.json`, and stay silent on `cat CLAUDE.md` and an
  Edit of `src_python/`. Activation found one defect: `gates_after_commit.py`
  looked for ruff only in its own tree's `.venv`, and a worktree has none, so
  every worktree commit came back red with "DID NOT RUN". It now falls back to
  the main checkout's `.venv`, then to PATH. A planted `import os` was
  reported as F401 with exit 2.
- **The data folders are not a channel** (Ola: "What I can't tolerate is two
  agents trying to communicate through the rasputin_folder or the temp
  folder"). This is now a rule in `.claude/REQUIRED-READING.md`. Write access to
  `../rasputin_data` and `../rasputin_scratch` is granted in Ola's user
  settings, not in the repo.

### R-B: a per-persona path guard is feasible

**The blocker in §5 is resolved in principle.** The Claude Code hooks
reference (code.claude.com/docs/en/hooks, "common input fields") says a hook
that fires inside a subagent receives `agent_type` (the agent's name, e.g.
`tester`) and `agent_id`. Both are absent in the main session. `PreToolUse`
inside a subagent fires the hooks from settings files. This has not been
observed live yet: step 1 of the build is a hook that logs `agent_type` in a
spawned `@tester`.

**Where it lives.** It should be one `PreToolUse` hook in
`.claude/settings.json` that looks up `agent_type` in a single table. The
alternative is hooks in each persona's frontmatter, which also works
according to the docs, but spreads the rules across six files and runs only
after workspace trust.

**The table, from the incidents in §1.2:**

| caller | blocked writes |
|---|---|
| `tester` | `src_python/`, `include/`, `bindings/`, `src/` (the 16b throwaway modules; the `geopackage.py` checkout) |
| `developer`, `perf` | `tests/` (`perf` also keeps `tools/bench.py` and `docs/benchmarks/`) |
| `architect` | anything outside `docs/` (the increment-22 prototype) |
| `orchestrator` | anything except its ledger (see below) |
| main session (no `agent_type`) | `ask`, not deny, on `src_python/`, `include/`, `bindings/`, `src/`, `tests/` |
| any subagent | `.md`/`.txt` under `../rasputin_data`, `../rasputin_scratch`, `/tmp`, a job `tmp/`: this enforces the rule against passing messages through those folders |

**Deny for subagents, ask for the main session.** An `ask` inside a subagent
waits for Ola, and Ola is often away while those runs go on. `deny` returns
the reason to the agent, which then hands the work back. That is the correct
outcome, and it needs no human.

**Limits.** The Edit, Write and NotebookEdit arm is exact. The Bash arm has
the same weakness as `guard_governance.py`: it reads the command as text, so a
path built from a variable escapes it. It would catch the recorded incidents
(`git checkout -- src_python/...`, heredocs into `src_python/`), but it is a
tripwire. The `.md`/`.txt` rule for the data folders is a tripwire too, since
a message could go in a `.json`.

**Cost.** One hook of about 150 lines, plus a test suite that plants each
forbidden write for each persona. It is code, so it goes through the loop:
`@architect` settles the table, `@tester` goes red, `@developer` makes it
green.

### The orchestrator ledger (Ola's idea, 2026-09-29)

Ola: "If it is the only one with permission to read and write a ledger, this
could be useful, also for self inspection." R-A removes the Write tool the
ledger needs. Frontmatter cannot give "Write, only to one file", but the R-B
hook can. The combined shape:

- `@orchestrator` gets `Write, Edit` back. The guard allows them only on the
  ledger path, for example `.claude/ledger/orchestrator.md`.
- Every other caller is denied Read, Edit and Write on that path, and the
  Bash arm catches the path in command text.
- The ledger has one step that writes it. Each dispatch appends persona, ask,
  expected file and outcome, and the orchestrator reads the ledger when it
  starts. A ledger nobody settles is worse than none
  (`docs/increments/README.md`, cost constraints). The step is what keeps
  this one from going stale.

Until the guard exists, R-A stands as applied. Restoring Write before the
guard exists would reopen what R-A closed.

**ASK OLA:** (1) Should the R-B guard go through the TDD loop as the next
increment? (2) Should the ledger be tracked (it survives across machines and
is auditable in git) or untracked (private to the orchestrator)? (3) May the
main session read the ledger, so it can relay the ledger to you, or is it
orchestrator-only?
