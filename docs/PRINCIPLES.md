# Principles

The rules this project works by.

Each entry is a rule, how to apply it, a pointer to where its incident is
recorded, and a status. The incident itself is not here.

**When you learn something, write it in the log** — `docs/retrospectives/`.
Only a retrospective promotes a log entry to a principle here. The file you
happen to have open is not the log, and neither is this one.

Designs live in `docs/increments/NN-*.md` and hold no findings.

**Origin** is a reference, not an account. Commits and files, no retelling.

**Last exercised** is set by the retrospective, not by whoever writes or edits a
principle. A principle not exercised in three increments is reviewed for
retirement at the next retrospective.

---

## A. Evidence

### A1 — Run the command before writing the claim

**Rule.** Run the check, then write the claim. Naming a check is not
performing one, and neither is being confident about it.
**Apply.** Any sentence of the form "X does Y" about a script, a suite, a gate
or a file — and any claim about what the user has asked for or ruled, which has
a command too. See A6.
**Origin.** `debd8d1`, corrected by `0fa05e3`.
**Status.** In force. Last exercised: 5c.

### A2 — Derive the probe set from the code as fixed, not the bug as found

**Rule.** Choose the probes from what the code accepts **now**, not from what
broke. Re-running the old probes tests the old code.
**Apply.** After any change to an input domain.
**Origin.** `fdbd532`.
**Status.** In force. Last exercised: 5b.

### A3 — Make the probe able to fail

**Rule.** Before trusting a pass, check the probe could have failed. A pass
that looks the same as a probe that never ran measured nothing.
**Apply.** Prefer a probe whose pass and its own absence look different. Time a
hang from outside the process.
**Origin.** `b0bf129`. Recurrence: 5b's ring search, `05b-noder-driver.md`,
"The bounded ring search".
**Status.** In force. Last exercised: 5b.

### A4 — Do not write "X is verified by Y" until Y has run against a broken X

**Rule.** As stated.
**Apply.** Before crediting any test, gate or probe with covering anything.
**Origin.** `7ece838`.
**Status.** In force. Last exercised: 5c.

### A5 — Four stale artifacts can make a run meaningless

**Rule.** Rebuild and reinstall before any `pytest` that measures C++.
`touch` after any `cp` restore. Read `cmake --build`'s exit status before you
read `ctest`'s summary. A clean-looking run proves nothing until you have.
**Apply.** Four artifacts go stale silently: `pytest` does not rebuild the
extension, `cmake --build` can skip a `cp` restore, `ctest` reports the previous
binary after a failed build, and Python imports a stale `.pyc` when mtime and
size are unchanged.
**Origin.** `.claude/REQUIRED-READING.md`, which carries the reproduction of
each.
**Status.** In force. Last exercised: 5c.

---

### A6 — A rule in a governance file is evidence someone wrote it

**Rule.** A rule in a repository document is evidence an agent wrote a sentence,
not evidence the user made a rule. Name the commit and the human turn before
citing it to the user or enforcing it in a gate.
**Apply.** The git author proves nothing: every commit here carries the user's
identity. The `Co-Authored-By` line and the transcript store are the evidence,
and `git log -S` finds the commit. A prohibition is never a design opinion for
this purpose, because provenance always has a command to run.
**Origin.** `docs/retrospectives/2026-09-24-unauthored-rules.md`, which measures
48 of 48 governance commits as agent-written.
**Status.** In force.

## B. Claims

### B1 — A claim names the object the code evaluates

**Rule.** Where there is nothing to run, ask one question: is this claim about
the same object the code evaluates?
**Apply.** Comments, design invariants, error messages, commit messages.
**Origin.** `1807e73`, `dd67a68`, `7473712`, `d86668f`.
**Status.** In force. Last exercised: 5c.

### B2 — Write the rule and the command, not the resolved value

**Rule.** When a claim's truth depends on the state of the tree, write what
resolves it, not what it resolved to. A resolved value is a citation with an
invisible expiry date.
**Apply.** Line numbers, commit hashes, counts, file lists. Prefer a section
name, a symbol, or a `grep`/`git grep` that returns the set.
**Origin.** `docs/increments/07-edge-properties.md`, its reconciliation section;
`docs/retrospectives/2026-09-22-increments-7-5b-5c.md`.
**Status.** In force. Last exercised: 5c.

### B3 — A correction is a change and can carry the next defect

**Rule.** Apply A1 and B1 to the correction itself, before committing it.
**Apply.** Especially when the correction is *about* accuracy, and in a
governance cleanup, where the mandate to delete false rules is what licenses
writing one.
**Origin.** `docs/retrospectives/2026-09-22-increments-7-5b-5c.md` for the
citation cases; `2026-09-24-unauthored-rules.md` for the chain of five commits
in which each repair carried the next defect.
**Status.** In force. Last exercised: 5c.

### B4 — An estimate is not a measurement and must not be reported as one

**Rule.** Label which is which, every time both appear.
**Apply.** Reconciliations, LOC tables, cost sections.
**Origin.** `docs/increments/05c-noder-wiring.md`, its reconciliation section.
**Status.** In force. Last exercised: 5c.

---

## C. The increment loop

### C1 — Design, red, green, review, in that order

**Rule.** `@architect` designs; `@tester` commits a failing suite; `@developer`
makes it pass; `@reviewer` audits. Detail in `docs/increments/README.md`.
**Status.** In force. Last exercised: 5c.

### C2 — The red commit touches no production file; the green commit touches no test file

**Rule.** As stated. Verified per commit, on master, after every merge.
**Apply.** `git show <sha> --name-only`. Merge commits, never squashes: a squash
destroys the trace irreversibly.
**Origin.** `docs/increments/README.md`.
**Status.** In force. Last exercised: 5c.

### C3 — A documentation defect found during an increment is fixed in that PR, or it is not recorded

**Rule.** As stated. No debt ledgers.
**Origin.** `docs/increments/README.md`, cost constraints.
**Status.** In force. Last exercised: 5c.

### C4 — The merge updates the roadmap row, in the same PR

**Rule.** As stated.
**Origin.** `docs/increments/README.md`, the merge section; PR #78.
**Status.** In force. Last exercised: 5c (first application).

---

## D. Working with agents

### D1 — One agent at a time

**Rule.** Dispatch serially. Disjoint files do not make parallel dispatch safe.
**Apply.** Queue the work and say what is queued. Prefer doing small mechanical
fixes inline over spawning an agent for each.
**Origin.** `docs/retrospectives/2026-09-22-increments-7-5b-5c.md`.
**Status.** In force. Last exercised: 5c.

### D2 — Verify a handback; do not relay it

**Rule.** Re-run the cheap claims in any report before acting on them or
repeating them.
**Apply.** Especially numbers, and especially when the report contradicts
something you believe.
**Origin.** `docs/retrospectives/2026-09-22-increments-7-5b-5c.md`.
**Status.** In force. Last exercised: 5c.

### D4 — Finish an open PR before starting new work

**Rule.** Merge it, or say why it is parked. Green and waiting is not finished,
and neither is opened. Report its checks until it settles.
**Apply.** Before branching anything new, `gh pr list --state open`. A PR of
yours in that list is the next task, ahead of whatever seems more interesting.
**Origin.** `docs/retrospectives/2026-09-22-increments-7-5b-5c.md`.
**Status.** In force. Last exercised: not yet.

### D3 — A compromise between two correct principles goes to `@orchestrator`

**Rule.** Factual errors route to the file's owner. A conflict between two rules
that are each correct is lifted, not settled in the main session and not handed
to a persona chosen ad hoc.
**Origin.** `docs/increments/05c-noder-wiring.md`, section 2.
**Status.** In force. Last exercised: 5c.

---

## E. Boundaries

### E1 — The working tree is yours; the remote is the user's

**Rule.** Editing, spawning, building, testing and committing run on an
instruction's momentum. `git push`, `gh pr create`, `gh pr merge`, any
force-push or history rewrite, and any change to permissions need a fresh yes.
**Apply.** One named act, one occurrence. Ask again if the tree has changed.
**Status.** In force. Last exercised: 5c.

### E2 — `@reviewer` runs before a branch's first push, and again before a push that follows findings

**Rule.** Not conditional on the branch containing code.
**Status.** In force. Last exercised: 5c.

### E3 — The I/O boundary

**Rule.** The C++ core never opens a file, sees a path or links a codec, and CRS
never crosses into it. `CLAUDE.md` §2 is the statement; this is the pointer.
**Origin.** `CLAUDE.md` §2; `docs/increments/05c-noder-wiring.md`, prior art.
**Status.** In force. Last exercised: 5c.

### E4 — Under 700 production lines per PR

**Rule.** `CLAUDE.md` §2 is the statement; this is the pointer.
**Status.** In force. Last exercised: 5c.
