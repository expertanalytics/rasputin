# Principles

The rules this project works by.

Each entry is a rule, how to apply it, a pointer to where its incident is
recorded, and a status. The incident itself is not here.

**When you learn something, report it in your handback** (`CLAUDE.md` §3,
*The main session dispatches*); it ends up in the log, `docs/retrospectives/`.
Only a retrospective promotes a log entry to a principle here. The file you
happen to have open is not the log, and neither is this one.

Designs live in `docs/increments/NN-*.md` and hold no findings.

**Origin** is a reference, not an account. Commits and files, no retelling.

---

## A. Evidence

### A6 — A rule in a governance file is evidence someone wrote it

**Rule.** A rule in a repository document is evidence an agent wrote a sentence,
not evidence the user made a rule. Name the commit and the human turn before
citing it to the user or enforcing it in a gate.
**Apply.** The git author proves nothing: every commit here carries the user's
identity. The `Co-Authored-By` line and the transcript store are the evidence,
and `git log -S` finds the commit. A prohibition is never a design opinion for
this purpose, because provenance always has a command to run.
**Origin.** `docs/retrospectives/2026-09-24-unauthored-rules.md`.
**Status.** In force.

## B. Claims

### B3 — A correction is a change and can carry the next defect

**Rule.** Check the correction itself as you would the claim it corrects
(`.claude/REQUIRED-READING.md`, *Claims*), before committing it.
**Apply.** Especially when the correction is *about* accuracy, and in a
governance cleanup, where the mandate to delete false rules is what licenses
writing one.
**Origin.** `docs/retrospectives/2026-09-22-increments-7-5b-5c.md` for the
citation cases; `2026-09-24-unauthored-rules.md` for the chain of five commits
in which each repair carried the next defect.
**Status.** In force.

### B4 — An estimate is not a measurement and must not be reported as one

**Rule.** Label which is which, every time both appear.
**Apply.** Reconciliations, LOC tables, cost sections.
**Origin.** `docs/increments/05c-noder-wiring.md`, its reconciliation section.
**Status.** In force.

---

## D. Working with agents

### D2 — Verify a handback; do not relay it

**Rule.** Re-run the cheap claims in any report before acting on them or
repeating them.
**Apply.** Especially numbers, and especially when the report contradicts
something you believe.
**Origin.** `docs/retrospectives/2026-09-22-increments-7-5b-5c.md`.
**Status.** In force.

### D3 — A compromise between two correct principles goes to `@orchestrator`

**Rule.** Factual errors route to the file's owner. A conflict between two rules
that are each correct is lifted, not settled in the main session and not handed
to a persona chosen ad hoc.
**Origin.** `docs/increments/05c-noder-wiring.md`, section 2.
**Status.** In force.
