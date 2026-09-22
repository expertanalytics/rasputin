# Principles

The rules this project works by. Normative only: each entry is a rule, how to
apply it, and a pointer to the incident that produced it. **The incident is not
here.** It is in the log, and a rule that cannot be stated without telling its
story is not yet a principle.

Three kinds of document, three homes, and nothing crosses:

| Kind | Home | Read when |
|---|---|---|
| Principles | this file | before acting |
| Logs | `docs/retrospectives/`, `docs/increments/*` audits | never, for rules |
| Designs | `docs/increments/NN-*.md` | working that increment |

**Status and exercise.** Every principle carries a status and the increment in
which it was last exercised. A principle not exercised in three increments is a
candidate for retirement, reviewed at the next retrospective. This is the only
defence against the register becoming the thing it replaced: this project has
already recorded that a ledger nobody settles is worse than none
(`docs/increments/README.md`, cost constraints).

**The exercise field is set by the retrospective, not by whoever writes or edits
a principle.** Today's values are a seed, not a measurement: they were assigned
by judgement when the register was drafted, and every one reads 5c because 5c
was in flight. That makes them worthless as evidence until a retrospective
replaces them with something checked. Stated here rather than left implicit,
because a resolved value asserted without a measurement is what B2 and B4 exist
to catch, and the register's own first draft committed it.

---

## A. Evidence

### A1 — Run the command before writing the claim

**Rule.** Naming a check is not performing one. Run it, then write it down.
**Apply.** Any sentence of the form "X does Y" about a script, a suite, a gate
or a file.
**Origin.** `debd8d1` shipped a one-word bug as an inherent limitation; `0fa05e3`
is the correction.
**Status.** In force. Last exercised: 5c.

### A2 — Derive the probe set from the code as fixed, not the bug as found

**Rule.** A fix that widens what the code accepts widens what can break it.
Re-running the old probes tests the old code.
**Apply.** After any change to an input domain.
**Origin.** `fdbd532` — five probes verified a glob widened from `*.md` to every
file; all five were UTF-8 because the bug had been.
**Status.** In force. Last exercised: 5b.

### A3 — Make the probe able to fail

**Rule.** If a pass looks the same as a probe that never ran, it measured
nothing.
**Apply.** Prefer a probe whose pass and its own absence look different. Time a
hang from outside the process.
**Origin.** `b0bf129` — a `signal.alarm` timeout raised inside the `except
OSError` it was meant to expose, reporting "no hang" in output produced by the
hang. Reoccurred in 5b: a strict filter admitted zero samples and would have
shipped `flipped 0 | SURVIVORS 0` as evidence.
**Status.** In force. Last exercised: 5b.

### A4 — Do not write "X is verified by Y" until Y has run against a broken X

**Rule.** As stated.
**Apply.** Before crediting any test, gate or probe with covering anything.
**Origin.** A comment credited a NaN fixture with killing a mutant it cannot
kill (`7ece838`).
**Status.** In force. Last exercised: 5c.

### A5 — Four stale artifacts can make a run meaningless

**Rule.** `pytest` does not rebuild the extension; `cmake --build` can skip a
`cp` restore; **`ctest` reports the previous binary after a failed build**;
Python imports a stale `.pyc` when mtime and size are unchanged.
**Apply.** Rebuild and reinstall before any `pytest` measuring C++. `touch`
after any `cp` restore. **Read `cmake --build`'s exit status before `ctest`'s
summary.**
**Origin.** All four met, not read about. Reproductions in
`.claude/REQUIRED-READING.md`.
**Status.** In force. Last exercised: 5c.

---

## B. Claims

### B1 — A claim names the object the code evaluates

**Rule.** Where there is nothing to run, ask one question: is this claim about
the same object the code evaluates?
**Apply.** Comments, design invariants, error messages, commit messages.
**Origin.** The dominant defect class. A §2 gate matched path-shaped keys
against bare CMake tokens (`1807e73`); tests bounded the compiler's FMA choice
rather than the kernel's error (`dd67a68`); an error message named a legal value
when the defect was the type.
**Status.** In force. Last exercised: 5c.

### B2 — Write the rule and the command, not the resolved value

**Rule.** When a claim's truth depends on the state of the tree, write what
resolves it, not what it resolved to. A resolved value is a citation with an
invisible expiry date.
**Apply.** Line numbers, commit hashes, counts, file lists. Prefer a section
name, a symbol, or a `grep`/`git grep` that returns the set.
**Origin.** Increment 7's review rounds found seven prose defects and zero
behavioural ones; four were in the text repairing the previous one. Two
sentences were deleted rather than corrected a third time, and those were the
first repairs that did not fail on landing.
**Status.** In force. Last exercised: 5c.

### B3 — A correction is a change and can carry the next defect

**Rule.** Apply A1 and B1 to the correction itself, before committing it.
**Apply.** Especially when the correction is *about* accuracy.
**Origin.** A citation with a line number written into a note explaining why
line numbers expire; "measured 18" where 18 was the estimate, inside a commit
fixing an estimate-versus-measurement confusion.
**Status.** In force. Last exercised: 5c.

### B4 — An estimate is not a measurement and must not be reported as one

**Rule.** Label which is which, every time both appear.
**Apply.** Reconciliations, LOC tables, cost sections.
**Origin.** 5b's C++ estimator was called "unbiased and imprecise" only after
5b came in at 1.52× and the pooled figure was re-derived with it in the sample.
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
**Origin.** A governance audit found every production file in the repo had
previously landed in the same commit as its test.
**Status.** In force. Last exercised: 5c.

### C3 — A documentation defect found during an increment is fixed in that PR, or it is not recorded

**Rule.** As stated. No debt ledgers.
**Origin.** Increments 2 and 3 each carried a debt section; between them they
cleared nothing, grew from three entries to five, and listed 5 of the 11 defects
an audit later found.
**Status.** In force. Last exercised: 5c.

### C4 — The merge updates the roadmap row, in the same PR

**Rule.** As stated.
**Origin.** Four increments shipped while `ROADMAP.md` said otherwise, because
nothing in the protocol referred to it.
**Status.** In force. Last exercised: 5c (first application).

---

## D. Working with agents

### D1 — One agent at a time

**Rule.** Dispatch serially. Disjoint files do not make parallel dispatch safe.
**Apply.** Queue the work and say what is queued. Prefer doing small mechanical
fixes inline over spawning an agent for each.
**Origin.** Three agents launched together were killed by the same session limit
within a minute; two were lost mid-task with nothing committed. Parallel edits
had previously produced four of six defects in one round.
**Status.** In force. Last exercised: 5c.

### D2 — Verify a handback; do not relay it

**Rule.** Re-run the cheap claims in any report before acting on them or
repeating them.
**Apply.** Especially numbers, and especially when the report contradicts
something you believe.
**Origin.** A reviewer's claim that a branch had deleted a call was relayed
unchecked; the call had been retyped, and the cited line was already stale
before the branch existed.
**Status.** In force. Last exercised: 5c.

### D3 — A compromise between two correct principles goes to `@orchestrator`

**Rule.** Factual errors route to the file's owner. A conflict between two rules
that are each correct is lifted, not settled in the main session and not handed
to a persona chosen ad hoc.
**Origin.** A conflict between "do not expose `--max-rounds`" and "the band
prints `describe()` verbatim" was lifted and dissolved rather than compromised:
the header text was simply wrong.
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
**Origin.** `legacy/bindings.cpp` took a proj4 string as a parameter and bound
CGAL's polygon types directly as the Python API.
**Status.** In force. Last exercised: 5c.

### E4 — Under 700 non-comment production lines per PR

**Rule.** `CLAUDE.md` §2 is the statement; this is the pointer.
**Status.** In force. Last exercised: 5c.
