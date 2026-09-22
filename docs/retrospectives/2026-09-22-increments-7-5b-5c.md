# Retrospective: increments 7, 5b and 5c

Dated 2026-09-22, covering `4834568..de3b049` — the edge property sets, the
noder, and wiring the noder through.

This is a log entry. It records what happened once. The rules it produced live
in `.claude/REQUIRED-READING.md` and `docs/increments/README.md`; where a rule
exists, this file points at it rather than restating it.

## What shipped

A road crossing a river produces a mesh with a constructed vertex at the
intersection. At a coarse spacing a road running along a river merges to one
edge carrying both properties by union. Before increment 5b, that input returned
`CdtStatus::NotNoded` and zero triangles.

## Cost, measured

| | commits | touch production | touch tests | touch prose |
|---|---|---|---|---|
| 7 | 26 | 8 | 9 | 16 |
| 5b | 25 | 4 | 10 | 11 |
| 5c | 23 | 4 | 7 | 12 |
| **total** | **74** | **16** | **26** | **39** |

Reproduce with `git rev-list --no-merges <range>` and classifying each commit by
the directories it touches.

**More than half of all commits maintain documents describing the code.** That
is the project's largest single cost and it was not being measured.

## The defect distribution, and what it implies

Across all three increments the review rounds found **zero behavioural defects**
and approximately **forty prose defects**. Not one blocking finding was in code.

Two readings, both true:

1. The TDD protocol works. Code defects were caught by the suites before review
   saw them — the `bool`-passes-`isinstance` hole, the `-fPIC` failure that is
   invisible on Mach-O, the guarantee-15 oracle that would have gone red on
   correct output. Review is the last line and the earlier lines held.
2. Review capacity is spent almost entirely on prose. That is a rational
   allocation only if the prose is load-bearing, and some of it was not.

## The dominant failure mode

Not wrong code. **A claim that names a different object than the code
evaluates.** Instances: an error message naming a legal value when the defect
was the type; a validation pattern credited with closing an injection hole it
cannot reach; a citation resolving to a blank line; a design section telling
readers three sentences had been corrected when none had.

The sharp form is **the correction that introduces the next defect**. Recorded
instances from this period:

- A citation with a line number, written into a note explaining why line numbers
  expire.
- "Measured 18 lines" where 18 was the estimate and 8 the measurement, inside a
  commit correcting an estimate-versus-measurement confusion.
- A stale line number quoted verbatim while being fixed, which created a new
  at-risk citation immediately.

What ended the sequence was not more care. It was the rule now in
`.claude/REQUIRED-READING.md`: when a claim's truth depends on the state of the
tree, write the rule and the command that resolves it, not the resolved value.
Two sentences were **deleted** rather than corrected a third time, and those
were the first repairs in the sequence that did not fail on landing.

## Four stale-artifact hazards, all met rather than read about

All four are in `.claude/REQUIRED-READING.md` with their reproductions.

1. `pytest` does not rebuild the C++ extension — `ninja` absent, editable
   auto-rebuild no-ops.
2. `cmake --build` can skip a `cp`-based restore when mtime and size match.
3. **`ctest` reports the previous binary after a failed build.** A failing suite
   reports green. Demonstrated: three targets failing gives 44 failures; a
   planted `int main(){return 0;}` at one registered path gives 40.
4. Python imports a stale `.pyc` when an edit preserves mtime and size.

Third is the most dangerous: the other three make a *passing* run meaningless,
that one makes a *failing* suite look green.

## What the process caught that a person would not have

- A test oracle that would have gone **red on correct output**. Latent:
  instrumenting both spellings found 0 divergences in ~1000 generated segments,
  so no seed would ever have reached it. Found by reading the rule against the
  code, not by running anything.
- `RingDegenerateAfterSnap` demoted to a self-check by bounded exhaustive
  search: two independently written programs, different enumerations,
  simplicity filters and lattices, `FLIPPED == grazed` and `SURVIVORS == 0` in
  every block.
- That ruling's hidden dependency on increment 2's implementation choice:
  `orientation<K>(A,B,C,B)` is `Clockwise` while `signed_area` of the same ring
  is identically zero. A shoelace spelling would reopen a status the increment
  had just closed, and nothing in the tree would catch it.

## What went wrong in how the work was run

- **Parallel dispatch.** Three agents launched together, all killed by the same
  session limit, two lost mid-task with nothing committed. Serial working had
  already been asked for. Recorded as a standing rule.
- **Relaying instead of verifying.** A reviewer's claim that the branch had
  deleted a call six documents argued from was passed on without checking. The
  call had been retyped, not removed, and the cited line was already stale on
  master before the branch existed.
- **A finished PR abandoned, and new work built on the tree without it.**
  PR #79 corrected `ROADMAP.md`'s intro, went green, and sat unmerged and
  unreported while this retrospective, the principle register and a governance
  test were written — all branched from a master that still carried the text
  #79 removed. The register's argument for separating principle from incident
  was drafted in a tree whose own roadmap opened with an anecdote about itself.

  Not caught by this retrospective. The user caught it, after it was written.
  A retrospective that misses a failure in progress is the strongest available
  evidence that the failure is structural rather than one of attention: nothing
  in the working sequence tracks an open PR to completion. "Push it and open the
  PR" was treated as the end of a task rather than the middle of one.

  Promoted to principle D4 from here, per the register's own rule that a
  retrospective is what promotes a log entry.
- **Routing decided at the wrong level.** A conflict between two correct
  principles was nearly handed to a persona chosen ad hoc. Sent to
  `@orchestrator` instead, which dissolved it rather than compromising: the
  header text was simply wrong and its own increment's design already said so.

## What worked and is worth keeping

Verifying handbacks by re-running them. This caught something in nearly every
round — a stale sanitizer count, a wrong mutant count, an estimate reported as a
measurement.

Agents refusing their briefs. A tester rewrote an assertion after finding it
protected something different from what its name claimed. An architect declined
to record a reconciliation it could not reproduce, and said so. Both were better
outcomes than compliance.

The red/green trace surviving every merge, verified on master after each: red
commits touch no production file, green commits touch no test file, red is an
ancestor of green.
