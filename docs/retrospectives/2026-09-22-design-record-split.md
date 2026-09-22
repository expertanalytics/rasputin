# Test: does splitting logs out of increment records help?

Dated 2026-09-22. A governance proposal, tested against two records, and mostly
refuted.

## The proposal

Increment records were hypothesised to be bloated with log material — incident
accounts, reconciliations, process ledgers — and would shrink if those moved to
`docs/retrospectives/`. `05d-corner-graze.md` was chosen as the test: small,
designed, unstarted.

## Result 1: a fresh design record has almost nothing to move

`05d-corner-graze.md` is 432 lines across thirteen sections. Classifying each:
the defect, the motivating case, the fix, the predicate, the cost to the shipped
header, the effect on 5b, degeneracy policy and prior art are all design. What
remains — "Files and LOC", "Risks" — is design input for scoping, not log.

**Nothing moves.** The split does not help a record that has not been through a
round.

## Result 2: the growth is real, and it is not log material

`05b-noder-driver.md`, measured with `git show <ref>:<path> | wc -l`:

| | lines |
|---|---|
| at its design commit `36e8023` | 1338 |
| at merge, on master | 2359 |

**76% growth during the round.** The eight sections that did not exist at design
time (`comm -13` over the two section lists) are:

- six rulings — the dense per-edge array, the merge at the node-id dedup,
  guarantee 15's oracle, the one-node chain, the bare-word prohibition, the
  renderer question
- one reconciliation
- one status note about a deferred question

Six of eight are **rulings**: they define behaviour, they were discovered during
the round, and they belong in the design record. They are not log material and
moving them would make the design incomplete.

## Result 3: the ledger sections are 16% of the largest record

Measured with `awk` over section boundaries in `05b-noder-driver.md`:

| Section | lines |
|---|---|
| Files and LOC (of which Reconciliation, 87) | 186 |
| Risks | 103 |
| Documentation this PR fixes | 96 |
| **total** | **385 of 2359** |

Of those, "Documentation this PR fixes" is about other files and is a process
ledger. Reconciliations are measurements of a finished thing. Risks and
estimates are design input.

**So the defensible move is about 180 lines of 2359 — under 8%.**

## What this means

The proposal's premise was wrong. Increment records are not large because log
material leaked into them. They are large because **the design genuinely grows
during the round**, by three quarters, as the red and review steps discover
rulings the design did not anticipate.

That is the protocol working. A ruling found at the red step is cheaper than one
found after merge, and the record is where it has to land.

The open question is a different one from the question this test was built to
answer: **is 2359 lines of design proportionate to 675 lines of code?** 3.5:1 may
be right for a correctness-critical geometry kernel and wrong for everything
else. Nothing in the project measures it.

## What to keep from the proposal

- Reconciliations and "documentation this PR fixes" move to the log.
- Rulings stay in the design record wherever in the round they are found.
- The principle register and `docs/retrospectives/` stand on their own merits;
  neither was justified by this hypothesis.
