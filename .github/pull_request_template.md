<!--
The three checks below are the ones no gate can make. Each exists because it
failed at least once here, silently, and CI stayed green throughout. CI covers
the rest; if you skip these, nothing covers them.
-->

## Summary



## Checks no gate makes

- [ ] **Red-step scaffolding is gone.** A TDD increment leaves comments saying headers "do not build yet — that is the intended red step". Three of those survived three merges in `tests/cpp/CMakeLists.txt`, each describing headers that existed by the time it merged.
- [ ] **Every prose claim this change touched is still true.** Read the docs the change affects against the code, not against the previous version of the docs. Eleven false or stale statements accumulated across six files before anyone looked — a README advertising a deleted program, three documents naming a test framework four suites explicitly decline, two listing three ctest targets when there were sixteen.
- [ ] **Actual LOC reconciled against the increment's estimate.** Measure it; there is no gate for this. Increment 3's design said "no split", named the seam to use if it overran, and named the likely cause. The overrun happened in exactly that place and the contingency never fired, because nobody re-measured. Non-comment production lines, per `CLAUDE.md` §2.

## Verification

<!-- Local runs are not CI. State what you ran. -->

- [ ] Release + `ctest`
- [ ] Debug + `-fsanitize=address,undefined -fno-sanitize-recover=all`
- [ ] `-ffp-contract=off` — required for anything touching floating-point geometry. Three tests once passed only because clang contracted a determinant to a single `fma`; the ubuntu leg would have gone red.
- [ ] Governance gates in `tools/`
- [ ] `gh pr checks` green — CI is authoritative

## Merge

- [ ] **Merge commit, not squash.** An increment's red commit must stay ahead of its green one in history; squashing destroys the only evidence the test-first claim is verifiable after the fact.
