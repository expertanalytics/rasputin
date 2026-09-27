# Retrospective: increments 14 to 20b (2026-09-27)

Called by Ola after 20 and 20b landed (PR #97). Ola's aims were to reassess how
well we are doing and how far we use state-of-the-art research. The main
session wrote the assessment, and Ola accepted it and added the last rule.

## Measured state

The 1 m benchmark is in the scratchpad at `bench1m/REPORT.md`. The scaling
sweep is at `scaling/raw.tsv`. Both runs were on battery.

- Quarter circle at 1 m: refine went from 3.29 s (16) to 0.24 s (20b). The
  worst angle went from 0.0117° to 0.396°, and max degree from 43 to 18.
- The whole tile at 1 m holds at about 0.26 s from 14b on. No quality measure
  got worse.
- Every run meets the tolerance, and from 14b on no run has a Delaunay
  violation. With 20 and 20b switched off, the output is byte-identical to
  18's.
- Scaling: refine speeds up at most 2.1× (1 to 20 threads, flat from about 7).
  Roughly half the single-thread time does not parallelise. The likely cause
  is the serial insert and flip phase, which has not been measured.

## Research

- Increments 14, 14b and 18 rebuilt greedy insertion with a per-triangle
  candidate heap and triangle scan-conversion. That is Garland and Heckbert
  1995, and nobody read it first.
- Increment 20 applies Chew/Ruppert refinement on purpose, but snapping to DEM
  nodes breaks its guarantees. It does not use off-centres (Üngör).
- Data-dependent triangulation was tried and parked. See
  `docs/research/data-dependent-triangulation.md`.
- What may be ours: an exact sup-norm guarantee on a DEM together with
  constraints and a domain polygon, deterministic parallel rounds, and
  tolerance-scaled constraint feet. A literature check is still needed before
  claiming any of it. Ola wants the option to publish kept open.

## What went wrong

The main session:
1. **Stated wrong facts with confidence**: the cause of the 0.0117° sliver;
   "tolerance > 0 triggers the Lawson bug", passed on without checking; C1-C3
   described wrongly; a DDT figure reported from a run still in progress.
2. **Pointed Ola to a report before it existed.**
3. **Wrote replies long enough to bury Ola's instructions**, then went silent.
4. **Produced unasked PNGs** and used jargon it had not defined.
5. **Designed the DDT flip rule badly** (local max error), against the
   literature it had just cited.

The process:
6. **Parallelism was a founding goal but was never measured** until after 20b.
7. **Research came after the work**, not before.
8. **No Delaunay oracle on refinement output until 20b.** The Lawson bug
   (c23583b) had been present since 14b.

## What worked

- A 14× speedup on the quarter circle.
- Clear quality gains from 20 and 20b.
- Deterministic output, green CI, and a tolerance guarantee that holds on
  every run.
- The oracle found a real bug, and the fix was proven by tests that go red
  without it.

## Rules from here

1. **Every design starts with a literature check** that cites the method it
   builds on and says what differs.
2. **Every increment's acceptance includes the 1 m benchmark and a scaling
   sweep**, run on AC power, comparable to the previous increment.
3. **Refinement property tests carry the constrained-Delaunay and tolerance
   oracles** on every path.
4. **The main session reports only finished artefacts and verified figures**,
   answers Ola's question first, and keeps it minimal.
5. **Every new round or increment opens with a structured recap and a roadmap
   reminder** (Ola). Ola works across days and at different times, and is not
   always fully present.
