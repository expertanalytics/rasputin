# Increment 16c acceptance on Bygdin: summary (@perf, 2026-09-29)

**Verdict: ACCEPTED** on the design's three pass criteria. AC power, Apple M1
Max, commit `f7d5f14`. Method, tables and scripts:
`2026-09-29/bygdin-landcover/README.md`.

Reduced Bygdin catchment (from 22), DTM10, CORINE 2018 Norway extract with
`--features-map corine`:

| criterion | 10 m | 1 m | pass |
|---|---|---|---|
| `vtkPolyDataReader` reads it; `land_cover_code` on every cell, 0 on every line | 92,053 cells | 1,159,508 cells | yes |
| seven codes, no code 0, each share within 0.01 pp of CORINE clipped | largest gap 0.000001 pp | the same | yes |
| design's centroid oracle, triangles checked / disagreements | 83,171 / 0 | 1,145,434 / 0 | yes |
| I1: unconstrained edges with differing codes | 0 of 116,398 | 0 of 1,705,117 | yes |

Shares: 333 42.34 %, 332 21.32 %, 322 17.00 %, 512 16.46 %, 335 2.44 %,
412 0.33 %, 142 0.10 %, equal to 22's table.

Time, not gated (median of 3): the `land cover` phase takes 0.071 s at 10 m
(83,171 triangles, total 2.08 s) and 1.296 s at 1 m (1,145,434 triangles,
total 3.89 s). At 1 m that is more than `refine` (0.559 s), which is the
condition the design's R2 sets for moving `regions` into the core.

No bench or thread sweep: 16c's commits change no C++ (design, Acceptance).
Pictures: `bygdin_landcover_oblique.png` (from the south-east, 2× vertical)
and `bygdin_landcover_top.png`. ParaView preset: `corine_natural.json`. The
manual ParaView import (design item 6) is left for Ola.
