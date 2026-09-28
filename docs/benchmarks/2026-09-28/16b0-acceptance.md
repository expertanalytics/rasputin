# Increment 16b-0 acceptance (@perf, 2026-09-28): summary

**Verdict: ACCEPTED.** Measured on battery against a battery baseline, back
to back. The branch is `4e83389`; the base is master `dc372a5`. Method,
tables and raw data: `16b0-acceptance/README.md`.

- **1 m benchmark and thread sweep** (`tools/bench.py`, three pairs): no
  change. Pooled over all 42 cells (each cell's median over the three base runs
  against its median over the three 16b-0 runs), the change is a median of
  -0.28 % (range -5.9 % to +2.8 %). bench.py said ACCEPTED on pairs 1 and 3. On pair 2 it
  flagged two tile cells (+6.9 % at 10 threads, +8.0 % at 16); those cells
  were +0.5 % to +2.1 % in the other pairs. The same build moves by up to
  +13.5 % between runs, and `refine` does not compile the changed header.
  The mesh sha256 is equal on both domains in all six runs. Ceiling:
  2.5-2.9x, as the base.
- **`node` on M3's CORINE squares.** Layout B at 48 km: **84.05 s -> 0.1435 s
  (586x)**, against the Acceptance section's 83.1 s. Layout A at
  12/24/36/48 km: 0.0041/0.0174/0.0364/0.0668 s (reference
  0.0041/0.0167/0.0364/0.0670). From 12 to 48 km the segments grow 13.8x;
  `node` grows 16-18x, against 185x on master. A 144 km block (870 k segments,
  layout B) nodes in 2.2 s. The noded output is identical between layouts A
  and B and between master and 16b-0.
- **The admitted worst case.** The ladder at 300 rungs takes 0.262 s
  (master: 210 s); at 1 200 rungs (2.9 M noded edges) it takes 17.3 s. The
  comb of parallel east-west lines is quadratic, as R8 admits: 4x the lines
  costs 14.4-15.6x, and 64 000 lines take 21.4 s. The split between the
  driver and the verifier has not been isolated.
- **The first CORINE baseline** (48 km square in `6603_4`, layout A, 16b-0).
  At 1 m: 6.75 M triangles, `node` 0.068 s, refine 5.48 s, peak RSS 2.35 GB.
  At 10 m, **the quality start makes the mesh 3.61x the featureless one**
  (480 961 triangles, against 212 255 with the quality start off and 133 379
  without features). That is 20c's input. `node` is no longer the cost.
- **Not measured**: a verifier-only split on the ladder and comb; master at
  96/144 km and on the large synthetic cases; any AC run.
