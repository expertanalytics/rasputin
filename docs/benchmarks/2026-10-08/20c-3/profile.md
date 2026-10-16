# 20c-3: where the land-cover clean-up time goes (Numedalslagen, one cProfile run)

**Method.** One run under `python -m cProfile` of the branch at `16143f4c`
(no source change since the `_core` build of `745bcb0d`; the clean-up is
Python), Numedalslagen with the defaults, the `NUM` arguments of
`scripts/run.sh`. Apple M1 Max, 10 cores, AC power (charged), threads 10, no
warm-up. Settings in the run record: repair 0.05 m, same-class merge on,
simplify 0 (off), outline rule 5 m, features already in EPSG:25833 (no
conversion). Under the profiler: total 14.87 s, `features clip` 9.44 s,
`features clip: clean-up` 9.33 s (62.7 %); the timed run without profiler had
8.78 s of 14.10 s. The profile file and mesh were in
`rasputin_scratch/perf-20c3-prof/` and are not kept; rerun the command above
to regenerate.

## Clean-up (`feature_input.py` `_clean`, 9.33 s), split

| part | calls | seconds | share of clean-up |
|---|---:|---:|---:|
| outline rule `snap_to_outline` (all of it) | 1 | 6.31 | 68 % |
|   `symmetric_difference(polygon, new)` (area-changed measure) | 21 | 2.20 | 24 % |
|   `union_all(rebuilt)` per class, and of `changed` | 26 | 2.17 | 23 % |
|   `_Outline.ring` Python loop over ring edges | 4 523 | 0.84 | 9 % |
|     of it, `STRtree` `dwithin` queries | 4 523 | 0.12 | 1 % |
|     of it, `_place` / `_cut` | 373 / 1 904 | 0.19 / 0.02 | 2 % |
|   `make_valid` of rebuilt shells | 8 | 0.62 | 7 % |
|   `_chains` | 373 | 0.21 | 2 % |
| `coverage_clean` (repair) | 1 | 1.18 | 13 % |
| `_add`: linework clipped to the domain (`intersection`) | 25 | 1.00 | 11 % |
| first clip to the read region (`intersection`, `_polygonal`) | 3 297 | about 0.75 | 8 % |
| same-class merge `coverage_union_all` | 25 | 0.12 | 1 % |
| `coverage_simplify` | 0 | 0 | off by default |
| moves into the computation CRS | 0 | 0 | none (source already EPSG:25833) |

The two `union_all` attributions and the two `intersection` ones are read
from the callers' cumulative times (shapely's decorator hides the direct
caller); they sum to the totals below.

## Top 15 by own time (whole run, 15.39 s under the profiler)

| # | call | calls | own s | share of run |
|---:|---|---:|---:|---:|
| 1 | shapely `symmetric_difference` | 21 | 2.196 | 14.3 % |
| 2 | shapely `union_all` | 26 | 2.172 | 14.1 % |
| 3 | shapely `intersection` | 11 307 | 1.542 | 10.0 % |
| 4 | `_core.refine` | 1 | 1.284 | 8.3 % |
| 5 | shapely `coverage_clean` | 1 | 1.177 | 7.6 % |
| 6 | shapely `make_valid` | 8 | 0.619 | 4.0 % |
| 7 | `_core.refine_strip` | 1 | 0.556 | 3.6 % |
| 8 | `_core.node` | 1 | 0.480 | 3.1 % |
| 9 | numpy `argsort` | 2 | 0.238 | 1.5 % |
| 10 | `feature_input._Outline.ring` | 4 523 | 0.226 | 1.5 % |
| 11 | `cli._masked_pairs` | 1 | 0.217 | 1.4 % |
| 12 | `stats.quality` | 1 | 0.192 | 1.2 % |
| 13 | `cli._chain_masks` | 1 | 0.179 | 1.2 % |
| 14 | tifffile `decode_other` | 1 152 | 0.162 | 1.1 % |
| 15 | `_core.triangulate` | 1 | 0.141 | 0.9 % |

## Top 15 by cumulative time (below the CLI wrappers)

| # | call | cumulative s | share of run |
|---:|---|---:|---:|
| 1 | `cli.mesh` | 15.09 | 98 % |
| 2 | `feature_input.open_features` / `source` | 9.71 | 63 % |
| 3 | `feature_input._take` | 9.46 | 61 % |
| 4 | `feature_input._clean` (the clean-up) | 9.33 | 61 % |
| 5 | `feature_input.snap_to_outline` | 6.31 | 41 % |
| 6 | `cli._dem_mesh` | 3.28 | 21 % |
| 7 | shapely `symmetric_difference` | 2.20 | 14 % |
| 8 | shapely `union_all` | 2.17 | 14 % |
| 9 | shapely `intersection` | 1.54 | 10 % |
| 10 | `_core.refine` | 1.28 | 8 % |
| 11 | `cli._open_dem` | 1.22 | 8 % |
| 12 | shapely `coverage_clean` | 1.18 | 8 % |
| 13 | `feature_input._add` | 1.00 | 6 % |
| 14 | `feature_input._Outline.ring` | 0.84 | 5 % |
| 15 | `mosaic.assemble` | 0.82 | 5 % |

## Reading

The outline rule is two thirds of the clean-up, and most of it is two
whole-class overlays per class: with the same-class merge on, each class is
one large polygon, and `snap_to_outline` rebuilds it with `union_all` and then
takes its `symmetric_difference` against the original only to measure the
area changed, though only rings near the outline moved. The per-edge Python
loop, the `STRtree` queries and the repair are each smaller. No fix is
proposed here.
