Basin (BHO level 2, ottobasin 76): 635,194.5 km² geodesic; OAS: 636,920 km².

## The piece: velhas76949, 11,667.6 km² (geodesic)

Medians of 5 timed runs per tolerance; one more ASCII run for quality and the final check.

| tolerance | triangles | vertices | vertices / grid nodes | triangles / km² | refine | process wall | max RSS | threads |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 8,699,303 | 4,353,814 | 33.6% | 745.6 | 5.24 s | 8.27 s | 3.11 GiB | 10 |
| 2 m | 4,584,784 | 2,296,371 | 17.7% | 393.0 | 2.68 s | 4.65 s | 1.69 GiB | 10 |
| 5 m | 1,576,681 | 792,105 | 6.1% | 135.1 | 0.90 s | 1.98 s | 0.94 GiB | 10 |
| 10 m | 649,845 | 328,614 | 2.5% | 55.7 | 0.42 s | 1.24 s | 0.78 GiB | 10 |
| 20 m | 264,833 | 136,077 | 1.1% | 22.7 | 0.22 s | 0.96 s | 0.69 GiB | 10 |
| 50 m | 74,841 | 41,076 | 0.3% | 6.4 | 0.09 s | 0.77 s | 0.62 GiB | 10 |

| tolerance | worst angle | max degree | Delaunay violations (checked, ambiguous) | within tolerance | control: grid nodes over tol | source nodes over tol, interior | source nodes over tol, boundary strip |
|---:|---:|---:|---:|---|---:|---:|---:|
| 1 m | 0.113° | 14 | 0 (13,044,793, 1,761,596) | True | 0 of 12,957,257 | 2,294,965 of 13,791,972 (16.64%; max 31.7 m) | 6,810 of 27,806 (24.49%; max 16.7 m) |
| 2 m | 0.113° | 14 | 0 (6,873,198, 589,150) | True | 0 of 12,957,257 | 760,657 of 13,791,972 (5.52%; max 31.7 m) | 2,883 of 27,806 (10.37%; max 15.8 m) |
| 5 m | 0.113° | 14 | 0 (2,361,258, 68,718) | True | 0 of 12,957,257 | 117,328 of 13,791,972 (0.85%; max 31.7 m) | 486 of 27,806 (1.75%; max 15.8 m) |
| 10 m | 0.113° | 14 | 0 (971,077, 9,013) | True | 0 of 12,957,257 | 24,383 of 13,791,972 (0.18%; max 31.7 m) | 77 of 27,806 (0.28%; max 14.1 m) |
| 20 m | 0.113° | 14 | 0 (393,590, 1,060) | True | 0 of 12,957,257 | 4,534 of 13,791,972 (0.03%; max 27.9 m) | 2 of 27,806 (0.01%; max 20.6 m) |
| 50 m | 0.113° | 14 | 0 (108,607, 84) | True | 0 of 12,957,257 | 514 of 13,791,972 (0.00%; max 56.6 m) | 0 of 27,806 (0.00%; max 36.0 m) |

Power: AC before and after every tolerance block: True.

Max RSS against triangles over the sweep, least squares: 0.56 GiB + 304 B per triangle (max residual 168 MiB).

## The basin, extrapolated from 200 random 10 km boxes (seed 1)

Mean density over the boxes times the basin's area; the interval is the bootstrap 95 % interval of the mean. Memory floor: the piece's per-triangle RSS slope (304 B) times the basin's triangles, plus the 8.6 GiB canvas of Q9. Refine time: the piece's refine time per triangle (10 threads), times the basin's triangles, i.e. linear.

| tolerance | triangles / km², mean (95 %) | median | p10-p90 | basin triangles (95 %) | basin vertices | refine, linear | memory floor |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 275.2 (244.2-307.7) | 216.9 | 64.7-592.3 | 174.8 M (155.1-195.5 M) | 87.4 M | 1.8 min | 58 GiB |
| 2 m | 132.5 (116.1-149.7) | 91.6 | 26.2-306.8 | 84.2 M (73.7-95.1 M) | 42.1 M | 0.8 min | 32 GiB |
| 5 m | 43.4 (37.6-49.8) | 29.6 | 6.7-102.0 | 27.6 M (23.9-31.6 M) | 13.8 M | 0.3 min | 16 GiB |
| 10 m | 17.4 (14.8-20.0) | 10.9 | 1.7-42.1 | 11.0 M (9.4-12.7 M) | 5.5 M | 0.1 min | 12 GiB |
| 20 m | 6.6 (5.6-7.7) | 3.8 | 0.1-16.1 | 4.2 M (3.5-4.9 M) | 2.1 M | 0.1 min | 10 GiB |
| 50 m | 1.5 (1.2-1.8) | 0.6 | 0.0-4.0 | 0.9 M (0.8-1.1 M) | 0.5 M | 0.0 min | 9 GiB |

The final check of Q6, measured: source-DEM nodes outside the tolerance after meshing the resampled grid, **interior only** (a box's 4-corner edges make its boundary strip unrepresentative of a real outline; the piece gives the strip). Basin count = mean per km² times the basin's area.

| tolerance | share of source nodes | per km², mean (95 %) | basin source nodes over tol | as a share of the basin's vertices |
|---:|---:|---:|---:|---:|
| 1 m | 4.56% | 51.9 (43.6-61.0) | 32.96 M | 38% |
| 2 m | 1.37% | 15.6 (12.8-18.7) | 9.92 M | 24% |
| 5 m | 0.20% | 2.3 (1.8-2.8) | 1.45 M | 10% |
| 10 m | 0.04% | 0.5 (0.4-0.6) | 0.30 M | 5% |
| 20 m | 0.01% | 0.1 (0.1-0.1) | 0.07 M | 3% |
| 50 m | 0.00% | 0.0 (0.0-0.0) | 0.01 M | 2% |

Relief inside a box (z max - z min of its grid): median 191 m, p10 54 m, p90 414 m.

## Check of the box method: 30 random 10 km boxes inside the piece against the piece meshed whole

| tolerance | piece, whole catchment, triangles / km² | boxes inside it, mean (95 %) | inside the interval |
|---:|---:|---:|---|
| 1 m | 745.6 | 786.8 (682.5-892.2) | True |
| 2 m | 393.0 | 416.8 (351.4-485.4) | True |
| 5 m | 135.1 | 144.8 (120.1-171.1) | True |
| 10 m | 55.7 | 60.4 (49.4-71.5) | True |
| 20 m | 22.7 | 24.5 (19.9-29.3) | True |
| 50 m | 6.4 | 5.8 (4.3-7.3) | True |

Relief inside a box: median 281 m.

## Refine thread scaling on the piece at 1 m

| threads | refine (median) | speed-up | process wall (median) |
|---:|---:|---:|---:|
| 1 | 10.95 s | 1.00 | 12.01 s |
| 2 | 7.62 s | 1.44 | 8.82 s |
| 4 | 6.06 s | 1.81 | 7.14 s |
| 6 | 5.51 s | 1.99 | 6.57 s |
| 8 | 5.25 s | 2.09 | 6.29 s |
| 10 | 5.22 s | 2.10 | 6.25 s |

Power: AC before and after: True.

