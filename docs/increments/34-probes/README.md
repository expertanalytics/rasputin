# Design probes for increment 34 (`@architect`, 2026-10-09)

Run by hand with the DTM10 tiles at `../rasputin_data/DTM10_UTM33_20260925`
and Copernicus GLO-30 at `../rasputin_data/germany_glo30` (the scripts name
the absolute paths). Python 3.14 in a scratch virtual environment with NumPy,
SciPy (the simulations only; not a rasputin dependency), tifffile, imagecodecs
(GLO-30's compressed tiles), shapely, pyproj, pydantic, typer and pytest;
`rasputin` is master at `9ba38490`: its `src_python/tin_engine` with the
`_core` that `@perf` built for increment 33's acceptance run, whose C++ sources
equal master's (`git diff ee4fb373 9ba38490 -- include src bindings
CMakeLists.txt` is empty). No C++ was built for this design. The Mac was on AC
power (`pmset -g batt`: "Now drawing from 'AC Power'").

Every slope below is section 3's (`slope_stats.py`, `horn`): Horn's 3 by 3
differences with a missing neighbour filled so that a plane is exact.
Round 1's probes used NumPy's edge padding instead (design review round 1,
B1); every output here was re-run with the stated rule.

The two windows, 10 by 10 km of DTM10 (1001 by 1001 nodes, EPSG:25833):
*Romsdalen*, x 118 000-128 000, y 6 938 000-6 948 000 (tile `6901_3`, the
valley walls below Trollveggen), and *Geilo-Ål*, x 127 000-137 000,
y 6 735 000-6 745 000 (tile `6701_3`, the north-west corner of increment 33's
section of Hallingdal).

1. `slope_stats.py`: the share of nodes at or above each angle by four slope
   measures, and Horn's slope with Gaussian noise added to z. Output, word for
   word (the Geilo-Ål line covers 33's whole section, x 127 000-150 000,
   y 6 720 000-6 745 000; the ns/node is NumPy's time for section 3's rule):

   ```
   romsdalen: 1001 x 1001 nodes, z 53..1797 m
     horn    15: 71.8% 20: 60.4% 25: 48.7% 30: 37.1% 35: 25.6% 40: 16.6% 45: 11.1%  (49.4 ns/node)
     central 15: 72.0% 20: 60.5% 25: 48.8% 30: 37.2% 35: 25.7% 40: 16.6% 45: 11.2%
     cellmax 15: 83.8% 20: 74.3% 25: 63.1% 30: 51.3% 35: 38.8% 40: 27.5% 45: 20.2%
     horn_s  15: 71.6% 20: 60.0% 25: 48.2% 30: 36.7% 35: 25.5% 40: 16.5% 45: 10.8%
     noise 0.1 m: 15: 71.8% 20: 60.4% 25: 48.7% 30: 37.1% 35: 25.6% 40: 16.6% 45: 11.1%  mean |change| 0.16 deg
     noise 0.3 m: 15: 71.9% 20: 60.4% 25: 48.7% 30: 37.1% 35: 25.7% 40: 16.6% 45: 11.1%  mean |change| 0.47 deg
     noise 1.0 m: 15: 72.2% 20: 60.7% 25: 49.0% 30: 37.3% 35: 26.0% 40: 16.8% 45: 11.2%  mean |change| 1.57 deg
     ground under 5 deg: 9.2%; p50/p90/p99 horn 24.4/46.3/66.2 deg
   geilo-al: 2501 x 2301 nodes, z 414..1243 m
     horn    15: 18.8% 20:  9.3% 25:  4.6% 30:  2.2% 35:  1.0% 40:  0.5% 45:  0.2%  (47.3 ns/node)
     central 15: 19.6% 20:  9.8% 25:  4.8% 30:  2.3% 35:  1.0% 40:  0.5% 45:  0.3%
     cellmax 15: 38.6% 20: 21.7% 25: 11.7% 30:  6.1% 35:  3.0% 40:  1.5% 45:  0.8%
     horn_s  15: 17.3% 20:  8.4% 25:  4.0% 30:  1.9% 35:  0.8% 40:  0.4% 45:  0.2%
     noise 0.1 m: 15: 18.9% 20:  9.3% 25:  4.6% 30:  2.2% 35:  1.0% 40:  0.5% 45:  0.2%  mean |change| 0.19 deg
     noise 0.3 m: 15: 19.0% 20:  9.4% 25:  4.6% 30:  2.2% 35:  1.0% 40:  0.5% 45:  0.2%  mean |change| 0.58 deg
     noise 1.0 m: 15: 20.5% 20: 10.0% 25:  4.9% 30:  2.3% 35:  1.0% 40:  0.5% 45:  0.3%  mean |change| 1.94 deg
     ground under 5 deg: 31.7%; p50/p90/p99 horn 7.8/19.5/34.7 deg
   ```

   Horn's differences are linear in z, so its slope noise scales with
   sigma / spacing: 1.0 m of noise on the 10 m grid stands for 0.1 m of noise
   on a 1 m grid (DTM1), which the probe did not have.

   1b. `border_check.py`: section 3's rule at the border and beside NoData,
   on planes (21 by 21 nodes, one NoData node at row 10, column 10), against
   round 1's stated rule (a missing neighbour takes the node's own z); then
   round 1's probe computation against section 3's rule on the Romsdalen
   window. Output, word for word:

   ```
   plane 40.0 deg towards 0.0 deg, 10 x 10 m:
     interior               round 1  40.000  section 3  40.000  class 80
     border                 round 1  32.183  section 3  40.000  class 80
     corner                 round 1  18.350  section 3  40.000  class 80
     next to NoData         round 1  32.183  section 3  40.000  class 80
     corner next to NoData  round 1  36.563  section 3  40.000  class 80
     section 3, every valid node: 40.000000000 .. 40.000000000 deg, classes [80, 81]
   plane 40.0 deg towards 30.0 deg, 10 x 10 m:
     interior               round 1  40.000  section 3  40.000  class 80
     border                 round 1  30.284  section 3  40.000  class 80
     corner                 round 1  12.947  section 3  40.000  class 80
     next to NoData         round 1  34.520  section 3  40.000  class 80
     corner next to NoData  round 1  39.581  section 3  40.000  class 81
     section 3, every valid node: 40.000000000 .. 40.000000000 deg, classes [80, 81]
   plane 37.3 deg towards 30.0 deg, 10 x 5 m:
     interior               round 1  37.300  section 3  37.300  class 75
     border                 round 1  27.932  section 3  37.300  class 75
     corner                 round 1  12.663  section 3  37.300  class 75
     next to NoData         round 1  31.982  section 3  37.300  class 75
     corner next to NoData  round 1  37.980  section 3  37.300  class 75
     section 3, every valid node: 37.300000000 .. 37.300000000 deg, classes [75]
   one-row grid, plane 40 deg towards 30 deg: 36.005 deg (the north-south part is lost)
   romsdalen: classes differ at 3565 nodes, all on the outer ring: True; ring 4000 nodes; at 30 deg or more on the ring: round 1 probe (numpy edge padding) 27.0%, section 3 41.0%, interior 38.2%
   ```

   A slope that is exactly on a class boundary (40 degrees is class 80) can
   come out one class higher by rounding (81): never lower.

   1c. `resample_check.py`: the slope on the resampled path. A window of
   GLO-30 N47 E011 (lon 11.0-11.5, lat 47.3-47.55: the Wetterstein and
   Karwendel) on its own grid, against a square grid in EPSG:25832 at the
   source's north-south spacing rounded to whole metres (as
   `target_grid.default_spacing`), bilinear from the source, over the window
   less 0.02 degrees on every side. Output, word for word:

   ```
   source window 900 x 1800 nodes, 20.96 m (east-west) by 30.88 m; target 785 x 1097 nodes at 31 m square
     at 15 deg or more: source  76.7%  target  76.3%
     at 25 deg or more: source  55.6%  target  54.4%
     at 30 deg or more: source  42.2%  target  41.0%
     at 35 deg or more: source  27.4%  target  26.5%
     at 40 deg or more: source  15.5%  target  15.2%
     at 45 deg or more: source   8.8%  target   8.6%
     median slope: source 27.2 deg, target 26.7 deg; p99 source 59.5, target 58.9
   ```

2. `greedy_sim.py`: greedy insertion to a slope-dependent tolerance, simulated
   in Python (the rules are in its docstring and in section 3 of the design),
   N = 2 m, F = 10 m. Commands: `python greedy_sim.py --size 1001`,
   `python greedy_sim.py --size 1001 --start 25 --end 35 --rules
   node,tri,tri+halo`, `python greedy_sim.py --case geilo-al --size 1001`
   (the 1001-node window from the corner of 33's section, the Geilo-Ål window
   above). Output, word for word:

   ```
   romsdalen 1001x1001 nodes; N 2.0 F 10.0 ramp 30.0..30.0 deg; steep (>= end) 37.1%
     uniform-F  triangles    56193 rounds  23 nodes over t(s(n)) 17.739% lazy-asked  50.5%  (8 s)
     uniform-N  triangles   412641 rounds  26 nodes over t(s(n))  0.000% lazy-asked  19.7%  (30 s)
     node       triangles   279337 rounds  24 nodes over t(s(n))  0.000% lazy-asked  22.8%  (20 s)
     tri        triangles   321830 rounds  25 nodes over t(s(n))  0.000% lazy-asked  21.4%  (23 s)
     tri+halo   triangles   344937 rounds  25 nodes over t(s(n))  0.000% lazy-asked  20.8%  (25 s)
     normal     triangles   271457 rounds  31 nodes over t(s(n))  1.032% lazy-asked  23.3%  (29 s)
   romsdalen 1001x1001 nodes; N 2.0 F 10.0 ramp 25.0..35.0 deg; steep (>= end) 25.6%
     node       triangles   253061 rounds  23 nodes over t(s(n))  0.000% lazy-asked  24.3%  (17 s)
     tri        triangles   302486 rounds  26 nodes over t(s(n))  0.000% lazy-asked  22.3%  (24 s)
     tri+halo   triangles   325818 rounds  25 nodes over t(s(n))  0.000% lazy-asked  21.5%  (24 s)
   geilo-al 1001x1001 nodes; N 2.0 F 10.0 ramp 30.0..30.0 deg; steep (>= end) 4.8%
     uniform-F  triangles    16460 rounds  21 nodes over t(s(n))  2.359% lazy-asked  60.6%  (6 s)
     uniform-N  triangles   197623 rounds  24 nodes over t(s(n))  0.000% lazy-asked  22.9%  (15 s)
     node       triangles    42016 rounds  20 nodes over t(s(n))  0.000% lazy-asked  40.2%  (6 s)
     tri        triangles    59394 rounds  23 nodes over t(s(n))  0.000% lazy-asked  35.1%  (8 s)
     tri+halo   triangles    73816 rounds  23 nodes over t(s(n))  0.000% lazy-asked  31.4%  (9 s)
     normal     triangles    39019 rounds  25 nodes over t(s(n))  0.449% lazy-asked  40.6%  (8 s)
   ```

   "nodes over t(s(n))" is the share of all nodes whose error exceeds the
   tolerance at their own slope when the rule has converged: the guarantee,
   which the `normal` rule (the triangle's own plane) does not keep.
   "lazy-asked" is the share of first evaluations of a triangle whose largest
   error lay between N and F, where increment 33's lazy test asks the policy.

3. `boxes.py`: rasputin master's own uniform meshes of the two windows, from
   the window's outline (`--domain`), F = 10 and N = 2. Output, word for word:

   ```
   romsdalen 10 ['38935 triangles. Every DEM node inside the mesh is within 10 m of it (largest difference 9.9992 m).']
   romsdalen 2 ['321985 triangles. Every DEM node inside the mesh is within 2 m of it (largest difference 2 m).']
   geilo-al 10 ['10502 triangles. Every DEM node inside the mesh is within 10 m of it (largest difference 9.996 m).']
   geilo-al 2 ['143783 triangles. Every DEM node inside the mesh is within 2 m of it (largest difference 2 m).']
   ```

   The simulation overcounts rasputin by 44 % and 28 % on Romsdalen, 57 % and
   37 % on Geilo-Ål (uniform F and N), so item 5 uses only its relative
   position between the uniform counts.

   The whole tile `6901_3` (2 550 km², `--dem` the tile, no domain), one run
   each with `--stats`: `--tolerance 10` 615 459 triangles, refine 0.375 s,
   scan 0.192 s, total 0.677 s, peak 0.47 GB; `--tolerance 2` 5 916 344
   triangles, refine 3.744 s, scan 0.891 s, total 6.321 s, peak 2.19 GB
   (`/usr/bin/time -l`, "peak memory footprint").

   The quick check's catchments at default flags, one run each with `--stats`
   (the commands of `docs/benchmarks/quick/cases.toml`, 10 threads):
   Numedalslagen, 17 882 by 15 667 nodes (2.80e8), 1 118 006 triangles, total
   11.12 s, decode 1.78 s, refine 1.37 s, scan 0.33 s, peak 2.60 GB; Lagan,
   target grid 3 638 by 4 777 nodes (1.74e7), 798 554 triangles, total
   21.56 s, refine 0.75 s, final check scan 0.34 s, peak 1.83 GB.

4. `fixture_figures.py`: the figures of section 9's fixtures V1 and V2.
   Command: `python fixture_figures.py PYTHON RASPUTIN_PY PKG`. Output, word
   for word:

   ```
   V1: 65 x 65 nodes, 10 x 10 m, z -3.0..204.4 m
     classes: >=15: 38.5% >=25: 36.6% >=30: 35.4% >=35: 35.4% >=40: 24.3%
     nodes held to N=2 by step 30: 1496 of 4225; by ramp 25..35, below F: 1545
     rasputin uniform 10: 28 triangles
     rasputin uniform 2: 796 triangles
     simulation uniform-F 32 triangles, 7 rounds, nodes over 9.964%
     simulation uniform-N 1142 triangles, 21 rounds, nodes over 0.000%
     simulation node      408 triangles, 11 rounds, nodes over 0.000%
     simulation tri+halo  594 triangles, 11 rounds, nodes over 0.000%
   V2: 33 x 41 nodes, 10 x 5 m, z -3.0..169.2 m
     classes: >=15: 51.2% >=25: 50.0% >=30: 48.9% >=35: 48.8% >=40: 32.0%
     nodes held to N=2 by step 30: 661 of 1353; by ramp 25..35, below F: 676
     rasputin uniform 10: 6 triangles
     rasputin uniform 2: 112 triangles
   ```

   4b. `mutant_fixtures.py`: the layouts that let tests 4 (b) and 6 kill M6
   and M3, and V1's line (its docstring says how each figure is counted).
   Command: `python mutant_fixtures.py PYTHON RASPUTIN_PY PKG`. Output, word
   for word:

   ```
   M6: --tolerance 2 on V1 cut by the edge (640.0, -300.0) to (150.0, -640.0): 339 vertices, 10 feet; DEM-node vertices on the wall within 5 m of the edge: 3, of them between eps(2) and eps(10) from it (footed under M6): 3
   M3: V1 cells whose corners straddle 30 deg: 129; check points (4 per cell, seed 7) in them nearest a corner below 30 deg: 246; of those with noise above N = 2 m: 136
   line: V1 nodes within 200 m of it 2871 of 4225; held tighter by the line than by the slope 1752, by the slope than by the line 1341
   ```

5. `estimate.py`: the estimates of sections 6 and 10 from items 2 and 3, and
   the whole tile's shares and NumPy's cost per node. Output, word for word:

   ```
   romsdalen step 30    node     position 0.626  rasputin estimate    216130
   romsdalen step 30    tri      position 0.745  rasputin estimate    249873
   romsdalen step 30    tri+halo position 0.810  rasputin estimate    268222
   romsdalen step 30    normal   position 0.604  rasputin estimate    209873
   romsdalen ramp 25-35 node     position 0.552  rasputin estimate    195265
   romsdalen ramp 25-35 tri      position 0.691  rasputin estimate    234513
   romsdalen ramp 25-35 tri+halo position 0.756  rasputin estimate    253040
   geilo-al  step 30    node     position 0.141  rasputin estimate     29303
   geilo-al  step 30    tri      position 0.237  rasputin estimate     42088
   geilo-al  step 30    tri+halo position 0.317  rasputin estimate     52699
   geilo-al  step 30    normal   position 0.125  rasputin estimate     27099
   tile 6901_3: 5051 x 5051 nodes; NumPy, one thread: section 3's rule 40.4 ns/node, plain Horn (the interior's arithmetic) 12.2 ns/node
     shares: >=25: 31.0% >=30: 22.3% >=35: 14.6%
     node rule, ramp 25-35: 1363235 to 3543164 triangles
   geilo-al window, lines with the slope: position 0.0398  rasputin estimate 18117 (lines alone 3687)
   ```

6. `lines_sim.py`: the slope together with increment 33's lines on the
   Geilo-Ål window, with the quick check's `geilo-al-ramp` flags (F 20, the
   lines' N 1, ramp 0 to 3000 m) and the slope at 2 m from 25 to 35 degrees.
   Command: `python lines_sim.py PYTHON RASPUTIN_PY PKG OUT_DIR` (after
   `boxes.py` has written the window's domain to `OUT_DIR`). Output, word for
   word:

   ```
   geilo-al window: nodes within 3000 m of the line 11.1%, at 35 deg or more 2.1%
     simulation uniform-1   473736 triangles
     simulation lines       5667 triangles
     simulation lines+slope 24288 triangles
     rasputin uniform-1   ['366418 triangles. Every DEM node inside the mesh is within 1 m of it (largest difference 1 m).']
     rasputin lines       ['3687 triangles. Every DEM node inside the mesh is within 20 m of it (largest difference 19.996 m).']
   ```
