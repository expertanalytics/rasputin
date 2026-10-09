# Design probes for increment 34 (`@architect`, 2026-10-09)

Run by hand with the DTM10 tiles at `../rasputin_data/DTM10_UTM33_20260925`
(the scripts name the absolute path). Python 3.14 in a scratch virtual
environment with NumPy, SciPy (the simulation only; not a rasputin
dependency), tifffile, shapely, pyproj, pydantic, typer and pytest;
`rasputin` is master at `9ba38490`: its `src_python/tin_engine` with the
`_core` that `@perf` built for increment 33's acceptance run, whose C++ sources
equal master's (`git diff ee4fb373 9ba38490 -- include src bindings
CMakeLists.txt` is empty). No C++ was built for this design. The Mac was on AC
power (`pmset -g batt`: "Now drawing from 'AC Power'").

The two windows, 10 by 10 km of DTM10 (1001 by 1001 nodes, EPSG:25833):
*Romsdalen*, x 118 000-128 000, y 6 938 000-6 948 000 (tile `6901_3`, the
valley walls below Trollveggen), and *Geilo-Ål*, x 127 000-137 000,
y 6 735 000-6 745 000 (tile `6701_3`, the north-west corner of increment 33's
section of Hallingdal).

1. `slope_stats.py`: the share of nodes at or above each angle by four slope
   measures, and Horn's slope with Gaussian noise added to z. Output, word for
   word (the Geilo-Ål line covers 33's whole section, x 127 000-150 000,
   y 6 720 000-6 745 000):

   ```
   romsdalen: 1001 x 1001 nodes, z 53..1797 m
     horn    15: 71.8% 20: 60.3% 25: 48.6% 30: 37.0% 35: 25.6% 40: 16.6% 45: 11.1%  (24.0 ns/node)
     central 15: 72.0% 20: 60.5% 25: 48.8% 30: 37.2% 35: 25.7% 40: 16.6% 45: 11.2%
     cellmax 15: 83.8% 20: 74.3% 25: 63.1% 30: 51.3% 35: 38.8% 40: 27.5% 45: 20.2%
     horn_s  15: 71.6% 20: 60.0% 25: 48.2% 30: 36.6% 35: 25.5% 40: 16.5% 45: 10.8%
     noise 0.1 m: 15: 71.8% 20: 60.3% 25: 48.6% 30: 37.0% 35: 25.6% 40: 16.6% 45: 11.1%  mean |change| 0.16 deg
     noise 0.3 m: 15: 71.8% 20: 60.4% 25: 48.6% 30: 37.1% 35: 25.6% 40: 16.6% 45: 11.1%  mean |change| 0.47 deg
     noise 1.0 m: 15: 72.2% 20: 60.7% 25: 48.9% 30: 37.3% 35: 25.9% 40: 16.8% 45: 11.2%  mean |change| 1.57 deg
     ground under 5 deg: 9.2%; p50/p90/p99 horn 24.4/46.3/66.2 deg
   geilo-al: 2501 x 2301 nodes, z 414..1243 m
     horn    15: 18.8% 20:  9.3% 25:  4.6% 30:  2.2% 35:  1.0% 40:  0.5% 45:  0.2%  (21.8 ns/node)
     central 15: 19.6% 20:  9.8% 25:  4.8% 30:  2.3% 35:  1.0% 40:  0.5% 45:  0.3%
     cellmax 15: 38.6% 20: 21.7% 25: 11.7% 30:  6.1% 35:  3.0% 40:  1.5% 45:  0.8%
     horn_s  15: 17.2% 20:  8.4% 25:  4.0% 30:  1.9% 35:  0.8% 40:  0.4% 45:  0.2%
     noise 0.1 m: 15: 18.8% 20:  9.3% 25:  4.6% 30:  2.2% 35:  1.0% 40:  0.5% 45:  0.2%  mean |change| 0.19 deg
     noise 0.3 m: 15: 19.0% 20:  9.4% 25:  4.6% 30:  2.2% 35:  1.0% 40:  0.5% 45:  0.2%  mean |change| 0.58 deg
     noise 1.0 m: 15: 20.5% 20: 10.0% 25:  4.9% 30:  2.3% 35:  1.0% 40:  0.5% 45:  0.3%  mean |change| 1.94 deg
     ground under 5 deg: 31.7%; p50/p90/p99 horn 7.8/19.5/34.7 deg
   ```

   Horn's differences are linear in z, so its slope noise scales with
   sigma / spacing: 1.0 m of noise on the 10 m grid stands for 0.1 m of noise
   on a 1 m grid (DTM1), which the probe did not have.

2. `greedy_sim.py`: greedy insertion to a slope-dependent tolerance, simulated
   in Python (the rules are in its docstring and in section 3 of the design),
   N = 2 m, F = 10 m. Commands: `python greedy_sim.py --size 1001`,
   `python greedy_sim.py --size 1001 --start 25 --end 35 --rules
   node,tri,tri+halo`, `python greedy_sim.py --case geilo-al --size 1001`
   (the 1001-node window from the corner of 33's section, the Geilo-Ål window
   above). Output, word for word:

   ```
   romsdalen 1001x1001 nodes; N 2.0 F 10.0 ramp 30.0..30.0 deg; steep (>= end) 37.0%
     uniform-F  triangles    56193 rounds  23 nodes over t(s(n)) 17.707% lazy-asked  50.5%  (8 s)
     uniform-N  triangles   412641 rounds  26 nodes over t(s(n))  0.000% lazy-asked  19.7%  (30 s)
     node       triangles   278199 rounds  24 nodes over t(s(n))  0.000% lazy-asked  23.1%  (20 s)
     tri        triangles   321685 rounds  25 nodes over t(s(n))  0.000% lazy-asked  21.4%  (23 s)
     tri+halo   triangles   344821 rounds  25 nodes over t(s(n))  0.000% lazy-asked  20.8%  (25 s)
     normal     triangles   271457 rounds  31 nodes over t(s(n))  1.019% lazy-asked  23.3%  (30 s)
   romsdalen 1001x1001 nodes; N 2.0 F 10.0 ramp 25.0..35.0 deg; steep (>= end) 25.6%
     node       triangles   253808 rounds  23 nodes over t(s(n))  0.000% lazy-asked  23.9%  (18 s)
     tri        triangles   302249 rounds  26 nodes over t(s(n))  0.000% lazy-asked  22.3%  (24 s)
     tri+halo   triangles   325719 rounds  25 nodes over t(s(n))  0.000% lazy-asked  21.5%  (24 s)
   geilo-al 1001x1001 nodes; N 2.0 F 10.0 ramp 30.0..30.0 deg; steep (>= end) 4.8%
     uniform-F  triangles    16460 rounds  21 nodes over t(s(n))  2.357% lazy-asked  60.6%  (6 s)
     uniform-N  triangles   197623 rounds  24 nodes over t(s(n))  0.000% lazy-asked  22.9%  (15 s)
     node       triangles    42174 rounds  19 nodes over t(s(n))  0.000% lazy-asked  40.1%  (6 s)
     tri        triangles    59364 rounds  23 nodes over t(s(n))  0.000% lazy-asked  35.1%  (8 s)
     tri+halo   triangles    73792 rounds  23 nodes over t(s(n))  0.000% lazy-asked  31.4%  (9 s)
     normal     triangles    39019 rounds  25 nodes over t(s(n))  0.448% lazy-asked  40.6%  (8 s)
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

4. `fixture_figures.py`: the figures of section 9's fixtures. Command:
   `python fixture_figures.py PYTHON RASPUTIN_PY PKG`. Output, word for word:

   ```
   V1: 65 x 65 nodes, 10 x 10 m, z -3.0..204.4 m
     classes: >=15: 38.5% >=25: 36.5% >=30: 35.4% >=35: 35.4% >=40: 24.3%
     nodes held to N=2 by step 30: 1495 of 4225; by ramp 25..35, below F: 1543
     rasputin uniform 10: 28 triangles
     rasputin uniform 2: 796 triangles
     simulation uniform-F 32 triangles, 7 rounds, nodes over 9.964%
     simulation uniform-N 1142 triangles, 21 rounds, nodes over 0.000%
     simulation node      408 triangles, 11 rounds, nodes over 0.000%
     simulation tri+halo  594 triangles, 11 rounds, nodes over 0.000%
   V2: 33 x 41 nodes, 10 x 5 m, z -3.0..169.2 m
     classes: >=15: 51.2% >=25: 48.1% >=30: 46.3% >=35: 46.3% >=40: 31.1%
     nodes held to N=2 by step 30: 627 of 1353; by ramp 25..35, below F: 651
     rasputin uniform 10: 6 triangles
     rasputin uniform 2: 112 triangles
   ```

5. `estimate.py`: the estimates of sections 6 and 10 from items 2 and 3, and
   the whole tile's shares and NumPy's cost of Horn per node. Output, word for
   word:

   ```
   romsdalen step 30    node     position 0.623  rasputin estimate    215227
   romsdalen step 30    tri      position 0.745  rasputin estimate    249758
   romsdalen step 30    tri+halo position 0.810  rasputin estimate    268130
   romsdalen step 30    normal   position 0.604  rasputin estimate    209873
   romsdalen ramp 25-35 node     position 0.554  rasputin estimate    195858
   romsdalen ramp 25-35 tri      position 0.690  rasputin estimate    234324
   romsdalen ramp 25-35 tri+halo position 0.756  rasputin estimate    252962
   geilo-al  step 30    node     position 0.142  rasputin estimate     29420
   geilo-al  step 30    tri      position 0.237  rasputin estimate     42066
   geilo-al  step 30    tri+halo position 0.316  rasputin estimate     52681
   geilo-al  step 30    normal   position 0.125  rasputin estimate     27099
   tile 6901_3: 5051 x 5051 nodes; NumPy Horn 15.5 ns/node, one thread
     shares: >=25: 31.0% >=30: 22.3% >=35: 14.6%
     node rule, ramp 25-35: 1367859 to 3554273 triangles
   ```
