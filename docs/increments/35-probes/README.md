# Design probes for increment 35 (`@architect`, 2026-10-11)

Run from the `cleanup-speed` worktree at master `c7dba32e` with its own
`.venv` (Python 3.14.7, shapely 2.2.0 on GEOS 3.14.1, pyproj 3.8.0 on PROJ
9.8.1, `_core` as `@perf` built it for the profile, C++ equal to master's).
Data: `/Users/skavhaug/projects/rasputin_data` (the paths are absolute in the
scripts). The quick check's `numedalslagen` and `lagan` arguments
(`docs/benchmarks/quick/cases.toml`), the CLI's clean-up defaults (repair 1 m,
merge on, band 50 m, outline rule 5 m). Timings are single runs and rank the
sub-steps; the first three scripts ran on battery (88 to 95 %), `mesh_check.py`
on AC power. No C++ was built.

1. `clean_probe.py CASE`: builds the `cover` list as `_Tally._take` does and
   runs `_Tally._clean`'s sub-steps one by one with a timer, each candidate
   beside master's call, outputs compared as WKB. Output, the lines that matter
   (Numedalslagen / Lagan):

   ```
   cover: 7932 polygons, 5417093 vertices before the clip   /  7866 polygons, 15888613 vertices
   C1 reprojection: per polygon 4.01 s, one call 4.08 s, same bits: True        (Lagan)
   clip to region: 0.68 s (7932 calls); C2: 3077 inside, 0.49 s, same WKB: False, same point sets: False
                   0.55 s (7866 calls); C2: 3468 inside, 0.42 s, same WKB: False, same point sets: False
   coverage_clean: 1.15 s / 0.86 s
   snap_to_outline: 2.25 s; inside it {'_Outline.ring': 0.56, '_chains': 0.13, 'make_valid@710': 0.13,
                   '_loops': 0.17, 'union_all@680': 0.56, 'make_valid@676': 0.29, 'union_all@687': 0.11}
                    2.13 s; inside it {'_Outline.ring': 0.42, '_chains': 0.06, 'make_valid@710': 0.14,
                   '_loops': 0.17, 'union_all@680': 0.36, 'make_valid@676': 0.51, 'union_all@687': 0.11}
   clip to domain: 0.81 s (25 calls); C6: 0.45 s, same WKB: False, same part sets: False
                   0.70 s (24 calls); C6: 0.41 s, same WKB: False, same part sets: False
   simplify_borders: 0.68 s / 0.64 s
   _outline_lines: 0.27 s (25 calls); C5: one outline tree costs 0.002 s, so about 0.06 s saved
                   0.71 s (24 calls); C5: one outline tree costs 0.010 s, so about 0.22 s saved
   _add lines: 2570; _untouched 0.02 s, 1944 pass; intersection of the rest 0.15 s; covers 0.02 s;
               C3 covers only where not untouched 0.02 s, same: True
               4255; _untouched 0.04 s, 3380 pass; intersection of the rest 0.62 s; covers 0.76 s;
               C3 covers only where not untouched 0.75 s, same: True
   C4: of 626 not untouched, 71 match the pattern (domain first 0.00 s, line first 0.09 s, agree: True);
       their intersection is the line itself: True; 75 of the rest split into parts
       of 875 not untouched, 850 match the pattern (domain first 0.01 s, line first 0.76 s, agree: True);
       their intersection is the line itself: True; 25 of the rest split into parts
   ```

   (The C3 row's time is the cost of `covers` on the lines that are *not*
   untouched, 626 and 875 of them; the saving is the difference once C4 has
   taken most of those too: measured end to end in `mesh_check.py`.)

2. `ring_order.py numedalslagen`: why C2 fails. Of the 3077 polygons the read
   region contains properly, `intersection(g, region)` returns the same point
   set as g (`equals` True) but every ring rewritten: the shell clockwise,
   starting at the input's last vertex and running backwards (3066 match that
   rule exactly), holes counter-clockwise and, in 11 polygons with holes,
   in another order or start. Output: `{'revrot': 3066, 'other': 11}`.

3. `stage_check.py CASE`: whether the next stage washes that order out, and
   the threaded transform. Output:

   ```
   C2 after coverage_clean (1.21 s): same WKB: False      (both cases)
   C6 after simplify_borders: same polygons: False; same lines: False   (both cases)
   C1' 15888613 points: one call 3.94 s, 2 threads 2.04 s, same bits: True
   C1' 15888613 points: one call 3.94 s, 4 threads 1.06 s, same bits: True
   C1' 15888613 points: one call 3.94 s, 8 threads 0.63 s, same bits: True
   ```

4. `mesh_check.py CASE DIR`: the four changes as monkeypatches on the installed
   `tin_engine`, `rasputin mesh` run as master and as patched in one process,
   AC power, one run each:

   ```
   numedalslagen master:  10.03 s wall, features clip 6.291, clean-up 6.219, features: 22 kept (15 cut ...), 4638 outside, 0 empty
   numedalslagen patched:  8.99 s wall, features clip 5.999, clean-up 5.936, the same features: line
   numedalslagen: meshes byte-identical: True  (sha256 349d0b56...)
   lagan master:   20.64 s wall, features clip 11.744, clean-up 7.547, features: 23 kept (0 cut ...), 4193 outside, 0 empty
   lagan patched:  16.28 s wall, features clip  6.742, clean-up 5.756, the same features: line
   lagan: meshes byte-identical: True  (sha256 6c110f9e...)
   ```

   An earlier run of the same script without the batch in `_take` (3.1 (a))
   left Lagan's reprojection at 4.05 s: no single polygon has the 500 000
   vertices that start the threads, so the batch across polygons is what makes
   them bite.
