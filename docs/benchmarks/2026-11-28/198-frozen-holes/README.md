# Increment 198 `--frozen-holes`: acceptance, untimed (2026-11-28)

Head `5e36c9389` (branch of increment 198) against master `949bc7d66`, design
`docs/increments/198-frozen-holes.md`, section 8. Machine: Apple M1 Max, 10 cores,
macOS 27.0.1, `_core` hardening "libc++ fast", on AC power (`pmset -g batt`: "Now drawing
from 'AC Power'", 100 %, charged).

**This run is untimed.** No speed gate is set for the flag (design section 7), and another
agent's untimed run shared the machine, so the wall times in `result.json` are not
benchmark figures. Section 8 also asks that, with the flag off, times stay within the
benchmark's noise. That check is not done here; the head runs the same code as master
when the flag is off (byte equality below), so no change is expected.

## Method

`compare.py` (this folder) runs every case once with master and once with the head, no
`--frozen-holes`. Each run writes a binary `.vtk`, and a second run writes a `.gpkg`, each
with `--record`. Both trees write to the same output path, because the `.gpkg` stores the
command line, output paths included. Each file is compared by sha256. The record JSON is
compared with its `command` entry dropped. To show the outputs are reproducible at all,
master was also run twice on `fx-square` and `tile`: equal.

Master's Python is `git archive 949bc7d66 src_python`, with the head's built `_core` copied
beside it. The increment changes no C++ (`git diff --stat 949bc7d66 5e36c9389` touches
only `cli.py`, `run_record.py`, tests and the design). `run_tree.py` runs it past the
venv's editable finder. The head runs as the venv's `rasputin`, after
`tools/rebuild_core.py`.

Cases: the design's fixture (section 5, built by `test_cli_mesh_frozen_holes.Inputs`:
a 101 x 101 node 10 m DEM, z = (x - 500)^2 / 1,000 m), and every case of
`docs/benchmarks/quick/cases.toml` (`19e410268`) with its thread counts. `tile` and
`quarter` are also the 1 m acceptance benchmark's two domains (`tools/bench.py`'s
defaults). None of the bench cases' domains has a hole (design section 7).

## Off: byte-equal to master

| case | `.vtk` sha256 (first 12) | `.gpkg` sha256 (first 12) | `.vtk` equal | `.gpkg` equal | records equal |
|---|---|---|---|---|---|
| fx-square (tolerance 0.1) | 8537a6a21cc7 | fbb0ae051613 | yes | yes | yes |
| fx-two-holes | 2989953fc144 | 54e05c02b436 | yes | yes | yes |
| fx-off-grid | 16efe5773751 | e96ca8a22c4d | yes | yes | yes |
| fx-no-hole | 7db0661ec54b | c3934c626b31 | yes | yes | yes |
| fx-square-budget (`--triangles 3000`) | 02d87a4438a9 | c1a376aece00 | yes | yes | yes |
| fx-square-line-through (x = 500.0004) | cf9645f75a8d | 406709d19460 | yes | yes | yes |
| fx-square-line-between (x = 512.3) | 24bb37d4901a | bbf215b5a09e | yes | yes | yes |
| fx-nodata (nodata at two hole corners) | 22a528a608ac | b9a5e2f1be74 | yes | yes | yes |
| tile, default threads | f5261521ad9d | 65c82803a986 | yes | yes | yes |
| tile, 1 thread | f5261521ad9d | 9193b16e7d4a | yes | yes | yes |
| quarter | de5bde173c21 | 2fbc05222346 | yes | yes | yes |
| numedalslagen | 260b7191c76e | cf18761eb624 | yes | yes | yes |
| lagan | 52aac9b27c10 | 613160fa9ab1 | yes | yes | yes |
| geilo-al-ramp | c9d7438fd63e | 40a8986bb08d | yes | yes | yes |
| romsdal-slope | e627393879cb | 75a1261c0a83 | yes | yes | yes |
| velhas | 9f9474cf0616 | efd3a8029c42 | yes | yes | yes |

The two tile `.gpkg` hashes differ from each other only because the stored command differs
(`--threads 1`). 16 of 16 cases are equal in all three files.

## On: the section-8 figures on the fixture (head only)

Each DEM node's error is recomputed from the written `.vtk`: its distance in z to the plane
of the triangle holding it. "On the rings" means within 1e-6 m of a hole edge. The dense
figure samples each written hole edge every 1 mm, the chord between its two written ends
against the bilinear DEM.

| run | hole vertices in / out | farthest move (m) | nodes off rings | largest error off rings (m) | nodes on rings | largest error on ring nodes (m) | dense 1 mm (m) | recorded `frozen_edges_max_error_m` | `frozen_holes` | `hole_vertices_snapped_max_m` |
|---|---|---|---|---|---|---|---|---|---|---|
| square, tolerance 0.1 | 10 / 10 | 0.0004 | 9,760 | 0.0999999829 | 80 | 0.6000000238 | 0.6000000238 | 0.6000000238 | 1 | 0.0004 |
| two holes (square, diamond) | 14 / 14 | 0.0004 | 9,659 | 0.0999999829 | 88 | 0.8999999762 | 0.8999999762 | 0.8999999762 | 2 | 0.0004 |
| off-grid hole | 4 / 4 | 0 | 9,802 | 0.0999999829 | 38 | 9.0499997139 | 9.0499997139 | 9.0499997139 | 1 | 0 |
| square, `--triangles 3000` | 10 / 10 | 0.0004 | 9,760 | 0.2000122070 | 80 | 0.6000000238 | 0.6000000238 | 0.6000000238 | 1 | 0.0004 |

Every hole vertex is kept, and none is added (out equals in). Off the rings every DEM node
is within the tolerance 0.1 m (the budget run has no tolerance; its 0.2000122 is its
`max_error_m`). The recorded figure equals the dense 1 mm sample to 10 decimals, and equals
the design's hand values 0.600, 0.900 and 9.05 to float32 rounding, so it is exact.

Mesh quality, flag off against on (reported, not gated: the option is off by default):

| fixture | triangles off / on | worst angle off / on (degrees) | median smallest angle off / on | max vertex degree off / on |
|---|---|---|---|---|
| square | 5,438 / 5,462 | 1.52 / 1.52 | 36.9 / 36.9 | 11 / 11 |
| two holes | 5,159 / 5,209 | 0.67 / 0.67 | 36.9 / 36.9 | 12 / 12 |
| off-grid hole | 5,686 / 5,670 | 1.33 / 0.70 | 36.9 / 36.9 | 10 / 16 |
| square, `--triangles 3000` | 2,702 / 2,232 | 1.74 / 1.62 | 33.7 / 35.8 | 9 / 9 |

The off-grid hole loses some quality next to the kept edges (worst angle 0.70, a vertex of
degree 16). This is the cost Triangle's `-Y` documents (design section 1.1).

## Regenerate

```bash
M=../rasputin_scratch/perf-198-master; mkdir -p $M
git archive 949bc7d66 src_python | tar -x -C $M
python3 tools/rebuild_core.py
cp .venv/lib/python3.14/site-packages/tin_engine/_core.*.so $M/src_python/tin_engine/
.venv/bin/python docs/benchmarks/2026-11-28/198-frozen-holes/compare.py \
    $M/src_python ../rasputin_scratch/perf-198-work docs/benchmarks/2026-11-28/198-frozen-holes/result.json
```

The bench cases need the data folder (`$RASPUTIN_DATA`, default `../rasputin_data`). Meshes
stay in the work folder, out of the repository. `result.json` holds every row, with each
run's exit code and wall time.

Verdict: **ACCEPTED** (flag off byte-equal on 16 of 16 cases; flag-on figures as designed).
Timing with the flag off is not measured in this untimed run.
