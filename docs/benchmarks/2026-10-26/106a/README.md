# 106a's polygon repair, clipped to the mesh's box: timing (@perf, 2026-10-26)

The fix timed: commit `993ed3a2` (`worktree-demo-fixes`, 106a on master
`79535812`) repairs a self-crossing land-cover or reservoir polygon with
`shapely.make_valid` only on its part inside the mesh points' box widened by
100 m (`_repaired`, `src_python/tin_engine/feature_input.py`). The run before
it, `ee1dda88`, repaired each polygon whole: Kultsjön's "land cover fractions"
step took 45.20 s against master's 10.43 s
(`docs/benchmarks/2026-10-26/106c/README.md` on `worktree-demo-fixes-c`,
`a9366418`).

## Method

Apple M1 Max (8 P + 2 E cores), 32 GiB, macOS 27.0.1, Python 3.14.7,
shapely 2.2.0 (GEOS 3.14.1). On AC ("Now drawing from 'AC Power'", read
before each run). Release `_core`, `--hardening on` (libc++ fast), built by
`bench.py`'s `build()` into each tree's `build-bench/pkg`: master from a
`git archive` of `79535812` in `rasputin_scratch/perf-106a-base/`, the head
in the worktree. Both `.so` files are the same bytes (sha256 `a353ac1e…`, as
in the 106c run); 106a has no C++, so no sanitizer run. Each run goes
through `bench.py _child --threads 0` (all 10 threads), master then head,
back to back, **one run each** (time-boxed; no median). `run.sh` holds the
commands; the times are the run's own phase clock (the `BENCH` line in each
`.log`). Desktop applications were open; nothing else of mine ran.

Kultsjön: `--dem glo30 --out-crs EPSG:3006 --domain
~/rasputin-output/opal-se2-2026-10-21/catchments/kultsjon.geojson
--land-cover clc2018 --triangles 6000`, no `--catchments`.
Nedre Sjodalsvatn (NVE gauge 2.13.0, 480.61 km², 1,663 vertices):
`--dem dtm10-norway --domain nedre-sjodalsvatn.geojson --land-cover
ar50-innlandet --reservoirs nve-reservoirs --triangles 6000`, the gap
analysis's command. Its outline was taken from the gauge catchments
(`~/.cache/rasputin/nve-gauge-catchments/features.gpkg`) and written as its
single Polygon (stored as a one-part MultiPolygon, which `--domain` refuses).

## Results

| catchment | build | land cover fractions (s) | run total (s) | polygons repaired | triangles | outputs |
|---|---|---:|---:|---|---:|---|
| Kultsjön | master `79535812` | 10.49 | 11.53 | (none) | 5,720 | |
| Kultsjön | head `993ed3a2` | **5.72** | 6.70 | 43 of 1,170 | 5,720 | equal to master's |
| Kultsjön | `ee1dda88` (106c run, earlier today) | 45.20 | 46.24 | 43 | 5,720 | |
| Nedre Sjodalsvatn | master | crashes: GEOS `TopologyException` ("side location conflict") at 2.6 s | | | | none |
| Nedre Sjodalsvatn | head | **2.49** | 3.20 | 32 of 2,843 | 5,940 | finishes |

"Equal to master's" (`same_output.py`): the tables `triangles` (5,720 rows),
`cells`, `land_cover_fractions` (10,298 rows), `land_cover_classes`, `lines`
and `layer_styles` are identical row for row, byte for byte; the `run` table
differs only in its `command` row (the package path and file names) and the
new `land_cover_repaired` row (43). Mesh quality is therefore master's.

Head against master, Kultsjön: the land-cover step is 4.77 s faster (0.55 of
master's). The 106c run's 34.8 s added by repairing whole polygons is gone.
Why the head beats master is thought to be the clip: the 43 polygons the
cut then meets are cut down to the box, so each overlay works on fewer
vertices. Not profiled. One run each, so a few tenths of a second are noise.

Nedre Sjodalsvatn finishes with 5,940 triangles, as the gap analysis found
with AR50 repaired by hand; master still throws.

## Files

- `run.sh`: the commands (`sh run.sh kultsjon`, `sh run.sh sjodal`).
- `same_output.py`: the table comparison, with its command in its docstring.
- `*.log`: each run's stderr, with the `BENCH` phase line and `time -p`.
- `*.stats.md`: each finished run's `--stats` output.
- `nedre-sjodalsvatn.geojson`: the domain used.

The meshes were written to `rasputin_scratch/perf-106a-out/` and removed;
`run.sh` regenerates them.

## The corrected repair, `9f33845e` (later the same day)

`9f33845e` replaces `make_valid` on the clipped polygon: it cuts the
polygon's rings to the same widened box, splits them where they cross, and
keeps the faces the odd-count rule (the one `make_valid` uses) puts inside.
Same machine, method and `run.sh`, unchanged; on AC before each run. Master
`79535812` was rebuilt from a fresh `git archive` (the earlier copy was
removed); both `.so` files are again sha256 `a353ac1e…`, so the C++ is the
same bytes and no sanitizer run is owed. Master then head, back to back, one
run each. The files are in `9f33845e/` (the head runs are `*-head.*`).

| catchment | build | land cover fractions (s) | run total (s) | polygons repaired | triangles | outputs |
|---|---|---:|---:|---|---:|---|
| Kultsjön | master `79535812` | 10.51 | 11.51 | (none) | 5,720 | |
| Kultsjön | head `9f33845e` | **5.61** | 6.58 | 43 | 5,720 | equal to master's |
| Kultsjön | `993ed3a2` (above) | 5.72 | 6.70 | 43 | 5,720 | equal to master's |
| Nedre Sjodalsvatn | master | crashes: GEOS `TopologyException` ("side location conflict") at 2.6 s | | | | none |
| Nedre Sjodalsvatn | head `9f33845e` | **2.50** | 3.20 | 32 | 5,940 | finishes |

Outputs (`same_output.py`, Kultsjön): `triangles`, `cells`,
`land_cover_fractions`, `land_cover_classes`, `lines` and `layer_styles` are
identical to master's row for row. In the `run` table the only value that
differs is the `command` row (the package path and output file names, which
name the build); the rest is the new `land_cover_repaired` row (43), which
shifts every later row's number by one, so `same_output.py` reports many
"rows differ" there: compared by name, every other run value is equal. Mesh
quality is therefore master's. Nedre Sjodalsvatn's `--stats` matches the
`993ed3a2` run's in every count; only phase times differ.

Head against master, Kultsjön: the land-cover step is 4.90 s faster (0.53 of
master's); against `993ed3a2` it is 0.11 s faster, within one run's noise.
