# A real piece of the São Francisco basin: triangle counts against tolerance (@perf, 2026-10-01)

ROADMAP "Order of work from 2026-09-30", item 2.1. Measured on AC power, Apple
M1 Max (8 performance + 2 efficiency cores, 32 GiB), macOS 27.0, a Release
`_core` built by `tools/bench.py`'s `build()` from commit `838796c`
(`build_and_machine.json` has the `.so`'s sha256 and the flags). No production
code was changed.

**The DEM is Copernicus GLO-30, not ANADEM.** ANADEM's host,
`metadados.snirh.gov.br`, answered HTTP 403 to every request from this machine
between 23:41 on 30 September and 02:11 on 1 October (UTC), including the
geonetwork front page; the earlier probes of `15-dem-mosaic.md` (B1-B9) reached
it. GLO-30 is ANADEM's parent (ANADEM is GLO-30 with the vegetation bias
removed), on AWS, readable by range requests. So every figure here is for a
**surface** model (canopy and roofs included). That ANADEM gives fewer
triangles at small tolerances is thought to be so (a terrain model is
smoother under forest), not measured. The scripts take `anadem` as the source
the moment the host answers (`prep_dem.py`, tile 23K; see "Not done").

## The piece, and why it

**BHO ottobasin 76949**: the upper and middle Rio das Velhas around Belo
Horizonte, 1,163 elementary catchments of BHO 2017 (5k) dissolved into one
polygon with 7,310 vertices and no holes, **11,667.6 km²** (geodesic; BHO's
own `NUAREACONT` sum is 11,667.561 km²). Chosen because:

- it lies wholly inside one ANADEM tile (23K) and one UTM zone (23S), so the
  same piece can be re-run on ANADEM without a mosaic;
- it is real terrain of the kind that drives triangle counts up: the
  Quadrilátero Ferrífero and the southern Espinhaço, heights 528-1,926 m, plus
  the Belo Horizonte urban area and lower ground to the north;
- a first 20 m run took 0.84 s, so a piece larger than "a few thousand km²"
  cost nothing in time; at 11.7 k km² it is 13.0 M source nodes, and the 1 m
  mesh still runs in under 10 s.

Because one steep piece cannot stand for the basin, the extrapolation uses a
second measurement: **200 random 10 km boxes over the whole basin**.

## Method

1. **Outline** (`fetch_bho.py outline`): ANA's ArcGIS REST service for BHO
   2017 5k drainage areas, `COBACIA LIKE '76949%'`, geometry in EPSG:4674,
   unioned with shapely. The basin outline is BHO level-2 ottobasin 76 (50k),
   object 60 of `SNIRH2016/Divisao_de_bacias` layer 6.
2. **DEM window** (`prep_dem.py fetch`): only the 1024² COG blocks under the
   domain's projected box (plus 4 target cells, mapped back to longitude and
   latitude and grown by 2 source cells) are range-read and stitched on GLO-30's
   1" lattice. Piece: 7,197 × 4,376 nodes from 45 blocks of 6 tiles, 80 s.
3. **Projection and resampling** (`prep_dem.py resample`), Q6's ruling done
   by hand: a square 30 m node grid in SIRGAS 2000 / UTM 23S (EPSG:31983),
   nodes on multiples of 30 m, z bilinear in the source's own
   (longitude, latitude) index space. Piece: 7,347 × 4,208 nodes, 0.99 s on 8
   threads. Written as a node-registered float32 GeoTIFF, which the existing
   reader accepts. PROJ takes WGS 84 to SIRGAS 2000 by "SIRGAS 2000 to WGS 84
   (1)", which moves (44°W, 19°S) by exactly 0.
4. **Mesh**: the existing pipeline, unchanged:
   `rasputin mesh --dem GRID --domain OUTLINE --tolerance T --binary`, the
   outline reprojected by increment 15b. Each run goes through
   `tools/bench.py _child` (the Release package; refine timed inside the
   process) under `/usr/bin/time -l`. Threads: the CLI default, 10
   (hardware concurrency, recorded in each `--stats` file).
5. **Quality**: one extra ASCII run per tolerance, judged by `tools/bench.py`'s
   `quality()` (worst angle, max degree, the constrained Delaunay check, exact
   where the float test is ambiguous).
6. **The final check of Q6, measured, not built**: the mesh is read linearly
   at every source (GLO-30) node inside the domain, projected, and compared with
   the source value; a node off by more than the tolerance is one the final
   check would find. Split into the **interior** and the **boundary strip**
   (within one grid spacing, 30 m, of the domain's edge; see "Surprises").
   Its **control** reads the mesh at the grid nodes it was refined against,
   where refine guarantees the tolerance.
7. **Basin sample** (`sample_boxes.py`): 200 box centres drawn uniformly by
   area over the basin (Albers equal-area, seed 1, 584 draws), a box kept only
   if wholly inside the basin; each box is a 10 km square in its centre's UTM
   zone (EPSG:31960 + zone), prepared and meshed exactly as the piece, once
   per tolerance. The basin figure is the mean density times the basin's
   area, with a bootstrap 95 % interval of the mean. A box centre is at least
   half a box (5 km) from the basin's edge, so the sample under-weights a band
   that wide along it, 5.6 % of the basin's area. A box's area is taken as
   100 km² in its UTM zone; its geodesic area differs by less than 0.2 %
   (0.9992-1.0016 over the 200 boxes). The interval is the sampling error of
   the mean only. It does not cover the unsampled edge band, the use of a
   surface model for a terrain model, or the box method's own bias (step 8).
8. **Check of the box method**: 30 random boxes inside the piece, against the
   piece meshed whole.

## Results

All tables are `analysis.md`, printed by `analyse.py` from the JSON in
`runs/`. Areas are geodesic on GRS80.

### The piece, 11,667.6 km² (medians of 5 timed runs)

| tolerance | triangles | vertices | vertices / grid nodes | triangles / km² | refine | process wall | max RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 10,431,955 | 5,220,159 | 40.3 % | 894.1 | 6.33 s | 9.67 s | 3.60 GiB |
| 2 m | 5,736,406 | 2,872,292 | 22.2 % | 491.7 | 3.36 s | 5.54 s | 2.18 GiB |
| 5 m | 1,993,312 | 1,000,480 | 7.7 % | 170.8 | 1.15 s | 2.38 s | 1.01 GiB |
| 10 m | 807,337 | 407,376 | 3.1 % | 69.2 | 0.50 s | 1.44 s | 0.79 GiB |
| 20 m | 314,507 | 160,916 | 1.2 % | 27.0 | 0.24 s | 1.03 s | 0.70 GiB |
| 50 m | 82,517 | 44,914 | 0.3 % | 7.1 | 0.10 s | 0.83 s | 0.62 GiB |

Threads 10; AC before and after every tolerance block. "Process wall" is the
whole child process (interpreter start, decode, refine, write, `--stats`).
Quality at every tolerance: worst angle 0.107-0.113°, max degree 13-14,
**0 constrained-Delaunay violations**, within tolerance, and the control at
**0 of 12,957,257** grid nodes over tolerance. The `--stats` phase timer puts
4.89-4.91 s of the 1 m refine's 6.32-6.37 s, 77 %, in "split + flip
(serial)" (five runs, `runs/velhas76949_glo30/logs/velhas76949_t1_run*.stats.md`).

### The basin, 635,194.5 km² (BHO level 2; OAS gives 636,920 km²)

From the 200 random boxes. "Refine, linear" is the piece's refine time per
triangle at 10 threads times the basin's triangles. "Memory floor" is the
piece's max-RSS slope over the sweep (least squares: 0.55 GiB + 310 B per
triangle, max residual 121 MiB) times the basin's triangles, plus Q9's 8.6 GiB
canvas; one process, the whole basin.

| tolerance | triangles / km², mean (95 %) | median | p10-p90 | **basin triangles** (95 %) | **basin vertices** | refine, linear | memory floor |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 m | 373.4 (338.1-410.9) | 315.6 | 101.3-757.0 | **237 M** (215-261 M) | 119 M | 2.4 min | 77 GiB |
| 2 m | 178.6 (158.4-199.6) | 138.5 | 40.3-391.5 | **114 M** (101-127 M) | 57 M | 1.1 min | 41 GiB |
| 5 m | 56.6 (49.2-64.6) | 39.6 | 8.4-139.3 | **36 M** (31-41 M) | 18 M | 0.3 min | 19 GiB |
| 10 m | 21.7 (18.5-24.9) | 13.6 | 2.2-53.4 | **14 M** (12-16 M) | 6.9 M | 0.1 min | 13 GiB |
| 20 m | 7.7 (6.5-9.1) | 4.3 | 0.1-19.9 | **4.9 M** (4.1-5.8 M) | 2.5 M | 0.1 min | 10 GiB |
| 50 m | 1.6 (1.3-2.0) | 0.7 | 0.0-4.6 | **1.0 M** (0.8-1.2 M) | 0.5 M | < 0.1 min | 9 GiB |

Extrapolating from the piece alone would give 568 M triangles at 1 m and 44 M
at 10 m: the piece's density is 2.4 times the sample's mean at 1 m, 3.2 times
at 10 m and 4.4 times at 50 m, because the piece is steeper
than the basin (median relief inside a box: 285 m in the piece, 192 m over the
basin).

**The box method holds on the piece**: at every tolerance the piece's own
density lies inside the 95 % interval of its 30 boxes (1 m: 894.1 against
934.2, 831.3-1,037.5; 10 m: 69.2 against 73.7, 60.7-86.9). Its boxes still
read 4.5-6.5 % above the whole-piece mesh at 1-20 m and 8.2 % below at 50 m.

### The final check of Q6, measured

Source nodes off the mesh by more than the tolerance after meshing the
resampled grid, interior only, over the basin's boxes:

| tolerance | share of source nodes | basin source nodes over tol | as a share of the basin's vertices |
|---:|---:|---:|---:|
| 1 m | 7.71 % | 52 M | 44 % |
| 2 m | 2.50 % | 17 M | 30 % |
| 5 m | 0.39 % | 2.6 M | 15 % |
| 10 m | 0.08 % | 0.53 M | 8 % |
| 20 m | 0.02 % | 0.12 M | 5 % |

On the piece: 22.30 % of interior source nodes at 1 m (max 27.1 m), 0.31 % at
10 m. A count of nodes over tolerance is not a count of insertions (one
insertion can fix its neighbours); it bounds the first pass. The resampled
grid itself, read at the piece's 13.0 M source nodes, is off by median
0.36 m, p99 3.14 m, max 27.8 m (`data_meta/*_check.json`, as B6 measured on
the Espinhaço window).

### Refine thread scaling, the piece at 1 m (medians of 3, interleaved)

| threads | refine | speed-up | process wall |
|---:|---:|---:|---:|
| 1 | 12.95 s | 1.00 | 14.12 s |
| 2 | 9.12 s | 1.42 | 10.31 s |
| 4 | 7.30 s | 1.77 | 8.49 s |
| 6 | 6.69 s | 1.93 | 7.87 s |
| 8 | 6.36 s | 2.03 | 7.55 s |
| 10 | 6.32 s | 2.05 | 7.50 s |

AC before and after. No baseline on this input exists; this is the first.

## Surprises

1. **Memory, not time, limits the basin.** Refine is fast (6.3 s for 10.4 M
   triangles; linear extrapolation 2.4 min for the basin at 1 m), but one
   process holding the basin needs at least ~77 GiB at 1 m and ~41 GiB at 2 m.
   At 5 m the floor is ~19 GiB, under the 32 GiB Mac, but a floor cannot show
   that it fits. The piece's own intercept (0.55 GiB for a 30.9 M-node grid,
   about 17 B per grid node once the ~66 MiB interpreter is taken off) is
   about four times the 4 B per node that the canvas term counts. At that rate
   the basin's 2.31 G-node canvas alone would be ~36 GiB. Whether 5 m fits is
   not measured.
2. **The final check is not a touch-up at small tolerances.** At 1 m it would
   find 52 M basin source nodes, 44 % of the mesh's vertex count; at 2 m 30 %;
   at 10 m 8 %. The resampling error alone (median 0.36 m, p99 3.1 m) is the
   size of a 1-2 m tolerance.
3. **The tolerance is unchecked in a strip along the domain's edge.**
   Increment 16 guarantees it "at every valid node inside the domain"; a
   sliver along an edge that holds no grid node is never checked. On a 10 km
   box with 4 corners, source nodes within 30 m of the edge were off by up to
   541.6 m at a 20 m tolerance (box 127; the interior's max there was
   27.3 m); point location was cross-checked by brute force. On the dense
   BHO outline of the piece the strip's worst is 17.9-38.5 m over the six
   tolerances, and 0-30 % of its 26,294 strip nodes are over tolerance. The
   final check of Q6 would see these nodes. For a projected DEM meshed
   directly (Norway) every DEM node inside the domain is a grid node, so the
   guarantee at nodes holds. The sliver between the edge and the outermost
   nodes is still unchecked. Measured against the DEM's own surface between
   nodes (bilinear, increment 16 R0), it is as far off as here. The cause is
   that a triangle whose closed node set is empty has error 0 and converges
   (`include/terrain/refinement/scan.hpp`), so box 127 at 20 m keeps its 4
   corners as its only boundary vertices.
4. **The basin lies in four ANADEM tiles, not six** (B5, bounding-box
   estimate): by BHO's outline and B3's band edges, 23K 158,153 km², 23L
   308,944 km², 24L 143,525 km², 24M 24,573 km²; nothing west of 47.65°W, so
   no 22K or 22L.
5. **At 1 m a 30 m surface model barely compresses**: the piece keeps 40 % of
   its grid nodes; over the basin, 119 M vertices against 688 M GLO-30 nodes
   (the boxes' mean of 1,083.9 source nodes per km² times the basin's area).
6. **Refine speeds up 2.05× at 10 threads** on this input; the serial
   split + flip is 77 % of refine time at 10 threads (the `--stats` timer;
   not a profile).

## Data and credits

The figures here are derived from **Copernicus GLO-30**, "produced using
Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space
GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all
rights reserved" (projected and resampled, so modified), and from the
**BHO 2017** catchment outlines of the Agência Nacional de Águas e Saneamento
Básico (ANA). Neither dataset is in the repository.

## Reproduce

From the repository root, the venv's python as `PY`, `D=../rasputin_data/sao_francisco_piece`,
`S=../rasputin_scratch/basin-piece`, `E=docs/benchmarks/2026-10-01/basin-piece`.
The data folder's `README.md` lists sources and licences; `data_meta/sha256.txt`
has the inputs' checksums.

```bash
$PY -c "import sys; sys.path.insert(0,'tools'); import bench; from pathlib import Path; \
print(bench.build(bench.make_runner(), Path('.').resolve()))"          # Release build-bench/pkg
$PY $E/fetch_bho.py outline $D 76949                                    # the piece's outline
curl -s -G "https://www.snirh.gov.br/arcgis/rest/services/SNIRH2016/Divisao_de_bacias/FeatureServer/6/query" \
  --data-urlencode objectIds=60 --data-urlencode outFields=* --data-urlencode outSR=4674 \
  --data-urlencode f=geojson -o $D/bho2017_level2_76_raw.geojson       # the basin's outline
O=$D/bho2017_5k_76949_outline_epsg4674.geojson
$PY $E/prep_dem.py fetch glo30 $O $D && $PY $E/prep_dem.py resample glo30 $O $D \
  && $PY $E/prep_dem.py check glo30 $O $D                               # window, grid, B6-style check
G=$D/derived/bho2017_5k_76949_glo30_epsg31983_30m.tif; W=$D/bho2017_5k_76949_glo30_window_epsg4326.tif
$PY $E/run_sweep.py velhas76949 $G $O $W $E/runs/velhas76949_glo30 $S --repeats 5   # ~12 min
$PY $E/scaling.py $G $O $E/runs/velhas76949_glo30 $S                    # ~6 min
$PY $E/sample_boxes.py glo30 $D/bho2017_level2_76_raw.geojson $D $E/runs/basin_boxes_glo30 $S --n 200   # ~20 min
$PY $E/sample_boxes.py glo30 $O $D/piece_boxes $E/runs/piece_boxes_glo30 $S --n 30
$PY $E/plant_check.py $G $O $W $S 20 15                                 # the control must fail when shifted
$PY $E/analyse.py $E/runs/velhas76949_glo30/results.json $O $D/bho2017_level2_76_raw.geojson \
  $E/runs/basin_boxes_glo30/boxes.json $E/runs/piece_boxes_glo30/boxes.json \
  $E/runs/velhas76949_glo30/scaling.json > $E/analysis.md
# Surprises 4: km² per ANADEM tile (6° zones, B3's band edges) and the west edge
$PY -c "
import json, sys, pyproj, shapely
from shapely.geometry import shape
g = shape(json.load(open(sys.argv[1]))['features'][0]['geometry'])
geod = pyproj.Geod(ellps='GRS80')
print('west edge', round(g.bounds[0], 3))
for z, (w, e) in {'22': (-54, -48), '23': (-48, -42), '24': (-42, -36), '25': (-36, -30)}.items():
    for b, (s, n) in {'K': (-24, -16.26), 'L': (-16.26, -8.13), 'M': (-8.13, 0)}.items():
        p = g.intersection(shapely.box(w, s, e, n))
        if not p.is_empty: print(z + b, round(abs(geod.geometry_area_perimeter(p)[0]) / 1e6), 'km2')
" $D/bho2017_level2_76_raw.geojson
# Determinism: the final runs against the superseded earlier run
$PY -c "
import json, sys
E, S = sys.argv[1] + '/runs/', sys.argv[2] + '/superseded_v1/'
for b in ('basin_boxes_glo30', 'piece_boxes_glo30'):
    old = {x['box']: x for x in json.load(open(S + b + '/boxes.json'))['boxes']}
    new = {x['box']: x for x in json.load(open(E + b + '/boxes.json'))['boxes']}
    key = lambda x, t: (x['tolerances'][t]['triangles'], x['tolerances'][t]['vertices'])
    print(b, len(old), 'boxes, differences', sum(key(o, t) != key(new[k], t) for k, o in old.items() for t in o['tolerances']))
o = json.load(open(S + 'velhas76949_glo30/results.json'))['runs']
n = json.load(open(E + 'velhas76949_glo30/results.json'))['runs']
print('piece: triangles and mesh sha256 equal at every tolerance:', all(
    a['timed'][0]['triangles'] == b['timed'][0]['triangles'] and a['quality']['mesh_sha256'] == b['quality']['mesh_sha256'] for a, b in zip(o, n)))
" $E $S
curl -s -o /dev/null -w "%{http_code}\n" -r 0-10 \
  https://metadados.snirh.gov.br/files/anadem_v1_tiles/anadem_v1_23K.tif   # 403 throughout this session
```

## Checks behind the claims

- **The final check's control can fail** (`runs/checks.txt`, `plant_check.py`):
  the piece at 20 m, as meshed, 0 grid nodes over tolerance and a max of
  19.999989197 m, equal to refine's own reported max error to 9 decimal places; the
  same mesh moved 15 m east, 22,548 over tolerance.
- **Point location** (`matplotlib.tri`) agreed with a brute-force barycentric
  search at the worst nodes of box 127 within 2e-8 m (run by hand, not
  committed).
- **Determinism**: every count and every quality-run mesh hash was identical
  in a full earlier run of the piece and both box samples
  (`runs/checks.txt`; the earlier run is kept in
  `../rasputin_scratch/basin-piece/superseded_v1/`). That earlier run is superseded only because it lacked
  the interior/strip split, and its two box samples had run concurrently
  sharing scratch mesh names (now unique per sample).
- **Script history**: `run_sweep.py`'s `final_check` read EPSG:31983 as fixed
  in its first version, which failed loudly on the first zone-24 box; it now
  reads the grid's own `ProjectedCSTypeGeoKey`. `scaling.json` was produced
  by the version before the interior/strip split, which does not touch the
  timed runs.

## Not done

- **ANADEM.** Re-run with `prep_dem.py fetch anadem ...` once the host
  answers; `fetch_anadem` reads only tile 23K, enough for the piece. The box
  sample would need tile selection by position (four tiles, overlaps of 8-87
  cells); not written.
- The final check's insertions (only its first-pass violators are counted),
  and a mesh of the whole basin in one run.
- Large outputs: meshes went to `../rasputin_scratch/basin-piece/` (the binary
  piece meshes are kept there; ASCII ones were deleted after reading). The
  BHO outline and the DEM grids stay in `../rasputin_data/sao_francisco_piece/`.

## Files

`fetch_bho.py`, `prep_dem.py`, `run_sweep.py`, `sample_boxes.py`,
`plant_check.py`, `scaling.py`, `analyse.py`: one-off measurement scripts, a
stopgap until rasputin fetches and projects its own data (15c, and the
basin's own inputs). `analysis.md`: every table. `runs/*/results.json`,
`boxes.json`, `scaling.json`: raw figures, including `pmset -g batt` before
and after each block. `runs/velhas76949_glo30/logs/`: each run's output and
`--stats`; the box samples' logs are in `logs.tar.gz`. `data_meta/`: the
window's and grid's metadata, the B6-style check and checksums.
`build_and_machine.json`: commit, build and machine.
