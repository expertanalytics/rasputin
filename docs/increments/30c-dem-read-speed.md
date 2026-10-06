# Increment 30c — the `decode` phase of `rasputin mesh`, made faster

Status: **built and accepted; the gate's base still to be recorded again at
the master merge `26a5d839` (section 6, "At the merge"); code review next.**
Red `bdf7b57d` (`@tester`'s pins beyond the design ruled in section 7: six
stand, one changes), green `0663efc9` (+3 net production lines,
`python3 tools/count_loc.py 6c729e97 0663efc9`; as built, section 3.3),
`@perf`'s acceptance ACCEPTED at `bc8d91eb` (section 9). Designed by
`@architect` 2026-10-06 on branch `worktree-dem-read-speed`, at `6c729e97`
(30b's approved head). Ola said "yes, build 30c" on 2026-10-06 (question 1,
section 13). Design review round 2 approved (`## Review`). 30b merged as
#200; master (with #200 to #203) was merged into this branch as `26a5d839`.

**What this is.** The third of three pull requests that remove the
bottlenecks `@perf` measured in `rasputin mesh` on two real catchments
(`docs/benchmarks/2026-10-06/bottlenecks/README.md`). Ola approved them on
2026-10-06 in this order: land cover (30a, merged as #199), the CORINE clip
(30b, `docs/increments/30b-clip-speed.md`), then reading the DEM (this one).
The ROADMAP row is 30.

**What changes for a user.** Nothing in any output: the assembled DEM, every
value in it, the seam report (`dem_seams`, the `--stats` table) and the mesh
stay byte-identical. The `decode` row gets about 40 % shorter: on AC power,
1.79 s to 1.06 s on Numedalslågen and 1.31 s to 0.82 s on Skiensvassdraget,
measured on a prototype through `rasputin mesh` (section 8). That is 0.7 s of
a 7.7 s run and 0.5 s of an 11.8 s run.

**Lean.** No mutation round (Ola's standing rule for lean rounds). No C++. Two
production files, about +3 net lines in the prototype.

## 1. Prior art: legacy and literature

### Literature

Nothing here is a method of its own. The two changes are:

- **Decoding independently compressed blocks in parallel.** A tiled TIFF
  compresses each tile on its own (TIFF 6.0, section 15, "Tiled Images"; a
  Cloud Optimized GeoTIFF is the same layout), so blocks decode in any order
  and on any thread. `decode_window` already does this on a pool of 4
  threads, each block into its own slice of one array (increment 23a-1,
  `docs/increments/23-basin-scale.md`, decided 7 and W2). This increment only
  sizes the pool to the machine. Recalled, not reread: the TIFF 6.0 section
  number.
- **Cheap tests before the expensive one** in the seam report: the same
  filter-then-test pattern 30b used (`docs/increments/30b-clip-speed.md`,
  section 1). Here the cheap tests (both values valid, gap at least the
  threshold) and the expensive one (the node inside the needed region) are
  all exact, and they are ANDed, so their order cannot change the set they
  select (section 4.2).

The facts this design relies on were read in the installed versions
(tifffile 2026.9.20, imagecodecs 2026.8.16, numpy 2.5.3, shapely 2.1.2 with
GEOS 3.13.1) or measured:

- **DTM10's tiles are LZW, no predictor, 512 × 512 blocks.** Read from a tile
  of `rasputin_data/DTM10_UTM33_20260925` with tifffile: shape 5051 × 5051,
  float32, compression 5, predictor 1, chunks (512, 512), 100 blocks.
  ANADEM (the São Francisco DEM, `rasputin_data/cache/anadem-v1`) is one
  Cloud Optimized GeoTIFF, Deflate (compression 8), no predictor, 512 × 512
  blocks, 8,700 of them in the cache.
- **`TiffPage.decode` decodes one block with no Python loop over pixels.**
  In `tifffile/tifffile.py` (2026.9.20), the closure `decode_other` (line
  7446) calls the codec with `out=size * dtypeitemsize`, so the size is known
  and the codec makes one pass, then reshapes.
- **imagecodecs' LZW decoder releases the GIL** (the lock that lets only one
  thread run Python at a time). `lzw_decode` in `imagecodecs/_imcd.pyx` wraps
  `imcd_lzw_decode` in `with nogil:` (read on GitHub's current source, not
  the 2026.8.16 tag, whose wheel ships only the compiled module). The
  measurement says the same thing on the installed version: on Numedalslågen
  the 16 windows take 3.07 s on 1 thread and 0.48 s on 8, 6.4 times faster
  (section 2), while the process CPU time of the whole `open_dem` call stays
  about the same: 3.84 s on 1 thread and 4.04 s on 8. A decoder holding the
  GIL could not scale like that. A rerun for design review round 1 (AC
  power, `time.process_time`, median of 3; the DEM code at master
  `ed125121` is the same as at `6c729e97`) gave 3.98 s and 4.51 s for the
  whole call, of which the windows alone are 3.21 s and 3.75 s.
- **Deflate scales the same way.** On 96 cached ANADEM blocks (one window of
  8 × 12 blocks), `decode_window` takes 0.194 s on 1 thread, 0.055 s on 4,
  0.031 s on 8 and 0.029 s on 10.
- **tifffile's own default is not used.** `TiffPage.asarray` would use
  `TIFF.MAXWORKERS`, half the cores unless `TIFFFILE_NUM_THREADS` is set
  (`tifffile.py` line 18440), but `decode_window` runs its own
  `ThreadPoolExecutor` over `page.decode` and never calls `asarray`, so that
  variable does not apply.

**Novelty: none claimed.** Nothing was searched.

### Legacy

```sh
$ git grep -l -i -E "thread|tifffile|imagecodecs|seam|decode" legacy-archive -- legacy
legacy-archive:legacy/rasputin/tin_repository.py
legacy-archive:legacy/rasputin/wfs_repository.py
```

Both hits are text decoding (`n.decode("utf-8")` at
`legacy/rasputin/tin_repository.py:56`, `.decode("utf-8")` at
`legacy/rasputin/wfs_repository.py:39`), not rasters. The legacy reader opens
one GeoTIFF whole with PIL (`legacy/rasputin/reader.py:409`, `Image.open`,
then `np.array(image)` at line 420), on one thread, with no mosaic and so no
seams. **Nothing is carried.**

## 2. What is slow, and why

`@perf`'s profile (section 4 of the bottlenecks README) was on battery, on
master `8199f30`. Measured again here, on **AC power** (`pmset -g batt`
before and after: AC), with a non-editable install of `6c729e97` in a scratch
venv, by a prototype script that times each step of `open_dem` in-process
(median of 5, the variants interleaved; canvas and seams compared with the
base in every run):

| step | Numedalslågen s | Skiensvassdraget s |
|---|---|---|
| read every tile's header (`footprints`, 254 tiles) | 0.084 | 0.084 |
| plan the mosaic on the domain (`_domain_plan`) | 0.103 | 0.075 |
| **decode the windows** (16 / 12 tiles, one after another, 4 threads each) | **0.873** | **0.649** |
| **`_covered`**: every overlap node tested against the needed region, only for the seam report | **0.409** | **0.243** |
| `_decide`, `_seam` | 0.054 | 0.035 |
| canvas allocation, filling and strips (the rest) | about 0.07 | about 0.04 |
| **`open_dem`, whole** | **1.594** | **1.130** |

The CLI's `decode` row is `open_dem` and nothing else
(`src_python/tin_engine/cli.py@6c729e97:1389-1393`), but it is higher: 1.79 s
and 1.31 s through `rasputin mesh` on AC (section 8). That is the first call
in a fresh process. Called three times in one process, `open_dem` took
1.75 s, then 1.61 s and 1.59 s (three processes, the same each time). A
cProfile of the first and second calls puts most of the difference in
`np.full` (the 1.1 GB NaN canvas: 0.091 s, then 0.019 s) and the copies into
the canvas: the first call touches memory fresh from the operating system.
`@perf`'s README left that difference unexplained. It is not attacked here.

**Decode threads.** The windows (16 on Numedalslågen, 12 on
Skiensvassdraget) on 1 to 16 threads (median of 3,
`open_dem` whole in the last column):

| threads | Numedalslågen: windows s | `open_dem` s | Skiensvassdraget: windows s | `open_dem` s |
|---|---|---|---|---|
| 1 | 3.073 | 3.823 | 2.285 | 2.767 |
| 2 | 1.610 | 2.336 | 1.205 | 1.687 |
| 4 (today) | 0.873 | 1.596 | 0.649 | 1.129 |
| 6 | 0.607 | 1.345 | 0.453 | 0.943 |
| 8 | 0.483 | 1.207 | 0.358 | 0.840 |
| 10 | 0.462 | 1.187 | 0.343 | 0.829 |
| 12, 16 | 0.462 | about 1.19 | 0.342 | about 0.83 |

The machine has 10 cores, 8 of them performance cores. Past 8 threads the
gain stops, and more threads than cores cost nothing measurable.

**`_covered`.** With a domain, the seam report counts only the overlap nodes
inside the needed region (Ola, 2026-09-28). Today every node of every overlap
is tested against that region with `shapely.intersects_xy`, 5.2 M and 3.1 M
nodes (`@perf`'s counts), and only then are the nodes the report can count
picked: both values valid and the gap at least `SEAM_THRESHOLD` (1 mm).
Those are 101,520 nodes on Numedalslågen and none on Skiensvassdraget.
Testing the region only on them costs about nothing.

**The most any change can save.** The `decode` row is 1.79 s and 1.31 s on
AC. The two changes above save about 0.74 s and 0.49 s of it (section 8).
What is left is about 1.06 s and 0.82 s:

- the LZW decoding itself, 3.07 s and 2.29 s of CPU time, at best about
  0.38 s and 0.29 s spread over 8 performance cores. This design gets it to
  0.46 s and 0.34 s, so at most about 0.1 s is left there;
- headers, plan, `_decide`, the canvas and the first touch of its memory:
  about 0.5 s and 0.4 s, which no change in this design touches
  (section 11).

So past this design, a further change could save a few tenths of a second
per catchment, no more.

## 3. The blueprint: what changes where

Two files: `src_python/tin_engine/io/cog.py` and
`src_python/tin_engine/mosaic.py`. No change to `open_dem`, `DemRequest`,
`assemble`'s signature, `Mosaic`, `Seam`, the repositories, `cli.py`, any
C++ or binding, and no new dependency.

```
cli._open_dem ──> open_dem(request)                                     unchanged
                    └─ assemble(plan, load, needed, load_window=...)     signature unchanged
                         ├─ per tile, in plan order: load_window(name, window)
                         │     └─ decode_window(source, meta, dtype, window, threads=None)
                         │           pool of (os.cpu_count() or 1) workers when threads is None   CHANGED (was a fixed 4)
                         ├─ per overlapping pair: _decide(...)                                    unchanged
                         └─ per overlapping pair: _seam(first, second, a, b, plan, box, needed)  CHANGED
                               valid in both, gap >= SEAM_THRESHOLD, then inside needed (only those nodes)
                    _covered                                                                     REMOVED
```

### 3.1 `decode_window`'s thread count (today `src_python/tin_engine/io/cog.py@6c729e97:110,146`)

- The keyword becomes `threads: int | None = None`. `None` means
  `os.cpu_count() or 1`, the rule `_open_reprojected` already uses for
  resampling (`src_python/tin_engine/dem_input.py@6c729e97:218`). An `int`
  is used as given, as today.
- The docstring says so, and says the count changes no value (W2).
- No caller passes `threads`, so both repositories (`TiffDemRepository` and
  `CacheRepository`) get the machine's count. The parameter stays, for a
  caller that wants to bound it.

This replaces half of increment 23a-1's decided 7 ("`threads=4`, fixed, not
read from the machine and not a flag"). The other half stands: it is not a
flag. The reason for "fixed" was that the output does not depend on the
count (W2), and that is still true; reading the machine costs no
determinism. `docs/increments/23-basin-scale.md` decided 7 gets a one-line
note pointing here, in this PR.

### 3.2 The seam report (today `src_python/tin_engine/mosaic.py@6c729e97:274-276,522-546`)

`_covered` is removed. `_seam` takes the pair's two whole strips, the overlap
`box` and `needed`, and selects the nodes it counts in this order:

1. `rows, cols = np.nonzero(valid_mask(a, nodata) & valid_mask(b, nodata))`;
2. `gaps = np.abs(a[rows, cols].astype(np.float64) - b[rows, cols].astype(np.float64))`,
   and keep `rows`, `cols`, `gaps` where `gaps >= SEAM_THRESHOLD`;
3. with `needed`: keep the gaps where
   `shapely.intersects_xy(needed, *plan.meta.node_xy(box.row0 + rows, box.col0 + cols))`;
4. none left: `None`; else `Seam(first, second, gaps.size, gaps.max(), median(gaps))`
   as today.

`assemble` calls it as
`_seam(a.name, b.name, strips[a.name, b.name], strips[b.name, a.name], plan, box, needed)`.
`needed` is still prepared once before the loop
(`src_python/tin_engine/mosaic.py@6c729e97:267`). `_seam` stays pure.

### 3.3 As built

Green `0663efc9` follows 3.1 and 3.2: +3 net production lines (19 added,
16 removed: `io/cog.py` +1, `mosaic.py` +2), the prototype's count in
section 8. Where the code says more than the design did, from `@developer`'s
notes:

- `threads=None` reads `os.cpu_count() or 1` at the call, and a given count
  is used as given
  (`src_python/tin_engine/io/cog.py@0663efc9:150`). So `threads=0` still
  raises `ValueError`, from `ThreadPoolExecutor` itself, as it did at the
  base (`src_python/tin_engine/io/cog.py@6c729e97:146`). No caller passes 0.
- `_covered` is gone, and `_seam`'s docstring now carries its note on where
  a node lies ("a node is `x_min + col * dx`, as in `_uncovered`",
  `src_python/tin_engine/mosaic.py@0663efc9:534-535`). `assemble` names the
  two strips `first` and `second` before the call
  (`src_python/tin_engine/mosaic.py@0663efc9:274-275`).

## 4. Why the output cannot change

### 4.1 The thread count

Each block is decoded by tifffile from its own bytes and written into its own
slice of one preallocated array; blocks never share a slice. So the array
does not depend on how many threads run or in which order they finish. That
is W2 (`tests/python/test_io_cog.py@6c729e97:237`, 1 thread against 8,
every block layout and codec of the fixtures), unchanged. The probe (section
6) shows the same on both catchments and every fixture.

### 4.2 The seam report

Today the counted set is `C ∩ V ∩ G`: `C` the overlap nodes `needed` covers,
`V` those valid in both tiles, `G` those whose gap is at least the threshold.
The new code computes `(V ∩ G) ∩ C`, the same set:

- each test is exact and decided per node, with the same inputs: `V` and the
  gap from the same strip values, the same float64 subtraction; `C` from the
  same `intersects_xy` on the same prepared geometry and the same node
  coordinates (`node_xy` is elementwise, `x_min + col * dx` on the same
  integer indices, so the floats are the same);
- `np.nonzero` lists nodes in row-major order and each filter keeps that
  order, so even the order of `gaps` is today's. The count, the largest gap
  and the median would not depend on order anyway.

Without `needed`, today's `kept` is `...` (every node), and the new code
skips step 3: the same set.

## 5. The blueprint against the assessment framework

- *Data and execution:* the thread count is execution, not data, so it is
  not put in `DemRequest`. It is read from the machine at the point of use,
  as `resample`'s already is.
- *State:* `_seam` stays pure; no new state. `decode_window` reads
  `os.cpu_count()`, which does not change during a run.
- *Dependency gravity:* none added.
- *Async-readiness:* unchanged in kind. `open_dem` is still blocking and an
  async caller still runs it in `asyncio.to_thread`. During decoding it now
  uses every core instead of 4. Section 10 says what that means for a server
  running several meshes at once.

## 6. The byte-identical gate

`docs/increments/30c-probes/dem_bytes.py` and its base output
`docs/increments/30c-probes/base_6c729e97.txt` are committed with this
design. The probe's docstring says how to run it. Both modes are the gate:

- **`fixtures`** runs the 14 suites that reach the DEM code in-process
  (`test_mosaic.py`, `test_mosaic_windowed.py`, `test_dem_input.py`,
  `test_dem_input_domain.py`, `test_io_cog.py`, `test_io_repository.py`,
  `test_io_cache.py`, and the `rasputin mesh` suites for DEMs, mosaics,
  domains, domain CRSs, the cache, geographic DEMs and plain output). It
  records every call of `assemble` and of `decode_window` wherever they are
  looked up (`mosaic`, `dem_input`, `catchment`, `io.cog`, `io.repository`),
  keyed by test id and call number. An `assemble` line holds the canvas's
  shape, dtype and a hash of its meta and bytes, and every seam's names,
  count, largest and median (floats as `repr`). A `decode_window` line holds
  the window and a hash of the tile. 1,178 lines at the base (416 `assemble`,
  762 `decode_window`), 86 of them with a non-empty seam report; pytest 709
  passed, 4 skipped.
- **`mesh`** runs `rasputin mesh` on Numedalslågen and Skiensvassdraget as
  `@perf`'s profile did, and records the same line for the assembled DEM and
  a hash of the `.vtk` written. The `.vtk` hashes are 30b's acceptance's
  (`34f7117e5e528e97`, `ab996019190166f9`).

The gate: **every `fixture ` and `mesh ` line of the base file appears,
unchanged, in the branch's run**; the run may also hold `fixture ` lines the
base lacks, but only those of the red suite's tests (section 7), and the
review lists them. With `run.txt` the branch's output and `base.txt` the base
file:

```bash
comm -23 <(grep -E '^(fixture|mesh) ' base.txt | sort) <(grep -E '^(fixture|mesh) ' run.txt | sort)   # must print nothing
comm -13 <(grep -E '^(fixture|mesh) ' base.txt | sort) <(grep -E '^(fixture|mesh) ' run.txt | sort)   # the allowed new lines, listed
```

**The base.** Recorded with a non-editable install of `6c729e97` in a
scratch venv outside the worktree (`git archive`, then `uv pip install
".[codecs]"`; pytest, pytest-asyncio and hypothesis added to that venv), on
AC power: `fixtures` twice and `mesh` twice, identical lines. **The probe
refuses a `tin_engine` outside a `site-packages` directory**, such as an
editable install or a loose copy of the source (30b's lesson). It cannot
tell which commit an installed package came from: a copy under any
`site-packages` passes. The probe's first line prints the package's path,
and the operator checks that it is the scratch venv installed from the
commit being measured. (The installed version is `0.2.0.dev0` at every
commit, and `direct_url.json` names only the directory installed from, so
neither is printed.) Run it from the repository root with the venv's own
`python`, so `tests/python` is the branch's.

**At the merge.** When 30b has merged and master is merged into this branch,
record the base again at that merge commit, the same way, as
`base_<merge>.txt`. 30b's code is already in `6c729e97`, so every `fixture`
and `mesh` line should equal this base's; if a line differs, the merge
changed something the probe reaches, and that is explained before
`@tester` starts (or before green, whichever comes later). The newer file is
then the gate's base.

**It can fail.** Plants in the prototype (section 8), each against the base
lines:

| plant | base `fixture` lines not matched (of 1,178) | pytest | Numedalslågen's `mesh` line |
|---|---|---|---|
| the needed-region test skipped (every valid gap counted) | 8 | 4 fail | differs |
| `gaps > SEAM_THRESHOLD` instead of `>=` | 2 | 2 fail | same |
| node coordinates without the overlap's origin (`node_xy(rows, cols)`) | 4 | 2 fail | differs |

**The prototype passes the gate**: all 1,178 `fixture` lines and both `mesh`
lines equal, `.vtk` hashes included.

**What the gate cannot see**: how many threads ran and how many nodes were
tested against the region. Both are speed, not output, so only the red tests
(R1, R2) and `@perf`'s timing guard them.

## 7. Tests for `@tester` (the red suite)

Behaviour does not change, so the red tests are the contract of the two
changes. The pins are green today, which is their point; they are committed
with the red ones, before any code. No mutation round. Both red tests were
tried in a scratch file against the base install (both fail: the pool is
made with 4 workers; 102 nodes tested where 2 count) and the prototype (both
pass).

**Red.**

- **R1. `decode_window` without `threads` uses the machine's cores**
  (`tests/python/test_io_cog.py`). With `os.cpu_count` patched to return 7
  and `cog.ThreadPoolExecutor` replaced by a spy that records its worker
  count and builds the real pool, a `decode_window` call that omits
  `threads` builds its pool with 7 workers; with `os.cpu_count` returning
  `None`, with 1. Any fixture page will do.
- **R2. The seam report tests the needed region only on nodes it could
  count** (`tests/python/test_mosaic.py`). Two overlapping tiles with a
  planted disagreement at a few overlap nodes, some inside `needed` and at
  least one outside, and many agreeing overlap nodes; `assemble(plan, load,
  needed)` with `mosaic.shapely.intersects_xy` wrapped by a spy that sums the
  points it is given (patch after `plan_mosaic`, which also calls it). The
  sum equals the number of planted nodes (valid in both, gap at least
  `SEAM_THRESHOLD`), and the seam counts only those inside.

**Pins (green today, must stay green).**

- **P1. The seam report's edge cases under a needed region**
  (`tests/python/test_mosaic.py`), against an oracle computed in the test the
  way today's code does it (every overlap node tested against `needed`, then
  both valid, then the gap). In one overlap: a qualifying node exactly on
  `needed`'s boundary (counted: closed); a gap exactly `SEAM_THRESHOLD` in
  float64 (counted; float64 tiles, for example 0.001 against 0.0) and one
  just below (not counted); NaN on one side; the NoData sentinel on one side;
  `+inf` against a finite value (counted, largest `inf`); `+inf` on both
  sides (`inf - inf` is NaN: not counted). `Seam`'s count, largest and median
  equal the oracle's, for an odd and an even count (the median of an even
  count is the mean of the middle two).
- **P2. The default thread count changes no byte**
  (`tests/python/test_io_cog.py`, next to W2): over W2's variants, a call
  without `threads` gives the same bytes as `threads=1`.

The existing tests (`TestW2Determinism`, the seam tests of `test_mosaic.py`
and `test_dem_input_domain.py`, `test_cli_mesh_mosaic.py`) stay as they are.

### What `@tester` pinned beyond this section, and the ruling

Ruled after the red commit `bdf7b57d` and before green, on the tests it added
(`tests/python/test_io_cog.py@bdf7b57d:253-304`,
`tests/python/test_mosaic.py@bdf7b57d:1544-1659`). Checked in a scratch copy
of the source with its own venv: at `bdf7b57d` R1's `seven_cores` and
`cores_unknown` and R2 fail, everything else in both files passes; with
sections 3.1 and 3.2 put in, both files pass. Each pin was then broken in
that copy (the plants named below).

`@architect`'s ruling, 2026-10-06: **six stand; pin 3 changes.**

1. *R1 patches `os.cpu_count`, so the code calls `os.cpu_count()` when it
   decodes.* Keep: 3.1 writes `os.cpu_count() or 1` and section 5 says it is
   read at the point of use. `from os import cpu_count`, and a core count
   read once at import, each fail `seven_cores` and `cores_unknown`.
2. *R1 watches `cog.ThreadPoolExecutor`: one pool per call, the count
   positional or `max_workers`.* Keep: the pool per call is today's and 3.1
   keeps it; `max_workers=` with a `thread_name_prefix` passes, a pool built
   through `concurrent.futures` directly fails all three cases.
3. *R1's third case: `threads=3` gives 3 workers.* Change: 3 is below the 7
   patched cores, so a code that caps `threads` at the core count passes it,
   and 3.1 says an `int` is used as given. `@tester` changes that case to
   `pytest.param(7, {"threads": 9}, 9, id="threads_given")`
   (`tests/python/test_io_cog.py@bdf7b57d:278`). Tried: it passes with 3.1,
   fails a cap at the cores and fails `threads` ignored, and passes at
   `bdf7b57d` as the 3 did.
4. *R2: the points given to `intersects_xy` during `assemble` sum to 4, in
   any number of calls.* Keep: that is R2. One call per node passes; testing
   the region before the gap filter (on the 54 valid nodes) fails. The spy
   counts the `x, y` form 3.2 step 3 writes; a single `(n, 2)` coordinate
   array counts as 8 and fails, so `@developer` uses 3.2's form.
5. *R2 and P1: a needed region with a hole, and the canvas is the whole
   grid.* Keep: the hole puts overlap nodes inside the region's box but
   outside the region, which a region test that is skipped would count (it
   fails R2 and both P1 cases). With the overlap at (0, 0), node
   coordinates without the overlap's origin pass here; the existing
   `TestSeamsInsideTheNeededRegion` in `test_dem_input_domain.py` fails on
   that plant, and so does the probe (section 6).
6. *P1 ignores the `inf - inf` RuntimeWarning.* Keep: the design does not
   say whether the subtraction warns. 3.2 as written and a version under
   `np.errstate(invalid="ignore")` both pass, also with
   `-W error::RuntimeWarning`.
7. *P1: float64 tiles, NoData -32767, the boundary node on the hole's
   edge.* Keep: 7's P1 asks for float64 (0.001 exactly) and a node on the
   region's boundary. `contains_xy` for `intersects_xy`, `>` for `>=`, the
   sentinel not treated as NoData, and the gap taken in float32 each fail
   both P1 cases.

## 8. Net production lines, and the prototype

The prototype is the installed `6c729e97` package in a copy of the base's
scratch venv, edited in place (no C++ build; removed after measuring). It
is sections 3.1 and 3.2 as written, with docstrings.
`python3 tools/count_loc.py` between its base and the edit: **+3 net** (19
added, 16 removed: `io/cog.py` +1, `mosaic.py` +2). Estimate for the PR:
**+3 to +10**.

Timings through `rasputin mesh --stats` on both catchments, base and
prototype alternated, three repeats each, **AC power** for all 24
`pmset -g batt` readings, each run printing which `tin_engine` it imported
(the base's 6 runs the base venv, the prototype's 6 the prototype's).
Medians:

| | Numedalslågen: `decode` s | total s | Skiensvassdraget: `decode` s | total s |
|---|---|---|---|---|
| base `6c729e97` | 1.794 | 7.711 | 1.307 | 11.765 |
| prototype | 1.057 | 6.939 | 0.815 | 11.280 |
| prototype / base | 0.589 | 0.900 | 0.624 | 0.959 |

Every run of a catchment wrote the same `.vtk` (`34f7117e5e528e97`,
`ab996019190166f9`). Per repeat, `decode`: base 1.794 / 1.788 / 1.822 and
1.323 / 1.307 / 1.289; prototype 1.057 / 1.111 / 1.046 and 0.815 / 0.815 /
0.811.

In-process (section 2's script, median of 5): `open_dem` 1.594 → 0.820 s
and 1.130 → 0.589 s with both changes on 10 threads; the thread change alone
1.181 and 0.821 s, the seam change alone 1.236 and 0.895 s.

A first timing run was void and was thrown away: the copied venv's
`bin/rasputin` script names the original venv's interpreter on its first
line, so the "prototype" runs ran the base code (they timed the same as the
base). The rerun called each venv's own `python` and printed
`tin_engine.__file__` per run. Section 9 asks `@perf` to do the same.

The suites the probe runs took 21.3 to 21.5 s on both the base and the
prototype (two runs each), so the extra threads cost the tests nothing.

## 9. Acceptance (`@perf`)

`tools/bench.py` and the thread sweep are not required: nothing under
`include/terrain/refinement/` or `include/terrain/mesh/` changes. `@perf`
runs, with the branch's own non-editable install and a base install of the
merge base, **the same power state for both, recorded with `pmset -g batt`
before and after each run**, and **each run printing the `tin_engine` path
it imported** (call each venv's own `python`; never a copied venv's entry
script):

1. **Byte-identical**: the probe's `mesh` mode on both catchments, every line
   equal to the gate's base (section 6), and its `fixtures` mode, every base
   line present and unchanged, the red suite's new lines the only additions.
2. **Time**: `rasputin mesh --stats` on both catchments, three repeats each,
   base and branch alternated, as
   `docs/benchmarks/2026-10-06/30b-clip/scripts/` does. Read the `decode`
   row's median. **Pass:** the branch's median is at most 0.7 of the base's
   on both catchments. Checked only at these two inputs (16 and 12 tiles,
   280 M and 200 M canvas nodes, on a 10-core machine); the prototype gave
   0.59 and 0.62.
3. The whole-run total and the other rows, recorded, not gated.

Evidence goes under `docs/benchmarks/<date>/30c-dem-read/`, with a README
whose numbers come from a script over the raw files (30b's
`scripts/summarize.py` is the model), not typed by hand.

**Result: ACCEPTED** (`@perf`, 2026-10-06, `bc8d91eb`; evidence
`docs/benchmarks/2026-10-06/30c-dem-read/README.md`, its tables written by
`scripts/summarize.py` from `raw/`). Base `6c729e97` against branch
`0663efc9`, both non-editable installs with the same `_core`, AC power for
all 24 `pmset -g batt` readings and no "Using Batt" line in `pmset -g log`
during any kept run, median of 3:

| catchment | decode, base → branch, s | branch / base | total, base → branch, s |
|---|---|---|---|
| Numedalslågen | 1.788 → 1.051 | **0.588** (pass) | 7.480 → 6.719 |
| Skiensvassdraget | 1.323 → 0.842 | **0.636** (pass) | 11.462 → 10.966 |

Each catchment's six runs wrote one `.vtk` sha256 (`34f7117e5e528e97`,
`ab996019190166f9`, `raw/stats/vtk_sha256.txt`), the hashes of section 6,
and one `dem_seams` value. Two runs were redone because `@perf`'s own shell
was busy beside them, none for power (`raw/discarded/README.txt`).

**The probe** (`raw/probe_compare.txt`): on both installs every line of
`base_6c729e97.txt` (1,178 `fixture`, 2 `mesh`) is present and unchanged,
and the 16 new `fixture` lines are all from the red suite's tests (R1's
three cases, P2's ten, P1's two, R2's one). On the branch pytest exited 0; on
the base the 3 red tests fail, as they should.

This acceptance compares `6c729e97` with `0663efc9`, both before the master
merge. The merge changes no DEM code (`git diff 0663efc9 26a5d839 --
src_python/tin_engine/io/cog.py src_python/tin_engine/mosaic.py
src_python/tin_engine/dem_input.py` is empty), but it does change the domain reading the probe's
`mesh` mode goes through (`domain.py`, `io/domain_file.py`, `io/geojson.py`).
Section 6's re-recorded base at the merge settles whether any probe line
moved.

## 10. Risks

- **Several meshes on one machine.** Each `decode_window` now starts up to
  `os.cpu_count()` threads. A server running several `open_dem` calls at once
  would run more decode threads than cores. The output is unchanged and the
  operating system shares the cores; only the decode's speed is shared. The
  `threads` parameter is the place to bound it if that server exists. Not
  measured.
- **Memory.** Each running thread holds one decoded block (512 × 512 float32,
  1 MB) and its compressed bytes: about 20 MB at 10 threads instead of about
  8 MB at 4, against a 1.1 GB canvas. Not measured beyond this arithmetic.
- **Efficiency cores.** `os.cpu_count()` counts the 2 efficiency cores of
  this machine too. Section 2's table shows 10 threads no slower than 8.
- **A tifffile or imagecodecs change** that held the GIL while decoding
  would make the extra threads useless but harmless. Nothing in CI would
  notice; `@perf`'s timings would.

## 11. Not in scope

- **Decoding several tiles at once.** The windows are decoded one tile after
  another, and the next tile waits for the slowest block of the last. Taking
  all blocks of all windows into one pool would get the windows from 0.46 s
  toward the 0.38 s floor of section 2 on Numedalslågen, but every decoded
  window would then be alive at once, which doubles the peak memory
  `assemble`'s docstring promises (the canvas plus about 3 tiles). Not worth
  0.1 s.
- **Reading 254 headers** (0.084 s) when 16 tiles are used: a header cache
  or a parallel read. Small.
- **The first touch of the canvas's memory** (about 0.15 s on the first call,
  section 2). Inherent in allocating a 1.1 GB array.

## 12. What 30c is worth

After 30a and 30b, `decode` is the largest phase but one on Numedalslågen
(1.79 s of 7.7 s, 23 %; `features clip` is 1.87 s; medians of section 8's base runs) and 11 % of the run on
Skiensvassdraget (1.31 s of 11.8 s). 30c takes 0.74 s and 0.49 s off, about
10 % and 4 % of the whole run. After it, `features clip`, `refine` and
`decode` lead on Numedalslågen (about 1.85, 1.19 and 1.06 s), and `refine`
and `land cover` on Skiensvassdraget.

**Against São Francisco** (the real target: about 635,000 km², ANADEM,
30 m): the seam change saves nothing there, because ANADEM is one file, so
there is no overlap and no seam report. The thread change applies: ANADEM's
Deflate blocks decode at about 0.58 ms each on 4 threads and 0.31 ms on 10
(section 1). For the 8,700 blocks in the cache, that is about 5.0 s and
2.7 s, so 30c saves about 2.3 s. Increment 23's arithmetic puts the rest of
that run in minutes (`docs/increments/23-basin-scale.md`, "The basin at 1 m,
as arithmetic": about 40 s of refining on 8 cores and about 24 s of
resampling). So on São Francisco 30c is worth a few per cent at most. The
96 blocks timed are one window; the rest of the cache was not timed.

In short: 30c is cheap (a few lines, a lean round) and visibly helps the two
Norwegian catchments, but it does little for São Francisco.

## 13. Questions for Ola

1. **Build 30c, or stop at 30b?** It saves about 0.7 s of a 7.7 s run on
   Numedalslågen and 0.5 s of 11.8 s on Skiensvassdraget, and about 2 s of a
   run of minutes on São Francisco, for about +3 to +10 lines and a lean
   round (tests, code, review, `@perf`'s acceptance). **Default: build it**,
   since you approved the order and the round is small.
   **Ola, 2026-10-06: "yes, build 30c".**

## 14. ROADMAP

Row 30's 30c entry now says what this design does: "30c reading the DEM
(decode threads from the machine's cores, the seam report's region test on
the nodes it counts only)", and its status says 30c is designed and waits
for question 1. Since Ola's yes, it says 30c is designed and in design
review.

## Review

**30c, design review, round 1, 2026-10-06.** Range `6c729e97..a3a7cfef`. Verdict: CHANGES REQUESTED. LOC 0. Blocking: (1) the 2-line note at `/Users/skavhaug/projects/rasputin/.claude/worktrees/dem-read-speed/docs/increments/23-basin-scale.md@a3a7cfef:1066-1067` shifts the self-citations at `@a3a7cfef:3281,3283` (`:3267,3271,3279` now point at 23b rounds 5 and 7, and at 23c-1 round 3); (2) `/Users/skavhaug/projects/rasputin/.claude/worktrees/dem-read-speed/docs/increments/30c-dem-read-speed.md@a3a7cfef:298-300` over-claims the refusal at `/Users/skavhaug/projects/rasputin/.claude/worktrees/dem-read-speed/docs/increments/30c-probes/dem_bytes.py@a3a7cfef:180` (a copy under any `site-packages` passes). Base reproduced exactly in a scratch install, the 4.2 same-set argument holds, and W2 plus my own LZW timing confirm the thread claim and GIL release.

**30c, design review, round 2, 2026-10-06.** Range `a3a7cfef..c2ec195b` (one commit, c2ec195b). Verdict: APPROVED. LOC 0. Both round-1 blockers are closed. (1) The note at `/Users/skavhaug/projects/rasputin/.claude/worktrees/dem-read-speed/docs/increments/23-basin-scale.md@c2ec195b:1064-1065` replaces the old two lines. The file is back to 3281 lines, the same as master `ed125121`. The self-citations, now at `:3279,3281`, again point at 23b round 6 (`:3267`), 23b round 8 (`:3271`) and the row-23 fix round 1 (`:3279`). (2) `/Users/skavhaug/projects/rasputin/.claude/worktrees/dem-read-speed/docs/increments/30c-dem-read-speed.md@c2ec195b:304-312` now describes what `/Users/skavhaug/projects/rasputin/.claude/worktrees/dem-read-speed/docs/increments/30c-probes/dem_bytes.py@c2ec195b:179-181` actually checks. Both suggestions were taken. The DEM code at `ed125121` is the same as at `6c729e97`. Nothing else changed. Not pushed; no CI.
