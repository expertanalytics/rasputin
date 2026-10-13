# Increment 38: threads and memory as settings

**Status: designed and review-approved after two rounds; red `4b9c0dfe` (`@tester`, 110 tests); green `b252a70d` (`@developer`, 96 net production lines by `python3 tools/count_loc.py 5d2e706d HEAD` against the estimate of 60, under the worst case of 120); test-side fixes `9de3212c` (`@tester`); section 12 is the as-built record; `@perf`'s quick check NO CHANGE, ACCEPTED; code review round 1 two prose blockers, fixed; round 2 APPROVED; next the push on Ola's yes;** by `@architect`, 2026-10-12, on master `5d2e706d` (branch
`worktree-threads`). The first fix from the generality audit
(`docs/increments/generality-audit.md`, section 6 item 1; rows D1, D2, D3, D5,
D6, D7). Ola's yes, 2026-10-12: "1, 2: yes" to a single pull request that adds
`--threads` (0 = every core, an environment variable as the default) to `mesh`
and `catchment`, passes it to every C++ call and Python pool, drops the fixed
cap of 8 threads on reprojection, prints the count used in the stats file, adds
`--memory-budget` honouring a container's limit, and deletes the unused
`decompose.py`. One PR, Python only, no C++ change. Full path (the diff is over
50 lines), no `@perf` acceptance run by the rule (no refine or mesh code is
touched), the quick check as the speed statement (section 9).

## 1. What Ola asked

Ola, 2026-10-12: "the number of cores should of course not be hardcoded. I'm a
bit shocked about the 'abstraction neglect' I've seen." Today every run uses
every core, one pool stops at 8 whatever the machine, the stats file prints
the core count as if it were the thread count, and the two memory caps read
the host's memory inside a container. Nothing lets a user, a shared server, a
laptop on battery, a CI job or an agent running several meshes at once say how
much of the machine a run may take.

## 2. Prior art: legacy and literature

**Legacy.** The legacy tree was serial and read no machine limit; nothing is
carried. The greps, run on this tree against the `legacy-archive` tag:

```
$ git grep -l -i -E 'cpu_count|hardware_concurrency|omp_get|num_threads|ThreadPool|multiprocessing' legacy-archive -- legacy
(nothing; exit 1)
$ git grep -l -i -E 'memory|rlimit|sysconf' legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp      # one hit: `#include <memory.h>`, not a limit
```

**Literature.** No novelty is claimed: "how many threads" is an engineering
convention, and the design follows the ones in use. What they do, each read
from its own documentation on 2026-10-12:

- *OpenMP* (specification 5.0, section 6.2, `OMP_NUM_THREADS`): "The value of
  this environment variable must be a list of positive integer values"; a
  value that "is not a positive integer", or more threads than the
  implementation supports, is "implementation defined". So 0 has no meaning in
  OpenMP, and a count above the machine is left to the implementation (GCC and
  LLVM honour it).
- *GNU coreutils `nproc`*: prints "the number of processing units available to
  the current process, which may be less than the number of online
  processors", so a `taskset` affinity mask is honoured; `OMP_NUM_THREADS` and
  `OMP_THREAD_LIMIT` "determine the minimum and maximum returned value".
- *Rayon* (`ThreadPoolBuilder::num_threads`): "If `num_threads` is 0, or you do
  not call this function, then the Rayon runtime will select the number of
  threads automatically"; the automatic count "is based on the
  `RAYON_NUM_THREADS` environment variable (if set), or the number of logical
  CPUs (otherwise)". Rayon's 0 means "choose for me", which is what the C++
  core already does (`include/terrain/parallel_util/chunks.hpp@5d2e706d:17`,
  "threads == 0 means hardware_concurrency"), and what `--threads 0` means here.
- *threadpoolctl* (`threadpool_limits`): "If None, this function does not do
  anything"; an integer "sets the maximum number of threads"; no zero
  convention; it reads no environment variable itself. Nothing to borrow
  beyond "None means leave it".
- *Python 3.13* (`os.process_cpu_count()`): "the number of logical CPU cores
  usable by the calling thread of the current process", so affinity counts;
  both it and `os.cpu_count()` "can be overridden through the new environment
  variable `PYTHON_CPU_COUNT` or the new command-line option `-X cpu_count`",
  "useful for users who need to limit CPU resources of a container system".
  The design uses it where present (this tree runs 3.14; `pyproject.toml`
  allows 3.12, where it falls back to `os.cpu_count()`).
- *The JDK* (JDK-8146115, `-XX:+UseContainerSupport`, on by default since
  OpenJDK 8u191 and 11): the CPU count and the memory limit are read from the
  cgroup file system, `-XX:ActiveProcessorCount` overrides the count,
  `-XX:MaxRAMPercentage` sizes the heap as a share of the container's memory.
  That is the shape here: a container-aware default, a user override, a share
  (half) as the memory default.
- *cgroup v2* (kernel `admin-guide/cgroup-v2`): `memory.max` is "a read-write
  single value file which exists on non-root cgroups. The default is 'max'",
  amounts in bytes; `cpu.max` is "$MAX $PERIOD", default "max 100000". In
  cgroup v1 the memory limit is `memory/memory.limit_in_bytes`, which reads as
  a very large number (about 2^63) when there is none.

What differs from all of them: nothing. The one choice they do not settle, a
count above the core count, is honoured (section 4.2), as OpenMP's
implementations and Rayon honour it.

## 3. Measured on this tree (`5d2e706d`)

The three inventories, by grep, so the reader can re-run them:

**Machine reads and pools.** `grep -rn -E 'cpu_count|ThreadPoolExecutor\(|sysconf' src_python/tin_engine`:

```
src_python/tin_engine/dem_input.py:218:    source, threads = TileWindows(mosaic.tile), os.cpu_count() or 1
src_python/tin_engine/target_grid.py:156:    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
src_python/tin_engine/crs.py:39:WORKERS = min(8, os.cpu_count() or 1)
src_python/tin_engine/crs.py:73:        with ThreadPoolExecutor(len(chunks)) as pool:
src_python/tin_engine/mosaic.py:57:    return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
src_python/tin_engine/cli.py:1590:            threads=os.cpu_count() if sizes.start_vertices is not None else None,
src_python/tin_engine/io/cog.py:124:    `threads` decoding workers, the machine's cores (`os.cpu_count() or 1`)
src_python/tin_engine/io/cog.py:203:    with ThreadPoolExecutor(threads if threads is not None else os.cpu_count() or 1) as pool:
```

Three pools: the reprojector (`src_python/tin_engine/crs.py@5d2e706d:73`, up
to `WORKERS` = min(8, cores)), the target grid's `_pool`
(`src_python/tin_engine/target_grid.py@5d2e706d:156`, serving `resample` at
:200 and `check_point_blocks` at :258, with the count from
`src_python/tin_engine/dem_input.py@5d2e706d:218`), and the block decoder
(`src_python/tin_engine/io/cog.py@5d2e706d:203`, serving `decode_window` :141
and `decode_class_window` :169, called without a count from
`src_python/tin_engine/io/repository.py@5d2e706d:148`, `:281` and `:285`). Four core-count reads
(`src_python/tin_engine/crs.py@5d2e706d:39`,
`src_python/tin_engine/dem_input.py@5d2e706d:218`,
`src_python/tin_engine/io/cog.py@5d2e706d:203`,
`src_python/tin_engine/cli.py@5d2e706d:1590`) and one memory read
(`src_python/tin_engine/mosaic.py@5d2e706d:57`), used as a cap at
`src_python/tin_engine/mosaic.py@5d2e706d:187` and
`src_python/tin_engine/catchment.py@5d2e706d:263`.

**`_core` entry points that take `threads`.** `grep -n 'threads' src_python/tin_engine/_core.pyi`
less the docstring lines: five, at `src_python/tin_engine/_core.pyi@5d2e706d:341`
(`SlopeTolerance.__init__`), `:450` (`refine`), `:562` (`refine_strip`), `:584`
(`refine_points`), `:764` (`class_areas`). Their call sites,
`grep -rn -E '\b(refine|refine_points|refine_strip|SlopeTolerance|class_areas)\(' src_python/tin_engine`
less the stub and the definitions:

```
src_python/tin_engine/cli.py:2137:                slope = SlopeTolerance(to_core(tile), *ramp)
src_python/tin_engine/cli.py:2139:        out = refine(
src_python/tin_engine/edge_strip.py:48:    out = refine_strip(
src_python/tin_engine/final_check.py:60:    out = refine_points(
src_python/tin_engine/land_cover_fractions.py:110:        parts.append(_core.class_areas(local, classes.array, row_areas(classes), threads=threads))
```

Five call sites. Four pass no `threads` (the C++ default 0 = every core); the
fifth passes `assign`'s `threads=0` default
(`src_python/tin_engine/land_cover_fractions.py@5d2e706d:94`), which no caller sets
(`src_python/tin_engine/cli.py@5d2e706d:1447`). The only route to another count today is the
benchmark's monkeypatch of `refine` (`tools/bench.py@5d2e706d:788-790`).

**Reprojector call sites.** `grep -rn 'reprojector(' src_python/tin_engine`
less the definition: 15, the grep as run on `5d2e706d`:

```
src_python/tin_engine/land_cover_fractions.py:67:    lonlat = reprojector(mesh_crs, meta.crs)(xy)
src_python/tin_engine/catchment_batch.py:180:    move = reprojector(stations_crs, segments_crs)
src_python/tin_engine/domain.py:49:        move = reprojector(source, target)
src_python/tin_engine/target_grid.py:136:    back = reprojector(grid.crs, meta.crs)
src_python/tin_engine/target_grid.py:181:        lonlat = reprojector(grid.crs, m.crs)(grid.xy(r0, r1))
src_python/tin_engine/target_grid.py:239:        move = reprojector(m.crs, grid.crs)
src_python/tin_engine/dem_input.py:212:    ((ax, ay),) = reprojector(target, meta.crs)([(c.x, c.y)])
src_python/tin_engine/catchment.py:222:    ((x, y),) = reprojector(request.seed_crs, dem_crs)([request.seed])
src_python/tin_engine/catchment.py:373:    ((px, py),) = reprojector(request.seed_crs, request.lakes_crs)([request.seed])
src_python/tin_engine/catchment.py:381:    move = reprojector(request.lakes_crs, dem_crs)
src_python/tin_engine/tolerance_field.py:83:    to_mesh = reprojector(found.crs, mesh_crs)
src_python/tin_engine/feature_input.py:183:        xy = reprojector(dem_crs, source_crs)(xy)
src_python/tin_engine/feature_input.py:339:        move = None if same_crs(src, self.dem) else reprojector(src, self.dem)
src_python/tin_engine/feature_input.py:563:        ((x, y),) = reprojector(point_crs, own)([point])
src_python/tin_engine/cli.py:2282:    ((x, y),) = reprojector(seed_crs, crs)([seed])
```

The pool opens at 500 000 points (`2 * CHUNK_POINTS`,
`src_python/tin_engine/crs.py@5d2e706d:70-73`). Three sites move arrays that
reach it in a real run: every feature vertex
(`src_python/tin_engine/feature_input.py@5d2e706d:339`, inside `_Tally._take`,
built only when the source's CRS differs from the DEM's), every mesh vertex
for the class fractions (`src_python/tin_engine/land_cover_fractions.py@5d2e706d:67`),
and the target grid's resample blocks
(`src_python/tin_engine/target_grid.py@5d2e706d:181`, up to `BLOCK_NODES` =
1 048 576 nodes a block, `src_python/tin_engine/target_grid.py@5d2e706d:36`
and `:163`, called inside `_pool`'s workers, so a block over the threshold
today opens a second pool of up to 8 threads inside each of the first pool's
workers). A fourth site runs inside `_pool`'s workers too,
`src_python/tin_engine/target_grid.py@5d2e706d:239` (`check_point_blocks`),
but moves a block's rim and then its valid nodes, at most `BLOCK` x `BLOCK` =
512 x 512 = 262 144 points (`:32`, `:238`, `:252`), under the threshold, so it
never nests a pool; it takes `threads=1` all the same (4.4). The other eleven
move a point (`dem_input.py:212`, `catchment.py:222`, `:373`,
`feature_input.py:563`, `cli.py:2282`), a station list
(`catchment_batch.py:180`), a domain outline (`domain.py:49`), the target
grid's rim (`target_grid.py:136`), the lakes (`catchment.py:381`), the
feature hull ring densified to 1 km (`feature_input.py:183`) or one tolerance
line's vertices at a time (`tolerance_field.py:83`, one `shapely.get_parts`
part per call): none reaches 500 000 points in a run, and all become serial by
the library default of 4.2 with no change in behaviour.

**This machine**, for the fixture figures below: 10 cores (`sysctl hw.ncpu`;
8 performance, 2 efficiency), 32 GiB (`hw.memsize` 34359738368), no cgroup;
`os.cpu_count()` and `os.process_cpu_count()` both 10 (Python 3.14.7).

**`decompose.py`**: 93 lines, about 55 counted (blank, comment and docstring
lines dropped by hand with the counter's rule); its three commits are 23c-1's
green steps; nothing imports it (`grep -rn decompose src_python tools`: the
module and one docstring, `src_python/tin_engine/io/mesh_index.py@5d2e706d:37`); its only importer is
`tests/python/test_decompose_partition.py` (263 lines).

## 4. The design

### 4.1 One place reads the machine

A new module `src_python/tin_engine/resources.py`, Pydantic and `os` only,
no `_core`, no typer. It is the one module allowed to read the machine; a test
greps for it (section 8, T1).

```python
class Resources(BaseModel):
    """What a run may use of the machine, resolved: never 0, never None."""
    model_config = ConfigDict(frozen=True)
    threads: int = Field(ge=1)          # the count every pool and _core call receives
    memory_budget: int = Field(ge=1)    # bytes one DEM window may take, flood or mosaic

    @classmethod
    def fallback(cls) -> Resources: ...  # one thread, half the machine's memory

def machine_cores(cgroup: Path = CGROUP) -> int: ...
def machine_memory(cgroup: Path = CGROUP) -> int: ...
def parse_bytes(text: str) -> int: ...
def resources(threads: int | None, memory_budget: str | None, environ: Mapping[str, str]) -> Resources: ...
```

- `machine_cores()`: `os.process_cpu_count()` where Python has it (3.13+;
  affinity and `PYTHON_CPU_COUNT` honoured), else `os.cpu_count()`, else 1;
  then the smaller of that and the cgroup v2 CPU quota, `ceil(MAX / PERIOD)`
  from `<cgroup>/cpu.max`. The memory rule applies to the file: a quota
  counts only when the file is readable and holds two words, the first a
  positive whole number (`max` is no quota) and the second a positive whole
  number (a zero or missing period is no quota, never a division); anything
  else is no quota. Never under 1.
- `machine_memory()`: `os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")`
  (what `mosaic.physical_memory` does today, which moves here), then the
  smaller of that and the cgroup limit: `<cgroup>/memory.max` (v2), else
  `<cgroup>/memory/memory.limit_in_bytes` (v1). A file that is absent,
  unreadable, reads `max`, or does not hold a positive whole number is no
  limit; v1's "no limit" value (about 2^63) loses to the physical figure by
  the `min` and needs no special case.
- `CGROUP = Path("/sys/fs/cgroup")`: the root as a container runtime mounts
  the container's own cgroup. A limit set on an ancestor cgroup of a process
  on a bare host (a systemd slice) is not read; stated as a limit in the
  module docstring. The scale assumed: one container, one cgroup root; checked
  against synthetic files only (T4), there being no container on this machine.
- `parse_bytes("8G")`: a whole number with an optional suffix K, M, G or T
  (either case), 1024-based; anything else, or 0, is a `ValueError` naming the
  text.
- `resources(threads, memory_budget, environ)`: the resolution rule of 4.2 and
  4.3; raises `ValueError` with plain words, which `cli.py` turns into a usage
  error with `param_hint` `--threads`, `RASPUTIN_THREADS` or `--memory-budget`.
  Pure: `environ` is passed in, as `cli.cache_root(option, environ)` already
  does for `RASPUTIN_DATA` (`src_python/tin_engine/cli.py@5d2e706d:1756`).

### 4.2 Threads: the rule

| Given | Count used |
|---|---|
| `--threads N`, N ≥ 1 | N, even above the core count |
| `--threads 0` | `machine_cores()` |
| no option, `RASPUTIN_THREADS=N` | as the option: N ≥ 1 is N, 0 is every core |
| no option, `RASPUTIN_THREADS` unset or empty | `machine_cores()` |
| `--threads -1` | refused by Typer (`min=0`) |
| `RASPUTIN_THREADS=abc` or `-1` | usage error: "RASPUTIN_THREADS=abc is not a whole number of threads (0 means every core)" |

Settled (the first two ruled by Ola, 2026-10-12, section 10):

- **A count above the core count is honoured, not clamped.** Reasons: OpenMP
  and Rayon honour it; the C++ honours it today
  (`include/terrain/parallel_util/chunks.hpp@5d2e706d:53` takes `min(threads, blocks)` workers, no core
  check); the benchmark's thread sweep runs to 20 on this 10-core machine
  (`tools/bench.py@5d2e706d:675`) and must keep meaning "20 threads"; and a
  user asking for more than the cores is making a measurement, not a mistake.
  The stats file shows the count, so the choice is visible.
- **`RASPUTIN_THREADS=0` means every core**, the same as `--threads 0`: one
  rule for both sources, Rayon's convention, and the C++'s own.
- **The C++ keeps taking a plain integer.** Resources is a Python model; the
  bindings and their default `threads = 0` are untouched (no C++ file in the
  diff). Python always passes a resolved count of 1 or more, so the C++ never
  sees 0 from the CLI; its 0 stays for direct callers and the C++ tests.
- **The pools' library defaults become one worker.** `reprojector(src, dst, *,
  threads=1)`, `decode_window(..., threads=1)`, `decode_class_window(...,
  threads=1)`, `_decode_into(..., threads)`: no module but `resources.py`
  reads the core count. A library caller that passes nothing gets a serial
  pool; the CLI passes its count everywhere a run's large arrays go (4.4). The
  reprojector's `WORKERS` constant goes; `CHUNK_POINTS` stays (a chunk of
  work, not a machine figure; the "about 60 ms" comment is reworded as one
  machine's figure).
- **No clamp on the chunked reprojector either**: with `threads=9` and
  2 250 000 points (9 chunks of 250 000), nine workers, where today's cap gives
  eight. With fewer chunks than threads, as many workers as chunks (today's
  rule, kept: `min(threads, len(points) // CHUNK_POINTS)`).

### 4.3 Memory: the rule

| Given | Budget |
|---|---|
| `--memory-budget 8G` (or `512M`, `1T`, a plain number of bytes) | that many bytes |
| no option | `machine_memory() // 2`: half the machine's memory, the container's limit counting as the machine's |
| `--memory-budget 0`, `8GB`, `-1`, `abc` | usage error naming the text |

The budget caps both DEM windows a run assembles: the mosaic canvas
(`src_python/tin_engine/mosaic.py@5d2e706d:187`, today "half the physical memory") and the flood
window (`src_python/tin_engine/catchment.py@5d2e706d:263`, the same words). `plan_mosaic` gains a
keyword `memory_budget: int | None = None` (None: `Resources.fallback()`'s);
the refusal reads "the window is R x C nodes, N bytes to flood, over the
memory budget of B bytes (set with --memory-budget; half the machine's memory
when it is not set)", and the mosaic's "the request is R x C nodes, N bytes
at float32, over the memory budget of B bytes (set with --memory-budget; half
the machine's memory when it is not set); narrow it with --bbox". Both
sentences are true whether or not a budget was given, which the earlier
wording ("half the machine's memory") was not. `mosaic.physical_memory` is
deleted; its one test (M5, `tests/python/test_mosaic.py@5d2e706d:658-666`)
moves to `resources`.

Ruled by Ola, 2026-10-12 (section 10): **`--memory-budget` goes on `mesh` as
well as `catchment` and `station-catchments`**, because `mesh`'s mosaic cap
reads the host's memory in a container too (the same D6 defect, cited in the
audit's D6 row), and the model is built once per command.

### 4.4 Where the value travels

```
cli.py  mesh / catchment / station-catchments
  │  threads: int | None = None   (--threads, min=0)
  │  memory_budget: str | None = None   (--memory-budget)
  │  res = resources(threads, memory_budget, os.environ)      ← built once, here only
  │
  ├─ _open_dem(..., res) ─ open_dem(request, resources=res)
  │     ├─ repository_for(sources, nodata, cached, threads=res.threads)
  │     │     TiffDemRepository / CacheRepository keep `threads`, pass it to
  │     │     decode_window / decode_class_window            (pool: io/cog.py)
  │     ├─ plan_mosaic(..., memory_budget=res.memory_budget)  (cap: mosaic.py)
  │     └─ _open_reprojected(..., res): resample(grid, source, res.threads),
  │           check_point_blocks(grid, source, domain, res.threads)   (pool: target_grid.py)
  │           and, inside those workers, reprojector(..., threads=1)  (no pool in a pool)
  ├─ _open_features(...): FeatureRequest(threads=res.threads) → _Tally → reprojector(threads)
  ├─ _dem_mesh(..., res):
  │     SlopeTolerance(view, *ramp, threads=res.threads)              (_core)
  │     refine(..., threads=res.threads)                               (_core)
  │     edge_strip.run(..., threads=res.threads) → refine_strip(threads=)   (_core)
  │     final_check.run(..., threads=res.threads) → refine_points(threads=) (_core)
  ├─ _class_fractions(..., res.threads): CacheRepository(threads=), assign(threads=)
  │     → index_space(..., threads) → reprojector(threads); class_areas(threads=)  (_core)
  ├─ record: threads=res.threads  → "Threads used" row in --record and --stats
  └─ _write_report(..., threads=res.threads) → Report.threads → "Threads used: N."

cli.py  catchment / station-catchments
  ├─ repository_for(dem, threads=res.threads)
  ├─ delineate(request, repository, resources=res)
  │     └─ _grow(..., budget=res.memory_budget): _plan → plan_mosaic(memory_budget=);
  │           `need > budget` refuses before the flood (today's check, the budget in place of half the memory)
  └─ run_batch(..., resources=res) → delineate(..., resources=res) per station
```

`delineate` and `open_dem` take `resources: Resources | None = None`; None is
`Resources.fallback()` (one thread, half the machine's memory), so every test
and library caller that passes nothing keeps working and reads no core count.
`Resources` travels whole through these two so a later field (say, an I/O
concurrency for `fetch`) changes no signature; the leaf functions take the one
integer they use.

The stats line changes from "Threads: 10 (hardware concurrency)."
(`src_python/tin_engine/stats.py@5d2e706d:242`) to "Threads used: 3." and is printed for every DEM
run, not only refined ones (`src_python/tin_engine/cli.py@5d2e706d:1590`'s condition goes). The run
record gains the entry `threads` with the wording "Threads used", an Inputs
row (`run_record.py`'s `WORDING` and `INPUTS`), set by the refined and the
stride records; the gallery's flat record has no DEM and no row.

### 4.5 What `--threads` does not reach, and why

- **`fetch` and `fetch-stations`**: I/O-bound; `fetch` already has
  `--connections` for its HTTP concurrency (`src_python/tin_engine/fetch/run.py@5d2e706d:95`) and
  decodes no block (it stores them; `tifffile` is used for page headers
  only). No option.
- **`draw`, `palette`, `version`, `mcp`**: no pool, no `_core` call that
  takes threads.
- **The MCP server (increment 37)**: its schemas are derived from the Typer
  parameters (`src_python/tin_engine/mcp/schema.py@5d2e706d:61-77`), so `threads` and
  `memory_budget` appear in the `mesh`, `catchment` and `station_catchments`
  tools with no server change (one line in `mcp/schema.py` keeps a zero lower bound, section 12), and an agent sets them per call; the child
  inherits the server's environment (`src_python/tin_engine/mcp/server.py@5d2e706d:111`,
  `env={**os.environ, **CHILD_ENV}`), so `RASPUTIN_THREADS` set where the
  server starts is the default for every call. One test pins the first (T9).
- **`tools/bench.py`**: its child still monkeypatches `refine`'s `threads`
  keyword (`tools/bench.py@5d2e706d:788-790`), which overrides the CLI's value
  for `refine` only, so the benchmark's timings are unchanged. Two things
  follow and are stated so nobody reads them as defects: the child's stats
  row reads "Threads used: 10." (the CLI's resolved count, every core on this
  machine) whatever the sweep's count, because the row reports what the CLI
  resolved, not what the monkeypatch passed; and the child's other `_core`
  calls and pools run on the CLI's count. Switching the child to `--threads`
  is a tooling idea for later (section 11).

### 4.6 `decompose.py` goes

`src_python/tin_engine/decompose.py` and `tests/python/test_decompose_partition.py`
are deleted; `src_python/tin_engine/io/mesh_index.py@5d2e706d:37`'s docstring says "as the partition
rule of increment 23c gave them" instead of naming the module;
`project_structure.md`'s `decompose.py` row (`project_structure.md@5d2e706d:170`)
goes and a `resources.py` row comes ("the one module that reads the machine:
cores, memory, cgroup limits; Pydantic and os only"). `PartitionRecord` in
`mesh_index.py` keeps its `memory_budget` field: it records what a cut run was
given, and is read by nothing that needs the deleted rule. Ola's Q4 in the
audit, default delete, taken.

### 4.7 Documentation in the same PR

README's "Quick example" paragraph gets one sentence on `--threads` and
`RASPUTIN_THREADS`; INSTALL's environment paragraph (`INSTALL.md@5d2e706d:188`,
where `RASPUTIN_DATA` is introduced) names `RASPUTIN_THREADS` beside it and
`--memory-budget` with the container rule. Not counted lines.

## 5. Boundaries checked

- *Data versus execution*: `Resources` is data, resolved once at the CLI; the
  algorithms take integers. The C++ sees a plain `unsigned`; nothing of the
  model, the environment or the cgroup crosses the binding (I/O boundary
  intact; no path reaches C++).
- *State*: no module-level machine read remains (`crs.WORKERS` at import time
  was one); `resources.py` has no state, only functions and a frozen model.
- *Dependency gravity*: none added; `os`, `pathlib`, `re`, Pydantic.
- *Async-readiness*: `catchment_batch` already runs `delineate` in
  `asyncio.to_thread`; a `Resources` per call lets an API worker give each job
  its share.
- *Testability*: `machine_cores(cgroup)` and `machine_memory(cgroup)` take the
  cgroup root, so a `tmp_path` with a synthetic `memory.max` tests the read
  without a container; `resources(...)` takes `environ`; every pool takes its
  count; no monkeypatch of `os` is needed anywhere in the new suite, and the
  catchment cap test (`tests/python/test_catchment.py@5d2e706d:296-304`) drops
  its two monkeypatches for a `Resources(memory_budget=...)`.

## 6. Degenerate inputs and refusals

- `--threads 0` on a machine whose count cannot be read (`os.cpu_count()`
  None): 1.
- `cpu.max` reading `max 100000`: no quota. `cpu.max` reading `50000 100000`:
  a quota of 0.5 CPU, `ceil` gives 1. `150000 100000`: 2. Malformed, by the
  memory rule, all no quota: `abc 100000` (first word not a number),
  `200000 0` (zero period), `200000` (one word), `200000 abc`, an empty file,
  an unreadable file.
- `memory.max` reading `max`: no limit; reading `0` or `-1` or text: ignored
  (no limit), never a budget of 0.
- `--memory-budget 1` (one byte): accepted by the parser; every window is then
  refused with the window's size and the budget in the message, which is the
  honest answer to that request.
- A window exactly at the budget is accepted (`need > budget` refuses, as
  today's `>`).
- `RASPUTIN_THREADS=" 4 "`: stripped, 4. `RASPUTIN_THREADS=""`: as unset.

## 7. Left out, with reasons

- A cgroup v1 CPU quota (`cpu/cpu.cfs_quota_us` and `cpu.cfs_period_us`): two
  files and a second format for hosts on the old hierarchy; the memory limit
  is read on both because the task names both and it is one file each.
- Ancestor cgroups on a bare host (walking `/proc/self/cgroup`): the
  container case, which is D6's, mounts the container's cgroup at the root.
- A per-phase thread count (decode on N, refine on M): one number is what Ola
  asked for; `Resources` can grow a field.
- `--threads` on `fetch` (4.5) and on the C++ default (the C++ is untouched).
- Reporting the machine's core count beside the count used ("3 of 10"): one
  more `Report` field for a figure the user can get from `nproc`; left out to
  keep the row plain (Ola's ruling, section 10).
- The quick check's `cases.toml` `threads` lists and `bench.py`'s child keep
  their own `--threads` plumbing (4.5).

## 8. Tests (`@tester`), red first

Fixtures: the micro DEMs of `tests/python/cli_driver.py` (`rough_dem(16)`, a
16 x 16-node GeoTIFF; the geographic one of `test_cli_mesh_geographic.py`,
`tests/fixtures/velhas/anadem_velhas.tif`, for the resampled path) and the
land-cover fraction fixture of `test_cli_mesh_land_cover_fractions.py`. Pools
are observed as `tests/python/test_io_cog.py@5d2e706d:268-299` observes one: the pool class
monkeypatched in the module under test, recording `max_workers`. `_core` calls
are observed through a spy that wraps the real function in the module that
calls it (`cli.refine`, `cli.SlopeTolerance`, `edge_strip.refine_strip`,
`final_check.refine_points`, `land_cover_fractions._core.class_areas`), records
`kwargs["threads"]`, and calls through, so the run still produces its mesh.
Figures on this machine (10 cores, 32 GiB) are given where a test's expectation
depends on the machine; the test computes them from `resources` itself, never
from a literal.

- **T1 one module reads the machine** (`tests/python/test_resources.py`):
  `grep`-equivalent in Python over `src_python/tin_engine/**/*.py`: the
  tokens `cpu_count`, `process_cpu_count`, `sched_getaffinity`, `sysconf`
  and `/sys/fs/cgroup` occur in `resources.py` only. Red today on
  `src_python/tin_engine/crs.py@5d2e706d:39`,
  `src_python/tin_engine/dem_input.py@5d2e706d:218`,
  `src_python/tin_engine/io/cog.py@5d2e706d:203`,
  `src_python/tin_engine/cli.py@5d2e706d:1590` and
  `src_python/tin_engine/mosaic.py@5d2e706d:57`.
- **T2 the threads rule** (section 4.2's table, one case each): `resources(None,
  None, {})` is `machine_cores()`; `(0, None, {})` likewise; `(3, None,
  {"RASPUTIN_THREADS": "7"})` is 3; `(None, None, {"RASPUTIN_THREADS": "7"})`
  is 7; `"0"` is `machine_cores()`; `" 4 "` is 4; `""` is as unset; `"abc"`
  and `"-1"` raise `ValueError` naming `RASPUTIN_THREADS` and the text; 64 on
  a 10-core machine is 64.
- **T3 `machine_cores` with a synthetic `cpu.max`**: `tmp_path/cpu.max` =
  `"200000 100000\n"` gives `min(2, count)` (2 here); `"max 100000"` and an
  absent file give the plain count (10 here); `"50000 100000"` gives 1; the
  malformed files of section 6 (`"abc 100000"`, `"200000 0"`, `"200000"`,
  `"200000 abc"`, an empty file) each give the plain count, never an
  exception.
- **T4 `machine_memory` with synthetic files**: `tmp_path/memory.max` =
  `"1073741824\n"` gives 1 GiB on this 32 GiB machine; `"max"` and an absent
  file give `os.sysconf` physical memory (34359738368 here); v1
  `tmp_path/memory/memory.limit_in_bytes` = `"9223372036854771712"` gives the
  physical figure, `"536870912"` gives 512 MiB; v2 wins when both exist.
- **T5 `parse_bytes`**: `"8G"` = 8589934592, `"512m"` = 536870912, `"1T"`,
  `"4096"`; `"0"`, `"8GB"`, `"-1"`, `"abc"`, `""` raise naming the text. The
  default budget is `machine_memory() // 2` (17179869184 here).
- **T6 the CLI option and the environment** (`test_cli_mesh_stats.py` or a new
  `test_cli_threads.py`): `mesh --dem rough16 --tolerance 1 --threads 3 --stats -`
  prints "Threads used: 3."; with `RASPUTIN_THREADS=5` in the runner's env and
  no option, "Threads used: 5."; with both, the option's 3; `--threads -1` is
  a usage error naming `--threads`; `RASPUTIN_THREADS=abc` is a usage error
  naming `RASPUTIN_THREADS`; a stride run (no `--tolerance`) prints the row
  too; `--record` carries `"threads": 3` with the wording "Threads used" among
  the inputs. `catchment --threads 3` and `station-catchments --threads 3`
  accept the option (the help lists it).
- **T7 every `_core` call site receives the value** (five tests, one per
  site, each a spy): with `--threads 3`, `refine` and `refine_strip` get
  `threads=3` on `rough_dem(16) --tolerance 1`; `refine_points` on the
  resampled Velhas run (`--out-crs` path, which is the `final_check` route);
  `SlopeTolerance` on a `--tolerance-slope` run; `class_areas` on the
  land-cover fractions fixture. Each test asserts the mesh is still written
  (the spy calls through).
- **T8 every pool receives the value** (five tests plus two defaults, one
  spy per large reprojector site of section 3): the
  block decoder's pool has `max_workers=3` on a `--threads 3` cached-source
  run (`test_cli_mesh_cache.py`'s fixture) and on a tile run; the target
  grid's `_pool` has `max_workers=3` on the Velhas `--out-crs` run, and the
  reprojector called inside its workers has `threads=1` (spy on
  `target_grid.reprojector`, which sees both the `:181` and the `:239`
  call); the feature reprojector built in `_Tally._take`
  (`src_python/tin_engine/feature_input.py@5d2e706d:339`) is called with
  `threads=3` on a `--threads 3` run of the reprojected-source case
  (`tests/python/test_cli_mesh_features.py@5d2e706d:328-338`, a GeoJSON in
  EPSG:4326 over the EPSG:25833 DEM, the one feature fixture whose CRS
  differs from the DEM's; spy on `feature_input.reprojector`, which also sees
  the hull call at `:183`, so the assertion picks the call whose source CRS
  is the feature file's); the class-fraction reprojector
  (`src_python/tin_engine/land_cover_fractions.py@5d2e706d:67`) is called
  with `threads=3` on a `--threads 3` run of the land-cover fractions fixture
  (spy on `land_cover_fractions.reprojector`; it is built on every fractions
  run, the class raster being EPSG:4326). These two are the sites a dropped
  `threads` would serialise unseen under the library default of one worker,
  which is why each has its own spy. The reprojector's pool has `max_workers=9` for
  `reprojector(a, b, threads=9)` over 2 250 000 points (36 MB of float64
  pairs; 9 chunks of `CHUNK_POINTS`), where today's `WORKERS` gives 8, and
  `max_workers=2` with `threads=9` over 500 000 points (two chunks).
  Defaults: `reprojector(a, b)` over 2 250 000 points opens no pool
  (`threads=1`: `one()` is called once); `decode_window` without `threads`
  builds its pool with 1 worker (replacing
  `tests/python/test_io_cog.py@5d2e706d:257-262`'s "the machine's cores").
  `tests/python/test_crs.py@5d2e706d:207-213`'s `WORKERS` assertion goes.
- **T9 the MCP schema carries the options**: the `mesh` entry of `derive(app)` has a model with
  fields `threads` (int or None, default None) and `memory_budget` (str or
  None), with one line changed in `mcp/schema.py` (a zero lower bound kept; section 12); `DerivedCommand.argv(model)`
  (`src_python/tin_engine/mcp/schema.py@5d2e706d:108`) with the entry's model
  built from `{"threads": 3}` plus the required arguments ends with
  `--threads 3`, and with `{"memory_budget": "8G"}` ends with
  `--memory-budget 8G`.
- **T10 the memory budget refuses a window above it**, on
  `tests/python/test_catchment.py@5d2e706d:296-304`'s fixture (`bowl()`
  through `repository_of`, `request(api)`). Its figures, measured by running
  `delineate` on this tree with `physical_memory` faked to twice each budget
  (so the cap is the budget), not read from the test's comment, which says
  48 x 48: the run floods three windows, 49 x 49, 113 x 95 and 224 x 175
  nodes (float32, 4 bytes a node for the mosaic, 6 to flood); after each
  flood the growth check plans a probe window at the base margin
  (`src_python/tin_engine/catchment.py@5d2e706d:279`), which `plan_mosaic`'s
  cap judges too. Measured:

  ```
  budget   8192  refused by the mosaic cap:  "the request is 49 x 49 nodes, 9604 bytes at float32, ..."
  budget   9604  refused by the flood cap:   "the window is 49 x 49 nodes, 14406 bytes to flood, ..."
  budget  13824  refused by the flood cap:   the same 14406
  budget  14405  refused by the flood cap:   the same 14406
  budget  14406  first window floods; refused by the mosaic cap on the probe window: "73 x 55 nodes, 16060 bytes"
  budget 235199  refused by the flood cap:   "the window is 224 x 175 nodes, 235200 bytes to flood, ..."
  budget 235200  floods, three windows
  ```

  The tests: `delineate(request, repository, resources=Resources(threads=1,
  memory_budget=12000))` raises `CatchmentError` whose text names
  "49 x 49", "14406 bytes to flood" and "12000" (the flood cap, not the
  mosaic cap, which 12 000 passes); `memory_budget=14406` raises naming
  "73 x 55" (so the first window was accepted at exactly its flood figure);
  `memory_budget=235199` raises naming "224 x 175" and "235200";
  `memory_budget=235200` floods and returns three windows (exactly at the
  budget is accepted). The test keeps the budgets as literals with the
  measured figures beside them, since they are the fixture's, not the
  machine's. The CLI: `catchment ... --memory-budget 8K` (8 192 bytes) is the
  mosaic refusal naming "49 x 49", "9604" and "8192" as a usage error;
  `--memory-budget 8KB` names the text. The mosaic cap on its own:
  `plan_mosaic(footprints, memory_budget=1)` on `test_mosaic.py`'s 2 x 2
  lattice refuses naming `--memory-budget`;
  `tests/python/test_mosaic.py@5d2e706d:658-666` (M5, `physical_memory`)
  becomes T4's physical case in `test_resources.py`.
- **T11 `decompose` is gone**: `importlib.util.find_spec("tin_engine.decompose")`
  is None; `project_structure.md` has no line containing `decompose.py` and
  one containing `resources.py`; `tests/python/test_decompose_partition.py`
  is deleted in the red commit. The layering table
  (`tests/python/test_layering.py@5d2e706d:44`) has a `"decompose": ""` row
  in L1, and its `test_the_table_names_every_module` (`:152-153`) fails on a
  row for a gone module (`:153`, "rows for modules that are gone"); the red
  commit drops that row and adds `"resources": ""` to L1 (`os` and Pydantic
  only, no first-party import), so the table test is red on `resources`
  being absent (`:152`, "modules with no row") and green once the module
  exists.
- **T12 wording**: `stats.render` with `threads=3` prints "Threads used: 3."
  and never "hardware concurrency" (`tests/python/test_stats.py@5d2e706d:246` and `:370`
  updated).

No mutation round: the invariant-critical surface is five keyword arguments
and two file reads, whose failure mode is a typo T7, T8 and T4 catch.

## 9. Estimate, speed, acceptance

**Lines** (production, `CLAUDE.md` section 2's count; added in the new
revision, removed in the old):

| Where | Added | Removed |
|---|---|---|
| `resources.py` (model, three readers, parser, resolution) | 55 | 0 |
| `cli.py` (two `Annotated` option aliases, three signatures, `resources()` once per command, threading through `_open_dem`, `_dem_mesh`, `_class_fractions`, `_write_report`, `catchment`, `station_catchments`, the record entry) | 32 | 4 |
| `dem_input.py`, `io/repository.py` (`resources` keyword, `threads` kept by two repositories, `repository_for`) | 16 | 2 |
| `io/cog.py`, `crs.py`, `target_grid.py` (defaults of 1, `WORKERS` gone, inner reprojector serial) | 8 | 6 |
| `feature_input.py`, `land_cover_fractions.py`, `edge_strip.py`, `final_check.py` | 12 | 2 |
| `mosaic.py`, `catchment.py`, `catchment_batch.py` (budget keyword, `physical_memory` gone, refusal wording) | 14 | 8 |
| `stats.py`, `run_record.py` | 3 | 1 |
| `decompose.py` deleted | 0 | 55 |
| **Total** | **140** | **78** |

Net about **60**; the audit said about 130 (210 added, 80 removed) before the
call sites were counted. Worst case at +39 % on the additions: 195 added, net
about 120. Both far under 700. One PR. Unmoved by review round 1: the
`cpu.max` rule (4.1) is one more condition in a reader already counted in
`resources.py`'s 55, the wording of two refusals changes no line count, and
the rest of the round's fixes are tests and prose.

**Speed, one line:** the refine path is touched only by passing a value it
already takes (today's implicit 0 resolves to the same `hardware_concurrency`
count that `machine_cores()` returns on this machine, 10 = 10), so every mesh
is bit-identical and no timing change is expected, and `python3
tools/bench_quick.py run` on the branch head against the stored baseline
confirms it; the one deliberate change in work scheduling, the reprojector
inside the target grid's workers running serial instead of opening a second
pool, can only shorten the `resample` phase of the Lagan GLO-30 case, and the
quick check's `decode`/`resample` rows are where it would show.

**Acceptance by the rule:** no file under `include/terrain/refinement/`,
`include/terrain/mesh/` or `bindings/` changes, so no `@perf` acceptance run;
the quick check is the gate.

## 10. Ola's rulings

Four questions were put to Ola with defaults; he ruled "Yes to all" on
2026-10-12, taking each default. No question is open.

1. **`mesh` gets `--memory-budget` too**, beside `catchment` and
   `station-catchments` (its mosaic cap reads the host's memory in a
   container as `catchment`'s does). Section 4.3.
2. **A thread count above the machine's cores is honoured, not clamped**
   (the benchmark's sweep to 20 on 10 cores, OpenMP's and Rayon's practice).
   Section 4.2.
3. **A container's CPU quota (`cpu.max`) counts toward "every core"**, as its
   memory limit counts toward the memory default. Section 4.1.
4. **The stats row reads "Threads used: 3."**, without the machine's count.
   Section 4.4; the longer row stays left out (section 7).

## 11. The ideas file

Entries this design takes up or overlaps: *Product review (b), generality
audit* (this is its first fix); *MCP server over the rasputin commands* (an
agent now sets `threads` per call through the derived schema, and
`RASPUTIN_THREADS` in the server's environment is the default; `mcp/schema.py`
changes by one line, section 12). Every other entry is untouched. New idea from this design,
for the handback: `tools/bench.py`'s child could pass `--threads` to the CLI
and drop its monkeypatch of `refine` (tooling, not counted).

## 12. As built (recorded by the main session from the handbacks)

Ola picked this as the audit's first fix on 2026-10-12 ("1, 2: yes") and
ruled the four design questions the same day ("Yes to all, but thousands of
sub-catchments is quite extensive", section 10). The design took two review
rounds; round 1's five blockers were figures and call sites, settled by
probing the fixture through a whole `delineate` run rather than its first
window.

**Red** `4b9c0dfe` (`@tester`): 110 tests, 85 in the new
`tests/python/test_resources.py` and `tests/python/test_cli_threads.py`
(the option and its environment default, the model reaching every call site
through one spy per site, the cgroup reads on synthetic files with the
malformed cases, the memory budget with the probed figures, the stats row),
25 in nine edited files; `test_decompose_partition.py` deleted. Every red
reason a missing thing: the module `tin_engine.resources`, the two options,
a `plan_mosaic` keyword, a `reprojector` keyword, the old wording, the 8- and
10-worker defaults, the layering row. Pins beyond the design, asserted:
`RASPUTIN_THREADS` "1.5" and "0x10" refused as not whole numbers; a negative
`cpu.max` quota means no quota; `parse_bytes` refuses "G" and "1.5G", reads
"2k" as 2048; the three commands' help names `RASPUTIN_THREADS`; `threads` is
a record entry and an Inputs row, not a mesh-file field; the flood and
mosaic refusals name `--memory-budget` and the figure. T11's assertion on
`project_structure.md` left out: a prose file, forbidden to tests by
`conftest.py`; the reviewer reads the rows. One citation into the deleted
test file (`docs/increments/23-basin-scale.md`, a 23c-1 review record) was
pinned at `5d2e706d` by the main session (`3b709b6f`).

**Green** `b252a70d` (`@developer`): `src_python/tin_engine/resources.py`
new (the `Resources` model: threads from `--threads`, 0 meaning all cores,
`RASPUTIN_THREADS` as the default when the option is absent, the machine's
cores otherwise with a well-formed `cpu.max` quota counted; the memory
budget from `--memory-budget` with `memory.max` read when present);
`--threads` and `--memory-budget` on `mesh`, `catchment` and
`station-catchments`; the count reaching all five `_core` entry points and
the three Python pools, the reprojector's cap of 8 gone and the pools'
library defaults set to one worker so only `resources.py` reads the machine;
both memory caps overridable by the budget; the stats row "Threads used: N."
and the `threads` record entry; `decompose.py` deleted; README, INSTALL and
`project_structure.md` updated. 96 net production lines (237 added, 141
removed; the overrun over 60 is `cli.py`'s option aliases and `_resources`
helper and the two repositories' `threads`). Gates clean: mypy, ruff,
prohibited deps, detria boundary, check_citations. Two fixes after the first
full run: the `--threads` metavar `N`, since `<int range>` made Rich truncate
long option names in `--help`; `_domain_plan`'s budget defaulting to None.
Pins beyond the design: `mcp/schema.py` changed by one line (`minimum: 0` no
longer dropped by `if v`), required by T9 although section 4.5 said no `mcp/`
change (a latent bug for any option with a zero bound; `@architect`'s
one-line correction of 4.5 is due); the `_core` wrappers `edge_strip.run` and
`final_check.run` default `threads=0` with the C++ meaning while the Python
pools default 1, not pinned by a test (`@tester` recommends 1, matching
section 4.2; default kept at 0, documented in the two docstrings; for the
code review to judge).

**Test-side fixes** `9de3212c` (`@tester`, seven files): the five layering
rows gain `resources` (each module grew by exactly `tin_engine.resources`,
checked with `importscan`); `as_16_opened_it` passes `resources` through;
the plain-output bit-identity test runs `--threads 1` and `--threads 2` to
the same paths (equal mesh bytes, records differing only in `threads` and
the echoed command); `plan_mosaic`'s parameter list gains `memory_budget`;
the feature-input fakes and the spied reprojector accept `threads=`; T9's
docstring no longer claims no `mcp/` change; no scaffolding left. Full
`pytest tests/python/` with the coverage gate: 7,422 passed, 18 skipped,
coverage 98.44 %, one failure beside the night's detached mesh batch
(`test_core_cdt.py`'s GIL ticker, 7 ticks against a bound of 3, load-bound;
5 of 5 alone, untouched by this increment). `ruff` clean on `tests/`.

**Speed:** no refine or mesh code touched. `@perf`'s quick check of this
branch at `8f9289c4` (2026-10-12 18:31 UTC, idle, AC, against the baseline
refreshed minutes before on master `f09b9863`; evidence in
`docs/benchmarks/2026-10-12/baseline-refresh/README.md`, evening section, on
`worktree-baseline-1012`): **NO CHANGE**, meshes identical on every case, no
phase outside its band; ACCEPTED.

**Lessons for the next retrospective:** a design that deletes a file lists
the unpinned citations into it first (`check_citations` would have shown the
23c-1 record); a red step that adds a keyword to a shared function widens
the suite's existing fakes of it in the same commit; the editable install's
loader keeps a deleted module's `find_spec` until reinstall; the GIL ticker
test's bound of 3 ticks is load-sensitive.

## Review

**Design review, round 1, 2026-10-12.** Range `5d2e706d..9e0d60c2`, 0 production lines, estimate 60 net with worst case 120, @reviewer: CHANGES REQUESTED, five blockers. (1) /Users/skavhaug/projects/rasputin/.claude/worktrees/threads/docs/increments/38-threads-and-memory.md@9e0d60c2:476-481: T10's figures are wrong; probed on the fixture, the first window is 49 x 49, the mosaic cap 9,604 bytes and the flood 14,406, so a budget of 8,192 is refused by the mosaic cap (`/Users/skavhaug/projects/rasputin/.claude/worktrees/threads/src_python/tin_engine/catchment.py@5d2e706d:394`) naming 9,604, not the flood, and 13,824 also refuses; fix: a budget in 9,604 to 14,405 for the flood case, 14,406 for "exactly at", the `8K` CLI case restated. (2) @9e0d60c2:485-488: T11 misses `/Users/skavhaug/projects/rasputin/.claude/worktrees/threads/tests/python/test_layering.py@5d2e706d:44` (the `decompose` row), whose line 152 fails on a gone module; the red commit drops that row and adds `resources`. (3) @9e0d60c2:457-470: T8 spies no reprojector at `/Users/skavhaug/projects/rasputin/.claude/worktrees/threads/src_python/tin_engine/feature_input.py@5d2e706d:339` or `/Users/skavhaug/projects/rasputin/.claude/worktrees/threads/src_python/tin_engine/land_cover_fractions.py@5d2e706d:67`, two of the three large sites; with the library default of one worker a dropped `threads` serialises them unseen; one spy per site. (4) @9e0d60c2:133-141: 15 reprojector callers, not 14 (`target_grid.py:136`, `tolerance_field.py:83`, `feature_input.py:183` missing; "other ten" is twelve); `target_grid.py:239`'s blocks are 512 x 512 = 262,144 points, under the 500,000 threshold, so only line 181 nests a pool. (5) @9e0d60c2:177-181: `cpu.max` with a non-numeric first word or a zero period is unspecified; apply the memory rule (no quota). Suggestions: the refusal text "half the machine's memory" (:254-256) is false when a budget was given; section 4.5 should say `bench.py` children's stats row still shows the CLI's count; T9's `argv` is `DerivedCommand.argv(model)` (`src_python/tin_engine/mcp/schema.py@5d2e706d:108`). Verified: all 37 pinned cites; the legacy greps; the I/O boundary, the pinned assumptions and the MCP claims; `check_citations` 46 at-risk lines, all `ROADMAP.md` rows 27 and 29, unmoved. Lesson: T10's figures came from a test comment, not a run; one probe script settled it.

**Design review, round 2, 2026-10-12.** Range `3ff48435..3018192a`, 0 production lines, estimate 60 net, @reviewer: APPROVED. Verified by running: T10's figures with the main venv (8,192 refused by the mosaic cap at 9,604; 9,604 to 14,405 by the flood cap at 14,406; 14,406 refused on the probe window 73 x 55 at 16,060; 235,199 refused naming 224 x 175; 235,200 floods three windows), every figure in the table matching; the `decompose` layering row and its two assertions; the three T8 spy sites and the EPSG:4326-over-25833 case; the caller grep reproducing 15 lines, `target_grid.py` under the 500,000 threshold; the `cpu.max` rule in 4.1 with its cases; section 10's four rulings stated as rulings with date and quote; check_citations exit 0. Idea, non-blocking: T3 omits section 6's unreadable-file case (chmod 000 in `tmp_path`, skipped as root).

**Code review, round 1, 2026-10-12.** Range `5d2e706d..f13d1370`, 96 net production lines (237 added, 141 removed; under 700), @reviewer: CHANGES REQUESTED, two prose blockers, the code sound; applied in this recording commit. (1) /Users/skavhaug/projects/rasputin/.claude/worktrees/threads/INSTALL.md@f13d1370:191-196: the new paragraph sat between a sentence ending in a colon and the code block it introduces, and said "Two more environment variables" naming one; moved below the block, "One more environment variable". (2) /Users/skavhaug/projects/rasputin/.claude/worktrees/threads/docs/increments/38-threads-and-memory.md@f13d1370:374, :552 and :671-672 said no `mcp/` change; the branch changes `mcp/schema.py` by one line; corrected. Pins accepted: the `mcp/schema.py` zero-bound fix; `edge_strip.run` and `final_check.run` defaulting to `threads=0` (section 4.2 keeps 0 for direct callers; the CLI always passes 1 or more). Verified: no scaffolding; the machine read only in `resources.py`; the model's rules, the cgroup reads and `parse_bytes` as sections 4.1 to 4.3 and 6; the count reaching all five C++ calls and three pools, the cap of 8 gone; stats row and record entry; both refusals naming `--memory-budget` and the figure; `decompose.py` gone; the C++ still takes a plain integer; check_citations exit 0; the four new suites 181 passed. Lesson: a paragraph inserted into a document breaks a "sentence, colon, code block" pattern unless the lines around it are read.

**Code review, round 2, 2026-10-12.** Range `f13d1370..7b6d4ddc`, @reviewer: APPROVED. Both round-1 fixes in place (the INSTALL paragraph after the code block, "One more environment variable"; the three design sentences naming the one-line `mcp/schema.py` change, matching the diff); round 1's citations resolve; nothing else changed; check_citations exit 0.

**Merge review, 2026-10-13.** Merge `3942af15` (origin/master with increment 41 into this branch, by `@developer`), @reviewer: APPROVED. The five `cli.py` conflict hunks keep both sides: `_open_features` takes the vocabulary then `threads` and builds `FeatureRequest` with both (/Users/skavhaug/projects/rasputin/.claude/worktrees/threads/src_python/tin_engine/cli.py@3942af15:1955), its one call passes them in that order (/Users/skavhaug/projects/rasputin/.claude/worktrees/threads/src_python/tin_engine/cli.py@3942af15:1168), and `_dem_mesh`'s one call passes `res.threads, vocabulary` in its signature's order; 450 targeted tests passed, 2 skipped (Ola's data).
