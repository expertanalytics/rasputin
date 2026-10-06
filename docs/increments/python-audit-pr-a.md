# PR A design: `audit-lattice` (F2, F9, F10's repository Protocol, F12's two `mosaic` edges)

Status: **code review approved in round 1 (1a84ea9); `@perf` accepted
(01751ec); -27 net production lines. Ready to push on Ola's yes once PR B
(`worktree-audit-crs`) has merged**. Ola's two questions (below) bind it;
their defaults are accept and no change.
Branch `worktree-audit-lattice` was written on PR B's approved head `618328b`
and now sits on B's merged head `ae493da` (B merged with master `fd64f8b`),
by a merge, not a rebase. A's diff outside `docs/` is unchanged by it: `git
diff ae493da HEAD` and `git diff 618328b 883c51f`, both outside `docs/`, are
the same, and `count_loc.py ae493da HEAD` is still -27.
Every `@618328b` citation below reads B's old head, which stays in this
branch's history and, through B's merge, in master's. The audit this designs from is
`docs/increments/python-audit.md` (sections 2, 5, 6 and 7); this file is
separate because several unpushed audit branches edit that one.

## What A takes, measured

A puts the node arithmetic of a DEM grid on `RasterMeta`, once, in the core's
spelling; moves `Bounds` and `TileFootprint` down to `io/models.py`; writes
the NoData rule once; adds `load_window` and `check` to the `DemRepository`
Protocol; and so removes `mosaic`'s two upward imports. Behaviour is
unchanged except where Ola ruled: **+-infinity is not NoData**
(`python-audit.md` section 7, ruling 1).

**About -28 net production lines, not the audit's about -100.** A scratch
prototype of this design (a clone of `618328b`, not committed, removed after)
measures -28 with `python3 tools/count_loc.py 618328b <prototype>`, with the
formatter run, mypy strict clean and ruff clean. Most sites go one line to
one line: the saving is the drift (one spelling of a node, one NoData rule),
not lines. The audit's figure counted the sites, not their rewrites.

Departures from the audit's shape (F2), each with its reason:

- **No `IndexWindow` methods** (`meet`, `within`, `clamp`,
  `from_inclusive`). `mosaic` is their only user, and it spells each once;
  moving them widens a layer-0 value type for one caller and saves no line.
  `mosaic`'s inclusive `(r0, c0, r1, c1)` tuples stay private to it.
- **No `cell_area`.** `delta_x * delta_y` cannot drift.
- **The three rounding rules stay three**, at their callers: nearest node
  (`round`, `np.round`), outward cover (`floor` low, `ceil` high), strictly
  inside (`ceil` low, `floor` high). They answer different questions; what
  becomes one is the fractional index they round.
- **`mosaic`'s lattice-relative arithmetic stays** (`_placed`, `_snapped`,
  `_covering`'s window, `_grid`): it counts from a lattice's reference node,
  not from one tile's corner, and is where the mosaic's global node indices
  live.
- **`burn._floor_node` stays**: it subtracts the search radius before
  dividing (`src_python/tin_engine/burn.py@618328b:71-72`); through
  `index_of` the sum would be formed first, which can move the quotient by
  an ulp, and `floor` or `ceil` can turn that into a node more or less in
  the search window.
- **`fetch/plan.plan_object`'s rule stays** ("a node strictly within one cell
  of the box", `src_python/tin_engine/fetch/plan.py@618328b:151-152`): it
  must fetch a superset of what `mosaic._snapped` reads (F2), and A does not
  touch the rule, so it cannot break that. Only its node coordinates go
  through `node_xy`.
- **`cli`'s two `--bbox` sites stay** (`src_python/tin_engine/cli.py@618328b:1268, 1382`):
  F1 gives `--bbox` to PR G; G then calls the `Bounds.of` added here.
- **`RasterMeta.node_box() -> tuple[float, float, float, float]`, not the
  audit's `node_bounds() -> Bounds`.** `Bounds` refuses a flat box
  (`x_min < x_max` and `y_min < y_max`). A one-row or one-column grid's
  node rectangle is flat, and `domain.check_extent` and `dem_input`'s union
  compute that flat rectangle today (the base probe's line
  `domain.check_extent | one row, domain outside` names `y 6600000.0 ..
  6600000.0`). A caller that wants a `Bounds` writes
  `Bounds.of(m.node_box())` (`catchment`), and so states that it refuses a
  flat grid.
- **F10 here is the repository Protocol and `dem_input`'s two `Any`s only.**
  `cli._placed` and `_reach_crs`'s `repository: Any`
  (`src_python/tin_engine/cli.py@618328b:1706, 1726`) are typed
  `DemRepository` by PR F (`worktree-audit-catchment`), whose design lists
  them (`docs/increments/python-audit.md@e2baa5f:781`) and whose green step
  wrote them (`src_python/tin_engine/cli.py@e2baa5f:1713, 1733`); A leaves
  them, so the two branches do not edit the same lines.
- **Added: `TargetGrid.node_box`**, the target grid's own rectangle, spelt
  twice in `target_grid.py` (`src_python/tin_engine/target_grid.py@618328b:117-118, 232-233`).

## Prior art: legacy and literature

*Literature.* None applies: A claims nothing new. The spelling it adopts is
the C++ core's, which the Python side must match bit for bit:
`RasterGeometry::node`
(`include/terrain/raster/geometry.hpp@618328b:69-72`, `x_min + col *
delta_x`, `y_max - row * delta_y`), the inverse in `CheckPoints::add`
(`include/terrain/refinement/check_points.hpp@618328b:82`), and
`Raster::is_nodata` (`include/terrain/raster/raster.hpp@618328b:63-66`: NaN
or the sentinel).

*Legacy.* `git grep -nE "is_nodata|isinf|isfinite|nodata" legacy-archive --
legacy` returns nothing: the legacy had no NoData rule. `git grep -nE "x_min
?\+|delta_x|dx ?\*|x0 ?\+" legacy-archive -- legacy` returns
`legacy/bindings.cpp` (151, 155, 157, 170), `legacy/rasputin/reader.py`
(302, 309, 321, 344, 354, 362, 364, 376, 377, 387, 392, 395, 399, 426) and
`legacy/rasputin/triangulate_dem.h` (337, 340). The one shape worth noting:
`legacy/rasputin/reader.py:309` puts `x_max` on the extents value type, as a method,
`self.x_min + self.delta_x*(self.shape[1] - 1)`, the direction A takes.
Nothing to carry across.

## The types (layer 0, `src_python/tin_engine/io/models.py`)

```python
class Bounds(BaseModel):          # moved unchanged from mosaic.py, plus `of`
    @classmethod
    def of(cls, box: tuple[float, float, float, float]) -> Self:
        """The box `(x_min, y_min, x_max, y_max)`, shapely's `bounds` order."""

class RasterMeta(BaseModel):      # four methods added; fields unchanged
    def node_xy[T: (float, npt.NDArray[Any])](self, rows: T, cols: T) -> tuple[T, T]:
        """Node `(row, col)`'s `(x, y)`, elementwise: `x_min + col * delta_x`,
        `y_max - row * delta_y`, the core's `RasterGeometry::node`."""
    def index_of[T: (float, npt.NDArray[Any])](self, x: T, y: T) -> tuple[T, T]:
        """The fractional `(row, col)` of `(x, y)`, unrounded:
        `(y_max - y) / delta_y`, `(x - x_min) / delta_x`."""
    def node_box(self) -> tuple[float, float, float, float]:
        """The node rectangle `(x_min, y_min, x_max, y_max)`; flat for one row
        or one column (not a `Bounds`, which refuses a flat box)."""
    def windowed(self, window: IndexWindow) -> RasterMeta:
        """This grid cut to `window`: the corner moved by whole cells."""

def valid_mask(values: npt.NDArray[Any], nodata: float | None) -> npt.NDArray[np.bool_]:
    """Not NoData. NoData is NaN or the sentinel; +-inf is data, as in the
    core's `Raster::is_nodata`."""

class TileFootprint(BaseModel):   # moved unchanged from io/repository.py
```

Rules the code must keep, because meshes depend on them bit for bit:

- `node_xy` is exactly `self.x_min + cols * self.delta_x, self.y_max - rows *
  self.delta_y`: no sum of steps, no reordering, no cast. `x` reads only
  `cols`, `y` only `rows`; there is no broadcasting between them, so a caller
  may pass a column range and a row band of different lengths
  (`reference.nodes_inside`, `fetch/plan.plan_object`). A caller passes what
  it passed before (Python ints, numpy integer arrays, float arrays), so the
  result's type is what it was (a `np.float64` stays one: a message that
  prints it would otherwise change).
- `node_box` is built from `node_xy(rows - 1, cols - 1)`;
  `windowed` from `node_xy(window.row0, window.col0)`, which is
  `io.cog.window_meta`'s expression.
- `valid_mask` adds no cast: each caller passes the array it compared before
  (`target_grid` its float64 copies, `mosaic` and `burn` the tile's dtype).
- The constrained type variable passes mypy strict on a probe file
  (`node_xy(1, 2)` is `tuple[float, float]`, `node_xy(np.arange(3),
  np.arange(4))` a pair of arrays) and in the prototype. If mypy refuses it
  in the real tree, `@overload` pairs are the fallback, about 8 lines more.

`TargetGrid.node_box() -> tuple[float, float, float, float]`, in
`target_grid.py`, keeps today's spelling: `x0, y1 = col0 * h, -row0 * h`,
then `x0 + (cols - 1) * h`, `y1 - (rows - 1) * h`.

`DemRepository` (`io/repository.py`) gains the two methods every caller
already uses, with the signatures both implementations have:

```python
def load_window(self, name: str, window: IndexWindow) -> DemTile: ...
def check(self, plan: MosaicPlan) -> None: ...   # MosaicPlan under TYPE_CHECKING, as now
```

`io/repository.py` imports `TileFootprint` from `.models` and drops it from
its `__all__`; `io/cog.py` loses `window_meta` and its `__all__` entry.

## The sites, as they change

Every node, index, rectangle and NoData expression in `src_python/`, at
`618328b`. "Same bits" means the new call evaluates the same floating-point
operations in the same order on the same operands.

| Site | Today | After |
|---|---|---|
| `src_python/tin_engine/mosaic.py@618328b:64-81` | `Bounds` | moved to `io/models.py`; `mosaic` imports it |
| `src_python/tin_engine/mosaic.py@618328b:30, 33-34` | `window_meta` from `io.cog`; `TileFootprint` from `io.repository` | both edges go; `mosaic` imports `io.models` only |
| `src_python/tin_engine/mosaic.py@618328b:265` | `assemble`'s docstring: "its meta must be `window_meta` of the listed one" | "its meta must be the listed one cut to that window, `placement.meta.windowed(placement.source)`" |
| `src_python/tin_engine/mosaic.py@618328b:376, 378` | node rectangle by hand | `m.node_box()[2]`, `[1]` |
| `src_python/tin_engine/mosaic.py@618328b:457-459, 466-467` | node coordinates | `meta.node_xy(rows, cols)`; the hit rows' first and last as one call |
| `src_python/tin_engine/mosaic.py@618328b:483` | `window_meta(placement.meta, placement.source)` | `placement.meta.windowed(placement.source)` |
| `src_python/tin_engine/mosaic.py@618328b:506-510` | `_valid` | deleted; `valid_mask` at its two callers |
| `src_python/tin_engine/mosaic.py@618328b:563` | node coordinates | `np.broadcast_arrays(*meta.node_xy(rows, cols))` |
| `src_python/tin_engine/io/cog.py@618328b:104-113, 163, 174` | `window_meta` | deleted; `meta.windowed(w)` |
| `src_python/tin_engine/io/repository.py@618328b:53-61, 64-74, 378` | `TileFootprint` here; Protocol of two methods | moved; four methods |
| `src_python/tin_engine/target_grid.py@618328b:29` | `Bounds` from `mosaic` | from `io.models`: the `target_grid -> mosaic` edge goes |
| `src_python/tin_engine/target_grid.py@618328b:117-118, 232-233` | the grid's rectangle, twice | `grid.node_box()` |
| `src_python/tin_engine/target_grid.py@618328b:153-155` | `_valid`, `isfinite`: +-inf is NoData | deleted; `valid_mask`: **+-inf is data (the ruling)** |
| `src_python/tin_engine/target_grid.py@618328b:179-180` | fractional index | `m.index_of(lonlat[:, 0], lonlat[:, 1])` |
| `src_python/tin_engine/target_grid.py@618328b:244, 251` | node coordinates | `m.node_xy(rr, cc)`, `m.node_xy(r0 + r, c0 + c)` |
| `src_python/tin_engine/burn.py@618328b:140` | NoData by hand | `valid_mask(raw, m.nodata)` |
| `src_python/tin_engine/burn.py@618328b:151, 169, 213, 223` | index and node coordinates | `m.index_of`, `m.node_xy` (`nx, ny = m.node_xy(*chain[placed])`) |
| `src_python/tin_engine/grid_domain.py@618328b:68-69` | the start nodes | `np.broadcast_arrays(*meta.node_xy(r, c))`; the docstring above names `node_xy` |
| `src_python/tin_engine/domain.py@618328b:146-147` | rectangle | `_, y_min, x_max, _ = meta.node_box()` |
| `src_python/tin_engine/dem_input.py@618328b:178-179` | union of rectangles | `boxes = [m.node_box() for m in metas]`, then the four `min`/`max` |
| `src_python/tin_engine/dem_input.py@618328b:190-192, 227` | `repository: Any`, `footprints: Any` | `DemRepository`, `Sequence[TileFootprint]` |
| `src_python/tin_engine/dem_input.py@618328b:236-238, 261` | `Bounds(**dict(zip(...)))`; rectangle | `Bounds.of(domain.polygon.bounds)`; `m.node_box()` |
| `src_python/tin_engine/catchment.py@618328b:384-385, 391-398` | index, three rounding rules | `m.index_of(...)`, each caller's rounding kept |
| `src_python/tin_engine/catchment.py@618328b:402, 331, 469, 475` | node coordinates | `m.node_xy(...)` |
| `src_python/tin_engine/catchment.py@618328b:407-417, 263, 317, 420-428` | `_extent`, `_bounds_of` | deleted; `Bounds.of(m.node_box())`; `_joined` reads `node_box` and two `node_xy` |
| `src_python/tin_engine/catchment_batch.py@618328b:159-164` | rectangle by hand | `[box(*f.meta.node_box()) for f in repository.footprints()]` |
| `src_python/tin_engine/reference.py@618328b:75-86` | index and node coordinates | `meta.index_of` for the window; `meta.node_xy(band, cols)` per row band |
| `src_python/tin_engine/cli.py@618328b:118, 1694-1701` | `Bounds` from `mosaic`; `_off_node` | from `io.models`; `meta.node_xy(*np.round(meta.index_of(...)))` against the input |
| `src_python/tin_engine/fetch/plan.py@618328b:29, 149-150` | `Bounds` from `mosaic`; node coordinates | from `io.models` (the edge goes); `meta.node_xy(np.arange(rows), np.arange(cols))` |
| `src_python/tin_engine/fetch/run.py@618328b:34` | `Bounds` from `mosaic` | from `io.models` (the edge goes) |

`src_python/tin_engine/fetch/plan.py@618328b:149-150` writes `meta.delta_x * np.arange(...)`, the
operands the other way round; IEEE multiplication commutes, so the bits are
the same, and the probe below agrees.

## The +-infinity ruling: what it changes

Only `target_grid`, the reprojected path. `mosaic` and `burn` already count
+-inf as data; the probe's lines for them do not change.

- **`resample`**: a target node whose 2 x 2 stencil holds an infinity is no
  longer NoData. Its value is the bilinear formula's: +-inf where the
  infinity has a positive weight; NaN where it has weight 0 (0 x inf), or
  where +inf meets -inf. NaN is NoData to the core, so a zero-weight
  infinity still makes that node NoData, where the core's own `bilinear`
  reads a node it lands on exactly and nothing else (increment 27). On the
  probe's 6 x 5 source with one +inf and no sentinel, the 11 x 9 half-cell
  grid goes from 16 NoData nodes (NaN) to 7 NaN and 9 at +inf; the 6 x 5
  on-node grid from 4 NaN to 3 NaN and 1 at +inf. With a sentinel, the
  NoData nodes were the sentinel and are now NaN: the same meaning to the
  core. Question 1.
- **`check_point_blocks`**: a source node at +-inf is now a check point with
  that `z`. The core's store drops it on `add`, counted as `outside`
  (`include/terrain/refinement/check_points.hpp@618328b:84`; probe: three
  points with `z` 1, inf and -inf store 1, `outside` 2). Nothing in Python
  reads `outside`, and `final_check.run` returns `store.size`, so the final
  check, its count and the mesh are unchanged.

## Layering table (`tests/python/test_layering.py`)

| Row | Today | After |
|---|---|---|
| `mosaic` | `io.cog io.models io.repository` | `io.models` |
| `target_grid` | `crs domain io.models mosaic` | `crs domain io.models` |
| `fetch.plan` | `crs domain fetch.http io.cog io.geotiff io.models mosaic` | drops `mosaic` |
| `fetch.run` | `crs fetch.http fetch.plan io.models io.repository mosaic sources tin_engine` | drops `mosaic` (one line: `crs fetch.http fetch.plan io.models io.repository sources tin_engine`) |
| `UPWARD` | `("mosaic", "io.cog")`, `("mosaic", "io.repository")` | both deleted |

`io.repository` keeps `mosaic` (`MosaicPlan`, under `TYPE_CHECKING`, layer
2 to layer 1: down). `cli`, `dem_input` and `catchment` keep `mosaic` for
`Seam`, `plan_mosaic` and the rest. The prototype's table passes
`test_layering.py`.

## The differential probe

`docs/increments/python-audit-probes/lattice_probe.py` reaches every
function A edits, directly or through its caller (18 entries:
`grid_domain.subsample`, `domain.check_extent`, `dem_input._past`,
`_reprojected`, `open_dem` (a box, two domains, and a domain reprojected to
UTM 32, on the committed seam fixture), `reference.nodes_inside`,
`catchment._seed_mask`, `catchment.delineate` (a lake, a reach, points, and
a valley whose window grows), `cli._off_node`, `fetch.plan.plan_object`,
`source_box`, `target_grid.source_region`, `resample`,
`check_point_blocks`, `burn.burn_reach`, `mosaic.plan_mosaic`, `assemble`
whole and windowed, on synthetic tiles and the committed DEM fixtures, and
`catchment_batch.run_batch`). Counted with a profiler hook at the base, the
private ones it reaches through them include `_joined` (11 calls),
`_covering` (18), `_domain_plan` (2), `_open_reprojected` (1), `_gauged`
(1) and `_loaded` (237). It runs over grids (10 m, non-dyadic, one row, one
column, one node, geographic, negative spacing, `None`, a dict) and values
(clean, NaN, +inf, -inf, both, the sentinel; float32 and float64; with and
without a sentinel). One line per pair: the value written exactly, or the
refusal. Its docstring says how to run it on a scratch copy; it imports only
the standard library, numpy, shapely and `tin_engine` (`check_prohibited_deps.py`
does not scan `docs/`). Its output at the base is committed beside it,
`lattice_probe-618328b.txt` (495 lines).

**The green step runs it on a scratch copy of its head and commits the
output beside the base's.** `diff` of the two is the record of changed
behaviour, and must be exactly these 88 lines, the prototype's:

1. **72 lines, the ruling**: `target_grid.resample` (48) and
   `target_grid.check_point_blocks` (24), every case whose values hold +inf,
   -inf, or both (3 value cases x 2 dtypes x 2 sentinels; 4 resample lines
   and 2 check-point lines each). The check filters on the case field (the
   second, `+inf`, `-inf` or `+inf beside -inf`), because every array
   outcome contains "inf" (`+inf=0`), so a whole-line `grep -v inf` passes
   any diff of `resample` or `check_point_blocks`: `diff <base> <head> | awk -F' [|] ' '/^> target_grid/ && $2 !~
   /^[+-]inf/' | wc -l` is 0. It can fail: on the base output, each line
   prefixed `> `, it counts 81 (the `target_grid` lines of the other cases),
   and the same filter inverted, on `resample` and `check_point_blocks`,
   counts the 72.
2. **16 lines, a grid of the wrong type** (`None` or a dict, which no
   command can pass): the `AttributeError` names the method it now calls
   (`node_box`, `index_of`, `node_xy`) instead of the field it read
   (`x_min`, `delta_x`): `domain.check_extent`, `dem_input._past` and
   `reference.nodes_inside` twice each, `cli._off_node` and
   `fetch.plan.plan_object` once, for each wrong type.

Any other changed line is a defect. The probe can fail: with `node_xy`
mutated to `(x_min / delta_x + cols) * delta_x`, an ulp off on non-dyadic
grids, 5 lines change (`grid_domain.subsample` on the non-dyadic and
geographic grids, `cli._off_node` on both). It can also fail on a mutated
`valid_mask`, but not on an ulp-level change to `index_of`: every caller
rounds its result, so the guard for `index_of` is
`tests/python/test_io_models.py`'s
`TestIndexOf::test_bit_for_bit_at_nodes_and_between_them`.

## Red tests (`@tester`, one commit, before any code)

Lean: no mutation round.

1. **`tests/python/test_io_models.py`** (new; red, nothing exists):
   `Bounds` and `TileFootprint` importable from `tin_engine.io.models`;
   `Bounds.of((x0, y0, x1, y1))` equals the keyword form and refuses what
   the constructor refuses; on the non-dyadic grid (origin 0.3, 0.7, spacing
   0.1, 7 x 6), `node_xy` equals `x_min + c * delta_x`, `y_max - r *
   delta_y` bit for bit for Python ints, numpy integer arrays and float
   arrays, and takes a row band and a column range of different lengths;
   `index_of` equals `((y_max - y) / delta_y, (x - x_min) / delta_x)` bit for
   bit, and `np.round(index_of(*node_xy(r, c)))` is `(r, c)` at every node;
   `node_box` on a 6 x 5, a one-row and a one-node grid (flat, no refusal);
   `windowed(w)` equals the `window_meta` expression with every other field
   unchanged; `valid_mask([1, nan, inf, -inf, s], s)` is `[T, F, T, T, F]`
   and with `None` `[T, F, T, T, T]`.
2. **The ruling, `tests/python/test_target_grid.py`** (red: today both are
   NoData or dropped): `resample` of a 6 x 5 source with +inf at one node,
   onto a same-CRS grid at half the spacing: a target node with that node at
   positive weight is +inf, one with it at weight 0 is NaN, one away from it
   unchanged; -inf mirrored. `check_point_blocks` on that source yields the
   +inf node, with `z` +inf.
3. **`DemRepository`** (red): `{"load_window", "check"} <=
   set(vars(DemRepository))` (CI runs Python 3.12, which lacks
   `typing.get_protocol_members`).
4. **`test_layering.py`**: the table above (red until the imports move).
5. **Re-point** (red until green: a test cannot import a name that does
   not exist yet, so each reads `Bounds`, `TileFootprint` or `windowed` at
   call time and fails on its own, not at collection; ruled below, pin 1):
   `tests/python/test_io_cog.py@618328b:131, 168-183, 233`
   and `tests/python/test_io_cache.py@618328b:190, 202` to `meta.windowed(w)`;
   `tests/python/test_mosaic_windowed.py@618328b:10, 31-32`'s fixture and its
   "HOW THIS FILE GOES RED" line; the imports of `Bounds` from `mosaic` and
   `TileFootprint` from `io.repository` to `io.models`
   (`tests/python/test_fetch_plan.py@618328b:55`,
   `tests/python/test_fetch_run.py@618328b:46`,
   `tests/python/test_io_cache.py@618328b:216, 240`,
   `tests/python/test_mosaic_windowed.py@618328b:26-27`,
   `tests/python/catchment_fixtures.py@618328b:33`,
   `tests/python/test_mosaic.py@618328b:81`); and the comment at
   `tests/cpp/unit/test_raster_node_at.cpp@618328b:78`, which cites
   `grid_domain.py:68-69`: pin it to `src_python/tin_engine/grid_domain.py@44fa7f5:68-69`.

On the prototype, the whole Python suite fails in exactly the `window_meta`
users (38 failed and 10 errors, in `test_io_cog.py`, `test_io_cache.py` and
`test_mosaic_windowed.py`) and passes elsewhere (5158 passed).

### Pinned by the red step (6cf4359), ruled

As `@tester` reported it: the base passes 5207; at 6cf4359, 203 failed, 148
errors and 4903 passed, each failure a missing name, a layering row or the
ruling; a scratch copy with the fix passes 5241, 0 failed. Each pin, ruled
by `@architect`:

1. **Re-pointed tests are red until green**, not "green before and after"
   as test 5 first said: they read `Bounds` and `TileFootprint` from
   `io.models` at call time (`importlib` inside `test_fetch_run.bounds`,
   `test_mosaic_windowed.footprints`, `test_fetch_plan.box_of` and
   `catchment_fixtures.MemoryRepository.footprints`, with the static imports
   under `TYPE_CHECKING`), so the catchment, catchment-batch, fetch and
   mosaic suites, `test_core_reduce` (5) and `test_dem_input_domain` (3)
   fail per test until green. **Accepted**: a test cannot import a name
   that does not exist yet, and a fallback to the old location would pin
   nothing; the collection stays whole. Test 5's wording is fixed above.
2. **One class each: `io.models.Bounds is mosaic.Bounds`,
   `io.models.TileFootprint is io.repository.TileFootprint`.** **Matches
   the design**: "moved", not copied. The old modules keep the name only as
   their own import, because they use it; neither re-exports it
   (`io/repository.py` drops it from `__all__`; `mosaic.py` has no
   `__all__`, and mypy strict does not re-export a plain import). So every
   `from tin_engine.mosaic import Bounds` in `src_python/` moves to
   `io.models`, including two the site table does not list:
   `src_python/tin_engine/dem_input.py@618328b:44` and the `mosaic` import
   block of `src_python/tin_engine/catchment.py@618328b:43-51`.
3. **`node_xy`**: Python ints give `float`, `np.int64` scalars give
   `np.float64`, `int32`, `int64` and `float64` arrays bit for bit.
   **Accepted**: the design's "the result's type is what it was".
4. **`index_of`**: Python floats give `float`; "unrounded" is checked with
   `approx(abs=1e-12)` on a unit grid, bit for bit elsewhere. **Accepted.**
5. **`node_box`**: four Python `float`s, including a 6 x 1 grid.
   **Accepted**: built from `node_xy` on Python ints, as the design says,
   so no numpy scalar enters.
6. **`valid_mask`**: a `np.bool_` array of the input's shape (2-D kept).
   **Accepted**: the design's signature.
7. **The infinity test**: 9 positive-weight and 7 zero-weight target nodes,
   listed by hand from the bilinear stencil; every other node byte-equal to
   the clean source's `resample`; `threads=1`; `check_point_blocks` checked
   at all 30 source nodes, the infinite one with `z` infinite.
   **Accepted**: the 9 and 7 are the design's counts on the 11 x 9 grid.
8. **`test_mosaic_windowed.py`**: `BOXES` holds tuples, turned into
   `Bounds` by `Bounds.of` at call time; the `window_meta` fixture keeps its
   name and returns `lambda meta, w: meta.windowed(w)`. **Accepted**: the
   smallest re-point; the fixture's docstring says what it is.
9. **Renames**: `test_window_meta_moves_the_corner_by_whole_cells` becomes
   `test_windowed_...`; `test_io_cache.py`'s now unused `cog` fixture is
   deleted. **Accepted.**
10. **The Protocol test** also checks that `TiffDemRepository` and
    `CacheRepository` have `load_window` and `check`. **Accepted**: green
    already, and it holds the design's "both implementations have" claim.

## Net production lines: about -28

Measured on the prototype (`count_loc.py 618328b <prototype>`), per file:

| File | Added | Removed | Net |
|---|---|---|---|
| `io/models.py` | 38 | 0 | +38 |
| `mosaic.py` | 11 | 36 | -25 |
| `catchment.py` | 20 | 35 | -15 |
| `io/cog.py` | 1 | 11 | -10 |
| `catchment_batch.py` | 1 | 6 | -5 |
| `target_grid.py` | 13 | 17 | -4 |
| `cli.py` | 4 | 8 | -4 |
| `dem_input.py` | 16 | 12 | +4 (the two signatures split over lines by the formatter, +5) |
| `io/repository.py` | 5 | 7 | -2 |
| `fetch/plan.py` | 2 | 4 | -2 |
| `domain.py`, `grid_domain.py`, `fetch/run.py` | | | -1 each |
| `burn.py`, `reference.py` | | | 0 each |
| **total** | 126 | 154 | **-28** |

No region is packed under `# fmt: skip`.

### As built (b3b38d2): -27

`python3 tools/count_loc.py 618328b b3b38d2`: 126 added, 153 removed,
**-27**. Every file matches the prototype's row above but three, each off by
one; the prototype was removed after measuring, so each is explained by
counting the code as written, not by a line-for-line comparison:

| File | Prototype | Built | The built file's lines |
|---|---|---|---|
| `io/models.py` | +38 | +37 (37/0) | `import math` 1; `Bounds` 17 (the moved class's 14, `of` 3); `node_xy` 2, `index_of` 2, `node_box` 3, `windowed` 4; `valid_mask` 3; `TileFootprint` 5 |
| `cli.py` | -4 (4/8) | -3 (4/7) | removed: the two import lines and `_off_node`'s five-line body (`src_python/tin_engine/cli.py@618328b:112, 118, 1696-1700`), every code line of the sites the table names; added: two import lines, two body lines |
| `target_grid.py` | -4 (13/17) | -3 (14/17) | added: the import 1; `TargetGrid.node_box` 4 (its `def`, `h = float(self.spacing)`, the corner, the `return`); `source_region` 2; `resample` 3; `check_point_blocks` 4 |

The total is one line from the design's -28, and every changed site is the
table's.

`@developer`'s pinned choices, each confirmed by `@architect` against the
expression it replaces at `618328b`:

1. **`reference.nodes_inside` builds each row band's `x` inside the loop**
   (`meta.node_xy(band, cols)` per band, `cols` built once). Same bits: `x`
   is `x_min + cols * delta_x` each time. The cost is one column vector per
   band, against the band's `meshgrid` and point-in-polygon tests of
   `band x cols` nodes. Confirmed.
2. **`catchment._joined` makes two `node_xy` calls**, `(row_max, col_min)`
   for the low corner and `(row_min, col_max)` for the high one. Each
   coordinate is the base's operand pair, and the margin is still added
   after. Confirmed.
3. **`mosaic._uncovered` makes one call, `node_xy(rows_hit[[0, -1]],
   cols_hit[[0, -1]])`.** The fancy index keeps an int64 array, so the four
   corners are `np.float64` as before, and the message prints them as
   before (the probe's `mosaic.plan_mosaic` lines are unchanged). Confirmed.
4. **`windowed` keeps `model_copy` without re-validation**, as `window_meta`
   did (`src_python/tin_engine/io/cog.py@618328b:104-113`), with its
   corner expression, now `node_xy(window.row0, window.col0)`. Confirmed:
   the design moves `window_meta`, it does not add a check.

**The probe at b3b38d2** (`lattice_probe-b3b38d2.txt`, run on a scratch copy
from `tools/scratch_copy.py`, `_core` copied in from a venv whose C++
matches; `git diff --quiet 618328b b3b38d2 -- include src bindings
CMakeLists.txt` exits 0). The base rerun at `618328b` reproduces
`lattice_probe-618328b.txt` byte for byte (`cmp` exits 0). `diff` of the two
outputs changes 88 lines: 72 `target_grid` lines whose case field is an
infinity (48 `resample`, 24 `check_point_blocks`), and 16 `AttributeError`
lines for `None` and a dict (`domain.check_extent`, `dem_input._past`,
`reference.nodes_inside` four each; `cli._off_node` and
`fetch.plan.plan_object` two each). The case-field filter counts 0 and its
inverse 72. Nothing else changed.

Both runs print numpy's `RuntimeWarning: invalid value encountered in
multiply` and `in add` from `resample`'s bilinear line
(`src_python/tin_engine/target_grid.py@b3b38d2:191`,
`src_python/tin_engine/target_grid.py@618328b:193`): the formula is
evaluated before the NoData mask, at both revisions, so a zero weight on an
infinity (0 x inf) or +inf meeting -inf warns. The warning is not new; what
the ruling changes is that the NaN it reports is now kept as the node's
value instead of being replaced by the fill. It is in Question 1.

## `@perf`

Owed: A touches what drives refine and mesh (`grid_domain.subsample` gives
the start vertices; `mosaic` and `target_grid` give the tile refine reads).
The gate is **byte-identical meshes**: `tools/bench.py` on the 1 m set,
A's head against `618328b` (`--tree`), back to back, every mesh file equal.
No timing claim is made, so none is owed beyond the runs. The bench's seam
(`tools/bench.py@618328b:750-761` swaps `refine` in `cli`'s namespace) does
not move. The evidence goes in `docs/benchmarks/<date of the run>/audit-pr-a/`
(the two run directories and a short `audit-pr-a.md` beside them with the
mesh hashes), as `docs/benchmarks/2026-10-04/` does for 20 and 23b.

The bench's domains (the whole tile and `--domain` files, in the DEM's own
CRS) reach at most `grid_domain`, `mosaic`, `domain` and `dem_input`, never
`target_grid` (the only change in behaviour), `burn`, `catchment` or
`fetch`: `tools/bench.py run` passes extra mesh arguments after `--`
(`tools/bench.py@618328b:636, 680`), and the acceptance run above passes
none, so no `--out-crs`.

**One reprojected run is added**, because `target_grid` holds A's only
change in behaviour and four of its sites (`resample`'s fractional index,
`check_point_blocks`' node coordinates, the grid's rectangle twice) are
rewritten: A's head against `618328b`, back to back,
`tools/bench.py run --label audit-pr-a-reprojected --tree <tree> --dem
../rasputin_data/sao_francisco_piece/bho2017_5k_76949_anadem_window_epsg4674.tif
--domain ../rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson
--tolerance 10 -- --out-crs EPSG:31983` (the Velhas piece of the São
Francisco basin, ANADEM at about 30 m in EPSG:4674, meshed in UTM 23S: the
domain and frame of `docs/benchmarks/2026-10-04/15f-3-acceptance/velhas.py`,
which reads the `anadem-v1` cache; here the window file beside the outline
in `../rasputin_data`, whose extent holds the outline's bounds).
Expected: **every mesh file byte-identical**. The window holds no infinity
and no NaN (its sidecar, `bho2017_5k_76949_anadem_window_epsg4674.json`,
says `nodata_nodes: 0`; `np.isinf` over the array counts 0), so the
ruling changes nothing there and any difference is a defect in a rewritten
site. `@perf` checks the command runs at `618328b` before relying on it;
if `bench.py` refuses the geographic file, it says so and the probe stays
the only gate for `target_grid`.

`burn`, `catchment`, `fetch` and the rest of `target_grid` are covered by
the differential probe alone: `open_dem` with a domain reprojected to UTM
32, `resample`, `check_point_blocks`, `burn_reach`, `delineate` and
`plan_object`, value for value. The probe is the gate for them; the bench
is the gate for the meshes refine builds.

**Accepted (01751ec):** `docs/benchmarks/2026-10-06/audit-pr-a/audit-pr-a.md`.
Meshes byte-identical on both gates: the 1 m set, tile and quarter (each 9
of 9 whole `.vtk` files equal: four base runs, four head runs and the head's
build run), and the reprojected Velhas run in EPSG:31983 (2 of 2). Quality
unchanged, as the hashes imply. Refine time on the 1 m set, four alternating
runs a side pooled, is within -2.5 to +2.6 % on every cell, against a
base-against-base spread of median 2.1 % and at most 5.3 %. Velhas ran once
a side (refine -2.0 to +1.3 %), so it has no noise floor of its own.

## Merging with B, C, D and F

The rule, as PR C's: whichever of A and another audit branch merges second
runs `git merge-tree --write-tree <its head> master` and resolves every file
it names. The heads move, so what follows is the prototype merged against the
heads it was measured at (F `e2baa5f`, C and D `b26beb8`), not a statement
about today's: D has since added rows to `test_layering.py`
(`git log b26beb8..worktree-audit-encoders -- tests/python/test_layering.py`),
and the merge-tree run at merge time is what counts.

- **F** (`worktree-audit-catchment`, `e2baa5f`): `catchment.py`'s import
  block. Merged: F's block, with B's `crs` line (`same_crs, single_crs`), the
  `io.models` line as `Bounds, DemTile, RasterMeta, TileFootprint`, and
  `io.repository` reduced to `DemRepository` (A drops `TileFootprint` from
  its `__all__`, so mypy strict refuses F's import of it from there). The rest
  of `catchment.py`, `catchment_batch.py`, `cli.py` and `test_layering.py`
  merge clean.
- **C** (`worktree-audit-geojson`, `b26beb8`) and **D**
  (`worktree-audit-encoders`, `b26beb8`): `test_layering.py`, two hunks.
  Merged: C's `fetch.nve` row with A's `fetch.plan` and `fetch.run` rows; in
  `UPWARD`, neither A's two lines nor C's `("chains", "feature_input")`.
  `domain.py`, `io/repository.py` and `cli.py` merge clean.
- **B**: A is built on it.
- **New since the design** (`git merge-tree --write-tree 1a84ea9 <head>`,
  2026-10-06): `project_structure.md` conflicts with C (`b26beb8`) and with
  D (`4725f12`), because A rewrote six of its rows after green; and
  `cli.py`'s `io.models` and `io.ply` import lines conflict with D, which
  adds `io.mesh_checks` there. Each resolves by keeping what both sides add:
  in `project_structure.md` the `dem_input`, `domain`, `crs` and
  `target_grid` rows' text, in `cli.py` the imported names.

## Citations this PR moves, pinned now

Found by running `python3 tools/check_citations.py --base 618328b` in the
prototype (37 at risk), then comparing each cited line at `618328b` with the
same line number in the prototype. 23 move. Pinned in this design's second
commit, each quotation re-read at its pin:

- to `44fa7f5` (master; each cited line reads the same there and at
  `618328b`): `docs/increments/15-dem-mosaic.md` line 616,
  `docs/increments/15e-memory-fixes.md` lines 53, 85, 93 and 95,
  `docs/increments/25-plain-output.md` lines 93 and 122, and
  `docs/increments/27-node-sampling.md` lines 114, 173, 174, 402 and 404;
- three review records whose citation had already moved before A, to the
  revision they reviewed: `docs/increments/29-nve-reference-catchments.md`
  lines 3088 and 3090 (`mosaic.py:385`, the `_mixed` refusal in
  `_covering`, to `e7b5923`, the end of round 4's range; round 5 changed no
  code) and line 3114 (`cli.py:1925-1928`, which the record itself reads "at
  `0588bfa`", to `0588bfa`).

None of these lines is edited by C, D, F or T1 (`git diff -U0 618328b...<branch>`
on the five files). Not pinned here: `tests/cpp/unit/test_raster_node_at.cpp`
line 78 (a test file: red test 5); `docs/increments/python-audit.md` line
1137 (`fetch/plan.py:118`, B's own text, true at `618328b`, not on master;
at A's head the cited `same_crs` line is
`src_python/tin_engine/fetch/plan.py@1a84ea9:117`: on the rebase
checklist, pinned when A is rebased onto master after B merges); and `python-audit.md`
lines 1271-1277, which list other files' citations in B's section 9 and sit
in the file every audit branch edits.

**Moved by the red step (6cf4359), pinned after green:**
`docs/increments/15e-memory-fixes.md` line 327 cited lines 399-400 and
514-515 of `test_target_grid.py` (the two "Went red at 9879805" comments);
the infinity tests inserted above them move them to 476 and 591. Both
citations are now pinned, `tests/python/test_target_grid.py@44fa7f5:399-400`
and `tests/python/test_target_grid.py@44fa7f5:514-515`, where they read as
cited.

`project_structure.md`'s rows for `io/models.py`, `io/cog.py`,
`io/repository.py`, `mosaic.py`, `target_grid.py` and `dem_input.py` are
rewritten against the code at b3b38d2.

## Constants

None added or changed. `mosaic.ALIGN_TOLERANCE` and `SEAM_THRESHOLD` stay
where they are, unchanged.

## Questions for Ola (defaults hold until he answers)

1. **A grid node next to an infinite DEM value, after reprojection.** Under
   your ruling, an infinite value is data. When the reprojected grid's node
   sits exactly on a source node whose neighbour is infinite, the
   interpolation multiplies that infinity by zero and gets "no value", so
   the node has no elevation, though the source node under it is finite.
   The direct (not reprojected) path reads that source node alone and gets
   its finite value. Such a DEM also makes numpy print a warning on the
   terminal during reprojection ("invalid value encountered in multiply");
   that warning is printed before PR A too, from the same line. Accept both
   for PR A, and treat them only if a real DEM ever holds an infinity?
   Default: accept, and leave the warning as it is until Question 2 is
   answered.
2. **Should a DEM tile that holds an infinite elevation be refused when it
   is read?** An infinite height is not terrain; today it is data (your
   ruling), and a mesh can carry it. Default: no change now; ask again if
   such a file turns up.

## Review

**PR A (`audit-lattice`), design review, round 1, 2026-10-06.** Range `618328b..62b5226`. Verdict: CHANGES REQUESTED. 0 production lines. The probe reproduces its base output byte for byte and fails under planted mutants (72 lines for the ruling's `_valid`, 3 for a `grid_domain` shift); all 15 re-pinned citations quote what they claim; the site table matches `git grep` at 618328b; the merge with F conflicts only in `catchment.py`'s import block. Blocking: line 268 must read 'no mutation round'; the site table lacks `src_python/tin_engine/mosaic.py@618328b:265`'s `window_meta` docstring.

**PR A (`audit-lattice`), design review, round 2, 2026-10-06.** Range `62b5226..0b20e1d`. Verdict: APPROVED. 0 production lines. Both round-1 blockers fixed ('Lean: no mutation round.'; the site-table row for `src_python/tin_engine/mosaic.py@618328b:265`). On `docs/increments/python-audit-probes/lattice_probe-618328b.txt@0b20e1d` the inf filter on the case field counts 81, inverted over `resample`/`check_point_blocks` 72. The `node_box` and F10 departures, the `@perf` claim about which paths the bench reaches, the probe docstring and the two corrected review records check out. `check_citations --base 618328b` exits 0.

**PR A (`audit-lattice`), code review, round 1, 2026-10-06.** Range `0b20e1d..1a84ea9`. Verdict: APPROVED. -27 net production lines (`count_loc.py 618328b 1a84ea9`: 126 added, 153 removed). Every rewrite matches the site table and the ten ruled pins. The probe rerun at 618328b and b3b38d2 reproduces both committed outputs byte for byte; their diff is 88 lines (72 infinity cases, 16 AttributeError), the case filter counts 0, and a planted `node_xy` mutant changes the 5 lines the design names. Red real at 6cf4359; suite green at the head (5254 passed, 17 skipped). mypy, ruff, gates and `check_citations --base 618328b` pass; 56 at-risk citations re-read.
