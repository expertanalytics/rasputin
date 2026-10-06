# PR A design: `audit-lattice` (F2, F9, F10's repository Protocol, F12's two `mosaic` edges)

Status: **design review round 1 answered (by `@architect`); next, `@reviewer`'s
design review round 2; no `@tester` step yet.** Branch
`worktree-audit-lattice`, stacked on PR B's approved head `618328b`
(`worktree-audit-crs`, not pushed); it is rebased onto master when B merges.
Every `@618328b` citation below reads B's head, which B's merge commit keeps
in master's history. The audit this designs from is
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
   any diff: `diff <base> <head> | awk -F' [|] ' '/^> target_grid/ && $2 !~
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
geographic grids, `cli._off_node` on both).

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
5. **Re-point** (green before and after, except the first three, which fail
   once `window_meta` goes): `tests/python/test_io_cog.py@618328b:131, 168-183, 233`
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
`fetch`, and the bench has no `--out-crs` to pass. Those four are covered
by the differential probe alone: `open_dem` with a domain reprojected to
UTM 32, `resample`, `check_point_blocks`, `burn_reach`, `delineate` and
`plan_object`, value for value. The probe is the gate for them; the bench
is the gate for the meshes refine builds from the DEM's own grid.

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
1137 (`fetch/plan.py:118`, B's own text, true at `618328b`, not on master:
pinned when A is rebased onto master after B merges); and `python-audit.md`
lines 1271-1277, which list other files' citations in B's section 9 and sit
in the file every audit branch edits.

`project_structure.md`'s `io/models.py` row ("Pydantic RasterMeta /
DemTile") is rewritten by `@architect` after the green commit, so it
describes the code as written (`Bounds`, `IndexWindow`, `TileFootprint`, the
node methods, `valid_mask`).

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
   its finite value. Accept this for PR A, and treat it only if a real DEM
   ever holds an infinity? Default: accept.
2. **Should a DEM tile that holds an infinite elevation be refused when it
   is read?** An infinite height is not terrain; today it is data (your
   ruling), and a mesh can carry it. Default: no change now; ask again if
   such a file turns up.

## Review

**PR A (`audit-lattice`), design review, round 1, 2026-10-06.** Range `618328b..62b5226`. Verdict: CHANGES REQUESTED. 0 production lines. The probe reproduces its base output byte for byte and fails under planted mutants (72 lines for the ruling's `_valid`, 3 for a `grid_domain` shift); all 15 re-pinned citations quote what they claim; the site table matches `git grep` at 618328b; the merge with F conflicts only in `catchment.py`'s import block. Blocking: line 268 must read 'no mutation round'; the site table lacks `src_python/tin_engine/mosaic.py@618328b:265`'s `window_meta` docstring.
