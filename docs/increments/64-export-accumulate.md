# Increment 64: `plugin_api` exports the flow accumulation

**Status: designed** (`@architect`, 2026-10-17). Short path; one public PR,
branch `worktree-export-accumulate`. 61 to 63 are taken (`git ls-tree` over
every ref of both remotes after `git fetch --all`; `ls` in every worktree).

## 1. Why

The private package's routing step (private increment 61, branch
`lsm-ddd-routing`, its §4.1; Ola's yes on 2026-10-16 to DDD-style routing)
needs each DEM node's receiver and drained count. The core computes both, but
a plugin may import only `tin_engine.plugin_api`, which lacks them.

## 2. Prior art: legacy and literature

- *Literature.* Unchanged method: Priority-Flood (Barnes, Lehman and Mulla 2014a,
  doi:10.1016/j.cageo.2013.04.024), as increments 22 and 29 cite it. Each node
  drains to the neighbour that flooded it (its *lowest* filled neighbour), not
  to the steepest as in O'Callaghan and Mark's D8 (1984); increment 22 names
  the departure, so private increment 61's "the core's D8" is loose. No novelty.
- *Legacy.* Nothing: `git grep -n -i -E "accumulat|flow_dir|d8" legacy-archive -- legacy`
  returned no file (exit 1, 2026-10-17).

## 3. What changes

- `src_python/tin_engine/plugin_api.py`: `from tin_engine.catchment_core import
  accumulate`, and `"accumulate"` in `__all__`. It re-exports the Python wrapper
  `catchment_core.accumulate(tile: DemTile) -> AccumulateOutcome`, the way it
  already re-exports `raster.core_args`, not the binding `_core.accumulate`,
  which takes a `RasterView`. The wrapper runs on a unit lattice: what it
  returns is in nodes and neighbour codes, not metres. It releases the GIL.
- `docs/plugins.md`'s table gains an `accumulate` row (this commit) promising
  `.count` and `.flow_to`, dtypes and codes; not `.reach` or the outcome type.
  The wrapper's docstring says the same, since that page makes a docstring
  "the full contract".

## 4. Tests (`@tester`, red)

- `tests/python/test_plugin_api.py`: `"catchment_core": ("accumulate",)` in
  `PROMISED`. Red now: `__all__` lacks the name, and the identity test finds
  no `plugin_api.accumulate`. No fixture (names and identity
  only); `tests/python/test_core_accumulate.py` already pins behaviour.
- `tests/python/test_layering.py`: `catchment_core` in `plugin_api`'s row.
  This permits the import; it is not a red pin.
- T6 (`tests/python/test_cli_help.py`, `TestThePluginsPage`) stays green: the
  table row above lands with this design.

## 5. Left out, and the estimate

Left out: a private binding. About 2 counted lines (`CLAUDE.md` §2); no C++, nothing under *Acceptance*.

## Review

- **Design review, round 1, 2026-10-17.** Commit `f5b46020`, 0 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session: the `accumulate` row in docs/plugins.md says which neighbour a node drains to and that the call releases the GIL. Held: the wrapper, not the binding, is the boundary; dtypes, codes, read-only views and the ValueError match the code.
