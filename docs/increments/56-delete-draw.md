# Increment 56: delete `rasputin draw` and the drawing package

**Status: designed** (`@architect`, 2026-10-16), one public PR, branch
`worktree-delete-draw`.

## 1. What Ola asked

Ola, 2026-10-16: "The best lines of codes are the ones you never built because
you understood they were not needed. The second best lines of codes are the
ones you deleted for the same reason :)". On the main session's proposal (delete
`rasputin draw`, the SVG renderer and their tests; move only the test shapes the
remaining tests need under `tests/`): "They all sound good to me."

Increment 54 (`docs/increments/54-open-core-cleanup.md`, section 7, question 1)
hid `draw` and left its deletion to this yes. Increment 54 also removed the
test-shape gallery from `rasputin mesh`, so `draw` is the gallery's last product
user.

## 2. Prior art: legacy and literature

Written before the design.

*Literature.* None applies: the increment adds no method and claims nothing
new. It removes a developer's picture tool.

*Legacy.* The legacy tree had no drawing tool to carry or to compare with. The
grep, at the tag:

```
$ git grep -l -i -E 'svg|matplotlib|def draw|plot' legacy-archive -- legacy
legacy-archive:legacy/tests/test_gml_repository.py
```

That one file imports `matplotlib` behind a `plot = False` switch in a test;
nothing is carried across.

## 3. Who uses the drawing package today

Run at `c5d4be39`, the branch's base:
`git grep -n -E 'tin_engine\.viz|viz\.|viz/|rasputin draw|GALLERY' -- src_python tests tools bindings include pyproject.toml`.

| User | What it uses | Fate |
|---|---|---|
| `src_python/tin_engine/viz/` (6 files, 1,211 lines) | itself | deleted |
| `cli.py`: `draw`, `Attempt`, `_triangulated`, `FIXTURE_CLOSED_ROLES`, `_PRECEDENCE`, `_BIT_OF`, `PROPERTY_STROKES`, `DEFAULT_LABEL_LIMIT`, the five `viz` imports, `app.command(hidden=True)(draw)` | the only product caller | deleted |
| `cli.py`: `ROLES`, `CORE_CLOSED_ROLES`, `_engine`, `_band`, `_snap_spacing`, `DEFAULT_SNAP_SPACING` | also used by `mesh` | kept; comments that cite the gallery or `viz` reworded |
| `pyproject.toml` line 195, `"src_python/tin_engine/viz"` | a path list | the entry deleted |
| `cli_kit.destination`: the `out is None` branch (a temporary `<name>.svg`) | reached only from `draw` | deleted, `out: Path` (section 7, question 1) |
| `features.py` module docstring and the `DEFAULT_VOCABULARY` comment | prose about `viz` | reworded; the name pattern stays (below) |
| `bindings/core.cpp:141`, `include/terrain/core/noded_pslg.hpp:44` | comments citing `viz/protocols.py`'s `PslgLike` | reworded, comments only |
| `tests/python/test_cli_draw.py`, `test_viz_protocols.py`, `test_viz_scene.py`, `test_viz_svg.py` (3,858 lines in all) | the drawing | deleted, except the obligations moved in section 5 |
| `tests/python/test_cli_mesh.py` | `GALLERY["road-crosses-river"]` through `cli._triangulated` | the one shape that moves (section 4) |
| `tests/python/test_cli_help.py` T4 | `draw --list` runs while hidden | replaced by the absence test |
| `tests/python/test_layering.py` | the `viz` rows and `viz.*` in the `cli` row | rows deleted |
| comments in `test_features.py`, `test_core_noding.py`, `test_core_cdt.py`, `test_io_ply.py`, `test_mcp_schema.py`, `tests/cpp/unit/test_noding_node.cpp`, `tests/cpp/unit/test_cdt_detria_backend.cpp` | prose naming `viz` or `rasputin draw road-crosses-river` | reworded to "the road-crosses-river shape" or cut |
| `testing.md` (`### viz`, the stylesheet exclusion), `project_structure.md` (the `viz/` tree, the binding and API paragraphs) | prose | the `viz` section and lines deleted; the binding paragraph says `cli.py` is the composition root |

**The "seven test files that load the shapes" are not gallery users.**
`test_cli_mesh_features.py` defines its own `GALLERY`, a list of GeoJSON
features, and `test_cli_class_maps.py`, `test_cli_features_cleanup.py`,
`test_cli_start_quality_gain.py`, `test_cli_threads.py` and
`test_cli_mesh_multi_features.py` import that one. None of them touches
`tin_engine.viz`; they do not change. No tool under `tools/` and no C++ test
reads the gallery: the C++ crossing (`test_noding_node.cpp:157-165`) is its own
copy of the numbers.

`EdgeProperty`'s name pattern `^[a-z][a-z0-9_]*$` was written for CSS class
tokens. It stays: the PLY header writes `feature_bit {bit} {name}`
(`src_python/tin_engine/io/ply.py:123`), space-separated, and the GeoPackage writer puts each name, quoted but not
escaped, into a `CREATE TABLE` statement as a column
(`src_python/tin_engine/io/geopackage_write.py:321`), so the pattern is now
what keeps a quote out of that SQL. The `features.py` docstring gives these
reasons instead.

## 4. The one shape that moves

Of the eleven gallery shapes only `road-crosses-river` is needed, by
`test_cli_mesh.py::TestTheChainMaskJoin`: a noded graph whose chains carry
feature bits, so that `cli._chain_masks` (kept, used by `mesh`) is checked
against an independent rebuild. It moves **inline into that test**, not into a
new helper module: eight vertices and three chains, built through
`cli._engine` (kept) instead of `cli._triangulated` (deleted). Its twin already
lives in `tests/python/test_core_noding.py::crossing_pslg`, at the same origin;
that one is a pytest fixture returning a `Pslg`, and the join needs a
`NodedPslg`, so the test writes the numbers again rather than chaining
fixtures across modules.

The shape, at origin `(430_000, 6_900_000)` (UTM 33N magnitudes, the old
`viz.fixtures.ORIGIN` and `test_core_noding.py`'s `EAST`, `NORTH`):

- outer ring `(0,0) (700,0) (700,700) (0,700)`, role `outer`, mask 0;
- road `(100,100)-(600,600)`, role `breakline`, mask 2;
- river `(100,600)-(600,100)`, role `breakline`, mask 1.

Measured at `DEFAULT_SNAP_SPACING` (1e-3 m), `_engine(..., delaunay=True)`:
`ok` True, status `Ok`, 9 noded vertices (the crossing at `(350,350)` added),
12 triangles, `_chain_masks` gives 8 pairs with values `{0: 4, 1: 2, 2: 2}`.
Command: the eight lines above through `cli._engine` and `cli._chain_masks`
with `/Users/skavhaug/projects/rasputin/.venv/bin/python` (master `9a7c312d`'s
build; `_engine` and `_chain_masks` are unchanged since).

## 5. Tests: what is deleted, what moves, what is new

All test work is `@tester`'s, in the red commit.

**Deleted whole:** `test_cli_draw.py`, `test_viz_protocols.py`,
`test_viz_scene.py`, `test_viz_svg.py`. What they pin is the picture, the
gallery shapes and `draw`'s options. The engine refusals the gallery showed
(validator rejection, backend `InvalidTopology`, the noder's
`CoordinateOutOfRange` and `NotConverged`) are pinned where they live:
`test_core_cdt.py`, `test_core_noding.py`, `test_cdt_detria_backend.cpp`,
`prop_noding_verifier_sweep.cpp`.

**Moved, because the code they guard stays** (from `test_cli_draw.py`):

- T1, `cli_kit.destination`'s boundary, into a new
  `tests/python/test_cli_kit.py` calling `destination` directly (no command,
  no compiled core): a path inside `out_parent` is returned resolved; `..` out
  of `out_parent` is refused with "outside"; an absolute path elsewhere is
  refused with "outside"; a symlink pointing back inside is refused with
  "symlink" and the target file is untouched; a missing directory is a
  `typer.BadParameter` naming it; the refusal's `param_hint` is the `hint`
  passed. Today only the "outside" and the symlink case are reached through
  `mesh` (`test_cli_mesh_plain_output.py::test_outside_out_parent`,
  `test_cli_mesh_stats.py::test_stats_through_a_symlink`).
- T2, `--snap-spacing` on `mesh`, in `test_cli_mesh.py` with the `bumpy`
  fixture (`rough_dem(16)`, 17 x 21 nodes near `(500_000, 6_600_000)`):
  `0`, `-1`, `nan`, `-inf` exit 2, the output names `snap-spacing`, and no
  `--out` file exists; `--snap-spacing=1e-12` exits non-zero and the output
  names `CoordinateOutOfRange` (the value reaches the engine). Measured on
  master `9a7c312d`: exit 2 and "snap-spacing must be finite and positive" for
  the first four; exit 2 and "has no mesh to write: CoordinateOutOfRange" for
  1e-12; 0.5 exits 0 with 640 triangles.
- T3, `ROLES` maps every `ChainRole` member and nothing else, into
  `test_cli_mesh.py` (measured True).
- T4, `TestTheChainMaskJoin` on the inline shape of section 4; the `attempt`
  fixture, `FIXTURE` and `FEATURED` go, and the module docstring stops
  pointing at `test_cli_draw.py`.

**Changed:**

- T5, absence, replacing `test_cli_help.py`'s
  `test_draw_is_hidden_and_still_runs` and `test_no_tool_is_offered_for_draw`:
  `rasputin draw` exits 2 with "No such command";
  `importlib.util.find_spec("tin_engine.viz")` is None; `cli` has none of
  `Attempt`, `_triangulated`, `FIXTURE_CLOSED_ROLES`, `PROPERTY_STROKES`,
  `DEFAULT_LABEL_LIMIT`. This is the test that is red before the green commit.
- T6, `test_layering.py`: the five `viz` rows and the `viz.*` names in the
  `cli` row deleted.
- `test_mcp_schema.py`: the comment above `TOOL_NAMES` stops saying `draw` is
  hidden; after the rebase (section 8) no `draw` is left in `TOOL_NAMES`.

Not added: a test for the `destination` `None` branch's absence; `mypy`
strict refuses `None` once the parameter is `Path`.

## 6. The green step

`@developer`: delete `src_python/tin_engine/viz/`; in `cli.py` delete the rows
of section 3, reword the module docstring and the comments on `ROLES`,
`CORE_CLOSED_ROLES`, `DEFAULT_SNAP_SPACING`, `_engine` and `_band` so none
names the gallery, `viz` or `svg.py`; delete the `pyproject.toml` entry; drop
`destination`'s `None` branch and type `out: Path`; reword `features.py` and
the two C++ comments. No test file. The prose files of section 3 change in the
same PR (they are not production code; any persona with the write right).

No refine or mesh code changes behaviour, so no `@perf` run: the C++ edits are
comments, and `mesh`'s path through `cli.py` loses no line it runs.

## 7. Questions for Ola

1. `plugin_api` promises `destination` and says a promised name changes "only
   in a release that says so". With `draw` gone, its `out=None` branch (write
   to a fresh temporary `<name>.svg`) has no caller; the one plugin passes a
   path. **Default: delete the branch and say so in the PR description.**

## 8. Size, overlap and order

**Counted size: about -717 net lines** (`tools/count_loc.py`): `viz/` -614,
`cli.py` -103; the `destination` branch -3 more. Measured on a scratch clone at
`c5d4be39` with the section 3 rows deleted (it parses, and `ruff --select F`
finds no unused name). Tests: about -3,860 lines deleted, about +120 moved or
new, not counted.

**Overlap with increment 55** (branch `worktree-mcp-datasets`, based on
`9a7c312d`, merges first): it edits `cli.py` (the `sources` import block next
to the `viz` imports, and the command registration at the end of the file),
`test_layering.py` (the same `cli` row this PR shortens), `test_mcp_schema.py`
(`TOOL_NAMES`, which lists `draw` on that branch because it predates increment
54), `project_structure.md`
and `ROADMAP.md` (its row 55 sits where row 56 goes). All are small textual
conflicts. Order: red and green may start now; the branch rebases onto master
once increment 55 merges, before `@reviewer`'s code review, and the rebase
takes increment 55's side and then deletes the `viz` names again.

Short path: no. The change touches two C++ files (comments) and is far from
small in lines moved, though negative in count.
