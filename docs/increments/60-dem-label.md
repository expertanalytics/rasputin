# Increment 60: delete the DEM's run name (`DemInput.label`)

**Status: designed** (`@architect`, 2026-10-16, on `ed0d313d`, increment 57's
head; merges after 57). Short path: Python only, under 50 net lines. One public
PR, branch `worktree-dem-label`.

## 1. What Ola asked

Ola, 2026-10-16: "The best lines of codes are the ones you never built because
you understood they were not needed. The second best lines of codes are the
ones you deleted for the same reason :)" and "We don't want to pile up legacy
so hard no to aliases."

Increment 57 removed `--out-parent`, the last reader of the run name `open_dem`
returns (`DemInput.label`: the directory's name, the first file's stem, or the
cached source's key). `@reviewer`, in 57's code review round 1, left its
deletion to a follow-up: `open_dem` and `repository_for` are both in
`tin_engine.plugin_api`, so two signatures a plugin sees change.

## 2. Prior art: legacy and literature

*Literature.* None applies: nothing is added and nothing new is claimed.

*Legacy.* Nothing to carry: the legacy tree has none of these names.

```
$ git grep -l -E 'DemInput|repository_for|open_dem' legacy-archive -- legacy
(no output, exit 1)
```

## 3. Who produces and who reads it

Run at `ed0d313d`: `grep -rn 'repository_for' src_python tests tools` and
`grep -rnE '\.label\b|label=' src_python tests`, less the unrelated
`label`s (class maps, `toml_data`, `bench.py`).

Producers, all in `src_python/tin_engine/dem_input.py`: the field
`DemInput.label`; `repository_for`'s second return value; `open_dem`, which
unpacks it and passes it on; `_open_reprojected`'s `label` parameter.

Production readers of the label: none. Callers of `repository_for` that
discard it with `[0]`: `cli.py:1671` (`_dem_repository`) and
`source_report.py:321`.

Test readers: `tests/python/test_dem_input.py:127,136,144` and
`tests/python/test_cli_mesh_cache.py:254`, each `assert ... .label ...`.
One more test line builds the pair: `tests/python/test_dem_input_domain.py:671`,
whose stand-in `repository_for` returns `(InMemory(tiles), "memory")`.

Outside the package: the private plugin `rasputin-lsm` calls
`repository_for(paths)[0]` once (`src/rasputin_lsm/cli.py:117` on
`private/master`). It pins rasputin at a commit, so nothing breaks until the
pin moves; the pin bump drops the `[0]`. Dated scripts under `docs/benchmarks/`
and `docs/increments/perf-audit-probes/` unpack the pair; they record what ran
at their commits and are left as they are.

## 4. The change

- `DemInput` loses `label`; its docstring loses "a name for the run (...)".
- `repository_for(...) -> TiffDemRepository | CacheRepository`: the
  repository alone. Its docstring drops "and the run's name ...".
- `open_dem` keeps its signature; its result has no `label`.
- `_open_reprojected` loses its `label` parameter.
- The two `[0]` callers drop the `[0]`.
- No alias, no deprecation shim (Ola: "hard no to aliases").
- `docs/plugins.md`, the `repository_for` row becomes:
  `| repository_for | The tile repository over a directory, files or a cached source. It no longer returns the run's name beside it, and the DemInput that open_dem returns no longer has label |`
  (with the code spans as the table writes them). This is the "release that
  says so" the page's import promise asks for, as the `destination` row does
  for 57.

## 5. Tests

Red (`@tester`), in `tests/python/test_dem_input.py`, no fixture beyond the
file's `micro_tiff` (one tile):

- `DemInput` has no field named `label` (`dataclasses.fields`);
- `repository_for((path,))` on that tile returns an object with
  `footprints()`, not a tuple.

Both fail at `ed0d313d`. The same commit deletes the four `.label` assertions
and changes the stand-in at `test_dem_input_domain.py:671` to return
`InMemory(tiles)`; those changes alone would also go red at `ed0d313d`, which
is fine, since green touches no test file.

## 6. Size

About -3 counted lines: `label: str` in `DemInput`, the `label: str,`
parameter and the `label=label,` argument go; every other touched line is
changed in place (return type, three returns, two call sites, the unpacking).
Checked by `python3 tools/count_loc.py <base> <head>` at review.

## 7. Left out

The dated scripts in `docs/` (above). The private plugin's `[0]`, done at its
next pin bump.
