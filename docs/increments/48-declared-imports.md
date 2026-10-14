# Increment 48 — Every third-party import is declared, and a bad `--z-unit` is a usage error

Status: **designed** (`@architect`, 2026-10-14), short path
(`docs/increments/README.md`, "Small changes: the short path"). **This
increment absorbs increment 47** (`--z-unit` usage error, designed on branch
`worktree-z-unit-usage`, commits `84b40250` and `069dd4f0`); the number 47 is
retired and not reused. Both defects have one cause and one fix (section 2),
so one PR carries the fix and both sets of tests; the main session took this
default while Ola was away. Every `file:line` is pinned to master `55b9f564`.

## 1. The asks

- `@perf` reported that a fresh install fails because
  `src_python/tin_engine/cli.py@55b9f564:58` does `import click` and
  `pyproject.toml` does not declare click.
- `rasputin mesh --z-unit furlong` prints a traceback and exits 1, where every
  other bad option value is a one-line usage error with exit 2.

Ola approved each as its own small PR on 2026-10-14 (as relayed by the main
session; no verbatim quote).

## 2. What the check found

**The missing import is real.** Reproduced in the main checkout's venv
(`.venv`, editable install, typer 0.27.2, no `mcp` extra):

```
$ /Users/skavhaug/projects/rasputin/.venv/bin/rasputin --help
  File ".../src_python/tin_engine/cli.py", line 58, in <module>
    import click
ModuleNotFoundError: No module named 'click'
```

typer 0.26.0 dropped click as a dependency and carries its own copy as
`typer._click` (typer 0.26.0 release notes, "Breaking Changes", PR #1774:
"Typer no longer depends on Click as a third party dependency, it vendors
(includes the source code of) Click. [...] Click-specific functionality is no
longer supported"). PyPI metadata: typer 0.25.0 requires `click>=8.2.1`;
0.26.0 and 0.27.0 do not. `pyproject.toml` says `typer>=0.12`, so pip picks
0.27.x. CI is green only because it installs `.[dev]`: `dev` holds
`mcp>=2.3`, which requires `uvicorn`, which requires `click>=7.0` (checked in
the `catalogue-data` worktree's venv: uvicorn is the only installed
distribution requiring click).

**The traceback has the same cause.** The one use of click is
`click_type=click.Choice(...)` on `Z_UNIT`
(`src_python/tin_engine/cli.py@55b9f564:229`), shared by `mesh`, `catchment`
and `station-catchments`. typer wraps the foreign `Choice` as a plain function
type, so the real-click `BadParameter` it raises is not the class typer's
handler catches. Probed by increment 47's design on master's code (typer
0.27.3, click 8.5.0):

| Command line | Exit | Last line |
|---|---|---|
| `mesh --z-unit furlong --out x.vtk` | 1 | `BadParameter: 'furlong' is not one of 'm', 'ft', 'us-ft'.` |
| `catchment --z-unit furlong` | 1 | same |
| `station-catchments --z-unit furlong` | 1 | same |
| `mesh --tolerance abc --out x.vtk` | 2 | `Invalid value for '--tolerance': 'abc' is not a valid float.` |

The same wrapping shows the option as `<function>` in `--help`, and gives the
MCP tool schema (`tin_engine/mcp/schema.py`, `_annotation`) any string for
`z_unit` instead of the three spellings, because the wrapped type has no
`choices`.

**So declaring click is the wrong fix**: the import would succeed and all
three faults would stay. The fix removes the import and declares the choice
the way typer supports.

**Every other third-party import.** An AST scan of all 70 `.py` files under
`src_python/tin_engine/` (every `import` and absolute `from` import,
including ones inside functions), less the standard library
(`sys.stdlib_module_names`) and `tin_engine` itself, finds ten top-level
names:

| import name | files | declared? |
|---|---|---|
| numpy | 38 | yes |
| pydantic | 27 | yes |
| shapely | 24 | yes |
| pyproj | 6 | yes |
| tifffile | 5 | yes |
| typer | 3 | yes |
| mcp | 1 (`mcp/server.py`) | yes, `mcp` extra |
| **click** | 1 (`cli.py`) | **no** |
| **mcp_types** | 2 (`mcp/schema.py`, `mcp/server.py`) | **no**; distribution `mcp-types`, pulled by `mcp` 2.3.0 as `mcp-types==2.3.0` |
| **anyio** | 1 (`mcp/server.py`) | **no**; pulled by `mcp` as `anyio>=4.9` (`>=4.10` on Python 3.14) |

`anyio` and `mcp-types` cannot fail today (the `mcp` extra always brings
them), but they are imported by name, so they are declared by name (deptry's
rule DEP003, section 3). Each import name maps to its distribution by PEP 503
normalisation (lower case, `_` to `-`); no alias table is needed.

## 3. Prior art: legacy and literature

*Literature.* This is packaging practice and the CLI library's own choice
type used as documented; nothing new is claimed.

- PEP 621, *Storing project metadata in pyproject.toml*: `dependencies` and
  `optional-dependencies` are the project's declared requirements.
- PEP 503, *Simple Repository API*, "Normalized Names": the normalisation the
  import gate uses to match an import name to a requirement.
- deptry (`https://deptry.com/rules-violations/`), rules DEP001 (an imported
  module has no declared package) and DEP003 (an imported module comes only
  through another package's dependency tree). Click here was a DEP003 case
  that became DEP001 the day typer 0.26 shipped. The gate below is those two
  rules for one package, as a short pytest rather than a new dev dependency.
- typer release notes (`https://typer.tiangolo.com/release-notes/`): 0.19.0,
  "Support `typing.Literal` to define a set of predefined choices"; 0.26.0,
  click vendored, click-specific functionality unsupported.

*Legacy.* Nothing to carry: the legacy tree has no command-line height unit,
no click, and no dependency declaration this increment could reuse.

```
$ git grep -l -iE 'z_unit|z-unit|click|vertical.?unit' legacy-archive -- legacy
legacy-archive:legacy/rasputin/reader.py
$ git grep -l -E "install_requires|requirements" legacy-archive -- legacy
(no output)
```

The one hit is `legacy/rasputin/reader.py@legacy-archive:86`,
`VerticalUnitsGeoKey = 4099`, a GeoTIFF key constant that
`tin_engine/io/geotiff.py` already reads.

## 4. Literal or enum: Literal

Increment 47 proposed `Literal["m", "ft", "us-ft"]` with `typer>=0.19`; this
file's first version proposed a `StrEnum` with `typer>=0.26`. **The Literal
is taken**, because:

- it keeps the value a plain `str`, as today, so `_z_unit`, `_z_words` and
  the MCP default handling see what they see now;
- typing `UNITS` with it makes mypy refuse a row whose spelling is not one of
  the three, where the enum needed a second hand-written list kept in step
  only by a test;
- its floor is sound. My reason for 0.26 was that an old typer can be paired
  with a click it does not parse with (measured: typer 0.12.0 with click 8.5.0
  answers `--z-unit ft` with "Got unexpected extra argument (ft)"). typer
  0.19.0 does not have that fault: measured on a toy app with the Literal
  option, typer 0.19.0 with click 8.5.0 and typer 0.26.0 with click
  uninstalled both accept `us-ft` (the function receives `'us-ft'`) and answer
  `yd` with exit 2 and "Invalid value for '--z-unit': 'yd' is not one of 'm',
  'ft', 'us-ft'."

## 5. What changes

1. `src_python/tin_engine/io/models.py`: beside `ZUnit`, the option's
   spellings `ZSpelling = Literal["m", "ft", "us-ft"]`, and `UNITS` typed
   `dict[int, tuple[ZUnit, float, ZSpelling, str]]`.
2. `src_python/tin_engine/cli.py`: `Z_UNIT` becomes
   `Annotated[ZSpelling | None, typer.Option("--z-unit", help=...)]` with no
   `click_type`; `import click` is removed; `_z_unit` takes
   `ZSpelling | None`. The help text is unchanged.
3. `pyproject.toml` (not counted): `typer>=0.19`; the `mcp` extra adds
   `anyio>=4.9` and `mcp-types>=2.3` (the bounds `mcp` 2.3.0 itself requires;
   the 3.14-only `anyio>=4.10` stays `mcp`'s to impose). `dev` is unchanged:
   its `mcp>=2.3` brings both. **Click is not declared**: after change 2 no
   first-party file imports it.

## 6. The tests that pin it (`@tester`, red first)

No fixture files: every CLI test is an argument list. click converts the given
`--z-unit` before it checks for missing required options (master's
`station-catchments --z-unit furlong`, with no `--dem`, already reaches the
conversion). Use the runner directly, not `refused()`, which appends `--out`
where `station-catchments` has `--out-dir`.

- **T1, the usage error** (parametrised over `mesh`, `catchment`,
  `station-catchments`): `--z-unit furlong` gives exit 2, the result's
  exception is `SystemExit`, and the output names `--z-unit`, `furlong`, `m`,
  `ft` and `us-ft` (Rich may wrap: compare with whitespace removed). Red on
  master: exit 1, the exception a `BadParameter`.
- **T2, the help** (same three): `--help` output, whitespace removed, contains
  `m|ft|us-ft` and not `<function>`. Red on master.
- **T3, the MCP schema** (in `tests/python/test_mcp_schema.py`, which skips
  without `mcp`): for `mesh`, `catchment` and `station_catchments`, the
  `z_unit` property admits exactly `m`, `ft`, `us-ft` and null. Red on master
  (any string).
- **T4, no drift**: the members of `ZSpelling`
  (`typing.get_args`) equal the spellings in `UNITS` (three: `m`, `ft`,
  `us-ft`). Red on master (no `ZSpelling`).
- **T5, the import gate** (new file `tests/python/test_declared_imports.py`,
  product suite, not `harness`): parse `pyproject.toml` with `tomllib`; take
  each requirement's name (the leading `[A-Za-z0-9._-]+`) and normalise it
  (PEP 503). Walk every `.py` under `src_python/tin_engine/` with `ast`
  (every `Import`, every `ImportFrom` with `level == 0`, inside functions and
  `TYPE_CHECKING` blocks too); take the top-level name, skip
  `sys.stdlib_module_names` and `tin_engine`, normalise it. A file under
  `tin_engine/mcp/` may use `dependencies` plus the `mcp` extra; every other
  file only `dependencies`. The assertion lists each offender as
  `path:line name`. Static: imports nothing, needs no installed extra.
  Fixture figures (master `55b9f564`): 70 files, ten third-party names
  (section 2's table); red with five offending lines: `cli.py:58 click`,
  `mcp/schema.py:20 mcp_types`, `mcp/server.py:31 anyio`,
  `mcp/server.py:32 mcp_types`, `mcp/server.py:34 anyio`. A core module
  importing `mcp` must also fail; `@tester` shows it once by a planted line
  in a scratch copy, not in the suite.
- **Happy path, kept**: `tests/python/test_cli_units.py`'s mesh with
  `--z-unit us-ft` (and the refusal of a feet DEM without it) stays green
  unchanged.

No mutation round: nothing here is invariant-critical.

## 7. Left out

- `tools/`, `tests/` and `docs/benchmarks/` imports: not shipped in the
  wheel. `tests/` uses `pytest`, `hypothesis` and `vtk`, all in extras.
- Unused declared dependencies (deptry DEP002): all six core dependencies and
  `imagecodecs` (loaded by tifffile, not imported) are in use.
- A CI job that installs without extras and runs `rasputin --help`: it would
  catch this class of fault at the cost of a C++ build per run; T5 catches it
  statically. Handed back as an idea.
- An upper bound on typer: a policy for Ola, not part of this fix.
- A general check that no option escapes as a crash: no other option is
  parsed through a foreign type (`grep -rn click_type src_python` returns
  only `cli.py:229`); options with a `callback=` raise `typer.BadParameter`
  and exit 2.
- From Ola's private ideas file, not taken up (as 47 noted): "`--stats`
  swallows the next option as its path" and "Two small fixes from 37's code
  review" (`schema._annotation` dropping a zero bound).

## 8. Estimate

Production lines (`CLAUDE.md` §2; TOML is not counted): about 0 net
(`import click` and the `click_type` line removed, one alias added to
`models.py`, three lines changed in place). Tests about 90 lines. Under the
short path's 50-line limit; no C++ file; nothing under *Acceptance*.

**Code review, round 1 (short path), 2026-10-14.** Range `55b9f564..d5d29945`, 0 net production lines (about 0 estimated), @reviewer: CHANGES REQUESTED. Blocking: `help_text` in `tests/python/cli_helpers.py` sets `COLUMNS=200` under `TERM=dumb`, which Rich ignores under `GITHUB_ACTIONS` (80 columns), so the three moved help tests fail on CI (reproduced with `GITHUB_ACTIONS=true`); fold it onto #258's `help_page` (`TERM=xterm`, `NO_COLOR=1`, `COLUMNS`) after #258 merges, as one helper. Accepted pins: T3's two checks, T5's 70-file floor, T1 and T2 in a new file. Held: short path; red before green, green touching no test; no scaffolding; prose matches code. Not pushed.
