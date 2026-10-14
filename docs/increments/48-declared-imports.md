# Increment 48 — Every third-party import is a declared dependency

Status: **designed** (`@architect`, 2026-10-14), short path
(`docs/increments/README.md`, "Small changes: the short path"). Shares its
code change with increment 47 (`--z-unit` usage errors); see section 6.

## 1. The ask

`@perf` reported that a fresh install fails because
`src_python/tin_engine/cli.py` line 58 does `import click` and
`pyproject.toml` does not declare click. Ola approved the fix as its own small
PR on 2026-10-14 (as relayed by the main session; no verbatim quote).

## 2. What the check found

**The failure is real.** Reproduced in the main checkout's venv
(`.venv`, editable install, typer 0.27.2, no `mcp` extra):

```
$ /Users/skavhaug/projects/rasputin/.venv/bin/rasputin --help
  File ".../src_python/tin_engine/cli.py", line 58, in <module>
    import click
ModuleNotFoundError: No module named 'click'
```

**Why click is missing.** typer 0.26.0 dropped click as a dependency and
vendors it as `typer._click` (typer 0.26.0 release notes, "Breaking
Changes", PR #1774: "Typer no longer depends on Click as a third party
dependency, it vendors (includes the source code of) Click. [...]
Click-specific functionality is no longer supported"). PyPI metadata:
typer 0.25.0 requires `click>=8.2.1`; 0.26.0 and 0.27.0 do not.
`pyproject.toml` says `typer>=0.12`, so pip picks 0.27.x.

**Why CI is green.** CI installs `.[dev]`, and `dev` holds `mcp>=2.3`, which
requires `uvicorn`, which requires `click>=7.0`. Checked in the
`catalogue-data` worktree's venv (mcp 2.3.0, click 8.5.0, typer 0.27.3) by
listing every installed distribution that requires click: only uvicorn. A
plain `pip install .` has no click.

**Declaring click is the wrong fix.** The one use of click is
`click_type=click.Choice(...)` on `--z-unit` (`cli.py` line 229). On typer
0.26 and later that passes a real-click type into typer's vendored click,
which the release notes say is no longer supported. Measured on a toy Typer
app in a scratch venv (typer 0.27.3, click 8.5.0) with the same
`click_type=click.Choice(["m", "ft"])` option:

| input | result |
|---|---|
| `--z-unit ft` | accepted, exit 0 |
| `--z-unit yd` | **traceback** (`click.exceptions.BadParameter` from real click, not caught by typer's vendored click), exit 1 |
| `--help` | metavar shown as `<function>` |

That traceback is increment 47's defect. With click declared, the import
succeeds and the traceback stays. The fix is to remove the import and
express the choice the way Typer itself supports: an enum.

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
| **mcp_types** | 2 (`mcp/schema.py`, `mcp/server.py`) | **no**; distribution `mcp-types`, pulled by `mcp==2.3.0` as `mcp-types==2.3.0` |
| **anyio** | 1 (`mcp/server.py`) | **no**; pulled by `mcp` as `anyio>=4.9` (`>=4.10` on Python 3.14) |

So beyond click, two direct imports of the `mcp` package's own dependencies
are undeclared. They cannot fail today (the `mcp` extra always brings them),
but they are imported by name, so they are declared by name (deptry's rule
DEP003, section 3). Each import name above maps to its distribution by
PEP 503 normalisation (lower case, `_` to `-`); no alias table is needed.

## 3. Prior art: legacy and literature

*Literature.* This is packaging practice, not a method; nothing new is
claimed.

- PEP 621, *Storing project metadata in pyproject.toml*: `dependencies` and
  `optional-dependencies` are the project's declared requirements; an
  installer installs nothing else on the project's behalf.
- PEP 503, *Simple Repository API*, section "Normalized Names": the name
  normalisation the test uses to match an import name to a requirement.
- deptry (`https://deptry.com/rules-violations/`), rules DEP001 (an imported
  module has no declared package) and DEP003 (an imported module comes only
  through another package's dependency tree). Click here is a DEP003 case that
  became DEP001 the day typer 0.26 shipped. The test below is those two rules
  for one package, written as a short pytest rather than a new dev
  dependency: deptry is a compiled tool, and one test over 70 files does not
  justify it.
- typer 0.26.0 release notes (`https://github.com/fastapi/typer/releases/tag/0.26.0`):
  click vendored, click-specific functionality unsupported. typer documents
  an `enum.Enum` annotation as its own way to declare a choice
  (`https://typer.tiangolo.com/tutorial/parameter-types/enum/`).

*Legacy.* Nothing to carry: the legacy tree neither uses click nor declares
dependencies in a form this increment could reuse.

```
$ git grep -l -i click legacy-archive -- legacy
(no output)
$ git grep -l -E "install_requires|requirements" legacy-archive -- legacy
(no output)
```

## 4. What changes

1. **`cli.py`: `--z-unit` becomes an enum option.** Remove `import click`.
   Add, beside `Z_UNIT`:

   ```python
   class ZUnitOption(enum.StrEnum):
       """``--z-unit``'s spellings, one per row of ``UNITS``."""
       M = "m"
       FT = "ft"
       US_FT = "us-ft"
   ```

   and `Z_UNIT = Annotated[ZUnitOption | None, typer.Option("--z-unit", help=...)]`
   with no `click_type`. Nothing downstream changes: `_z_unit(option: str | None)`
   and `_z_words` take a `str`, a `StrEnum` member is one, `==` against
   `"ft"` holds and `str()` gives the spelling. The MCP schema
   (`mcp/schema.py`, `_annotation`) reads `param_type.choices` by duck typing,
   so the tool's `z_unit` property stays `Literal["m", "ft", "us-ft"]`.
   Measured on the toy app with this enum: typer 0.26.0 with click
   uninstalled, 0.25.1 with click 8.5.0, and 0.12.0 with click 8.1.8 all
   accept `ft` (the function receives `ZUnitOption.FT`, equal to `"ft"`) and
   answer `yd` with "Invalid value for '--z-unit': 'yd' is not one of 'm',
   'ft', 'us-ft'." (exit 2 measured on 0.26.0 and 0.25.1); typer's `choices` is
   `('m', 'ft', 'us-ft')`; `--help` shows `<m|ft|us-ft>`.
   The members are written out, not built from `UNITS`, because mypy strict
   rejects a functional `Enum(...)` call with a computed second argument;
   test 2 keeps the two in step.
2. **`pyproject.toml`: the floor of typer becomes `typer>=0.26`.** Below
   0.26 typer requires click with no upper bound, and the resolver can pair an
   old typer with a click it does not work with: typer 0.12.0 with click 8.5.0
   (what pip picks for `typer==0.12.0` today) answers the toy's
   `--z-unit ft` with "Got unexpected extra argument (ft)". From 0.26 typer
   carries its own click, so the floor is the first version whose click
   cannot drift.
3. **`pyproject.toml`: the `mcp` extra declares `anyio>=4.9` and
   `mcp-types>=2.3`.** The bounds are the ones `mcp` 2.3.0 itself requires
   (the 3.14-only `anyio>=4.10` stays `mcp`'s to impose). `dev` is unchanged:
   it lists `mcp>=2.3`, which brings both.

Click is **not** declared: after change 1 no first-party file imports it.

## 5. The tests that pin it

New file `tests/python/test_declared_imports.py`, product suite (not
`harness`: it reads `src_python/`).

1. **Every third-party import is declared.** Parse `pyproject.toml` with
   `tomllib`; take each requirement's name (the leading
   `[A-Za-z0-9._-]+`) and normalise it (PEP 503). Walk every `.py` under
   `src_python/tin_engine/` with `ast` (every `Import`, every `ImportFrom`
   with `level == 0`, inside functions and `TYPE_CHECKING` blocks too); take
   the top-level name, skip `sys.stdlib_module_names` and `tin_engine`,
   normalise it. A file under `tin_engine/mcp/` may use `dependencies` plus
   the `mcp` extra; every other file only `dependencies`. The assertion lists
   each offender as `path:line name`. Static: it imports nothing and needs
   no installed extra, so it runs the same in every CI leg.
   Fixture figures (master `55b9f564`): 70 files, ten third-party names
   (section 2's table); on master it fails with five offending lines:
   `cli.py:58 click`, `mcp/schema.py:20 mcp_types`, `mcp/server.py:31 anyio`,
   `mcp/server.py:32 mcp_types`, `mcp/server.py:34 anyio`. After change 3 only
   the first remains; after change 1, none. A core module importing `mcp`
   would also fail, which is the second rule's point; `@tester` shows it by a
   planted line in a copy, not in the suite.
2. **The enum matches the unit table.** `{m.value for m in ZUnitOption} ==
   {row[2] for row in UNITS.values()}` (three spellings: `m`, `ft`, `us-ft`).
   Fails red on `ImportError` until `ZUnitOption` exists.
3. **A wrong `--z-unit` is a usage error.** `CliRunner().invoke(app,
   ["mesh", ..., "--z-unit", "yd"])` exits 2, says "Invalid value for
   '--z-unit'", and `result.exception` is `None` or `SystemExit`. The other
   arguments are `@tester`'s to choose so that parsing reaches `--z-unit`
   (click converts every option before the command body runs, so no DEM is
   read); the red step records the argv used. If increment 47 lands first
   with this pin, 48 drops it (section 6).

No mutation round: no invariant-critical suite here.

## 6. Overlap with increment 47

Increment 47 (designed at the same time on `worktree-z-unit-usage`) is about
`rasputin mesh --z-unit <bad value>` printing a traceback instead of a usage
error. Section 2 shows that traceback and this missing import have **one
cause** (a real-click type handed to typer's vendored click) and **one fix**
(change 1). Whichever increment lands second finds its red test already
green. The main session sequences them: either one PR carries change 1 with
both increments' tests, or 48 lands first and 47 keeps only the pins 48 does
not have (for instance the same check on the two other commands taking
`--z-unit`, `catchment` and `station-catchments`; `Z_UNIT` is used at
`cli.py` lines 967, 2654 and 2917).

## 7. Left out

- `tools/`, `tests/` and `docs/benchmarks/` imports: not shipped in the wheel.
  `tests/` uses `pytest`, `hypothesis` and `vtk`, all in extras.
- Unused declared dependencies (deptry DEP002): all six core dependencies and
  `imagecodecs` (loaded by tifffile, not imported) are in use; no check added.
- A fresh-venv install test in CI (`pip install .` without extras, then
  `rasputin --help`): it would catch this class at the cost of a C++ build per
  run; test 1 catches it statically. Handed back as an idea.
- Upper bounds on typer: typer is 0.x and has broken its API in a minor
  release once already, but a cap is a policy for Ola, not part of this fix.

## 8. Estimate

Production lines (`CLAUDE.md` §2; TOML is not counted): `cli.py` about +3 net
(`import enum` +1, `import click` -1, the class 4, the `click_type` line -1).
Tests about 60 lines. Well under the short path's 50-line limit; no C++ file;
nothing under *Acceptance*.
