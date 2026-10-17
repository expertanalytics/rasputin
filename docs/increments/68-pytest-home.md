# Increment 68: the test session never reads the user's config or home

**Status: designed** (`@architect`, 2026-10-18, on `9c80b7de`). Short path:
0 counted lines (every changed file is under `tests/`), no C++, nothing the
*Acceptance* section covers. Touches no governed file: not `pyproject.toml`,
not `.claude/`, not `tools/`.

**Plan.** Short path (`docs/increments/README.md`, *Small changes: the short
path*): `@tester` commits the red suite, `@developer` the green step, and one
`@reviewer` round then judges the design, the red step and the green step
together. Citations below are pinned to `9c80b7de` unless they name another
commit.

## 1. The defect

On 2026-10-17 a pytest run of increment 66's branch gave 7843 passed and 50
errors; every error came from Ola's `~/.config/rasputin/sources.toml`, written
in a source schema that branch did not know. With an empty
`XDG_CONFIG_HOME` the same files passed
(`docs/increments/66-coincident-written-vertices.md@9c80b7de:409`).

The suite already has autouse fixtures giving each test an empty
`XDG_CONFIG_HOME` and `XDG_CACHE_HOME`
(`tests/python/conftest.py@9c80b7de:122` and `:138`, on 66's branch since
`f82dd00c`). They did not help, because **they are function-scoped**: pytest
sets up a module-, class- or session-scoped fixture before the
function-scoped ones of the first test that needs it, so such a fixture sees
the environment pytest was started with. Fixtures of that kind run the CLI,
for example `full` (`tests/python/test_cli_station_catchments.py@9c80b7de:110`),
`f2` and `f2_twice` (`tests/python/test_cli_mesh_gpkg.py@9c80b7de:122` and
`:140`) and `velhas` (`tests/python/test_cli_catchment_geographic.py@9c80b7de:80`),
and the CLI reads the user's sources file through `os.environ`. When the fixture
fails, every test that uses it is an error.

Reproduced on `9c80b7de`, with a planted sources file whose one source has
an unknown `kind`:

```
$ XDG_CONFIG_HOME=<dir holding rasputin/sources.toml, kind = "nonsense-kind"> \
    pytest -q --no-cov tests/python/test_cli_station_catchments.py \
    tests/python/test_cli_mesh_gpkg.py tests/python/test_cli_catchment_geographic.py
35 passed, 32 errors in 33.66s
$ XDG_CONFIG_HOME=<empty dir> pytest -q --no-cov <same three files>
67 passed in 10.34s
```

## 2. Prior art: legacy and literature

*Literature.* No published method; this is test hygiene. The practice is the
XDG Base Directory Specification (freedesktop.org, version 0.8), which the
package already follows (`xdg_base`, `src_python/tin_engine/sources.py@9c80b7de:316`),
and pytest's own documentation on fixture scope order ("Higher-scoped
fixtures are executed first", *Fixture instantiation order*). Setting the
environment in `pytest_configure` follows pytest's hook order: it runs once,
before collection and before any fixture. Nothing new is claimed.

*Legacy.* Nothing to carry: the legacy tree has no user-path reads and no
test environment set-up.

```
$ git grep -lnE 'XDG_CONFIG|XDG_CACHE|Path\.home|expanduser|conftest' legacy-archive -- legacy
(no output, exit 1)
$ git grep -l 'monkeypatch' legacy-archive -- legacy
(no output)
```

## 3. Every place the package reads a user path or variable

Found with `grep -rnE "~/|'~|\"~|environ\[|os\.environ|environ\.get|getenv|Path\.home|expanduser|XDG_" --include='*.py' src_python`
on `9c80b7de`. Each read takes its variables from the mapping passed in, and
every caller passes `os.environ`, so setting `os.environ` covers them all.

| Read | Where | Variables |
|---|---|---|
| user sources file | `user_source_candidates`, `xdg_base`, `data_root` (`src_python/tin_engine/sources.py@9c80b7de:325`, `:316`, `:309`); called from `_known_sources` (`src_python/tin_engine/cli.py@9c80b7de:1632`, four commands) and `src_python/tin_engine/source_report.py@9c80b7de:358` | `RASPUTIN_DATA`, `XDG_CONFIG_HOME`, `HOME` |
| tile cache root | `cache_root` (`src_python/tin_engine/cli.py@9c80b7de:1619`, via `xdg_base`) | `RASPUTIN_DATA`, `XDG_CACHE_HOME`, `HOME` |
| `~` in a path a catalogue file names | `_beside` (`src_python/tin_engine/sources.py@9c80b7de:360`, `expanduser`) | `HOME` |
| MCP output root | `build_server` (`src_python/tin_engine/mcp/server.py@9c80b7de:288`, `Path.home()`); reported at `src_python/tin_engine/source_report.py@9c80b7de:396` | `HOME`, `RASPUTIN_MCP_OUT_ROOT` |
| thread count | `resources` (`src_python/tin_engine/resources.py@9c80b7de:102`), from `src_python/tin_engine/cli_kit.py@9c80b7de:60` | `RASPUTIN_THREADS` |
| MCP child processes | `src_python/tin_engine/mcp/server.py@9c80b7de:141` passes `os.environ` on; `tests/python/mcp_driver.py@9c80b7de:86` and `:94` start servers with `dict(os.environ)` | all of the above, inherited |

The same search over `tests/`, `tools/` and `.github/`
(`grep -rnoE 'RASPUTIN_[A-Z_]+' tests tools .github .claude/hooks src_python`)
finds one more environment variable that is read, and it is a **test
control, kept**: `RASPUTIN_INSTALL_TRAP`, set by CI's "Hardening install
trap" step (`.github/workflows/main.yaml@9c80b7de:344`) and read in a
`skipif` at collection time (`tests/python/test_hardening.py@9c80b7de:255`).
Removing it would skip the trap tests and leave that step green having
checked nothing. The other names the search finds are not environment reads:
`RASPUTIN_HARDENING`, `RASPUTIN_EXPECT_HARDENING`, `RASPUTIN_BUILD_PYTHON` and
`RASPUTIN_BUILD_TESTS` are CMake options (`tools/new_worktree.py@9c80b7de:111`),
`RASPUTIN_T6_SENTINEL` is a CMake cache entry a test writes
(`tests/python/test_hardening.py@9c80b7de:167`), and `RASPUTIN_UNRELATED` is
a variable a test sets itself (`tests/python/test_count_loc.py@9c80b7de:832`).

`tools/session_state.py@9c80b7de:57` builds `TRANSCRIPTS` from `Path.home()`
when imported; harness tests that import it are covered by the same change.
Installed plugin packages (`entry_points`) are not a user path and are out of
scope.

## 4. The change

**Where:** `tests/python/conftest.py` only. **Not** in pytest settings:
pytest has no built-in setting for environment variables (that needs the
`pytest-env` plugin, a new dependency, declared in `pyproject.toml`, which is
governed). A session-scoped autouse fixture is also not enough: it runs
after collection, so a read at import or collection time would escape it.

1. A `pytest_configure` hook, the first thing the session does:
   - makes one temporary directory (`tempfile.mkdtemp`, prefix
     `rasputin-test-home-`) holding an empty `home/`;
   - sets `HOME` to it, and `XDG_CONFIG_HOME`, `XDG_CACHE_HOME`,
     `XDG_DATA_HOME`, `XDG_STATE_HOME` to `home/.config`, `home/.cache`,
     `home/.local/share`, `home/.local/state` (not created: an absent
     directory is what a new user has);
   - removes exactly three variables, the package's user settings:
     `RASPUTIN_DATA`, `RASPUTIN_THREADS` and `RASPUTIN_MCP_OUT_ROOT`. Not
     every `RASPUTIN_` name: `RASPUTIN_INSTALL_TRAP` (section 3) must reach
     the trap tests. A user setting added to the package later is added to
     this list in the same change.
2. A `pytest_unconfigure` hook removes the directory. It does not put the
   old values back: the process ends with the session, and no test runs
   pytest in-process (`grep -rn 'pytest.main\|runpytest\|inline_run'
   tests/python` finds nothing), so there is no reader of restored values.
3. `_own_xdg_cache` and `_own_xdg_config` stay: they still give each test a
   fresh directory, which tests that write the cache rely on. Their
   docstrings say the session hook covers fixtures of wider scope.
4. Nothing in the hook imports `tin_engine`: the harness job runs this
   conftest without the package.

Moving `HOME` for the whole session is safe on `9c80b7de`. With `HOME` an
empty directory and `XDG_CONFIG_HOME`, `XDG_CACHE_HOME` and `RASPUTIN_DATA`
unset, `pytest -m harness tests/python` gave 2452 passed, and the product
tests gave one failure, `test_open_core.py::test_no_private_module_imports[tin_engine.snow]`.
That one also fails with the real `HOME`: the main checkout, whose venv ran
the probe, holds a leftover `src_python/tin_engine/snow/` with only
`__pycache__` in it. No test wrote anything into the empty `HOME`.

A test that needs its own values keeps setting them with `monkeypatch`, as
the tests in 8 files that set `HOME` or delete `XDG_CONFIG_HOME` do today
(`grep -rlE 'setenv\("HOME"|delenv\("XDG_CONFIG_HOME"' tests/python`).
Real-data tests find `../rasputin_data` by path
(`tests/python/gpkg_fixtures.py@9c80b7de:35`), not
through `RASPUTIN_DATA`, so removing the variable skips none of them.

## 5. Red tests (`@tester`)

A new file, `tests/python/test_session_home.py`.

1. **`test_a_module_fixture_sees_no_user_sources_file`** (product test). A
   module-scoped fixture records `os.environ["HOME"]`,
   `user_source_candidates(os.environ)` and the three removed variables. The
   test asserts `HOME` differs from the account's home
   (`pwd.getpwuid(os.getuid()).pw_dir`), no candidate exists and none is
   under the account's home, and none of `RASPUTIN_DATA`, `RASPUTIN_THREADS`,
   `RASPUTIN_MCP_OUT_ROOT` is set. Red on any machine today, CI included: the
   module fixture sees the real `HOME`.
2. **`test_the_install_trap_switch_reaches_the_tests`** (`mark.harness`).
   An inner pytest run (as in test 3) with `RASPUTIN_INSTALL_TRAP=1` in the
   outer environment; the inner test asserts it still reads `"1"`. Green
   today, and red against a hook that removes every `RASPUTIN_` name: the
   guard on the reviewer's finding.
3. **`test_a_planted_user_config_does_not_reach_a_module_fixture`**
   (`mark.harness`, imports nothing from `tin_engine`). Runs an inner pytest
   in a subprocess, as `tests/python/test_harness_collection.py@9c80b7de:75` does (`-p conftest`,
   `tests/python` on `PYTHONPATH`), with an outer environment holding a
   planted `HOME` whose `.config/rasputin/sources.toml` and an
   `XDG_CONFIG_HOME` whose `rasputin/sources.toml` are both unreadable TOML,
   and `RASPUTIN_DATA` pointing at a third such file's folder. The inner
   test's module-scoped fixture fails if any of the three files is under its
   `HOME`, `XDG_CONFIG_HOME` or `RASPUTIN_DATA`. Asserts the inner run exits
   0. Red today: the inner run exits 1 with an error. This is the incident's
   shape, planted.
4. **`test_a_monkeypatched_home_still_wins`** (product test, green today
   and after: a guard on the existing behaviour). A test that sets `HOME`
   and deletes `XDG_CONFIG_HOME` with `monkeypatch` gets
   `user_source_candidates` under its own `HOME`, not the session's.

Fixture figures: the planted files are three bytes of invalid TOML each;
the inner run is one test file of one test, with `-o addopts=` (no
coverage), as the harness collection test's inner runs are.

## 6. Left out

- Resolving `../rasputin_data` from the main checkout, so real-data tests
  run in a worktree: a separate idea ("Data-dependent tests never run in a
  worktree"). If it is built through `RASPUTIN_DATA`, it must read the
  variable before this hook removes it.
- Any change to the package: its reads are already pure in the mapping they
  are given.
- `pytest-env` or any setting in `pyproject.toml`.

## 7. Estimate

0 counted lines (`python3 tools/count_loc.py` counts nothing under `tests/`).
About 30 lines in `conftest.py` and 90 in the new test file.
