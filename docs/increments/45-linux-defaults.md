# Increment 45: a cache without setup, general result names, the benchmark on Linux

Status: **designed** (`@architect`, 2026-10-13), not yet reviewed. Written
against master `b1353075`, to be built after increment 44 merges (section 3).
Defaults taken while Ola was away are section 10's questions.

## 1. What Ola asked

The generality audit's fix 6 (`docs/increments/generality-audit.md`, section 6,
item 6: "Cache default, names, tooling on Linux (G1, C5, I1, E1)"); Ola
approved queuing it on 2026-10-13, after fix 4 (increment 44). The audit's
question, in his words of 2026-10-12: "see if we're making a system that is
good only on my computer or datasets, and not generally applicable everywhere
on the globe".

The four rows, as the audit states them:

- **G1** (degrades): the tile cache is `--cache` or `$RASPUTIN_DATA/cache`;
  without either, `fetch` and `mesh --dem <source>` refuse. A first run
  anywhere, and an agent through the MCP server (rasputin's tool interface for
  outside AI agents), must set an environment variable first.
- **C5** (cosmetic): the batch's result fields `nve_in_ours` and `ours_in_nve`
  carry the Norwegian water authority's name (NVE), though the comparison
  works on any reference polygon set.
- **I1** (degrades): `tools/bench.py` and `tools/bench_quick.py` read the
  machine and the power state with macOS commands only (`pmset`, `sysctl`,
  `sw_vers`), so no benchmark runs on Linux.
- **E1** (cosmetic): four of the quick check's six cases read Ola's data
  folder, and a missing input stops the whole check; another machine has no
  check at all.

The speed question (`docs/increments/perf-quick-check.md`, "When it runs"):
**no effect on mesh speed.** Nothing here changes what `rasputin mesh` computes;
it changes where a cache is looked for, two result field names and the
benchmark tools.

## 2. Prior art: legacy and literature

*Literature.*

- **The cache directory** follows the XDG Base Directory Specification,
  version 0.8 (freedesktop.org, read 2026-10-13): "If $XDG_CACHE_HOME is either
  not set or empty, a default equal to $HOME/.cache should be used", and "If an
  implementation encounters a relative path in any of these variables it should
  consider the path invalid and ignore it." Both sentences are taken as they
  stand. The specification is written for Linux desktops; macOS's own
  convention is `~/Library/Caches`. Tools differ on this Mac (measured:
  `ls -d ~/.cache/* ~/Library/Caches/pip ~/Library/Caches/pypoetry` lists
  `~/.cache/gh`, `~/.cache/uv`, `~/Library/Caches/pip` and
  `~/Library/Caches/pypoetry`); uv's documentation states
  "`$XDG_CACHE_HOME/uv` or `$HOME/.cache/uv` on Unix". The design takes the XDG
  rule on every platform, as uv and gh do (question Q1). The `platformdirs`
  package would give the per-platform answer; it is not a dependency
  (`grep -n platformdirs pyproject.toml` is empty) and six lines do not justify one.
- **The power state on Linux** reads `/sys/class/power_supply/<supply>/`, the
  kernel's documented interface (`Documentation/ABI/testing/sysfs-class-power`):
  `type` is one of "Battery", "UPS", "Mains", "USB", "Wireless"; `status` one
  of "Unknown", "Charging", "Discharging", "Not charging", "Full"; `online`
  0 (offline), 1 or 2 (online); `capacity` an integer percentage. A `scope`
  attribute ("Unknown", "System", "Device", added in 2011 by Jeremy
  Fitzhardinge's "power_supply: add SCOPE attribute to power supplies") marks a
  peripheral's battery, a wireless mouse's for example, as "Device". UPower and
  Waybar skip "Device" supplies and treat a missing scope as "System"; this
  design does the same. Departure: UPower also reports a UPS; this design
  ignores UPS entries (a UPS running on its battery is a power cut, not a
  laptop setting; section 10, Q3 records the no-battery default).
- **Cores and memory on Linux** come from the kernel's CPU topology files
  (`Documentation/admin-guide/cputopology.rst`: `topology/physical_package_id`,
  `topology/core_id` per `cpuN`), `/proc/cpuinfo`'s `model name` and
  `processor` lines, `/proc/meminfo`'s `MemTotal` (in kB), and
  `/etc/os-release`'s `PRETTY_NAME` (systemd's `os-release(5)`).

No novelty is claimed.

*Legacy.* Nothing. The commands, run on this tree:

```
$ git grep -l -i -E 'XDG|expanduser|\.cache|power_supply|cpuinfo|pmset|appdirs|platformdirs' legacy-archive -- legacy
(no output; exit 1)
$ git grep -l -i -E 'nve_in|ours_in|reference' legacy-archive -- legacy
legacy-archive:legacy/bindings.cpp
```

The one hit is pybind11's `py::return_value_policy::reference_internal`
(`legacy/bindings.cpp` lines 208, 330, 354, 355), not a result field. The
legacy tree had no tile cache, no reference comparison and no benchmark tool.

## 3. How this meets increment 44

Increment 44 (`docs/increments/44-catalogue-as-data.md` on branch
`worktree-catalogue-data`, in design review) touches the same function:

- 44 adds `cli.data_root(environ) -> Path | None`, "the one reader of
  `RASPUTIN_DATA`", and makes `cache_root` derive from it; the user catalogue
  file is `<data root>/sources.toml`, read only when `RASPUTIN_DATA` is set.
- 45 adds one step **after** `data_root`: with no `--cache` and no data root,
  the cache is the XDG one. `data_root` is unchanged, so the user catalogue file
  is still read only from `$RASPUTIN_DATA` (Q4 asks whether it should also have
  a default place).
- 45 is built on master after 44 merges, so `@developer` changes 44's
  `cache_root`, not today's. If 44 changes the three "no cache" refusals into
  one, 45 rewords that one; otherwise 45 makes them one constant (section 5.1).
- Both add a ROADMAP row after row 43; whichever merges second resolves a
  one-line conflict in `ROADMAP.md`.

## 4. Design

### 4.1 The cache default (G1)

`cache_root(option, environ)` keeps its signature and returns, in order:

1. `--cache`, when given;
2. `<data root>/cache`, when `RASPUTIN_DATA` is set and not empty (44's
   `data_root`);
3. `$XDG_CACHE_HOME/rasputin`, when `XDG_CACHE_HOME` is set to an absolute path;
4. `$HOME/.cache/rasputin`, when `HOME` is set to an absolute path;
5. `None` otherwise: the refusal, as today.

It stays a pure function of its two arguments (no `Path.home()`, which also
consults the password database), so every row of the rule is a unit test.
`RASPUTIN_DATA` wins over XDG because it is rasputin's own, explicit setting;
Ola's setup, with `RASPUTIN_DATA` set, sees no change.

The directory need not exist: the cache writer creates it
(`src_python/tin_engine/io/repository.py@b1353075:366`,
`self.directory.mkdir(parents=True, exist_ok=True)`), and a `mesh --dem
<source>` against a missing directory already says the source is not in the
cache and names the `fetch` that gets it.

Because the default now writes into the home directory without being asked,
`fetch` says where: its summary line gains `into <directory>` (the source's
directory under the cache root). This changes an existing line, adding none.

The three refusals (`mesh --dem <source>`, `mesh --land-cover <source>`,
`fetch`) become one message, held once:
`no cache directory: pass --cache, or set RASPUTIN_DATA (the cache is then
$RASPUTIN_DATA/cache) or XDG_CACHE_HOME (then $XDG_CACHE_HOME/rasputin)`.
The two `--cache` help texts become "Default: $RASPUTIN_DATA/cache, else
$XDG_CACHE_HOME/rasputin, else ~/.cache/rasputin." `INSTALL.md`'s cache
sentence (line 188) says the same.

### 4.2 The result names (C5)

Renamed everywhere in `src_python/tin_engine`, with no alias for the old name:

| today | after |
|---|---|
| `nve_in_ours` | `reference_in_ours` |
| `ours_in_nve` | `ours_in_reference` |

The fields are `reference.Agreement`'s and `reference.Summary`'s, the
`_MEASURES` tuple, `catchment_batch.StationResult`'s and the key tuple that
copies them (`src_python/tin_engine/catchment_batch.py@b1353075:255`), so the
batch's CSV header and summary JSON keys change with them. The stderr line
(`src_python/tin_engine/cli.py@b1353075:2735-2739`) becomes "reference in ours
98.2 %, ours in reference 97.6 %, area ratio 1.003". The module docstring of
`reference.py` says "a reference polygon" where it says "NVE's polygon".

Not renamed: `nve_area_km2`, the area NVE's station list states for a station.
It is a column of the station file format (`io/station_set.py`,
`fetch/nve.py`, the shipped `data/nve_hrd_2025.csv`), not a comparison result;
renaming it changes a file users write (Q5). The committed evidence under
`docs/benchmarks/2026-10-05/` keeps the old names: its scripts read its own
committed CSV and JSON files, which keep their old headers, and no product
code reads a batch's output back (`grep -rn -E 'DictReader|read_csv|json.load'
src_python/tin_engine` finds only the NVE fetch, the cache manifest and the MCP
server's record argument).

### 4.3 Machine and power on Linux (I1)

One protocol in `tools/bench.py`, two implementations, one factory:

```
Runner ──> MacHost(runner) ── sysctl, sw_vers, pmset -g batt ──┐
                                                               ├─> Host.machine() -> Machine
root "/" ──> LinuxHost(root) ── reads text files only ─────────┘    Host.power()   -> Power
make_host(runner, system=sys.platform): "darwin" -> MacHost, "linux" -> LinuxHost,
                                         else BenchError (exit 3, nothing written)
```

`MacHost` is today's `_machine` and `parse_pmset(runner.run(PMSET))`, moved,
unchanged in what it runs. `LinuxHost` starts no process; every file it reads
is under `root`, so a test hands it a directory of fixture text.

**Machine on Linux:**

- `cpu_brand`: the first `model name` value in `<root>/proc/cpuinfo`, else
  the first `Hardware` value, else `"unknown"` (ARM servers print neither).
- `p_cores`: the number of distinct `(physical_package_id, core_id)` pairs over
  `<root>/sys/devices/system/cpu/cpu[0-9]*/topology/`; when that directory is
  absent, the number of `processor` lines in `/proc/cpuinfo`. `e_cores` is 0:
  Intel's hybrid parts are not split into performance and efficiency cores
  (the two counts only identify a machine for comparison, which one number
  does on one machine).
- `memory_bytes`: `MemTotal` in `<root>/proc/meminfo`, times 1024.
- `os_version`: `PRETTY_NAME` in `<root>/etc/os-release`, unquoted, else
  `"Linux"`.

In a container, `/proc/cpuinfo` and the topology describe the host; that is
the machine's identity, which is what the record needs (the threads a run uses
are increment 38's business).

**`Machine` gains `system`** (`"macOS"` or `"Linux"`, default `"macOS"`) and
renames `macos` to `os_version`, which still reads the old key
(`Field(validation_alias=AliasChoices("os_version", "macos"))`). Measured
(`docs/increments/45-probes/probe45.out`): all 174 stored `run.json` files
under `docs/benchmarks/` and the quick baseline load under the new model; the
baseline loads as `system="macOS"`, `os_version="27.0"`, and its dump has the
same keys a new Mac run's has, so `bench_quick._mismatch` sees no spurious
difference. The evidence README's method line says `{system} {os_version}`
where it says `macOS {macos}`, and "read before and after (raw text in
run.json)" where it names `pmset -g batt`.

**Power on Linux**, over the supplies in `<root>/sys/class/power_supply/`
whose `scope` is not `Device` (a missing `scope` counts as `System`):

1. any `Mains`, `USB` or `Wireless` supply with `online` not `0`: `ac`;
2. else any `Battery` with `status` `Discharging`: `battery`;
3. else no mains-type supply and no battery at all: `ac` (a server, a
   desktop, a CI runner or a virtual machine runs on mains; Q3);
4. else every battery `Charging`, `Full` or `Not charging`: `ac` (external
   power, on a laptop that lists no adapter);
5. else `unknown`.

`percent` is the first battery's `capacity` (supplies in name order), else
`None`. `raw` is one line per supply read, `name: key=value, ...` with the
attributes sorted, or `no power supply under /sys/class/power_supply`, so
`run.json` keeps the evidence as it keeps `pmset`'s text today. `UPS` entries
are read into `raw` and decide nothing. `combine_power` is unchanged.

`bench.run` and `bench_quick.run` call `make_host(runner)` once and use
`host.machine()` and `host.power()` where they call `_machine` and `PMSET`
today; `bench_quick` keeps its own `make_host = bench.make_host`, as it keeps
`make_runner`. `bench.build` is already portable (`_core*.so`, CMake).

### 4.4 The quick check on another machine (E1)

- **A case whose input is absent is skipped**, not a refusal: the run prints
  `SKIPPED: <case>: no input <path>` for each and judges the rest. Only the
  present cases' inputs are fingerprinted. A run with every case skipped cannot
  happen while the in-tree cases exist; it is not special-cased.
- **One baseline file per machine and power state**:
  `docs/benchmarks/quick/baseline-<cpu>-<power>.json`, where `<cpu>` is the
  CPU name lower-cased with every run of other characters than `a-z0-9` made one
  `-`, ends stripped. Measured: "Apple M1 Max" gives `apple-m1-max`; "Intel(R)
  Xeon(R) Platinum 8370C CPU @ 2.80GHz" gives
  `intel-r-xeon-r-platinum-8370c-cpu-2-80ghz`. So a second machine's
  `--save-baseline` cannot overwrite Ola's file. The machine is read before the
  stored baseline is looked up (today it is read after the build). `_mismatch`
  is unchanged and still compares the machine and every input.
- **A public case on another continent**, `velhas`: the committed ANADEM
  extract of the upper Rio das Velhas (`tests/fixtures/velhas/anadem_velhas.tif`,
  EPSG:4674, degrees, southern hemisphere), resampled onto a 30 m grid in
  EPSG:31983, `--tolerance 0.5`, threads default, 1 warm-up + 5 runs:

  ```
  args = ["mesh", "--dem", "$TREE/tests/fixtures/velhas/anadem_velhas.tif",
          "--out-crs", "EPSG:31983", "--bbox", "592000", "7845500", "599000", "7854000",
          "--tolerance", "0.5"]
  inputs = ["$TREE/tests/fixtures/velhas/anadem_velhas.tif"]
  ```

  It is the only case on a geographic DEM, the resampled target grid and a
  southern UTM zone; the two other in-tree cases are the Norwegian DTM10 tile.
  Measured on this Mac (AC, 10 threads, the worktree's `_core`, not a Release
  build from `bench.build`; `docs/increments/45-probes/velhas-stats.md`): a
  288 × 239 grid, 99,236 triangles, largest height error 0.49999833 m, clock
  total 0.126 s, process 0.39-0.40 s over three runs; decode the largest phase
  at 28.4 %, so no phase reaches the 40 % hotspot line. It is small: the
  judge's 0.05 s floor lets it catch a slowdown of about 40 % or more, and its
  worth is the code path, not the timing. The extract's whole box at 30 m is
  this case; the CLI has no option for a finer target spacing.
- `--bbox`'s help says "in the DEM's CRS"; with `--out-crs` the box is in the
  target CRS (the case above passes metres for a DEM in degrees, and the
  refusal without it says "a box in the target CRS"). The help says so: "in the
  DEM's CRS, or in --out-crs's when given". A documentation defect fixed here,
  per `docs/increments/README.md`.

Adding a case and renaming the file means the stored baseline no longer
applies: after green, `@perf` runs `python tools/bench_quick.py run
--save-baseline` once on AC, on this branch (nothing it measures changed, so
the branch is master for this purpose), which writes
`baseline-apple-m1-max-ac.json`; `baseline-ac.json` is deleted in the same
commit. This is not `@perf`'s acceptance (no refine or mesh code is touched).

### 4.5 Rule text

Two rule files name `pmset -g batt` as how the power state is read:
`docs/increments/README.md` line 101 and `.claude/agents/perf.md`. Each gains
"on macOS; `/sys/class/power_supply` on Linux". `docs/benchmarks/bench-py.md`
(the tool's design) and `docs/increments/perf-quick-check.md` (the baseline's
name, the cases table with the `velhas` row, the skip) are updated to match.

## 5. Interfaces and types

### 5.1 `src_python/tin_engine/cli.py`

```python
NO_CACHE: Final = "no cache directory: pass --cache, or set RASPUTIN_DATA (...) or XDG_CACHE_HOME (...)"
def cache_root(option: Path | None, environ: Mapping[str, str]) -> Path | None: ...
    # --cache; <data root>/cache; $XDG_CACHE_HOME/rasputin (absolute only);
    # $HOME/.cache/rasputin (absolute only); None
```

### 5.2 `reference.py`, `catchment_batch.py`

Field renames only (4.2); no type changes.

### 5.3 `tools/bench.py`

```python
class Machine(BaseModel):
    cpu_brand: str
    p_cores: int
    e_cores: int
    memory_bytes: int
    system: Literal["macOS", "Linux"] = "macOS"
    os_version: str = Field(validation_alias=AliasChoices("os_version", "macos"))
    python: str; numpy: str; shapely, geos, pyproj, proj: str | None = None

class Host(Protocol):
    def machine(self) -> Machine: ...
    def power(self) -> Power: ...

class MacHost:
    def __init__(self, runner: Runner) -> None: ...
class LinuxHost:
    def __init__(self, root: Path = Path("/")) -> None: ...

def make_host(runner: Runner, system: str = sys.platform) -> Host: ...  # BenchError otherwise
```

`parse_pmset` and `PMSET` stay (MacHost uses them; their tests stand).

### 5.4 `tools/bench_quick.py`

```python
make_host = bench.make_host
def baseline_path(machine: Machine, power: str) -> Path: ...  # QUICK / f"baseline-{slug}-{power}.json"
```

## 6. Tests `@tester` writes first

Figures come from `docs/increments/45-probes/probe45.py`, run as
`.venv/bin/python -I probe45.py <tree>` (output `probe45.out`); its sketches
of the rules are the probe, not the product.

**Isolation first.** With a default in the home directory, a test that clears
`RASPUTIN_DATA` would reach the real `~/.cache/rasputin`. An autouse fixture
in `tests/python/conftest.py` sets `XDG_CACHE_HOME` to a fresh `tmp_path`
directory for every test. The existing tests that expect the refusal
(`tests/python/test_cli_fetch.py@b1353075:280-286`,
`tests/python/test_cli_mesh_cache.py@b1353075:325-326`, `:357-362`,
`:368-383`, `tests/python/test_cli_mesh_land_cover_fractions.py@b1353075:343-357`)
also delete `XDG_CACHE_HOME` and `HOME`, and assert the new message.

**T1, the cache rule** (unit, `cache_root` with a mapping; the probe's nine
rows):

| `--cache` | environment | result |
|---|---|---|
| `/given` | `RASPUTIN_DATA=/data`, `XDG_CACHE_HOME=/xdg`, `HOME=/home/u` | `/given` |
| none | the same | `/data/cache` |
| none | `RASPUTIN_DATA=""`, `XDG_CACHE_HOME=/xdg`, `HOME=/home/u` | `/xdg/rasputin` |
| none | `XDG_CACHE_HOME=/xdg`, `HOME=/home/u` | `/xdg/rasputin` |
| none | `XDG_CACHE_HOME=rel/cache`, `HOME=/home/u` | `/home/u/.cache/rasputin` |
| none | `XDG_CACHE_HOME=""`, `HOME=/home/u` | `/home/u/.cache/rasputin` |
| none | `HOME=/home/u` | `/home/u/.cache/rasputin` |
| none | `HOME=""` | `None` |
| none | nothing | `None` |

`TestW8CacheRoot.test_neither_is_none` becomes the last row.

**T2, the default through the CLI.** With only `XDG_CACHE_HOME=<tmp>` set:
`fetch <key> --bbox ...` against `test_cli_fetch.py`'s fake remote writes the
tiles under `<tmp>/rasputin/<key>/` and its summary line ends with that
directory; `mesh --dem <key>` against a cache written at `<tmp>/rasputin`
exits 0. With `RASPUTIN_DATA`, `XDG_CACHE_HOME` and `HOME` all unset, each of
the three commands refuses as a usage error whose text holds `--cache`,
`RASPUTIN_DATA` and `XDG_CACHE_HOME`. `mesh --help` and `fetch --help` name
`~/.cache/rasputin`.

**T3, the result names.** The existing assertions on `nve_in_ours` and
`ours_in_nve` in `tests/python/test_reference.py` (24 lines) and
`tests/python/test_catchment_batch.py` (10) move to the new names unchanged in
value; the batch's CSV header holds `reference_in_ours,ours_in_reference` in
the old place; the stderr line holds "reference in ours" and "ours in
reference"; and a source scan finds neither old name anywhere under
`src_python/tin_engine` (`nve_area_km2` is allowed).

**T4, Linux machine** (fixture trees under `tmp_path`, built as the probe
builds them):

- x86, 8 `processor` entries, topology pairs (0,0) (0,0) (0,1) (0,1) (1,0)
  (1,0) (1,1) (1,1), `MemTotal: 16374624 kB`, `PRETTY_NAME="Ubuntu 24.04.3
  LTS"`: `cpu_brand` "Intel(R) Xeon(R) Platinum 8370C CPU @ 2.80GHz",
  `p_cores` 4, `e_cores` 0, `memory_bytes` 16767614976, `system` "Linux",
  `os_version` "Ubuntu 24.04.3 LTS".
- ARM, 4 `processor` entries with `CPU implementer`/`CPU part` only and no
  topology directory: `cpu_brand` "unknown", `p_cores` 4.
- No `/etc/os-release`: `os_version` "Linux".
- No subprocess is started: the test's runner raises on any call.

**T5, Linux power** (fixture trees; the probe's six rows, plus the raw text):

| supplies | state | percent |
|---|---|---|
| none | `ac` | None |
| `AC` Mains online 1; `BAT0` Charging 81; `hidpp_battery_0` scope Device, Discharging 40 | `ac` | 81 |
| `AC` Mains online 0; `BAT0` Discharging 64 | `battery` | 64 |
| `BAT0` Full 100, no mains entry | `ac` | 100 |
| `AC` Mains online 0; `BAT0` Unknown 50 | `unknown` | 50 |
| `ups0` UPS Discharging only | `ac` | None |

In the second row `raw` is `AC: online=1, type=Mains` and `BAT0:
capacity=81, status=Charging, type=Battery`, one per line, and holds no
`hidpp_battery_0`. With none, `raw` is `no power supply under
/sys/class/power_supply`.

**T6, the factory.** `make_host(runner, "darwin")` is a `MacHost`,
`make_host(runner, "linux")` a `LinuxHost` rooted at `/`, and
`make_host(runner, "win32")` raises `BenchError` naming `win32`; `bench.py run`
and `bench_quick.py run` with that factory exit 3 and write nothing.
The existing run tests (`tests/python/test_bench.py@b1353075:998`, `use_runner`,
`tests/python/test_bench_quick.py@b1353075:655`) pin `make_host` to `MacHost`
over their fake runner, so they mean the same on CI's Linux as on the Mac;
without it they would read the CI runner's own `/sys`.

**T7, the record format.** Every `run.json` under `docs/benchmarks/` (174
today) and `docs/benchmarks/quick/baseline-ac.json` load as `RunRecord` and
`QuickRecord` (the baseline: `system` "macOS", `os_version` "27.0"); the
evidence README's method line holds `Linux Ubuntu 24.04.3 LTS` for a Linux
machine and no `pmset`. (The suite may read them: `tools/ci_changes.py` classes only Markdown as
prose.)

**T8, the real kernel** (`skipif sys.platform != "linux"`; runs on CI's
Ubuntu legs, never on the Mac): `LinuxHost().machine()` has `p_cores >= 1`,
`memory_bytes > 0`, `system == "Linux"`; `LinuxHost().power().state` is one of
`ac`, `battery`, `unknown`.

**T9, the quick check.**

- With one case's input absent (the existing two-case fixture,
  `tests/python/test_bench_quick.py@b1353075:640`, its second input deleted):
  the output holds `SKIPPED: beta: no input <path>`, alpha runs and is
  judged, the exit code is the verdict's and not 3, and `inputs` holds alpha's
  only.
- `baseline_path` for "Apple M1 Max" and `ac` is
  `QUICK/baseline-apple-m1-max-ac.json`; for the Xeon name above,
  `baseline-intel-r-xeon-r-platinum-8370c-cpu-2-80ghz-ac.json`;
  `--save-baseline` writes there and the next run reads it (the existing
  tests at `:671`, `:698`, `:716`, `:796` move to the new name).
- The shipped cases (`test_the_shipped_cases_are_sections_3_table`) gain
  `("velhas", [0], 1, 5)` at the end, tolerance 0.5, every input under
  `$TREE` and present in the tree.

Mutation testing: not required (no invariant-critical suite; value rules
and tooling).

## 7. How the Linux provider is tested from this Mac, and what CI covers

This Mac has no Linux (no Docker, Podman, Colima, Lima, OrbStack or
Multipass: `which` finds none). The Linux readers take a root directory, so
T4 and T5 run here on fixture trees of `/proc` and `/sys` text written by the
tests, and on CI as well. The fixture text follows the kernel's documented
formats (section 2), not a capture from a real machine. CI's Python legs run
on `ubuntu-latest` (`.github/workflows/main.yaml`) and run `tests/python/test_bench.py`
and `test_bench_quick.py` (they carry no `harness` mark), so T8 is the one
read of a real kernel's files. Not covered anywhere: a full `bench.py run` or
`bench_quick.py run` on Linux, which needs a Linux machine with a C++
toolchain and an hour; nobody has one today.

## 8. Left out

- Splitting Intel hybrid parts into performance and efficiency cores
  (`/sys/devices/cpu_core/cpus`, `/sys/devices/cpu_atom/cpus`).
- Windows (the audit's I2: not a declared platform).
- A default place for 44's user catalogue file (Q4).
- Renaming `nve_area_km2` and the station file's column (Q5).
- The data-folder layout (`../rasputin_data`, audit row G2): the ideas
  file's "Structure for rasputin_data and rasputin_scratch".

## 9. Size

Counted per `CLAUDE.md` section 2 (Python under `src_python/` and `tools/`;
TOML, Markdown and tests not counted):

| part | net lines |
|---|---|
| `cache_root` rule, one refusal constant (three sites lose a line each), help and `fetch` line reworded | +4 |
| C5 renames (changed lines, not added) | 0 |
| `Machine` fields, `Host`, `MacHost` (moved), `make_host` | +18 |
| `LinuxHost.machine` | +20 |
| `LinuxHost.power` | +22 |
| `bench.run` and README line | +2 |
| `bench_quick`: skip, `baseline_path`, host, order | +12 |
| **total** | **about 80** (worst case 120) |

One PR, far under 700. The audit estimated about 120.

## 10. Defaults taken while Ola was away

Each is a question for Ola, with the default the design takes.

- **Q1** Where should rasputin keep downloaded tiles when you have not said?
  Default: `~/.cache/rasputin` (or `$XDG_CACHE_HOME/rasputin`) on both Linux
  and macOS, as uv and gh do on your Mac; the alternative on macOS is
  `~/Library/Caches/rasputin`, which pip uses. `RASPUTIN_DATA` and `--cache`
  still win, so your own setup does not change.
- **Q2** Today, an agent that runs `rasputin fetch` without `RASPUTIN_DATA`
  is refused; after this, it downloads into your home directory (possibly
  gigabytes). Is that acceptable? Default: yes, and `fetch` now prints the
  directory it wrote to.
- **Q3** On a Linux machine with no battery and no power adapter listed (a
  server, a CI runner, a virtual machine), should the benchmark record the
  power as "AC"? Default: yes, so such machines can compare runs; "unknown"
  would make every run there incomparable.
- **Q4** Should your own source list (increment 44's `sources.toml`) also be
  found without `RASPUTIN_DATA`, for example in `~/.config/rasputin/`?
  Default: no, not in this increment.
- **Q5** Should the station list's `nve_area_km2` column (the area NVE states
  for each gauging station) also get a general name? Default: no; it is a
  column of the station file you and NVE's service write, and renaming it
  breaks those files.
- **Q6** Should the quick check keep one stored baseline per machine, named
  after its processor, so another machine never overwrites yours? Default: yes.
