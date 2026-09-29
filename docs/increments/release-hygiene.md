# Release hygiene — README, INSTALL, and relicensing to MIT

Status: **planned.** Written by `@architect` on 2026-09-29, on branch
`release-docs-licence` off master `340e7a8` (16b and 16c merged, increment 22
not). This branch carries the rewritten `README.md`, the new `INSTALL.md` and
this plan. Nothing below is built yet; every item that changes code, tests,
`LICENSE`, `CLAUDE.md` or `.claude/` is left to the persona named for it.

**The ask** (Ola, 2026-09-29): "we should now rewrite the README, dropping the
references to CGAL and the legacy code. We also need an updated INSTALL, such
that a new user can easily install and run. The license should be permissive,
like MIT or similar, unless we are forced to use anything else. In the latter
case, we might have to change dependencies."

**Ownership.** Ola: "I am Expert Analytics. Owner, Founder, and CEO." Ola may
therefore relicense Ola's own code and the company's.

**Answer in one line.** Nothing forces a copyleft licence: every dependency is
permissive or is weak copyleft used unmodified, and outside `legacy/` the tree
is Ola's except for a handful of boilerplate lines. MIT works once `legacy/`
leaves the tree and three gates that read it are changed.

## Prior art: legacy and literature

This is not an algorithmic increment: it cites no method and claims nothing
new. Its subject is the legacy tree itself, which it removes from the working
tree; history keeps it (section 3).

## 1. Where the licence stands today

- `LICENSE` is the GPL v3 text, added in the repository's first commit
  (`f2f6722`, "Initial commit", Ola, 2018-10-26; the whole diff is `LICENSE`,
  674 lines). The pipeline then was built on CGAL, whose triangulation
  packages are GPL; that is most likely why.
- `pyproject.toml` says MIT: the classifier
  `License :: OSI Approved :: MIT License`, next to
  `license = { file = "LICENSE" }`, which points at the GPL text.
- The built wheel carries the contradiction. A `pip install .` from this branch
  writes `License: GNU GENERAL PUBLIC LICENSE` into the wheel's `METADATA` and
  ships the GPL text as its only licence file
  (`<venv>/lib/python3.13/site-packages/rasputin-0.2.0.dev0.dist-info/licenses/LICENSE`),
  while the classifier says MIT.
- GitHub reports the repository as GPL v3
  (`gh repo view expertanalytics/rasputin --json licenseInfo`).
- CGAL has left the built tree: `python3 tools/check_prohibited_deps.py`
  passes, and it checks imports, includes, declared dependencies and build
  directives outside `legacy/`.

## 2. Who wrote what

Checked on master `340e7a8` with, for every tracked text file outside
`legacy/`, `lib/`, `tests/fixtures/` and `docs/benchmarks/`:

```sh
git blame --line-porcelain -- "$f" | awk '/^author /{a=substr($0,8); c[a]++}
    END {for (k in c) if (k != "Ola Skavhaug") print k, c[k]}'
```

Every line is Ola's except these.

| file | lines, author | what the lines are | action |
|---|---|---|---|
| `CMakeLists.txt` | 7, Magne Nordaas (2019) | `set(CMAKE_CXX_STANDARD 20)`, `set(CMAKE_CXX_STANDARD_REQUIRED ON)`, two `endif()`, three blank lines | none: standard CMake idiom with no expression of its own |
| `.github/workflows/main.yaml` | 12, Sigmund Slang (2020) | `name: CI`, the `on:` push/pull-request block for `master`, `jobs:`, `runs-on: ubuntu-latest`, `steps:`, blank lines | none: GitHub Actions boilerplate |
| `.gitignore` | 4, Sigmund Slang; 1, Vinzenz Gregor Eck | `.DS_Store`, `.env`, `lib/boo*`, `lib/geometry`, `lib/build_*` | delete the three `lib/` lines, which ignore directories that no longer exist (`lib/` holds only `detria/`); keep the other two |
| `README.md` (master) | Magne Nordaas 13, Vinzenz Gregor Eck 6, Olga Silantyeva 8, Stian Lågstad 3 | install notes, blank lines, the two publication references | rewritten on this branch; the publication references are bibliographic facts and stay, reformatted |

`lib/detria/` is third-party (section 5).

**Committed test data.** `git log --follow` shows Sigmund Slang committed these
in March 2020 (`8e30af4`, `5eb4828`, `d920254`); a grep outside `legacy/`
(`git grep -n -I "<name>" -- . ':!legacy'`) shows who reads them:

| file | read by | action |
|---|---|---|
| `tests/fixtures/materials/material.yaml` | nothing (Ola created it in `cc75908`; Sigmund fixed a path in `b4ac0a8`) | delete |
| `tests/fixtures/textures/green_blue.jpg`, `green_red.jpg` | only `material.yaml` | delete |
| `tests/fixtures/tin_archive/ingoya.h5`, `ingoya.xdmf` | only each other | delete |
| `tests/fixtures/corine/0000_4326_corine2018_4e6064_GML.xsd` | nothing: the GML names it in `xsi:schemaLocation`, but `src_python/tin_engine/io/gml.py` never reads a schema | delete, and drop ".xsd" from entry 2 of `tests/fixtures/corine/NOTICE`. The schema is OGR's generated output, not anyone's work, so keeping it would also be fine; it is deleted because nothing reads it |
| `tests/fixtures/corine/0000_4326_corine2018_4e6064_GML.gml` | `tests/python/gpkg_fixtures.py` (`LEGACY_GML`), `test_io_gml.py`, `test_feature_input.py`, `test_cli_mesh_features.py` | keep. It is CORINE data (third-party, credited in `corine/NOTICE`), which Sigmund only committed; Ola's ruling on 16b's Q6 kept it |
| `tests/fixtures/dem_archive/7908_3_10m_z33.tif` | the benchmark tile, and many tests | keep. It is Kartverket DTM10 data, which Sigmund only committed. It has **no credit in the tree**: `dtm10/extract.py` credits its own windows (© Kartverket, CC BY 4.0), but not this tile. Credit it in `NOTICE.md` (section 5) |

**ASK OLA (A1).** Were Magne Nordaas, Sigmund Slang, Vinzenz Gregor Eck, Stian
Lågstad and Olga Silantyeva working for Expert Analytics when they contributed?
If so the company holds those lines anyway. The plan does not depend on the
answer, because the lines it keeps are boilerplate and it deletes the rest, but
the answer belongs in the record. (This is an engineering reading, not legal
advice.)

**A later risk, not this PR's.** The radiation design (branch
`design-terrain-radiation`, TR1) ports `legacy/rasputin/solar_position.h`.
`git blame -C -C` gives that header 681 lines by Ola and 71 by Magne Nordaas.
TR1's design should either rewrite Magne's lines or record Magne's (or the
company's) permission, and should say where the solar-position algorithm itself
came from (it is described there as an SPA implementation checked against NREL's
test values; NREL's own SPA C code has a licence of its own).

## 3. Removing `legacy/`

`legacy/` is 30 tracked files, 5,224 lines, the CGAL-era pipeline by several
authors. It keeps the GPL question open for as long as it ships in the tree, so
it goes; history keeps it. **Before** it is removed, tag the last commit that
holds it (proposed name `legacy-archive`), so a citation or a porting session
can read it back with `git show legacy-archive:legacy/rasputin/reader.py`.
Pushing that tag is a publishing act and needs Ola's yes.

What breaks, found by deleting `legacy/` in a scratch clone of master
`340e7a8` and running every gate and the C++ build:

1. **`tools/check_legacy_imports.py`** exits 1: "legacy/rasputin does not
   exist". Delete the tool, CI's "Archive integrity" step in
   `.github/workflows/main.yaml`, its line in
   `.claude/hooks/gates_after_commit.py` (`GATES`) and its line in `CLAUDE.md`
   section 4.
2. **The C++ build fails.** `tests/cpp/unit/test_solar_position.cpp` includes
   `<rasputin/solar_position.h>` through
   `target_include_directories(test_solar_position PRIVATE "${CMAKE_SOURCE_DIR}/legacy")`
   in `tests/cpp/CMakeLists.txt`: "fatal error: 'rasputin/solar_position.h'
   file not found". Delete the test file and its CMake block. It characterises
   legacy code that nothing in the new tree uses; TR1 will write its own suite
   against the tag.
3. **`tools/check_citations.py`** exits 1 with 173 broken citations. The
   prose cites `legacy/` by line 147 distinct times, in 21 files
   (`git grep -l -E "legacy/[A-Za-z0-9_./-]+:[0-9]+" -- ':!legacy'`), mostly
   the increment records' "Legacy" sections; the tool scans the 19 under
   `docs/` (its default `--paths` are `docs`, `.claude` and `CLAUDE.md`), and
   `project_structure.md` and `parallel_refinement.md` hold the other two. Rewriting 147 citations would
   destroy the records' evidence, so instead the tool resolves a citation whose
   path starts with `legacy/` against the tag
   (`git cat-file -p legacy-archive:<path>`) and reports it broken if the tag
   or the file is missing or the line is past the end. CI's governance job then
   needs the tag: `fetch-depth: 0` on its checkout, or
   `git fetch --depth=1 origin tag legacy-archive` before the step.
4. **`tools/check_prohibited_deps.py`** still passes, but its docstring, its
   `skip` tuple and its failure message describe a `legacy/` exemption that no
   longer has an object. Remove them, and the sentence "`legacy/` is exempt."
   from `CLAUDE.md` section 2.
5. **Configuration and prose that name the directory** (none fails, all go
   stale):
   - `pyproject.toml`: `"legacy"` in ruff's `exclude`;
   - `CLAUDE.md` section 4: "`legacy/` is excluded" after `ruff check .`;
   - `.github/workflows/main.yaml`: the comment above the removed step;
   - `docs/increments/README.md` step 1, the "Legacy" half: its grep runs on
     the tag (`git grep -n <pattern> legacy-archive -- legacy`);
   - `.claude/agents/migration-expert.md` (reads legacy from the tag), and
     `.claude/agents/orchestrator.md` and `architect.md` where they say
     "legacy";
   - `project_structure.md`: the `legacy/` entry in the layout, "What gets
     deleted, eventually", "Existing files to integrate", and the paragraph on
     `legacy/rasputin/` under "Python API surface";
   - `tests/fixtures/corine/NOTICE`: cites `legacy/tests/test_gml_repository.py`;
     make it name the tag;
   - code comments that cite `legacy/` paths, which `check_citations.py` does
     not scan (it reads `.md` only): `include/terrain/core/edge_properties.hpp`,
     `include/terrain/raster/{geometry,raster,sample}.hpp`,
     `tests/python/landcover_fixtures.py`, `test_features.py`,
     `test_always_xy.py`. Prefix each path with the tag.

**ASK OLA (A2).** Does `@migration-expert` stay, reading from the tag, or
retire with `legacy/`? This plan assumes it stays.

## 4. Dependency licences, checked

Each licence below was read from the installed package's own metadata
(`License-Expression`, `License` or classifier) or licence files in its
`dist-info`, in this repository's `.venv` (macOS arm64, Python 3.14) or, for
the build tools, from their `dist-info` in Homebrew and the uv cache. Versions
are the installed ones. The runtime closure was walked from `pyproject.toml`'s
`dependencies` through every `Requires-Dist` without an extra marker.

| package | version | licence | role | bundled native code, and its licence |
|---|---|---|---|---|
| detria | `8aa25f3` | WTFPL or MIT; **MIT elected** (`lib/detria/README.md`) | vendored, compiled into `_core` | — |
| pybind11 | 3.1.0 | BSD-3-Clause | build; its headers compile into `_core` | — |
| scikit-build-core | 1.1.0 | Apache-2.0 | build only | — |
| Catch2 | v3.6.0 | BSL-1.0 | C++ tests only, fetched by CMake | — |
| numpy | 2.5.3 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 | runtime | OpenBLAS; libgfortran, libgcc: GPL-3 **with the GCC Runtime Library Exception**, which puts no terms on programs using them |
| pydantic, pydantic-core | 2.13.5, 2.46.5 | MIT | runtime | — |
| shapely | 2.1.2 | BSD-3-Clause | runtime | **GEOS: LGPL-2.1** (`dist-info/licenses/LICENSE_GEOS`), as a separate shared library |
| pyproj | 3.8.0 | MIT | runtime | PROJ (MIT, `LICENSE_proj`); also libtiff, libjpeg, liblzma, libwebp, libzstd, whose licence files the wheel does not ship (below) |
| typer | 0.27.2 | MIT | runtime | — |
| tifffile | 2026.9.20 | BSD-3-Clause | runtime | — |
| annotated-types, annotated-doc, typing-inspection, rich, markdown-it-py, mdurl | — | MIT | runtime (via pydantic, typer) | — |
| typing-extensions | 4.16.0 | PSF-2.0 | runtime | — |
| pygments | 2.21.0 | BSD-2-Clause | runtime (via rich) | — |
| shellingham | 1.5.4 | ISC | runtime (via typer) | — |
| certifi | 2026.7.22 | **MPL-2.0** | runtime (via pyproj) | — |
| imagecodecs | 2026.8.16 | BSD-3-Clause | extra `codecs` | 48 libraries; `grep -l -i "GNU GENERAL PUBLIC\|GNU LESSER\|GNU LIBRARY"` over its `licenses/` directory finds none |
| vtk | 9.7.0 | BSD | extra `viewer` | pulls matplotlib (its own permissive licence), pillow (MIT-CMU), contourpy (BSD-3), fonttools, pyparsing, six (MIT), and others; none copyleft by metadata |
| pytest, mypy, ruff | 9.1.1, 2.3.1, 0.16.8 | MIT | extra `dev` | — |
| hypothesis | 6.168.1 | **MPL-2.0** | extra `dev` | — |

**Does anything force a copyleft licence? No.**

- GEOS is LGPL-2.1, but rasputin only imports shapely, which loads GEOS as a
  shared library; rasputin neither modifies nor distributes it. The LGPL puts
  duties on whoever distributes GEOS, and permits use from a program under any
  licence. If Expert Analytics ever ships a bundle that includes shapely's
  GEOS (a frozen desktop application, say), that bundle must keep GEOS
  replaceable and offer its source; that is a note for such a bundle, not a
  reason to change dependencies.
- MPL-2.0 (certifi, hypothesis) is copyleft per file, on changes to those
  files only. Rasputin changes neither.
- The GCC Runtime Library Exception exists so that the GPL on libgfortran and
  libgcc does not reach programs that use them.

**Not verified.** Only the macOS arm64 wheels were inspected. The Linux
(manylinux) wheels of shapely, pyproj, numpy and imagecodecs bundle their own
native libraries; `@reviewer`'s Linux install should repeat the grep over their
`dist-info/licenses` and `*.libs/`. pyproj's bundled libtiff, libjpeg and
friends are, by their upstream projects, under permissive licences, but the
wheel does not ship their licence files, so this plan has not read them.

**What rasputin must do for what it redistributes.** The wheel contains
compiled detria (MIT) and pybind11 (BSD-3-Clause) code, and both licences
require their notice in copies, binary ones included. Today the wheel ships
only `LICENSE`. The fix is a `NOTICE.md` at the root (section 5):
scikit-build-core's default `wheel.license-files` is
`["LICEN[CS]E*", "COPYING*", "NOTICE*", "AUTHORS*"]` (its `skbuild_model.py`),
so `NOTICE.md` lands in the wheel's `licenses/` with no configuration.
`@reviewer` checks that it does.

## 5. `NOTICE.md`: third-party code and data credits

A new root file. Its content, to be written by `@developer` (or the main
session, as it is prose):

- **Rasputin itself**: MIT, © Expert Analytics (the holder string as in
  `LICENSE`).
- **Vendored and compiled in**:
  - detria, © Kimbatt, `lib/detria/`, commit `8aa25f3`, MIT by election; the
    MIT text is `lib/detria/LICENSE-MIT.txt`, reproduced in full;
  - pybind11, © Wenzel Jakob, BSD-3-Clause, compiled into `tin_engine._core`;
    its licence text reproduced in full.
- **Fetched for tests, not redistributed**: Catch2 v3.6.0, BSL-1.0.
- **Data in `tests/fixtures/`**:
  - Kartverket DTM10: `dem_archive/7908_3_10m_z33.tif` and `dtm10/`,
    © Kartverket, CC BY 4.0 (as `dtm10/extract.py` and
    `docs/increments/15-dem-mosaic.md` record for the archive);
  - CORINE Land Cover 2018: `corine/`, © European Union, Copernicus Land
    Monitoring Service 2018, EEA; the full attribution and the modifications
    are in `tests/fixtures/corine/NOTICE`, which stays and is referred to.
- **Data in `docs/benchmarks/`**: pictures and domains derived from DTM10 and
  CORINE carry the same two credits. NVE: **no NVE data is committed on
  master**, and nothing names NVE (`git grep -n -i -w NVE -- ':!legacy'` is
  empty). Increment
  22 compares its Bygdin catchment with NVE's published area; if 22 commits any
  NVE geometry, it adds NVE's credit to this file (NVE publishes under NLOD or
  CC BY 4.0 depending on the dataset; **not verified**).
- **Runtime dependencies** are not redistributed (pip installs them, each with
  its own licence); the file names them and points at section 4 of this plan.

## 6. The change to `LICENSE`

Replace the GPL text with the MIT text, with the copyright line

```text
Copyright (c) 2018-2026 Expert Analytics AS
```

**ASK OLA (A3).** The holder string: the company's registered name (is it
"Expert Analytics AS"?), and whether the years start in 2018 (the first commit)
or are left out. `pyproject.toml`'s `authors` stays Ola; the README says
"developed by Expert Analytics", and its Licence section becomes "MIT; see
LICENSE" in the same PR.

`pyproject.toml` keeps `license = { file = "LICENSE" }` and the MIT
classifier, which then agree. (PEP 639's `license = "MIT"` would be neater, but
it forbids the licence classifier and needs a check that scikit-build-core 0.9,
the declared minimum, accepts it; not worth it in this PR.)

## 7. Who does what, in order

One PR, on this branch. One agent at a time.

| step | who | what | files |
|---|---|---|---|
| 0 | main session, **Ola's yes** | create and push the tag `legacy-archive` on master's tip | — |
| 1 | `@tester` (red) | new `tests/python/test_check_citations.py`: a `legacy/` citation resolves through the tag; one past the end of the tagged file is broken; with the tag absent it is reported, not passed. Delete `tests/cpp/unit/test_solar_position.cpp` and its block in `tests/cpp/CMakeLists.txt`. Delete the unused fixtures (section 2). The prose-only `legacy/` mentions in `tests/python/*.py` get the tag prefix | `tests/` |
| 2 | `@developer` (green), each `tools/check_*` edit behind the guard's prompt, so **Ola's yes** | `check_citations.py` resolves `legacy/` through the tag; `check_prohibited_deps.py` loses its exemption; delete `check_legacy_imports.py`; workflow: remove "Archive integrity", fetch the tag in the governance job; `pyproject.toml` ruff `exclude`; `.gitignore`'s three stale lines; tag prefix on the `include/` comments | `tools/`, `.github/`, `pyproject.toml`, `.gitignore`, `include/` (comments) |
| 3 | main session, **Ola's yes** | `git rm -r legacy`; `LICENSE` to MIT with Ola's holder string; `CLAUDE.md` sections 2 and 4; `.claude/hooks/gates_after_commit.py`; `.claude/agents/{migration-expert,orchestrator,architect}.md`; `docs/increments/README.md` step 1 | guarded files |
| 4 | `@architect` or main session | `NOTICE.md`; `project_structure.md`; `tests/fixtures/corine/NOTICE`; README's licence line; a ROADMAP row | docs |
| 5 | `@reviewer` | CI green; INSTALL end to end on a clean Linux and a clean macOS; the built wheel's `dist-info/licenses/` holds `LICENSE` (MIT) and `NOTICE.md`, and its `METADATA` no longer says GPL; the Linux-wheel licence grep of section 4; `python3 tools/check_citations.py` with the at-risk list re-read; after merge, `gh repo view --json licenseInfo` says MIT | — |

**Size.** Production lines by `CLAUDE.md` section 2's count: about 40 in
`check_citations.py`, 10-15 across the workflow, `pyproject.toml` and
`.gitignore`; well under 700. Removed: `tools/check_legacy_imports.py`
(116 lines) and `legacy/` (30 files, 5,224 lines). Tests: about 60 new
lines. Docs: `NOTICE.md` about 80 lines including the two reproduced licence
texts; `LICENSE` 21 lines.

## 8. Open for Ola

- **A1** Contributors' employment (section 2); for the record only.
- **A2** `@migration-expert` stays or retires (section 3).
- **A3** The copyright holder string and years (section 6).
- **A4** The tag name `legacy-archive`, and the yes to push it (section 3).
