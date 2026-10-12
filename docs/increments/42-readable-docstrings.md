# Increment 42 — Readable Python docstrings and comments

Status: **designed** (`@architect`, 2026-10-13, off master `dd33a9ce`), not yet
reviewed. Six PRs (section 8). Every PR but the checker's changes no code: its
docstrings and comments are rewritten, and an AST comparison (section 6) shows
that no code token moved.

## 1. Ola's ask

Ola, 2026-10-12: "So main, I've just had a bit of looking at the Python code.
The documentation in the files heavily references increments, fixes, etc. The
amount of cross referencing is very hard to work with for a human." He approved
this as the first half of his product review (the second half, the generality
audit, is `docs/increments/generality-audit.md`): rewrite the docstrings and
comments of `src_python/tin_engine/` so they say what the code does and why, in
plain words. The history stays in git and `docs/increments/`. A rule for new
code follows.

The reader this is written for: someone who opens `mosaic.py` and has never seen
an increment file.

## 2. Prior art: legacy and literature

*Literature.* The rule follows four published sources and claims nothing new.

- PEP 257, *Docstring Conventions* (Goodger and van Rossum, 2001): a docstring
  states what the object does, its arguments, its return value, its side
  effects and the exceptions it raises. It does not ask for history. Rasputin's
  docstrings keep their house style (a noun phrase where PEP 257 asks for a
  command); this increment changes what they say, not their grammar.
- Google Python Style Guide, section 3.8, *Comments and Docstrings*: comments
  explain the tricky parts and the reasons, and "never describe the code".
- R. C. Martin, *Clean Code* (2008), chapter 4, *Comments*: "Journal comments"
  (a changelog at the top of a file) and "Attributions and bylines" (who
  decided) are listed as bad comments, because source control keeps both.
  `io/repository.py`'s module docstring is a journal comment in this sense
  (example 1, section 5).
- J. Ousterhout, *A Philosophy of Software Design* (2018), chapter 13 (comments
  describe what is not obvious from the code) and section 16.3, *Comments belong
  in the code, not the commit log*. This is the limit of "leave the history to
  git": **the reason for a design choice stays in the code**; only the label of
  the decision (its increment, its ruling number, who made it) leaves. A rewrite
  that deletes a reason along with its label departs from Ousterhout and is a
  defect of the sweep.

Comments that go stale as code changes are a measured problem (Wen, Nagy,
Bavota and Lanza, *A Large-Scale Empirical Study on Code-Comment
Inconsistencies*, ICPC 2019). A label such as "increment 15a" is worse than a
stale fact: it is correct but tells the reader nothing without a second file.

Existing linters check docstring form, not content: ruff's `D` rules
(pydocstyle) check shape and punctuation; `ERA` flags commented-out code; `FIX`
and `TD` flag TODO and FIXME markers. None knows a project's labels, so the
checker of section 6 is a project-specific pattern list. No novelty is claimed.

*Legacy.* The legacy tree has docstrings in Sphinx `:param:` style and no
history references:

```
$ git grep -l -iE 'increment [0-9]' legacy-archive -- legacy | wc -l
       0
```

Nothing is carried across.

## 3. Measured on master `dd33a9ce`

All figures come from `docs/increments/42-probes/history_probe.py`. Its `scan`
output is committed as `docs/increments/42-probes/scan_output.txt`, one line
per match.

**The text in scope.** 65 Python files, holding 3,034 lines of module, class and
function docstrings and comments. These are counted by `ast` (docstring spans)
and `tokenize` (comment tokens). A bare string after an assignment is not a
docstring to `ast` and is not counted.

**History labels.** These are found by the probe's six patterns (section 6),
which run over docstrings and comments only:

```
$ python3 docs/increments/42-probes/history_probe.py scan src_python/tin_engine \
    > docs/increments/42-probes/scan_output.txt
TOTAL 601 {'increment': 147, 'review': 28, 'label': 156, 'ruling-id': 216, 'person': 54}
```

That is 601 matches on 404 lines in 62 of the 65 files. There are no commit
hashes (`hash`: 0). The only hex runs in the package are the WGS 84 constants
in `io/geopackage_write.py`. The matches per file:

```
cli.py 99, feature_input.py 43, mosaic.py 38, io/repository.py 36, dem_input.py 31,
io/geotiff.py 22, io/ply.py 21, target_grid.py 20, io/models.py 20, stats.py 15,
io/vtk_legacy.py 15, sources.py 12, run_record.py 12, final_check.py 12,
features.py 12, crs.py 12, io/__init__.py 9, decompose.py 9, chains.py 9,
io/geopackage.py 8, edge_strip.py 8, domain.py 8, io/cog.py 7, viz/fixtures.py 6,
landcover.py 6, io/gml.py 6, grid_domain.py 6, fetch/plan.py 6, palettes.py 5,
land_cover_fractions.py 5, fetch/run.py 5, catchment.py 5, and 30 files with 1 to 4
```

(Paths under `src_python/tin_engine/`; `cut -f1 docs/increments/42-probes/scan_output.txt`
gives the full list.)

The broader words the brief named were measured too: "fixed in", "used to",
"no longer", "previously", "formerly" (1 match in total), and "was", "now",
"since" (31 lines). Of those 31 lines, 26 are ordinary English ("Seconds since
the clock was made"). Of the 5 that are history, 4 also carry a label the
patterns catch. These words are therefore left out of the checker: they are
for the reviewer to judge, not for a pattern.

**Precision of the patterns.** Every `ruling-id` line without another history
marker was read (100 lines). All of them are design-file IDs (`R5`, `D2`, `Q11`,
`J6`, `B14`). The draft's two false positives, `1e-6` and `1e-4`, are now
excluded as scientific notation, and so are labels inside a
`docs/increments/` path (`15c-geographic-dem.md`). No false positive is left
among the 601 matches that a full read of `scan_output.txt` found.

**A sample of 30, classified.** Every twelfth line of the broader scan (352
lines, including the was/now/since words) was taken and judged by hand:

| Class | Lines | What happens to it |
|---|---|---|
| Pure history: the label adds nothing the sentence lacks | 19 | Delete the label, keep the sentence (`(increment 37, section 6.4)`, `26b:`, `(PR 4, "Lake gauges")`) |
| A reason disguised as history | 7 | Rewrite as the reason. Examples: "because 6b-i deduplicated them" becomes "because `build_scene` deduplicates them"; "that is 16b's bit 7" becomes "that is the feature bit `land_cover`"; "the Q5 refusal" becomes the refusal's condition |
| A pointer a maintainer needs | 1 | `mcp/server.py`'s module docstring pointing at the MCP design. Keep it, in the heading form of section 4 |
| Not history (ordinary "was", a number) | 3 | Untouched |
| Text a test or tool reads | 0 in the sample | See below |

The checker's patterns flag 26 of the 27 history lines in the sample. The one
they miss is `src_python/tin_engine/catchment_core.py@dd33a9ce:1`, "(python-audit.md, section 11)": a bare
file name with no label. Other misses seen while browsing are of the same
kind: a bare increment number in parentheses, `(32)`
(`src_python/tin_engine/feature_input.py@dd33a9ce:8`), or followed by a section, `35, section 3.2`
(`src_python/tin_engine/feature_input.py@dd33a9ce:252`), `# 34, 4.3:` (`src_python/tin_engine/run_record.py@dd33a9ce:250`). The reviewer
catches these. A pattern for a bare number would flag every parenthesised
count.

**What reads this text.**

- *Tests.* `tests/python/test_io_repository.py:461` reads
  `src_python/tin_engine/io/__init__.py`'s docstring and asserts it contains
  `repository.py` and does not contain "Nothing here opens a file.".
  `tests/python/test_cli_mesh_landcover.py:377` asserts `rasputin mesh --help`
  names `land_cover_code` and `rasputin palette`. Typer prints the `mesh`
  function's docstring as its help, so that test reads `cli.py`'s docstring.
  The other `--help` tests (`test_cli_draw.py`,
  `test_cli_constraint_feet.py`, `test_cli_mesh_land_cover_fractions.py`)
  assert option names, which come from `help=` strings, not from docstrings.
- *The product.* Typer turns the docstring of each of `cli.py`'s nine commands
  into `rasputin <command> --help`. `src_python/tin_engine/mcp/schema.py@dd33a9ce:97` turns its first paragraph
  into the MCP tool description. So users see these docstrings, and today they
  print increment numbers: `mesh` has 8 (`increment 12`, `14`, `16`, `15a`,
  `17`, `16c`, `26b`, `25`), and `palette`, `fetch`, `fetch-stations`,
  `catchment` and `station-catchments` have 1 or 2 each. Two help texts that
  are code, not docstrings, print one each: `src_python/tin_engine/cli.py@dd33a9ce:280` (`mcp`, "(increment
  37)") and `src_python/tin_engine/cli.py@dd33a9ce:2395` (`--river`, "increment 29").
- *Tools.* `tools/check_citations.py` scans comments and docstrings in
  `src_python/` for citations and checks that they resolve. No other tool reads
  this text (`grep -n "src_python\|__doc__\|get_docstring" tools/*.py`: the
  hits are directory lists, `argparse` descriptions of the tools' own
  docstrings, and `bench.py` listing files).
- *Citations into the package.* 134 unpinned line citations, in 27 files,
  point into `src_python/` (counted with `check_citations.py`'s own `CITATION`
  pattern and `resolve`, read-only). 49 of them are in
  `docs/benchmarks/2026-10-07/30d-buffer/README.md`. The targets with the most
  are `feature_input.py` (30), `target_grid.py` (25), `dem_input.py` (19) and
  `cli.py` (13). A rewrite moves line numbers. So each of these 134 becomes
  "at risk" on the PR that touches its target, and a citation past the
  shortened file's end becomes "broken", which fails CI. Section 7 deals with
  them.

## 4. The rule

Proposed text. Where it lives is Ola's call (section 10). Until then the sweep
applies it as written here.

> **Docstrings and comments describe the code as it is, for a reader who has
> never seen the project's history.**
>
> - A docstring says what the object does: what it takes, what it returns,
>   what it refuses, and the reason for anything a careful reader would
>   otherwise question. A comment says why a line is as it is, when the code
>   cannot say it.
> - They never carry history: no increment numbers or labels (`increment 29`,
>   `16c`, `23a-2`), no IDs from a design file (`R5`, `D2`, `Q11`), no review
>   or PR names (`audit PR C`, `round 3`), no record of who decided
>   (`Ola's ruling`), no commit hashes, no account of what the code used to
>   do. Git and `docs/increments/` keep that.
> - A reason that came from a decision is written as the reason, never as the
>   decision's label. A value that was chosen rather than derived says so:
>   "a chosen threshold: below 1 mm two tiles are taken to agree".
> - A pointer to a design file is kept only where the derivation or
>   measurement behind a value is too long to repeat. It names the file and a
>   heading, in the form `tools/check_citations.py` verifies:
>   `` `docs/increments/32-landcover-simplify.md`, *5.3 The band in plan, not over the DEM* ``. At most one
>   per module docstring. None in functions.
> - Tool directives (`# type: ignore[...]`, `# noqa: ...`, `# fmt: skip`) are
>   code: they are kept exactly as they are.
>
> Checked by `tools/check_code_history.py`.

## 5. Three before-and-after examples

The "after" text is the shape the sweep aims for. Its facts are checked by the
sweep's writer and its reviewer against the code, not taken from this file.

**Example 1: a changelog at the top of a module**
(`src_python/tin_engine/io/repository.py@dd33a9ce:1-27`, pure history around true
statements). Before (an excerpt):

```
"""A DEM held as tiles on disk: headers first, pixels on request (increment 15a).

`docs/increments/15-dem-mosaic.md` R1-R3 and Ola's Q3 ruling: this is the one
module in `io/` that opens files, and it opens them read-only (`"rb"`). ...

Increment 23a-1 adds windowed loads (`load_window`, only the blocks a window
meets) and the tile cache's read side (`CacheRepository`): a source fetched
by `rasputin fetch` (23a-2) is read from `<cache>/<source>/` alone, offline.
Increment 23a-2 adds the cache's write side (`CacheWriter`), ...
Increment 16b (R2, R3) adds :func:`open_geopackage`: SQLite cannot read from a
Python stream, ...
Increment 29 adds :func:`read_json`, ...
Audit PR C adds :func:`read_text`, ...
Increment 36 adds :func:`create_sqlite`, ...
Increment 26a adds the class raster's reads from the cache ...
"""
```

After:

```
"""A DEM held as tiles on disk: headers first, pixels on request.

The one module in `io/` that opens files; the readers and writers elsewhere in
`io/` take bytes or an open connection from here. Reads open files read-only
(`"rb"`); the two writers are `CacheWriter` and `create_sqlite`.

`TiffDemRepository` stores, it does not choose: which tiles there are
(`footprints`, headers only), one tile's pixels (`load`), or only the blocks a
window meets (`load_window`). Choosing and stitching tiles is
`tin_engine.mosaic`'s, the same arithmetic for any storage.

`CacheRepository` reads what `rasputin fetch` copied to `<cache>/<source>/`,
offline: DEM tiles, a land-cover class raster (`load_class_window`,
`class_strips`) and the day each object was first fetched. `CacheWriter` is the
one writer of that directory.

Other files opened here for other modules: `open_geopackage` (SQLite cannot
read a Python stream, so a GeoPackage is read through a read-only connection),
`create_sqlite` (the GeoPackage writer's file), `read_json` (station and river
files) and `read_text` (a WKT domain).
"""
```

**Example 2: a function contract written as a list of rulings**
(`src_python/tin_engine/cli.py@dd33a9ce:2059-2077`, `_dem_mesh`). Before:

```
    Without ``tolerance`` this is increment 12's R6: z sampled bilinearly at
    the stride grid. With it, increment 14's R9: the stride grid is the start
    mesh, refined against the DEM's nodes; with ``domain`` (already in the
    DEM's CRS, 15b; ``domain_name`` is its file's), increment 16's R3, the
    polygon's rings are, and ``features``' lines (16b). ``min_angle`` > 0
    improves the start's angles first (increment 20), ``gain`` >= 0 with
    20c's soft criterion; ``feet`` inserts
    constraint feet (increment 20b). With ``grid`` and its ``checks`` (15c-2),
    the refined mesh is checked against the source's nodes (D5), after the
    tile is dropped (15e, fix 3). ...
    ``clock`` gets R5's phases.
```

After:

```
    Without ``tolerance``, z is sampled bilinearly at every ``stride``-th
    node. With it, those nodes are only the start mesh, refined at DEM nodes
    until every triangle is within ``tolerance``. With ``domain`` (already in
    the DEM's CRS; ``domain_name`` names its file) the start mesh is the
    polygon's rings and the lines of ``features`` instead.
    ``min_angle`` > 0 adds points to improve the start mesh's angles first, and
    ``gain`` >= 0 makes each such point pay for itself in angle gained.
    ``feet`` puts a point that falls near a constraint onto the constraint.
    ``lines`` and ``steep`` tighten the tolerance near named lines and on steep
    ground, in all three refinement calls.
    With ``grid`` and ``checks`` (a DEM resampled onto a target grid), the
    refined mesh is then checked against the source DEM's own nodes. The tile is
    released first, so the check never holds both in memory.
    ``clock`` gets the phase timings.
```

**Example 3: a comment whose reason hides behind a ruling**
(`src_python/tin_engine/features.py@dd33a9ce:191-194`). Before:

```
# Increment 23c (question A as ruled): a cut run's piece files mark a seam,
# the edge a partition line adds, on bit 9, the lowest the default leaves
# free. A vocabulary of its own, so the uncut path, its fingerprint and the
# classes feature_input builds from the default stay as they are.
```

After:

```
# A cut run's piece files mark a seam (an edge a partition line adds) on bit 9,
# the lowest bit the default vocabulary leaves free. It is a vocabulary of its
# own so that an uncut run's vocabulary, fingerprint and feature classes are
# exactly the default's.
```

## 6. The checker: `tools/check_code_history.py`

The sweep alone fixes today's tree. Without a gate, new code brings the labels
back: the measured count grew with every increment. So: **a checker, as a CI
gate** (default D1).

*What it does.* For each Python file in its scope it reads the comments
(`tokenize`) and the module, class and function docstrings (`ast`), never
string literals, and reports `path:line: rule: matched text` for each match of:

| Rule | Pattern (the probe's) | Today |
|---|---|---|
| `increment` | `increment` + a number and optional letter | 147 |
| `label` | a number and a letter `a`-`f`, optional `-2` or `-ii` (`16c`, `23a-2`, `6b-ii`); not inside a path, not scientific notation (`1e-6`), not a format spec | 156 |
| `ruling-id` | a design ID: one of `R D Q B C S U J M F T` and a number (`R5`, `Q11`, `R4.3`); `D8` excepted, because D8 is also the hydrology flow-direction scheme (the probe does not except it: 3 of its 216 are `D8`, all design IDs today) | 216 |
| `review` | `audit PR`, `PR 4`, `PR C`, `review round`, `round 3` | 28 |
| `person` | `Ola`, `ruled`, `ruling` | 54 |
| `hash` | a backticked hex run of 7-12 characters | 0 |

Exit 1 when anything is reported, 0 otherwise. There is no inline escape: on a
measured 0 false positives in 601, a false positive is reworded, not
silenced.

*Scope.* The scope is a list of paths in `pyproject.toml`, under
`[tool.check_code_history] paths`, read with `tomllib`. It is kept outside the
tool because `tools/check_*` is a governed path (`guard_governance.py` asks
before every write), and each sweep PR's red commit extends the scope. The list
starts empty in PR 1 and ends as `["src_python/tin_engine"]` in PR 6. A path
given on the command line overrides it.

*Same-code mode.* `--same-code BASE` checks, for every `.py` file that differs
between `BASE` and the working tree, that the code is unchanged. The base
version is read with `git show BASE:path`. Two things are compared:

1. `ast.dump` of each version, after removing every statement that is a bare
   string constant (docstrings included). `ast.dump` leaves out line and column
   numbers by default, so reflowed text and moved lines do not count. Comments
   are not in the AST at all.
2. The tool directives in comments (`type: ignore[...]`, `noqa: <codes>`,
   `fmt: ...`, `pragma ...`), in order. The text after a `noqa` code is prose
   and is not compared.

`ast.dump` alone is **not** enough, because a dropped `# noqa` or
`# type: ignore` changes what a gate does without changing the AST. Part 2
closes that gap. Probes, run with the prototype on `io/repository.py` and
`catchment_batch.py`:

| Change | Result |
|---|---|
| a docstring rewritten | same code |
| a comment rewritten | same code |
| one string literal changed (`"rb"` to `"r+b"`) | CODE DIFFERS |
| `# noqa: N818` dropped | CODE DIFFERS |
| only the reason after `# noqa: N818` reworded | same code |
| every package file against itself | 0 differ |

What same-code mode does not show, and how each is covered:

- A docstring's text reaches the user through `__doc__`: the CLI's `--help`
  and the MCP tool descriptions. PR 5 commits the nine `--help` outputs before
  and after into the review, and the suite's help tests run.
- Line numbers move: section 7.
- A bare string removed from a function whose body is only that string leaves
  an empty body. Python refuses that, so the writer must add `...`, and
  same-code mode reports it. Such a docstring is shortened, never deleted.

*Size.* About 90 counted lines: the scan about 35, same-code mode about 30, the
scope, arguments and output about 25. These are production lines under
`CLAUDE.md` §2 (a Python file under `tools/`).

## 7. Citations into the package: PR 0

Each of the 134 unpinned line citations into `src_python/` is pinned to the
commit that last changed the citing line (`git blame --porcelain -L n,n`). That
is the tree its author was reading: a citation of line 120 of a package file
becomes the same path and line pinned to that commit, `<path>@<commit>:120`. `check_citations.py` resolves a
pin with `git cat-file`, so no later rewrite can break or endanger it. This is
done before any sweep PR, so that no sweep PR's review has to re-read 134
quotations. The change is to docs, tests and two benchmark scripts' comments,
in 27 files and 134 lines, and it is mechanical. A citation that
was already stale stays as stale as it was: pinning neither fixes nor hides
it.

## 8. The PR split

Under `CLAUDE.md` §2, only lines that are not comments or docstrings count. So
a sweep PR counts 0 production lines, except PR 5's two help strings, and the
limit is how much a reviewer can read with care. The changed-line figures are
estimates: about twice the flagged lines (rewrapping a paragraph touches its
neighbours), plus whole module docstrings in the changelog style of example 1.

| PR | Scope | Files | Flagged lines | Lines changed (est.) | Counted lines |
|---|---|---|---|---|---|
| 0 | Pin the 134 citations into the package (section 7) | 27 docs, tests and benchmark files | — | about 134 | 0 |
| 1 | `tools/check_code_history.py`, its tests, scope empty, CI step | 2 + workflow + pyproject | — | — | about 90 |
| 2 | The catchment family (`catchment*.py`, `gauge`, `burn`, `outline`, `grow`, `hydrography`, `reference`, `sensitivity`), `fetch/`, `mcp/`, `viz/`, `sources.py`, `__init__.py`, `__main__.py` | 27 | 61 | about 130 | 0 |
| 3 | `io/` | 16 | 103 | about 250 | 0 |
| 4 | DEM input and checks: `dem_input`, `mosaic`, `crs`, `target_grid`, `domain`, `grid_domain`, `raster`, `elevation`, `final_check`, `edge_strip`, `tolerance_field`, `decompose` | 12 | 100 | about 220 | 0 |
| 5 | `cli.py`, and its two help strings (`src_python/tin_engine/cli.py@dd33a9ce:280`, `src_python/tin_engine/cli.py@dd33a9ce:2395`) | 1 | 69 | about 160 | 2 |
| 6 | Features and reports: `chains`, `features`, `feature_input`, `landcover`, `palettes`, `land_cover_fractions`, `border_simplify`, `run_record`, `stats`; the scope becomes the whole package | 9 | 71 | about 160 | 0 |

*Order.* PR 2 goes first because it is small, so Ola sees the rule applied
before the large PRs. `cli.py` comes late because the in-flight code PRs
(threads and memory as settings, increment 38, and the generality-audit fixes)
touch it most. A sweep PR on a package waits for any open PR that touches the
same files (`gh pr list` and `git diff --name-only` per PR), so the sweep never
causes a conflict in someone else's branch.

*Steps of a sweep PR* (the short path of `docs/increments/README.md`: under 50
counted lines and no C++):

1. Red commit: the PR's package is added to the scope in `pyproject.toml`; the
   checker fails, and its output is the worklist.
2. Green commit(s) by `@developer`: docstrings and comments rewritten to the
   rule. No test file is touched. A docstring found to be false while it is
   rewritten is corrected in the same PR and listed in the PR description.
   Example: `mesh`'s first line, "Write a mesh as legacy VTK or as PLY", leaves
   out the GeoPackage output; it is the MCP tool's description too.
3. Checks: the checker on the scope (exit 0), `--same-code` against the red
   commit (exit 0; in PR 5, against the commit with the two help strings), the
   full suite, `mypy`, `ruff`, `check_citations.py`.
4. `@reviewer`: reads the diff as a reader new to the project. The question for
   every rewritten reason: is it about the object the code evaluates
   (`docs/PRINCIPLES.md`, B1)? A rewrite that dropped a reason along with its
   label is a finding.

PR 1 takes the full path: `@tester` writes the checker's failing tests first.
No mutation round: the checker is not an invariant-critical suite.

## 9. Tests for the checker (`@tester`, PR 1)

Every fixture is a module text written inline in the test. The figures below
are the counts each fixture must produce.

1. *Each rule fires, in both places.* One module of 14 lines carries one
   instance of each of the six rules in a docstring and again in a comment: 12
   reports, each naming its rule and line.
2. *Code is not prose.* `help="increment 12"` in a call, `X = "R5 16c"`, and
   an f-string holding `Ola`: 0 reports.
3. *The exclusions.* `1e-6`, `1e-4`, `f"{x:.0f}"` inside a comment,
   `` `docs/increments/15c-geographic-dem.md` ``, `D8 flow directions`: 0
   reports.
4. *The pointer form passes.*
   `` `docs/increments/32-landcover-simplify.md`, *5.3 The band in plan, not over the DEM* `` in a module
   docstring: 0 reports.
5. *Exit status and output.* A clean file gives exit 0 and no output; the file
   from test 1 gives exit 1 and 12 lines of the form `path:line: rule: text`.
6. *Scope.* With `[tool.check_code_history] paths = []`, the dirty file is not
   checked (exit 0). With its directory listed, it is (exit 1). A listed path
   that does not exist is an error naming it.
7. *Same-code mode*, in a temporary git repository with one commit of a
   20-line module: a docstring rewrite gives same; a comment rewrite gives
   same; a string literal changed gives "differs", naming the file; a dropped
   `# noqa: N818` gives "differs"; the reason after the `noqa` code reworded
   gives same; a new file gives "differs" (there is nothing to compare it
   with).

## 10. Defaults taken while Ola was away, and the question

- **D1.** A checker, gated in CI (section 6), rather than the rule alone.
- **D2.** The scope is `src_python/tin_engine/`. Out of scope: `tests/` (717
  label matches; tests pin design IDs on purpose, so a test docstring "R5:
  ..." names the requirement it checks), `tools/`, and the C++ headers
  (`include/`, 84 matches by the same count). The C++ headers are handed back
  as an idea.
- **D3.** The 134 citations are pinned in PR 0, each to the commit that last
  changed its citing line.
- **D4.** Design IDs are deleted and the rule they name is written out. No ID
  becomes a pointer, except where section 4 keeps one.
- **D5.** Sweep PRs take the short path, with "package added to the scope" as
  the red commit.
- **D6.** The order of section 8.
- **D7.** False docstrings found during the sweep are fixed in that PR and
  listed (`docs/increments/README.md`, *Cost constraints*: a documentation
  defect is fixed in the increment's PR).

**The question for Ola: where does the rule live?** Each place is a rule file,
so writing it needs Ola's go. The choices:

- (a) **the Python skill file, `.claude/skills/python-development/SKILL.md`**
  (default): every agent loads it before writing, judging or planning Python;
  it is 35 lines today and the rule adds about 12;
- (b) `@developer`'s persona file, `.claude/agents/developer.md`, read only by
  the agent that writes the code;
- (c) a new `CONTRIBUTING.md` at the root, for human contributors as well, with
  a one-line pointer from the skill.

## 11. Left out

- Docstrings in `tests/` and `tools/` and C++ comments (D2).
- Docstring grammar and format (PEP 257 imperative mood, Sphinx fields):
  unchanged.
- Docstrings added where none exist. The sweep rewrites what is there, and
  same-code mode would report an added bare string only as "same" (it is
  removed before comparing), so adding docstrings is allowed but is not this
  increment's work.

## 12. Ideas file

This takes up part (a) of "Product review: readable code docs, and general
applicability beyond Ola's machine and data". Part (b), the generality audit,
is a separate record. It overlaps this one in one place: comments that name a
data set as if it were the only one (`src_python/tin_engine/mosaic.py@dd33a9ce:32`, "metres for DTM10") are
rewritten in general terms by the sweep.

## Estimate

PR 1: about 90 counted lines. The other PRs: 0, except PR 5's 2. About 920
changed lines of docstrings and comments over five sweep PRs, and 134
pinned citations in PR 0.

## Review

**Design review, round 1, 2026-10-13.** Range `dd33a9ce..61460bef`, 0 counted lines, @reviewer: CHANGES REQUESTED, two blockers. (1) /Users/skavhaug/projects/rasputin/.claude/worktrees/docstrings/docs/increments/42-readable-docstrings.md@61460bef:348-354 claims the directive check closes the gap, but it matches `noqa` only directly after `#`; file-level `# ruff: noqa` (as in `src_python/tin_engine/stats.py@dd33a9ce:12`), `# mypy:`, `# pyright:` are missed, and adding `# ruff: noqa` gives "same code" (probed); match `noqa` anywhere in a comment, add the `ruff:`, `mypy:`, `pyright:` and `isort:` prefixes, and a test 7 fixture. (2) /Users/skavhaug/projects/rasputin/.claude/worktrees/docstrings/docs/increments/42-readable-docstrings.md@61460bef:191-199 sets limits on design-file pointers (heading form, at most one per module docstring, none in functions) and says "Checked by", but no rule checks them; 35 of the 65 lines with `.md` pointers are flagged by no pattern (`src_python/tin_engine/crs.py@dd33a9ce:101-102` is a pointer in a function docstring); add a `pointer` rule with a test, or say the reviewer checks pointers and correct the sample paragraph. Suggestions: fold PR 0 into PR 1 or pin per sweep (it touches `tests/` and runs full CI anyway); `help=` strings are exempt so labels can return to `--help`; `--river` should be `--rivers`; same-code mode looks only in `body`, so a bare string under `else`, `finally` or `except` gives a false "differs". Verified: the probe reproduces (601 labels, 404 lines, 62 of 65 files, 3,034 lines); 134 unpinned citations in 27 files; `mesh` carries 8 labels and two `help=` strings carry labels; `tools/check_*` is governed; same-code mode catches a changed default and a dropped `# noqa: N818`; check_citations exit 0. Question for Ola, default yes: PR 1 also lists the checker among `CLAUDE.md` section 4's gates (a rule-file edit). Lesson: probe a directive check with file-level forms as well as inline ones.
