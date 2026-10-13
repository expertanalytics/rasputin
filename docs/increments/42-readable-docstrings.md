# Increment 42 — Readable Python docstrings and comments

Status: **PR 1 merged (#246); PR 2 built and review-approved** (red `c75f90ee`, green `b1b43199`, review wordings `c1e41b77`; 0 net counted lines), next its push and CI; PRs 3 to 6 to follow (2026-10-13). Six PRs (section 8). Every PR but the
checker's changes no code: its docstrings and comments are rewritten, and an
AST comparison (section 6) shows that no code token moved.

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

**History labels.** These are found by the probe's seven rules (section 6),
which run over docstrings, comments and the CLI's `help=` strings only:

```
$ python3 docs/increments/42-probes/history_probe.py scan src_python/tin_engine \
    > docs/increments/42-probes/scan_output.txt
TOTAL 668 {'increment': 149, 'pointer': 65, 'review': 28, 'label': 156, 'ruling-id': 216, 'person': 54}
```

That is 668 matches on 441 lines in 63 of the 65 files. Without the `pointer`
rule and the `help=` strings it is 601 matches on 404 lines in 62 files: the
`pointer` rule adds 35 lines no other rule flags, and the two `help=` strings
add 2 `increment` matches. There are no commit hashes (`hash`: 0). The only hex
runs in the package are the WGS 84 constants in `io/geopackage_write.py`. The
matches per file:

```
cli.py 104, feature_input.py 44, mosaic.py 39, io/repository.py 37, dem_input.py 32,
io/geotiff.py 24, io/models.py 21, io/ply.py 21, target_grid.py 21, io/vtk_legacy.py 16,
stats.py 16, crs.py 15, sources.py 14, features.py 13, run_record.py 13, final_check.py 12,
decompose.py 11, io/__init__.py 11, domain.py 10, chains.py 9, io/cog.py 9,
io/geopackage.py 9, edge_strip.py 8, fetch/plan.py 7, io/gml.py 7, landcover.py 7,
catchment.py 6, fetch/run.py 6, grid_domain.py 6, land_cover_fractions.py 6, palettes.py 6,
viz/fixtures.py 6, viz/svg.py 6, catchment_batch.py 5, fetch/http.py 5, gauge.py 5,
io/domain_file.py 5, io/mesh_index.py 5, raster.py 5, reference.py 5, sensitivity.py 5,
and 22 files with 1 to 4
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
among the 601 label matches that a full read of `scan_output.txt` found. The
65 `pointer` matches and the 2 `help=` matches were read too: all true.

**A sample of 30, classified.** Every twelfth line of the broader scan (352
lines, including the was/now/since words) was taken and judged by hand:

| Class | Lines | What happens to it |
|---|---|---|
| Pure history: the label adds nothing the sentence lacks | 19 | Delete the label, keep the sentence (`(increment 37, section 6.4)`, `26b:`, `(PR 4, "Lake gauges")`) |
| A reason disguised as history | 7 | Rewrite as the reason. Examples: "because 6b-i deduplicated them" becomes "because `build_scene` deduplicates them"; "that is 16b's bit 7" becomes "that is the feature bit `land_cover`"; "the Q5 refusal" becomes the refusal's condition |
| A pointer a maintainer needs | 1 | `mcp/server.py`'s module docstring pointing at the MCP design. Keep it, in the heading form of section 4 |
| Not history (ordinary "was", a number) | 3 | Untouched |
| Text a test or tool reads | 0 in the sample | See below |

The checker's rules flag all 27 history lines in the sample. The one the label
rules miss, `src_python/tin_engine/catchment_core.py@dd33a9ce:1`, "(python-audit.md, section 11)", a
bare file name with no label, is flagged by the `pointer` rule (a `.md` name in
a function's or class's docstring, or in a comment). The misses seen while
browsing are bare increment numbers: in parentheses, `(32)`
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
  37)") and `src_python/tin_engine/cli.py@dd33a9ce:2395` (`--rivers`, "increment 29"). The checker
  therefore reads `help=` strings too (section 6), so labels cannot come back
  into `--help` through them.
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

The text, with Ola's rulings of 2026-10-13 (section 10): it lives in the
Python skill file, and outside bug-report links are allowed as below.

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
>   `` `docs/increments/32-landcover-simplify.md`, *5.3 The band in plan, not over the DEM* ``, on one
>   line. At most one per module docstring. None in classes, functions,
>   comments or help text.
> - The same holds for the help text the CLI prints (`help=` strings).
> - A link to an outside bug report (an issue in NumPy, PROJ, GitHub) is
>   allowed where the code works around that bug, and only with the reason
>   written beside it, so the comment still reads if the link dies.
> - Tool directives (`# type: ignore[...]`, `# noqa: ...`, `# ruff: noqa`,
>   `# mypy: ...`, `# fmt: skip`) are code: they are kept exactly as they are.
>
> Checked by `tools/check_code_history.py`, except two things the reviewer
> judges: whether a sentence still describes the code, and history written
> without a label ("used to", a bare increment number).

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
(`tokenize`), the module, class and function docstrings (`ast`), and the string
value of every `help=`, `short_help=` or `epilog=` keyword argument (`ast`; the
constant parts of an f-string, reported at the line where the value starts).
It reads no other string literal. It reports `path:line: rule: matched text`
for each match of:

| Rule | Pattern (the probe's) | Today |
|---|---|---|
| `increment` | `increment` + a number and optional letter | 149 (2 in `help=` strings) |
| `label` | a number and a letter `a`-`f`, optional `-2` or `-ii` (`16c`, `23a-2`, `6b-ii`); not inside a path, not scientific notation (`1e-6`), not a format spec; a label that ends a sentence (`as in 16c.`) counts, a full stop followed by a word character (`16c.txt`) does not | 156 |
| `ruling-id` | a design ID: one of `R D Q B C S U J M F T` and a number (`R5`, `Q11`, `R4.3`); `D8` excepted, because D8 is also the hydrology flow-direction scheme (the probe does not except it: 3 of its 216 are `D8`, all design IDs today) | 216 |
| `review` | `audit PR`, `PR 4`, `PR C`, `review round`, `round 3` | 28 |
| `person` | `Ola`, `ruled`, `ruling` | 54 |
| `hash` | a backticked hex run of 7-12 characters | 0 |
| `pointer` | a `.md` file name (`check_citations.py`'s `MD_PATH`) anywhere but the one allowed place: the first one in the module docstring, in the heading form `` `<path>.md`, *<heading>* `` on one line (`check_citations.py`'s `HEADING_CITATION`, which also checks the heading exists) | 65, on 65 lines in 51 files; 35 of those lines no other rule flags |

*The label rule at a full stop* (decided after the green commit). As built,
the `label` pattern ends in the lookahead `(?![\w.]|-\d)`, so a label followed
by a full stop is missed: "as in 16c." gives no report. It is changed to
`(?!\w|\.\w|-\d)`: a full stop that ends a sentence no longer hides a label,
and one that begins a file suffix (`16c.txt`) or a decimal still does. Measured
with the changed pattern against the checker's own `prose` on this branch: 0
added matches in `src_python/tin_engine` and `tools/`, 9 in `tests/`, all
true labels (7 of them after the word "increment", which the `increment` rule
already flags; the 2 new ones are "Python after 5c." and "before 15c-2.").
No false positive. The probe keeps its old pattern, as the record of the
measurements above.

Exit 1 when anything is reported, 0 otherwise. There is no inline escape: on a
measured 0 false positives in 668, a false positive is reworded, not
silenced. (Every `pointer` hit today is a true one: none of the 65 is in the
heading form, so each breaks the rule as written.)

*Why a `pointer` rule rather than the reviewer.* The rule of section 4 sets
three limits on pointers, and 35 of the 65 lines holding one carry no label
another rule would flag (`src_python/tin_engine/crs.py@dd33a9ce:101-102` is a pointer in a function
docstring). Left to the reviewer, a limit the checker does not see would return
the way the labels did. The rule costs about 10 lines and has a measured 0
false positives. It copies `check_citations.py`'s two patterns rather than
importing them, so neither gate depends on the other.

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
   string constant (docstrings included) from every statement list: `body`,
   `orelse` and `finalbody`, so a bare string under `else`, `except`, `finally`
   or a loop's `else` is removed like one under `def`. `ast.dump` leaves out
   line and column numbers by default, so reflowed text and moved lines do not
   count. Comments are not in the AST at all.
2. The tool directives in comments, in order. A directive is matched anywhere
   in a comment, not only right after `#`: `noqa` with its codes,
   `type: ignore[...]`, `fmt: ...`, `pragma ...`, and a comment's text from a
   `ruff:`, `mypy:`, `pyright:` or `isort:` prefix to its end (so
   `# ruff: noqa: RUF001`, as in `src_python/tin_engine/stats.py@dd33a9ce:12`, is compared whole). The
   text after a `noqa` code list is prose and is not compared. A prose comment
   that merely contains the word `noqa` counts as a directive too; that errs
   towards "differs", the safe side.

`ast.dump` alone is **not** enough, because a dropped `# noqa` or
`# type: ignore` changes what a gate does without changing the AST. Part 2
closes that gap. Probes, run with the prototype on `io/repository.py`,
`catchment_batch.py`, `stats.py` and a 12-line `try`/`for` module:

| Change | Result |
|---|---|
| a docstring rewritten (`io/repository.py` line 1) | same code |
| two comments rewritten (`io/repository.py` lines 187 and 352) | same code |
| one string literal changed (`"rb"` to `"r+b"`) | CODE DIFFERS |
| `# noqa: N818` dropped | CODE DIFFERS |
| only the reason after `# noqa: N818` reworded | same code |
| `# ruff: noqa` added at the top of `io/repository.py` | CODE DIFFERS (the round-1 prototype said "same code") |
| `# mypy: ignore-errors`, `# pyright: basic` or `# isort: skip_file` added | CODE DIFFERS, each |
| `# ruff: noqa: RUF001` dropped from `stats.py`, or a code added to it | CODE DIFFERS, each |
| a bare string added under `except`, `else`, `finally` and a `for`'s `else` | same code (the round-1 prototype said "differs") |
| a function's only docstring replaced by `...` | CODE DIFFERS |
| every package file against itself | 0 differ |

What same-code mode does not show, and how each is covered:

- A docstring's text reaches the user through `__doc__`: the CLI's `--help`
  and the MCP tool descriptions. PR 5 commits the nine `--help` outputs before
  and after into the review, and the suite's help tests run.
- Line numbers move: section 7.
- A changed `help=` string is a changed string literal, so same-code mode
  reports it: PR 5's two are made in a commit of their own, and the sweep is
  compared against that commit.
- A bare string removed from a function whose body is only that string leaves
  an empty body. Python refuses that, so the writer must add `...`, and
  same-code mode reports it. Such a docstring is shortened, never deleted.

*Size.* About 105 counted lines: the scan about 35, the `help=` strings and the
`pointer` rule about 15, same-code mode about 30, the scope, arguments and
output about 25. These are production lines under
`CLAUDE.md` §2 (a Python file under `tools/`).

## 7. Citations into the package: pinned in PR 1

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

*Where.* The pins are a commit of PR 1. The design had them first, by
`@tester`, ahead of the red commit; as built they came after green, by
`@architect` for the docs and benchmark files (section 13). The pinning script is committed as
`docs/increments/42-probes/pin_citations.py`, so the reviewer can rerun it and
see no further change. The round-1 design had them as a PR 0 of their own.
Folding them in saves one PR, one CI run and one review round, at no cost: PR
1 touches `tests/` and runs the full CI anyway, and it merges before every
sweep PR, which is the one order the pins need. Pinning per sweep PR was the
other choice; it was not taken because it spreads the same 134 mechanical
lines over five reviews instead of one.

A citation into the package added by another PR after PR 1 merges shows up in
the sweep PR's `check_citations.py` run, as at risk or broken; that sweep PR
pins it the same way.

## 8. The PR split

Under `CLAUDE.md` §2, only lines that are not comments or docstrings count. So
a sweep PR counts 0 production lines, except PR 5's two help strings, and the
limit is how much a reviewer can read with care. The changed-line figures are
estimates: about twice the flagged lines (rewrapping a paragraph touches its
neighbours), plus whole module docstrings in the changelog style of example 1.

| PR | Scope | Files | Flagged lines | Lines changed (est.) | Counted lines |
|---|---|---|---|---|---|
| 1 | The 134 citations into the package pinned (section 7, first commit); `tools/check_code_history.py`, its tests, scope empty, CI step | 27 docs, tests and benchmark files; 2 + workflow + pyproject | — | about 134 pins | about 105 |
| 2 | The catchment family (`catchment*.py`, `gauge`, `burn`, `outline`, `grow`, `hydrography`, `reference`, `sensitivity`), `fetch/`, `mcp/`, `viz/`, `sources.py`, `__init__.py`, `__main__.py` | 27 | 79 | about 160 | 0 |
| 3 | `io/` | 16 | 112 | about 270 | 0 |
| 4 | DEM input and checks: `dem_input`, `mosaic`, `crs`, `target_grid`, `domain`, `grid_domain`, `raster`, `elevation`, `final_check`, `edge_strip`, `tolerance_field`, `decompose` | 12 | 103 | about 230 | 0 |
| 5 | `cli.py`, and its two help strings (`src_python/tin_engine/cli.py@dd33a9ce:280`, `src_python/tin_engine/cli.py@dd33a9ce:2395`) | 1 | 74 | about 170 | 2 |
| 6 | Features and reports: `chains`, `features`, `feature_input`, `landcover`, `palettes`, `land_cover_fractions`, `border_simplify`, `run_record`, `stats`; the scope becomes the whole package | 9 | 73 | about 165 | 0 |

Flagged lines are the lines of `scan_output.txt` in each PR's files; they sum
to 441.

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

1. *Each rule fires, in both places.* One module carries one instance of each
   of the seven rules in a function docstring (seven lines) and again in seven
   comments: 14 reports, each naming its rule and line.
2. *Code is not prose; help text is.* `X = "R5 16c"`, `print("increment 12")`
   and an f-string holding `Ola`: 0 reports. `typer.Option(help="increment
   12")`, and `@app.command(help="Serve." f"{x}" " (increment 37)")` written
   over two lines: 2 reports, each at the line where the `help=` value starts.
3. *The exclusions.* `1e-6`, `1e-4`, `f"{x:.0f}"` inside a comment,
   `` `docs/increments/15c-geographic-dem.md` ``, `D8 flow directions`: 0
   reports. *Follow-up (section 6, the full stop):* a comment
   `# as in 16c.` gives 1 `label` report; a comment `# read 16c.txt` and a
   comment `# 3f.5` give 0.
4. *Pointers.* `` `docs/increments/32-landcover-simplify.md`, *5.3 The band in plan, not over the DEM* `` on
   one line of a module docstring: 0 reports. Add, in the same module, a
   second pointer in the heading form in the module docstring, one there not
   in the heading form (`` See `x.md` section 3. ``), one in a function
   docstring, one in a comment and one in a `help=` string: 5 `pointer`
   reports, at those five lines. (The prototype gives 4 on the first four; it
   has no `help=` case in that fixture.)
5. *Exit status and output.* A clean file gives exit 0 and no output; the file
   from test 1 gives exit 1 and 14 lines of the form `path:line: rule: text`.
6. *Scope.* With `[tool.check_code_history] paths = []`, the dirty file is not
   checked (exit 0). With its directory listed, it is (exit 1). A listed path
   that does not exist is an error naming it.
7. *Same-code mode*, in a temporary git repository with one commit of a
   20-line module: a docstring rewrite gives same; a comment rewrite gives
   same; a string literal changed gives "differs", naming the file; a dropped
   `# noqa: N818` gives "differs"; the reason after the `noqa` code reworded
   gives same; a new file gives "differs" (there is nothing to compare it
   with). Directives anywhere in a comment, each its own case: `# ruff: noqa`
   added as the first line gives "differs"; the module's
   `# ruff: noqa: RUF001` dropped gives "differs"; `# mypy: ignore-errors`,
   `# pyright: basic` and `# isort: skip_file` added give "differs" each.
   Bare strings beyond `body`: a module holding `try`/`except`/`else`/
   `finally` and a `for` loop with an `else`, with a bare string added at the
   top of each of those four blocks, gives same.

## 10. Defaults taken while Ola was away, and the questions

- **D1.** A checker, gated in CI (section 6), rather than the rule alone.
- **D2.** The scope is `src_python/tin_engine/`. Out of scope: `tests/` (717
  label matches; tests pin design IDs on purpose, so a test docstring "R5:
  ..." names the requirement it checks), `tools/`, and the C++ headers
  (`include/`, 84 matches by the same count). The C++ headers are handed back
  as an idea.
- **D3.** The 134 citations are pinned in PR 1's first commit, each to the
  commit that last changed its citing line (section 7).
- **D4.** Design IDs are deleted and the rule they name is written out. No ID
  becomes a pointer, except where section 4 keeps one.
- **D5.** Sweep PRs take the short path, with "package added to the scope" as
  the red commit.
- **D6.** The order of section 8.
- **D7.** False docstrings found during the sweep are fixed in that PR and
  listed (`docs/increments/README.md`, *Cost constraints*: a documentation
  defect is fixed in the increment's PR).
- **D8.** The pointer limits are checked by the checker's `pointer` rule, not
  left to the reviewer (section 6), and the checker reads `help=` strings.

**Question 1 for Ola: where does the rule live?** Each place is a rule file,
so writing it needs Ola's go. The choices:

- (a) **the Python skill file, `.claude/skills/python-development/SKILL.md`**
  (default): every agent loads it before writing, judging or planning Python;
  it is 35 lines today and the rule adds about 12;
- (b) `@developer`'s persona file, `.claude/agents/developer.md`, read only by
  the agent that writes the code;
- (c) a new `CONTRIBUTING.md` at the root, for human contributors as well, with
  a one-line pointer from the skill.

**Question 2 for Ola: should PR 1 also add the new checker to the list of
gates in `CLAUDE.md`, section 4?** That list is what agents run before they
commit, so a gate missing from it is found only by CI. `CLAUDE.md` is a rule
file, so this needs Ola's go. Default: yes, one line under *Governance gates*,
`python3 tools/check_code_history.py`.

**Ola's rulings, 2026-10-13.** Question 1: (a), the Python skill file, with one
exception added to the rule, a link to an outside bug report where the code
works around that bug, with the reason beside it (section 4). Ola, on being
asked: "Good topic to discuss further, in human generated code, what is
considered state-of-the-art here?"; then, to the default with that exception:
"yes to 1 - 3". Question 2: yes, the checker joins `CLAUDE.md`'s gates in PR 1.

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

## 13. As built (PR 1)

- *Red*, `769ae82a` (`@tester`): 38 tests in
  `tests/python/test_check_code_history.py`, section 9's seven cases plus the
  upper-case `NOQA` case the round-2 review asked for.
- *Green*, `539cf514` (`@developer`): `tools/check_code_history.py`, 128
  counted lines against the estimate of about 105
  (`python3 tools/count_loc.py 769ae82a 539cf514`); the scope
  `[tool.check_code_history] paths = []` in `pyproject.toml`; a CI step. Two
  choices beyond the design: `D8` is excepted from `ruling-id` (section 6 says
  so, the probe does not), and same-code mode reads
  `git diff --name-only BASE -- '*.py'`, so an untracked file is not compared.
- *The rule and the gate*: section 4's text in
  `.claude/skills/python-development/SKILL.md`, section 6; the checker under
  *Governance gates* in `CLAUDE.md` section 4.
- *The pins*: `docs/increments/42-probes/pin_citations.py` found 176 unpinned
  citations into `src_python/` on this branch, not 134: the 42 more are in
  `docs/increments/38-threads-and-memory.md`, which merged after `dd33a9ce`.
  163 are pinned, in 22 docs and benchmark files. 13 are left, outside the
  docs write limit: 10 in `tests/python/` (`test_edge_strip.py`,
  `test_feature_input.py`, `test_features.py`, `test_grow.py`,
  `test_viz_svg.py`), and 3 in `src_python/tin_engine/features.py` lines
  42-45, which cite `viz/svg.py`. Those 13 were pinned by `@tester`
  (`ae34c441`, the 10 in tests) and `@developer` (`8c47988f`, the 3 in
  `features.py`); the script now lists none.
- *The label rule's full stop* (section 6): red `1719aae5` (`@tester`, four
  tests), green `8c47988f` (`@developer`, the lookahead `(?!\w|\.\w|-\d)`);
  44 tests pass.
- *Code review round 1*: the CI step moves `.github/workflows/main.yaml` down
  six lines, so the six citations of its lines 261 and 306 in
  `37-mcp-server.md` (two), `cpp-audit.md`, `h11-ci-path-filter.md` and
  `h16-harness-fixes.md` (two) are pinned, each to a revision where the line
  holds what it says: the two in `37-mcp-server.md` to master `ce8c364f`
  (line 261, the `.[dev]` install); `h11-ci-path-filter.md` and both in
  `h16-harness-fixes.md` to `95b68da` (lines 306-307, the job `CI result`);
  `cpp-audit.md`'s line 552 to `1035690`.

## Estimate

PR 1: about 105 counted lines, and 134 pinned citations. The other PRs: 0,
except PR 5's 2. About 995 changed lines of docstrings and comments over five
sweep PRs.

## Review

**Design review, round 1, 2026-10-13.** Range `dd33a9ce..61460bef`, 0 counted lines, @reviewer: CHANGES REQUESTED, two blockers. (1) /Users/skavhaug/projects/rasputin/.claude/worktrees/docstrings/docs/increments/42-readable-docstrings.md@61460bef:348-354 claims the directive check closes the gap, but it matches `noqa` only directly after `#`; file-level `# ruff: noqa` (as in `src_python/tin_engine/stats.py@dd33a9ce:12`), `# mypy:`, `# pyright:` are missed, and adding `# ruff: noqa` gives "same code" (probed); match `noqa` anywhere in a comment, add the `ruff:`, `mypy:`, `pyright:` and `isort:` prefixes, and a test 7 fixture. (2) /Users/skavhaug/projects/rasputin/.claude/worktrees/docstrings/docs/increments/42-readable-docstrings.md@61460bef:191-199 sets limits on design-file pointers (heading form, at most one per module docstring, none in functions) and says "Checked by", but no rule checks them; 35 of the 65 lines with `.md` pointers are flagged by no pattern (`src_python/tin_engine/crs.py@dd33a9ce:101-102` is a pointer in a function docstring); add a `pointer` rule with a test, or say the reviewer checks pointers and correct the sample paragraph. Suggestions: fold PR 0 into PR 1 or pin per sweep (it touches `tests/` and runs full CI anyway); `help=` strings are exempt so labels can return to `--help`; `--river` should be `--rivers`; same-code mode looks only in `body`, so a bare string under `else`, `finally` or `except` gives a false "differs". Verified: the probe reproduces (601 labels, 404 lines, 62 of 65 files, 3,034 lines); 134 unpinned citations in 27 files; `mesh` carries 8 labels and two `help=` strings carry labels; `tools/check_*` is governed; same-code mode catches a changed default and a dropped `# noqa: N818`; check_citations exit 0. Question for Ola, default yes: PR 1 also lists the checker among `CLAUDE.md` section 4's gates (a rule-file edit). Lesson: probe a directive check with file-level forms as well as inline ones.

**Design review, round 2, 2026-10-13.** Range `423223aa..aa70c016`, 0 counted lines, @reviewer: APPROVED. Verified by rerunning the probe: 668 matches on 441 lines in 63 files, 65 pointer hits, identical to `scan_output.txt`; three pointer hits read, all real; a file-level ruff directive, a changed ruff code and a dropped `type: ignore`, `pragma` or `fmt: skip` all give "differs"; bare strings under `except`, `else`, `finally`, a loop's `else` and a `match` case give "same", a real change under them still "differs"; `--rivers` fixed; check_citations exit 0. Suggestion for the build: ruff also obeys an upper-case `NOQA`; make the directive match case-insensitive (/Users/skavhaug/projects/rasputin/.claude/worktrees/docstrings/docs/increments/42-readable-docstrings.md@aa70c016:376-383) and add a test 7 fixture.

**Code review, round 1, 2026-10-13.** Range `origin/master..8c47988f`, 128 counted lines, @reviewer: CHANGES REQUESTED, two blockers. (1) The new CI step moves `.github/workflows/main.yaml` down 6 lines, so 6 unpinned citations into it now land on the wrong line. (2) The increment file's Status and section 13 are stale. Checked and holding: the SKILL.md rule matches section 4 with the bug-report bullet and no history, the CLAUDE.md gate line, all 176 pins against `git blame`, the empty scope green in a clean environment. Ideas: the `hash` rule misses `path@hash:line` pins (three in `features.py`, for PR 6 to remove); test 94's message "does not exist yet" reads like the red step.

**Code review, round 2, 2026-10-13.** Range `origin/master..7acac12a`, 128 counted lines, @reviewer: CHANGES REQUESTED: h11:451, h16:1444/:1505 and cpp-audit:552 pinned to ce8c364f, where line 306 is not `CI result`; pin h11 and h16 to 95b68da, cpp-audit to 1035690; fix section 13's "written against" claim. Idea not taken: the reviewer's seventh mention, in cpp-audit.md's list of citations left unpinned, stays as written, since pinning it would contradict that list.

**Code review, round 3, 2026-10-13.** Range `origin/master..c148ed63`, 128 counted lines, @reviewer: CHANGES REQUESTED. The pins at h11:451, h16:1444 and :1505 (`95b68da`, lines 306-307 are `CI result`) and cpp-audit:552 (`1035690`) hold, and section 13 is right; agreed on leaving cpp-audit:664 unpinned. Blocker: the Status line still says "next code review round 2". check_citations exits 0 both ways; CI has not run (no PR yet). The blocker is a one-line wording fix the reviewer named, applied by the main session with the suggestion ("(unpinned then)" at h16:1444); by Ola's ruling of 2026-10-13 on such fixes, the branch counts as reviewed.

**PR 2 code review, round 1, 2026-10-13.** Range `origin/master..b1b43199`, 0 net counted lines, @reviewer: APPROVED WITH WORDINGS (catchment.py:16 "downstream along the chain"; nve.py LAKE_ENVELOPE_HALF "cells up to 33 m"; gauge.py:113 prefix comment; gauge.py:33 U-floor reason). Applied by @developer in `c1e41b77`.
