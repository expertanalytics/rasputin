# Increment 31: binary output by default

**Status:** implemented, as built in section 9; next: code review.
Designed by `@architect` on `2764ef71`, design approved in review round 2.
Red `e3fe391a` and `9008811c` (`@tester`, tests only; its four pins ruled
in Rulings, P1-P4), green `29bc5ad3` (`@developer`). Question 1 is still
open with Ola; what was built is its default, both formats binary. The
worktree's `.venv` holds a plain copy of `tin_engine` and a copied `_core`
(no editable install), so each changed `src_python/tin_engine` file is copied
into `.venv/lib/python3.*/site-packages/tin_engine/` before running `pytest`.
One small PR, Python only. Not refine or mesh code, so no `@perf` acceptance
run (section 7).

## 1. What Ola asked

The performance review (PR #208, `docs/increments/perf-audit.md` on branch
`worktree-perf-audit`, finding F6 and question 3) found that the default
`.vtk` is text, written one number at a time in Python: 0.97 s of
`write: encode` on the Rio das Velhas piece, and an estimated several minutes
at basin scale (10⁸ triangles). It asked: keep text as the default and speed
up its writer later, or make binary the default? Its proposed default was to
keep text.

Ola, 2026-10-07: "yes og 1,2, change 3 to default binary." So `rasputin mesh`
writes binary unless asked for text.

This reverses part of increment 13's ruling 8 and its U2 (a)
(`docs/increments/13-bundled-mesh.md`, ruling 8 and U2), which made text the
default "because the user reads the output". Text stays one flag away.

## 2. Prior art: legacy and literature

*Literature.* Nothing to cite beyond the format itself: the legacy VTK format
(Kitware, *The VTK User's Guide*, "VTK File Formats") defines `BINARY` as
big-endian packed values, which `io/vtk_legacy.py` already writes; PLY's
`binary_little_endian` is likewise already written by `io/ply.py`. No new
method, no novelty claim.

*Legacy.* Nothing. The archived tree has no mesh writer of either format:

    $ git grep -l -i -e "vtk" -e "\.ply" -e "binary" legacy-archive -- 'legacy/*'
    (no files; the tag holds 30 files under legacy/)

## 3. How it works today

- `src_python/tin_engine/cli.py@2764ef71:632-634`: one Typer flag pair for both
  formats, `--binary/--ascii`, a `bool` named `binary`, default `False`, help
  "Packed records, or text `head` can read."
- `src_python/tin_engine/cli.py@2764ef71:1022-1052`: the composition root passes it explicitly to
  each writer: `write_vtk(..., binary=binary)` and `write_ply(..., ascii=not
  binary)` for the surface and the edge file.
- The writers (`io/vtk_legacy.py`, `io/ply.py`; #205 moved the encoders into
  `io/`) keep their own keyword defaults, `binary=False` and `ascii=True`. The
  CLI never relies on them.

## 4. The change

**D1. Flip the CLI default; keep the flag pair.** `binary: ... = True`. The
flag stays `--binary/--ascii`: it is already the Typer boolean-pair idiom of
`cli.py`'s other pair, `--delaunay/--no-delaunay`, `--ascii` keeps
every existing script and test that names it working, and `--binary` becomes a
harmless no-op spelling. No `--vtk-format` option: increment 13 ruling 8
rejected a format flag that could contradict the suffix, and an encoding
choice is a yes/no.

**D2. One default for both formats.** The pair covers `.vtk` and `.ply`, so
`.ply` becomes binary by default too. That keeps increment 13's U2 principle,
"one command should not have two defaults", and PLY's text writer has the
same per-number cost. (Question 1 below.)

**D3. Help text in plain words** (Ola's plain-output rule). Proposed:

    help="Binary files, small and fast; --ascii writes text you can read."

Typer prints the default itself: today's `rasputin mesh --help` ends the
line with `[default: ascii]`, and after the flip it reads `[default: binary]`.

The `mesh` command's own description is printed on the same help page, and
says the old default in words (`src_python/tin_engine/cli.py@2764ef71:755`,
"Text is the default for both formats."). That sentence becomes:

    Both formats are binary by default; ``--ascii`` writes text you can read.

(Ruling P1; red test 3 checks the old sentence is gone.)

**D4. The writers' own defaults stay.** `write_vtk(binary=False)` and
`write_ply(ascii=True)` are library defaults the CLI overrides explicitly;
flipping them would churn `test_io_vtk_legacy.py` and `test_io_ply.py` (their
`write()` helpers rely on them) for no user-visible gain. Their module
docstrings get one sentence each saying `rasputin mesh` writes binary unless
`--ascii` (increment 31), so the "ASCII by default" ruling text there is not
read as the command's behaviour.

**D5. Docs.** The `mesh` docstring sentence of D3, and `README.md:82` ("Both formats are text by default; `--binary`
writes packed records.") becomes "Both formats are binary by default; `--ascii`
writes text you can read with `head`." Nothing in `INSTALL.md` or
`project_structure.md` names the default (`git grep -n -i -e ascii -e
binary -- INSTALL.md project_structure.md` finds only `checked_ascii`,
`escaped_ascii` and a binary stream).

## 5. What else reads the output, and what depends on text

Checked on `2764ef71`; each reader reads both encodings or passes the flag
explicitly.

| Reader | Encoding | Effect |
|---|---|---|
| `tests/python/vtkread.py` `read_vtk` (used by `cli_driver.mesh_to_vtk`, `recordread`, `landcover_fixtures`, every CLI suite) | both (`tests/python/vtkread.py@2764ef71:207-209`) | none |
| `tests/python/plyread.py` | both (`tests/python/plyread.py@2764ef71:170-173`) | none |
| `tests/python/test_io_vtk_readback.py` (CI's viewer step, real `vtk`) | passes `--binary` or `--ascii` explicitly | none |
| `tools/bench.py` | timing runs pass `--binary`, the quality run passes `--ascii` and reads it with `read_vtk_ascii` (`tools/bench.py@2764ef71:342`, `:607`, `:616`) | none; stored baselines unchanged |
| Golden digests (`test_refine_golden.digest`, used by `test_cli_constraint_feet.py`) | hash the in-memory refine arrays, not file bytes | none |
| Two byte-for-byte golden tests: `tests/python/test_refine_golden.py` `test_the_kartverket_stride_vtk_is_unchanged_by_node_sampling` (`@2764ef71:202-206`, SHA-256 of the whole file against `GOLDEN_STRIDE_VTK`) and `tests/python/test_cli_mesh_refine.py` `TestWithoutTolerance.test_the_uniform_mesh_is_unchanged` (`@2764ef71:212-219`, SHA-256 from `POINTS` on against `INCREMENT_12_LARGER_STRIDE_2`) | **text, by leaving the flag out** | fail after the flip; the red commit passes `--ascii` in both runs and keeps both hashes (section 6) |
| `tests/python/test_cli_mesh_landcover.py` `test_every_cell_carries_a_code_in_both_encodings` (`@2764ef71:289-303`) | **text, by leaving the flag out** (`["--binary"] if binary else []`) | does not fail, but its `ascii` case would quietly write binary; the red commit passes `--ascii` there (section 6) |
| The one-off probes under `docs/increments/26-probes/` (`docs/increments/26-probes/fractions_probe.py@2764ef71:52-61` `read_mesh`, which `water_probe.py` imports) | both: VTK's own `vtkPolyDataReader` | none |
| `@perf`'s byte-identity probes of increments 30a-30c and the Python audit (`docs/increments/30*-probes/*_bytes.py`, `python-audit-probes/encoder_bytes.py`) and the stats scripts of `docs/benchmarks/2026-10-06/` | pass `--binary` (and the audit also `--ascii`) explicitly | none |
| Ola's converter `../rasputin_data/lands/scripts/vtk_to_gpkg.py` (outside the repository) | both: reads line 3 for `BINARY` | none |
| ParaView, QGIS | both | none |
| The run record and `--stats` | record the command line as typed and no encoding field (`git grep -n -i -e binary -e encoding -- src_python/tin_engine` finds only `cli.py`) | no wording change |

Two stored byte-for-byte hashes depend on the text default, and one test
picks text by leaving the flag out; all three get `--ascii` in the red commit.
The hashes themselves stay: `GOLDEN_STRIDE_VTK` "must not change"
(`tests/python/test_refine_golden.py@2764ef71:170-190`), and a hash of the text
file still pins the same mesh. A future comparison across this increment's
merge (an old commit against a new one) must pass `--binary` or `--ascii`
explicitly, as every script above already does.

**Tests that assert the old default** (to be inverted, section 6):
`tests/python/test_cli_mesh_vtk.py` `TestEncoding.test_vtk_is_ascii_by_default`
and `test_ply_is_ascii_by_default`; `tests/python/test_cli_mesh.py`
`test_both_files_are_ascii_by_default`. Tests that name `--ascii` explicitly
(`TestTheAsciiFlag`, `test_ascii_is_the_explicit_spelling`, `test_cli_fetch.py`)
stay as they are. Their stale wording changes with them: `TestEncoding`'s
docstring (`tests/python/test_cli_mesh_vtk.py@2764ef71:103`, "text by default
(U2 (a))") and the module docstring of `tests/python/test_cli_mesh.py`
(`@2764ef71:3-4`), which describe the old default.

## 6. Tests `@tester` writes red first

1. Invert the three default tests above: a `.vtk` written with no encoding
   flag has `BINARY` on line 3; both `.ply` files are `binary_little_endian`.
   Rename them `..._binary_by_default` and note the reversal of increment 13
   U2 (a) in the comment, as the existing comment does for increment 10.
2. `--ascii` still writes text for both formats (exists; keep).
3. `rasputin mesh --help` shows the default as binary (Typer's
   `[default: binary]` marker on the `--binary / --ascii` line), so the help
   and the behaviour cannot drift apart. Strip ANSI and line-wrap before
   matching, as `test_cli_mesh_geographic.py` does.

4. Amendments in the same red commit, since `@developer` may not touch tests
   (section 5's table):
   - `test_the_kartverket_stride_vtk_is_unchanged_by_node_sampling` and
     `TestWithoutTolerance.test_the_uniform_mesh_is_unchanged` pass `--ascii`;
     their hashes stay as recorded.
   - `test_every_cell_carries_a_code_in_both_encodings` builds its flag as
     `["--binary"] if binary else ["--ascii"]`.
   - The `TestEncoding` docstring and `test_cli_mesh.py`'s module docstring
     say binary by default (increment 31).

Before committing red, `@tester` runs the whole Python suite once against a
local, uncommitted flip of the default in its own worktree (then reverts it).
Expected there: 5 failures, the three default tests of item 1 and the two
golden hashes of item 4; any other failure is a test that silently relied on
text output and is amended in the red commit. A flip run finds only tests
that fail, so `@tester` also greps for tests that choose an encoding by
leaving the flag out and makes each one name its flag:

    git grep -n -e '\["--binary"\] if' -e 'if binary else \[\]' -- tests

On `2764ef71` this finds only `test_cli_mesh_landcover.py:293`. No mutation
round: a one-line default is not an invariant-critical suite.

*Note for the next default change:* grep the command's own docstring too, not
only the option's help and the README; Typer prints the docstring in
`--help`. This design missed `src_python/tin_engine/cli.py@2764ef71:755`, and `@tester` found it
(ruling P1). For this repository:

    git grep -n -i -e 'is the default' -e 'by default' -- src_python README.md

## 7. Cost, `@perf`, LOC

- **Production lines:** about +1 net (`cli.py`: the default and the help
  string; ruff may wrap the option onto one more line). Docstrings and
  `README.md` do not count. Far under the 700 limit.
- **`@perf`: not needed.** The diff touches neither `include/terrain/refinement/`
  nor `include/terrain/mesh/` nor what drives them (`docs/increments/README.md`,
  "Acceptance"); it changes which existing encoder the CLI picks. `bench.py`
  passes the encoding explicitly, so no baseline moves. The speed-up for a
  run that omits the flag is the gap F6 measured between the two encoders,
  not a new claim of this increment.
- **Dependencies:** none added.

## 8. Questions for Ola

1. **Should `.ply` become binary by default too, or only `.vtk`?** One flag
   covers both today. *Default: both, so the command keeps one default
   (increment 13's principle).* A "no" reopens D1, D3 and red test 3: the
   flag would become three-state (`bool | None`, unset meaning "by suffix":
   binary for `.vtk`, text for `.ply`), Typer could no longer print a single
   `[default: ...]` marker, so the help states the two defaults in words and
   test 3 matches that sentence. A few more production lines, still far under
   the limit.
   On the red commit `e3fe391a`, a "no" changes these tests:
   `tests/python/test_cli_mesh_vtk.py@e3fe391a:129-138`
   (`test_ply_is_binary_by_default` goes back to text),
   `tests/python/test_cli_mesh.py@e3fe391a:107-112`
   (`test_both_files_are_binary_by_default`, the `.ply` pair, goes back to
   text), and `tests/python/test_cli_mesh_vtk.py@e3fe391a:140-154`
   (`test_the_help_shows_binary_as_the_default` matches the two-default
   sentence instead of one `[default: ...]` marker).
   *Still open with Ola after the green step.* What was built is the
   default, one flag with both formats binary (section 9); a "no" now also
   reopens the green code, `src_python/tin_engine/cli.py@29bc5ad3:632-638`
   and the `mesh` description at `:759-760`.

## 9. As built

Green `29bc5ad3`, on red `e3fe391a` and its citation fix `9008811c`.

- **D1, D2.** `src_python/tin_engine/cli.py@29bc5ad3:632-638`: the
  `--binary/--ascii` pair now defaults to `True`, for `.vtk` and both `.ply`
  files alike. The composition root is unchanged; it already passes the flag
  to each writer.
- **D3.** The help text is D3's sentence, and `COLUMNS=100 rasputin mesh
  --help` ends the option's entry with `[default: binary]`. The `mesh`
  description says Ruling P1's sentence, wrapped over two docstring lines
  (`src_python/tin_engine/cli.py@29bc5ad3:759-760`).
- **D4.** The writers keep `binary=False` and `ascii=True`; one sentence each
  says `rasputin mesh` writes binary unless given `--ascii`
  (`src_python/tin_engine/io/vtk_legacy.py@29bc5ad3:14-15`,
  `src_python/tin_engine/io/ply.py@29bc5ad3:14-15`).
- **D5.** `README.md`'s output section says both formats are binary by
  default and `--ascii` writes text you can read with `head`.
- **Production lines: 4 net** (6 added, 2 removed;
  `python3 tools/count_loc.py 2764ef71 29bc5ad3`), against the estimate of
  about 1. The overrun is formatting: with the longer help string the
  option's inner line would be 118 characters, over ruff's 100, so
  `ruff format` spreads the `binary` parameter over seven lines where it had
  three. No logic beyond the flipped default.
- **Tests:** the whole suite passes on `29bc5ad3`
  (`.venv/bin/python -m pytest -q`: 5690 passed, 17 skipped).

## Rulings

Pins of the red step `e3fe391a` (`@tester`), ruled by `@architect`,
2026-10-07. All four stand.

- **P1. The `mesh` description loses "Text is the default".** Upheld. The
  help test (`tests/python/test_cli_mesh_vtk.py@e3fe391a:140-154`) also asserts
  that sentence (`src_python/tin_engine/cli.py@2764ef71:755`) is gone from
  `--help`; D3 and D5 missed it. The new wording is set in D3:
  "Both formats are binary by default; ``--ascii`` writes text you can read."
  `@developer` writes exactly that sentence. A docstring line, so no
  production-line cost.
- **P2. The help test sets `COLUMNS=100`.** Upheld: without it the width of
  the caller's terminal decides whether Rich cuts the option names, and the
  test would pass or fail by terminal. 100 columns is a test pin, not a product
  setting.
- **P3. The `.vtk` default test also checks that the raw third line is
  `BINARY`,** beside `read_vtk`'s encoding. Upheld: it checks the bytes
  directly, as section 6 item 1 states it, so it does not depend on the reader.
- **P4. A third stale docstring updated:** the module docstring of
  `tests/python/test_cli_mesh_vtk.py@e3fe391a:1-7` (lines 3-7 describe the
  pair's default). Upheld; section 5 should have listed it with the other two.

## Review

Design review round 1 (@reviewer, 2026-10-07): CHANGES REQUESTED on 2764ef71..ad91b5dd (docs only, 0 net production lines by tools/count_loc.py): two byte-for-byte golden tests hash the default text output, contrary to section 5 (/Users/skavhaug/projects/rasputin/tests/python/test_refine_golden.py@2764ef71:185, /Users/skavhaug/projects/rasputin/tests/python/test_cli_mesh_refine.py@2764ef71:134), and the ascii case of /Users/skavhaug/projects/rasputin/tests/python/test_cli_mesh_landcover.py@2764ef71:293 would quietly become a second binary case.

Design review round 2 (@reviewer, 2026-10-07): APPROVED on ad91b5dd..09cbee96 (docs only, 0 net production lines by tools/count_loc.py): round 1's two findings and four suggestions are answered in sections 5, 6 and 8; round 1's first two citations were wrong and are corrected here: the golden tests are at /Users/skavhaug/projects/rasputin/tests/python/test_refine_golden.py@2764ef71:202-206 and /Users/skavhaug/projects/rasputin/tests/python/test_cli_mesh_refine.py@2764ef71:212-219, and the land-cover citation /Users/skavhaug/projects/rasputin/tests/python/test_cli_mesh_landcover.py@2764ef71:293 stands.
