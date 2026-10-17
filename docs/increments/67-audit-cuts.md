# Increment 67: the code-audit cuts

**Status: designed** (`@architect`, 2026-10-17, on `df923398`). Three pull
requests. P1 (Python and tests) and P2 (C++) touch no governed file and may run
while Ola is away; P3 changes the gates and rule files and waits for his
prompt. Each PR goes through the full path: design, red, green, test run,
review.

## 1. What Ola asked

Ola, 2026-10-17: "yes to the audit cuts tonight". The audit is section 3.2 of
`docs/retrospectives/2026-10-17-overdue-check-night-merges-audit.md`
(branch `worktree-retro-1017`, commit `0e1616a2`): eight cuts, A to H. His
standing rules apply: delete unneeded code; no aliases on a rename; rule files
hold rules and checks, no history.

## 2. Prior art: legacy and literature

*Literature.* Two named refactorings, nothing new claimed. Cut H is *Extract
Function* with *Replace Magic Literal* (Fowler, *Refactoring*, 2nd ed., 2018,
chapters 6 and 9): the lambda becomes named functions and the coded 0/1/2
becomes an enum, with behaviour unchanged. Cut G follows the flag-argument
advice (Martin, *Clean Code*, 2008, ch. 3, "Flag Arguments") as ruff's
`flake8-boolean-trap` rules apply it: the call site names what each boolean
means. No departure from either.

*Legacy.* Nothing to carry: the legacy tree has none of the subjects.

```
$ git grep -lnE 'as_written|written_at|plural\(|check_code_history' legacy-archive -- legacy
(no output, exit 1)
$ git grep -c 'docs/increments' legacy-archive -- legacy
(no output, exit 0)
```

## 3. Each cut, checked against `df923398`

Every count below was rerun on `df923398`; where it differs from the audit
(made on `e5b9b044`), this file's number is the one to use.

| Cut | Checked with | Found | PR |
|---|---|---|---|
| A, history in the 14 modules outside the history gate | `python3 tools/check_code_history.py src_python/tin_engine` | 116 hits, 12 modules (`feature_input.py` 42, `stats.py` 16, `run_record.py` 13, `chains.py` 9, `palettes.py` 8, `features.py` 8, `landcover.py` 7, `land_cover_fractions.py` 6, `resources.py`, `geodesy.py`, `border_simplify.py` 2 each, `legends.py` 1). `source_report.py` and `toml_data.py` are clean. The covered modules give 0 | P1 |
| B, history in C++ comments | the gate's own `RULES` and pointer pattern applied to `//` comments in `include/`, `bindings/`, `src/`; and the audit's regex, `git grep -nP '//.*(\b[0-9]{2}[a-z]?:\|\b[RNDQB][0-9]+\|step [0-9]\|docs/increments/)'` | the gate's rules: 190 lines in 43 files (109 design IDs, 76 pointers, 50 labels, 29 increment numbers). The audit's regex: 160 lines. Its extra hits are IDs the gate's pattern lacks (`N1` to `N19`, `L5`, `I1`, `FE6`), design "step N" references and line-leading `65:`, `34:` labels | P2 |
| C, the rename alias | `src_python/tin_engine/cli.py:66-68`, `tests/python/test_plugin_api.py:35,49-52` | 3 aliases (`_resources`, `_destination`, `_words`) with 18 call sites in `cli.py`, pinned by one parametrised test. Nothing else uses the old names | P1 |
| D, absence tests | `git grep -nE 'assert not hasattr\(\|find_spec\("tin_engine'` and the audit's list | more than the audit's 12; the full list is section 4.4 | P1 |
| E, a word built by hand | `src_python/tin_engine/stats.py:176-177` | 1; `run_record.plural` gives the same words | P1 |
| F, nine copies of a `git` runner | the audit's list | all in `tools/` and `.claude/hooks/`, which are governed | P3 |
| G, positional booleans | `.venv/bin/ruff check --select FBT003 .` | 39: 8 in `src_python`, 16 in `tests`, 6 in `tools`, 1 in `.claude/hooks`, 8 in `docs/increments/*-probes`. Also `source_report.py:175` multiplies a tuple by a boolean, which no rule flags | P1 (src_python, tests); P3 (the rest, and the rule) |
| H, the 78-line lambda | `include/terrain/refinement/refine_points.hpp:397-474` (`as_written`) | 1 | P2 |

## 4. P1: Python and tests (may run while Ola is away)

Files: `src_python/tin_engine/` and `tests/python/` only. No C++, so no build.

### 4.1 Red step (`@tester`)

1. **A's pin**, a new harness-marked test in `tests/python/` (it imports
   nothing from `tin_engine`): `tools/check_code_history.py src_python/tin_engine`
   exits 0 and prints nothing. Red today: 116 lines, exit 1. P3 makes the gate
   read the whole package and deletes this test.
2. **E's pin**: in `tests/python/test_layering.py:34`, the `stats` row becomes
   `"run_record"`. Both are layer 0, so the edge is allowed. Red today: `stats`
   imports nothing. The wording is already pinned by
   `tests/python/test_stats.py:225` ("0 holes") and `:397-415` ("1 hole",
   "2 holes").
3. **C**: delete `MOVED` and
   `test_cli_calls_the_moved_helpers_under_their_old_names`
   (`test_plugin_api.py:35,49-52`). No replacement test: a test that the old
   names are gone is an absence test, the kind D deletes.
4. **D**: the deletions in 4.4, in **a commit of their own**, so that one can
   be reverted if Ola wants the scanner of increment 58 kept (section 7).
5. **G in tests**: keyword arguments at the 16 FBT003 hits in `tests/python/`,
   for example `_engine(xy, chains, delaunay=True, spacing=SNAP)`. At
   `test_station_set.py:412`, `_typed(False)` passes the value under test, so
   it becomes `_typed(value=False)`.

### 4.2 Green step (`@developer`)

- **A**: rewrite the 116 hits by the rule in the `python-development` skill
  (*Docstrings and Comments*): say what the code does and why, with no labels.
  A pointer to a design file stays only where the derivation is too long to
  repeat: at most one per module docstring, in the form `` `<path>.md`, *<heading>* ``.
  Prose only, in a commit of its own made first. Check, at that commit:
  `python3 tools/check_code_history.py --same-code <that commit>^`, which must
  print nothing.
- **C**: `from tin_engine import cli_kit` replaces the three aliased imports,
  and the call sites become `cli_kit.cli_resources(...)`,
  `cli_kit.destination(...)` and `cli_kit.words(...)`. The calls are
  qualified, not imported by their bare names, because `cli.py` has local
  variables named `words` (`:1260`, `:1737`, `:1856`) that would shadow the
  imported function.
- **E**: `stats.py` imports `plural` from `run_record` and writes
  `f"... (1 ring, {plural(holes, 'hole', 'holes')})"`. The local `plural`
  variable goes.
- **G**: keyword arguments at the 8 hits in `src_python`
  (`catchment.py:250,359`, `chains.py:53`, `feature_input.py:515,526`,
  `mosaic.py:377`, `source_report.py:257,329`), and an `if` in place of the
  tuple-times-boolean at `source_report.py:175`. No signature changes: making
  the parameters keyword-only (ruff `FBT001`/`FBT002`, 31 hits) was not part
  of the audit.
  `.venv/bin/ruff check --select FBT003 src_python tests` must print nothing.

### 4.3 What guards P1

The full `pytest` (including the 85 % coverage floor), `mypy`, which reports a
missed rename as an undefined name and a misspelt keyword as an unknown
argument, and `ruff`. By cut: A, the pin and `--same-code`; C, the
`test_cli_*` suites, which go through every call site's command; E,
`test_stats.py` and the layering row; G, the suites of each touched module
(`test_catchment*.py`, `test_feature_input.py`, `test_mosaic*.py`,
`test_source_fits.py`, `test_chains*.py`).

**No `@perf` run.** P1 changes no refine or mesh code. The `cli.py` edits rename
call sites and change no step of the meshing path. `@reviewer` may rule
otherwise.

### 4.4 D: what goes and what stays

The rule: a test goes if it only checks that a deleted name, module or
spelling is still absent. A test stays if it checks a behaviour a user or a
plugin meets, or a boundary rule that a future addition could break.

Deleted whole:

- `test_cli_help.py:224` `test_the_help_names_no_gallery`, `:271`
  `test_the_drawing_package_is_not_importable`, and `:284`
  `test_the_cli_keeps_no_drawing_name` (5 cases);
- `test_reference.py:546` `test_no_product_source_names_the_old_overlap_fields`;
- `test_domain.py:301` `test_the_must_match_rule_is_gone`;
- `test_io_gml.py:228` `test_a_document_has_no_completeness_flag`;
- `test_resources.py:287` `test_t11_decompose_is_gone`;
- `test_legend_name.py:1-148`: the repository scanner
  (`test_no_old_spelling_outside_the_history_records`, `:106`) and its six
  self-tests (`TestTheRules`, `:113-148`). The file keeps its three refusal
  tests (`:151-169`) and a one-line docstring saying what they check.

The absence lines deleted from tests that check something else as well:

- `test_core_refine.py:186`;
- `test_crs.py:227`, with its docstring sentence;
- `test_feature_input.py:744`, and delete the whole test
  `test_read_lake_polygons_replaces_read_lakes`: what is left checks only that a
  name exists, which the next test calls anyway;
- `test_gauge.py:191-192`, with the test renamed for what is left (the
  `U_FLOOR_CELLS` value);
- `test_io_geojson.py:151-154`, the `domain` import and its loop, and the
  docstring 'Moved here from `domain.py`…' goes;
- `test_station_set.py:244,246`, and the docstring 'Audit PR C, red test 5: no
  alias keeps the old name' goes.

Kept, with the reason:

- the old commands and options are refused with a usage error, which a user's
  old script meets: `test_cli_help.py:215` (`--flat`, `--crs`), `:264`
  (`draw`), `test_mcp_vector_sources.py:325` (`--out-parent`),
  `test_legend_name.py:151-169`;
- a signature plugins call: `test_cli_kit.py:36`;
- boundary rules: `test_core_noding.py:474,555` (bindings left out on
  purpose), `:862` (no composite entry point on `_core`),
  `test_features.py:387-388` (`mask` and `names` stay methods);
- `test_open_core.py`, whole. The audit lists its `:66`, but the file guards
  the line between the public and the private repository against a private
  branch merged by mistake (its module docstring). That can happen at any
  merge; it is not deleted code staying deleted.

## 5. P2: C++ comments and the `as_written` function (may run while Ola is away, after 66 merges)

Files: `include/`, `bindings/`, `src/`, `tests/cpp/`, and the B pin in
`tests/python/`. It starts from master **after increment 66 has merged**.
Increment 66 (`worktree-coincident-vertices`) changes `refine.hpp`,
`constraint_points.hpp`, `strip_scan.hpp` and
`prop_refinement_written_delaunay.cpp`, all of which P2 touches. P2 needs a C++
build and `@perf`, so it runs with no other C++ build or `@perf` run beside it.

### 5.1 B: the comment sweep

- **The rule for C++ comments** is the Python one (the `python-development`
  skill, *Docstrings and Comments*) applied to `//` and `/* */` comments: no
  increment numbers or labels, no design-file IDs (including the `N`, `L`,
  `I` and `FE` series the gate's pattern lacks), no "step N" pointing into a
  design file, no `65:` line labels. A comment that numbers the steps of its
  own algorithm (`intersect.hpp:171-173`) is not history and stays. At most
  one design-file pointer per file, in the file's opening comment, in heading
  form (`` `<path>.md`, *<heading>* ``); `tools/check_citations.py` already
  checks C++ comments.
- **Red** (`@tester`): a harness-marked test in `tests/python/` that imports
  `RULES`, `SCIENTIFIC` and `HEADING_CITATION` from `tools/check_code_history.py`
  and applies them to the comment text of every `.hpp`, `.cpp` and `.h` under
  `include/`, `bindings/` and `src/`. It allows only one heading-form pointer
  per file, in the file's opening comment. It must expect no hits; red today
  on 190 lines. P3 moves this into the gate and deletes the test. The extra IDs
  and the "step N" references are found by the audit's regex in section 3,
  judged by `@reviewer`, not pinned: the gate's pattern would flag legitimate
  Python text (P3).
- **Green** (`@developer`): comments only. `@reviewer` checks that only
  comments changed: with `//` and `/* */` comments stripped by a scratch
  script, base and head are the same text in every changed file.
- `tests/cpp/` comments (394 hits by the audit's regex) are out of scope:
  the rule covers production code, as the Python gate does.

### 5.2 H: `as_written` as named functions

`point_loop`'s lambda (`refine_points.hpp:397-474`) uses no template
parameter of `point_loop`, so it becomes non-template `inline` functions in
`terrain::refinement::detail`, unit-testable without running the loop:

```cpp
enum class WrittenKind : std::uint8_t { Faithful, Candidate, Flat };

struct WrittenPoints {
    std::vector<mesh::MeshVertex> at;  // each vertex as written, (x, -y)
    std::vector<char> faithful;        // 1 where the written point is an exact translate of its grid point
};

struct AsWrittenPass {
    bool changed = false;                    // m was replaced by the legalised mesh
    std::size_t flat = 0;                    // triangles flat or turned over as written
    std::size_t flips = 0;                   // flips of the pass, when changed
    std::optional<std::size_t> left_illegal; // set when the result was refused on the grid; m unchanged
};

[[nodiscard]] WrittenPoints written_points(const mesh::LatticeMesh& m, const IndexedMesh2& start,
                                           const raster::RasterGeometry& g, const mesh::LatticeFrame& frame);
[[nodiscard]] std::vector<WrittenKind> written_kinds(const mesh::LatticeMesh& m, const WrittenPoints& w);
[[nodiscard]] AsWrittenPass legalise_as_written(mesh::LatticeMesh& m, WrittenPoints w,
                                                std::uint32_t frozen_mask, std::vector<char>& touched);
```

The call site keeps today's assignments exactly: `out.flat_as_written = pass.flat`
on every call; `out.illegal_as_written = *pass.left_illegal` only when it is
set; `out.flips` and `out.flips_as_written` grow by `pass.flips` only when
`changed`; `touched` is marked by the legalisation as now, also when the
result is refused. The names in the copied mesh (`ctris`, `cbits`, `cmasks`,
`slot`) become `kept_triangles`, `kept_sides`, `kept_masks` and `loop_slot`;
`sides` becomes a named function `constraint_sides`; `exact_sum` becomes a
named helper. No other behaviour changes.

- **Red** (`@tester`): a unit test file in `tests/cpp/unit/`, registered in
  `tests/cpp/CMakeLists.txt`, that fails to compile until the functions exist.
  Fixture: T3's quad from `tests/cpp/property/prop_refinement_written_delaunay.cpp:310-340`
  (`RasterGeometry{0.0, 30.0, 3.0, 3.0, 11, 11}`; A, B, C collinear as
  written, orientation exactly 0; counter-clockwise on the grid, about
  9.25e-20 as that test states). Fed through `detail::to_lattice`
  (`refine.hpp:176`) as `refine_points` does, the expected result is what T3
  asserts through `refine_points`: one `Flat` triangle from `written_kinds`,
  and `flat == 1` from `legalise_as_written` with the mesh left as given. T3
  runs the whole loop, so `@tester` confirms these values against the current
  code, through a scratch copy of the lambda, before committing them. The T5
  fixture (`:363-380`, multiples of 1/64 near (500000, 6600000)) may serve as
  a second case under the same check.
- **Not invariant-critical**: no mutation round. It is a refactor guarded
  by the suites below.

### 5.3 What guards P2

`cmake --build` with `-Werror` (exit status read before `ctest`); `ctest`,
and among it `prop_refinement_written_delaunay.cpp` T1 to T5; `pytest`, among
it `test_refine_golden.py` (output hashes) and `test_core_refine_points.py`;
and the B pin.

**`@perf` applies** (refine code changes): `tools/bench.py`, the 1 m benchmark
and the thread sweep. Both bench hashes must be unchanged, because the
behaviour is the same bit for bit. Time is within the noise of the previous
run, with the power state recorded.

## 6. P3: the gates and rule files (waits for Ola's prompt)

Every file here is governed except `pyproject.toml` (`tools/check_*` and the
listed tools, `.claude/hooks/`, `.claude/agents/`; see
`.claude/hooks/guard_governance.py@4ad42e45:58-84`). After P1 and P2 the tree is already clean, so each gate goes
in green.

1. **A's gate**: `[tool.check_code_history].paths` becomes `src_python/tin_engine`,
   with no exclusions, so a new module is covered by default. P1's pin test
   is deleted.
2. **B's gate**: `check_code_history.py` also reads `//` and `/* */` comments in
   `.hpp`, `.cpp` and `.h` under `include/`, `bindings/` and `src/`, tests first
   in `test_check_code_history.py`. P2's pin test is deleted. Taken up with it:
   the ideas-file entry *Docstring checker gaps* (a capitalised "Ruling", K
   and W IDs, a bare "§N") and the C++ ID series of 5.1, each with a planted
   case and a check that it does not flag the Python package as it stands.
3. **G's gate**: `FBT003` added to `[tool.ruff.lint] select`, with a
   per-file ignore for `docs/**` (the probe scripts are kept as evidence). The 7
   hits in `tools/` and `.claude/hooks/` get keyword arguments.
4. **F**: one `git` runner in `tools/harness_mode.py`, used by the nine
   copies. Low priority; may be left out.
5. **The reviewer's checklist**: `.claude/agents/reviewer.md` gains one
   line for what no gate sees: rename aliases, absence tests, words built
   from variables, long lambdas with coded values. Rules only.

## 7. Questions and defaults

- *`test_legend_name.py`'s scanner.* When Ola approved increment 58 he said
  "Keep the guard tight"; today he approved the audit, which lists the
  scanner for deletion. Ola ruled (2026-10-17): "We don't keep tests that are
  not relevant anymore." Delete it, in its own commit (4.1 item 4), and
  keep the three refusal tests. Reverting that one commit keeps it.

## 8. Size

Counted lines (`CLAUDE.md` §2; comments, docstrings and `tests/` do not count):

| PR | Estimate | Why |
|---|---|---|
| P1 | about +3 (−2 to +10) | A is prose (0); C −2 imports, call sites rewritten in place; E about 0; G +0 to +8 where keyword calls wrap |
| P2 | about +15 (+5 to +30) | B is comments (0); H adds three signatures, two structs and an enum where the lambda had captures |
| P3 | about +20 (−10 to +50) | the C++ comment reader in `check_code_history.py` +40; F about −30 if done; `pyproject.toml` is not counted |

All three are far under the limit. In the diff, P1 is mostly test deletions
(about −300 lines under `tests/`), and P2 is about 230 comment lines in about
45 files.

## 9. Order and conflicts

- P1 and the in-flight increment 63 (`worktree-fetch-vectors`) both change
  `source_report.py` near lines 257 and 329. Whichever merges second merges
  master in (`tools/merge_master.py`); the conflict is in four call lines.
  P1 also touches `cli.py`, as 63 and 23c-2 (`worktree-23c`) do, but not
  the same lines (`git diff` of 63's `cli.py` has no alias line).
- P2 waits for 66 (section 5). `worktree-23c` changes `refine_points.hpp`
  around the same loop; its own writer merges P2 in when it next merges master.

## Review

- **Design review, round 1, 2026-10-17.** Commit `55acdec4`, 0 counted lines, @reviewer: APPROVED WITH WORDINGS, applied by the main session: the history check runs at the prose-only commit; three test cuts finished (a whole test and two docstrings); `pyproject.toml` is not governed. The scanner question: Ola ruled delete. Idea kept: a later sweep of 2771 history labels in `tests/python`.
